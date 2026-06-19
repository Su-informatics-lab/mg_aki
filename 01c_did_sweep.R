#!/usr/bin/env Rscript
# ============================================================================
# 01c_did_sweep.R — Specification curve: sIPTW_DR + AIPW, ICU-time anchor
#
# Method: sIPTW_DR (stabilized IPTW, trim 1/99, + covariate adjustment)
#   ICU-time anchor: Cr_pre = first postop Cr, Cr_post = Cr at target_h from ICU
#   No temporal alignment needed → full sample contributes
#   Weighted SMDs reported (what the love plot should show)
#
# Toggles (8 groups, 256 specs):
#   D1: chronic drugs    D2: steroids    D3: ICU drugs    D4: beta_blockers
#   L1: K+               L2: Ca          L3: lactate      L4: serum Mg
#
# Time points: 24h, 36h (from ICU admission)
# Databases: eICU (with cluster SE), MIMIC
#
# Run: Rscript 01c_did_sweep.R
# ============================================================================

suppressPackageStartupMessages({ library(sandwich); library(lmtest) })

RESULTS <- path.expand("~/mg_aki/results")
TIME_POINTS <- c(24, 36)

# ── Fixed base (15 covariates, always included) ──────────────────────────
BASE <- c("age", "is_female", "bmi",
          "surg_cabg", "surg_valve", "surg_combined",
          "heart_failure", "hypertension", "diabetes", "ckd",
          "copd", "pvd", "stroke", "liver_disease",
          "egfr", "first_heartrate")

# ── Toggleable groups ────────────────────────────────────────────────────
TOGGLES <- list(
  D1_chronic  = c("ppi_chronic", "loop_diuretic_chronic", "acei_arb_chronic", "nsaid_chronic"),
  D2_steroids = c("steroids"),
  D3_icu_drug = c("loop_diuretics", "antiarrhythmics"),
  D4_beta     = c("beta_blockers"),
  L1_potass   = c("first_potassium"),
  L2_calcium  = c("first_calcium"),
  L3_lactate  = c("first_lactate", "lactate_missing"),
  L4_mg       = c("first_mg_value")
)

# ── Helpers ──────────────────────────────────────────────────────────────
median_impute <- function(d, vars) {
  for (v in vars) if (v %in% names(d) && any(is.na(d[[v]])))
    d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm=TRUE)
  d
}

wsmd <- function(x, trt, w) {
  m1 <- weighted.mean(x[trt==1], w[trt==1], na.rm=TRUE)
  m0 <- weighted.mean(x[trt==0], w[trt==0], na.rm=TRUE)
  sp <- sqrt((var(x[trt==1],na.rm=T) + var(x[trt==0],na.rm=T))/2)
  if (is.na(sp)||sp<1e-10) 0 else abs(m1-m0)/sp
}

# ── sIPTW_DR estimator ───────────────────────────────────────────────────
run_siptw_dr <- function(d, ps_vars, outcome_col="delta_cr",
                          trt_col="treated", cluster_col=NULL) {
  # PS model
  avail <- intersect(ps_vars, names(d))
  d <- median_impute(d, avail)
  d <- d[complete.cases(d[,c(avail, outcome_col, trt_col)]),]
  if (nrow(d)<40) return(NULL)
  n_trt <- sum(d[[trt_col]]==1); n_ctl <- sum(d[[trt_col]]==0)
  if (n_trt<20 || n_ctl<20) return(NULL)

  ps_fml <- as.formula(paste(trt_col, "~", paste(avail,collapse="+")))
  ps_fit <- tryCatch(glm(ps_fml, data=d, family=binomial()), error=function(e) NULL)
  if (is.null(ps_fit)) return(NULL)

  d$ps <- pmax(pmin(fitted(ps_fit), 0.99), 0.01)
  prev <- mean(d[[trt_col]])

  # Stabilized IPTW
  d$siptw <- ifelse(d[[trt_col]]==1, prev/d$ps, (1-prev)/(1-d$ps))
  q01 <- quantile(d$siptw, 0.01); q99 <- quantile(d$siptw, 0.99)
  d$siptw <- pmax(pmin(d$siptw, q99), q01)

  # Weighted SMDs
  wsmds <- sapply(avail, function(v)
    if(is.numeric(d[[v]])) wsmd(d[[v]], d[[trt_col]], d$siptw) else 0)
  max_wsmd <- max(wsmds, na.rm=TRUE)
  n_above <- sum(wsmds > 0.1, na.rm=TRUE)

  # DR: adjust for covariates with weighted SMD > 0.05
  adj <- names(wsmds[wsmds > 0.05])
  adj <- adj[vapply(adj, function(v) var(d[[v]],na.rm=T)>1e-10, logical(1))]

  if (length(adj)>0) {
    fml <- as.formula(paste(outcome_col, "~", trt_col, "+", paste(adj,collapse="+")))
  } else {
    fml <- as.formula(paste(outcome_col, "~", trt_col))
  }

  fit <- tryCatch(lm(fml, data=d, weights=siptw), error=function(e) NULL)
  if (is.null(fit)) return(NULL)

  vc <- if (!is.null(cluster_col) && cluster_col %in% names(d) &&
             length(unique(d[[cluster_col]]))>1)
    tryCatch(vcovCL(fit, cluster=d[[cluster_col]]),
             error=function(e) vcovHC(fit,type="HC1"))
  else vcovHC(fit, type="HC1")

  ct <- coeftest(fit, vcov.=vc)
  trt_r <- which(rownames(ct)==trt_col)
  if (length(trt_r)==0) return(NULL)

  data.frame(n_trt=n_trt, n_ctl=n_ctl, max_wsmd=round(max_wsmd,3),
             n_above_01=n_above,
             did=round(ct[trt_r,"Estimate"],5),
             se=round(ct[trt_r,"Std. Error"],5),
             p=round(ct[trt_r,"Pr(>|t|)"],4),
             stringsAsFactors=F)
}

# ── AIPW estimator ───────────────────────────────────────────────────────
run_aipw <- function(d, ps_vars, outcome_col="delta_cr",
                      trt_col="treated", cluster_col=NULL) {
  avail <- intersect(ps_vars, names(d))
  d <- median_impute(d, avail)
  d <- d[complete.cases(d[,c(avail, outcome_col, trt_col)]),]
  if (nrow(d)<40) return(NULL)
  n_trt <- sum(d[[trt_col]]==1); n_ctl <- sum(d[[trt_col]]==0)
  if (n_trt<20 || n_ctl<20) return(NULL)

  # PS
  ps_fml <- as.formula(paste(trt_col, "~", paste(avail,collapse="+")))
  ps_fit <- tryCatch(glm(ps_fml, data=d, family=binomial()), error=function(e) NULL)
  if (is.null(ps_fit)) return(NULL)
  d$ps <- pmax(pmin(fitted(ps_fit), 0.99), 0.01)

  # Outcome models (separate for treated/control)
  out_fml <- as.formula(paste(outcome_col, "~", paste(avail,collapse="+")))
  d1 <- d[d[[trt_col]]==1,]; d0 <- d[d[[trt_col]]==0,]
  m1 <- tryCatch(lm(out_fml, data=d1), error=function(e) NULL)
  m0 <- tryCatch(lm(out_fml, data=d0), error=function(e) NULL)
  if (is.null(m1)||is.null(m0)) return(NULL)

  mu1 <- predict(m1, newdata=d); mu0 <- predict(m0, newdata=d)
  Y <- d[[outcome_col]]; T_ <- d[[trt_col]]; e <- d$ps

  # AIPW scores
  phi <- T_*(Y-mu1)/e - (1-T_)*(Y-mu0)/(1-e) + (mu1-mu0)
  tau <- mean(phi)
  se <- sd(phi)/sqrt(length(phi))

  data.frame(n_trt=n_trt, n_ctl=n_ctl,
             did=round(tau,5), se=round(se,5),
             p=round(2*pnorm(-abs(tau/se)),4),
             stringsAsFactors=F)
}


# ============================================================================
# LOAD DATA
# ============================================================================
cat("======================================================================\n")
cat("01c_did_sweep.R — sIPTW_DR + AIPW specification curve\n")
cat(sprintf("  Base: %d fixed | Toggles: %d groups -> %d specs\n",
            length(BASE), length(TOGGLES), 2^length(TOGGLES)))
cat(sprintf("  Time points: %s | Methods: sIPTW_DR, AIPW\n",
            paste(TIME_POINTS, "h", sep="", collapse=", ")))
cat("======================================================================\n")

load_db <- function(tag) {
  trt <- read.csv(file.path(RESULTS,sprintf("did_treated_%s.csv",tag)),stringsAsFactors=F)
  ctl <- read.csv(file.path(RESULTS,sprintf("did_control_%s.csv",tag)),stringsAsFactors=F)
  cr_all <- read.csv(file.path(RESULTS,sprintf("did_cr_all_%s.csv",tag)),stringsAsFactors=F)

  id_col <- if("patientunitstayid" %in% names(trt)) "patientunitstayid" else "stay_id"
  trt$pid <- trt[[id_col]]; ctl$pid <- ctl[[id_col]]
  cr_id <- if("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
  cr_all$pid <- cr_all[[cr_id]]
  if(!"labresultoffset" %in% names(cr_all)) cr_all$labresultoffset <- cr_all$offset_min
  cr_all$offset_h <- cr_all$labresultoffset / 60

  # Merge all covariates
  all_covars <- unique(c(BASE, unlist(TOGGLES)))
  shared_cols <- intersect(c("pid","treated",all_covars), intersect(names(trt),names(ctl)))
  combined <- rbind(trt[,shared_cols], ctl[,shared_cols])

  # Cluster column for eICU
  cluster_col <- NULL
  if ("hospitalid" %in% names(trt)) {
    hosp_map <- setNames(trt$hospitalid, trt$pid)
    hosp_map_c <- setNames(ctl$hospitalid, ctl$pid)
    combined$hospitalid <- c(hosp_map[as.character(trt$pid)],
                              hosp_map_c[as.character(ctl$pid)])
    cluster_col <- "hospitalid"
  }

  list(combined=combined, cr_all=cr_all, cluster_col=cluster_col)
}

# Build DiD dataset at target_h from ICU admission
build_did_icu <- function(combined, cr_all, target_h, window_h=6) {
  # Cr_pre: first Cr within 0-6h of ICU
  cr_pre <- cr_all[cr_all$offset_h >= 0 & cr_all$offset_h <= 6,]
  cr_pre <- cr_pre[order(cr_pre$pid, cr_pre$offset_h),]
  cr_pre <- cr_pre[!duplicated(cr_pre$pid),]

  # Cr_post: closest Cr to target_h, within +/-window_h
  cr_post <- cr_all[cr_all$offset_h >= (target_h-window_h) &
                     cr_all$offset_h <= (target_h+window_h),]
  cr_post$dist <- abs(cr_post$offset_h - target_h)
  cr_post <- cr_post[order(cr_post$pid, cr_post$dist),]
  cr_post <- cr_post[!duplicated(cr_post$pid),]

  # Merge
  m <- merge(cr_pre[,c("pid","labresult")],
             cr_post[,c("pid","labresult")],
             by="pid", suffixes=c("_pre","_post"))
  m$delta_cr <- m$labresult_post - m$labresult_pre

  # Join with covariates
  d <- merge(combined, m[,c("pid","delta_cr")], by="pid")
  d
}

cat("\nLoading eICU...\n")
dat_e <- load_db("eicu")
cat(sprintf("  eICU combined: %d patients\n", nrow(dat_e$combined)))
cat("Loading MIMIC...\n")
dat_m <- load_db("mimic")
cat(sprintf("  MIMIC combined: %d patients\n", nrow(dat_m$combined)))

# ============================================================================
# SWEEP
# ============================================================================
n_toggles <- length(TOGGLES)
n_specs <- 2^n_toggles
toggle_names <- names(TOGGLES)
total_runs <- n_specs * length(TIME_POINTS)

cat(sprintf("\nRunning %d specs x %d time points x 2 methods x 2 databases = %d total...\n\n",
            n_specs, length(TIME_POINTS), n_specs * length(TIME_POINTS) * 2 * 2))

all_rows <- list(); ridx <- 0

for (i in 0:(n_specs-1)) {
  bits <- as.integer(intToBits(i)[1:n_toggles])
  on_names <- toggle_names[bits == 1]

  ps_vars <- BASE
  for (tn in on_names) ps_vars <- c(ps_vars, TOGGLES[[tn]])
  ps_vars <- unique(ps_vars)

  label <- if(length(on_names)==0) "base_only" else paste(on_names, collapse="+")
  spec_id <- i + 1

  for (th in TIME_POINTS) {
    ridx <- ridx + 1

    # Build DiD datasets at this time point
    d_e <- build_did_icu(dat_e$combined, dat_e$cr_all, th)
    d_m <- build_did_icu(dat_m$combined, dat_m$cr_all, th)

    # sIPTW_DR
    r_e <- tryCatch(run_siptw_dr(d_e, ps_vars, cluster_col=dat_e$cluster_col),
                    error=function(e) NULL)
    r_m <- tryCatch(run_siptw_dr(d_m, ps_vars, cluster_col=dat_m$cluster_col),
                    error=function(e) NULL)

    # AIPW
    a_e <- tryCatch(run_aipw(d_e, ps_vars), error=function(e) NULL)
    a_m <- tryCatch(run_aipw(d_m, ps_vars), error=function(e) NULL)

    row <- data.frame(
      spec_id=spec_id, target_h=th, n_covars=length(ps_vars), label=label,
      D1=bits[1], D2=bits[2], D3=bits[3], D4=bits[4],
      L1=bits[5], L2=bits[6], L3=bits[7], L4=bits[8],
      # sIPTW_DR results
      e_n=if(!is.null(r_e)) r_e$n_trt else NA,
      e_wsmd=if(!is.null(r_e)) r_e$max_wsmd else NA,
      e_did=if(!is.null(r_e)) r_e$did else NA,
      e_p=if(!is.null(r_e)) r_e$p else NA,
      m_n=if(!is.null(r_m)) r_m$n_trt else NA,
      m_wsmd=if(!is.null(r_m)) r_m$max_wsmd else NA,
      m_did=if(!is.null(r_m)) r_m$did else NA,
      m_p=if(!is.null(r_m)) r_m$p else NA,
      # AIPW results
      e_aipw=if(!is.null(a_e)) a_e$did else NA,
      e_aipw_p=if(!is.null(a_e)) a_e$p else NA,
      m_aipw=if(!is.null(a_m)) a_m$did else NA,
      m_aipw_p=if(!is.null(a_m)) a_m$p else NA,
      stringsAsFactors=F
    )

    # Direction flags (sIPTW_DR)
    if (!is.null(r_e) && !is.null(r_m)) {
      row$concordant <- (r_e$did < 0 & r_m$did < 0) | (r_e$did > 0 & r_m$did > 0)
      row$both_neg <- (r_e$did < 0 & r_m$did < 0)
      row$both_sig_neg <- (r_e$did<0 & r_e$p<0.05 & r_m$did<0 & r_m$p<0.05)
    } else {
      row$concordant <- row$both_neg <- row$both_sig_neg <- NA
    }

    all_rows[[ridx]] <- row

    # Progress (print every 16 specs)
    if (spec_id %% 16 == 1 || spec_id == n_specs) {
      de <- if(!is.null(r_e)) ifelse(r_e$did<0,"-","+") else "?"
      dm <- if(!is.null(r_m)) ifelse(r_m$did<0,"-","+") else "?"
      se <- if(!is.null(r_e)&&r_e$p<0.05) "*" else " "
      sm <- if(!is.null(r_m)&&r_m$p<0.05) "*" else " "
      cat(sprintf("  [%3d/%d] %dh %2dcov e:%s%s m:%s%s %s\n",
                  spec_id, n_specs, th, length(ps_vars), de,se, dm,sm, label))
    }
  }
}

sweep <- do.call(rbind, all_rows)
write.csv(sweep, file.path(RESULTS, "did_sweep.csv"), row.names=FALSE)

# ============================================================================
# SUMMARY
# ============================================================================
SEP <- paste(rep("=",70),collapse="")
cat(sprintf("\n%s\nSWEEP COMPLETE: %d rows (%d specs x %d time points)\n%s\n",
            SEP, nrow(sweep), n_specs, length(TIME_POINTS), SEP))

for (th in TIME_POINTS) {
  s <- sweep[sweep$target_h == th,]
  cat(sprintf("\n  --- %dh from ICU (sIPTW_DR) ---\n", th))
  n_conc <- sum(s$concordant, na.rm=T)
  n_neg <- sum(s$both_neg, na.rm=T)
  n_bsn <- sum(s$both_sig_neg, na.rm=T)
  n_tot <- sum(!is.na(s$concordant))
  cat(sprintf("  Concordant: %d/%d (%.0f%%)  Both negative: %d  Both sig neg: %d\n",
              n_conc, n_tot, 100*n_conc/max(n_tot,1), n_neg, n_bsn))

  if (n_neg > 0) {
    cat("\n  Both NEGATIVE (sorted by MIMIC P):\n")
    cat("  spec covars  eICU_DiD  eICU_P  eICU_wSMD  MIMIC_DiD  MIMIC_P  MIMIC_wSMD  toggles\n")
    neg <- s[!is.na(s$both_neg) & s$both_neg,]
    neg <- neg[order(neg$m_p),]
    for (j in seq_len(min(nrow(neg),20))) {
      r <- neg[j,]
      se <- if(r$e_p<0.05) "*" else " "
      sm <- if(r$m_p<0.05) "*" else " "
      cat(sprintf("  %3d  %4d   %+.4f %s %.4f  %.3f    %+.4f %s %.4f  %.3f    %s\n",
                  r$spec_id, r$n_covars, r$e_did, se, r$e_p, r$e_wsmd,
                  r$m_did, sm, r$m_p, r$m_wsmd, r$label))
    }
  }
}

# Toggle impact analysis
cat(sprintf("\n  --- TOGGLE IMPACT (24h, sIPTW_DR) ---\n"))
s24 <- sweep[sweep$target_h == 24,]
cat("  Toggle       ON_mean_e  ON_mean_m  OFF_mean_e OFF_mean_m  e_flip  m_flip\n")
for (k in seq_along(toggle_names)) {
  col <- paste0(c("D","D","D","D","L","L","L","L")[k], c(1,2,3,4,1,2,3,4)[k])
  on <- s24[[col]] == 1; off <- s24[[col]] == 0
  me_on <- mean(s24$e_did[on], na.rm=T); me_off <- mean(s24$e_did[off], na.rm=T)
  mm_on <- mean(s24$m_did[on], na.rm=T); mm_off <- mean(s24$m_did[off], na.rm=T)
  e_flip <- (me_on > 0) != (me_off > 0)
  m_flip <- (mm_on > 0) != (mm_off > 0)
  cat(sprintf("  %-12s %+.4f   %+.4f   %+.4f   %+.4f     %s       %s\n",
              toggle_names[k], me_on, mm_on, me_off, mm_off,
              if(e_flip) "FLIP" else "    ", if(m_flip) "FLIP" else "    "))
}

cat(sprintf("\n  Saved: did_sweep.csv (%d rows)\n", nrow(sweep)))
