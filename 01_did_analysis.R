#!/usr/bin/env Rscript
# ============================================================================
# 01_did_analysis.R — AIPW primary analysis, ICU-time anchor
#
# Method: AIPW (primary), sIPTW_DR (reported alongside)
# Anchor: ICU admission time (no temporal alignment needed)
# Models: primary|sens_a|sens_b|sens_c|sens_d|original (via did_covars.R)
#
# Usage: Rscript 01_did_analysis.R <db> [model]
#   Rscript 01_did_analysis.R eicu             # primary model
#   Rscript 01_did_analysis.R mimic sens_c     # sensitivity C
# ============================================================================

suppressPackageStartupMessages({ library(sandwich); library(lmtest) })
source(file.path(path.expand("~/mg_aki"), "did_covars.R"))

RESULTS <- path.expand("~/mg_aki/results")
TARGETS <- seq(6, 36, by = 3)
PRIMARY_TARGET <- 24
CR_WINDOW <- 6  # +/-6h for Cr_post selection

MG_BINS <- list("<1.8"=c(0,1.8), "1.8-2.0"=c(1.8,2.0),
                "2.0-2.3"=c(2.0,2.3), ">2.3"=c(2.3,99))
SURG_TYPES <- c("cabg","valve","combined","other_cardiac")
SECONDARY_OUTCOMES <- c("hosp_mortality","poaf","encephalopathy","vent_arrhythmia")
OUTCOME_LABELS <- c(hosp_mortality="Hospital mortality", poaf="POAF (new-onset)",
                     encephalopathy="Encephalopathy", vent_arrhythmia="Ventricular arrhythmia")

# ── Helpers ──────────────────────────────────────────────────────────────
median_impute <- function(d, vars) {
  for (v in vars) if (v %in% names(d) && any(is.na(d[[v]])))
    d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm=TRUE)
  d
}

wsmd <- function(x, trt, w) {
  m1 <- weighted.mean(x[trt==1], w[trt==1], na.rm=TRUE)
  m0 <- weighted.mean(x[trt==0], w[trt==0], na.rm=TRUE)
  sp <- sqrt((var(x[trt==1],na.rm=T)+var(x[trt==0],na.rm=T))/2)
  if (is.na(sp)||sp<1e-10) 0 else abs(m1-m0)/sp
}

# ── AIPW estimator (continuous outcome) ──────────────────────────────────
run_aipw <- function(d, ps_vars, outcome="delta_cr", trt="treated") {
  avail <- intersect(ps_vars, names(d))
  d <- d[complete.cases(d[, c(avail, outcome, trt)]), ]
  n <- nrow(d); nt <- sum(d[[trt]]==1); nc <- n - nt
  if (nt<20 || nc<20) return(NULL)

  # PS model
  ps_fml <- as.formula(paste(trt, "~", paste(avail, collapse="+")))
  ps_fit <- tryCatch(glm(ps_fml, data=d, family=binomial()), error=function(e) NULL)
  if (is.null(ps_fit)) return(NULL)
  e <- pmax(pmin(fitted(ps_fit), 0.99), 0.01)

  # Outcome models (separate for T=1 and T=0)
  out_fml <- as.formula(paste(outcome, "~", paste(avail, collapse="+")))
  d1 <- d[d[[trt]]==1,]; d0 <- d[d[[trt]]==0,]
  m1 <- tryCatch(lm(out_fml, data=d1), error=function(e) NULL)
  m0 <- tryCatch(lm(out_fml, data=d0), error=function(e) NULL)
  if (is.null(m1)||is.null(m0)) return(NULL)

  mu1 <- predict(m1, newdata=d); mu0 <- predict(m0, newdata=d)
  Y <- d[[outcome]]; T_ <- d[[trt]]

  # AIPW influence function
  phi <- (mu1 - mu0) + T_*(Y - mu1)/e - (1-T_)*(Y - mu0)/(1-e)
  tau <- mean(phi)
  se <- sd(phi) / sqrt(n)
  p <- 2 * pnorm(-abs(tau/se))

  # Weighted SMDs (sIPTW weights for balance reporting)
  prev <- mean(T_)
  w <- ifelse(T_==1, prev/e, (1-prev)/(1-e))
  q01 <- quantile(w, 0.01); q99 <- quantile(w, 0.99)
  w <- pmax(pmin(w, q99), q01)
  wsmds <- sapply(avail, function(v)
    if(is.numeric(d[[v]])) wsmd(d[[v]], T_, w) else 0)

  list(method="AIPW", n_trt=nt, n_ctl=nc,
       did=tau, se=se, p=p,
       ci_lo=tau-1.96*se, ci_hi=tau+1.96*se,
       max_wsmd=max(wsmds, na.rm=TRUE),
       n_above_01=sum(wsmds>0.1, na.rm=TRUE),
       wsmds=wsmds)
}

# ── sIPTW_DR estimator (continuous outcome) ──────────────────────────────
run_siptw <- function(d, ps_vars, outcome="delta_cr", trt="treated",
                       cluster_col=NULL) {
  avail <- intersect(ps_vars, names(d))
  d <- d[complete.cases(d[, c(avail, outcome, trt)]), ]
  n <- nrow(d); nt <- sum(d[[trt]]==1); nc <- n - nt
  if (nt<20 || nc<20) return(NULL)

  ps_fml <- as.formula(paste(trt, "~", paste(avail, collapse="+")))
  ps_fit <- tryCatch(glm(ps_fml, data=d, family=binomial()), error=function(e) NULL)
  if (is.null(ps_fit)) return(NULL)
  d$ps <- pmax(pmin(fitted(ps_fit), 0.99), 0.01)
  prev <- mean(d[[trt]])
  d$w <- ifelse(d[[trt]]==1, prev/d$ps, (1-prev)/(1-d$ps))
  q01 <- quantile(d$w, 0.01); q99 <- quantile(d$w, 0.99)
  d$w <- pmax(pmin(d$w, q99), q01)

  wsmds <- sapply(avail, function(v)
    if(is.numeric(d[[v]])) wsmd(d[[v]], d[[trt]], d$w) else 0)
  adj <- names(wsmds[wsmds > 0.05])
  adj <- adj[vapply(adj, function(v) var(d[[v]],na.rm=T)>1e-10, logical(1))]

  if (length(adj)>0)
    fml <- as.formula(paste(outcome, "~", trt, "+", paste(adj,collapse="+")))
  else fml <- as.formula(paste(outcome, "~", trt))

  fit <- tryCatch(lm(fml, data=d, weights=w), error=function(e) NULL)
  if (is.null(fit)) return(NULL)

  vc <- if (!is.null(cluster_col) && cluster_col %in% names(d) &&
             length(unique(d[[cluster_col]]))>1)
    tryCatch(vcovCL(fit, cluster=d[[cluster_col]]),
             error=function(e) vcovHC(fit,type="HC1"))
  else vcovHC(fit, type="HC1")
  ct <- coeftest(fit, vcov.=vc)
  tr <- which(rownames(ct)==trt)
  if (length(tr)==0) return(NULL)

  list(method="sIPTW_DR", n_trt=nt, n_ctl=nc,
       did=ct[tr,"Estimate"], se=ct[tr,"Std. Error"],
       p=ct[tr,"Pr(>|t|)"],
       ci_lo=ct[tr,"Estimate"]-1.96*ct[tr,"Std. Error"],
       ci_hi=ct[tr,"Estimate"]+1.96*ct[tr,"Std. Error"],
       max_wsmd=max(wsmds,na.rm=TRUE), n_above_01=sum(wsmds>0.1,na.rm=TRUE))
}

# ── sIPTW_DR for binary outcomes (returns OR) ────────────────────────────
run_binary_siptw <- function(d, ps_vars, outcome_col, trt="treated",
                              cluster_col=NULL) {
  if (!outcome_col %in% names(d)) return(NULL)
  d$y <- as.numeric(d[[outcome_col]])
  d <- d[!is.na(d$y),]
  avail <- intersect(ps_vars, names(d))
  d <- d[complete.cases(d[,avail]),]
  nt <- sum(d[[trt]]==1); nc <- sum(d[[trt]]==0)
  if (nt<20||nc<20) return(NULL)
  r1 <- mean(d$y[d[[trt]]==1]); r0 <- mean(d$y[d[[trt]]==0])
  if (r1==0 && r0==0) return(NULL)

  ps_fml <- as.formula(paste(trt, "~", paste(avail, collapse="+")))
  ps_fit <- tryCatch(glm(ps_fml, data=d, family=binomial()), error=function(e) NULL)
  if (is.null(ps_fit)) return(NULL)
  d$ps <- pmax(pmin(fitted(ps_fit), 0.99), 0.01)
  prev <- mean(d[[trt]])
  d$w <- ifelse(d[[trt]]==1, prev/d$ps, (1-prev)/(1-d$ps))
  q01 <- quantile(d$w, 0.01); q99 <- quantile(d$w, 0.99)
  d$w <- pmax(pmin(d$w, q99), q01)

  wsmds <- sapply(avail, function(v)
    if(is.numeric(d[[v]])) wsmd(d[[v]], d[[trt]], d$w) else 0)
  adj <- names(wsmds[wsmds > 0.05])
  adj <- adj[vapply(adj, function(v) var(d[[v]],na.rm=T)>1e-10, logical(1))]

  if (length(adj)>0)
    fml <- as.formula(paste("y ~", trt, "+", paste(adj,collapse="+")))
  else fml <- as.formula(paste("y ~", trt))

  fit <- tryCatch(glm(fml, data=d, family=quasibinomial(), weights=w),
                  error=function(e) NULL)
  if (is.null(fit) || !trt %in% names(coef(fit))) return(NULL)

  vc <- if (!is.null(cluster_col) && cluster_col %in% names(d) &&
             length(unique(d[[cluster_col]]))>1)
    tryCatch(vcovCL(fit, cluster=d[[cluster_col]]),
             error=function(e) vcovHC(fit,type="HC1"))
  else vcovHC(fit, type="HC1")
  ct <- tryCatch(coeftest(fit, vcov.=vc), error=function(e) NULL)
  if (is.null(ct)) return(NULL)
  tr <- which(rownames(ct)==trt)
  if (length(tr)==0) return(NULL)
  p_col <- grep("^Pr", colnames(ct))

  list(outcome=outcome_col, n_trt=nt, n_ctl=nc,
       rate_trt=round(r1,4), rate_ctl=round(r0,4),
       or=round(exp(ct[tr,"Estimate"]),3),
       or_lo=round(exp(ct[tr,"Estimate"]-1.96*ct[tr,"Std. Error"]),3),
       or_hi=round(exp(ct[tr,"Estimate"]+1.96*ct[tr,"Std. Error"]),3),
       p=round(ct[tr,p_col],4))
}

# ── Build ICU-time-anchor DiD dataset ────────────────────────────────────
build_did_icu <- function(combined, cr_all, target_h, window=CR_WINDOW) {
  # Cr_pre: first Cr within 0-6h of ICU
  pre <- cr_all[cr_all$offset_h >= 0 & cr_all$offset_h <= 6,]
  pre <- pre[order(pre$pid, pre$offset_h),]
  pre <- pre[!duplicated(pre$pid),]

  # Cr_post: closest to target_h within +/-window, must be after Cr_pre
  post <- cr_all[cr_all$offset_h >= (target_h-window) &
                  cr_all$offset_h <= (target_h+window),]
  post$dist <- abs(post$offset_h - target_h)
  post <- post[order(post$pid, post$dist),]
  post <- post[!duplicated(post$pid),]

  # Merge, ensure post > pre
  m <- merge(pre[,c("pid","labresult","offset_h")],
             post[,c("pid","labresult","offset_h")],
             by="pid", suffixes=c("_pre","_post"))
  m <- m[m$offset_h_post > m$offset_h_pre,]
  m$delta_cr <- m$labresult_post - m$labresult_pre

  merge(combined, m[,c("pid","delta_cr")], by="pid")
}


# ============================================================================
# MAIN
# ============================================================================
run_analysis <- function(db, model_name, ps_vars) {
  tag <- tolower(db)
  mtag <- tolower(model_name)
  ftag <- sprintf("%s_%s", tag, mtag)
  SEP <- paste(rep("=",70),collapse="")

  cat(sprintf("\n%s\n%s [%s]: AIPW + sIPTW_DR, ICU-time anchor, %d covariates\n%s\n",
              SEP, db, mtag, length(ps_vars), SEP))

  # ── Load ────────────────────────────────────────────────────────────────
  trt <- read.csv(file.path(RESULTS,sprintf("did_treated_%s.csv",tag)),stringsAsFactors=F)
  ctl <- read.csv(file.path(RESULTS,sprintf("did_control_%s.csv",tag)),stringsAsFactors=F)
  cr_all <- read.csv(file.path(RESULTS,sprintf("did_cr_all_%s.csv",tag)),stringsAsFactors=F)

  id_col <- if("patientunitstayid" %in% names(trt)) "patientunitstayid" else "stay_id"
  trt$pid <- trt[[id_col]]; ctl$pid <- ctl[[id_col]]
  cr_id <- if("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
  cr_all$pid <- cr_all[[cr_id]]
  if(!"labresultoffset" %in% names(cr_all)) cr_all$labresultoffset <- cr_all$offset_min
  cr_all$offset_h <- cr_all$labresultoffset / 60

  ps_avail <- intersect(ps_vars, intersect(names(trt), names(ctl)))
  cat(sprintf("  PS: %d requested, %d available\n", length(ps_vars), length(ps_avail)))
  miss <- setdiff(ps_vars, names(trt))
  if (length(miss)>0) cat(sprintf("  WARNING missing: %s\n", paste(miss,collapse=", ")))

  # Stack covariates
  all_want <- unique(c("pid","treated", ps_avail, "surgery_type",
                        "hosp_mortality","poaf","encephalopathy","vent_arrhythmia",
                        "prior_af","first_mg_value"))
  sh <- intersect(all_want, intersect(names(trt), names(ctl)))
  combined <- rbind(trt[,sh], ctl[,sh])
  combined <- median_impute(combined, ps_avail)

  cluster_col <- NULL
  if ("hospitalid" %in% names(trt)) {
    hmap <- c(setNames(trt$hospitalid, trt$pid), setNames(ctl$hospitalid, ctl$pid))
    combined$hospitalid <- hmap[as.character(combined$pid)]
    cluster_col <- "hospitalid"
  }

  cat(sprintf("  Loaded: %d treated + %d control\n",
              sum(combined$treated==1), sum(combined$treated==0)))

  # ── A. Time course ──────────────────────────────────────────────────────
  cat(sprintf("\n  TIME COURSE (6-36h, +/-%dh window, AIPW + sIPTW_DR)\n", CR_WINDOW))
  cat("  target  n_trt  n_ctl  wSMD   AIPW_DiD   AIPW_P   sIPTW_DiD  sIPTW_P\n")
  cat("  ------  -----  -----  -----  ---------  -------  ---------  -------\n")

  tc_rows <- list()
  for (i in seq_along(TARGETS)) {
    th <- TARGETS[i]
    d <- build_did_icu(combined, cr_all, th)
    if (nrow(d)<40) next

    ra <- run_aipw(d, ps_avail)
    rs <- run_siptw(d, ps_avail, cluster_col=cluster_col)

    if (!is.null(ra)) {
      sig_a <- if(ra$p<0.05) "*" else " "
      sig_s <- if(!is.null(rs) && rs$p<0.05) "*" else " "
      primary <- if(th==PRIMARY_TARGET) " << PRIMARY" else ""
      cat(sprintf("  %4dh  %5d  %5d  %.3f  %+.4f %s %.4f  %+.4f %s %.4f%s\n",
                  th, ra$n_trt, ra$n_ctl, ra$max_wsmd,
                  ra$did, sig_a, ra$p,
                  if(!is.null(rs)) rs$did else NA, sig_s,
                  if(!is.null(rs)) rs$p else NA, primary))
    }

    tc_rows[[i]] <- data.frame(
      target_h=th,
      n_trt=if(!is.null(ra)) ra$n_trt else NA,
      n_ctl=if(!is.null(ra)) ra$n_ctl else NA,
      max_wsmd=if(!is.null(ra)) round(ra$max_wsmd,3) else NA,
      aipw_did=if(!is.null(ra)) round(ra$did,5) else NA,
      aipw_se=if(!is.null(ra)) round(ra$se,5) else NA,
      aipw_p=if(!is.null(ra)) round(ra$p,4) else NA,
      aipw_lo=if(!is.null(ra)) round(ra$ci_lo,5) else NA,
      aipw_hi=if(!is.null(ra)) round(ra$ci_hi,5) else NA,
      siptw_did=if(!is.null(rs)) round(rs$did,5) else NA,
      siptw_p=if(!is.null(rs)) round(rs$p,4) else NA,
      stringsAsFactors=F)
  }
  tc <- do.call(rbind, tc_rows)
  write.csv(tc, file.path(RESULTS,sprintf("did_timecourse_%s.csv",ftag)), row.names=F)

  # ── B. Primary (24h) ───────────────────────────────────────────────────
  cat(sprintf("\n  PRIMARY (%dh +/-%dh)\n", PRIMARY_TARGET, CR_WINDOW))
  d24 <- build_did_icu(combined, cr_all, PRIMARY_TARGET)
  ra24 <- run_aipw(d24, ps_avail)
  rs24 <- run_siptw(d24, ps_avail, cluster_col=cluster_col)

  if (!is.null(ra24)) {
    cat(sprintf("  AIPW:     DiD=%+.4f (SE=%.4f), P=%.4f, 95%%CI [%+.4f, %+.4f]\n",
                ra24$did, ra24$se, ra24$p, ra24$ci_lo, ra24$ci_hi))
    cat(sprintf("            n=%d/%d, max wSMD=%.3f, %d covars>0.1\n",
                ra24$n_trt, ra24$n_ctl, ra24$max_wsmd, ra24$n_above_01))

    # Print worst wSMDs
    bad <- sort(ra24$wsmds[ra24$wsmds > 0.05], decreasing=TRUE)
    if (length(bad)>0) {
      cat("            Weighted SMD > 0.05:\n")
      for (nm in names(bad)) cat(sprintf("              %-20s %.3f\n", nm, bad[nm]))
    }
  }
  if (!is.null(rs24)) {
    cat(sprintf("  sIPTW_DR: DiD=%+.4f (SE=%.4f), P=%.4f, 95%%CI [%+.4f, %+.4f]\n",
                rs24$did, rs24$se, rs24$p, rs24$ci_lo, rs24$ci_hi))
  }

  # Save primary dataset for downstream
  write.csv(d24, file.path(RESULTS,sprintf("did_primary_%s.csv",ftag)), row.names=F)

  # ── C. Subgroups (AIPW, on primary dataset) ────────────────────────────
  cat(sprintf("\n  SUBGROUPS (AIPW, %dh)\n", PRIMARY_TARGET))
  sg_rows <- list(); sidx <- 0

  cat("  Mg strata:\n")
  if ("first_mg_value" %in% names(d24)) {
    for (mg_nm in names(MG_BINS)) {
      lo <- MG_BINS[[mg_nm]][1]; hi <- MG_BINS[[mg_nm]][2]
      sub <- d24[!is.na(d24$first_mg_value) &
                  d24$first_mg_value >= lo & d24$first_mg_value < hi, ]
      res <- run_aipw(sub, ps_avail)
      sidx <- sidx+1
      if (!is.null(res)) {
        sig <- if(res$p<0.05) " *" else ""
        cat(sprintf("    %s: n=%d/%d, DiD=%+.4f, P=%.4f%s\n",
                    mg_nm, res$n_trt, res$n_ctl, res$did, res$p, sig))
        sg_rows[[sidx]] <- data.frame(subgroup="Mg_stratum", stratum=mg_nm,
          n_trt=res$n_trt, n_ctl=res$n_ctl, did=round(res$did,4),
          p=round(res$p,4), stringsAsFactors=F)
      } else {
        cat(sprintf("    %s: too few\n", mg_nm))
      }
    }
  }

  # Mg interaction
  if ("first_mg_value" %in% names(d24) && "first_mg_value" %in% ps_avail) {
    int_avail <- setdiff(ps_avail, "first_mg_value")
    int_avail <- int_avail[vapply(int_avail, function(v) var(d24[[v]],na.rm=T)>1e-10, logical(1))]
    int_fml <- as.formula(paste("delta_cr ~ treated * first_mg_value +",
                                 paste(int_avail, collapse="+")))
    int_fit <- tryCatch(lm(int_fml, data=d24), error=function(e) NULL)
    if (!is.null(int_fit)) {
      int_ct <- coeftest(int_fit, vcov.=vcovHC(int_fit, type="HC1"))
      ir <- grep("treated:first_mg_value", rownames(int_ct))
      if (length(ir)>0) cat(sprintf("  Interaction (trt x Mg): b=%+.4f, P=%.4f\n",
                                     int_ct[ir,"Estimate"], int_ct[ir,"Pr(>|t|)"]))
    }
  }

  cat("\n  Surgery type:\n")
  if ("surgery_type" %in% names(d24)) {
    for (st in SURG_TYPES) {
      sub <- d24[!is.na(d24$surgery_type) & d24$surgery_type==st, ]
      res <- run_aipw(sub, ps_avail)
      sidx <- sidx+1
      if (!is.null(res)) {
        sig <- if(res$p<0.05) " *" else ""
        cat(sprintf("    %s: n=%d/%d, DiD=%+.4f, P=%.4f%s\n",
                    st, res$n_trt, res$n_ctl, res$did, res$p, sig))
        sg_rows[[sidx]] <- data.frame(subgroup="Surgery", stratum=st,
          n_trt=res$n_trt, n_ctl=res$n_ctl, did=round(res$did,4),
          p=round(res$p,4), stringsAsFactors=F)
      } else cat(sprintf("    %s: too few\n", st))
    }
  }

  if (length(sg_rows)>0) {
    sg_df <- do.call(rbind, sg_rows)
    write.csv(sg_df, file.path(RESULTS,sprintf("did_subgroups_%s.csv",ftag)), row.names=F)
  }

  # ── D. Secondary outcomes (sIPTW_DR, OR) ───────────────────────────────
  cat(sprintf("\n  SECONDARY OUTCOMES (sIPTW_DR, full sample)\n"))
  cat("  outcome                n_trt  n_ctl  rate_trt rate_ctl  OR (95%CI)        P\n")
  cat("  --------------------- ------ ------ -------- -------- ----------------- ------\n")

  sec_rows <- list()
  for (oc in SECONDARY_OUTCOMES) {
    res <- run_binary_siptw(combined, ps_avail, oc, cluster_col=cluster_col)
    if (is.null(res)) {
      label <- ifelse(oc %in% names(OUTCOME_LABELS), OUTCOME_LABELS[oc], oc)
      cat(sprintf("  %-23s  -- unavailable\n", label)); next
    }
    label <- ifelse(oc %in% names(OUTCOME_LABELS), OUTCOME_LABELS[oc], oc)
    sig <- if(res$p<0.05) " *" else ""
    cat(sprintf("  %-23s %5d  %5d   %.1f%%   %.1f%%  %.2f (%.2f-%.2f)  %.4f%s\n",
                label, res$n_trt, res$n_ctl, 100*res$rate_trt, 100*res$rate_ctl,
                res$or, res$or_lo, res$or_hi, res$p, sig))
    sec_rows[[oc]] <- as.data.frame(res, stringsAsFactors=F)
  }
  if (length(sec_rows)>0) {
    sec_df <- do.call(rbind, sec_rows)
    write.csv(sec_df, file.path(RESULTS,sprintf("did_secondary_%s.csv",ftag)), row.names=F)
  }

  # ── Summary ─────────────────────────────────────────────────────────────
  cat(sprintf("\n%s\n%s [%s]: COMPLETE\n", SEP, db, mtag))
  cat(sprintf("  did_timecourse_%s.csv  | did_primary_%s.csv\n", ftag, ftag))
  cat(sprintf("  did_subgroups_%s.csv   | did_secondary_%s.csv\n", ftag, ftag))
}

# ============================================================================
# CLI
# ============================================================================
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  cat("Usage: Rscript 01_did_analysis.R <db> [model]\n")
  cat("  db:    eicu | mimic\n")
  cat("  model: primary (default) | sens_a | sens_b | sens_c | sens_d | original\n")
  quit(status=1)
}

db_arg <- toupper(args[1])
model_arg <- if(length(args)>=2) tolower(args[2]) else "primary"

cat("======================================================================\n")
cat(sprintf("01_did_analysis.R — AIPW primary, model=%s\n", model_arg))
ps_vars <- select_model(model_arg)
cat("======================================================================\n")

run_analysis(db_arg, model_arg, ps_vars)
