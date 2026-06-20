#!/usr/bin/env Rscript
# ============================================================================
# 01_did_analysis.R — PSM primary + AIPW sensitivity
#
# Primary: PSM (1:1 and 1:4), ICU-time anchor, 21 covariates
# Sensitivity: AIPW (reported alongside)
# No sIPTW, no OW.
#
# Usage: Rscript 01_did_analysis.R eicu [model]
#        Rscript 01_did_analysis.R mimic
# ============================================================================

suppressPackageStartupMessages({
  library(MatchIt); library(sandwich); library(lmtest)
})
source(file.path(path.expand("~/mg_aki"), "did_covars.R"))

RESULTS    <- path.expand("~/mg_aki/results")
CALIPER    <- 0.2
TARGETS    <- seq(6, 36, by = 3)
PRIMARY_T  <- 24
WINDOW     <- 6

MG_BINS <- list("<1.8"=c(0,1.8),"1.8-2.0"=c(1.8,2.0),
                "2.0-2.3"=c(2.0,2.3),">2.3"=c(2.3,99))
SURG_TYPES <- c("cabg","valve","combined","other_cardiac")
SEC_OC <- c("hosp_mortality","poaf","encephalopathy","vent_arrhythmia")
SEC_LB <- c(hosp_mortality="Hospital mortality",poaf="POAF (new-onset)",
            encephalopathy="Encephalopathy",vent_arrhythmia="Ventricular arrhythmia")

# ── Helpers ──────────────────────────────────────────────────────────────
median_impute <- function(d, vars) {
  for (v in vars) if (v %in% names(d) && any(is.na(d[[v]])))
    d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm=TRUE)
  d
}

compute_smds <- function(d, vars) {
  sapply(vars, function(v) {
    if (!v %in% names(d)) return(NA)
    x1 <- d[[v]][d$treated==1]; x0 <- d[[v]][d$treated==0]
    sp <- sqrt((var(x1,na.rm=T)+var(x0,na.rm=T))/2)
    if (is.na(sp)||sp<1e-10) NA else abs(mean(x1,na.rm=T)-mean(x0,na.rm=T))/sp
  })
}

# Doubly-robust DiD in matched set
did_dr <- function(df, ps_vars, pair_col="subclass") {
  nt <- sum(df$treated==1); nc <- sum(df$treated==0)
  if (nt<10||nc<10) return(NULL)
  smds <- compute_smds(df, ps_vars)
  adj <- names(smds[!is.na(smds) & smds>0.05])
  adj <- intersect(adj, names(df))
  adj <- adj[vapply(adj, function(v) var(df[[v]],na.rm=T)>1e-10, logical(1))]

  get_vc <- function(fit) {
    if (pair_col %in% names(df) && length(unique(df[[pair_col]]))>1)
      tryCatch(vcovCL(fit, cluster=df[[pair_col]]),
               error=function(e) vcovHC(fit,type="HC1"))
    else vcovHC(fit,type="HC1")
  }

  fit0 <- lm(delta_cr ~ treated, data=df)
  ct0 <- coeftest(fit0, vcov.=get_vc(fit0))
  if (length(adj)>0) {
    fml <- as.formula(paste("delta_cr ~ treated +", paste(adj,collapse="+")))
    fit1 <- tryCatch(lm(fml,data=df), error=function(e) NULL)
    if (!is.null(fit1)) ct1 <- coeftest(fit1, vcov.=get_vc(fit1))
    else ct1 <- ct0
  } else ct1 <- ct0

  list(n_trt=nt, n_ctl=nc,
       did_unadj=ct0["treated","Estimate"], p_unadj=ct0["treated","Pr(>|t|)"],
       did_adj=ct1["treated","Estimate"], se_adj=ct1["treated","Std. Error"],
       p_adj=ct1["treated","Pr(>|t|)"],
       ci_lo=ct1["treated","Estimate"]-1.96*ct1["treated","Std. Error"],
       ci_hi=ct1["treated","Estimate"]+1.96*ct1["treated","Std. Error"],
       n_adjust=length(adj))
}

# AIPW (for sensitivity)
run_aipw <- function(d, ps_vars, outcome="delta_cr") {
  avail <- intersect(ps_vars, names(d))
  d <- d[complete.cases(d[, c(avail, outcome, "treated")]), ]
  nt <- sum(d$treated==1); nc <- sum(d$treated==0)
  if (nt<15||nc<15) return(NULL)
  ps_fml <- as.formula(paste("treated ~", paste(avail,collapse="+")))
  ps_fit <- tryCatch(glm(ps_fml, data=d, family=binomial()), error=function(e) NULL)
  if (is.null(ps_fit)) return(NULL)
  e <- pmax(pmin(fitted(ps_fit),0.99),0.01)
  out_fml <- as.formula(paste(outcome,"~",paste(avail,collapse="+")))
  m1 <- tryCatch(lm(out_fml,data=d[d$treated==1,]),error=function(e) NULL)
  m0 <- tryCatch(lm(out_fml,data=d[d$treated==0,]),error=function(e) NULL)
  if (is.null(m1)||is.null(m0)) return(NULL)
  mu1 <- predict(m1,newdata=d); mu0 <- predict(m0,newdata=d)
  Y <- d[[outcome]]; T_ <- d$treated
  phi <- (mu1-mu0)+T_*(Y-mu1)/e-(1-T_)*(Y-mu0)/(1-e)
  tau <- mean(phi); se <- sd(phi)/sqrt(nrow(d))
  list(did=tau, se=se, p=2*pnorm(-abs(tau/se)),
       ci_lo=tau-1.96*se, ci_hi=tau+1.96*se, n_trt=nt, n_ctl=nc)
}

# Build ICU-time ΔCr dataset
build_dcr <- function(combined, cr_all, target_h) {
  pre <- cr_all[cr_all$offset_h>=0 & cr_all$offset_h<=6,]
  pre <- pre[order(pre$pid,pre$offset_h),]; pre <- pre[!duplicated(pre$pid),]
  post <- cr_all[cr_all$offset_h>=(target_h-WINDOW) & cr_all$offset_h<=(target_h+WINDOW),]
  post$dist <- abs(post$offset_h - target_h)
  post <- post[order(post$pid,post$dist),]; post <- post[!duplicated(post$pid),]
  m <- merge(pre[,c("pid","labresult","offset_h")],
             post[,c("pid","labresult","offset_h")],
             by="pid",suffixes=c("_pre","_post"))
  m <- m[m$offset_h_post > m$offset_h_pre,]
  m$cr_pre <- m$labresult_pre; m$cr_post <- m$labresult_post
  m$delta_cr <- m$cr_post - m$cr_pre
  m$aki <- as.integer(m$delta_cr >= 0.3)
  merge(combined, m[,c("pid","cr_pre","cr_post","delta_cr","aki")], by="pid")
}


# ============================================================================
run_analysis <- function(db, model_name, ps_vars) {
  tag <- tolower(db); mtag <- tolower(model_name)
  ftag <- sprintf("%s_%s", tag, mtag)
  SEP <- paste(rep("=",70),collapse="")

  cat(sprintf("\n%s\n%s [%s]: PSM + AIPW, ICU-time anchor, %d covariates\n%s\n",
              SEP, db, mtag, length(ps_vars), SEP))

  # ── Load ────────────────────────────────────────────────────────────
  trt <- read.csv(file.path(RESULTS,sprintf("did_treated_%s.csv",tag)),stringsAsFactors=F)
  ctl <- read.csv(file.path(RESULTS,sprintf("did_control_%s.csv",tag)),stringsAsFactors=F)
  cr_all <- read.csv(file.path(RESULTS,sprintf("did_cr_all_%s.csv",tag)),stringsAsFactors=F)

  id_col <- if("patientunitstayid" %in% names(trt)) "patientunitstayid" else "stay_id"
  trt$pid <- trt[[id_col]]; ctl$pid <- ctl[[id_col]]
  cr_id <- if("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
  cr_all$pid <- cr_all[[cr_id]]
  if(!"labresultoffset" %in% names(cr_all)) cr_all$labresultoffset <- cr_all$offset_min
  cr_all$offset_h <- cr_all$labresultoffset / 60

  all_want <- unique(c("pid","treated",ps_vars,"surgery_type","first_mg_value",
                        "hosp_mortality","poaf","encephalopathy","vent_arrhythmia"))
  sh <- intersect(all_want, intersect(names(trt),names(ctl)))
  combined <- rbind(trt[,sh], ctl[,sh])
  combined <- median_impute(combined, ps_vars)
  ps_avail <- intersect(ps_vars, names(combined))

  cat(sprintf("  Loaded: %d treated + %d control\n",
              sum(combined$treated==1), sum(combined$treated==0)))

  # ── A. PSM (1:1 and 1:4) ───────────────────────────────────────────
  for (ratio in c(1, 4)) {
    cat(sprintf("\n  ── PSM 1:%d (caliper=%.1f, replace=T) ──\n", ratio, CALIPER))
    ps_fml <- as.formula(paste("treated ~", paste(ps_avail,collapse="+")))
    m <- suppressWarnings(matchit(ps_fml, data=combined, method="nearest",
                                   distance="glm", ratio=ratio,
                                   caliper=CALIPER, replace=TRUE))
    md <- match.data(m)
    smds_raw <- compute_smds(combined, ps_avail)
    smds_match <- compute_smds(md, ps_avail)

    cat(sprintf("  Raw max SMD: %.3f | Matched max SMD: %.3f | n>0.1: %d/%d\n",
                max(smds_raw,na.rm=T), max(smds_match,na.rm=T),
                sum(smds_match>0.1,na.rm=T), length(ps_avail)))
    worst <- sort(smds_match[smds_match>0.1], decreasing=TRUE)
    if (length(worst)>0) {
      cat("  Covariates with matched SMD > 0.1:\n")
      for (nm in names(worst)) cat(sprintf("    %-20s %.3f\n", nm, worst[nm]))
    }

    # Save balance for love plot
    bal <- data.frame(covariate=names(smds_raw),
                       raw_smd=round(as.numeric(smds_raw),4),
                       matched_smd=round(as.numeric(smds_match),4),
                       stringsAsFactors=F)
    write.csv(bal, file.path(RESULTS,sprintf("did_balance_%s_%d.csv",ftag,ratio)), row.names=F)

    n_trt_matched <- sum(md$treated==1); n_ctl_matched <- sum(md$treated==0)
    cat(sprintf("  Matched: %d treated, %d control\n", n_trt_matched, n_ctl_matched))

    # ── Time course (PSM, ICU-time anchor) ────────────────────────────
    cat(sprintf("\n  TIME COURSE (1:%d, ICU-time anchor)\n", ratio))
    cat("  target  n_trt  n_ctl  PSM_DiD     P       AIPW_DiD    P\n")
    cat("  ------  -----  -----  --------  ------  ----------  ------\n")

    tc_rows <- list()
    for (ti in seq_along(TARGETS)) {
      th <- TARGETS[ti]
      # Build ΔCr for matched patients
      dcr <- build_dcr(md, cr_all, th)
      if (nrow(dcr)<40) next

      res_psm <- did_dr(dcr, ps_avail)
      res_aipw <- run_aipw(dcr, ps_avail)  # AIPW on matched set

      if (!is.null(res_psm)) {
        sig_p <- if(res_psm$p_adj<0.05) "*" else " "
        sig_a <- if(!is.null(res_aipw) && res_aipw$p<0.05) "*" else " "
        primary <- if(th==PRIMARY_T) " << PRIMARY" else ""
        cat(sprintf("  %4dh  %5d  %5d  %+.4f %s %.4f  %+.4f %s %.4f%s\n",
                    th, res_psm$n_trt, res_psm$n_ctl,
                    res_psm$did_adj, sig_p, res_psm$p_adj,
                    if(!is.null(res_aipw)) res_aipw$did else NA, sig_a,
                    if(!is.null(res_aipw)) res_aipw$p else NA, primary))
      }

      tc_rows[[ti]] <- data.frame(
        target_h=th, ratio=ratio,
        psm_n_trt=if(!is.null(res_psm)) res_psm$n_trt else NA,
        psm_n_ctl=if(!is.null(res_psm)) res_psm$n_ctl else NA,
        psm_did=if(!is.null(res_psm)) round(res_psm$did_adj,5) else NA,
        psm_se=if(!is.null(res_psm)) round(res_psm$se_adj,5) else NA,
        psm_p=if(!is.null(res_psm)) round(res_psm$p_adj,4) else NA,
        psm_lo=if(!is.null(res_psm)) round(res_psm$ci_lo,5) else NA,
        psm_hi=if(!is.null(res_psm)) round(res_psm$ci_hi,5) else NA,
        aipw_did=if(!is.null(res_aipw)) round(res_aipw$did,5) else NA,
        aipw_p=if(!is.null(res_aipw)) round(res_aipw$p,4) else NA,
        stringsAsFactors=F)
    }
    tc <- do.call(rbind, tc_rows)
    write.csv(tc, file.path(RESULTS,sprintf("did_timecourse_%s_r%d.csv",ftag,ratio)), row.names=F)

    # ── Primary 24h detail ────────────────────────────────────────────
    if (ratio == 1) {
      cat(sprintf("\n  PRIMARY (%dh, 1:1 matched)\n", PRIMARY_T))
      d24 <- build_dcr(md, cr_all, PRIMARY_T)
      res24 <- did_dr(d24, ps_avail)
      aipw24 <- run_aipw(d24, ps_avail)

      if (!is.null(res24)) {
        cat(sprintf("  PSM-DR:  DiD=%+.4f (SE=%.4f), P=%.4f, 95%%CI [%+.4f,%+.4f]\n",
                    res24$did_adj, res24$se_adj, res24$p_adj, res24$ci_lo, res24$ci_hi))
        cat(sprintf("           n=%d/%d, %d covars adjusted\n",
                    res24$n_trt, res24$n_ctl, res24$n_adjust))
      }
      if (!is.null(aipw24)) {
        cat(sprintf("  AIPW:    DiD=%+.4f (SE=%.4f), P=%.4f, 95%%CI [%+.4f,%+.4f]\n",
                    aipw24$did, aipw24$se, aipw24$p, aipw24$ci_lo, aipw24$ci_hi))
      }

      # Save matched set
      write.csv(d24, file.path(RESULTS,sprintf("did_matched_%s_24h.csv",ftag)), row.names=F)

      # ── Subgroups ──────────────────────────────────────────────────
      cat("\n  SUBGROUPS (PSM-DR, 24h)\n")
      if ("first_mg_value" %in% names(d24)) {
        cat("  Mg strata:\n")
        for (mg_nm in names(MG_BINS)) {
          lo <- MG_BINS[[mg_nm]][1]; hi <- MG_BINS[[mg_nm]][2]
          sub <- d24[!is.na(d24$first_mg_value) & d24$first_mg_value>=lo & d24$first_mg_value<hi,]
          r <- if(sum(sub$treated==1)>=10 && sum(sub$treated==0)>=10) did_dr(sub,ps_avail) else NULL
          if(!is.null(r)) {
            sig <- if(r$p_adj<0.05) " *" else ""
            cat(sprintf("    %s: n=%d/%d, DiD=%+.4f, P=%.4f%s\n",
                        mg_nm,r$n_trt,r$n_ctl,r$did_adj,r$p_adj,sig))
          }
        }
      }
      if ("surgery_type" %in% names(d24)) {
        cat("  Surgery type:\n")
        for (st in SURG_TYPES) {
          sub <- d24[!is.na(d24$surgery_type) & d24$surgery_type==st,]
          r <- if(sum(sub$treated==1)>=10 && sum(sub$treated==0)>=10) did_dr(sub,ps_avail) else NULL
          if(!is.null(r)) {
            sig <- if(r$p_adj<0.05) " *" else ""
            cat(sprintf("    %s: n=%d/%d, DiD=%+.4f, P=%.4f%s\n",
                        st,r$n_trt,r$n_ctl,r$did_adj,r$p_adj,sig))
          }
        }
      }

      # ── Secondary outcomes (matched logistic) ──────────────────────
      cat("\n  SECONDARY OUTCOMES (matched set, logistic)\n")
      cat("  outcome                n_trt  n_ctl  rate_trt rate_ctl  OR (95%CI)         P\n")
      cat("  --------------------- ------ ------ -------- -------- -----------------  ------\n")
      for (oc in SEC_OC) {
        if (!oc %in% names(md)) next
        md$y <- as.numeric(md[[oc]]); md_c <- md[!is.na(md$y),]
        nt <- sum(md_c$treated==1); nc <- sum(md_c$treated==0)
        if (nt<20||nc<20) next
        r1 <- mean(md_c$y[md_c$treated==1]); r0 <- mean(md_c$y[md_c$treated==0])
        fit <- tryCatch(glm(y~treated, data=md_c, family=quasibinomial()),
                        error=function(e) NULL)
        if (is.null(fit)||!"treated" %in% names(coef(fit))) next
        ct <- tryCatch(coeftest(fit, vcov.=vcovHC(fit,type="HC1")),
                       error=function(e) NULL)
        if (is.null(ct)) next
        tr <- which(rownames(ct)=="treated")
        if (length(tr)==0) next
        or <- exp(ct[tr,"Estimate"])
        or_lo <- exp(ct[tr,"Estimate"]-1.96*ct[tr,"Std. Error"])
        or_hi <- exp(ct[tr,"Estimate"]+1.96*ct[tr,"Std. Error"])
        p <- ct[tr,grep("^Pr",colnames(ct))]
        sig <- if(p<0.05) " *" else ""
        label <- ifelse(oc %in% names(SEC_LB), SEC_LB[oc], oc)
        cat(sprintf("  %-23s %5d  %5d   %.1f%%   %.1f%%  %.2f (%.2f-%.2f)  %.4f%s\n",
                    label,nt,nc,100*r1,100*r0,or,or_lo,or_hi,p,sig))
      }
    }
  }

  # ── AIPW on FULL sample (sensitivity) ──────────────────────────────
  cat(sprintf("\n  ── AIPW on full sample (sensitivity) ──\n"))
  d24_full <- build_dcr(combined, cr_all, PRIMARY_T)
  aipw_full <- run_aipw(d24_full, ps_avail)
  if (!is.null(aipw_full)) {
    cat(sprintf("  AIPW (full): DiD=%+.4f (SE=%.4f), P=%.4f, 95%%CI [%+.4f,%+.4f]\n",
                aipw_full$did, aipw_full$se, aipw_full$p,
                aipw_full$ci_lo, aipw_full$ci_hi))
    cat(sprintf("              n=%d/%d\n", aipw_full$n_trt, aipw_full$n_ctl))
  }

  cat(sprintf("\n%s\n%s [%s]: COMPLETE\n%s\n", SEP, db, mtag, SEP))
}

# ============================================================================
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  cat("Usage: Rscript 01_did_analysis.R <db> [model]\n"); quit(status=1)
}
db_arg <- toupper(args[1])
model_arg <- if(length(args)>=2) tolower(args[2]) else "primary"

cat("======================================================================\n")
cat(sprintf("01_did_analysis.R — PSM (1:1 + 1:4) + AIPW, model=%s\n", model_arg))
ps_vars <- select_model(model_arg)
cat("======================================================================\n")

run_analysis(db_arg, model_arg, ps_vars)
