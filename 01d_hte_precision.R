#!/usr/bin/env Rscript
# ============================================================================
# 01d_hte_precision.R — Heterogeneous Treatment Effect / Precision Medicine
#
# Q: Which patients benefit most from IV Mg? Where is it null? Harmful?
#
# Analyses:
#   A. Risk-stratified AIPW: build AKI risk score from controls → quintiles
#   B. Progressive enrichment: exclude low-risk → watch signal amplify
#   C. Sliding-window AIPW: treatment effect as f(eGFR), f(Cr_pre), f(age)
#   D. Multi-dimensional enrichment: combined risk factors
#
# All use AIPW estimator with ICU-time anchor.
#
# Run: Rscript 01d_hte_precision.R eicu
#      Rscript 01d_hte_precision.R mimic
# ============================================================================

source(file.path(path.expand("~/mg_aki"), "did_covars.R"))

RESULTS   <- path.expand("~/mg_aki/results")
PS_COVARS <- PS_PRIMARY
TARGET_H  <- 24; WINDOW_H <- 6

median_impute <- function(d, vars) {
  for (v in vars) if (v %in% names(d) && any(is.na(d[[v]])))
    d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm=TRUE)
  d
}

# ── AIPW (continuous ΔCr) ───────────────────────────────────────────────
aipw <- function(d, ps_vars, outcome="delta_cr") {
  avail <- intersect(ps_vars, names(d))
  d <- d[complete.cases(d[, c(avail, outcome, "treated")]), ]
  nt <- sum(d$treated==1); nc <- sum(d$treated==0)
  if (nt<15||nc<15) return(NULL)

  ps_fml <- as.formula(paste("treated ~", paste(avail, collapse="+")))
  ps_fit <- tryCatch(glm(ps_fml, data=d, family=binomial()), error=function(e) NULL)
  if (is.null(ps_fit)) return(NULL)
  e <- pmax(pmin(fitted(ps_fit), 0.99), 0.01)

  out_fml <- as.formula(paste(outcome, "~", paste(avail, collapse="+")))
  m1 <- tryCatch(lm(out_fml, data=d[d$treated==1,]), error=function(e) NULL)
  m0 <- tryCatch(lm(out_fml, data=d[d$treated==0,]), error=function(e) NULL)
  if (is.null(m1)||is.null(m0)) return(NULL)

  mu1 <- predict(m1, newdata=d); mu0 <- predict(m0, newdata=d)
  Y <- d[[outcome]]; T_ <- d$treated
  phi <- (mu1-mu0) + T_*(Y-mu1)/e - (1-T_)*(Y-mu0)/(1-e)
  tau <- mean(phi); se <- sd(phi)/sqrt(nrow(d))

  list(tau=tau, se=se, p=2*pnorm(-abs(tau/se)),
       ci_lo=tau-1.96*se, ci_hi=tau+1.96*se,
       n_trt=nt, n_ctl=nc, n=nrow(d))
}

# ── AIPW (binary AKI) ──────────────────────────────────────────────────
aipw_bin <- function(d, ps_vars, outcome="aki") {
  avail <- intersect(ps_vars, names(d))
  d <- d[complete.cases(d[, c(avail, outcome, "treated")]), ]
  nt <- sum(d$treated==1); nc <- sum(d$treated==0)
  if (nt<15||nc<15) return(NULL)
  r1 <- mean(d[[outcome]][d$treated==1]); r0 <- mean(d[[outcome]][d$treated==0])

  ps_fml <- as.formula(paste("treated ~", paste(avail, collapse="+")))
  ps_fit <- tryCatch(glm(ps_fml, data=d, family=binomial()), error=function(e) NULL)
  if (is.null(ps_fit)) return(NULL)
  e <- pmax(pmin(fitted(ps_fit), 0.99), 0.01)

  out_fml <- as.formula(paste(outcome, "~", paste(avail, collapse="+")))
  m1 <- tryCatch(lm(out_fml, data=d[d$treated==1,]), error=function(e) NULL)
  m0 <- tryCatch(lm(out_fml, data=d[d$treated==0,]), error=function(e) NULL)
  if (is.null(m1)||is.null(m0)) return(NULL)

  mu1 <- predict(m1, newdata=d); mu0 <- predict(m0, newdata=d)
  Y <- d[[outcome]]; T_ <- d$treated
  phi <- (mu1-mu0) + T_*(Y-mu1)/e - (1-T_)*(Y-mu0)/(1-e)
  rd <- mean(phi); se <- sd(phi)/sqrt(nrow(d))

  list(rd=rd, se=se, p=2*pnorm(-abs(rd/se)),
       ci_lo=rd-1.96*se, ci_hi=rd+1.96*se,
       n_trt=nt, n_ctl=nc, rate_trt=r1, rate_ctl=r0)
}


# ============================================================================
run_hte <- function(db) {
  tag <- tolower(db)
  SEP <- paste(rep("=",70), collapse="")
  cat(sprintf("\n%s\n%s: HTE / Precision Medicine Analysis\n%s\n", SEP, db, SEP))

  # ── Load & build ────────────────────────────────────────────────────
  trt <- read.csv(file.path(RESULTS,sprintf("did_treated_%s.csv",tag)),stringsAsFactors=F)
  ctl <- read.csv(file.path(RESULTS,sprintf("did_control_%s.csv",tag)),stringsAsFactors=F)
  cr_all <- read.csv(file.path(RESULTS,sprintf("did_cr_all_%s.csv",tag)),stringsAsFactors=F)

  id_col <- if("patientunitstayid" %in% names(trt)) "patientunitstayid" else "stay_id"
  trt$pid <- trt[[id_col]]; ctl$pid <- ctl[[id_col]]
  cr_id <- if("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
  cr_all$pid <- cr_all[[cr_id]]
  if(!"labresultoffset" %in% names(cr_all)) cr_all$labresultoffset <- cr_all$offset_min
  cr_all$offset_h <- cr_all$labresultoffset / 60

  all_want <- unique(c("pid","treated",PS_COVARS,"surgery_type","first_mg_value"))
  sh <- intersect(all_want, intersect(names(trt),names(ctl)))
  full <- rbind(trt[,sh], ctl[,sh])
  full <- median_impute(full, PS_COVARS)

  # Build 24h ΔCr
  pre <- cr_all[cr_all$offset_h>=0 & cr_all$offset_h<=6,]
  pre <- pre[order(pre$pid, pre$offset_h),]; pre <- pre[!duplicated(pre$pid),]
  post <- cr_all[cr_all$offset_h>=(TARGET_H-WINDOW_H) & cr_all$offset_h<=(TARGET_H+WINDOW_H),]
  post$dist <- abs(post$offset_h - TARGET_H)
  post <- post[order(post$pid, post$dist),]; post <- post[!duplicated(post$pid),]
  m <- merge(pre[,c("pid","labresult","offset_h")],
             post[,c("pid","labresult","offset_h")],
             by="pid", suffixes=c("_pre","_post"))
  m <- m[m$offset_h_post > m$offset_h_pre,]
  m$cr_pre <- m$labresult_pre; m$cr_post <- m$labresult_post
  m$delta_cr <- m$cr_post - m$cr_pre
  m$aki <- as.integer(m$delta_cr >= 0.3)

  df <- merge(full, m[,c("pid","cr_pre","cr_post","delta_cr","aki")], by="pid")
  cat(sprintf("  Dataset: %d treated + %d control = %d total\n",
              sum(df$treated==1), sum(df$treated==0), nrow(df)))

  # ══════════════════════════════════════════════════════════════════════
  # A. RISK-STRATIFIED AIPW
  # ══════════════════════════════════════════════════════════════════════
  cat(sprintf("\n%s\nA. Risk-stratified AIPW (AKI risk score from controls)\n%s\n",
              paste(rep("-",60),collapse=""), paste(rep("-",60),collapse="")))

  # Build risk score from controls
  controls <- df[df$treated==0,]
  risk_vars <- intersect(PS_COVARS, names(controls))
  risk_fml <- as.formula(paste("aki ~", paste(risk_vars, collapse="+")))
  risk_mod <- glm(risk_fml, data=controls, family=binomial())

  # Predict for all patients
  df$risk_score <- predict(risk_mod, newdata=df, type="response")
  df$risk_q <- cut(df$risk_score,
                    breaks=quantile(df$risk_score, probs=seq(0,1,0.2), na.rm=T),
                    labels=c("Q1 (lowest)","Q2","Q3","Q4","Q5 (highest)"),
                    include.lowest=TRUE)

  cat("  AKI risk score (predicted from controls):\n")
  cat(sprintf("    Overall AKI rate: treated %.1f%%, control %.1f%%\n",
              100*mean(df$aki[df$treated==1]), 100*mean(df$aki[df$treated==0])))

  cat("\n  Risk quintile    n_trt  n_ctl  AKI_trt  AKI_ctl  AIPW_ΔCr     P       AIPW_RD_AKI    P\n")
  cat("  ──────────────  ─────  ─────  ──────  ───────  ──────────  ──────  ──────────────  ──────\n")

  risk_rows <- list()
  for (q in levels(df$risk_q)) {
    sub <- df[df$risk_q == q, ]
    r_cr <- aipw(sub, PS_COVARS, "delta_cr")
    r_aki <- aipw_bin(sub, PS_COVARS, "aki")

    if (!is.null(r_cr) && !is.null(r_aki)) {
      sig_cr <- if(r_cr$p<0.05) "*" else " "
      sig_aki <- if(r_aki$p<0.05) "*" else " "
      cat(sprintf("  %-16s %5d  %5d   %5.1f%%  %5.1f%%   %+.4f %s  %.4f  %+5.1f%% [%+.1f,%+.1f] %s %.4f\n",
                  q, r_cr$n_trt, r_cr$n_ctl,
                  100*r_aki$rate_trt, 100*r_aki$rate_ctl,
                  r_cr$tau, sig_cr, r_cr$p,
                  100*r_aki$rd, 100*r_aki$ci_lo, 100*r_aki$ci_hi, sig_aki, r_aki$p))

      risk_rows[[q]] <- data.frame(
        quintile=q, n_trt=r_cr$n_trt, n_ctl=r_cr$n_ctl,
        aki_trt=round(100*r_aki$rate_trt,1), aki_ctl=round(100*r_aki$rate_ctl,1),
        aipw_dcr=round(r_cr$tau,4), dcr_p=round(r_cr$p,4),
        aipw_rd=round(100*r_aki$rd,2), rd_p=round(r_aki$p,4),
        stringsAsFactors=F)
    }
  }

  # Risk score interaction test
  ps_avail <- intersect(PS_COVARS, names(df))
  int_fml <- as.formula(paste("delta_cr ~ treated * risk_score +",
                               paste(ps_avail, collapse="+")))
  int_fit <- lm(int_fml, data=df)
  int_ct <- coef(summary(int_fit))
  ir <- grep("treated:risk_score", rownames(int_ct))
  if (length(ir)>0)
    cat(sprintf("\n  Interaction (treated × risk_score): β=%+.3f, P=%.4f\n",
                int_ct[ir,"Estimate"], int_ct[ir,"Pr(>|t|)"]))

  # ══════════════════════════════════════════════════════════════════════
  # B. PROGRESSIVE ENRICHMENT
  # ══════════════════════════════════════════════════════════════════════
  cat(sprintf("\n%s\nB. Progressive enrichment (exclude low-risk patients)\n%s\n",
              paste(rep("-",60),collapse=""), paste(rep("-",60),collapse="")))

  enrichment_steps <- list(
    list(name="All patients", filter=expression(TRUE)),
    list(name="Excl eGFR>=90", filter=expression(egfr < 90)),
    list(name="Excl eGFR>=90 + Cr<=0.7", filter=expression(egfr < 90 & cr_pre > 0.7)),
    list(name="Excl eGFR>=90 + Cr<=0.7 + no-DM-CKD",
         filter=expression(egfr < 90 & cr_pre > 0.7 & (diabetes==1 | ckd==1))),
    list(name="eGFR<60 only", filter=expression(egfr < 60)),
    list(name="eGFR<60 + DM or CKD", filter=expression(egfr < 60 & (diabetes==1 | ckd==1))),
    list(name="eGFR<45 only", filter=expression(egfr < 45)),
    list(name="Cr_pre > 1.2", filter=expression(cr_pre > 1.2)),
    list(name="Cr_pre > 1.2 + eGFR<60", filter=expression(cr_pre > 1.2 & egfr < 60)),
    list(name="Risk Q4+Q5 only", filter=expression(risk_score >= quantile(risk_score, 0.6)))
  )

  cat("  Population                         n_trt  n_ctl  AIPW_ΔCr       P       AKI_RD(%)     P\n")
  cat("  ─────────────────────────────────  ─────  ─────  ──────────  ──────  ───────────── ──────\n")

  enrich_rows <- list()
  for (i in seq_along(enrichment_steps)) {
    es <- enrichment_steps[[i]]
    sub <- tryCatch(df[eval(es$filter, df), ], error=function(e) df[0,])
    if (nrow(sub)<40) {
      cat(sprintf("  %-35s  -- too few (%d)\n", es$name, nrow(sub))); next
    }
    r_cr <- aipw(sub, PS_COVARS, "delta_cr")
    r_aki <- aipw_bin(sub, PS_COVARS, "aki")

    if (!is.null(r_cr)) {
      sig_cr <- if(r_cr$p<0.05) "*" else " "
      sig_aki <- if(!is.null(r_aki) && r_aki$p<0.05) "*" else " "
      aki_str <- if(!is.null(r_aki)) sprintf("%+5.1f%%  %s %.4f", 100*r_aki$rd, sig_aki, r_aki$p) else "  —"
      cat(sprintf("  %-35s %5d  %5d  %+.4f %s  %.4f  %s\n",
                  es$name, r_cr$n_trt, r_cr$n_ctl,
                  r_cr$tau, sig_cr, r_cr$p, aki_str))
      enrich_rows[[i]] <- data.frame(step=es$name, n_trt=r_cr$n_trt, n_ctl=r_cr$n_ctl,
        aipw_dcr=round(r_cr$tau,4), dcr_p=round(r_cr$p,4),
        aki_rd=if(!is.null(r_aki)) round(100*r_aki$rd,2) else NA,
        aki_p=if(!is.null(r_aki)) round(r_aki$p,4) else NA, stringsAsFactors=F)
    }
  }

  # ══════════════════════════════════════════════════════════════════════
  # C. SLIDING WINDOW: treatment effect as f(covariate)
  # ══════════════════════════════════════════════════════════════════════
  cat(sprintf("\n%s\nC. Sliding window AIPW: treatment effect vs continuous covariates\n%s\n",
              paste(rep("-",60),collapse=""), paste(rep("-",60),collapse="")))

  slide_vars <- list(
    list(var="egfr", label="eGFR", breaks=seq(20, 100, 10)),
    list(var="cr_pre", label="Baseline Cr", breaks=seq(0.5, 2.5, 0.2)),
    list(var="age", label="Age", breaks=seq(40, 85, 5)),
    list(var="first_mg_value", label="Serum Mg", breaks=seq(1.2, 2.8, 0.2))
  )

  all_slide <- list()
  for (sv in slide_vars) {
    if (!sv$var %in% names(df)) next
    cat(sprintf("\n  %s:\n", sv$label))
    cat("    midpoint   n_trt  n_ctl  AIPW_ΔCr      P\n")
    cat("    ────────  ─────  ─────  ──────────  ──────\n")

    brks <- sv$breaks
    width <- brks[2] - brks[1]
    for (mid in brks) {
      lo <- mid - width; hi <- mid + width
      sub <- df[!is.na(df[[sv$var]]) & df[[sv$var]] >= lo & df[[sv$var]] < hi, ]
      if (nrow(sub)<60) next
      r <- aipw(sub, PS_COVARS, "delta_cr")
      if (!is.null(r)) {
        sig <- if(r$p<0.05) "*" else " "
        cat(sprintf("    %7.1f   %5d  %5d  %+.4f %s  %.4f\n",
                    mid, r$n_trt, r$n_ctl, r$tau, sig, r$p))
        all_slide[[length(all_slide)+1]] <- data.frame(
          var=sv$var, label=sv$label, midpoint=mid,
          n_trt=r$n_trt, n_ctl=r$n_ctl,
          tau=round(r$tau,5), se=round(r$se,5), p=round(r$p,4),
          ci_lo=round(r$ci_lo,5), ci_hi=round(r$ci_hi,5), stringsAsFactors=F)
      }
    }
  }

  # ══════════════════════════════════════════════════════════════════════
  # D. MULTI-DIMENSIONAL: Combined high-risk phenotype
  # ══════════════════════════════════════════════════════════════════════
  cat(sprintf("\n%s\nD. Multi-dimensional enrichment: combined phenotypes\n%s\n",
              paste(rep("-",60),collapse=""), paste(rep("-",60),collapse="")))

  phenotypes <- list(
    list(name="High-risk renal (eGFR<60 OR Cr>1.2)",
         filter=expression(egfr<60 | cr_pre>1.2)),
    list(name="Metabolic syndrome (DM + BMI>=30)",
         filter=expression(diabetes==1 & bmi>=30)),
    list(name="Cardiac burden (HF + CABG/combined)",
         filter=expression(heart_failure==1 & (surg_cabg==1|surg_combined==1))),
    list(name="Low Mg + high-risk renal",
         filter=expression(first_mg_value<1.8 & (egfr<60|cr_pre>1.2))),
    list(name="Normal Mg + high-risk renal (prophylactic benefit?)",
         filter=expression(first_mg_value>=1.8 & first_mg_value<2.3 & (egfr<60|cr_pre>1.2))),
    list(name="Elderly + CKD (age>=65 + eGFR<60)",
         filter=expression(age>=65 & egfr<60)),
    list(name="Young + preserved renal (age<55 + eGFR>=90)",
         filter=expression(age<55 & egfr>=90)),
    list(name="CABG + eGFR<90",
         filter=expression(surg_cabg==1 & egfr<90))
  )

  cat("  Phenotype                                    n_trt  n_ctl  AIPW_ΔCr       P       AKI_RD      P\n")
  cat("  ───────────────────────────────────────────  ─────  ─────  ──────────  ──────  ────────── ──────\n")

  pheno_rows <- list()
  for (ph in phenotypes) {
    sub <- tryCatch(df[eval(ph$filter, df), ], error=function(e) df[0,])
    if (nrow(sub)<40) {
      cat(sprintf("  %-45s  -- too few (%d)\n", ph$name, nrow(sub))); next
    }
    r_cr <- aipw(sub, PS_COVARS, "delta_cr")
    r_aki <- aipw_bin(sub, PS_COVARS, "aki")
    if (!is.null(r_cr)) {
      sig_cr <- if(r_cr$p<0.05) "*" else " "
      sig_aki <- if(!is.null(r_aki) && r_aki$p<0.05) "*" else " "
      aki_str <- if(!is.null(r_aki)) sprintf("%+5.1f%% %s %.4f", 100*r_aki$rd, sig_aki, r_aki$p) else "  —"
      cat(sprintf("  %-45s %5d  %5d  %+.4f %s  %.4f  %s\n",
                  ph$name, r_cr$n_trt, r_cr$n_ctl, r_cr$tau, sig_cr, r_cr$p, aki_str))
      pheno_rows[[length(pheno_rows)+1]] <- data.frame(
        phenotype=ph$name, n_trt=r_cr$n_trt, n_ctl=r_cr$n_ctl,
        aipw_dcr=round(r_cr$tau,4), dcr_p=round(r_cr$p,4),
        aki_rd=if(!is.null(r_aki)) round(100*r_aki$rd,2) else NA,
        aki_p=if(!is.null(r_aki)) round(r_aki$p,4) else NA, stringsAsFactors=F)
    }
  }

  # ── Save all results ────────────────────────────────────────────────
  if (length(risk_rows)>0)
    write.csv(do.call(rbind, risk_rows),
              file.path(RESULTS, sprintf("hte_risk_strata_%s.csv",tag)), row.names=F)
  if (length(enrich_rows)>0)
    write.csv(do.call(rbind, enrich_rows),
              file.path(RESULTS, sprintf("hte_enrichment_%s.csv",tag)), row.names=F)
  if (length(all_slide)>0)
    write.csv(do.call(rbind, all_slide),
              file.path(RESULTS, sprintf("hte_sliding_%s.csv",tag)), row.names=F)
  if (length(pheno_rows)>0)
    write.csv(do.call(rbind, pheno_rows),
              file.path(RESULTS, sprintf("hte_phenotypes_%s.csv",tag)), row.names=F)

  cat(sprintf("\n  Saved: hte_risk_strata_%s.csv, hte_enrichment_%s.csv\n", tag, tag))
  cat(sprintf("         hte_sliding_%s.csv, hte_phenotypes_%s.csv\n", tag, tag))
}

# ============================================================================
cat("======================================================================\n")
cat("01d_hte_precision.R — Precision Medicine / HTE Analysis\n")
cat("  A: Risk-stratified AIPW (quintiles)\n")
cat("  B: Progressive enrichment\n")
cat("  C: Sliding-window AIPW (continuous modifiers)\n")
cat("  D: Multi-dimensional phenotypes\n")
cat("======================================================================\n")

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) { cat("Usage: Rscript 01d_hte_precision.R eicu|mimic\n"); quit(status=1) }
for (a in args) run_hte(toupper(a))
