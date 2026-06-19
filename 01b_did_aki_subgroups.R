#!/usr/bin/env Rscript
# ============================================================================
# 01b_did_aki_subgroups.R — AKI binary endpoint sweep (AIPW, ICU-time anchor)
#
# PS model (always included):
#   Demographics (3):   age, sex, BMI
#   Surgery type (3):   CABG, valve, combined
#   Comorbidities (8):  HF, HTN, DM, CKD, COPD, PVD, stroke, liver
#   Renal function (1): eGFR
#   Vitals (1):         heart rate
#   Labs (5):           K+, Ca, lactate, lactate_missing, Mg
#   ─────────────────── Total: 21 covariates (primary model)
#
# Method: AIPW risk difference (primary) + sIPTW_DR OR (secondary)
# Anchor: ICU admission time (Cr_pre: 0-6h, Cr_post: 18-30h = 24h +/-6h)
#
# Run: Rscript 01b_did_aki_subgroups.R eicu
#      Rscript 01b_did_aki_subgroups.R mimic
# ============================================================================

suppressPackageStartupMessages({ library(sandwich); library(lmtest) })
source(file.path(path.expand("~/mg_aki"), "did_covars.R"))

RESULTS    <- path.expand("~/mg_aki/results")
PS_COVARS  <- PS_PRIMARY   # 21 covariates from did_covars.R
TARGET_H   <- 24
WINDOW_H   <- 6

# ── Helpers ──────────────────────────────────────────────────────────────
median_impute <- function(d, vars) {
  for (v in vars) if (v %in% names(d) && any(is.na(d[[v]])))
    d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm=TRUE)
  d
}

# AIPW for binary outcome → risk difference
run_aipw_binary <- function(d, ps_vars, outcome="aki") {
  avail <- intersect(ps_vars, names(d))
  d <- d[complete.cases(d[, c(avail, outcome, "treated")]), ]
  nt <- sum(d$treated==1); nc <- sum(d$treated==0)
  if (nt<15 || nc<15) return(NULL)
  r1 <- mean(d[[outcome]][d$treated==1]); r0 <- mean(d[[outcome]][d$treated==0])
  ev1 <- sum(d[[outcome]][d$treated==1]); ev0 <- sum(d[[outcome]][d$treated==0])
  if (ev1+ev0 < 5) return(NULL)

  ps_fml <- as.formula(paste("treated ~", paste(avail, collapse="+")))
  ps_fit <- tryCatch(glm(ps_fml, data=d, family=binomial()), error=function(e) NULL)
  if (is.null(ps_fit)) return(NULL)
  e <- pmax(pmin(fitted(ps_fit), 0.99), 0.01)

  out_fml <- as.formula(paste(outcome, "~", paste(avail, collapse="+")))
  d1 <- d[d$treated==1,]; d0 <- d[d$treated==0,]
  m1 <- tryCatch(lm(out_fml, data=d1), error=function(e) NULL)
  m0 <- tryCatch(lm(out_fml, data=d0), error=function(e) NULL)
  if (is.null(m1)||is.null(m0)) return(NULL)

  mu1 <- predict(m1, newdata=d); mu0 <- predict(m0, newdata=d)
  Y <- d[[outcome]]; T_ <- d$treated
  phi <- (mu1-mu0) + T_*(Y-mu1)/e - (1-T_)*(Y-mu0)/(1-e)
  rd <- mean(phi); se <- sd(phi)/sqrt(nrow(d))
  p <- 2*pnorm(-abs(rd/se))

  list(n_trt=nt, n_ctl=nc, events_trt=ev1, events_ctl=ev0,
       rate_trt=r1, rate_ctl=r0, rd=rd, rd_se=se, rd_p=p,
       nnt=if(rd<0 && abs(rd)>0.001) round(1/abs(rd)) else NA)
}

# sIPTW_DR for binary → OR
run_siptw_binary <- function(d, ps_vars, outcome="aki") {
  avail <- intersect(ps_vars, names(d))
  d <- d[complete.cases(d[, c(avail, outcome, "treated")]), ]
  nt <- sum(d$treated==1); nc <- sum(d$treated==0)
  if (nt<15||nc<15) return(NULL)
  ev1 <- sum(d[[outcome]][d$treated==1]); ev0 <- sum(d[[outcome]][d$treated==0])
  if (ev1+ev0<5) return(NULL)

  ps_fml <- as.formula(paste("treated ~", paste(avail, collapse="+")))
  ps_fit <- tryCatch(glm(ps_fml, data=d, family=binomial()), error=function(e) NULL)
  if (is.null(ps_fit)) return(NULL)
  d$ps <- pmax(pmin(fitted(ps_fit),0.99),0.01)
  prev <- mean(d$treated)
  d$w <- ifelse(d$treated==1, prev/d$ps, (1-prev)/(1-d$ps))
  q01 <- quantile(d$w,0.01); q99 <- quantile(d$w,0.99)
  d$w <- pmax(pmin(d$w,q99),q01)

  d$y <- d[[outcome]]
  fit <- tryCatch(glm(y ~ treated, data=d, family=quasibinomial(), weights=w),
                  error=function(e) NULL)
  if (is.null(fit)||!"treated" %in% names(coef(fit))) return(NULL)
  vc <- tryCatch(vcovHC(fit,type="HC1"), error=function(e) vcov(fit))
  ct <- tryCatch(coeftest(fit,vcov.=vc), error=function(e) NULL)
  if (is.null(ct)) return(NULL)
  tr <- which(rownames(ct)=="treated")
  if (length(tr)==0) return(NULL)
  p_col <- grep("^Pr",colnames(ct))

  list(or=exp(ct[tr,"Estimate"]),
       or_lo=exp(ct[tr,"Estimate"]-1.96*ct[tr,"Std. Error"]),
       or_hi=exp(ct[tr,"Estimate"]+1.96*ct[tr,"Std. Error"]),
       or_p=ct[tr,p_col])
}

# ============================================================================
run_sweep <- function(db) {
  tag <- tolower(db)
  SEP <- paste(rep("=",70),collapse="")
  cat(sprintf("\n%s\n%s: AKI Subgroup Sweep (AIPW, ICU-time anchor, %d PS covariates)\n%s\n",
              SEP, db, length(PS_COVARS), SEP))

  # Load raw data
  trt <- read.csv(file.path(RESULTS,sprintf("did_treated_%s.csv",tag)),stringsAsFactors=F)
  ctl <- read.csv(file.path(RESULTS,sprintf("did_control_%s.csv",tag)),stringsAsFactors=F)
  cr_all <- read.csv(file.path(RESULTS,sprintf("did_cr_all_%s.csv",tag)),stringsAsFactors=F)

  id_col <- if("patientunitstayid" %in% names(trt)) "patientunitstayid" else "stay_id"
  trt$pid <- trt[[id_col]]; ctl$pid <- ctl[[id_col]]
  cr_id <- if("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
  cr_all$pid <- cr_all[[cr_id]]
  if(!"labresultoffset" %in% names(cr_all)) cr_all$labresultoffset <- cr_all$offset_min
  cr_all$offset_h <- cr_all$labresultoffset / 60

  # Stack covariates
  all_want <- unique(c("pid","treated",PS_COVARS,"surgery_type","first_mg_value"))
  sh <- intersect(all_want, intersect(names(trt),names(ctl)))
  combined <- rbind(trt[,sh], ctl[,sh])
  combined <- median_impute(combined, PS_COVARS)

  # Build ICU-time-anchored Cr
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

  df <- merge(combined, m[,c("pid","cr_pre","cr_post","delta_cr")], by="pid")
  cat(sprintf("  Built dataset: %d patients (%d treated, %d control)\n",
              nrow(df), sum(df$treated==1), sum(df$treated==0)))

  # Define KDIGO AKI stages
  df$aki_kdigo1 <- as.integer(df$delta_cr >= 0.3)
  df$aki_kdigo2 <- as.integer(df$cr_pre>0 & df$cr_post >= 2.0*df$cr_pre)
  df$aki_kdigo3 <- as.integer((df$cr_pre>0 & df$cr_post >= 3.0*df$cr_pre) | df$cr_post >= 4.0)

  cat("  AKI rates (overall):\n")
  for (oc in c("aki_kdigo1","aki_kdigo2","aki_kdigo3")) {
    r1 <- 100*mean(df[[oc]][df$treated==1],na.rm=T)
    r0 <- 100*mean(df[[oc]][df$treated==0],na.rm=T)
    cat(sprintf("    %s: treated %.1f%%, control %.1f%%\n", oc, r1, r0))
  }

  # Define subgroups
  subgroups <- list(
    list(name="Overall", filter=NULL),
    list(name="Age < 65", filter=expression(age < 65)),
    list(name="Age >= 65", filter=expression(age >= 65)),
    list(name="Age >= 75", filter=expression(age >= 75)),
    list(name="Female", filter=expression(is_female == 1)),
    list(name="Male", filter=expression(is_female == 0)),
    list(name="eGFR < 45", filter=expression(egfr < 45)),
    list(name="eGFR 45-60", filter=expression(egfr >= 45 & egfr < 60)),
    list(name="eGFR 60-90", filter=expression(egfr >= 60 & egfr < 90)),
    list(name="eGFR >= 90", filter=expression(egfr >= 90)),
    list(name="Cr_pre <= 1.0", filter=expression(cr_pre <= 1.0)),
    list(name="Cr_pre > 1.0", filter=expression(cr_pre > 1.0)),
    list(name="Cr_pre > 1.5", filter=expression(cr_pre > 1.5)),
    list(name="Diabetes", filter=expression(diabetes == 1)),
    list(name="No diabetes", filter=expression(diabetes == 0)),
    list(name="CKD", filter=expression(ckd == 1)),
    list(name="No CKD", filter=expression(ckd == 0)),
    list(name="Heart failure", filter=expression(heart_failure == 1)),
    list(name="No HF", filter=expression(heart_failure == 0)),
    list(name="CABG", filter=expression(surg_cabg == 1)),
    list(name="Valve", filter=expression(surg_valve == 1)),
    list(name="Mg < 1.8", filter=expression(first_mg_value < 1.8)),
    list(name="Mg 1.8-2.0", filter=expression(first_mg_value >= 1.8 & first_mg_value < 2.0)),
    list(name="Mg >= 2.0", filter=expression(first_mg_value >= 2.0)),
    list(name="BMI < 25", filter=expression(bmi < 25)),
    list(name="BMI 25-30", filter=expression(bmi >= 25 & bmi < 30)),
    list(name="BMI >= 30", filter=expression(bmi >= 30))
  )

  # Sweep
  aki_eps <- c("aki_kdigo1","aki_kdigo2","aki_kdigo3")
  aki_lab <- c(aki_kdigo1="KDIGO>=1", aki_kdigo2="KDIGO>=2", aki_kdigo3="KDIGO>=3")

  all_rows <- list(); ridx <- 0
  for (sg in subgroups) {
    sub <- if(is.null(sg$filter)) df else tryCatch(df[eval(sg$filter, df),], error=function(e) df[0,])
    sub <- sub[complete.cases(sub[,c("treated","delta_cr","cr_pre","cr_post")]),]

    for (ep in aki_eps) {
      ridx <- ridx+1
      ra <- tryCatch(run_aipw_binary(sub, PS_COVARS, ep), error=function(e) NULL)
      rs <- tryCatch(run_siptw_binary(sub, PS_COVARS, ep), error=function(e) NULL)

      all_rows[[ridx]] <- data.frame(
        subgroup=sg$name, endpoint=ep,
        n_trt=if(!is.null(ra)) ra$n_trt else sum(sub$treated==1),
        n_ctl=if(!is.null(ra)) ra$n_ctl else sum(sub$treated==0),
        events_trt=if(!is.null(ra)) ra$events_trt else NA,
        events_ctl=if(!is.null(ra)) ra$events_ctl else NA,
        rate_trt=if(!is.null(ra)) round(100*ra$rate_trt,1) else NA,
        rate_ctl=if(!is.null(ra)) round(100*ra$rate_ctl,1) else NA,
        aipw_rd=if(!is.null(ra)) round(100*ra$rd,1) else NA,
        aipw_p=if(!is.null(ra)) round(ra$rd_p,4) else NA,
        or=if(!is.null(rs)) round(rs$or,2) else NA,
        or_lo=if(!is.null(rs)) round(rs$or_lo,2) else NA,
        or_hi=if(!is.null(rs)) round(rs$or_hi,2) else NA,
        or_p=if(!is.null(rs)) round(rs$or_p,4) else NA,
        nnt=if(!is.null(ra)) ra$nnt else NA,
        stringsAsFactors=F)
    }
  }

  results <- do.call(rbind, all_rows)
  write.csv(results, file.path(RESULTS,sprintf("did_aki_subgroups_%s.csv",tag)), row.names=F)

  # Print KDIGO >= 1
  cat(sprintf("\n%s\nKDIGO >= Stage 1 (dCr >= 0.3 mg/dL)\n%s\n", SEP, SEP))
  cat("  subgroup              n_trt  n_ctl  AKI_trt AKI_ctl  AIPW_RD   AIPW_P  OR (95%CI)          OR_P    NNT\n")
  cat("  ───────────────────  ─────  ─────  ─────── ───────  ───────  ───────  ──────────────────  ──────  ────\n")

  k1 <- results[results$endpoint=="aki_kdigo1",]
  for (i in seq_len(nrow(k1))) {
    r <- k1[i,]
    if (is.na(r$aipw_rd)) {
      cat(sprintf("  %-21s  %5d  %5d      —       —       —        —          —               —     —\n",
                  r$subgroup, r$n_trt, r$n_ctl)); next
    }
    sig_a <- if(!is.na(r$aipw_p) && r$aipw_p<0.05) "*" else " "
    sig_o <- if(!is.na(r$or_p) && r$or_p<0.05) "*" else " "
    nnt_s <- if(!is.na(r$nnt)) sprintf("%4d",r$nnt) else "   —"
    or_s <- if(!is.na(r$or)) sprintf("%.2f (%.2f-%.2f)", r$or, r$or_lo, r$or_hi) else "       —        "
    cat(sprintf("  %-21s  %5d  %5d   %5.1f%%  %5.1f%%  %+5.1f%%  %.4f%s %s  %.4f%s %s\n",
                r$subgroup, r$n_trt, r$n_ctl,
                r$rate_trt, r$rate_ctl, r$aipw_rd, r$aipw_p, sig_a,
                or_s, r$or_p, sig_o, nnt_s))
  }

  # Notable KDIGO >= 2/3
  for (ep in c("aki_kdigo2","aki_kdigo3")) {
    sub <- results[results$endpoint==ep & !is.na(results$aipw_p) & results$aipw_p<0.1,]
    if (nrow(sub)>0) {
      cat(sprintf("\n  %s — notable (AIPW P<0.1):\n", aki_lab[ep]))
      for (i in seq_len(nrow(sub))) {
        r <- sub[i,]
        sig <- if(r$aipw_p<0.05) " *" else ""
        cat(sprintf("    %s: %.1f%% vs %.1f%%, RD=%+.1f%%, P=%.4f%s\n",
                    r$subgroup, r$rate_trt, r$rate_ctl, r$aipw_rd, r$aipw_p, sig))
      }
    }
  }

  cat(sprintf("\n  Saved: did_aki_subgroups_%s.csv (%d rows)\n", tag, nrow(results)))
}

# ============================================================================
cat("======================================================================\n")
cat("01b_did_aki_subgroups.R — AIPW + sIPTW_DR, ICU-time anchor\n")
cat(sprintf("  PS model (21 covariates):\n"))
cat(sprintf("    Demographics:    age, sex, BMI\n"))
cat(sprintf("    Surgery:         CABG, valve, combined\n"))
cat(sprintf("    Comorbidities:   HF, HTN, DM, CKD, COPD, PVD, stroke, liver\n"))
cat(sprintf("    Renal function:  eGFR\n"))
cat(sprintf("    Vitals:          heart rate\n"))
cat(sprintf("    Labs:            K+, Ca, lactate, lactate_missing, Mg\n"))
cat("  Endpoints: KDIGO >= 1, 2, 3\n")
cat("======================================================================\n")

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) { cat("Usage: Rscript 01b_did_aki_subgroups.R eicu|mimic\n"); quit(status=1) }
for (a in args) run_sweep(toupper(a))
