#!/usr/bin/env Rscript
# ============================================================================
# 01b_did_aki_subgroups.R — Downstream subgroup analysis (sIPTW_DR OR)
#
# Approach: compute IPTW weights ONCE on full sample, then stratify.
#   AKI outcomes: 24h ICU-time-anchor subset (Cr_pre + Cr_post)
#   Other outcomes: full sample (mortality, POAF, encephalopathy, vent arr)
# Reports OR + interaction P per subgroup × outcome. Ref groups explicit.
#
# Run: Rscript 01b_did_aki_subgroups.R eicu
#      Rscript 01b_did_aki_subgroups.R mimic
# ============================================================================

suppressPackageStartupMessages({ library(sandwich); library(lmtest) })
source(file.path(path.expand("~/mg_aki"), "did_covars.R"))

RESULTS   <- path.expand("~/mg_aki/results")
PS_COVARS <- PS_PRIMARY
TARGET_H  <- 24; WINDOW_H <- 6

# ── Helpers ──────────────────────────────────────────────────────────────
median_impute <- function(d, vars) {
  for (v in vars) if (v %in% names(d) && any(is.na(d[[v]])))
    d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm=TRUE)
  d
}

# Weighted OR from pre-computed IPTW
weighted_or <- function(d, outcome, wt_col="w") {
  d$y <- as.numeric(d[[outcome]]); d <- d[!is.na(d$y),]
  nt <- sum(d$treated==1); nc <- sum(d$treated==0)
  if (nt<15||nc<15) return(NULL)
  ev1 <- sum(d$y[d$treated==1]); ev0 <- sum(d$y[d$treated==0])
  if (ev1+ev0<3) return(NULL)
  r1 <- mean(d$y[d$treated==1]); r0 <- mean(d$y[d$treated==0])

  fit <- tryCatch(glm(y ~ treated, data=d, family=quasibinomial(), weights=d[[wt_col]]),
                  error=function(e) NULL)
  if (is.null(fit)||!"treated" %in% names(coef(fit))) return(NULL)
  vc <- tryCatch(vcovHC(fit, type="HC1"), error=function(e) vcov(fit))
  ct <- tryCatch(coeftest(fit, vcov.=vc), error=function(e) NULL)
  if (is.null(ct)) return(NULL)
  tr <- which(rownames(ct)=="treated")
  if (length(tr)==0) return(NULL)
  p_col <- grep("^Pr", colnames(ct))

  list(n_trt=nt, n_ctl=nc, ev_trt=ev1, ev_ctl=ev0,
       rate_trt=r1, rate_ctl=r0,
       or=exp(ct[tr,"Estimate"]),
       or_lo=exp(ct[tr,"Estimate"]-1.96*ct[tr,"Std. Error"]),
       or_hi=exp(ct[tr,"Estimate"]+1.96*ct[tr,"Std. Error"]),
       p=ct[tr, p_col])
}

# Interaction P: treated × continuous/binary subgroup variable
interaction_p <- function(d, outcome, sg_var, wt_col="w") {
  d$y <- as.numeric(d[[outcome]]); d <- d[!is.na(d$y) & !is.na(d[[sg_var]]),]
  if (nrow(d)<40) return(NA)
  fml <- as.formula(paste("y ~ treated *", sg_var))
  fit <- tryCatch(glm(fml, data=d, family=quasibinomial(), weights=d[[wt_col]]),
                  error=function(e) NULL)
  if (is.null(fit)) return(NA)
  vc <- tryCatch(vcovHC(fit, type="HC1"), error=function(e) vcov(fit))
  ct <- tryCatch(coeftest(fit, vcov.=vc), error=function(e) NULL)
  if (is.null(ct)) return(NA)
  ir <- grep(paste0("treated:",sg_var), rownames(ct))
  if (length(ir)==0) ir <- grep(paste0(sg_var,":treated"), rownames(ct))
  if (length(ir)==0) return(NA)
  p_col <- grep("^Pr", colnames(ct))
  ct[ir[1], p_col]
}


# ============================================================================
run_sweep <- function(db) {
  tag <- tolower(db)
  SEP <- paste(rep("=",70), collapse="")
  cat(sprintf("\n%s\n%s: Subgroup Analysis (sIPTW_DR downstream, %d PS covariates)\n%s\n",
              SEP, db, length(PS_COVARS), SEP))

  # ── Load & stack ─────────────────────────────────────────────────────
  trt <- read.csv(file.path(RESULTS,sprintf("did_treated_%s.csv",tag)),stringsAsFactors=F)
  ctl <- read.csv(file.path(RESULTS,sprintf("did_control_%s.csv",tag)),stringsAsFactors=F)
  cr_all <- read.csv(file.path(RESULTS,sprintf("did_cr_all_%s.csv",tag)),stringsAsFactors=F)

  id_col <- if("patientunitstayid" %in% names(trt)) "patientunitstayid" else "stay_id"
  trt$pid <- trt[[id_col]]; ctl$pid <- ctl[[id_col]]
  cr_id <- if("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
  cr_all$pid <- cr_all[[cr_id]]
  if(!"labresultoffset" %in% names(cr_all)) cr_all$labresultoffset <- cr_all$offset_min
  cr_all$offset_h <- cr_all$labresultoffset / 60

  all_want <- unique(c("pid","treated",PS_COVARS,"surgery_type","first_mg_value",
                        "hosp_mortality","poaf","encephalopathy","vent_arrhythmia","prior_af"))
  sh <- intersect(all_want, intersect(names(trt),names(ctl)))
  full <- rbind(trt[,sh], ctl[,sh])
  full <- median_impute(full, PS_COVARS)
  ps_avail <- intersect(PS_COVARS, names(full))

  # ── Compute IPTW weights ONCE ───────────────────────────────────────
  ps_fml <- as.formula(paste("treated ~", paste(ps_avail, collapse="+")))
  ps_fit <- glm(ps_fml, data=full, family=binomial())
  full$ps <- pmax(pmin(fitted(ps_fit), 0.99), 0.01)
  prev <- mean(full$treated)
  full$w <- ifelse(full$treated==1, prev/full$ps, (1-prev)/(1-full$ps))
  q01 <- quantile(full$w, 0.01); q99 <- quantile(full$w, 0.99)
  full$w <- pmax(pmin(full$w, q99), q01)

  cat(sprintf("  Full sample: %d treated + %d control\n",
              sum(full$treated==1), sum(full$treated==0)))

  # ── Build AKI dataset (24h ICU-time anchor) ─────────────────────────
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

  aki_df <- merge(full, m[,c("pid","cr_pre","cr_post","delta_cr")], by="pid")
  aki_df$aki_kdigo1 <- as.integer(aki_df$delta_cr >= 0.3)
  aki_df$aki_kdigo2 <- as.integer(aki_df$cr_pre>0 & aki_df$cr_post >= 2.0*aki_df$cr_pre)
  aki_df$aki_kdigo3 <- as.integer((aki_df$cr_pre>0 & aki_df$cr_post >= 3.0*aki_df$cr_pre) | aki_df$cr_post >= 4.0)

  cat(sprintf("  AKI subset (24h): %d treated + %d control\n",
              sum(aki_df$treated==1), sum(aki_df$treated==0)))

  # ── Outcomes ────────────────────────────────────────────────────────
  # AKI outcomes use aki_df; binary outcomes use full
  outcomes <- list(
    list(name="AKI KDIGO>=1", col="aki_kdigo1", data="aki"),
    list(name="AKI KDIGO>=2", col="aki_kdigo2", data="aki"),
    list(name="AKI KDIGO>=3", col="aki_kdigo3", data="aki"),
    list(name="Hospital mortality", col="hosp_mortality", data="full"),
    list(name="POAF (new-onset)", col="poaf", data="full"),
    list(name="Encephalopathy", col="encephalopathy", data="full"),
    list(name="Vent arrhythmia", col="vent_arrhythmia", data="full")
  )

  # ── Subgroup definitions (var, levels, ref) ─────────────────────────
  subgroup_defs <- list(
    list(var="age_cat", label="Age", ref="<65",
         build=function(d){ d$age_cat <- ifelse(d$age<65,"<65",">=65"); d }),
    list(var="is_female", label="Sex", ref="0",
         build=function(d){ d$is_female <- as.character(d$is_female); d }),
    list(var="egfr_cat", label="eGFR", ref=">=90",
         build=function(d){
           d$egfr_cat <- cut(d$egfr, c(-Inf,45,60,90,Inf), labels=c("<45","45-59","60-89",">=90"), right=FALSE)
           d$egfr_cat <- as.character(d$egfr_cat); d }),
    list(var="cr_cat", label="Baseline Cr", ref="<=1.0",
         build=function(d){ d$cr_cat <- ifelse(d$cr_pre<=1.0,"<=1.0",">1.0"); d }),
    list(var="mg_cat", label="Serum Mg", ref="1.8-2.0",
         build=function(d){
           d$mg_cat <- cut(d$first_mg_value, c(-Inf,1.8,2.0,2.3,Inf),
                           labels=c("<1.8","1.8-2.0","2.0-2.3",">2.3"), right=FALSE)
           d$mg_cat <- as.character(d$mg_cat); d }),
    list(var="surgery_type", label="Surgery", ref="other_cardiac",
         build=function(d){ d }),
    list(var="diabetes", label="Diabetes", ref="0",
         build=function(d){ d$diabetes <- as.character(d$diabetes); d }),
    list(var="ckd", label="CKD", ref="0",
         build=function(d){ d$ckd <- as.character(d$ckd); d }),
    list(var="heart_failure", label="Heart failure", ref="0",
         build=function(d){ d$heart_failure <- as.character(d$heart_failure); d }),
    list(var="bmi_cat", label="BMI", ref="25-30",
         build=function(d){
           d$bmi_cat <- cut(d$bmi, c(-Inf,25,30,Inf), labels=c("<25","25-30",">=30"), right=FALSE)
           d$bmi_cat <- as.character(d$bmi_cat); d })
  )

  # ── Run sweep ──────────────────────────────────────────────────────
  all_rows <- list(); ridx <- 0

  for (oc in outcomes) {
    d <- if(oc$data=="aki") aki_df else full
    if (!oc$col %in% names(d)) next

    cat(sprintf("\n%s\n%s\n%s\n", paste(rep("-",60),collapse=""), oc$name,
                paste(rep("-",60),collapse="")))
    cat("  subgroup              level       n_trt  n_ctl  rate_trt rate_ctl  OR (95%CI)          P     ref  int_P\n")
    cat("  ───────────────────  ──────────  ─────  ─────  ─────── ────────  ──────────────────  ──────  ───  ─────\n")

    # Overall
    res <- weighted_or(d, oc$col)
    if (!is.null(res)) {
      ridx <- ridx+1
      sig <- if(res$p<0.05) " *" else ""
      cat(sprintf("  %-21s %-10s  %5d  %5d   %5.1f%%   %5.1f%%  %.2f (%.2f-%.2f)  %.4f%s\n",
                  "Overall","—",res$n_trt,res$n_ctl,100*res$rate_trt,100*res$rate_ctl,
                  res$or,res$or_lo,res$or_hi,res$p,sig))
      all_rows[[ridx]] <- data.frame(outcome=oc$name, subgroup="Overall", level="—",
        n_trt=res$n_trt, n_ctl=res$n_ctl, rate_trt=round(100*res$rate_trt,1),
        rate_ctl=round(100*res$rate_ctl,1), or=round(res$or,2),
        or_lo=round(res$or_lo,2), or_hi=round(res$or_hi,2),
        p=round(res$p,4), ref="", int_p=NA, stringsAsFactors=F)
    }

    # Per subgroup
    for (sg in subgroup_defs) {
      # cr_cat only for AKI outcomes
      if (sg$var=="cr_cat" && oc$data!="aki") next

      d2 <- tryCatch(sg$build(d), error=function(e) NULL)
      if (is.null(d2) || !sg$var %in% names(d2)) next
      d2 <- d2[!is.na(d2[[sg$var]]),]

      levels <- sort(unique(d2[[sg$var]]))
      if (length(levels)<2) next

      # Interaction P (treated × subgroup as factor)
      d2$sg_fac <- factor(d2[[sg$var]], levels=c(sg$ref, setdiff(levels, sg$ref)))
      int_p <- interaction_p(d2, oc$col, "sg_fac")
      int_str <- if(!is.na(int_p)) sprintf("%.3f",int_p) else "  —"

      for (lv in levels) {
        sub <- d2[d2[[sg$var]]==lv,]
        res <- weighted_or(sub, oc$col)
        is_ref <- (lv == sg$ref)
        ref_str <- if(is_ref) "ref" else ""
        ridx <- ridx+1

        if (!is.null(res)) {
          sig <- if(res$p<0.05) " *" else ""
          # Print int_p only on first level
          ip_str <- if(lv==levels[1]) int_str else "     "
          cat(sprintf("  %-21s %-10s  %5d  %5d   %5.1f%%   %5.1f%%  %.2f (%.2f-%.2f)  %.4f%s %s  %s\n",
                      sg$label, lv, res$n_trt, res$n_ctl,
                      100*res$rate_trt, 100*res$rate_ctl,
                      res$or, res$or_lo, res$or_hi, res$p, sig, ref_str, ip_str))
          all_rows[[ridx]] <- data.frame(outcome=oc$name, subgroup=sg$label, level=lv,
            n_trt=res$n_trt, n_ctl=res$n_ctl, rate_trt=round(100*res$rate_trt,1),
            rate_ctl=round(100*res$rate_ctl,1), or=round(res$or,2),
            or_lo=round(res$or_lo,2), or_hi=round(res$or_hi,2),
            p=round(res$p,4), ref=ref_str,
            int_p=if(lv==levels[1] && !is.na(int_p)) round(int_p,4) else NA,
            stringsAsFactors=F)
        } else {
          ip_str <- if(lv==levels[1]) int_str else "     "
          cat(sprintf("  %-21s %-10s  %5d  %5d      —        —          —               — %s  %s\n",
                      sg$label, lv, sum(sub$treated==1), sum(sub$treated==0), ref_str, ip_str))
          all_rows[[ridx]] <- data.frame(outcome=oc$name, subgroup=sg$label, level=lv,
            n_trt=sum(sub$treated==1), n_ctl=sum(sub$treated==0),
            rate_trt=NA, rate_ctl=NA, or=NA, or_lo=NA, or_hi=NA, p=NA,
            ref=ref_str, int_p=NA, stringsAsFactors=F)
        }
      }
    }
  }

  results <- do.call(rbind, all_rows)
  write.csv(results, file.path(RESULTS,sprintf("did_subgroups_full_%s.csv",tag)), row.names=F)
  cat(sprintf("\n  Saved: did_subgroups_full_%s.csv (%d rows)\n", tag, nrow(results)))
}

# ============================================================================
cat("======================================================================\n")
cat("01b — Downstream subgroup analysis (sIPTW_DR)\n")
cat("  IPTW weights computed ONCE on full sample, then stratified\n")
cat("  Outcomes: AKI KDIGO 1/2/3 + mortality + POAF + encephalopathy + vent arr\n")
cat("  Subgroups: age, sex, eGFR, Cr, Mg, surgery, DM, CKD, HF, BMI\n")
cat("  Each subgroup: OR + interaction P, with explicit reference group\n")
cat("======================================================================\n")

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) { cat("Usage: Rscript 01b_did_aki_subgroups.R eicu|mimic\n"); quit(status=1) }
for (a in args) run_sweep(toupper(a))
