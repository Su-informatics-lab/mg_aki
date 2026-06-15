#!/usr/bin/env Rscript
# ============================================================================
# 08b_hospital_re.R — Hospital Random Effects (The Hospital Test)
#
# The >2.3 stratum signal in eICU could reflect hospital-level practice
# variation: hospitals that supplement everyone AND have better AKI outcomes
# for unrelated reasons. This script tests that by:
#
#   1. Mixed-effects logistic with hospital random intercept (glmer)
#   2. ICC: how much AKI variation is between vs within hospitals?
#   3. Hospital exposure variation: is supplementation clustered?
#   4. Narrow sub-bands within >2.3 (2.3-2.6, 2.6-3.0, >3.0)
#   5. Within-hospital meta-analysis (hospitals with ≥5 per arm)
#
# If the treatment OR survives hospital RE → hospital confounding unlikely
# If the treatment OR → null under hospital RE → hospital confounding
#
# Run: Rscript 08b_hospital_re.R
# ============================================================================

suppressPackageStartupMessages({
  library(lme4); library(sandwich); library(lmtest)
})
RESULTS <- path.expand("~/mg_aki/results")

# ── Standardize ──────────────────────────────────────────────────────────
stdz <- function(d) {
  rmap <- c(mg_supplementation="mg_supp", age_num="age",
    baseline_cr="baseline_creatinine", baseline_egfr="egfr",
    hx_chf="heart_failure", hx_hypertension="hypertension",
    hx_diabetes="diabetes", hx_ckd="ckd", hx_copd="copd",
    hx_pvd="pvd", hx_stroke="stroke", hx_liver="liver_disease",
    nephrotox_loop_diuretic="loop_diuretics", nephrotox_nsaid="nsaids",
    nephrotox_acei_arb="acei_arb", nephrotox_ppi="ppi",
    has_betablocker="beta_blockers", has_steroid="steroids",
    preop_antiarrhythmic="antiarrhythmics",
    first_k_value="first_potassium", first_ca_value="first_calcium",
    first_hr="first_heartrate", has_vasopressor="vasopressor_6h")
  for (old in names(rmap)) {
    new <- rmap[[old]]
    if (old %in% names(d) && !new %in% names(d)) names(d)[names(d)==old] <- new
  }
  if (is.character(d$age)) { d$age <- suppressWarnings(as.numeric(d$age)); d$age[is.na(d$age)] <- 90 }
  if ("surgery_type" %in% names(d)) {
    d$surg_cabg <- as.integer(d$surgery_type=="cabg")
    d$surg_valve <- as.integer(d$surgery_type=="valve")
    d$surg_combined <- as.integer(d$surgery_type=="combined")
  }
  if ("first_lactate" %in% names(d)) {
    d$lactate_missing <- as.integer(is.na(d$first_lactate))
    d$first_lactate[is.na(d$first_lactate)] <- median(d$first_lactate, na.rm=TRUE)
  }
  for (v in c("bmi","first_heartrate","first_calcium","first_potassium"))
    if (v %in% names(d) && any(is.na(d[[v]])))
      d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm=TRUE)
  d
}

# PS covariates (without first_mg_value — stratified on it)
ps_covars <- c("age","is_female","bmi",
  "surg_cabg","surg_valve","surg_combined",
  "heart_failure","hypertension","diabetes","ckd",
  "copd","pvd","stroke","liver_disease",
  "baseline_creatinine","egfr",
  "loop_diuretics","nsaids","acei_arb","ppi",
  "beta_blockers","steroids","antiarrhythmics",
  "first_potassium","first_calcium","first_heartrate",
  "vasopressor_6h")

# ============================================================================
# LOAD DATA
# ============================================================================
cat(strrep("=",65), "\n")
cat("08b: HOSPITAL RANDOM EFFECTS (The Hospital Test)\n")
cat(strrep("=",65), "\n")

dat <- stdz(read.csv(file.path(RESULTS, "01_analysis_a_cohort.csv"), stringsAsFactors=FALSE))
if (!"hospitalid" %in% names(dat)) stop("hospitalid not found in eICU cohort")

avail_covars <- intersect(ps_covars, names(dat))
if ("first_lactate" %in% names(dat)) avail_covars <- c(avail_covars, "first_lactate", "lactate_missing")
dat <- dat[complete.cases(dat[, avail_covars]),]

# Mg strata
dat$mg_stratum <- cut(dat$first_mg_value,
  breaks=c(0, 1.8, 2.0, 2.3, Inf),
  labels=c("<1.8","1.8-2.0","2.0-2.3",">2.3"),
  right=FALSE, include.lowest=TRUE)

cat(sprintf("  N=%d, hospitals=%d\n", nrow(dat), length(unique(dat$hospitalid))))

results <- list()

# ============================================================================
# 1. ICC: How much AKI variation is between hospitals?
# ============================================================================
cat(sprintf("\n%s\n1. INTRACLASS CORRELATION (AKI ~ 1 + (1|hospital))\n%s\n",
            strrep("-",65), strrep("-",65)))

icc_fit <- glmer(aki_kdigo1 ~ 1 + (1 | hospitalid), data=dat, family=binomial,
                 nAGQ=1, control=glmerControl(optimizer="bobyqa"))
vc <- as.data.frame(VarCorr(icc_fit))
sigma2_h <- vc$vcov[1]
icc <- sigma2_h / (sigma2_h + pi^2/3)
cat(sprintf("  Hospital variance (logit scale): %.4f\n", sigma2_h))
cat(sprintf("  ICC: %.3f (%.1f%% of AKI variation is between hospitals)\n", icc, 100*icc))
cat(sprintf("  Interpretation: %s\n",
    ifelse(icc > 0.05, "SUBSTANTIAL hospital-level clustering — hospital confounding plausible",
           "Modest clustering — hospital confounding less concerning")))

# ============================================================================
# 2. Hospital exposure variation: is Mg supplementation clustered?
# ============================================================================
cat(sprintf("\n%s\n2. HOSPITAL EXPOSURE VARIATION\n%s\n", strrep("-",65), strrep("-",65)))

hosp_stats <- aggregate(mg_supp ~ hospitalid, data=dat, FUN=function(x) {
  c(n=length(x), trt=sum(x), rate=mean(x))
})
hosp_stats <- do.call(data.frame, hosp_stats)
names(hosp_stats) <- c("hospitalid","n","n_trt","trt_rate")
hosp_stats <- hosp_stats[hosp_stats$n >= 10,]

cat(sprintf("  Hospitals with >=10 patients: %d\n", nrow(hosp_stats)))
cat(sprintf("  Supplementation rate across hospitals:\n"))
cat(sprintf("    Median: %.1f%%  IQR: [%.1f%%, %.1f%%]\n",
    100*median(hosp_stats$trt_rate), 100*quantile(hosp_stats$trt_rate, 0.25),
    100*quantile(hosp_stats$trt_rate, 0.75)))
cat(sprintf("    Range:  %.1f%% – %.1f%%\n",
    100*min(hosp_stats$trt_rate), 100*max(hosp_stats$trt_rate)))
n_zero <- sum(hosp_stats$trt_rate == 0)
n_high <- sum(hosp_stats$trt_rate > 0.30)
cat(sprintf("    Zero supplementation: %d hospitals (%.0f%%)\n",
    n_zero, 100*n_zero/nrow(hosp_stats)))
cat(sprintf("    >30%% supplementation: %d hospitals (%.0f%%)\n",
    n_high, 100*n_high/nrow(hosp_stats)))

# Correlation: hospital supplementation rate vs hospital AKI rate
hosp_aki <- aggregate(aki_kdigo1 ~ hospitalid, data=dat, FUN=mean)
names(hosp_aki) <- c("hospitalid","aki_rate")
hosp_merged <- merge(hosp_stats, hosp_aki)
hosp_merged <- hosp_merged[hosp_merged$n >= 20,]
if (nrow(hosp_merged) >= 10) {
  cr <- cor.test(hosp_merged$trt_rate, hosp_merged$aki_rate)
  cat(sprintf("\n  Hospital-level correlation (supplementation rate vs AKI rate):\n"))
  cat(sprintf("    r = %.3f, P = %.4f (among %d hospitals with >=20 pts)\n",
      cr$estimate, cr$p.value, nrow(hosp_merged)))
  cat(sprintf("    %s\n",
      ifelse(cr$estimate < -0.15 && cr$p.value < 0.05,
             "NEGATIVE correlation: high-supplementation hospitals have lower AKI → SUSPICIOUS",
             ifelse(cr$estimate > 0.15 && cr$p.value < 0.05,
                    "POSITIVE correlation: high-supplementation hospitals have higher AKI",
                    "No significant correlation — hospital confounding less likely at this level"))))
}

# ============================================================================
# 3. MIXED-EFFECTS MODELS (each Mg stratum + overall)
# ============================================================================
cat(sprintf("\n%s\n3. MIXED-EFFECTS LOGISTIC REGRESSION\n   aki_kdigo1 ~ mg_supp + covariates + (1|hospitalid)\n%s\n",
            strrep("-",65), strrep("-",65)))

cov_str <- paste(avail_covars, collapse=" + ")

for (stratum in c("Overall", "<1.8", "1.8-2.0", "2.0-2.3", ">2.3")) {
  d_s <- if (stratum == "Overall") dat else dat[dat$mg_stratum == stratum,]
  n_trt <- sum(d_s$mg_supp); n_ctrl <- nrow(d_s) - n_trt
  n_hosp <- length(unique(d_s$hospitalid))

  if (n_trt < 15 || n_ctrl < 15 || n_hosp < 5) {
    cat(sprintf("  %-10s SKIP (trt=%d, ctrl=%d, hosp=%d)\n", stratum, n_trt, n_ctrl, n_hosp))
    next
  }

  # Drop hospitals with <2 patients (glmer needs variance within clusters)
  hosp_n <- table(d_s$hospitalid)
  d_s <- d_s[d_s$hospitalid %in% names(hosp_n[hosp_n >= 2]),]

  cat(sprintf("\n  ── %s (N=%d, trt=%d, hosp=%d) ──\n", stratum, nrow(d_s), sum(d_s$mg_supp),
              length(unique(d_s$hospitalid))))

  # 3a. Standard logistic (no hospital adjustment) — reference
  fml_fixed <- as.formula(paste("aki_kdigo1 ~ mg_supp +", cov_str))
  tryCatch({
    fit_fixed <- glm(fml_fixed, data=d_s, family=binomial())
    ct_f <- coeftest(fit_fixed, vcov.=vcovCL(fit_fixed, cluster=d_s$hospitalid))
    or_f <- exp(ct_f[2,1]); lo_f <- exp(ct_f[2,1]-1.96*ct_f[2,2])
    hi_f <- exp(ct_f[2,1]+1.96*ct_f[2,2]); p_f <- 2*pnorm(-abs(ct_f[2,1]/ct_f[2,2]))
    cat(sprintf("    Fixed (cluster-robust):  OR %.3f (%.3f–%.3f) P=%.4f\n",
                or_f, lo_f, hi_f, p_f))
    results[[length(results)+1]] <- data.frame(
      stratum=stratum, model="fixed_clusterSE",
      n=nrow(d_s), n_trt=sum(d_s$mg_supp), n_hosp=length(unique(d_s$hospitalid)),
      or=round(or_f,3), lo=round(lo_f,3), hi=round(hi_f,3), p=round(p_f,4))
  }, error=function(e) cat(sprintf("    Fixed model failed: %s\n", e$message)))

  # 3b. Mixed-effects logistic with hospital random intercept
  fml_re <- as.formula(paste("aki_kdigo1 ~ mg_supp +", cov_str, "+ (1 | hospitalid)"))
  tryCatch({
    fit_re <- glmer(fml_re, data=d_s, family=binomial, nAGQ=1,
                    control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=50000)))
    sf <- summary(fit_re)
    ct_re <- sf$coefficients
    tr <- which(rownames(ct_re) == "mg_supp")
    or_re <- exp(ct_re[tr,1]); lo_re <- exp(ct_re[tr,1]-1.96*ct_re[tr,2])
    hi_re <- exp(ct_re[tr,1]+1.96*ct_re[tr,2])
    p_re <- 2*pnorm(-abs(ct_re[tr,1]/ct_re[tr,2]))
    sig <- ifelse(p_re < 0.05, " *", "")

    # Hospital variance in this stratum
    vc_s <- as.data.frame(VarCorr(fit_re))
    sigma2_s <- vc_s$vcov[1]
    icc_s <- sigma2_s / (sigma2_s + pi^2/3)

    cat(sprintf("    Hospital RE (glmer):     OR %.3f (%.3f–%.3f) P=%.4f%s\n",
                or_re, lo_re, hi_re, p_re, sig))
    cat(sprintf("    Hospital variance: %.4f  ICC: %.3f\n", sigma2_s, icc_s))

    # Compare with fixed
    if (exists("or_f")) {
      change <- (or_re - or_f) / or_f * 100
      cat(sprintf("    Change from fixed: %+.1f%% (%s)\n", change,
          ifelse(abs(change) < 10, "MINIMAL — hospital RE barely changes estimate",
                 ifelse(or_re > or_f, "ATTENUATED toward null — hospital confounding present",
                        "STRENGTHENED — hospital RE makes estimate more protective"))))
    }

    results[[length(results)+1]] <- data.frame(
      stratum=stratum, model="hospital_RE",
      n=nrow(d_s), n_trt=sum(d_s$mg_supp), n_hosp=length(unique(d_s$hospitalid)),
      or=round(or_re,3), lo=round(lo_re,3), hi=round(hi_re,3), p=round(p_re,4))
  }, error=function(e) cat(sprintf("    Hospital RE FAILED: %s\n", e$message)))
}

# ============================================================================
# 4. NARROW SUB-BANDS WITHIN >2.3
# ============================================================================
cat(sprintf("\n%s\n4. NARROW SUB-BANDS WITHIN >2.3 (where is the signal?)\n%s\n",
            strrep("-",65), strrep("-",65)))

d_hi <- dat[dat$mg_stratum == ">2.3",]
d_hi$mg_subband <- cut(d_hi$first_mg_value,
  breaks=c(2.3, 2.6, 3.0, Inf),
  labels=c("2.3-2.6","2.6-3.0",">3.0"),
  right=FALSE, include.lowest=TRUE)

for (sb in c("2.3-2.6","2.6-3.0",">3.0")) {
  d_sb <- d_hi[d_hi$mg_subband == sb,]
  n_trt <- sum(d_sb$mg_supp); n_ctrl <- nrow(d_sb) - n_trt
  if (n_trt < 10 || n_ctrl < 10) {
    cat(sprintf("  %-10s N=%d trt=%d — SKIP\n", sb, nrow(d_sb), n_trt))
    next
  }
  fml <- as.formula(paste("aki_kdigo1 ~ mg_supp +", cov_str))
  tryCatch({
    fit <- glm(fml, data=d_sb, family=binomial())
    ct <- coeftest(fit, vcov.=vcovHC(fit, type="HC1"))
    or <- exp(ct[2,1]); lo <- exp(ct[2,1]-1.96*ct[2,2])
    hi <- exp(ct[2,1]+1.96*ct[2,2]); p <- 2*pnorm(-abs(ct[2,1]/ct[2,2]))
    sig <- ifelse(p < 0.05, " *", "")
    cat(sprintf("  %-10s N=%4d trt=%3d  OR %.3f (%.3f–%.3f) P=%.4f  AKI: trt=%.0f%% ctrl=%.0f%%%s\n",
        sb, nrow(d_sb), n_trt, or, lo, hi, p,
        100*mean(d_sb$aki_kdigo1[d_sb$mg_supp==1]),
        100*mean(d_sb$aki_kdigo1[d_sb$mg_supp==0]), sig))
    results[[length(results)+1]] <- data.frame(
      stratum=paste0(">2.3:",sb), model="fixed_subband",
      n=nrow(d_sb), n_trt=n_trt, n_hosp=length(unique(d_sb$hospitalid)),
      or=round(or,3), lo=round(lo,3), hi=round(hi,3), p=round(p,4))
  }, error=function(e) cat(sprintf("  %-10s FAILED: %s\n", sb, e$message)))
}

# ============================================================================
# 5. WITHIN-HOSPITAL META-ANALYSIS (>2.3 stratum)
# ============================================================================
cat(sprintf("\n%s\n5. WITHIN-HOSPITAL META-ANALYSIS (>2.3 stratum)\n%s\n",
            strrep("-",65), strrep("-",65)))

d_hi <- dat[dat$mg_stratum == ">2.3",]
hosp_counts <- aggregate(mg_supp ~ hospitalid, data=d_hi,
                         FUN=function(x) c(n=length(x), trt=sum(x), ctrl=sum(x==0)))
hosp_counts <- do.call(data.frame, hosp_counts)
names(hosp_counts) <- c("hospitalid","n","n_trt","n_ctrl")
eligible_hosp <- hosp_counts$hospitalid[hosp_counts$n_trt >= 5 & hosp_counts$n_ctrl >= 5]
cat(sprintf("  Hospitals with >=5 per arm in >2.3 stratum: %d\n", length(eligible_hosp)))

if (length(eligible_hosp) >= 3) {
  hosp_ors <- list()
  for (h in eligible_hosp) {
    d_h <- d_hi[d_hi$hospitalid == h,]
    # Simple adjusted model within hospital (fewer covariates to avoid separation)
    simple_covars <- intersect(c("age","is_female","baseline_creatinine","surg_valve",
                                  "surg_combined","vasopressor_6h"), names(d_h))
    # Check for separation
    if (length(unique(d_h$aki_kdigo1)) < 2) next
    if (length(unique(d_h$mg_supp)) < 2) next
    fml_h <- as.formula(paste("aki_kdigo1 ~ mg_supp +", paste(simple_covars, collapse="+")))
    tryCatch({
      fit_h <- glm(fml_h, data=d_h, family=binomial())
      ct_h <- coef(summary(fit_h))
      if ("mg_supp" %in% rownames(ct_h)) {
        hosp_ors[[length(hosp_ors)+1]] <- data.frame(
          hospitalid=h, n=nrow(d_h), n_trt=sum(d_h$mg_supp),
          logOR=ct_h["mg_supp",1], se=ct_h["mg_supp",2])
      }
    }, error=function(e) {})
  }

  if (length(hosp_ors) >= 3) {
    ho <- do.call(rbind, hosp_ors)
    # Remove extreme estimates (likely separation artifacts)
    ho <- ho[abs(ho$logOR) < 5 & ho$se < 5,]
    cat(sprintf("  Usable hospital-level estimates: %d\n", nrow(ho)))

    # Fixed-effects meta across hospitals
    w <- 1 / ho$se^2
    pool_logOR <- sum(w * ho$logOR) / sum(w)
    pool_se <- sqrt(1 / sum(w))
    pool_or <- exp(pool_logOR)
    pool_lo <- exp(pool_logOR - 1.96*pool_se)
    pool_hi <- exp(pool_logOR + 1.96*pool_se)
    pool_p <- 2*pnorm(-abs(pool_logOR/pool_se))

    Q <- sum(w * (ho$logOR - pool_logOR)^2)
    df <- nrow(ho) - 1
    I2 <- max(0, (Q-df)/Q)*100

    cat(sprintf("  Within-hospital pooled OR: %.3f (%.3f–%.3f) P=%.4f  I²=%.0f%%\n",
                pool_or, pool_lo, pool_hi, pool_p, I2))
    cat(sprintf("  Interpretation: %s\n",
        ifelse(pool_or < 0.85 && pool_p < 0.10,
               "Protective signal persists WITHIN hospitals → NOT hospital confounding",
               ifelse(pool_or > 0.95,
                      "Signal VANISHES within hospitals → likely hospital confounding",
                      "Attenuated but directional — inconclusive"))))

    results[[length(results)+1]] <- data.frame(
      stratum=">2.3", model="within_hospital_meta",
      n=sum(ho$n), n_trt=sum(ho$n_trt), n_hosp=nrow(ho),
      or=round(pool_or,3), lo=round(pool_lo,3), hi=round(pool_hi,3), p=round(pool_p,4))
  }
} else {
  cat("  Insufficient hospitals for within-hospital meta-analysis\n")
}

# ============================================================================
# SAVE + VERDICT
# ============================================================================
all_res <- do.call(rbind, results)
outpath <- file.path(RESULTS, "08b_hospital_re.csv")
write.csv(all_res, outpath, row.names=FALSE)

cat(sprintf("\n%s\nVERDICT SUMMARY\n%s\n", strrep("=",65), strrep("=",65)))
cat("  Compare these columns:\n")
cat("    fixed_clusterSE = original analysis (cluster-robust SEs)\n")
cat("    hospital_RE     = mixed model with hospital random intercept\n\n")

for (s in unique(all_res$stratum)) {
  fixed <- all_res[all_res$stratum==s & all_res$model=="fixed_clusterSE",]
  re    <- all_res[all_res$stratum==s & all_res$model=="hospital_RE",]
  if (nrow(fixed)>0 && nrow(re)>0) {
    cat(sprintf("  %-12s Fixed: OR %.3f P=%.4f  →  Hospital RE: OR %.3f P=%.4f  [Δ%+.0f%%]\n",
        s, fixed$or, fixed$p, re$or, re$p, (re$or-fixed$or)/fixed$or*100))
  }
}

cat(sprintf("\n  If all Δ are small (<10%%): hospital confounding NOT driving the signal\n"))
cat(sprintf("  If >2.3 Δ is large (>20%%): hospital confounding IS the explanation\n"))

cat(sprintf("\n✓ Saved: %s (%d rows)\n", outpath, nrow(all_res)))
