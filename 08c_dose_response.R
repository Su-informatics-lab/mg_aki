#!/usr/bin/env Rscript
# ============================================================================
# 08c_dose_response.R — MIMIC Dose-Response Analysis
#
# MIMIC has mg_total_dose (total mg of MgSO4 administered within 6h).
# This is a within-treated analysis: among patients who received Mg,
# does higher dose predict lower AKI?
#
# This sidesteps treated-vs-untreated confounding entirely.
# If a dose-response exists, it strongly supports a causal/threshold effect.
#
# Also: combine with serum Mg data to show dose → achieved level → AKI
#       (the full causal chain)
#
# Run: Rscript 08c_dose_response.R
# ============================================================================

suppressPackageStartupMessages({
  library(sandwich); library(lmtest)
})
RESULTS <- path.expand("~/mg_aki/results")

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
  if ("surgery_type" %in% names(d)) {
    d$surg_cabg <- as.integer(d$surgery_type=="cabg")
    d$surg_valve <- as.integer(d$surgery_type=="valve")
    d$surg_combined <- as.integer(d$surgery_type=="combined")
  }
  for (v in c("bmi","first_heartrate","first_calcium","first_potassium"))
    if (v %in% names(d) && any(is.na(d[[v]])))
      d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm=TRUE)
  d
}

cat(strrep("=", 65), "\n")
cat("08c: MIMIC DOSE-RESPONSE ANALYSIS\n")
cat(strrep("=", 65), "\n")

dat <- stdz(read.csv(file.path(RESULTS, "04_mimic_cohort.csv"), stringsAsFactors=FALSE))

if (!"mg_total_dose" %in% names(dat)) stop("mg_total_dose not found in MIMIC cohort")

# ── Focus on supplemented patients ───────────────────────────────────
trt <- dat[dat$mg_supp == 1 & !is.na(dat$mg_total_dose) & dat$mg_total_dose > 0, ]
cat(sprintf("\n  Supplemented patients with dose data: %d\n", nrow(trt)))
cat(sprintf("  Dose (mg MgSO4): median=%.0f, IQR=[%.0f, %.0f], range=[%.0f, %.0f]\n",
    median(trt$mg_total_dose), quantile(trt$mg_total_dose, 0.25),
    quantile(trt$mg_total_dose, 0.75), min(trt$mg_total_dose), max(trt$mg_total_dose)))
cat(sprintf("  AKI rate in supplemented: %.1f%%\n", 100*mean(trt$aki_kdigo1)))

# ── Dose tertiles ────────────────────────────────────────────────────
trt$dose_tertile <- cut(trt$mg_total_dose,
  breaks=quantile(trt$mg_total_dose, c(0, 1/3, 2/3, 1), na.rm=TRUE),
  labels=c("T1 (low)", "T2 (mid)", "T3 (high)"),
  include.lowest=TRUE)

cat(sprintf("\n%s\n1. DOSE TERTILE ANALYSIS\n%s\n", strrep("-",65), strrep("-",65)))
for (tert in c("T1 (low)", "T2 (mid)", "T3 (high)")) {
  sub <- trt[trt$dose_tertile == tert, ]
  cat(sprintf("  %-12s N=%3d  dose: %.0f–%.0f mg  AKI=%.1f%%  Mg: %.2f±%.2f\n",
      tert, nrow(sub),
      min(sub$mg_total_dose), max(sub$mg_total_dose),
      100*mean(sub$aki_kdigo1),
      mean(sub$first_mg_value), sd(sub$first_mg_value)))
}

# ── Dose-response trend test (adjusted) ──────────────────────────────
cat(sprintf("\n%s\n2. ADJUSTED DOSE-RESPONSE (continuous dose → AKI)\n%s\n",
            strrep("-",65), strrep("-",65)))

adj_vars <- intersect(c("age","is_female","bmi",
  "surg_cabg","surg_valve","surg_combined",
  "heart_failure","hypertension","diabetes","ckd","copd","pvd","stroke","liver_disease",
  "baseline_creatinine","egfr",
  "loop_diuretics","nsaids","acei_arb","ppi",
  "beta_blockers","steroids","antiarrhythmics",
  "first_potassium","first_calcium","first_heartrate",
  "vasopressor_6h","first_mg_value"), names(trt))

# Standardize dose for interpretability (per 1g = 1000mg)
trt$dose_g <- trt$mg_total_dose / 1000

# Continuous dose
fml_cont <- as.formula(paste("aki_kdigo1 ~ dose_g +", paste(adj_vars, collapse="+")))
trt_cc <- trt[complete.cases(trt[, c("dose_g", adj_vars, "aki_kdigo1")]),]
cat(sprintf("  Complete cases: %d\n", nrow(trt_cc)))

tryCatch({
  fit_cont <- glm(fml_cont, data=trt_cc, family=binomial())
  ct <- coeftest(fit_cont, vcov.=vcovHC(fit_cont, type="HC1"))
  dose_row <- which(rownames(ct) == "dose_g")
  or_per_g <- exp(ct[dose_row, 1])
  lo <- exp(ct[dose_row, 1] - 1.96*ct[dose_row, 2])
  hi <- exp(ct[dose_row, 1] + 1.96*ct[dose_row, 2])
  p <- 2*pnorm(-abs(ct[dose_row, 1]/ct[dose_row, 2]))
  cat(sprintf("  OR per 1g MgSO4 increase: %.3f (%.3f–%.3f) P=%.4f\n", or_per_g, lo, hi, p))
  cat(sprintf("  Interpretation: %s\n",
      ifelse(or_per_g < 0.9 && p < 0.10,
             "PROTECTIVE dose-response — higher dose, lower AKI",
             ifelse(or_per_g > 1.1,
                    "HARMFUL dose-response — higher dose, more AKI (confounding by severity?)",
                    "No significant dose-response"))))
}, error=function(e) cat(sprintf("  FAILED: %s\n", e$message)))

# Tertile contrasts (T3 vs T1, adjusted)
cat(sprintf("\n  Tertile contrasts:\n"))
trt_cc$dose_t3 <- as.integer(trt_cc$dose_tertile == "T3 (high)")
trt_t1t3 <- trt_cc[trt_cc$dose_tertile %in% c("T1 (low)", "T3 (high)"),]
if (nrow(trt_t1t3) > 30) {
  fml_t <- as.formula(paste("aki_kdigo1 ~ dose_t3 +", paste(adj_vars, collapse="+")))
  tryCatch({
    fit_t <- glm(fml_t, data=trt_t1t3, family=binomial())
    ct_t <- coeftest(fit_t, vcov.=vcovHC(fit_t, type="HC1"))
    or_t <- exp(ct_t[2,1]); lo_t <- exp(ct_t[2,1]-1.96*ct_t[2,2])
    hi_t <- exp(ct_t[2,1]+1.96*ct_t[2,2]); p_t <- 2*pnorm(-abs(ct_t[2,1]/ct_t[2,2]))
    cat(sprintf("  T3 (high) vs T1 (low): OR %.3f (%.3f–%.3f) P=%.4f\n", or_t, lo_t, hi_t, p_t))
  }, error=function(e) cat(sprintf("  FAILED: %s\n", e$message)))
}

# ── Dose → follow-up serum Mg (the pharmacokinetic link) ─────────────
cat(sprintf("\n%s\n3. DOSE → ACHIEVED SERUM Mg (pharmacokinetic check)\n%s\n",
            strrep("-",65), strrep("-",65)))

if ("followup_mg_value" %in% names(trt)) {
  trt_fu <- trt[!is.na(trt$followup_mg_value),]
  cat(sprintf("  Patients with follow-up Mg (6-48h): %d\n", nrow(trt_fu)))

  for (tert in c("T1 (low)", "T2 (mid)", "T3 (high)")) {
    sub <- trt_fu[trt_fu$dose_tertile == tert,]
    if (nrow(sub) > 5) {
      cat(sprintf("  %-12s N=%3d  follow-up Mg: %.2f±%.2f  delta_Mg: %+.2f\n",
          tert, nrow(sub),
          mean(sub$followup_mg_value), sd(sub$followup_mg_value),
          mean(sub$delta_mg, na.rm=TRUE)))
    }
  }

  # Does higher dose → higher follow-up Mg?
  if (nrow(trt_fu) > 30) {
    cr_dose_mg <- cor.test(trt_fu$mg_total_dose, trt_fu$followup_mg_value)
    cat(sprintf("\n  Dose–followup Mg correlation: r=%.3f P=%.4f\n",
        cr_dose_mg$estimate, cr_dose_mg$p.value))
  }
} else {
  cat("  No follow-up Mg data available\n")
}

# ── Follow-up serum Mg → AKI (the threshold test) ───────────────────
cat(sprintf("\n%s\n4. ACHIEVED SERUM Mg → AKI (the threshold test)\n%s\n",
            strrep("-",65), strrep("-",65)))

if ("followup_mg_value" %in% names(trt)) {
  trt_fu <- trt[!is.na(trt$followup_mg_value),]
  if (nrow(trt_fu) > 30) {
    # Achieved Mg tertiles
    trt_fu$fu_mg_cat <- cut(trt_fu$followup_mg_value,
      breaks=c(0, 2.0, 2.3, Inf),
      labels=c("<2.0", "2.0-2.3", ">2.3"))

    cat("  AKI by achieved follow-up Mg level (within supplemented):\n")
    for (cat_label in c("<2.0", "2.0-2.3", ">2.3")) {
      sub <- trt_fu[trt_fu$fu_mg_cat == cat_label,]
      if (nrow(sub) > 5) {
        cat(sprintf("    %-10s N=%3d  AKI=%.1f%%  Mg: %.2f±%.2f\n",
            cat_label, nrow(sub), 100*mean(sub$aki_kdigo1),
            mean(sub$followup_mg_value), sd(sub$followup_mg_value)))
      }
    }

    # Adjusted: follow-up Mg → AKI
    adj_fu <- intersect(adj_vars, names(trt_fu))
    fml_fu <- as.formula(paste("aki_kdigo1 ~ followup_mg_value +",
                                paste(adj_fu, collapse="+")))
    trt_fu_cc <- trt_fu[complete.cases(trt_fu[, c("followup_mg_value", adj_fu)]),]
    tryCatch({
      fit_fu <- glm(fml_fu, data=trt_fu_cc, family=binomial())
      ct_fu <- coeftest(fit_fu, vcov.=vcovHC(fit_fu, type="HC1"))
      fu_row <- which(rownames(ct_fu) == "followup_mg_value")
      or_fu <- exp(ct_fu[fu_row,1]); lo_fu <- exp(ct_fu[fu_row,1]-1.96*ct_fu[fu_row,2])
      hi_fu <- exp(ct_fu[fu_row,1]+1.96*ct_fu[fu_row,2])
      p_fu <- 2*pnorm(-abs(ct_fu[fu_row,1]/ct_fu[fu_row,2]))
      cat(sprintf("\n  Adjusted OR per 1 mg/dL follow-up Mg: %.3f (%.3f–%.3f) P=%.4f\n",
          or_fu, lo_fu, hi_fu, p_fu))
      cat(sprintf("  Interpretation: %s\n",
          ifelse(or_fu < 0.8 && p_fu < 0.10,
                 "Higher achieved Mg → LOWER AKI (supports threshold hypothesis)",
                 "No clear achieved-level effect within supplemented patients")))
    }, error=function(e) cat(sprintf("  FAILED: %s\n", e$message)))
  }
}

# ── Combined: dose + baseline Mg interaction ─────────────────────────
cat(sprintf("\n%s\n5. DOSE × BASELINE Mg INTERACTION\n%s\n",
            strrep("-",65), strrep("-",65)))
cat("  Does dose matter more when baseline Mg is already high?\n")

trt_cc$mg_high <- as.integer(trt_cc$first_mg_value >= 2.0)
fml_int <- as.formula(paste("aki_kdigo1 ~ dose_g * mg_high +",
                             paste(adj_vars[adj_vars != "first_mg_value"], collapse="+")))
tryCatch({
  fit_int <- glm(fml_int, data=trt_cc, family=binomial())
  ct_int <- coeftest(fit_int, vcov.=vcovHC(fit_int, type="HC1"))
  int_row <- grep("dose_g:mg_high", rownames(ct_int))
  if (length(int_row) > 0) {
    int_or <- exp(ct_int[int_row,1])
    int_p <- ct_int[int_row,4]
    cat(sprintf("  dose_g × baseline_Mg_high interaction: OR %.3f P=%.4f\n", int_or, int_p))
    cat(sprintf("  %s\n",
        ifelse(int_or < 0.7 && int_p < 0.10,
               "Dose effect STRONGER when baseline Mg is high → threshold model supported",
               "No significant interaction")))
  }
}, error=function(e) cat(sprintf("  FAILED: %s\n", e$message)))

# ── Also: ALL-PATIENT dose analysis (including untreated as dose=0) ──
cat(sprintf("\n%s\n6. ALL-PATIENT ANALYSIS: Dose as continuous (0 for untreated)\n%s\n",
            strrep("-",65), strrep("-",65)))

dat$dose_g_all <- ifelse(dat$mg_supp == 1, dat$mg_total_dose / 1000, 0)
dat_cc <- dat[complete.cases(dat[, c("dose_g_all", adj_vars, "aki_kdigo1")]),]
cat(sprintf("  N=%d (untreated=%d, supplemented=%d)\n",
    nrow(dat_cc), sum(dat_cc$dose_g_all==0), sum(dat_cc$dose_g_all>0)))

fml_all <- as.formula(paste("aki_kdigo1 ~ dose_g_all +", paste(adj_vars, collapse="+")))
tryCatch({
  fit_all <- glm(fml_all, data=dat_cc, family=binomial())
  ct_all <- coeftest(fit_all, vcov.=vcovHC(fit_all, type="HC1"))
  dose_row <- which(rownames(ct_all) == "dose_g_all")
  or_all <- exp(ct_all[dose_row,1]); lo_all <- exp(ct_all[dose_row,1]-1.96*ct_all[dose_row,2])
  hi_all <- exp(ct_all[dose_row,1]+1.96*ct_all[dose_row,2])
  p_all <- 2*pnorm(-abs(ct_all[dose_row,1]/ct_all[dose_row,2]))
  cat(sprintf("  OR per 1g MgSO4: %.3f (%.3f–%.3f) P=%.4f\n", or_all, lo_all, hi_all, p_all))
}, error=function(e) cat(sprintf("  FAILED: %s\n", e$message)))

# ── Save ─────────────────────────────────────────────────────────────
cat(sprintf("\n%s\nSUMMARY\n%s\n", strrep("=",65), strrep("=",65)))
cat("  Key question: within supplemented patients, does higher dose → lower AKI?\n")
cat("  If yes → causal dose-response supports pharmacological effect\n")
cat("  If no  → binary supplementation signal may reflect confounding\n")
cat("\n  Secondary: does dose → higher achieved Mg → lower AKI?\n")
cat("  (the full pharmacokinetic → pharmacodynamic chain)\n")
cat(sprintf("\n✓ 08c_dose_response.R COMPLETE\n"))
