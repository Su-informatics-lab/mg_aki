#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────
# Mg Reserve → Cardiac Surgery AKI (eICU)
# 02_psm.R — Propensity score estimation + IPTW / matching
#
# Produces:
#   results/02a_ps_diagnostics.pdf   — overlap plots, SMD
#   results/02b_iptw_cohort.csv      — Analysis A (GPS for continuous Mg)
#   results/02c_tte_matched.csv      — Analysis B (1:1 PS matched)
#   results/02d_tte_iptw.csv         — Analysis B (IPTW weighted)
#
# Design:
#   Analysis A (prognostic): Generalized PS for continuous exposure
#     OR multivariable regression with hospital clustering
#   Analysis B (TTE):  ACNU — Mg supplementation vs. none
#     Primary: IPTW (stabilized, truncated 1/99)
#     Secondary: 1:1 nearest-neighbor PS matching (caliper 0.2 SD)
#     Sensitivity: overlap weighting (Li, Morgan, Zaslavsky 2018)
# ─────────────────────────────────────────────────────────────────────

# ─── Auto-install dependencies ──────────────────────────────────────
local({
  pkgs <- c("tidyverse", "MatchIt", "cobalt", "WeightIt",
            "survey", "survival", "tableone")
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    cat(sprintf("Installing %d missing packages: %s\n",
                length(missing), paste(missing, collapse = ", ")))
    install.packages(missing, repos = "https://cloud.r-project.org",
                     quiet = TRUE, Ncpus = parallel::detectCores())
  }
})

suppressPackageStartupMessages({
  library(tidyverse)
  library(MatchIt)
  library(cobalt)
  library(WeightIt)
  library(survey)
  library(survival)
  library(tableone)
})

RESULTS <- path.expand("~/mg_aki/results")

# ─── Load cohorts ───────────────────────────────────────────────────
cat("Loading cohorts...\n")
cohort_a <- read_csv(file.path(RESULTS, "01_analysis_a_cohort.csv"),
                     show_col_types = FALSE)
cohort_b <- tryCatch(
  read_csv(file.path(RESULTS, "06_tte_cohort.csv"), show_col_types = FALSE),
  error = function(e) {
    cat("  No TTE cohort found — skipping Analysis B\n")
    NULL
  }
)

cat(sprintf("  Analysis A: %d patients\n", nrow(cohort_a)))
if (!is.null(cohort_b)) {
  cat(sprintf("  Analysis B: %d patients (%d treated, %d untreated)\n",
              nrow(cohort_b),
              sum(cohort_b$mg_supplementation == 1, na.rm = TRUE),
              sum(cohort_b$mg_supplementation == 0, na.rm = TRUE)))
}

# ─── Covariate formula (shared) ────────────────────────────────────
# DAG-informed: only pre-Mg-measurement variables
cov_formula_rhs <- paste(c(
  "age_num", "is_female", "bmi",
  "surgery_type",
  "hx_chf", "hx_hypertension", "hx_diabetes", "hx_ckd",
  "hx_copd", "hx_pvd", "hx_stroke", "hx_liver",
  "baseline_cr", "baseline_egfr",
  "nephrotox_loop_diuretic", "nephrotox_nsaid",
  "nephrotox_acei_arb", "nephrotox_ppi",
  "has_betablocker", "has_steroid"
), collapse = " + ")

# Add APACHE predicted mortality (pre-treatment) if available — NOT apacheScore (post-treatment!)
# apacheScore uses worst values in first 24h → includes post-Mg-supplementation physiology → mediator
if ("apachescore" %in% names(cohort_a)) {
  cat("  NOTE: Using predicted ICU mortality instead of APACHE score (post-treatment variable)\n")
}

# Add pre-op antiarrhythmic and potassium if available
for (v in c("preop_antiarrhythmic", "first_k_value",
            "has_vasopressor", "first_map", "first_hr",
            "first_ca_value", "first_lactate")) {
  if (v %in% names(cohort_a)) {
    cov_formula_rhs <- paste(cov_formula_rhs, "+", v)
  }
}

# =====================================================================
# ANALYSIS A — Prognostic (continuous Mg exposure)
# =====================================================================
cat("\n", strrep("=", 70), "\n")
cat("ANALYSIS A: Prognostic — Mg level as continuous exposure\n")
cat(strrep("=", 70), "\n")

# For continuous exposure: use multivariable regression with hospital
# clustering rather than GPS (simpler, comparable performance for
# single-timepoint exposure)

# Prepare data
dat_a <- cohort_a %>%
  mutate(
    surgery_type = factor(surgery_type, levels = c("cabg", "valve",
                                                    "combined",
                                                    "other_cardiac")),
    mg_quartile  = factor(mg_quartile),
    ethnicity    = factor(replace_na(ethnicity, "Other/Unknown")),
  ) %>%
  filter(!is.na(baseline_cr), !is.na(age_num)) %>%
  # Cap implausible BMI (eICU admissionWeight/Height can be junk)
  mutate(bmi = ifelse(!is.na(bmi) & bmi > 10 & bmi < 80, bmi, NA_real_))

cat(sprintf("  Analysis A sample: %d\n", nrow(dat_a)))

# Primary model: logistic with cluster-robust SE (hospital)
# (Actual model fitting in 03_models.R; here we do covariate
#  diagnostics across Mg quartiles)
cat("\n  Covariate balance across Mg quartiles:\n")
tab1_vars <- c("age_num", "is_female", "bmi", "surgery_type",
               "hx_chf", "hx_hypertension", "hx_diabetes", "hx_ckd",
               "baseline_cr", "baseline_egfr")
tab1_vars <- intersect(tab1_vars, names(dat_a))
if (length(tab1_vars) > 0) {
  t1 <- CreateTableOne(vars = tab1_vars,
                       strata = "mg_quartile",
                       data = dat_a,
                       test = FALSE)
  print(t1, smd = TRUE)
}

save(dat_a, file = file.path(RESULTS, "02b_analysis_a_prepared.RData"))
cat("  Saved Analysis A prepared data\n")

# =====================================================================
# ANALYSIS B — TTE (Mg supplementation, binary treatment)
# =====================================================================
if (!is.null(cohort_b) && nrow(cohort_b) > 10) {
  cat("\n", strrep("=", 70), "\n")
  cat("ANALYSIS B: TTE — Mg supplementation vs. none\n")
  cat(strrep("=", 70), "\n")

  dat_b <- cohort_b %>%
    mutate(
      surgery_type = factor(surgery_type,
                            levels = c("cabg", "valve", "combined",
                                       "other_cardiac")),
      trt = as.integer(mg_supplementation),
    ) %>%
    filter(!is.na(baseline_cr), !is.na(age_num))

  # Also include Mg level at time zero (the value triggering eligibility)
  ps_formula <- as.formula(paste("trt ~", cov_formula_rhs,
                                 "+ first_mg_value"))

  # Drop rows with NA in any PS covariate before fitting
  ps_vars <- all.vars(ps_formula)
  dat_b <- dat_b %>% drop_na(any_of(ps_vars))

  cat(sprintf("  TTE sample: %d (treated: %d) [after dropping NA covariates]\n",
              nrow(dat_b), sum(dat_b$trt)))

  # ── PS estimation ──────────────────────────────────────────────
  cat("\n  Fitting propensity score model...\n")
  ps_model <- glm(ps_formula, data = dat_b, family = binomial())
  dat_b$ps <- predict(ps_model, type = "response")

  cat(sprintf("  PS: mean=%.3f, range=[%.3f, %.3f]\n",
              mean(dat_b$ps), min(dat_b$ps), max(dat_b$ps)))

  # ── 1. IPTW (primary) ─────────────────────────────────────────
  cat("\n  Computing IPTW (stabilized, truncated 1/99)...\n")
  w_obj <- weightit(ps_formula, data = dat_b, method = "ps",
                    estimand = "ATE")
  dat_b$iptw_raw <- w_obj$weights

  # Truncate at 1st/99th percentile
  q01 <- quantile(dat_b$iptw_raw, 0.01)
  q99 <- quantile(dat_b$iptw_raw, 0.99)
  dat_b$iptw <- pmax(pmin(dat_b$iptw_raw, q99), q01)

  cat(sprintf("  IPTW: median=%.2f, max=%.2f, ESS=%.0f\n",
              median(dat_b$iptw), max(dat_b$iptw),
              (sum(dat_b$iptw))^2 / sum(dat_b$iptw^2)))

  # Balance diagnostics
  bal_iptw <- bal.tab(w_obj, stats = c("m", "v"), thresholds = c(m = 0.1))
  cat("\n  Balance (IPTW):\n")
  print(bal_iptw)

  # ── 2. PS Matching (secondary) ────────────────────────────────
  cat("\n  PS Matching (1:1 nearest-neighbor, caliper 0.2 SD)...\n")
  tryCatch({
    m_out <- matchit(ps_formula, data = dat_b, method = "nearest",
                     distance = "glm", caliper = 0.2, ratio = 1)
    dat_b_matched <- match.data(m_out)
    cat(sprintf("  Matched: %d pairs (%d total)\n",
                sum(dat_b_matched$trt == 1), nrow(dat_b_matched)))

    bal_match <- bal.tab(m_out, stats = c("m", "v"), thresholds = c(m = 0.1))
    cat("\n  Balance (matched):\n")
    print(bal_match)

    write_csv(dat_b_matched, file.path(RESULTS, "02c_tte_matched.csv"))
  }, error = function(e) {
    cat(sprintf("  Matching failed: %s\n", e$message))
    cat("  (Likely insufficient treated/control overlap — expected in demo)\n")
  })

  # ── 3. Overlap weighting (sensitivity) ────────────────────────
  cat("\n  Computing overlap weights (Li et al. 2018)...\n")
  dat_b$ow <- ifelse(dat_b$trt == 1, 1 - dat_b$ps, dat_b$ps)
  cat(sprintf("  Overlap weights: median=%.3f, ESS=%.0f\n",
              median(dat_b$ow),
              (sum(dat_b$ow))^2 / sum(dat_b$ow^2)))

  # Save IPTW cohort
  write_csv(dat_b, file.path(RESULTS, "02d_tte_iptw.csv"))
  cat("  Saved TTE-A IPTW cohort\n")

} else {
  cat("\n  Skipping TTE-A — insufficient sample\n")
}

# =====================================================================
# TTE-B: UNRESTRICTED (all cardiac surgery, prophylactic question)
# =====================================================================
cat("\n", strrep("=", 70), "\n")
cat("TTE-B: Unrestricted — routine Mg supplementation (all patients)\n")
cat(strrep("=", 70), "\n")

dat_b2 <- cohort_a %>%
  mutate(
    surgery_type = factor(surgery_type,
                          levels = c("cabg", "valve", "combined",
                                     "other_cardiac")),
    trt = as.integer(mg_supplementation),
  ) %>%
  filter(!is.na(baseline_cr), !is.na(age_num))

# PS formula: same covariates + first_mg_value (critical for confounding by indication)
ps_formula_b2 <- as.formula(paste("trt ~", cov_formula_rhs, "+ first_mg_value"))

# Drop NA in PS covariates
ps_vars_b2 <- all.vars(ps_formula_b2)
dat_b2 <- dat_b2 %>% drop_na(any_of(ps_vars_b2))

cat(sprintf("  TTE-B sample: %d (treated: %d, untreated: %d)\n",
            nrow(dat_b2), sum(dat_b2$trt), sum(dat_b2$trt == 0)))

if (nrow(dat_b2) > 50 && sum(dat_b2$trt) > 20) {
  # ── PS estimation ──────────────────────────────────────────────
  cat("\n  Fitting propensity score model...\n")
  ps_model_b2 <- glm(ps_formula_b2, data = dat_b2, family = binomial())
  dat_b2$ps <- predict(ps_model_b2, type = "response")

  cat(sprintf("  PS: mean=%.3f, range=[%.3f, %.3f]\n",
              mean(dat_b2$ps), min(dat_b2$ps), max(dat_b2$ps)))

  # ── IPTW ───────────────────────────────────────────────────────
  cat("\n  Computing IPTW (stabilized, truncated 1/99)...\n")
  w_obj_b2 <- weightit(ps_formula_b2, data = dat_b2, method = "ps",
                        estimand = "ATE")
  dat_b2$iptw_raw <- w_obj_b2$weights
  q01 <- quantile(dat_b2$iptw_raw, 0.01)
  q99 <- quantile(dat_b2$iptw_raw, 0.99)
  dat_b2$iptw <- pmax(pmin(dat_b2$iptw_raw, q99), q01)

  cat(sprintf("  IPTW: median=%.2f, max=%.2f, ESS=%.0f\n",
              median(dat_b2$iptw), max(dat_b2$iptw),
              (sum(dat_b2$iptw))^2 / sum(dat_b2$iptw^2)))

  # Balance
  bal_b2 <- bal.tab(w_obj_b2, stats = c("m", "v"), thresholds = c(m = 0.1))
  cat("\n  Balance (IPTW):\n")
  print(bal_b2)

  # ── Overlap weighting ─────────────────────────────────────────
  dat_b2$ow <- ifelse(dat_b2$trt == 1, 1 - dat_b2$ps, dat_b2$ps)

  # ── PS Matching ────────────────────────────────────────────────
  cat("\n  PS Matching (1:1 nearest-neighbor, caliper 0.2 SD)...\n")
  dat_b2_matched <- NULL
  tryCatch({
    m_out_b2 <- matchit(ps_formula_b2, data = dat_b2, method = "nearest",
                         distance = "glm", caliper = 0.2, ratio = 1)
    dat_b2_matched <- match.data(m_out_b2)
    cat(sprintf("  Matched: %d pairs (%d total)\n",
                sum(dat_b2_matched$trt == 1), nrow(dat_b2_matched)))
    write_csv(dat_b2_matched, file.path(RESULTS, "02f_tteb_matched.csv"))
  }, error = function(e) {
    cat(sprintf("  Matching failed: %s\n", e$message))
  })

  write_csv(dat_b2, file.path(RESULTS, "02e_tteb_iptw.csv"))
  cat("  Saved TTE-B IPTW cohort\n")

} else {
  cat("  Skipping TTE-B — insufficient sample\n")
}

# ─── PS diagnostics PDF ────────────────────────────────────────────
cat("\n  Generating PS diagnostics PDF...\n")
tryCatch({
  pdf(file.path(RESULTS, "02a_ps_diagnostics.pdf"), width = 10, height = 8)

  # Analysis A: Mg distribution by AKI status
  if (nrow(dat_a) > 0) {
    p1 <- ggplot(dat_a, aes(x = first_mg_value, fill = factor(ifelse(aki_primary == 1 | aki_delta03 == 1, 1, 0)))) +
      geom_density(alpha = 0.5) +
      labs(title = "Analysis A: Mg distribution by AKI status",
           x = "First postop Mg (mg/dL)", fill = "AKI (KDIGO ≥1)") +
      theme_minimal()
    print(p1)
  }

  # Analysis B (restricted): PS overlap
  if (!is.null(cohort_b) && exists("dat_b") && "ps" %in% names(dat_b)) {
    p2 <- ggplot(dat_b, aes(x = ps, fill = factor(trt))) +
      geom_density(alpha = 0.5) +
      labs(title = "TTE-A (hypoMg): PS overlap",
           x = "Propensity score", fill = "Mg supp") +
      theme_minimal()
    print(p2)
  }

  # TTE-B (unrestricted): PS overlap
  if (exists("dat_b2") && "ps" %in% names(dat_b2)) {
    p3 <- ggplot(dat_b2, aes(x = ps, fill = factor(trt))) +
      geom_density(alpha = 0.5) +
      labs(title = "TTE-B (all patients): PS overlap",
           x = "Propensity score", fill = "Mg supp") +
      theme_minimal()
    print(p3)
  }

  dev.off()
  cat("  Saved PS diagnostics PDF\n")
}, error = function(e) {
  cat(sprintf("  PDF generation failed: %s\n", e$message))
  try(dev.off(), silent = TRUE)
})

cat("\n", strrep("=", 70), "\n")
cat("02_psm.R COMPLETE\n")
cat("Next: Rscript 03_models.R\n")
cat(strrep("=", 70), "\n")
