#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────
# 02_psm.R — Propensity score estimation + IPTW / matching
#
# Produces:
#   results/02a_ps_diagnostics.pdf
#   results/02b_analysis_a_prepared.RData
#   results/02c_hypo_matched.csv     — hypoMg 1:1 PS matched
#   results/02d_hypo_iptw.csv        — hypoMg IPTW weighted
#   results/02e_all_iptw.csv         — all-patient IPTW weighted
#   results/02f_all_matched.csv      — all-patient 1:1 PS matched
# ─────────────────────────────────────────────────────────────────────

local({
  pkgs <- c("tidyverse", "MatchIt", "cobalt", "WeightIt",
            "survey", "survival", "tableone")
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    install.packages(missing, repos = "https://cloud.r-project.org",
                     quiet = TRUE, Ncpus = parallel::detectCores())
  }
})

suppressPackageStartupMessages({
  library(tidyverse); library(MatchIt); library(cobalt)
  library(WeightIt); library(survey); library(survival); library(tableone)
})

RESULTS <- path.expand("~/mg_aki/results")

cat("Loading cohorts...\n")
cohort_a <- read_csv(file.path(RESULTS, "01_analysis_a_cohort.csv"), show_col_types = FALSE)
cohort_hypo <- tryCatch(
  read_csv(file.path(RESULTS, "01b_hypo_cohort.csv"), show_col_types = FALSE),
  error = function(e) { cat("  No hypoMg cohort found\n"); NULL })

cat(sprintf("  Prognostic: %d patients\n", nrow(cohort_a)))
if (!is.null(cohort_hypo))
  cat(sprintf("  TTE hypoMg: %d patients (%d treated, %d untreated)\n",
              nrow(cohort_hypo),
              sum(cohort_hypo$mg_supplementation == 1, na.rm = TRUE),
              sum(cohort_hypo$mg_supplementation == 0, na.rm = TRUE)))

# ─── Covariate formula ─────────────────────────────────────────────
cov_formula_rhs <- paste(c(
  "age_num", "is_female", "bmi", "surgery_type",
  "hx_chf", "hx_hypertension", "hx_diabetes", "hx_ckd",
  "hx_copd", "hx_pvd", "hx_stroke", "hx_liver",
  "baseline_cr", "baseline_egfr",
  "nephrotox_loop_diuretic", "nephrotox_nsaid",
  "nephrotox_acei_arb", "nephrotox_ppi",
  "has_betablocker", "has_steroid"
), collapse = " + ")

if ("apachescore" %in% names(cohort_a))
  cat("  NOTE: Using predicted ICU mortality instead of APACHE score (post-treatment variable)\n")

for (v in c("preop_antiarrhythmic", "first_k_value",
            "has_vasopressor", "first_hr", "first_ca_value"))
  if (v %in% names(cohort_a)) cov_formula_rhs <- paste(cov_formula_rhs, "+", v)

for (v in c("first_k_value", "first_hr", "first_ca_value")) {
  if (v %in% names(cohort_a) && any(is.na(cohort_a[[v]]))) {
    med_val <- median(cohort_a[[v]], na.rm = TRUE)
    n_imp <- sum(is.na(cohort_a[[v]]))
    cohort_a[[v]][is.na(cohort_a[[v]])] <- med_val
    cat(sprintf("  Imputed %s: %d NA -> median %.2f\n", v, n_imp, med_val))
  }
}

# =====================================================================
# PROGNOSTIC
# =====================================================================
cat("\n", strrep("=", 70), "\n")
cat("PROGNOSTIC: Mg level as continuous exposure\n")
cat(strrep("=", 70), "\n")

dat_a <- cohort_a %>%
  mutate(surgery_type = factor(surgery_type, levels = c("cabg","valve","combined","other_cardiac")),
         mg_quartile = factor(mg_quartile),
         ethnicity = factor(replace_na(ethnicity, "Other/Unknown"))) %>%
  filter(!is.na(baseline_cr), !is.na(age_num)) %>%
  mutate(bmi = ifelse(!is.na(bmi) & bmi > 10 & bmi < 80, bmi, NA_real_))

cat(sprintf("  Prognostic sample: %d\n", nrow(dat_a)))

tab1_vars <- intersect(c("age_num","is_female","bmi","surgery_type",
  "hx_chf","hx_hypertension","hx_diabetes","hx_ckd",
  "baseline_cr","baseline_egfr"), names(dat_a))
if (length(tab1_vars) > 0) {
  cat("\n  Covariate balance across Mg quartiles:\n")
  print(CreateTableOne(vars=tab1_vars, strata="mg_quartile", data=dat_a, test=FALSE), smd=TRUE)
}

save(dat_a, file = file.path(RESULTS, "02b_analysis_a_prepared.RData"))
cat("  Saved prognostic prepared data\n")

# =====================================================================
# TTE HYPOMG (Mg < 2.0)
# =====================================================================
if (!is.null(cohort_hypo) && nrow(cohort_hypo) > 10) {
  cat("\n", strrep("=", 70), "\n")
  cat("TTE HYPO: Mg supplementation vs. none (Mg < 2.0)\n")
  cat(strrep("=", 70), "\n")

  dat_hypo <- cohort_hypo %>%
    mutate(surgery_type = factor(surgery_type, levels = c("cabg","valve","combined","other_cardiac")),
           trt = as.integer(mg_supplementation)) %>%
    filter(!is.na(baseline_cr), !is.na(age_num))

  for (v in c("first_k_value","first_hr","first_ca_value"))
    if (v %in% names(dat_hypo) && any(is.na(dat_hypo[[v]])))
      dat_hypo[[v]][is.na(dat_hypo[[v]])] <- median(dat_hypo[[v]], na.rm=TRUE)

  ps_formula <- as.formula(paste("trt ~", cov_formula_rhs, "+ first_mg_value"))
  dat_hypo <- dat_hypo %>% drop_na(any_of(all.vars(ps_formula)))

  cat(sprintf("  TTE hypoMg sample: %d (treated: %d)\n", nrow(dat_hypo), sum(dat_hypo$trt)))

  cat("\n  Fitting propensity score model...\n")
  ps_model <- glm(ps_formula, data = dat_hypo, family = binomial())
  dat_hypo$ps <- predict(ps_model, type = "response")
  cat(sprintf("  PS: mean=%.3f, range=[%.3f, %.3f]\n",
              mean(dat_hypo$ps), min(dat_hypo$ps), max(dat_hypo$ps)))

  cat("\n  Computing IPTW (stabilized, truncated 1/99)...\n")
  w_obj <- weightit(ps_formula, data=dat_hypo, method="ps", estimand="ATE")
  dat_hypo$iptw_raw <- w_obj$weights
  q01 <- quantile(dat_hypo$iptw_raw, 0.01); q99 <- quantile(dat_hypo$iptw_raw, 0.99)
  dat_hypo$iptw <- pmax(pmin(dat_hypo$iptw_raw, q99), q01)
  cat(sprintf("  IPTW: median=%.2f, max=%.2f, ESS=%.0f\n",
              median(dat_hypo$iptw), max(dat_hypo$iptw),
              (sum(dat_hypo$iptw))^2/sum(dat_hypo$iptw^2)))
  cat("\n  Balance (IPTW):\n"); print(bal.tab(w_obj, stats=c("m","v"), thresholds=c(m=0.1)))

  cat("\n  PS Matching (1:1 nearest-neighbor, caliper 0.2 SD)...\n")
  tryCatch({
    m_out <- matchit(ps_formula, data=dat_hypo, method="nearest", distance="glm", caliper=0.2, ratio=1)
    dat_hypo_matched <- match.data(m_out)
    cat(sprintf("  Matched: %d pairs (%d total)\n", sum(dat_hypo_matched$trt==1), nrow(dat_hypo_matched)))
    cat("\n  Balance (matched):\n"); print(bal.tab(m_out, stats=c("m","v"), thresholds=c(m=0.1)))
    write_csv(dat_hypo_matched, file.path(RESULTS, "02c_hypo_matched.csv"))
  }, error = function(e) cat(sprintf("  Matching failed: %s\n", e$message)))

  cat("\n  Computing overlap weights (Li et al. 2018)...\n")
  dat_hypo$ow <- ifelse(dat_hypo$trt==1, 1-dat_hypo$ps, dat_hypo$ps)
  cat(sprintf("  Overlap weights: median=%.3f, ESS=%.0f\n",
              median(dat_hypo$ow), (sum(dat_hypo$ow))^2/sum(dat_hypo$ow^2)))

  write_csv(dat_hypo, file.path(RESULTS, "02d_hypo_iptw.csv"))
  cat("  Saved TTE hypoMg IPTW cohort\n")
} else {
  cat("\n  Skipping TTE hypoMg — insufficient sample\n")
}

# =====================================================================
# TTE ALL PATIENTS
# =====================================================================
cat("\n", strrep("=", 70), "\n")
cat("TTE ALL: Routine Mg supplementation (all patients)\n")
cat(strrep("=", 70), "\n")

dat_all <- cohort_a %>%
  mutate(surgery_type = factor(surgery_type, levels = c("cabg","valve","combined","other_cardiac")),
         trt = as.integer(mg_supplementation)) %>%
  filter(!is.na(baseline_cr), !is.na(age_num))

for (v in c("first_k_value","first_hr","first_ca_value"))
  if (v %in% names(dat_all) && any(is.na(dat_all[[v]])))
    dat_all[[v]][is.na(dat_all[[v]])] <- median(dat_all[[v]], na.rm=TRUE)

ps_formula_all <- as.formula(paste("trt ~", cov_formula_rhs, "+ first_mg_value"))
dat_all <- dat_all %>% drop_na(any_of(all.vars(ps_formula_all)))

cat(sprintf("  TTE all sample: %d (treated: %d, untreated: %d)\n",
            nrow(dat_all), sum(dat_all$trt), sum(dat_all$trt==0)))

if (nrow(dat_all) > 50 && sum(dat_all$trt) > 20) {
  cat("\n  Fitting propensity score model...\n")
  ps_model_all <- glm(ps_formula_all, data=dat_all, family=binomial())
  dat_all$ps <- predict(ps_model_all, type="response")
  cat(sprintf("  PS: mean=%.3f, range=[%.3f, %.3f]\n",
              mean(dat_all$ps), min(dat_all$ps), max(dat_all$ps)))

  cat("\n  Computing IPTW (stabilized, truncated 1/99)...\n")
  w_obj_all <- weightit(ps_formula_all, data=dat_all, method="ps", estimand="ATE")
  dat_all$iptw_raw <- w_obj_all$weights
  q01 <- quantile(dat_all$iptw_raw, 0.01); q99 <- quantile(dat_all$iptw_raw, 0.99)
  dat_all$iptw <- pmax(pmin(dat_all$iptw_raw, q99), q01)
  cat(sprintf("  IPTW: median=%.2f, max=%.2f, ESS=%.0f\n",
              median(dat_all$iptw), max(dat_all$iptw),
              (sum(dat_all$iptw))^2/sum(dat_all$iptw^2)))
  cat("\n  Balance (IPTW):\n"); print(bal.tab(w_obj_all, stats=c("m","v"), thresholds=c(m=0.1)))

  dat_all$ow <- ifelse(dat_all$trt==1, 1-dat_all$ps, dat_all$ps)

  cat("\n  PS Matching (1:1 nearest-neighbor, caliper 0.2 SD)...\n")
  tryCatch({
    m_out_all <- matchit(ps_formula_all, data=dat_all, method="nearest",
                         distance="glm", caliper=0.2, ratio=1)
    dat_all_matched <- match.data(m_out_all)
    cat(sprintf("  Matched: %d pairs (%d total)\n", sum(dat_all_matched$trt==1), nrow(dat_all_matched)))
    write_csv(dat_all_matched, file.path(RESULTS, "02f_all_matched.csv"))
  }, error = function(e) cat(sprintf("  Matching failed: %s\n", e$message)))

  write_csv(dat_all, file.path(RESULTS, "02e_all_iptw.csv"))
  cat("  Saved TTE all IPTW cohort\n")
} else {
  cat("  Skipping TTE all — insufficient sample\n")
}

# ─── PS diagnostics PDF ────────────────────────────────────────────
cat("\n  Generating PS diagnostics PDF...\n")
tryCatch({
  pdf(file.path(RESULTS, "02a_ps_diagnostics.pdf"), width=10, height=8)
  if (nrow(dat_a) > 0)
    print(ggplot(dat_a, aes(x=first_mg_value,
      fill=factor(ifelse(aki_primary==1|aki_delta03==1,1,0)))) +
      geom_density(alpha=0.5) + labs(title="Prognostic: Mg by AKI",
        x="First postop Mg (mg/dL)", fill="AKI (KDIGO >=1)") + theme_minimal())
  if (exists("dat_hypo") && "ps" %in% names(dat_hypo))
    print(ggplot(dat_hypo, aes(x=ps, fill=factor(trt))) + geom_density(alpha=0.5) +
      labs(title="TTE hypoMg: PS overlap", x="Propensity score", fill="Mg supp") + theme_minimal())
  if (exists("dat_all") && "ps" %in% names(dat_all))
    print(ggplot(dat_all, aes(x=ps, fill=factor(trt))) + geom_density(alpha=0.5) +
      labs(title="TTE all: PS overlap", x="Propensity score", fill="Mg supp") + theme_minimal())
  dev.off()
  cat("  Saved PS diagnostics PDF\n")
}, error = function(e) { cat(sprintf("  PDF failed: %s\n", e$message)); try(dev.off(), silent=TRUE) })

cat("\n", strrep("=", 70), "\n")
cat("02_psm.R COMPLETE\nNext: Rscript 03_models.R\n")
cat(strrep("=", 70), "\n")
