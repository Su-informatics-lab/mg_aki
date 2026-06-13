#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────
# Mg Reserve → Cardiac Surgery AKI (eICU)
# 03_models.R — Main regression models + sensitivity analyses
#
# Analysis A (prognostic):
#   Primary: mixed-effects logistic (Mg continuous + quartiles),
#            hospital random intercept, cluster-robust SE
#   Secondary: Cox proportional hazards for time-to-AKI
#
# Analysis B (TTE — Mg supplementation):
#   Primary: IPTW-weighted logistic for 7-day AKI
#   Secondary: IPTW-weighted Cox, PS-matched logistic
#   Sensitivity: overlap-weighted, E-value
#
# Per TTE skill Tier 1: E-value for unmeasured confounding
# Per TTE skill Tier 2: negative control outcomes
# ─────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(survival)
  library(survey)
  library(sandwich)
  library(lmtest)
  library(EValue)
  library(broom)
  library(splines)
})

RESULTS <- Sys.getenv("RESULTS", unset = "~/mg_aki/results")

# ─── Load data ──────────────────────────────────────────────────────
cat("Loading prepared data...\n")
load(file.path(RESULTS, "02b_analysis_a_prepared.RData"))  # dat_a
dat_b <- tryCatch(
  read_csv(file.path(RESULTS, "02d_tte_iptw.csv"), show_col_types = FALSE),
  error = function(e) NULL
)
dat_b_m <- tryCatch(
  read_csv(file.path(RESULTS, "02c_tte_matched.csv"), show_col_types = FALSE),
  error = function(e) NULL
)

cov_rhs <- paste(c(
  "age_num", "is_female", "bmi",
  "surgery_type",
  "hx_chf", "hx_hypertension", "hx_diabetes", "hx_ckd",
  "hx_copd", "hx_pvd", "hx_stroke",
  "baseline_cr", "baseline_egfr",
  "nephrotox_loop_diuretic", "nephrotox_nsaid",
  "nephrotox_acei_arb"
), collapse = " + ")
if ("apacheScore" %in% names(dat_a)) cov_rhs <- paste(cov_rhs, "+ apacheScore")

results_list <- list()

# =====================================================================
# ANALYSIS A: Prognostic — Mg level → AKI
# =====================================================================
cat("\n", strrep("=", 70), "\n")
cat("ANALYSIS A: Mg level → AKI (prognostic)\n")
cat(strrep("=", 70), "\n")

# ── A1: Continuous Mg (per 1 mg/dL decrease) ─────────────────────
cat("\n  A1: Continuous Mg (per 1 mg/dL decrease)...\n")
tryCatch({
  # Negate Mg so positive coefficient = lower Mg → higher AKI
  dat_a$mg_neg <- -dat_a$first_mg_value

  if (length(unique(dat_a$hospitalID)) > 1) {
    # Mixed-effects with hospital random intercept
    f_a1 <- as.formula(paste("aki_primary ~ mg_neg +", cov_rhs,
                             "+ (1 | hospitalID)"))
    m_a1 <- glmer(f_a1, data = dat_a, family = binomial(),
                  control = glmerControl(optimizer = "bobyqa"))
  } else {
    f_a1 <- as.formula(paste("aki_primary ~ mg_neg +", cov_rhs))
    m_a1 <- glm(f_a1, data = dat_a, family = binomial())
  }

  s_a1 <- tidy(m_a1, conf.int = TRUE, exponentiate = TRUE)
  cat("  Key result (OR per 1 mg/dL DECREASE in Mg):\n")
  mg_row <- s_a1 %>% filter(term == "mg_neg")
  if (nrow(mg_row) > 0) {
    cat(sprintf("    OR = %.2f (95%% CI: %.2f–%.2f), p = %.4f\n",
                mg_row$estimate, mg_row$conf.low, mg_row$conf.high,
                mg_row$p.value))
    results_list[["A1_continuous"]] <- mg_row
  }
  cat("  Full model:\n")
  print(s_a1 %>% select(term, estimate, conf.low, conf.high, p.value), n = 30)
}, error = function(e) cat(sprintf("  A1 failed: %s\n", e$message)))

# ── A2: Mg quartiles (Q4=highest as reference) ──────────────────
cat("\n  A2: Mg quartiles (Q4 = reference)...\n")
tryCatch({
  dat_a$mg_q <- relevel(factor(dat_a$mg_quartile), ref = "Q4")

  f_a2 <- as.formula(paste("aki_primary ~ mg_q +", cov_rhs))
  m_a2 <- glm(f_a2, data = dat_a, family = binomial())
  # Cluster-robust SE by hospital
  if (length(unique(dat_a$hospitalID)) > 1) {
    m_a2_robust <- coeftest(m_a2, vcov = vcovCL(m_a2, cluster = dat_a$hospitalID))
    cat("  Quartile results (cluster-robust SE):\n")
    print(m_a2_robust[grep("mg_q", rownames(m_a2_robust)), ])
  }
  s_a2 <- tidy(m_a2, conf.int = TRUE, exponentiate = TRUE)
  q_rows <- s_a2 %>% filter(grepl("mg_q", term))
  cat("\n  Quartile ORs:\n")
  print(q_rows %>% select(term, estimate, conf.low, conf.high, p.value))
  results_list[["A2_quartiles"]] <- q_rows
}, error = function(e) cat(sprintf("  A2 failed: %s\n", e$message)))

# ── A3: Restricted cubic spline (dose-response) ─────────────────
cat("\n  A3: RCS dose-response...\n")
tryCatch({
  f_a3 <- as.formula(paste("aki_primary ~ ns(first_mg_value, df=4) +", cov_rhs))
  m_a3 <- glm(f_a3, data = dat_a, family = binomial())
  cat("  Spline model AIC:", AIC(m_a3), "\n")
  results_list[["A3_spline_aic"]] <- AIC(m_a3)
}, error = function(e) cat(sprintf("  A3 failed: %s\n", e$message)))

# ── A4: Cox for time-to-AKI ─────────────────────────────────────
cat("\n  A4: Cox proportional hazards (time-to-AKI)...\n")
tryCatch({
  dat_a_surv <- dat_a %>%
    filter(!is.na(time_to_event_hours), time_to_event_hours > 0)

  if (nrow(dat_a_surv) > 20) {
    f_a4 <- as.formula(paste(
      "Surv(time_to_event_hours, aki_primary) ~ mg_neg +", cov_rhs
    ))
    if (length(unique(dat_a_surv$hospitalID)) > 1) {
      m_a4 <- coxph(f_a4, data = dat_a_surv,
                     cluster = hospitalID)
    } else {
      m_a4 <- coxph(f_a4, data = dat_a_surv)
    }
    s_a4 <- tidy(m_a4, conf.int = TRUE, exponentiate = TRUE)
    mg_hr <- s_a4 %>% filter(term == "mg_neg")
    if (nrow(mg_hr) > 0) {
      cat(sprintf("    HR per 1 mg/dL decrease = %.2f (%.2f–%.2f), p = %.4f\n",
                  mg_hr$estimate, mg_hr$conf.low, mg_hr$conf.high,
                  mg_hr$p.value))
    }
    results_list[["A4_cox"]] <- mg_hr
  }
}, error = function(e) cat(sprintf("  A4 failed: %s\n", e$message)))

# ── A_sens: Eadon baseline sensitivity ──────────────────────────
cat("\n  A_sens: Eadon baseline (lowest 48h) sensitivity...\n")
tryCatch({
  if ("baseline_cr_eadon" %in% names(dat_a)) {
    dat_a_eadon <- dat_a %>%
      filter(!is.na(baseline_cr_eadon)) %>%
      mutate(baseline_cr_orig = baseline_cr,
             baseline_cr = baseline_cr_eadon)
    # Re-derive AKI with Eadon baseline (simplified — ratio only)
    dat_a_eadon <- dat_a_eadon %>%
      mutate(max_cr_ratio_eadon = max_followup_cr / baseline_cr,
             aki_eadon = as.integer(max_cr_ratio_eadon >= 1.5))

    f_eadon <- as.formula(paste("aki_eadon ~ mg_neg +", cov_rhs))
    m_eadon <- glm(f_eadon, data = dat_a_eadon, family = binomial())
    s_eadon <- tidy(m_eadon, conf.int = TRUE, exponentiate = TRUE)
    mg_eadon <- s_eadon %>% filter(term == "mg_neg")
    if (nrow(mg_eadon) > 0) {
      cat(sprintf("    OR (Eadon baseline) = %.2f (%.2f–%.2f)\n",
                  mg_eadon$estimate, mg_eadon$conf.low, mg_eadon$conf.high))
    }
    results_list[["A_sens_eadon"]] <- mg_eadon
  }
}, error = function(e) cat(sprintf("  Eadon sensitivity failed: %s\n", e$message)))

# =====================================================================
# ANALYSIS B: TTE — Mg supplementation → AKI
# =====================================================================
if (!is.null(dat_b) && nrow(dat_b) > 10) {
  cat("\n", strrep("=", 70), "\n")
  cat("ANALYSIS B: TTE — Mg supplementation → AKI\n")
  cat(strrep("=", 70), "\n")

  # ── B1: IPTW-weighted logistic (primary) ───────────────────────
  cat("\n  B1: IPTW-weighted logistic...\n")
  tryCatch({
    if ("hospitalID" %in% names(dat_b)) {
      des <- svydesign(ids = ~hospitalID, weights = ~iptw, data = dat_b)
    } else {
      des <- svydesign(ids = ~1, weights = ~iptw, data = dat_b)
    }
    m_b1 <- svyglm(aki_primary ~ trt, design = des,
                    family = quasibinomial())
    s_b1 <- tidy(m_b1, conf.int = TRUE, exponentiate = TRUE)
    trt_row <- s_b1 %>% filter(term == "trt")
    if (nrow(trt_row) > 0) {
      cat(sprintf("    IPTW OR = %.2f (95%% CI: %.2f–%.2f), p = %.4f\n",
                  trt_row$estimate, trt_row$conf.low, trt_row$conf.high,
                  trt_row$p.value))
      results_list[["B1_iptw"]] <- trt_row
    }
  }, error = function(e) cat(sprintf("  B1 failed: %s\n", e$message)))

  # ── B2: IPTW-weighted Cox ─────────────────────────────────────
  cat("\n  B2: IPTW-weighted Cox...\n")
  tryCatch({
    dat_b_surv <- dat_b %>%
      filter(!is.na(time_to_event_hours), time_to_event_hours > 0)
    if (nrow(dat_b_surv) > 10) {
      m_b2 <- coxph(Surv(time_to_event_hours, aki_primary) ~ trt,
                     weights = iptw, data = dat_b_surv,
                     robust = TRUE)
      s_b2 <- tidy(m_b2, conf.int = TRUE, exponentiate = TRUE)
      cat(sprintf("    IPTW HR = %.2f (%.2f–%.2f)\n",
                  s_b2$estimate[1], s_b2$conf.low[1], s_b2$conf.high[1]))
      results_list[["B2_iptw_cox"]] <- s_b2
    }
  }, error = function(e) cat(sprintf("  B2 failed: %s\n", e$message)))

  # ── B3: Overlap-weighted (sensitivity) ─────────────────────────
  cat("\n  B3: Overlap-weighted logistic...\n")
  tryCatch({
    des_ow <- svydesign(ids = ~1, weights = ~ow, data = dat_b)
    m_b3 <- svyglm(aki_primary ~ trt, design = des_ow,
                    family = quasibinomial())
    s_b3 <- tidy(m_b3, conf.int = TRUE, exponentiate = TRUE)
    trt_ow <- s_b3 %>% filter(term == "trt")
    if (nrow(trt_ow) > 0) {
      cat(sprintf("    OW OR = %.2f (%.2f–%.2f)\n",
                  trt_ow$estimate, trt_ow$conf.low, trt_ow$conf.high))
      results_list[["B3_overlap"]] <- trt_ow
    }
  }, error = function(e) cat(sprintf("  B3 failed: %s\n", e$message)))

  # ── B4: PS-matched logistic ────────────────────────────────────
  if (!is.null(dat_b_m) && nrow(dat_b_m) > 4) {
    cat("\n  B4: PS-matched logistic...\n")
    tryCatch({
      m_b4 <- glm(aki_primary ~ trt, data = dat_b_m,
                   family = binomial(), weights = weights)
      s_b4 <- tidy(m_b4, conf.int = TRUE, exponentiate = TRUE)
      trt_m <- s_b4 %>% filter(term == "trt")
      if (nrow(trt_m) > 0) {
        cat(sprintf("    Matched OR = %.2f (%.2f–%.2f)\n",
                    trt_m$estimate, trt_m$conf.low, trt_m$conf.high))
        results_list[["B4_matched"]] <- trt_m
      }
    }, error = function(e) cat(sprintf("  B4 failed: %s\n", e$message)))
  }

  # ── E-value (TTE Tier 1 — unmeasured confounding) ─────────────
  cat("\n  E-value (unmeasured confounding sensitivity)...\n")
  tryCatch({
    if (exists("trt_row") && nrow(trt_row) > 0) {
      ev <- evalues.OR(est = trt_row$estimate,
                       lo  = trt_row$conf.low,
                       hi  = trt_row$conf.high,
                       rare = (mean(dat_b$aki_primary) < 0.15))
      cat("  E-value results:\n")
      print(ev)
      results_list[["B_evalue"]] <- ev
    }
  }, error = function(e) cat(sprintf("  E-value failed: %s\n", e$message)))
}

# =====================================================================
# SAVE RESULTS SUMMARY
# =====================================================================
cat("\n", strrep("=", 70), "\n")
cat("RESULTS SUMMARY\n")
cat(strrep("=", 70), "\n")

summary_rows <- list()
for (nm in names(results_list)) {
  obj <- results_list[[nm]]
  if (is.data.frame(obj) && "estimate" %in% names(obj)) {
    for (i in seq_len(nrow(obj))) {
      summary_rows[[length(summary_rows) + 1]] <- tibble(
        model = nm,
        term = obj$term[i],
        estimate = obj$estimate[i],
        conf.low = ifelse("conf.low" %in% names(obj), obj$conf.low[i], NA),
        conf.high = ifelse("conf.high" %in% names(obj), obj$conf.high[i], NA),
        p.value = ifelse("p.value" %in% names(obj), obj$p.value[i], NA),
      )
    }
  }
}

if (length(summary_rows) > 0) {
  summary_df <- bind_rows(summary_rows)
  write_csv(summary_df, file.path(RESULTS, "03_results_summary.csv"))
  cat("\n")
  print(summary_df, n = 50)
  cat(sprintf("\n  Saved: %s\n", file.path(RESULTS, "03_results_summary.csv")))
} else {
  cat("  No results to summarize (likely demo data too small)\n")
}

cat("\n03_models.R COMPLETE\n")
