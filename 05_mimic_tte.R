#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────
# 05_mimic_tte.R — MIMIC-IV External Validation
# Replicates TTE-B (unrestricted) from eICU in MIMIC-IV
# Focus: Mg supplementation → AKI (KDIGO ≥1, 48h window, dose-response)
# ─────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyverse); library(survey); library(broom)
  library(WeightIt); library(cobalt); library(MatchIt)
})

RESULTS <- path.expand("~/mg_aki/results")
cat("Loading MIMIC-IV cohort...\n")
dat <- read_csv(file.path(RESULTS, "04_mimic_cohort.csv"), show_col_types = FALSE)
cat(sprintf("  N = %d, treated = %d (%.1f%%)\n",
            nrow(dat), sum(dat$mg_supplementation), mean(dat$mg_supplementation)*100))

# ── Prepare ─────────────────────────────────────────────────────────
dat <- dat %>%
  mutate(
    trt = as.integer(mg_supplementation),
    aki_kdigo1 = as.integer(aki_kdigo1),
    aki_primary = as.integer(aki_primary),
    aki_primary_48h = as.integer(aki_primary_48h),
    aki_stage2 = as.integer(aki_stage2),
    aki_stage3 = as.integer(aki_stage3),
  )

# Median-impute missing covariates
for (v in c("first_k_value", "first_ca_value", "first_mg_value")) {
  if (v %in% names(dat) && any(is.na(dat[[v]]))) {
    med <- median(dat[[v]], na.rm = TRUE)
    n_imp <- sum(is.na(dat[[v]]))
    dat[[v]][is.na(dat[[v]])] <- med
    cat(sprintf("  Imputed %s: %d NA -> %.2f\n", v, n_imp, med))
  }
}

# ── PS model (simplified — match eICU covariates where available) ──
# MIMIC has: age, sex, baseline_cr, first_mg, first_k, first_ca, mortality
# Missing vs eICU: surgery_type, comorbidities, BMI, nephrotoxins, vasopressors, HR
# This is a SIMPLIFIED validation — fewer covariates but independent database

ps_vars <- c("age", "is_female", "baseline_cr", "first_mg_value")
for (v in c("first_k_value", "first_ca_value")) {
  if (v %in% names(dat)) ps_vars <- c(ps_vars, v)
}

ps_formula <- as.formula(paste("trt ~", paste(ps_vars, collapse = " + ")))
cat(sprintf("\n  PS formula: %s\n", deparse(ps_formula)))

dat_clean <- dat %>% drop_na(any_of(c(ps_vars, "aki_kdigo1")))
cat(sprintf("  After dropping NA: %d (treated: %d)\n",
            nrow(dat_clean), sum(dat_clean$trt)))

# ── IPTW ────────────────────────────────────────────────────────────
cat("\n  Fitting PS model + IPTW...\n")
w <- weightit(ps_formula, data = dat_clean, method = "ps", estimand = "ATE")
dat_clean$iptw_raw <- w$weights
q01 <- quantile(dat_clean$iptw_raw, 0.01)
q99 <- quantile(dat_clean$iptw_raw, 0.99)
dat_clean$iptw <- pmax(pmin(dat_clean$iptw_raw, q99), q01)

cat(sprintf("  IPTW: median=%.2f, max=%.2f, ESS=%.0f\n",
            median(dat_clean$iptw), max(dat_clean$iptw),
            (sum(dat_clean$iptw))^2 / sum(dat_clean$iptw^2)))

# Balance
bal <- bal.tab(w, stats = c("m"), thresholds = c(m = 0.1))
cat("\n  Balance:\n")
print(bal)

# ── PS Matching ─────────────────────────────────────────────────────
cat("\n  PS Matching (1:1)...\n")
dat_matched <- NULL
tryCatch({
  m_out <- matchit(ps_formula, data = dat_clean, method = "nearest",
                   distance = "glm", caliper = 0.2, ratio = 1)
  dat_matched <- match.data(m_out)
  cat(sprintf("  Matched: %d pairs\n", sum(dat_matched$trt)))
}, error = function(e) cat(sprintf("  Matching failed: %s\n", e$message)))

# ── Helper ──────────────────────────────────────────────────────────
run_iptw <- function(d, outcome, label) {
  if (!outcome %in% names(d)) return(NULL)
  d <- d %>% filter(!is.na(.data[[outcome]]))
  nev <- sum(d[[outcome]])
  if (nev < 5) return(NULL)
  tryCatch({
    des <- svydesign(ids = ~1, weights = ~iptw, data = d)
    m <- svyglm(as.formula(paste(outcome, "~ trt")),
                design = des, family = quasibinomial())
    s <- tidy(m, conf.int = TRUE, exponentiate = TRUE) %>% filter(term == "trt")
    if (nrow(s) > 0) { s$n_events <- nev; s$label <- label }
    s
  }, error = function(e) NULL)
}

ptbl <- function(rows, title) {
  cat(sprintf("\n  %s\n", title))
  cat(sprintf("  %-40s %6s %18s %8s\n", "Outcome", "Events", "OR (95% CI)", "p"))
  cat("  ", strrep("-", 76), "\n")
  for (r in rows) if (!is.null(r) && nrow(r) > 0)
    cat(sprintf("  %-40s %6d %6.2f (%.2f-%.2f) %8.4f\n",
                r$label, r$n_events, r$estimate, r$conf.low, r$conf.high, r$p.value))
}

# =====================================================================
cat("\n", strrep("=", 70), "\n")
cat("MIMIC-IV EXTERNAL VALIDATION RESULTS\n")
cat(strrep("=", 70), "\n")

# ── Pipeline check: Serum Mg elevation ──────────────────────────────
cat("\n  Pipeline check: Serum Mg elevation\n")
# Use follow-up Mg (not implemented in ETL yet — use crude comparison)
trt_mg <- dat_clean %>% filter(trt == 1) %>% pull(first_mg_value)
ctrl_mg <- dat_clean %>% filter(trt == 0) %>% pull(first_mg_value)
cat(sprintf("    Treated baseline Mg: %.2f ± %.2f\n", mean(trt_mg), sd(trt_mg)))
cat(sprintf("    Untreated baseline Mg: %.2f ± %.2f\n", mean(ctrl_mg), sd(ctrl_mg)))

# ── AKI severity-stratified ─────────────────────────────────────────
aki_oc <- c(aki_kdigo1="KDIGO >=1", aki_primary="Ratio >=1.5x",
            aki_primary_48h="Ratio >=1.5x (48h)",
            aki_stage2="Stage >=2", aki_stage3="Stage >=3")
aki_r <- lapply(names(aki_oc), function(o) run_iptw(dat_clean, o, aki_oc[[o]]))
ptbl(aki_r, "AKI: IPTW (primary)")

# ── PS Matched ──────────────────────────────────────────────────────
if (!is.null(dat_matched)) {
  cat("\n  AKI: PS Matched\n")
  for (ov in c("aki_kdigo1", "aki_primary", "aki_primary_48h")) {
    if (ov %in% names(dat_matched)) {
      tryCatch({
        m <- glm(as.formula(paste(ov, "~ trt")), data = dat_matched,
                 family = binomial(), weights = weights)
        s <- tidy(m, conf.int = TRUE, exponentiate = TRUE) %>% filter(term == "trt")
        if (nrow(s) > 0)
          cat(sprintf("    %s: OR %.2f (%.2f-%.2f), p=%.4f\n",
                      aki_oc[[ov]], s$estimate, s$conf.low, s$conf.high, s$p.value))
      }, error = function(e) NULL)
    }
  }
}

# ── Hospital mortality ──────────────────────────────────────────────
r_mort <- run_iptw(dat_clean, "hosp_mortality", "Hospital mortality")
ptbl(list(r_mort), "Mortality")

# ── Dose-response (MIMIC unique!) ───────────────────────────────────
if ("mg_total_dose" %in% names(dat_clean)) {
  cat("\n  Dose-response (MIMIC unique — eICU has no dose data):\n")
  treated <- dat_clean %>% filter(trt == 1, mg_total_dose > 0)
  if (nrow(treated) > 50) {
    # Tertiles of dose
    treated$dose_q <- ntile(treated$mg_total_dose, 3)
    cat(sprintf("    Dose tertiles: Q1=%.0f, Q2=%.0f, Q3=%.0f (median units)\n",
                median(treated$mg_total_dose[treated$dose_q == 1]),
                median(treated$mg_total_dose[treated$dose_q == 2]),
                median(treated$mg_total_dose[treated$dose_q == 3])))

    # AKI rate by dose tertile
    for (q in 1:3) {
      sub <- treated %>% filter(dose_q == q)
      cat(sprintf("    Q%d (n=%d): AKI KDIGO1 = %.1f%%, ratio1.5x = %.1f%%\n",
                  q, nrow(sub),
                  mean(sub$aki_kdigo1)*100, mean(sub$aki_primary)*100))
    }

    # Trend test
    tryCatch({
      m <- glm(aki_kdigo1 ~ mg_total_dose, data = treated, family = binomial())
      s <- tidy(m, conf.int = TRUE, exponentiate = TRUE) %>%
        filter(term == "mg_total_dose")
      if (nrow(s) > 0)
        cat(sprintf("    Dose-response (continuous): OR %.4f per unit, p=%.4f\n",
                    s$estimate, s$p.value))
    }, error = function(e) NULL)
  }
}

# ── Comparison table ────────────────────────────────────────────────
cat("\n", strrep("=", 70), "\n")
cat("COMPARISON: eICU vs MIMIC-IV\n")
cat(strrep("=", 70), "\n")
cat("                              eICU TTE-B    MIMIC-IV\n")
cat(sprintf("  N                           %8d     %8d\n", 7924, nrow(dat_clean)))
cat(sprintf("  Treated                     %8d     %8d\n", 1104, sum(dat_clean$trt)))
cat(sprintf("  AKI KDIGO >=1 rate          %7.1f%%    %7.1f%%\n",
            22.2, mean(dat_clean$aki_kdigo1)*100))

# Print eICU results for comparison
cat("\n  Key findings comparison:\n")
cat("  Outcome              eICU OR (p)       MIMIC OR (p)\n")
cat("  ", strrep("-", 56), "\n")
for (r in aki_r) {
  if (!is.null(r) && nrow(r) > 0) {
    cat(sprintf("  %-22s                    %.2f (%.4f)\n",
                r$label, r$estimate, r$p.value))
  }
}

cat("\n05_mimic_tte.R COMPLETE\n")
