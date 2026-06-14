#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────
# 05_mimic_tte.R — MIMIC-IV External Validation (Full PS Model)
# Mirrors eICU 02_psm.R + 03_models.R in one script
#
# Fixes vs prior version:
#   - POAF removed (same-admission ICD can't distinguish pre/new-onset)
#   - Results saved to 05_mimic_results_summary.csv (matching eICU format)
#   - Hardcoded eICU comparison table removed (06_meta.R does this)
# ─────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyverse); library(survey); library(broom)
  library(WeightIt); library(cobalt); library(MatchIt)
})

RESULTS <- path.expand("~/mg_aki/results")
dat <- read_csv(file.path(RESULTS, "04_mimic_cohort.csv"), show_col_types = FALSE)
cat(sprintf("MIMIC-IV cohort: N=%d, treated=%d (%.1f%%)\n",
            nrow(dat), sum(dat$mg_supplementation), mean(dat$mg_supplementation)*100))

dat$trt <- as.integer(dat$mg_supplementation)
dat$surgery_type <- factor(dat$surgery_type,
                           levels = c("cabg","valve","combined","other_cardiac"))

results_list <- list()

# ── PS formula (matching eICU's covariates) ─────────────────────────
ps_vars <- c(
  "age", "is_female",
  "surgery_type",
  "hx_chf", "hx_hypertension", "hx_diabetes", "hx_ckd",
  "hx_copd", "hx_pvd", "hx_stroke", "hx_liver",
  "baseline_cr", "baseline_egfr",
  "nephrotox_loop_diuretic", "nephrotox_nsaid",
  "nephrotox_acei_arb", "nephrotox_ppi",
  "has_betablocker", "has_steroid", "preop_antiarrhythmic",
  "first_k_value", "first_ca_value",
  "has_vasopressor",
  "first_mg_value"
)
if ("bmi" %in% names(dat)) ps_vars <- c(ps_vars, "bmi")
if ("first_hr" %in% names(dat)) ps_vars <- c(ps_vars, "first_hr")
ps_vars <- ps_vars[ps_vars %in% names(dat)]
cat(sprintf("PS covariates: %d matched\n", length(ps_vars)))

# Median-impute continuous covariates
for (v in c("first_k_value", "first_ca_value", "first_hr",
            "bmi", "baseline_egfr")) {
  if (v %in% names(dat) && any(is.na(dat[[v]]))) {
    med <- median(dat[[v]], na.rm = TRUE)
    n_imp <- sum(is.na(dat[[v]]))
    dat[[v]][is.na(dat[[v]])] <- med
    cat(sprintf("  Imputed %s: %d NA -> %.1f\n", v, n_imp, med))
  }
}

ps_formula <- as.formula(paste("trt ~", paste(ps_vars, collapse = " + ")))
dat_clean <- dat %>% drop_na(any_of(c(ps_vars, "aki_kdigo1")))
cat(sprintf("After NA drop: %d (treated: %d)\n", nrow(dat_clean), sum(dat_clean$trt)))

# ── IPTW ────────────────────────────────────────────────────────────
cat("\nFitting PS + IPTW...\n")
w <- weightit(ps_formula, data = dat_clean, method = "ps", estimand = "ATE")
dat_clean$iptw_raw <- w$weights
q01 <- quantile(dat_clean$iptw_raw, 0.01)
q99 <- quantile(dat_clean$iptw_raw, 0.99)
dat_clean$iptw <- pmax(pmin(dat_clean$iptw_raw, q99), q01)
dat_clean$ow <- ifelse(dat_clean$trt == 1,
                       1 - predict(glm(ps_formula, data=dat_clean, family=binomial()), type="response"),
                       predict(glm(ps_formula, data=dat_clean, family=binomial()), type="response"))

cat(sprintf("IPTW: median=%.2f, max=%.2f, ESS=%.0f\n",
            median(dat_clean$iptw), max(dat_clean$iptw),
            (sum(dat_clean$iptw))^2 / sum(dat_clean$iptw^2)))

bal <- bal.tab(w, stats = "m", thresholds = c(m = 0.1))
cat("\nBalance:\n")
print(bal)

# ── PS Matching ─────────────────────────────────────────────────────
cat("\nPS Matching (1:1)...\n")
dat_matched <- NULL
tryCatch({
  m_out <- matchit(ps_formula, data = dat_clean, method = "nearest",
                   distance = "glm", caliper = 0.2, ratio = 1)
  dat_matched <- match.data(m_out)
  cat(sprintf("Matched: %d pairs\n", sum(dat_matched$trt)))
  bal_m <- bal.tab(m_out, stats = "m", thresholds = c(m = 0.1))
  cat("\nMatched balance:\n")
  print(bal_m)
}, error = function(e) cat(sprintf("Matching failed: %s\n", e$message)))

# ── Helpers ─────────────────────────────────────────────────────────
iptw_or <- function(d, outcome, label, wt = "iptw") {
  if (!outcome %in% names(d)) return(NULL)
  d <- d %>% filter(!is.na(.data[[outcome]]))
  nev <- sum(d[[outcome]])
  if (nev < 5) return(NULL)
  tryCatch({
    des <- svydesign(ids = ~1, weights = as.formula(paste0("~",wt)), data = d)
    m <- svyglm(as.formula(paste(outcome, "~ trt")), design = des,
                family = quasibinomial())
    s <- tidy(m, conf.int = TRUE, exponentiate = TRUE) %>% filter(term == "trt")
    if (nrow(s) > 0) { s$n_events <- nev; s$label <- label; s$n_total <- nrow(d) }
    s
  }, error = function(e) NULL)
}

ptbl <- function(rows, title) {
  cat(sprintf("\n  %s\n", title))
  cat(sprintf("  %-42s %6s %18s %8s\n", "Outcome", "Events", "OR (95% CI)", "p"))
  cat("  ", strrep("-", 76), "\n")
  for (r in rows) if (!is.null(r) && nrow(r) > 0)
    cat(sprintf("  %-42s %6d %6.2f (%.2f-%.2f) %8.4f\n",
                r$label, r$n_events, r$estimate, r$conf.low, r$conf.high, r$p.value))
}

store <- function(key, r) {
  if (!is.null(r) && nrow(r) > 0) results_list[[key]] <<- r
}

# =====================================================================
cat("\n", strrep("=", 70), "\n")
cat("MIMIC-IV TTE-B RESULTS (Full PS Model)\n")
cat(strrep("=", 70), "\n")

# ── Pipeline check ──────────────────────────────────────────────────
cat("\n  Pipeline check: Serum Mg elevation\n")
if ("delta_mg" %in% names(dat_clean)) {
  d_mg <- dat_clean %>% filter(!is.na(delta_mg))
  if (nrow(d_mg) > 20) {
    tryCatch({
      des <- svydesign(ids=~1, weights=~iptw, data=d_mg)
      m <- svyglm(delta_mg ~ trt, design=des)
      s <- tidy(m, conf.int=TRUE) %>% filter(term=="trt")
      cat(sprintf("    IPTW diff = %.3f (%.3f to %.3f), p=%.4f\n",
                  s$estimate, s$conf.low, s$conf.high, s$p.value))
      results_list[["TTEB_mg_elevation"]] <- tibble(
        model="TTEB_mg_elevation", term="trt", estimate=s$estimate,
        conf.low=s$conf.low, conf.high=s$conf.high, p.value=s$p.value,
        n_events=nrow(d_mg), label="Serum Mg elevation (pipeline)")
    }, error = function(e) cat("    Failed\n"))
  }
}

# ── VT/VF (boundary analysis) ──────────────────────────────────────
r_vt <- iptw_or(dat_clean, "vent_arrhythmia", "VT/VF (positive control)")
store("TTEB_vent_arrhythmia", r_vt)
if (!is.null(r_vt))
  cat(sprintf("\n  VT/VF: OR %.2f (%.2f-%.2f), p=%.4f  [n=%d]\n",
              r_vt$estimate, r_vt$conf.low, r_vt$conf.high, r_vt$p.value, r_vt$n_events))

# ── AKI severity-stratified ─────────────────────────────────────────
aki_oc <- c(aki_delta03="Delta >=0.3", aki_kdigo1="KDIGO >=1",
            aki_primary="Ratio >=1.5x", aki_stage2="Stage >=2",
            aki_stage3="Stage >=3")
aki_r <- lapply(names(aki_oc), function(o) iptw_or(dat_clean, o, aki_oc[[o]]))
ptbl(aki_r, "AKI: IPTW")
for (r in aki_r) store(paste0("TTEB_aki_", r$label), r)

# Time-windowed
tw_oc <- c(aki_primary_24h="AKI 1.5x <=24h", aki_primary_48h="AKI 1.5x <=48h",
           aki_primary_72h="AKI 1.5x <=72h")
tw_r <- lapply(names(tw_oc), function(o) iptw_or(dat_clean, o, tw_oc[[o]]))
ptbl(tw_r, "AKI: Time-windowed")
for (r in tw_r) store(paste0("TTEB_tw_", r$label), r)

# Peak Cr ratio
if ("peak_cr_ratio" %in% names(dat_clean)) {
  d_pcr <- dat_clean %>% filter(!is.na(peak_cr_ratio))
  tryCatch({
    des <- svydesign(ids=~1, weights=~iptw, data=d_pcr)
    m <- svyglm(peak_cr_ratio ~ trt, design=des)
    s <- tidy(m, conf.int=TRUE) %>% filter(term=="trt")
    cat(sprintf("\n  Peak Cr ratio diff = %.3f (%.3f to %.3f), p=%.4f\n",
                s$estimate, s$conf.low, s$conf.high, s$p.value))
    results_list[["TTEB_peak_cr"]] <- tibble(
      model="TTEB_peak_cr", term="trt", estimate=s$estimate,
      conf.low=s$conf.low, conf.high=s$conf.high, p.value=s$p.value,
      n_events=nrow(d_pcr), label="Peak Cr ratio diff")
  }, error = function(e) NULL)
}

# ── PS Matched results ──────────────────────────────────────────────
if (!is.null(dat_matched)) {
  cat("\n  AKI: PS Matched\n")
  for (ov in names(aki_oc)) {
    if (ov %in% names(dat_matched)) tryCatch({
      m <- glm(as.formula(paste(ov, "~ trt")), data=dat_matched,
               family=binomial(), weights=weights)
      s <- tidy(m, conf.int=TRUE, exponentiate=TRUE) %>% filter(term=="trt")
      if (nrow(s) > 0) {
        cat(sprintf("    %s: OR %.2f (%.2f-%.2f), p=%.4f\n",
                    aki_oc[[ov]], s$estimate, s$conf.low, s$conf.high, s$p.value))
        s$n_events <- sum(dat_matched[[ov]], na.rm=TRUE)
        s$label <- aki_oc[[ov]]
        s$n_total <- nrow(dat_matched)
        results_list[[paste0("TTEB_matched_", aki_oc[[ov]])]] <- s
      }
    }, error = function(e) NULL)
  }
}

# ── OW ──────────────────────────────────────────────────────────────
r_ow <- iptw_or(dat_clean, "aki_primary", "OW ratio >=1.5x", wt="ow")
store("TTEB_ow", r_ow)
if (!is.null(r_ow))
  cat(sprintf("\n  OW: OR %.2f (%.2f-%.2f), p=%.4f\n",
              r_ow$estimate, r_ow$conf.low, r_ow$conf.high, r_ow$p.value))

# ── Mortality ───────────────────────────────────────────────────────
r_mort <- iptw_or(dat_clean, "hosp_mortality", "Hospital mortality")
store("TTEB_sec_Hospital mortality", r_mort)
ptbl(list(r_mort), "Mortality")

# ── Negative controls ───────────────────────────────────────────────
nc_r <- list(
  iptw_or(dat_clean, "nc_fracture", "Fracture"),
  iptw_or(dat_clean, "nc_uti", "UTI")
)
ptbl(nc_r, "Negative controls")
for (r in nc_r) store(paste0("TTEB_nc_", r$label), r)

# ── Neuro ───────────────────────────────────────────────────────────
neuro_r <- list(
  iptw_or(dat_clean, "neuro_delirium", "Delirium"),
  iptw_or(dat_clean, "neuro_encephalopathy", "Encephalopathy")
)
ptbl(neuro_r, "Neuro (exploratory)")
for (r in neuro_r) store(paste0("TTEB_neuro_", r$label), r)

# ── Dose-response (MIMIC unique) ────────────────────────────────────
if ("mg_total_dose" %in% names(dat_clean)) {
  cat("\n  Dose-response (MIMIC unique):\n")
  treated <- dat_clean %>% filter(trt == 1, mg_total_dose > 0)
  if (nrow(treated) > 50) {
    treated$dose_q <- ntile(treated$mg_total_dose, 3)
    for (q in 1:3) {
      sub <- treated %>% filter(dose_q == q)
      cat(sprintf("    Q%d (n=%d, median dose=%.0f): AKI KDIGO1=%.1f%%, 1.5x=%.1f%%\n",
                  q, nrow(sub), median(sub$mg_total_dose),
                  mean(sub$aki_kdigo1)*100, mean(sub$aki_primary)*100))
    }
    tryCatch({
      m <- glm(aki_kdigo1 ~ mg_total_dose, data=treated, family=binomial())
      s <- tidy(m, conf.int=TRUE, exponentiate=TRUE) %>% filter(term=="mg_total_dose")
      if (nrow(s)>0) cat(sprintf("    Trend: OR %.4f per unit, p=%.4f\n", s$estimate, s$p.value))
    }, error = function(e) NULL)
  }
}


# =====================================================================
cat("\n", strrep("=", 70), "\n")
cat("TTE-A: HYPOMAGNESEMIA ONLY (Mg < 2.0)\n")
cat(strrep("=", 70), "\n")

dat_a <- dat_clean %>% filter(first_mg_value < 2.0)
cat(sprintf("TTE-A: N=%d (treated: %d)\n", nrow(dat_a), sum(dat_a$trt)))

if (sum(dat_a$trt) >= 30) {
  tryCatch({
    w_a <- weightit(ps_formula, data = dat_a, method = "ps", estimand = "ATE")
    dat_a$iptw_raw <- w_a$weights
    q01a <- quantile(dat_a$iptw_raw, 0.01); q99a <- quantile(dat_a$iptw_raw, 0.99)
    dat_a$iptw <- pmax(pmin(dat_a$iptw_raw, q99a), q01a)

    # AKI
    aki_a <- lapply(names(aki_oc), function(o) iptw_or(dat_a, o, aki_oc[[o]]))
    ptbl(aki_a, "TTE-A AKI")
    for (r in aki_a) store(paste0("TTEA_aki_", r$label), r)

    # Time-windowed
    tw_a <- lapply(names(tw_oc), function(o) iptw_or(dat_a, o, tw_oc[[o]]))
    ptbl(tw_a, "TTE-A Time-windowed AKI")
    for (r in tw_a) store(paste0("TTEA_tw_", r$label), r)

    # Mortality
    r_mort_a <- iptw_or(dat_a, "hosp_mortality", "Hospital mortality")
    store("TTEA_sec_Hospital mortality", r_mort_a)
    ptbl(list(r_mort_a), "TTE-A Mortality")

    # Negative controls
    nc_a <- list(
      iptw_or(dat_a, "nc_fracture", "Fracture"),
      iptw_or(dat_a, "nc_uti", "UTI")
    )
    ptbl(nc_a, "TTE-A Negative controls")
    for (r in nc_a) store(paste0("TTEA_nc_", r$label), r)

    # PS matching TTE-A
    tryCatch({
      m_a <- matchit(ps_formula, data = dat_a, method = "nearest",
                     distance = "glm", caliper = 0.2, ratio = 1)
      dat_ma <- match.data(m_a)
      cat(sprintf("\n  TTE-A matched: %d pairs\n", sum(dat_ma$trt)))
      for (ov in c("aki_kdigo1", "aki_primary", "aki_stage2", "hosp_mortality")) {
        if (ov %in% names(dat_ma)) tryCatch({
          m <- glm(as.formula(paste(ov, "~ trt")), data = dat_ma,
                   family = binomial(), weights = weights)
          s <- tidy(m, conf.int=TRUE, exponentiate=TRUE) %>% filter(term=="trt")
          lbl <- if (ov=="hosp_mortality") "Hospital mortality" else aki_oc[[ov]]
          if (nrow(s)>0) {
            cat(sprintf("    %s: OR %.2f (%.2f-%.2f), p=%.4f\n",
                        lbl, s$estimate, s$conf.low, s$conf.high, s$p.value))
            s$n_events <- sum(dat_ma[[ov]], na.rm=TRUE)
            s$label <- lbl; s$n_total <- nrow(dat_ma)
            results_list[[paste0("TTEA_matched_", lbl)]] <- s
          }
        }, error = function(e) NULL)
      }
    }, error = function(e) cat(sprintf("  TTE-A matching failed: %s\n", e$message)))

  }, error = function(e) cat(sprintf("TTE-A failed: %s\n", e$message)))
}


# =====================================================================
cat("\n", strrep("=", 70), "\n")
cat("SURGERY-TYPE INTERACTION (Cardioplegia Hypothesis)\n")
cat(strrep("=", 70), "\n")

tryCatch({
  dat_clean$surg_complex <- ifelse(dat_clean$surgery_type %in% c("combined","valve"),
                                    "complex", "simple")
  cov_rhs_mimic <- paste(setdiff(ps_vars, "first_mg_value"), collapse = " + ")
  m_int <- glm(as.formula(paste("aki_kdigo1 ~ first_mg_value * surg_complex +", cov_rhs_mimic)),
               data = dat_clean, family = binomial())
  s_int <- tidy(m_int, conf.int = TRUE, exponentiate = TRUE)
  main <- s_int %>% filter(term == "first_mg_value")
  inter <- s_int %>% filter(grepl(":", term))
  if (nrow(main) > 0)
    cat(sprintf("  Mg effect (simple): OR %.2f (%.2f-%.2f), p=%.4f\n",
                main$estimate, main$conf.low, main$conf.high, main$p.value))
  if (nrow(inter) > 0)
    cat(sprintf("  Interaction (complex×Mg): OR %.2f (%.2f-%.2f), p=%.4f\n",
                inter$estimate, inter$conf.low, inter$conf.high, inter$p.value))
  for (st in c("simple", "complex")) {
    dsub <- dat_clean %>% filter(surg_complex == st)
    msub <- glm(as.formula(paste("aki_kdigo1 ~ first_mg_value +", cov_rhs_mimic)),
                data = dsub, family = binomial())
    ssub <- tidy(msub, conf.int=TRUE, exponentiate=TRUE) %>% filter(term=="first_mg_value")
    if (nrow(ssub) > 0)
      cat(sprintf("  Stratified %s (n=%d): OR %.2f (%.2f-%.2f), p=%.4f\n",
                  st, nrow(dsub), ssub$estimate, ssub$conf.low, ssub$conf.high, ssub$p.value))
  }
}, error = function(e) cat(sprintf("  Interaction failed: %s\n", e$message)))


# =====================================================================
cat("\n", strrep("=", 70), "\n")
cat("E-VALUES\n")
cat(strrep("=", 70), "\n")

tryCatch({
  if (!require(EValue, quietly=TRUE)) install.packages("EValue", lib="~/R/libs", repos="https://cloud.r-project.org")
  library(EValue)

  evals <- list(
    c("Mortality (IPTW)", 0.65, 0.45, 0.95),
    c("AKI KDIGO>=1 (IPTW)", 0.91, 0.74, 1.13),
    c("AKI Ratio 1.5x (IPTW)", 0.90, 0.71, 1.13),
    c("AKI Stage>=3 (matched)", 0.69, 0.40, 1.18),
    c("Encephalopathy (IPTW)", 0.47, 0.26, 0.86)
  )
  for (ev in evals) {
    tryCatch({
      or <- as.numeric(ev[2]); lo <- as.numeric(ev[3]); hi <- as.numeric(ev[4])
      res <- evalues.OR(or, lo, hi, rare = TRUE)
      cat(sprintf("  %s: OR=%.2f, E-value=%.2f (CI bound=%.2f)\n",
                  ev[1], or, res["E-values","point"], res["E-values","lower"]))
    }, error = function(e) cat(sprintf("  %s: failed\n", ev[1])))
  }
}, error = function(e) cat(sprintf("E-value computation failed: %s\n", e$message)))


# =====================================================================
# SAVE RESULTS (matching 03_results_summary.csv format)
# =====================================================================
cat("\n", strrep("=", 70), "\n")
cat("SAVING MIMIC RESULTS\n")
cat(strrep("=", 70), "\n")

rows <- list()
for (nm in names(results_list)) {
  obj <- results_list[[nm]]
  if (is.data.frame(obj) && "estimate" %in% names(obj))
    for (i in seq_len(nrow(obj)))
      rows[[length(rows)+1]] <- tibble(
        model = nm,
        term = if ("label" %in% names(obj)) obj$label[i] else obj$term[i],
        estimate = obj$estimate[i],
        conf.low = if ("conf.low" %in% names(obj)) obj$conf.low[i] else NA,
        conf.high = if ("conf.high" %in% names(obj)) obj$conf.high[i] else NA,
        p.value = if ("p.value" %in% names(obj)) obj$p.value[i] else NA,
        n_events = if ("n_events" %in% names(obj)) obj$n_events[i] else NA)
}
if (length(rows) > 0) {
  df <- bind_rows(rows)
  out <- file.path(RESULTS, "05_mimic_results_summary.csv")
  write_csv(df, out)
  print(df, n = 100)
  cat(sprintf("\nSaved: %s (%d estimates)\n", out, nrow(df)))
} else {
  cat("WARNING: No results to save\n")
}

cat("\n05_mimic_tte.R COMPLETE\n")
