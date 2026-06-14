#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────
# 06_meta.R — Fixed-Effects Meta-Analysis: eICU + MIMIC-IV
# Reads results from both 03_results_summary.csv and
# 05_mimic_results_summary.csv (no more hardcoded values)
# ─────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({ library(tidyverse) })

RESULTS <- path.expand("~/mg_aki/results")

cat(rep("=", 70), "\n", sep = "")
cat("META-ANALYSIS: eICU + MIMIC-IV\n")
cat(rep("=", 70), "\n", sep = "")

# ── Load results from both databases ────────────────────────────────
eicu_path  <- file.path(RESULTS, "03_results_summary.csv")
mimic_path <- file.path(RESULTS, "05_mimic_results_summary.csv")

if (!file.exists(eicu_path))  stop("Missing: ", eicu_path, "\n  Run: Rscript 03_models.R")
if (!file.exists(mimic_path)) stop("Missing: ", mimic_path, "\n  Run: Rscript 05_mimic_tte.R")

eicu  <- read_csv(eicu_path,  show_col_types = FALSE)
mimic <- read_csv(mimic_path, show_col_types = FALSE)

cat(sprintf("  eICU:  %d estimates from %s\n", nrow(eicu), eicu_path))
cat(sprintf("  MIMIC: %d estimates from %s\n", nrow(mimic), mimic_path))

# ── Helper: extract OR/CI/n from results df by model name ──────────
get_est <- function(df, model_name) {
  r <- df %>% filter(model == model_name)
  if (nrow(r) == 0) return(NULL)
  list(or = r$estimate[1], lo = r$conf.low[1], hi = r$conf.high[1],
       n = r$n_events[1])
}

# ── Helper: inverse-variance fixed-effects meta on log-OR scale ────
meta_or <- function(label, or1, lo1, hi1, n1, or2, lo2, hi2, n2) {
  logOR1 <- log(or1); se1 <- (log(hi1) - log(lo1)) / (2 * 1.96)
  logOR2 <- log(or2); se2 <- (log(hi2) - log(lo2)) / (2 * 1.96)

  w1 <- 1/se1^2; w2 <- 1/se2^2
  pooled_logOR <- (w1*logOR1 + w2*logOR2) / (w1 + w2)
  pooled_se <- sqrt(1 / (w1 + w2))
  pooled_OR <- exp(pooled_logOR)
  pooled_lo <- exp(pooled_logOR - 1.96 * pooled_se)
  pooled_hi <- exp(pooled_logOR + 1.96 * pooled_se)
  pooled_z <- pooled_logOR / pooled_se
  pooled_p <- 2 * pnorm(-abs(pooled_z))

  Q <- w1 * (logOR1 - pooled_logOR)^2 + w2 * (logOR2 - pooled_logOR)^2
  I2 <- max(0, (Q - 1) / Q * 100)
  Q_p <- pchisq(Q, df = 1, lower.tail = FALSE)

  cat(sprintf("\n  %-30s\n", label))
  cat(sprintf("    eICU:   OR %.2f (%.2f-%.2f), events=%d\n", or1, lo1, hi1, n1))
  cat(sprintf("    MIMIC:  OR %.2f (%.2f-%.2f), events=%d\n", or2, lo2, hi2, n2))
  cat(sprintf("    Pooled: OR %.2f (%.2f-%.2f), p=%.4f\n", pooled_OR, pooled_lo, pooled_hi, pooled_p))
  cat(sprintf("    I²=%.0f%%, Q p=%.3f\n", I2, Q_p))

  data.frame(outcome = label,
             eicu_or = or1, eicu_lo = lo1, eicu_hi = hi1,
             mimic_or = or2, mimic_lo = lo2, mimic_hi = hi2,
             pooled_or = pooled_OR, pooled_lo = pooled_lo,
             pooled_hi = pooled_hi, pooled_p = pooled_p,
             I2 = I2, Q_p = Q_p)
}

# ── Define outcome pairs ───────────────────────────────────────────
# Each pair: (label, eICU model name, MIMIC model name)
# Both CSVs use the same model naming convention (TTEB_*, TTEA_*)
pairs <- list(
  # AKI severity-stratified (TTE-B)
  list("KDIGO Stage >=1",   "TTEB_aki_KDIGO >=1",   "TTEB_aki_KDIGO >=1"),
  list("Ratio >=1.5x",      "TTEB_aki_Ratio >=1.5x", "TTEB_aki_Ratio >=1.5x"),
  list("Delta >=0.3",       "TTEB_aki_Delta >=0.3",  "TTEB_aki_Delta >=0.3"),
  list("Stage >=2",         "TTEB_aki_Stage >=2",    "TTEB_aki_Stage >=2"),
  list("Stage >=3",         "TTEB_aki_Stage >=3",    "TTEB_aki_Stage >=3"),
  # Time-windowed
  list("AKI 1.5x <=48h",   "TTEB_tw_AKI 1.5x <=48h", "TTEB_tw_AKI 1.5x <=48h"),
  list("AKI 1.5x <=72h",   "TTEB_tw_AKI 1.5x <=72h", "TTEB_tw_AKI 1.5x <=72h"),
  # Mortality (TTE-B)
  list("Hospital mortality (TTE-B)", "TTEB_sec_Hospital mortality", "TTEB_sec_Hospital mortality"),
  # Mortality (eICU TTE-A + MIMIC TTE-B) — the individually-significant pool
  list("Hospital mortality (eICU TTE-A + MIMIC TTE-B)", "TTEA_sec_Hospital mortality", "TTEB_sec_Hospital mortality"),
  # Negative controls
  list("Fracture (negative control)", "TTEB_nc_Fracture", "TTEB_nc_Fracture"),
  # Neuro
  list("Encephalopathy",    "TTEB_neuro_Encephalopathy", "TTEB_neuro_Encephalopathy")
)

# ── Run meta-analyses ───────────────────────────────────────────────
results <- list()
skipped <- character(0)

for (p in pairs) {
  label  <- p[[1]]
  e_key  <- p[[2]]
  m_key  <- p[[3]]
  e <- get_est(eicu, e_key)
  m <- get_est(mimic, m_key)

  if (is.null(e)) { skipped <- c(skipped, paste0(label, " (missing eICU: ", e_key, ")")); next }
  if (is.null(m)) { skipped <- c(skipped, paste0(label, " (missing MIMIC: ", m_key, ")")); next }

  results[[length(results)+1]] <- meta_or(label, e$or, e$lo, e$hi, e$n,
                                                  m$or, m$lo, m$hi, m$n)
}

if (length(skipped) > 0) {
  cat("\n  Skipped (missing data):\n")
  for (s in skipped) cat(sprintf("    - %s\n", s))
}

# ── Summary table ───────────────────────────────────────────────────
cat("\n", rep("=", 70), "\n", sep = "")
cat("POOLED RESULTS SUMMARY\n")
cat(rep("=", 70), "\n", sep = "")

if (length(results) > 0) {
  res_df <- bind_rows(results)
  cat(sprintf("\n  %-36s %8s %8s %18s %8s %5s\n",
              "Outcome", "eICU OR", "MIMIC OR", "Pooled OR (95%CI)", "p", "I²"))
  cat("  ", strrep("-", 84), "\n")
  for (i in seq_len(nrow(res_df))) {
    r <- res_df[i, ]
    sig <- ifelse(r$pooled_p < 0.05, " *", "")
    cat(sprintf("  %-36s %8.2f %8.2f %6.2f (%.2f-%.2f) %8.4f %4.0f%%%s\n",
                r$outcome, r$eicu_or, r$mimic_or,
                r$pooled_or, r$pooled_lo, r$pooled_hi,
                r$pooled_p, r$I2, sig))
  }

  out <- file.path(RESULTS, "06_meta_results.csv")
  write_csv(res_df, out)
  cat(sprintf("\nSaved: %s (%d pooled estimates)\n", out, nrow(res_df)))
} else {
  cat("WARNING: No pairs could be pooled\n")
}

cat("\n06_meta.R COMPLETE\n")
