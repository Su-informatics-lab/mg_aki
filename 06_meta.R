#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────
# 06_meta.R — Fixed-Effects Meta-Analysis: eICU + MIMIC-IV
# ─────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({ library(tidyverse) })

RESULTS <- path.expand("~/mg_aki/results")

cat(rep("=", 70), "\n", sep = "")
cat("META-ANALYSIS: eICU + MIMIC-IV\n")
cat(rep("=", 70), "\n", sep = "")

eicu_path  <- file.path(RESULTS, "03_results_summary.csv")
mimic_path <- file.path(RESULTS, "05_mimic_results_summary.csv")
if (!file.exists(eicu_path))  stop("Missing: ", eicu_path)
if (!file.exists(mimic_path)) stop("Missing: ", mimic_path)

eicu  <- read_csv(eicu_path,  show_col_types = FALSE)
mimic <- read_csv(mimic_path, show_col_types = FALSE)
cat(sprintf("  eICU:  %d estimates\n  MIMIC: %d estimates\n", nrow(eicu), nrow(mimic)))

get_est <- function(df, model_name) {
  r <- df %>% filter(model == model_name)
  if (nrow(r) == 0) return(NULL)
  list(or = r$estimate[1], lo = r$conf.low[1], hi = r$conf.high[1], n = r$n_events[1])
}

meta_or <- function(label, or1, lo1, hi1, n1, or2, lo2, hi2, n2) {
  logOR1 <- log(or1); se1 <- (log(hi1) - log(lo1)) / (2 * 1.96)
  logOR2 <- log(or2); se2 <- (log(hi2) - log(lo2)) / (2 * 1.96)
  w1 <- 1/se1^2; w2 <- 1/se2^2
  pooled_logOR <- (w1*logOR1 + w2*logOR2) / (w1 + w2)
  pooled_se <- sqrt(1 / (w1 + w2))
  pooled_OR <- exp(pooled_logOR)
  pooled_lo <- exp(pooled_logOR - 1.96 * pooled_se)
  pooled_hi <- exp(pooled_logOR + 1.96 * pooled_se)
  pooled_p <- 2 * pnorm(-abs(pooled_logOR / pooled_se))
  Q <- w1 * (logOR1 - pooled_logOR)^2 + w2 * (logOR2 - pooled_logOR)^2
  I2 <- max(0, (Q - 1) / Q * 100)

  cat(sprintf("\n  %-40s\n", label))
  cat(sprintf("    eICU:   OR %.2f (%.2f-%.2f), events=%d\n", or1, lo1, hi1, n1))
  cat(sprintf("    MIMIC:  OR %.2f (%.2f-%.2f), events=%d\n", or2, lo2, hi2, n2))
  cat(sprintf("    Pooled: OR %.2f (%.2f-%.2f), p=%.4f, I²=%.0f%%\n",
              pooled_OR, pooled_lo, pooled_hi, pooled_p, I2))

  data.frame(outcome = label,
             eicu_or = or1, eicu_lo = lo1, eicu_hi = hi1,
             mimic_or = or2, mimic_lo = lo2, mimic_hi = hi2,
             pooled_or = pooled_OR, pooled_lo = pooled_lo,
             pooled_hi = pooled_hi, pooled_p = pooled_p, I2 = I2)
}

# ── Outcome pairs: (label, eICU model, MIMIC model) ────────────────
pairs <- list(
  list("KDIGO Stage >=1",        "ALL_aki_KDIGO >=1",        "ALL_aki_KDIGO >=1"),
  list("Ratio >=1.5x",           "ALL_aki_Ratio >=1.5x",     "ALL_aki_Ratio >=1.5x"),
  list("Delta >=0.3",            "ALL_aki_Delta >=0.3",       "ALL_aki_Delta >=0.3"),
  list("Stage >=2",              "ALL_aki_Stage >=2",         "ALL_aki_Stage >=2"),
  list("Stage >=3",              "ALL_aki_Stage >=3",         "ALL_aki_Stage >=3"),
  list("AKI 1.5x <=48h",        "ALL_tw_AKI 1.5x <=48h",    "ALL_tw_AKI 1.5x <=48h"),
  list("AKI 1.5x <=72h",        "ALL_tw_AKI 1.5x <=72h",    "ALL_tw_AKI 1.5x <=72h"),
  list("Mortality (exploratory)","ALL_sec_Hospital mortality", "ALL_sec_Hospital mortality"),
  list("AKI KDIGO >=1 (AC)",     "AC_aki_KDIGO >=1",          "AC_aki_KDIGO >=1"),
  list("Fracture NC (AC)",       "AC_nc_Fracture",            "AC_nc_Fracture"),
  list("Fracture (neg ctrl)",    "ALL_nc_Fracture",           "ALL_nc_Fracture"),
  list("Encephalopathy",         "ALL_neuro_Encephalopathy",  "ALL_neuro_Encephalopathy")
)

results <- list()
skipped <- character(0)
for (p in pairs) {
  e <- get_est(eicu, p[[2]]); m <- get_est(mimic, p[[3]])
  if (is.null(e)) { skipped <- c(skipped, paste0(p[[1]], " (eICU: ", p[[2]], ")")); next }
  if (is.null(m)) { skipped <- c(skipped, paste0(p[[1]], " (MIMIC: ", p[[3]], ")")); next }
  results[[length(results)+1]] <- meta_or(p[[1]], e$or, e$lo, e$hi, e$n,
                                                  m$or, m$lo, m$hi, m$n)
}
if (length(skipped) > 0) {
  cat("\n  Skipped:\n")
  for (s in skipped) cat(sprintf("    - %s\n", s))
}

# ── Summary ─────────────────────────────────────────────────────────
cat("\n", rep("=", 70), "\n", sep = "")
cat("POOLED RESULTS\n")
cat(rep("=", 70), "\n", sep = "")

if (length(results) > 0) {
  res_df <- bind_rows(results)
  cat(sprintf("\n  %-40s %7s %7s %18s %8s %5s\n",
              "Outcome", "eICU", "MIMIC", "Pooled (95%CI)", "p", "I²"))
  cat("  ", strrep("-", 86), "\n")
  for (i in seq_len(nrow(res_df))) {
    r <- res_df[i, ]
    sig <- ifelse(r$pooled_p < 0.05, " *", "")
    cat(sprintf("  %-40s %7.2f %7.2f %6.2f (%.2f-%.2f) %8.4f %4.0f%%%s\n",
                r$outcome, r$eicu_or, r$mimic_or,
                r$pooled_or, r$pooled_lo, r$pooled_hi,
                r$pooled_p, r$I2, sig))
  }
  write_csv(res_df, file.path(RESULTS, "06_meta_results.csv"))
  cat(sprintf("\nSaved: %s (%d pooled estimates)\n",
              file.path(RESULTS, "06_meta_results.csv"), nrow(res_df)))
}
cat("\n06_meta.R COMPLETE\n")
