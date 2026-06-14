#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────
# 06_meta.R — Fixed-Effects Meta-Analysis: eICU + MIMIC-IV
# Pools treatment effects across two independent databases
# ─────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyverse)
})

RESULTS <- path.expand("~/mg_aki/results")

cat("=" , rep("=", 69), "\n", sep = "")
cat("META-ANALYSIS: eICU + MIMIC-IV\n")
cat(rep("=", 70), "\n", sep = "")

# ── eICU TTE-B results (from 03_results_summary.csv) ────────────────
eicu <- read_csv(file.path(RESULTS, "03_results_summary.csv"), show_col_types = FALSE)

# ── Hardcoded MIMIC results (from 05_mimic_tte.R output) ────────────
# Update these after running 05_mimic_tte.R with full model
# Format: outcome, OR, lower, upper, n_events, n_total

# These will be filled from the actual output — placeholder structure
cat("\nReading eICU results...\n")

# Helper: inverse-variance fixed-effects meta on log-OR scale
meta_or <- function(label, or1, lo1, hi1, n1, or2, lo2, hi2, n2) {
  logOR1 <- log(or1); se1 <- (log(hi1) - log(lo1)) / (2 * 1.96)
  logOR2 <- log(or2); se2 <- (log(hi2) - log(lo2)) / (2 * 1.96)

  # Fixed-effects inverse-variance
  w1 <- 1/se1^2; w2 <- 1/se2^2
  pooled_logOR <- (w1*logOR1 + w2*logOR2) / (w1 + w2)
  pooled_se <- sqrt(1 / (w1 + w2))
  pooled_OR <- exp(pooled_logOR)
  pooled_lo <- exp(pooled_logOR - 1.96 * pooled_se)
  pooled_hi <- exp(pooled_logOR + 1.96 * pooled_se)
  pooled_z <- pooled_logOR / pooled_se
  pooled_p <- 2 * pnorm(-abs(pooled_z))

  # Heterogeneity
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

# ── Extract eICU estimates ──────────────────────────────────────────
get_eicu <- function(model_name) {
  r <- eicu %>% filter(model == model_name)
  if (nrow(r) == 0) return(NULL)
  list(or = r$estimate[1], lo = r$conf.low[1], hi = r$conf.high[1],
       n = r$n_events[1])
}

# ── Run meta-analyses ───────────────────────────────────────────────
cat("\n", rep("-", 70), "\n", sep = "")
cat("AKI OUTCOMES (TTE-B IPTW)\n")
cat(rep("-", 70), "\n", sep = "")

results <- list()

# KDIGO >=1
e <- get_eicu("TTEB_aki_KDIGO >=1")
results[[1]] <- meta_or("KDIGO Stage >=1",
  e$or, e$lo, e$hi, e$n,
  0.91, 0.74, 1.13, 1227)

# Ratio >=1.5x
e <- get_eicu("TTEB_aki_Ratio >=1.5x")
results[[2]] <- meta_or("Ratio >=1.5x",
  e$or, e$lo, e$hi, e$n,
  0.90, 0.71, 1.13, 908)

# Delta >=0.3
e <- get_eicu("TTEB_aki_Delta >=0.3")
results[[3]] <- meta_or("Delta >=0.3",
  e$or, e$lo, e$hi, e$n,
  0.94, 0.75, 1.19, 943)

# Stage >=2
e <- get_eicu("TTEB_aki_Stage >=2")
results[[4]] <- meta_or("Stage >=2",
  e$or, e$lo, e$hi, e$n,
  0.87, 0.63, 1.20, 446)

# Stage >=3
e <- get_eicu("TTEB_aki_Stage >=3")
results[[5]] <- meta_or("Stage >=3",
  e$or, e$lo, e$hi, e$n,
  0.85, 0.50, 1.45, 159)

cat("\n", rep("-", 70), "\n", sep = "")
cat("TIME-WINDOWED AKI\n")
cat(rep("-", 70), "\n", sep = "")

# 48h
e <- get_eicu("TTEB_tw_AKI 1.5x <=48h")
results[[6]] <- meta_or("AKI 1.5x <=48h",
  e$or, e$lo, e$hi, e$n,
  0.95, 0.70, 1.27, 483)

# 72h
e <- get_eicu("TTEB_tw_AKI 1.5x <=72h")
results[[7]] <- meta_or("AKI 1.5x <=72h",
  e$or, e$lo, e$hi, e$n,
  0.95, 0.72, 1.25, 577)

cat("\n", rep("-", 70), "\n", sep = "")
cat("MORTALITY\n")
cat(rep("-", 70), "\n", sep = "")

# Hospital mortality — eICU TTE-B
e <- get_eicu("TTEB_sec_Hospital mortality")
if (!is.null(e)) {
  results[[8]] <- meta_or("Hospital mortality (TTE-B)",
    e$or, e$lo, e$hi, e$n,
    0.65, 0.45, 0.95, 353)
}

# Hospital mortality — eICU TTE-A vs MIMIC full
e_a <- eicu %>% filter(model == "TTEA_sec_Hospital mortality")
if (nrow(e_a) > 0) {
  results[[9]] <- meta_or("Hospital mortality (eICU TTE-A + MIMIC)",
    e_a$estimate[1], e_a$conf.low[1], e_a$conf.high[1], e_a$n_events[1],
    0.65, 0.45, 0.95, 353)
}

cat("\n", rep("-", 70), "\n", sep = "")
cat("NEGATIVE CONTROLS\n")
cat(rep("-", 70), "\n", sep = "")

e <- get_eicu("TTEB_nc_Fracture")
results[[10]] <- meta_or("Fracture (negative control)",
  e$or, e$lo, e$hi, e$n,
  1.01, 0.56, 1.80, 100)

cat("\n", rep("-", 70), "\n", sep = "")
cat("NEURO (exploratory)\n")
cat(rep("-", 70), "\n", sep = "")

e <- get_eicu("TTEB_neuro_Encephalopathy")
results[[11]] <- meta_or("Encephalopathy",
  e$or, e$lo, e$hi, e$n,
  0.47, 0.26, 0.86, 184)

# ── Summary table ───────────────────────────────────────────────────
cat("\n", rep("=", 70), "\n", sep = "")
cat("POOLED RESULTS SUMMARY\n")
cat(rep("=", 70), "\n", sep = "")

res_df <- bind_rows(results)
cat(sprintf("\n  %-30s %8s %8s %18s %8s %5s\n",
            "Outcome", "eICU OR", "MIMIC OR", "Pooled OR (95%CI)", "p", "I²"))
cat("  ", strrep("-", 80), "\n")
for (i in seq_len(nrow(res_df))) {
  r <- res_df[i, ]
  sig <- ifelse(r$pooled_p < 0.05, " *", "")
  cat(sprintf("  %-30s %8.2f %8.2f %6.2f (%.2f-%.2f) %8.4f %4.0f%%%s\n",
              r$outcome, r$eicu_or, r$mimic_or,
              r$pooled_or, r$pooled_lo, r$pooled_hi,
              r$pooled_p, r$I2, sig))
}

# Save
out <- file.path(RESULTS, "06_meta_results.csv")
write_csv(res_df, out)
cat(sprintf("\nSaved: %s\n", out))

cat("\n06_meta.R COMPLETE\n")
