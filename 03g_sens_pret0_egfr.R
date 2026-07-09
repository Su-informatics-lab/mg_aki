#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════
# 03g_sens_pret0_egfr.R — Sensitivity: pre-supplementation eGFR
#
# Clinical rationale (Yan): clinicians deciding whether to give Mg
# consider the CURRENT kidney function at the time of Mg, not the
# admission-baseline eGFR. This sensitivity re-stratifies the
# eGFR-specific HTE analysis using the LAST Cr before T₀ (ICU-only,
# no hospital-admission fallback) to compute a "pre-supplementation
# eGFR", then compares with the primary analysis (baseline eGFR).
#
# Reads:
#   did_hte_data_{db}.csv   — matched pair outcomes
#   did_labs_all_{db}.csv   — raw lab timings (offset_h)
#   did_all_{db}.csv        — mg_offset_h for T₀
#
# Outputs:
#   sens_pret0_egfr_{db}.csv — side-by-side comparison
#
# Usage:
#   Rscript 03g_sens_pret0_egfr.R mimic
#   Rscript 03g_sens_pret0_egfr.R eicu
# ══════════════════════════════════════════════════════════════════
suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
db   <- if (length(args) >= 1) tolower(args[1]) else "mimic"
cat(sprintf("\n══ 03g_sens_pret0_egfr.R [%s] ══\n", db))

RESULTS <- path.expand("~/mg_aki/results")

# ── CKD-EPI 2021 (race-free) ─────────────────────────────────
ckd_epi_2021 <- function(cr, age, is_female) {
  # cr in mg/dL, age in years, is_female = 0/1
  kappa <- ifelse(is_female == 1, 0.7, 0.9)
  alpha <- ifelse(is_female == 1, -0.241, -0.302)
  sex_mult <- ifelse(is_female == 1, 1.012, 1.0)
  142 * pmin(cr / kappa, 1)^alpha *
      pmax(cr / kappa, 1)^(-1.200) *
      0.9938^age * sex_mult
}

# ── Load data ─────────────────────────────────────────────────
hte <- fread(file.path(RESULTS, sprintf("did_hte_data_%s.csv", db)))
hte[, pid := as.character(pid)]

all_pts <- fread(file.path(RESULTS, sprintf("did_all_%s.csv", db)))
all_pts[, pid := as.character(pid)]

labs_file <- file.path(RESULTS, sprintf("did_labs_all_%s.csv", db))
labs <- fread(labs_file)
pid_col <- if ("patientunitstayid" %in% names(labs)) "patientunitstayid" else "stay_id"
labs[, pid := as.character(get(pid_col))]

# Get mg_offset_h (T₀) for each patient
mg_map <- all_pts[, .(pid, mg_offset_h)]

# ── Compute two eGFR values per patient ───────────────────────
cr_labs <- labs[lab_name == "creatinine" & value >= 0.1 & value <= 25.0]
cr_labs <- merge(cr_labs, mg_map, by = "pid", all.x = TRUE)

# 1. FIRST ICU Cr (earliest post-op measurement) → "baseline" eGFR
cr_first <- cr_labs[offset_h >= 0, .SD[which.min(offset_h)], by = pid]
setnames(cr_first, "value", "cr_first")

# 2. LAST Cr before T₀ (ICU only, no fallback) → "pre-supplementation" eGFR
#    Treated: offset_h < mg_offset_h AND offset_h >= 0 (ICU only)
#    Control: mg_offset_h is NA → all ICU Cr qualify → take last
cr_pre_t0 <- cr_labs[offset_h >= 0 & (is.na(mg_offset_h) | offset_h < mg_offset_h)]
cr_pre_t0 <- cr_pre_t0[, .SD[which.max(offset_h)], by = pid]
setnames(cr_pre_t0, "value", "cr_pre_t0")

# Merge into HTE data
hte <- merge(hte, cr_first[, .(pid, cr_first)], by = "pid", all.x = TRUE)
hte <- merge(hte, cr_pre_t0[, .(pid, cr_pre_t0)], by = "pid", all.x = TRUE)

# Compute eGFR from each
hte[, egfr_first  := ckd_epi_2021(cr_first,  age, is_female)]
hte[, egfr_pre_t0 := ckd_epi_2021(cr_pre_t0, age, is_female)]

cat(sprintf("  Patients with both eGFR values: %d / %d (%.1f%%)\n",
            sum(!is.na(hte$egfr_first) & !is.na(hte$egfr_pre_t0)),
            nrow(hte),
            100 * sum(!is.na(hte$egfr_first) & !is.na(hte$egfr_pre_t0)) / nrow(hte)))

# ── Compare eGFR distributions ────────────────────────────────
cat("\n  eGFR comparison (matched patients):\n")
cat(sprintf("    Baseline (first Cr):    median %.1f, IQR [%.1f-%.1f]\n",
            median(hte$egfr_first, na.rm = TRUE),
            quantile(hte$egfr_first, 0.25, na.rm = TRUE),
            quantile(hte$egfr_first, 0.75, na.rm = TRUE)))
cat(sprintf("    Pre-T₀ (last Cr < T₀): median %.1f, IQR [%.1f-%.1f]\n",
            median(hte$egfr_pre_t0, na.rm = TRUE),
            quantile(hte$egfr_pre_t0, 0.25, na.rm = TRUE),
            quantile(hte$egfr_pre_t0, 0.75, na.rm = TRUE)))
cat(sprintf("    Pearson r: %.3f\n",
            cor(hte$egfr_first, hte$egfr_pre_t0, use = "complete.obs")))

# Reclassification table
egfr_cut <- function(x) cut(x, breaks = c(0, 30, 45, 60, 90, Inf),
                            labels = c("<30", "30-44", "45-59", "60-89", ">=90"),
                            right = FALSE)
hte[, strat_first  := egfr_cut(egfr_first)]
hte[, strat_pre_t0 := egfr_cut(egfr_pre_t0)]

cat("\n  Reclassification matrix:\n")
reclass <- table(Baseline = hte$strat_first, `Pre-T0` = hte$strat_pre_t0, useNA = "ifany")
print(reclass)
pct_same <- sum(hte$strat_first == hte$strat_pre_t0, na.rm = TRUE) /
            sum(!is.na(hte$strat_first) & !is.na(hte$strat_pre_t0))
cat(sprintf("\n  Same stratum: %.1f%%\n", 100 * pct_same))

# ── Run eGFR-stratified OR for both definitions ──────────────
run_strat <- function(dat, egfr_col, outcome_col = "aki_48h") {
  dat$egfr_strat <- egfr_cut(dat[[egfr_col]])
  results <- list()

  for (s in c("Overall", levels(dat$egfr_strat))) {
    if (s == "Overall") {
      sub <- dat
    } else {
      sub <- dat[egfr_strat == s]
    }
    sub <- sub[!is.na(get(outcome_col))]
    n <- nrow(sub)
    n1 <- sum(sub$treated == 1)
    n0 <- sum(sub$treated == 0)
    if (n1 < 10 || n0 < 10) {
      results[[s]] <- data.table(stratum = s, n = n, or = NA_real_,
                                  or_lo = NA_real_, or_hi = NA_real_, p = NA_real_)
      next
    }
    fit <- tryCatch(
      glm(as.formula(paste(outcome_col, "~ treated")),
          data = sub, family = binomial),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      results[[s]] <- data.table(stratum = s, n = n, or = NA_real_,
                                  or_lo = NA_real_, or_hi = NA_real_, p = NA_real_)
      next
    }
    ci <- suppressMessages(confint(fit))
    or <- exp(coef(fit)["treated"])
    lo <- exp(ci["treated", 1])
    hi <- exp(ci["treated", 2])
    pv <- summary(fit)$coefficients["treated", 4]
    results[[s]] <- data.table(stratum = s, n = n, or = or, or_lo = lo,
                                or_hi = hi, p = pv)
  }
  rbindlist(results)
}

outcomes <- c("aki_48h", "aki_7d")
all_results <- list()

for (oc in outcomes) {
  if (!(oc %in% names(hte))) {
    cat(sprintf("  WARN: %s not in data, skipping\n", oc))
    next
  }

  r_base <- run_strat(hte, "egfr", oc)
  r_base[, egfr_def := "Primary (baseline eGFR)"]

  r_first <- run_strat(hte, "egfr_first", oc)
  r_first[, egfr_def := "Sens: first ICU Cr"]

  r_pret0 <- run_strat(hte, "egfr_pre_t0", oc)
  r_pret0[, egfr_def := "Sens: last Cr before T0"]

  combined <- rbind(r_base, r_first, r_pret0)
  combined[, outcome := oc]
  all_results[[oc]] <- combined

  cat(sprintf("\n  %s — eGFR-stratified ORs:\n", oc))
  cat(sprintf("  %-30s %6s %8s %17s %8s\n", "Definition / Stratum", "n", "OR", "95% CI", "P"))
  cat("  ", strrep("-", 72), "\n", sep = "")
  for (i in seq_len(nrow(combined))) {
    r <- combined[i]
    label <- sprintf("%s / %s", r$egfr_def, r$stratum)
    if (is.na(r$or)) {
      cat(sprintf("  %-30s %6d %8s\n", label, r$n, "—"))
    } else {
      ci <- sprintf("(%.2f-%.2f)", r$or_lo, r$or_hi)
      sig <- if (!is.na(r$p) && r$p < 0.05) "*" else ""
      cat(sprintf("  %-30s %6d %8.2f %17s %7.4f%s\n",
                  label, r$n, r$or, ci, r$p, sig))
    }
  }
}

# ── Save ──────────────────────────────────────────────────────
out <- rbindlist(all_results)
outfile <- file.path(RESULTS, sprintf("sens_pret0_egfr_%s.csv", db))
fwrite(out, outfile)
cat(sprintf("\n  Saved: %s\n", outfile))

# ── Concordance summary ───────────────────────────────────────
cat("\n  Concordance (direction of OR vs 1.0):\n")
for (oc in outcomes) {
  if (!(oc %in% out$outcome)) next
  base <- out[outcome == oc & egfr_def == "Primary (baseline eGFR)" & stratum != "Overall"]
  sens <- out[outcome == oc & egfr_def == "Sens: last Cr before T0" & stratum != "Overall"]
  m <- merge(base[, .(stratum, or_base = or)],
             sens[, .(stratum, or_sens = or)], by = "stratum")
  m <- m[!is.na(or_base) & !is.na(or_sens)]
  same_dir <- sum((m$or_base > 1) == (m$or_sens > 1))
  cat(sprintf("    %s: %d/%d strata concordant in direction\n",
              oc, same_dir, nrow(m)))
}

cat("\n══ DONE ══\n")
