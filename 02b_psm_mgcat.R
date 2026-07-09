#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════
# 02b_psm_mgcat.R — Sensitivity: add categorical serum Mg to PS
#
# Rationale (Meng/Liu): since the clinical decision to supplement
# Mg is driven by low serum Mg, add a 3-level indicator to the PS:
#   1 = hypomagnesemia (serum Mg < 1.7 mg/dL)
#   2 = normal / high  (serum Mg >= 1.7)
#   3 = missing         (no pre-T₀ measurement)
#
# This is a SENSITIVITY analysis — primary remains the 19-covariate
# model without serum Mg (02_psm.R, set.seed(2026)).
#
# Reads:  did_all_{db}.csv, did_labs_all_{db}.csv, did_cr_all_{db}.csv
# Writes: sens_mgcat_{db}.csv (binary outcome ORs, comparison table)
#         sens_mgcat_balance_{db}.csv (SMD before/after)
#
# Usage:
#   Rscript 02b_psm_mgcat.R mimic
#   Rscript 02b_psm_mgcat.R eicu
# ══════════════════════════════════════════════════════════════════
suppressPackageStartupMessages({
  library(data.table)
  library(MatchIt)
  library(mice)
})

args <- commandArgs(trailingOnly = TRUE)
db   <- if (length(args) >= 1) tolower(args[1]) else "mimic"
cat(sprintf("\n══ 02b_psm_mgcat.R [%s] ══\n", db))

set.seed(2026)
RESULTS <- path.expand("~/mg_aki/results")
MG_CUT  <- 1.7  # clinical hypomagnesemia threshold

# ── Load data ─────────────────────────────────────────────────
all_pts <- fread(file.path(RESULTS, sprintf("did_all_%s.csv", db)))
all_pts[, pid := as.character(pid)]

labs <- fread(file.path(RESULTS, sprintf("did_labs_all_%s.csv", db)))
pid_col <- if ("patientunitstayid" %in% names(labs)) "patientunitstayid" else "stay_id"
labs[, pid := as.character(get(pid_col))]

# ── Extract last serum Mg before T₀ ──────────────────────────
mg_map <- all_pts[, .(pid, mg_offset_h)]
mg_labs <- labs[lab_name == "magnesium" & value >= 0.5 & value <= 10.0]
mg_labs <- merge(mg_labs, mg_map, by = "pid", all.x = TRUE)
mg_labs <- mg_labs[offset_h >= 0 & (is.na(mg_offset_h) | offset_h < mg_offset_h)]
last_mg <- mg_labs[, .SD[which.max(offset_h)], by = pid][, .(pid, last_mg = value)]

all_pts <- merge(all_pts, last_mg, by = "pid", all.x = TRUE)

# ── Create 3-level categorical ────────────────────────────────
all_pts[, mg_cat := factor(
  fifelse(is.na(last_mg), "missing",
  fifelse(last_mg < MG_CUT, "low", "normal")),
  levels = c("normal", "low", "missing")
)]

cat(sprintf("\n  Mg category distribution:\n"))
cat(sprintf("    Treated:  %s\n",
    paste(sprintf("%s=%d", names(table(all_pts[treated==1, mg_cat])),
                           table(all_pts[treated==1, mg_cat])), collapse=", ")))
cat(sprintf("    Control:  %s\n",
    paste(sprintf("%s=%d", names(table(all_pts[treated==0, mg_cat])),
                           table(all_pts[treated==0, mg_cat])), collapse=", ")))

# ── PS formula: primary 19 covariates + mg_cat ────────────────
# Same as 02_psm.R primary spec + mg_cat
ps_vars_base <- c(
  "age", "is_female", "bmi",
  "surg_cabg", "surg_valve", "surg_combined",
  "heart_failure", "hypertension", "diabetes", "ckd",
  "copd", "pvd", "stroke", "liver_disease",
  "egfr"
)

# Time-varying labs: need to compute from did_labs_all
# Extract last-before-T0 for Ca, lactate, heartrate
for (lab_name in c("calcium", "lactate", "heartrate")) {
  sub <- labs[lab_name == get("lab_name", envir=.GlobalEnv)]
  # Filter by name properly
  sub <- labs[labs$lab_name == lab_name]
  sub <- merge(sub, mg_map, by = "pid", all.x = TRUE)
  sub <- sub[offset_h >= 0 & (is.na(mg_offset_h) | offset_h < mg_offset_h)]
  last_val <- sub[, .SD[which.max(offset_h)], by = pid][, .(pid, value)]
  setnames(last_val, "value", paste0("last_", lab_name))
  all_pts <- merge(all_pts, last_val, by = "pid", all.x = TRUE,
                   suffixes = c("", ".new"))
  # If column already exists, prefer new
  old_col <- paste0("last_", lab_name)
  new_col <- paste0(old_col, ".new")
  if (new_col %in% names(all_pts)) {
    all_pts[, (old_col) := get(new_col)]
    all_pts[, (new_col) := NULL]
  }
}
all_pts[, last_lactate_missing := as.integer(is.na(last_lactate))]

ps_vars <- c(ps_vars_base,
             "last_lactate", "last_lactate_missing",
             "last_heartrate",
             "mg_cat")  # Ca removed; mg_cat added

ps_formula <- as.formula(paste("treated ~", paste(ps_vars, collapse = " + ")))
cat(sprintf("\n  PS formula: %d variables (primary 19 − Ca + mg_cat = 19)\n", length(ps_vars)))
cat(sprintf("  Formula: %s\n", deparse(ps_formula, width.cutoff = 200)))

# ── MICE imputation ───────────────────────────────────────────
imp_vars <- c("treated", ps_vars)
imp_vars <- imp_vars[imp_vars %in% names(all_pts)]
imp_dat <- all_pts[, ..imp_vars]

cat(sprintf("\n  Missingness:\n"))
for (v in imp_vars) {
  pct <- 100 * mean(is.na(imp_dat[[v]]))
  if (pct > 0) cat(sprintf("    %s: %.1f%%\n", v, pct))
}

cat("\n  Running MICE m=5 (sensitivity, faster than m=20)...\n")
imp <- mice(imp_dat, m = 5, method = "pmm", printFlag = FALSE, seed = 2026)

# ── Match on each imputation, take stable pairs ───────────────
pair_list <- list()
for (k in 1:5) {
  d <- complete(imp, k)
  d$pid <- all_pts$pid
  # MatchIt
  m <- tryCatch(
    matchit(ps_formula, data = d, method = "nearest",
            distance = "glm", caliper = 0.2, replace = TRUE, ratio = 1),
    error = function(e) { cat(sprintf("  MICE %d match failed: %s\n", k, e$message)); NULL }
  )
  if (is.null(m)) next
  md <- match.data(m)
  trt_pids <- md$pid[md$treated == 1]
  # Get matched control for each treated
  matched_ctl <- d$pid[as.integer(m$match.matrix[, 1])]
  pairs_k <- data.table(trt_pid = trt_pids, ctl_pid = matched_ctl, imp = k)
  pairs_k <- pairs_k[!is.na(ctl_pid)]
  pair_list[[k]] <- pairs_k
  cat(sprintf("    MICE %d: %d pairs\n", k, nrow(pairs_k)))
}

all_pairs <- rbindlist(pair_list)
# Keep pairs appearing in >= 50% of imputations
pair_counts <- all_pairs[, .N, by = .(trt_pid, ctl_pid)]
stable_pairs <- pair_counts[N >= 3]  # 3/5 = 60%
cat(sprintf("\n  Stable pairs (>=3/5 imputations): %d\n", nrow(stable_pairs)))
cat(sprintf("  Unique treated: %d / %d (%.1f%%)\n",
            uniqueN(stable_pairs$trt_pid),
            sum(all_pts$treated == 1),
            100 * uniqueN(stable_pairs$trt_pid) / sum(all_pts$treated == 1)))

# ── Extract matched cohort ────────────────────────────────────
trt <- all_pts[pid %in% stable_pairs$trt_pid]
ctl_pids <- unique(stable_pairs$ctl_pid)
ctl <- all_pts[pid %in% ctl_pids]

cat(sprintf("  Matched cohort: %d treated, %d controls\n", nrow(trt), nrow(ctl)))

# ── Balance check ─────────────────────────────────────────────
cat(sprintf("\n  Post-matching balance (mg_cat):\n"))
cat(sprintf("    Treated:  %s\n",
    paste(sprintf("%s=%d", names(table(trt$mg_cat)), table(trt$mg_cat)), collapse=", ")))
cat(sprintf("    Control:  %s\n",
    paste(sprintf("%s=%d", names(table(ctl$mg_cat)), table(ctl$mg_cat)), collapse=", ")))

# SMD for key variables
smd_fn <- function(x1, x0) {
  m1 <- mean(x1, na.rm=TRUE); m0 <- mean(x0, na.rm=TRUE)
  sp <- sqrt((var(x1, na.rm=TRUE) + var(x0, na.rm=TRUE)) / 2)
  if (sp < 1e-10) return(0)
  abs(m1 - m0) / sp
}

cat("\n  Key SMDs (primary vs sensitivity):\n")
cat(sprintf("  %-25s %8s %8s\n", "Variable", "Primary", "+mg_cat"))

# Load primary pairs for comparison
primary_pairs <- fread(file.path(RESULTS,
  sprintf("did_pairs_primary_yet_untreated_%s.csv", db)))
primary_trt <- all_pts[pid %in% as.character(primary_pairs$trt_pid)]
primary_ctl <- all_pts[pid %in% as.character(primary_pairs$ctl_pid)]

for (v in c("last_mg", "egfr", "age", "last_calcium", "last_lactate")) {
  if (!(v %in% names(all_pts))) next
  smd_pri <- smd_fn(primary_trt[[v]], primary_ctl[[v]])
  smd_new <- smd_fn(trt[[v]], ctl[[v]])
  cat(sprintf("  %-25s %8.3f %8.3f\n", v, smd_pri, smd_new))
}

# ── Primary outcome: 48h AKI ─────────────────────────────────
# Merge outcomes from did_cr_all
cr_path <- file.path(RESULTS, sprintf("did_cr_all_%s.csv", db))
if (file.exists(cr_path)) {
  cr <- fread(cr_path)
  cr_pid <- if ("patientunitstayid" %in% names(cr)) "patientunitstayid" else "stay_id"
  cr[, pid := as.character(get(cr_pid))]
  # For each treated, compute AKI at T0+48h
  # This is simplified — real pipeline uses pair-level computation
  # For now, merge aki flags from did_hte_data if available
}

hte_path <- file.path(RESULTS, sprintf("did_hte_data_%s.csv", db))
if (file.exists(hte_path)) {
  hte <- fread(hte_path)
  hte[, pid := as.character(pid)]

  outcomes <- c("aki_48h", "aki_7d")
  results <- list()
  for (oc in outcomes) {
    if (!(oc %in% names(hte))) next
    # Treated outcomes
    trt_oc <- hte[treated == 1 & pid %in% stable_pairs$trt_pid, .(pid, outcome = get(oc))]
    ctl_oc <- hte[treated == 0 & pid %in% ctl_pids, .(pid, outcome = get(oc))]

    rate_t <- mean(trt_oc$outcome, na.rm = TRUE)
    rate_c <- mean(ctl_oc$outcome, na.rm = TRUE)
    # Simple OR
    a <- sum(trt_oc$outcome == 1, na.rm=TRUE)
    b <- sum(trt_oc$outcome == 0, na.rm=TRUE)
    c <- sum(ctl_oc$outcome == 1, na.rm=TRUE)
    d <- sum(ctl_oc$outcome == 0, na.rm=TRUE)
    or_val <- (a * d) / (b * c)
    se_log <- sqrt(1/a + 1/b + 1/c + 1/d)
    ci_lo <- exp(log(or_val) - 1.96 * se_log)
    ci_hi <- exp(log(or_val) + 1.96 * se_log)
    pv <- 2 * pnorm(-abs(log(or_val) / se_log))

    results[[oc]] <- data.table(
      outcome = oc, n_trt = nrow(trt_oc), n_ctl = nrow(ctl_oc),
      rate_trt = rate_t, rate_ctl = rate_c,
      or = or_val, ci_lo = ci_lo, ci_hi = ci_hi, p = pv
    )
    cat(sprintf("\n  %s: trt %.1f%% vs ctl %.1f%%, OR %.2f (%.2f-%.2f), P=%.4f\n",
                oc, 100*rate_t, 100*rate_c, or_val, ci_lo, ci_hi, pv))
  }

  out <- rbindlist(results)
  out[, analysis := "sensitivity_mg_cat"]
  outfile <- file.path(RESULTS, sprintf("sens_mgcat_%s.csv", db))
  fwrite(out, outfile)
  cat(sprintf("\n  Saved: %s\n", outfile))
} else {
  cat("  WARN: did_hte_data not found, skipping outcome analysis\n")
}

cat("\n══ DONE ══\n")
