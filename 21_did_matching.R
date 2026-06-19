#!/usr/bin/env Rscript
# ============================================================================
# 21_did_matching.R — 1:4 PSM + temporal Cr alignment for DiD
#
# For each database:
#   1. Combine treated + control, estimate PS on 28 covariates
#   2. 1:4 nearest-neighbor matching (caliper 0.2 SD logit-PS, with replacement)
#   3. Transfer treated patient's Cr measurement times to matched controls
#   4. Find control's closest Cr to those times (±4h tolerance)
#   5. Compute ΔCr for controls, report match quality
#
# Inputs:
#   results/20_did_treated_{db}.csv
#   results/20_did_control_{db}.csv
#   results/20_did_cr_all_{db}.csv
#
# Outputs:
#   results/21_did_matched_{db}.csv     — matched dataset (treated + controls)
#   results/21_did_match_quality_{db}.csv — SMD + temporal gap report
#
# Run:  Rscript 21_did_matching.R            # both
#       Rscript 21_did_matching.R eicu       # eICU only
#       Rscript 21_did_matching.R mimic      # MIMIC only
# ============================================================================

suppressPackageStartupMessages({
  library(MatchIt)
  library(tableone)
})

RESULTS <- path.expand("~/mg_aki/results")
TOLERANCE_H <- 4    # ±4h temporal tolerance for Cr alignment
PRIMARY_WINDOW <- "6_24h"

# ── PS covariates (28 — no vasopressor/transfusion) ─────────────────────
PS_COVARS <- c(
  "age", "is_female", "bmi",
  "surg_cabg", "surg_valve", "surg_combined",
  "heart_failure", "hypertension", "diabetes", "ckd",
  "copd", "pvd", "stroke", "liver_disease",
  "egfr",
  "loop_diuretics", "nsaids", "acei_arb", "ppi",
  "beta_blockers", "steroids", "antiarrhythmics",
  "first_potassium", "first_calcium", "first_heartrate",
  "first_mg_value", "first_lactate", "lactate_missing"
)

# ── Helpers ──────────────────────────────────────────────────────────────
median_impute <- function(d, vars) {
  for (v in vars) {
    if (v %in% names(d) && any(is.na(d[[v]]))) {
      med <- median(d[[v]], na.rm = TRUE)
      n_miss <- sum(is.na(d[[v]]))
      d[[v]][is.na(d[[v]])] <- med
      cat(sprintf("    Imputed %s: %d missing → median %.2f\n", v, n_miss, med))
    }
  }
  d
}

smd_table <- function(d, vars, trt_var = "treated") {
  # Compute standardized mean differences
  out <- data.frame(variable = vars, smd_raw = NA_real_, smd_matched = NA_real_,
                    stringsAsFactors = FALSE)
  for (i in seq_along(vars)) {
    v <- vars[i]
    if (!v %in% names(d)) next
    x1 <- d[[v]][d[[trt_var]] == 1]
    x0 <- d[[v]][d[[trt_var]] == 0]
    sp <- sqrt((var(x1, na.rm = TRUE) + var(x0, na.rm = TRUE)) / 2)
    if (is.na(sp) || sp < 1e-10) next
    out$smd_raw[i] <- abs(mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) / sp
  }
  out
}

# ============================================================================
# MAIN: run_matching() — one database at a time
# ============================================================================
run_matching <- function(db) {
  tag <- tolower(db)
  SEP <- paste(rep("=", 70), collapse = "")
  cat(sprintf("\n%s\n%s: 1:4 PSM + Temporal Cr Alignment\n%s\n", SEP, db, SEP))

  # ── Load data ────────────────────────────────────────────────────────────
  trt_path <- file.path(RESULTS, sprintf("20_did_treated_%s.csv", tag))
  ctl_path <- file.path(RESULTS, sprintf("20_did_control_%s.csv", tag))
  cr_path  <- file.path(RESULTS, sprintf("20_did_cr_all_%s.csv", tag))

  if (!file.exists(trt_path)) { cat("  File not found:", trt_path, "\n"); return(NULL) }

  trt <- read.csv(trt_path, stringsAsFactors = FALSE)
  ctl <- read.csv(ctl_path, stringsAsFactors = FALSE)
  cr_all <- read.csv(cr_path, stringsAsFactors = FALSE)

  cat(sprintf("  Treated: %d  Controls: %d  Cr measurements: %d\n",
              nrow(trt), nrow(ctl), nrow(cr_all)))

  # ── Harmonize ID column ──────────────────────────────────────────────────
  # eICU uses patientunitstayid, MIMIC uses stay_id
  id_col <- if ("patientunitstayid" %in% names(trt)) "patientunitstayid" else "stay_id"
  trt$pid <- trt[[id_col]]
  ctl$pid <- ctl[[id_col]]

  # Same for Cr_all
  cr_id_col <- if ("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
  cr_all$pid <- cr_all[[cr_id_col]]

  # Cr_all offset column: labresultoffset (minutes)
  if (!"labresultoffset" %in% names(cr_all) && "offset_min" %in% names(cr_all))
    cr_all$labresultoffset <- cr_all$offset_min

  # ── Align columns for stacking ──────────────────────────────────────────
  # Control needs egfr computed from first_postop_cr
  if (!"egfr" %in% names(ctl) && "first_postop_cr" %in% names(ctl)) {
    # Already computed in ETL
  }

  # Stack treated + control
  common_cols <- intersect(names(trt), names(ctl))
  # Ensure key columns exist
  needed <- c("pid", "treated", PS_COVARS)
  missing_trt <- setdiff(needed, names(trt))
  missing_ctl <- setdiff(needed, names(ctl))
  if (length(missing_trt) > 0)
    cat(sprintf("  ⚠ Missing in treated: %s\n", paste(missing_trt, collapse = ", ")))
  if (length(missing_ctl) > 0)
    cat(sprintf("  ⚠ Missing in control: %s\n", paste(missing_ctl, collapse = ", ")))

  # Use available covariates only
  ps_vars <- intersect(PS_COVARS, intersect(names(trt), names(ctl)))
  cat(sprintf("  PS covariates available: %d / %d\n", length(ps_vars), length(PS_COVARS)))

  stack_cols <- unique(c("pid", "treated", ps_vars))
  stack_cols <- intersect(stack_cols, common_cols)
  combined <- rbind(trt[, stack_cols], ctl[, stack_cols])
  rownames(combined) <- seq_len(nrow(combined))
  cat(sprintf("  Combined: %d rows (%d treated + %d control)\n",
              nrow(combined), sum(combined$treated == 1), sum(combined$treated == 0)))

  # ── Median imputation ───────────────────────────────────────────────────
  cat("\n  Imputation:\n")
  combined <- median_impute(combined, ps_vars)

  # ── PS model + matching ─────────────────────────────────────────────────
  cat("\n  Propensity score matching...\n")
  ps_formula <- as.formula(paste("treated ~", paste(ps_vars, collapse = " + ")))

  # MatchIt: 1:4, nearest neighbor, caliper 0.2 SD, with replacement
  m <- matchit(ps_formula, data = combined, method = "nearest",
               distance = "glm", ratio = 4, caliper = 0.2,
               replace = TRUE)

  cat(sprintf("  Matched: %d treated, %d control uses\n",
              sum(m$treat), sum(!is.na(m$subclass) & m$treat == 0)))

  # Extract matched data
  md <- match.data(m)
  n_trt_matched <- sum(md$treated == 1)
  n_ctl_matched <- sum(md$treated == 0)
  cat(sprintf("  Matched dataset: %d treated, %d control rows\n",
              n_trt_matched, n_ctl_matched))

  # Unmatched treated
  n_unmatched <- sum(combined$treated == 1) - n_trt_matched
  if (n_unmatched > 0)
    cat(sprintf("  ⚠ Unmatched treated (outside caliper): %d\n", n_unmatched))

  # ── Covariate balance ───────────────────────────────────────────────────
  cat("\n  Covariate balance (SMD):\n")
  bal <- data.frame(variable = ps_vars, smd_raw = NA_real_,
                    smd_matched = NA_real_, stringsAsFactors = FALSE)
  for (i in seq_along(ps_vars)) {
    v <- ps_vars[i]
    if (!v %in% names(combined)) next
    # Raw
    x1 <- combined[[v]][combined$treated == 1]
    x0 <- combined[[v]][combined$treated == 0]
    sp <- sqrt((var(x1, na.rm = TRUE) + var(x0, na.rm = TRUE)) / 2)
    if (!is.na(sp) && sp > 1e-10)
      bal$smd_raw[i] <- abs(mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) / sp
    # Matched
    x1m <- md[[v]][md$treated == 1]
    x0m <- md[[v]][md$treated == 0]
    spm <- sqrt((var(x1m, na.rm = TRUE) + var(x0m, na.rm = TRUE)) / 2)
    if (!is.na(spm) && spm > 1e-10)
      bal$smd_matched[i] <- abs(mean(x1m, na.rm = TRUE) - mean(x0m, na.rm = TRUE)) / spm
  }
  bal <- bal[!is.na(bal$smd_raw), ]
  cat(sprintf("    Max raw SMD:     %.3f (%s)\n",
              max(bal$smd_raw, na.rm = TRUE),
              bal$variable[which.max(bal$smd_raw)]))
  cat(sprintf("    Max matched SMD: %.3f (%s)\n",
              max(bal$smd_matched, na.rm = TRUE),
              bal$variable[which.max(bal$smd_matched)]))
  cat(sprintf("    SMD > 0.10 after matching: %d / %d\n",
              sum(bal$smd_matched > 0.10, na.rm = TRUE), nrow(bal)))
  cat(sprintf("    SMD > 0.05 after matching: %d / %d\n",
              sum(bal$smd_matched > 0.05, na.rm = TRUE), nrow(bal)))

  # ════════════════════════════════════════════════════════════════════════
  # TEMPORAL Cr ALIGNMENT
  # ════════════════════════════════════════════════════════════════════════
  cat(sprintf("\n%s\nTemporal Cr alignment (tolerance: ±%dh)\n%s\n",
              paste(rep("─", 50), collapse = ""), TOLERANCE_H,
              paste(rep("─", 50), collapse = "")))

  # Get match pairs: for each treated patient, which controls were matched?
  # MatchIt subclass links them
  trt_pids_matched <- md$pid[md$treated == 1]

  # Build explicit pair table from MatchIt internals
  match_matrix <- m$match.matrix  # rows = treated, cols = 1:ratio
  trt_indices <- as.integer(rownames(match_matrix))
  pairs <- list()

  for (i in seq_len(nrow(match_matrix))) {
    trt_idx <- trt_indices[i]
    trt_id <- combined$pid[trt_idx]
    for (j in seq_len(ncol(match_matrix))) {
      ctl_idx <- match_matrix[i, j]
      if (!is.na(ctl_idx)) {
        ctl_id <- combined$pid[as.integer(ctl_idx)]
        pairs[[length(pairs) + 1]] <- data.frame(
          trt_pid = trt_id, ctl_pid = ctl_id,
          stringsAsFactors = FALSE)
      }
    }
  }
  pairs_df <- do.call(rbind, pairs)
  cat(sprintf("  Total match pairs: %d\n", nrow(pairs_df)))
  cat(sprintf("  Unique treated: %d, unique controls: %d\n",
              length(unique(pairs_df$trt_pid)),
              length(unique(pairs_df$ctl_pid))))

  # ── For each window, do temporal alignment ──────────────────────────────
  windows <- c("6_24h", "6_48h", "0_24h")
  cr_offset_cols <- paste0("cr_post_offset_", windows)
  tolerance_min <- TOLERANCE_H * 60

  # Pre-index Cr measurements by pid for fast lookup
  cr_list <- split(cr_all[, c("pid", "labresult", "labresultoffset")],
                   cr_all$pid)

  results_by_window <- list()

  for (wname in windows) {
    cat(sprintf("\n── Window: %s%s ──\n", wname,
                if (wname == PRIMARY_WINDOW) " (PRIMARY)" else ""))

    offset_col <- paste0("cr_post_offset_", wname)
    cr_post_col <- paste0("cr_post_", wname)
    delta_col <- paste0("delta_cr_", wname)

    # Get treated patients with valid Cr_post for this window
    trt_valid <- trt[!is.na(trt[[cr_post_col]]),
                     c("pid", "cr_pre", "cr_pre_offset_min",
                       cr_post_col, offset_col, delta_col)]
    names(trt_valid)[names(trt_valid) == cr_post_col] <- "cr_post"
    names(trt_valid)[names(trt_valid) == offset_col] <- "cr_post_offset_min"
    names(trt_valid)[names(trt_valid) == delta_col] <- "delta_cr"
    trt_valid$treated <- 1
    trt_valid$match_pair_id <- seq_len(nrow(trt_valid))

    cat(sprintf("  Treated with Cr_post: %d\n", nrow(trt_valid)))

    # For each pair, find control's Cr at treated patient's time points
    n_valid <- 0
    n_pre_fail <- 0
    n_post_fail <- 0
    n_no_cr <- 0
    ctl_rows <- list()
    temporal_gaps <- list()

    pairs_for_window <- pairs_df[pairs_df$trt_pid %in% trt_valid$pid, ]
    cat(sprintf("  Match pairs to align: %d\n", nrow(pairs_for_window)))

    for (r in seq_len(nrow(pairs_for_window))) {
      tpid <- pairs_for_window$trt_pid[r]
      cpid <- pairs_for_window$ctl_pid[r]

      # Get treated patient's time points
      tidx <- which(trt_valid$pid == tpid)[1]
      if (is.na(tidx)) next
      t_pre_min <- trt_valid$cr_pre_offset_min[tidx]
      t_post_min <- trt_valid$cr_post_offset_min[tidx]
      trt_pair_id <- trt_valid$match_pair_id[tidx]

      # Get control's Cr measurements
      ctl_cr <- cr_list[[as.character(cpid)]]
      if (is.null(ctl_cr) || nrow(ctl_cr) < 2) {
        n_no_cr <- n_no_cr + 1
        next
      }

      # Find closest Cr to t_pre_min
      pre_diffs <- abs(ctl_cr$labresultoffset - t_pre_min)
      pre_idx <- which.min(pre_diffs)
      pre_gap_min <- pre_diffs[pre_idx]

      if (pre_gap_min > tolerance_min) {
        n_pre_fail <- n_pre_fail + 1
        next
      }

      # Find closest Cr to t_post_min (must be DIFFERENT from pre measurement
      # and must be AFTER the pre measurement)
      post_candidates <- ctl_cr[ctl_cr$labresultoffset > ctl_cr$labresultoffset[pre_idx], ]
      if (nrow(post_candidates) == 0) {
        n_post_fail <- n_post_fail + 1
        next
      }
      post_diffs <- abs(post_candidates$labresultoffset - t_post_min)
      post_idx <- which.min(post_diffs)
      post_gap_min <- post_diffs[post_idx]

      if (post_gap_min > tolerance_min) {
        n_post_fail <- n_post_fail + 1
        next
      }

      # Valid match
      cr_pre_ctl <- ctl_cr$labresult[pre_idx]
      cr_post_ctl <- post_candidates$labresult[post_idx]
      cr_pre_offset_ctl <- ctl_cr$labresultoffset[pre_idx]
      cr_post_offset_ctl <- post_candidates$labresultoffset[post_idx]

      n_valid <- n_valid + 1
      ctl_rows[[n_valid]] <- data.frame(
        pid = cpid,
        cr_pre = cr_pre_ctl,
        cr_pre_offset_min = cr_pre_offset_ctl,
        cr_post = cr_post_ctl,
        cr_post_offset_min = cr_post_offset_ctl,
        delta_cr = cr_post_ctl - cr_pre_ctl,
        treated = 0,
        match_pair_id = trt_pair_id,
        stringsAsFactors = FALSE)

      temporal_gaps[[n_valid]] <- data.frame(
        pre_gap_h = pre_gap_min / 60,
        post_gap_h = post_gap_min / 60,
        stringsAsFactors = FALSE)
    }

    cat(sprintf("  Valid temporal matches: %d / %d (%.1f%%)\n",
                n_valid, nrow(pairs_for_window),
                100 * n_valid / max(nrow(pairs_for_window), 1)))
    cat(sprintf("  Failed: no Cr=%d, pre-tolerance=%d, post-tolerance=%d\n",
                n_no_cr, n_pre_fail, n_post_fail))

    if (n_valid == 0) {
      cat("  ⚠ No valid matches — skipping window\n")
      next
    }

    ctl_matched <- do.call(rbind, ctl_rows)
    gaps <- do.call(rbind, temporal_gaps)

    cat(sprintf("  Temporal gaps (h): pre median=%.1f, post median=%.1f\n",
                median(gaps$pre_gap_h), median(gaps$post_gap_h)))
    cat(sprintf("    Pre  P25/P75: %.1f / %.1f\n",
                quantile(gaps$pre_gap_h, 0.25), quantile(gaps$pre_gap_h, 0.75)))
    cat(sprintf("    Post P25/P75: %.1f / %.1f\n",
                quantile(gaps$post_gap_h, 0.25), quantile(gaps$post_gap_h, 0.75)))

    # Effective match ratio
    n_trt_with_match <- length(unique(ctl_matched$match_pair_id))
    matches_per_trt <- table(ctl_matched$match_pair_id)
    cat(sprintf("  Treated with ≥1 temporal match: %d / %d\n",
                n_trt_with_match, nrow(trt_valid)))
    cat(sprintf("  Matches per treated: mean=%.1f, median=%d, range=[%d, %d]\n",
                mean(matches_per_trt), median(matches_per_trt),
                min(matches_per_trt), max(matches_per_trt)))

    # Combine treated + matched controls
    matched_dataset <- rbind(
      trt_valid[trt_valid$match_pair_id %in% unique(ctl_matched$match_pair_id), ],
      ctl_matched
    )

    # Merge covariates back
    covar_cols <- c("pid", ps_vars, "surgery_type", "hosp_mortality")
    covar_cols <- intersect(covar_cols, common_cols)
    matched_dataset <- merge(matched_dataset, combined[, covar_cols],
                             by = "pid", all.x = TRUE)

    # ΔCr summary
    cat(sprintf("\n  ΔCr summary:\n"))
    for (grp in c(1, 0)) {
      sub <- matched_dataset$delta_cr[matched_dataset$treated == grp]
      label <- if (grp == 1) "Treated" else "Control"
      cat(sprintf("    %s: n=%d, mean=%.4f, median=%.3f, SD=%.3f\n",
                  label, length(sub), mean(sub, na.rm=TRUE),
                  median(sub, na.rm=TRUE), sd(sub, na.rm=TRUE)))
    }

    # Crude DiD estimate
    did_raw <- mean(matched_dataset$delta_cr[matched_dataset$treated == 1], na.rm = TRUE) -
               mean(matched_dataset$delta_cr[matched_dataset$treated == 0], na.rm = TRUE)
    cat(sprintf("  Crude DiD (ΔΔCr): %.4f mg/dL\n", did_raw))

    # Quick t-test
    tt <- t.test(delta_cr ~ treated, data = matched_dataset)
    cat(sprintf("  t-test: diff=%.4f, P=%.4f, 95%% CI [%.4f, %.4f]\n",
                -diff(tt$estimate), tt$p.value, -tt$conf.int[2], -tt$conf.int[1]))

    # Save
    out_path <- file.path(RESULTS, sprintf("21_did_matched_%s_%s.csv", tag, wname))
    write.csv(matched_dataset, out_path, row.names = FALSE)
    cat(sprintf("  Saved: %s (%d rows)\n", basename(out_path), nrow(matched_dataset)))

    results_by_window[[wname]] <- list(
      n_trt = n_trt_with_match,
      n_ctl = n_valid,
      did_raw = did_raw,
      p = tt$p.value,
      ci = c(-tt$conf.int[2], -tt$conf.int[1])
    )
  }

  # ── Save balance table ─────────────────────────────────────────────────
  bal_path <- file.path(RESULTS, sprintf("21_did_match_quality_%s.csv", tag))
  write.csv(bal, bal_path, row.names = FALSE)
  cat(sprintf("\n  Balance table: %s\n", basename(bal_path)))

  # ── Summary ────────────────────────────────────────────────────────────
  cat(sprintf("\n%s\n%s: MATCHING SUMMARY\n%s\n", SEP, db, SEP))
  for (wname in names(results_by_window)) {
    r <- results_by_window[[wname]]
    primary <- if (wname == PRIMARY_WINDOW) " ◀ PRIMARY" else ""
    cat(sprintf("  %s: %d treated, %d control matches, DiD=%.4f, P=%.4f%s\n",
                wname, r$n_trt, r$n_ctl, r$did_raw, r$p, primary))
  }

  return(results_by_window)
}

# ============================================================================
# ENTRY POINT
# ============================================================================
cat("======================================================================\n")
cat("21_did_matching.R — 1:4 PSM + Temporal Cr Alignment\n")
cat(sprintf("  Tolerance: ±%dh | Primary window: %s\n", TOLERANCE_H, PRIMARY_WINDOW))
cat("======================================================================\n")

args <- commandArgs(trailingOnly = TRUE)
run_all <- length(args) == 0

if (run_all || "eicu" %in% tolower(args))
  run_matching("eICU")

if (run_all || "mimic" %in% tolower(args))
  run_matching("MIMIC")

cat("\n======================================================================\n")
cat("NEXT: Rscript 22_did_analysis.R — formal DiD with clustered SEs,\n")
cat("      subgroup analyses, Mg-stratified effects\n")
cat("======================================================================\n")
