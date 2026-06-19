#!/usr/bin/env Rscript
# ============================================================================
# 21b_did_sweep.R — Parameter sweep for PSM + temporal alignment
#
# Sweeps across matching and temporal parameters to understand tradeoffs:
#   Phase 1: PS matching quality (caliper × ratio × replace)
#   Phase 2: Temporal alignment rate (tolerance for each PS config)
#
# Output: results/21b_sweep_{db}.csv — one row per configuration
#
# Run:  Rscript 21b_did_sweep.R eicu
#       Rscript 21b_did_sweep.R mimic
# ============================================================================

suppressPackageStartupMessages(library(MatchIt))

RESULTS <- path.expand("~/mg_aki/results")
PRIMARY_WINDOW <- "6_24h"

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

# ── Sweep grid ───────────────────────────────────────────────────────────
CALIPER    <- 0.2                         # fixed
RATIOS     <- c(1, 2, 4)
TOLERANCES <- c(2, 4, 6, 8, 12)          # hours

# ── Helpers ──────────────────────────────────────────────────────────────
median_impute <- function(d, vars) {
  for (v in vars)
    if (v %in% names(d) && any(is.na(d[[v]])))
      d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm = TRUE)
  d
}

compute_smd <- function(d, vars, trt_var = "treated") {
  smds <- numeric(length(vars))
  for (i in seq_along(vars)) {
    v <- vars[i]
    if (!v %in% names(d)) { smds[i] <- NA; next }
    x1 <- d[[v]][d[[trt_var]] == 1]
    x0 <- d[[v]][d[[trt_var]] == 0]
    sp <- sqrt((var(x1, na.rm = TRUE) + var(x0, na.rm = TRUE)) / 2)
    smds[i] <- if (!is.na(sp) && sp > 1e-10)
      abs(mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) / sp else NA
  }
  names(smds) <- vars
  smds
}

# ============================================================================
run_sweep <- function(db) {
  tag <- tolower(db)
  cat(sprintf("\n%s\n%s: Parameter Sweep\n%s\n",
              paste(rep("=", 70), collapse = ""), db,
              paste(rep("=", 70), collapse = "")))

  # ── Load ─────────────────────────────────────────────────────────────────
  trt <- read.csv(file.path(RESULTS, sprintf("20_did_treated_%s.csv", tag)),
                  stringsAsFactors = FALSE)
  ctl <- read.csv(file.path(RESULTS, sprintf("20_did_control_%s.csv", tag)),
                  stringsAsFactors = FALSE)
  cr_all <- read.csv(file.path(RESULTS, sprintf("20_did_cr_all_%s.csv", tag)),
                     stringsAsFactors = FALSE)

  id_col <- if ("patientunitstayid" %in% names(trt)) "patientunitstayid" else "stay_id"
  trt$pid <- trt[[id_col]]; ctl$pid <- ctl[[id_col]]
  cr_id <- if ("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
  cr_all$pid <- cr_all[[cr_id]]
  if (!"labresultoffset" %in% names(cr_all) && "offset_min" %in% names(cr_all))
    cr_all$labresultoffset <- cr_all$offset_min

  ps_vars <- intersect(PS_COVARS, intersect(names(trt), names(ctl)))
  stack_cols <- intersect(unique(c("pid", "treated", ps_vars)),
                          intersect(names(trt), names(ctl)))

  combined <- rbind(trt[, stack_cols], ctl[, stack_cols])
  rownames(combined) <- seq_len(nrow(combined))
  combined <- median_impute(combined, ps_vars)

  cat(sprintf("  %d treated + %d control = %d total\n",
              sum(combined$treated == 1), sum(combined$treated == 0), nrow(combined)))
  cat(sprintf("  PS covariates: %d\n", length(ps_vars)))

  # Raw SMDs
  raw_smds <- compute_smd(combined, ps_vars)
  cat(sprintf("  Raw max SMD: %.3f (%s)\n",
              max(raw_smds, na.rm = TRUE), names(which.max(raw_smds))))

  # Pre-index Cr by pid
  cr_list <- split(cr_all[, c("pid", "labresult", "labresultoffset")], cr_all$pid)

  # Primary window Cr_post info for treated
  offset_col <- paste0("cr_post_offset_", PRIMARY_WINDOW)
  cr_post_col <- paste0("cr_post_", PRIMARY_WINDOW)
  trt_valid <- trt[!is.na(trt[[cr_post_col]]),
                   c("pid", "cr_pre", "cr_pre_offset_min", cr_post_col, offset_col)]
  names(trt_valid)[4:5] <- c("cr_post", "cr_post_offset_min")
  cat(sprintf("  Treated with Cr_post (%s): %d\n", PRIMARY_WINDOW, nrow(trt_valid)))

  ps_formula <- as.formula(paste("treated ~", paste(ps_vars, collapse = " + ")))

  # ── Sweep: ratio × tolerance ─────────────────────────────────────────
  results <- list()
  idx <- 0
  n_configs <- length(RATIOS) * length(TOLERANCES)
  cat(sprintf("\n  Sweeping %d ratio configs × %d tolerances = %d total...\n",
              length(RATIOS), length(TOLERANCES), n_configs))

  for (ratio in RATIOS) {

    config_label <- sprintf("cal=0.2_r=%d_rep=T", ratio)

    # PS matching (fixed: caliper=0.2, replace=TRUE)
    m <- tryCatch(
      matchit(ps_formula, data = combined, method = "nearest",
              distance = "glm", ratio = ratio, caliper = CALIPER,
              replace = TRUE),
      error = function(e) NULL,
      warning = function(w) {
        suppressWarnings(
          matchit(ps_formula, data = combined, method = "nearest",
                  distance = "glm", ratio = ratio, caliper = CALIPER,
                  replace = TRUE))
      })
    if (is.null(m)) {
      cat(sprintf("    %s: FAILED\n", config_label)); next
    }

    md <- match.data(m)
    n_trt_m <- sum(md$treated == 1)
    n_ctl_m <- sum(md$treated == 0)
    if (n_trt_m == 0 || n_ctl_m == 0) {
      cat(sprintf("    %s: empty match\n", config_label)); next
    }

    # Balance
    matched_smds <- compute_smd(md, ps_vars)
    max_smd <- max(matched_smds, na.rm = TRUE)
    mean_smd <- mean(matched_smds, na.rm = TRUE)
    n_above_10 <- sum(matched_smds > 0.10, na.rm = TRUE)

    # Extract match pairs
    mm <- m$match.matrix
    if (is.null(mm)) next
    trt_idx <- as.integer(rownames(mm))
    pairs <- list()
    for (i in seq_len(nrow(mm))) {
      for (j in seq_len(ncol(mm))) {
        ci <- mm[i, j]
        if (!is.na(ci))
          pairs[[length(pairs) + 1]] <- c(combined$pid[trt_idx[i]],
                                            combined$pid[as.integer(ci)])
      }
    }
    if (length(pairs) == 0) next
    pairs_mat <- do.call(rbind, pairs)
    colnames(pairs_mat) <- c("trt_pid", "ctl_pid")
    pairs_df <- as.data.frame(pairs_mat)
    pairs_valid <- pairs_df[pairs_df$trt_pid %in% trt_valid$pid, ]

    cat(sprintf("    %s: n_trt=%d, n_ctl=%d, max_SMD=%.3f, n_smd>0.1=%d, pairs=%d\n",
                config_label, n_trt_m, n_ctl_m, max_smd, n_above_10, nrow(pairs_valid)))

    # Sweep temporal tolerances for this ratio
    for (tol_h in TOLERANCES) {
      tol_min <- tol_h * 60
      n_valid <- 0; n_pre_fail <- 0; n_post_fail <- 0; n_no_cr <- 0
      did_ctl_sum <- 0; did_ctl_n <- 0
      matched_trt_pids <- c()
      pre_gaps <- c(); post_gaps <- c()

      for (r in seq_len(nrow(pairs_valid))) {
        tpid <- pairs_valid$trt_pid[r]
        cpid <- pairs_valid$ctl_pid[r]

        tidx <- which(trt_valid$pid == tpid)[1]
        if (is.na(tidx)) next
        t_pre <- trt_valid$cr_pre_offset_min[tidx]
        t_post <- trt_valid$cr_post_offset_min[tidx]

        ctl_cr <- cr_list[[as.character(cpid)]]
        if (is.null(ctl_cr) || nrow(ctl_cr) < 2) {
          n_no_cr <- n_no_cr + 1; next
        }

        # Pre alignment
        pd <- abs(ctl_cr$labresultoffset - t_pre)
        pi <- which.min(pd)
        if (pd[pi] > tol_min) { n_pre_fail <- n_pre_fail + 1; next }

        # Post alignment (must be after pre)
        pc <- ctl_cr[ctl_cr$labresultoffset > ctl_cr$labresultoffset[pi], ]
        if (nrow(pc) == 0) { n_post_fail <- n_post_fail + 1; next }
        pod <- abs(pc$labresultoffset - t_post)
        poi <- which.min(pod)
        if (pod[poi] > tol_min) { n_post_fail <- n_post_fail + 1; next }

        n_valid <- n_valid + 1
        matched_trt_pids <- c(matched_trt_pids, tpid)
        pre_gaps <- c(pre_gaps, pd[pi] / 60)
        post_gaps <- c(post_gaps, pod[poi] / 60)
        did_ctl_sum <- did_ctl_sum + (pc$labresult[poi] - ctl_cr$labresult[pi])
        did_ctl_n <- did_ctl_n + 1
      }

      trt_matched <- trt_valid[trt_valid$pid %in% unique(matched_trt_pids), ]
      n_trt_temporal <- nrow(trt_matched)
      did_trt <- if (n_trt_temporal > 0)
        mean(trt_matched$cr_post - trt_matched$cr_pre) else NA
      did_ctl <- if (did_ctl_n > 0) did_ctl_sum / did_ctl_n else NA
      did_est <- if (!is.na(did_trt) && !is.na(did_ctl)) did_trt - did_ctl else NA

      idx <- idx + 1
      results[[idx]] <- data.frame(
        ratio = ratio,
        tolerance_h = tol_h,
        n_trt_ps_matched = n_trt_m,
        n_ctl_ps_matched = n_ctl_m,
        max_smd = round(max_smd, 3),
        mean_smd = round(mean_smd, 3),
        n_smd_above_10 = n_above_10,
        n_pairs = nrow(pairs_valid),
        n_temporal_valid = n_valid,
        temporal_rate_pct = round(100 * n_valid / max(nrow(pairs_valid), 1), 1),
        n_pre_fail = n_pre_fail,
        n_post_fail = n_post_fail,
        n_trt_with_match = n_trt_temporal,
        median_pre_gap_h = if (length(pre_gaps) > 0) round(median(pre_gaps), 1) else NA,
        median_post_gap_h = if (length(post_gaps) > 0) round(median(post_gaps), 1) else NA,
        did_estimate = if (!is.na(did_est)) round(did_est, 4) else NA,
        stringsAsFactors = FALSE
      )
    }
  }

  # ── Compile and save ───────────────────────────────────────────────────
  sweep <- do.call(rbind, results)

  out_path <- file.path(RESULTS, sprintf("21b_sweep_%s.csv", tag))
  write.csv(sweep, out_path, row.names = FALSE)
  cat(sprintf("\n  Saved: %s (%d configurations)\n", basename(out_path), nrow(sweep)))

  # ── Print top configurations ───────────────────────────────────────────
  cat(sprintf("\n%s\nRESULTS (caliper=0.2, replace=TRUE, %s window)\n%s\n",
              paste(rep("─", 60), collapse = ""), PRIMARY_WINDOW,
              paste(rep("─", 60), collapse = "")))

  # Full table
  cat("\n  ratio  tol_h  ps_trt  max_smd  smd>0.1  temporal  trt_matched  DiD\n")
  cat("  ───── ────── ────── ──────── ──────── ──────── ─────────── ─────────\n")
  for (i in seq_len(nrow(sweep))) {
    r <- sweep[i, ]
    cat(sprintf("  %3d    %3dh   %5d    %.3f      %2d     %5d     %5d     %s\n",
                r$ratio, r$tolerance_h, r$n_trt_ps_matched,
                r$max_smd, r$n_smd_above_10, r$n_temporal_valid,
                r$n_trt_with_match,
                if (is.na(r$did_estimate)) "   NA" else sprintf("%+.4f", r$did_estimate)))
  }

  # Highlight: best balance configs with decent yield
  cat("\n  Recommended configs:\n")
  for (th in TOLERANCES) {
    sub <- sweep[sweep$tolerance_h == th & !is.na(sweep$max_smd), ]
    if (nrow(sub) == 0) next
    best <- sub[which.min(sub$max_smd), ]
    most <- sub[which.max(sub$n_trt_with_match), ]
    cat(sprintf("    ±%2dh: best balance r=%d (SMD=%.3f, n=%d) | most yield r=%d (SMD=%.3f, n=%d)\n",
                th, best$ratio, best$max_smd, best$n_trt_with_match,
                most$ratio, most$max_smd, most$n_trt_with_match))
  }

  return(sweep)
}

# ============================================================================
# ENTRY
# ============================================================================
cat("======================================================================\n")
cat("21b_did_sweep.R — Ratio × Tolerance sweep (caliper=0.2, replace=T)\n")
cat(sprintf("  Ratios: %s | Tolerances: ±%s h | Window: %s\n",
            paste(RATIOS, collapse = ", "), paste(TOLERANCES, collapse = ", "),
            PRIMARY_WINDOW))
cat(sprintf("  Total configs: %d\n", length(RATIOS) * length(TOLERANCES)))
cat("======================================================================\n")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  cat("Usage: Rscript 21b_did_sweep.R eicu|mimic\n")
  quit(status = 1)
}

for (a in args) run_sweep(toupper(a))

cat("\nDone. Review the CSV and the tradeoff table above.\n")
cat("Pick the config, then update 21_did_matching.R with those parameters.\n")
