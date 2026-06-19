#!/usr/bin/env Rscript
# ============================================================================
# 21_did_final.R — Production DiD: matching + doubly robust analysis
#
# Design:
#   PS matching (caliper=0.2, replace=TRUE) → temporal Cr alignment →
#   doubly robust DiD (regression-adjusted within matched sample)
#
# Configurations:
#   Primary:       r=1, ±6h, first Cr 6-24h
#   Sensitivity 1: r=1, ±4h, first Cr 6-24h  (tighter temporal)
#   Sensitivity 2: r=2, ±6h, first Cr 6-24h  (more power)
#   Sensitivity 3: r=4, ±4h, first Cr 6-24h  (max power, tight temporal)
#
# Cr_post strategies (each run through primary config):
#   first:  first Cr 6-24h after IV Mg (primary)
#   peak:   max Cr 6-48h after IV Mg
#   at_24h: Cr closest to 24h after IV Mg
#
# Run:  Rscript 21_did_final.R eicu
#       Rscript 21_did_final.R mimic
# ============================================================================

suppressPackageStartupMessages({
  library(MatchIt)
  library(sandwich)
  library(lmtest)
})

RESULTS <- path.expand("~/mg_aki/results")
CALIPER <- 0.2

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
  for (v in vars)
    if (v %in% names(d) && any(is.na(d[[v]])))
      d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm = TRUE)
  d
}

compute_smds <- function(d, vars) {
  sapply(vars, function(v) {
    if (!v %in% names(d)) return(NA)
    x1 <- d[[v]][d$treated == 1]; x0 <- d[[v]][d$treated == 0]
    sp <- sqrt((var(x1, na.rm=T) + var(x0, na.rm=T)) / 2)
    if (is.na(sp) || sp < 1e-10) NA else abs(mean(x1,na.rm=T) - mean(x0,na.rm=T)) / sp
  })
}

# Doubly robust DiD: OLS within matched sample, cluster SEs by pair
did_robust <- function(matched_df, ps_vars, pair_col = "match_pair_id") {
  # Identify covariates with SMD > 0.05 for adjustment
  smds <- compute_smds(matched_df, ps_vars)
  adjust_vars <- names(smds[!is.na(smds) & smds > 0.05])

  # Unadjusted
  fit0 <- lm(delta_cr ~ treated, data = matched_df)
  cl0 <- if (length(unique(matched_df[[pair_col]])) > 1)
    vcovCL(fit0, cluster = matched_df[[pair_col]]) else vcovHC(fit0, type = "HC1")
  ct0 <- coeftest(fit0, vcov. = cl0)

  # Adjusted (doubly robust)
  avail <- intersect(adjust_vars, names(matched_df))
  if (length(avail) > 0) {
    fml <- as.formula(paste("delta_cr ~ treated +", paste(avail, collapse = " + ")))
    fit1 <- lm(fml, data = matched_df)
    cl1 <- if (length(unique(matched_df[[pair_col]])) > 1)
      vcovCL(fit1, cluster = matched_df[[pair_col]]) else vcovHC(fit1, type = "HC1")
    ct1 <- coeftest(fit1, vcov. = cl1)
  } else {
    ct1 <- ct0
    avail <- character(0)
  }

  list(
    n = nrow(matched_df),
    n_trt = sum(matched_df$treated == 1),
    n_ctl = sum(matched_df$treated == 0),
    n_adjust = length(avail),
    # Unadjusted
    did_unadj = ct0["treated", "Estimate"],
    se_unadj = ct0["treated", "Std. Error"],
    p_unadj = ct0["treated", "Pr(>|t|)"],
    # Doubly robust
    did_adj = ct1["treated", "Estimate"],
    se_adj = ct1["treated", "Std. Error"],
    p_adj = ct1["treated", "Pr(>|t|)"],
    ci_lo = ct1["treated", "Estimate"] - 1.96 * ct1["treated", "Std. Error"],
    ci_hi = ct1["treated", "Estimate"] + 1.96 * ct1["treated", "Std. Error"],
    adjust_vars = paste(avail, collapse = ", ")
  )
}

# ============================================================================
run_did <- function(db) {
  tag <- tolower(db)
  SEP <- paste(rep("=", 70), collapse = "")
  cat(sprintf("\n%s\n%s: DiD Matching + Doubly Robust Analysis\n%s\n", SEP, db, SEP))

  # ── Load ─────────────────────────────────────────────────────────────────
  trt <- read.csv(file.path(RESULTS, sprintf("20_did_treated_%s.csv", tag)), stringsAsFactors=F)
  ctl <- read.csv(file.path(RESULTS, sprintf("20_did_control_%s.csv", tag)), stringsAsFactors=F)
  cr_all <- read.csv(file.path(RESULTS, sprintf("20_did_cr_all_%s.csv", tag)), stringsAsFactors=F)

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

  cat(sprintf("  %d treated + %d control\n", sum(combined$treated==1), sum(combined$treated==0)))

  # Pre-index Cr
  cr_list <- split(cr_all[, c("pid","labresult","labresultoffset")], cr_all$pid)

  # ── Build Cr_post strategies for treated ─────────────────────────────────
  cat("\n── Cr_post strategies for treated ──\n")

  # All post-IV-Mg Cr for treated patients
  cr_trt <- cr_all[cr_all$pid %in% trt$pid, ]
  mg_times <- setNames(trt$mg_offset_min, trt$pid)
  cr_trt$mg_off <- mg_times[as.character(cr_trt$pid)]
  cr_trt$post_h <- (cr_trt$labresultoffset - cr_trt$mg_off) / 60

  strategies <- list()

  # Strategy 1: FIRST Cr 6-24h after IV Mg
  s1 <- cr_trt[cr_trt$post_h >= 6 & cr_trt$post_h <= 24, ]
  s1 <- s1[order(s1$pid, s1$labresultoffset), ]
  s1 <- s1[!duplicated(s1$pid), ]
  s1$strategy <- "first_6_24h"
  strategies[["first_6_24h"]] <- s1
  cat(sprintf("  first_6_24h: %d treated\n", nrow(s1)))

  # Strategy 2: PEAK Cr 6-48h after IV Mg
  s2 <- cr_trt[cr_trt$post_h >= 6 & cr_trt$post_h <= 48, ]
  s2 <- s2[order(s2$pid, -s2$labresult), ]
  s2 <- s2[!duplicated(s2$pid), ]
  s2$strategy <- "peak_6_48h"
  strategies[["peak_6_48h"]] <- s2
  cat(sprintf("  peak_6_48h:  %d treated\n", nrow(s2)))

  # Strategy 3: CLOSEST to 24h after IV Mg
  cr_trt$dist_24h <- abs(cr_trt$post_h - 24)
  s3 <- cr_trt[cr_trt$post_h >= 0 & cr_trt$post_h <= 48, ]
  s3 <- s3[order(s3$pid, s3$dist_24h), ]
  s3 <- s3[!duplicated(s3$pid), ]
  s3$strategy <- "closest_24h"
  strategies[["closest_24h"]] <- s3
  cat(sprintf("  closest_24h: %d treated\n", nrow(s3)))

  # Add cr_pre to each strategy
  cr_pre_map <- setNames(trt$cr_pre, trt$pid)
  cr_pre_off_map <- setNames(trt$cr_pre_offset_min, trt$pid)
  for (nm in names(strategies)) {
    strategies[[nm]]$cr_pre <- cr_pre_map[as.character(strategies[[nm]]$pid)]
    strategies[[nm]]$cr_pre_offset_min <- cr_pre_off_map[as.character(strategies[[nm]]$pid)]
    strategies[[nm]]$cr_post <- strategies[[nm]]$labresult
    strategies[[nm]]$cr_post_offset_min <- strategies[[nm]]$labresultoffset
    strategies[[nm]]$delta_cr <- strategies[[nm]]$cr_post - strategies[[nm]]$cr_pre
  }

  # ── Matching configurations ──────────────────────────────────────────────
  configs <- list(
    list(label="primary",  ratio=1, tol_h=6, cr_strategy="first_6_24h"),
    list(label="sens_t4",  ratio=1, tol_h=4, cr_strategy="first_6_24h"),
    list(label="sens_r2",  ratio=2, tol_h=6, cr_strategy="first_6_24h"),
    list(label="sens_r4t4",ratio=4, tol_h=4, cr_strategy="first_6_24h"),
    # Cr strategy sensitivities (use primary matching config)
    list(label="sens_peak", ratio=1, tol_h=6, cr_strategy="peak_6_48h"),
    list(label="sens_at24", ratio=1, tol_h=6, cr_strategy="closest_24h")
  )

  ps_formula <- as.formula(paste("treated ~", paste(ps_vars, collapse = " + ")))

  # Cache PS matches by ratio (avoid re-running MatchIt)
  match_cache <- list()
  all_results <- list()

  for (cfg in configs) {
    cat(sprintf("\n%s\n  Config: %s (r=%d, ±%dh, %s)\n%s\n",
                paste(rep("─",60), collapse=""), cfg$label,
                cfg$ratio, cfg$tol_h, cfg$cr_strategy,
                paste(rep("─",60), collapse="")))

    # ── PS matching (cached by ratio) ──────────────────────────────────
    r_key <- as.character(cfg$ratio)
    if (is.null(match_cache[[r_key]])) {
      cat("  Running MatchIt...\n")
      m <- suppressWarnings(matchit(ps_formula, data=combined, method="nearest",
                                     distance="glm", ratio=cfg$ratio,
                                     caliper=CALIPER, replace=TRUE))
      md <- match.data(m)

      # Balance
      smds <- compute_smds(md, ps_vars)
      cat(sprintf("  PS matched: %d trt, %d ctl, max SMD=%.3f, n>0.1=%d\n",
                  sum(md$treated==1), sum(md$treated==0),
                  max(smds,na.rm=T), sum(smds>0.1,na.rm=T)))

      # Extract pairs
      mm <- m$match.matrix
      trt_idx <- as.integer(rownames(mm))
      pairs <- list()
      for (i in seq_len(nrow(mm)))
        for (j in seq_len(ncol(mm))) {
          ci <- mm[i,j]
          if (!is.na(ci))
            pairs[[length(pairs)+1]] <- c(combined$pid[trt_idx[i]],
                                           combined$pid[as.integer(ci)])
        }
      pairs_df <- as.data.frame(do.call(rbind, pairs))
      names(pairs_df) <- c("trt_pid","ctl_pid")

      match_cache[[r_key]] <- list(m=m, md=md, smds=smds, pairs=pairs_df)
    } else {
      cat("  Using cached match (r=", r_key, ")\n")
    }

    mc <- match_cache[[r_key]]

    # ── Get treated Cr_post for this strategy ──────────────────────────
    strat <- strategies[[cfg$cr_strategy]]
    trt_valid <- strat[, c("pid","cr_pre","cr_pre_offset_min","cr_post",
                            "cr_post_offset_min","delta_cr")]
    trt_valid$treated <- 1

    pairs_valid <- mc$pairs[mc$pairs$trt_pid %in% trt_valid$pid, ]
    cat(sprintf("  Treated with Cr_post: %d, pairs: %d\n",
                nrow(trt_valid), nrow(pairs_valid)))

    # ── Temporal alignment ─────────────────────────────────────────────
    tol_min <- cfg$tol_h * 60
    n_valid <- 0; n_pre_fail <- 0; n_post_fail <- 0; n_no_cr <- 0
    ctl_rows <- list()
    pre_gaps <- c(); post_gaps <- c()

    for (r in seq_len(nrow(pairs_valid))) {
      tpid <- pairs_valid$trt_pid[r]
      cpid <- pairs_valid$ctl_pid[r]
      tidx <- which(trt_valid$pid == tpid)[1]
      if (is.na(tidx)) next
      t_pre <- trt_valid$cr_pre_offset_min[tidx]
      t_post <- trt_valid$cr_post_offset_min[tidx]
      pair_id <- tidx

      ctl_cr <- cr_list[[as.character(cpid)]]
      if (is.null(ctl_cr) || nrow(ctl_cr) < 2) { n_no_cr <- n_no_cr+1; next }

      pd <- abs(ctl_cr$labresultoffset - t_pre)
      pi <- which.min(pd)
      if (pd[pi] > tol_min) { n_pre_fail <- n_pre_fail+1; next }

      pc <- ctl_cr[ctl_cr$labresultoffset > ctl_cr$labresultoffset[pi], ]
      if (nrow(pc) == 0) { n_post_fail <- n_post_fail+1; next }
      pod <- abs(pc$labresultoffset - t_post)
      poi <- which.min(pod)
      if (pod[poi] > tol_min) { n_post_fail <- n_post_fail+1; next }

      n_valid <- n_valid + 1
      pre_gaps <- c(pre_gaps, pd[pi]/60)
      post_gaps <- c(post_gaps, pod[poi]/60)
      ctl_rows[[n_valid]] <- data.frame(
        pid=cpid, cr_pre=ctl_cr$labresult[pi],
        cr_pre_offset_min=ctl_cr$labresultoffset[pi],
        cr_post=pc$labresult[poi],
        cr_post_offset_min=pc$labresultoffset[poi],
        delta_cr=pc$labresult[poi] - ctl_cr$labresult[pi],
        treated=0, match_pair_id=pair_id, stringsAsFactors=F)
    }

    if (n_valid == 0) {
      cat("  ⚠ No valid temporal matches\n"); next
    }

    ctl_matched <- do.call(rbind, ctl_rows)
    trt_with_match <- trt_valid[trt_valid$pid %in%
      trt_valid$pid[unique(ctl_matched$match_pair_id)], ]
    trt_with_match$match_pair_id <- match(trt_with_match$pid,
                                           trt_valid$pid[unique(ctl_matched$match_pair_id)])
    # Fix: use sequential pair IDs matching the ctl_matched
    valid_pairs <- sort(unique(ctl_matched$match_pair_id))
    trt_out <- trt_valid[valid_pairs, ]
    trt_out$match_pair_id <- seq_along(valid_pairs)
    ctl_matched$match_pair_id <- match(ctl_matched$match_pair_id, valid_pairs)

    matched_df <- rbind(trt_out, ctl_matched)

    cat(sprintf("  Temporal: %d valid / %d pairs (%.1f%%)\n",
                n_valid, nrow(pairs_valid), 100*n_valid/max(nrow(pairs_valid),1)))
    cat(sprintf("  Treated with ≥1 match: %d\n", nrow(trt_out)))
    cat(sprintf("  Pre gap median: %.1fh | Post gap median: %.1fh\n",
                median(pre_gaps), median(post_gaps)))

    # ── Merge covariates for doubly robust ─────────────────────────────
    covar_want <- c("pid", ps_vars, "surgery_type", "hosp_mortality")
    covar_trt <- trt[, intersect(covar_want, names(trt)), drop=F]
    covar_ctl <- ctl[, intersect(covar_want, names(ctl)), drop=F]
    shared <- intersect(names(covar_trt), names(covar_ctl))
    covar_all <- rbind(covar_trt[,shared], covar_ctl[,shared])
    covar_all <- covar_all[!duplicated(covar_all$pid), ]
    matched_df <- merge(matched_df, covar_all, by="pid", all.x=TRUE, suffixes=c("",".cov"))
    matched_df <- median_impute(matched_df, ps_vars)

    # ── Doubly robust DiD ──────────────────────────────────────────────
    cat("\n  ── Doubly robust DiD ──\n")
    res <- did_robust(matched_df, ps_vars)

    cat(sprintf("  Unadjusted:     DiD = %+.4f (SE=%.4f, P=%.4f)\n",
                res$did_unadj, res$se_unadj, res$p_unadj))
    cat(sprintf("  Doubly robust:  DiD = %+.4f (SE=%.4f, P=%.4f)\n",
                res$did_adj, res$se_adj, res$p_adj))
    cat(sprintf("                  95%% CI [%+.4f, %+.4f]\n", res$ci_lo, res$ci_hi))
    cat(sprintf("  Adjusted for %d covariates (SMD>0.05)\n", res$n_adjust))

    # Interpretation
    if (res$p_adj < 0.05) {
      dir <- if (res$did_adj < 0) "PROTECTIVE" else "HARMFUL"
      cat(sprintf("  → %s: IV Mg associated with %.3f mg/dL less Cr rise (P=%.4f)\n",
                  dir, abs(res$did_adj), res$p_adj))
    } else {
      cat(sprintf("  → Not significant at α=0.05\n"))
    }

    # Save matched dataset
    out_path <- file.path(RESULTS, sprintf("21_matched_%s_%s.csv", tag, cfg$label))
    write.csv(matched_df, out_path, row.names=FALSE)
    cat(sprintf("  Saved: %s (%d rows)\n", basename(out_path), nrow(matched_df)))

    # Store result
    all_results[[cfg$label]] <- data.frame(
      config = cfg$label,
      ratio = cfg$ratio,
      tolerance_h = cfg$tol_h,
      cr_strategy = cfg$cr_strategy,
      n_trt = res$n_trt,
      n_ctl = res$n_ctl,
      did_unadj = round(res$did_unadj, 4),
      p_unadj = round(res$p_unadj, 4),
      did_adj = round(res$did_adj, 4),
      se_adj = round(res$se_adj, 4),
      p_adj = round(res$p_adj, 4),
      ci_lo = round(res$ci_lo, 4),
      ci_hi = round(res$ci_hi, 4),
      n_adjust = res$n_adjust,
      stringsAsFactors = FALSE
    )
  }

  # ── Summary table ──────────────────────────────────────────────────────
  cat(sprintf("\n%s\n%s: RESULTS SUMMARY\n%s\n", SEP, db, SEP))
  res_df <- do.call(rbind, all_results)
  write.csv(res_df, file.path(RESULTS, sprintf("21_did_results_%s.csv", tag)), row.names=F)

  cat("\n  config         r  tol  cr_strategy    n_trt  n_ctl  DiD_unadj  DiD_adj    P_adj    95% CI\n")
  cat("  ───────────── ── ──── ───────────── ────── ────── ───────── ──────── ──────── ──────────────\n")
  for (i in seq_len(nrow(res_df))) {
    r <- res_df[i,]
    sig <- if (r$p_adj < 0.05) " *" else ""
    cat(sprintf("  %-14s %d  ±%dh  %-13s %5d  %5d   %+.4f  %+.4f   %.4f  [%+.4f,%+.4f]%s\n",
                r$config, r$ratio, r$tolerance_h, r$cr_strategy,
                r$n_trt, r$n_ctl, r$did_unadj, r$did_adj, r$p_adj,
                r$ci_lo, r$ci_hi, sig))
  }

  cat(sprintf("\n  Key: DiD < 0 = IV Mg protective (less Cr rise)\n"))
  cat(sprintf("  Doubly robust adjusts for covariates with SMD > 0.05 after matching\n"))
  cat(sprintf("  Cluster-robust SEs by match pair\n"))
  cat(sprintf("\n  Saved: 21_did_results_%s.csv\n", tag))

  return(res_df)
}

# ============================================================================
cat("======================================================================\n")
cat("21_did_final.R — DiD: PSM + temporal alignment + doubly robust\n")
cat("  Primary: r=1, ±6h, first Cr 6-24h\n")
cat("  Cr strategies: first_6_24h, peak_6_48h, closest_24h\n")
cat("======================================================================\n")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) { cat("Usage: Rscript 21_did_final.R eicu|mimic\n"); quit(status=1) }
for (a in args) run_did(toupper(a))

cat("\n======================================================================\n")
cat("Done. Review the summary table and matched CSVs.\n")
cat("Next: subgroup analyses (surgery type, Mg strata) on the primary matched set.\n")
cat("======================================================================\n")
