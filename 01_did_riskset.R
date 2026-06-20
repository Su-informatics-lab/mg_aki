#!/usr/bin/env Rscript
# ============================================================================
# 01_did_riskset.R — Risk-Set PSM with Temporal Alignment (v5 Design)
#
# For each treated patient (first IV Mg at t_mg):
#   1. Build risk set at t_mg (alive + in ICU + no Mg + no AKI + has Cr)
#   2. PSM 1:1 WITH REPLACEMENT (global PS, caliper 0.2 SD)
#   3. ΔCr = Cr(t_mg + 48h) - Cr_pre(t_mg), both groups
#   4. Pool with cluster SE on reused controls
#
# Usage: Rscript 01_did_riskset.R eicu
#        Rscript 01_did_riskset.R mimic
# ============================================================================

suppressPackageStartupMessages({ library(sandwich); library(lmtest); library(mice) })
source(path.expand("~/mg_aki/did_covars.R"))

RESULTS   <- path.expand("~/mg_aki/results")
PRIMARY_H <- 48       # primary outcome: ΔCr at t_mg + 48h
CR_WINDOW <- 12       # ±12h tolerance for Cr lookup
CALIPER_SD <- 0.2     # PS caliper in SD units
M_IMP     <- 20       # MICE imputations
TARGETS   <- c(6, 12, 18, 24, 30, 36, 42, 48)  # time course

# ── Helpers ────────────────────────────────────────────────────────────────

# Find Cr closest to target_h for a patient, from pre-sorted Cr list
find_cr <- function(cr_pt, target_h, window_h = CR_WINDOW) {
  if (is.null(cr_pt) || nrow(cr_pt) == 0) return(c(value = NA, time = NA))
  cand <- cr_pt[cr_pt$offset_h >= (target_h - window_h) &
                cr_pt$offset_h <= (target_h + window_h), ]
  if (nrow(cand) == 0) return(c(value = NA, time = NA))
  best <- cand[which.min(abs(cand$offset_h - target_h)), ]
  c(value = best$labresult, time = best$offset_h)
}

# Find last Cr before time t
find_cr_pre <- function(cr_pt, t_h) {
  if (is.null(cr_pt) || nrow(cr_pt) == 0) return(c(value = NA, time = NA))
  cand <- cr_pt[cr_pt$offset_h >= 0 & cr_pt$offset_h < t_h, ]
  if (nrow(cand) == 0) return(c(value = NA, time = NA))
  best <- cand[which.max(cand$offset_h), ]  # latest before t
  c(value = best$labresult, time = best$offset_h)
}

# Check if AKI exists at time t (relative to first Cr)
has_aki_at <- function(cr_pt, t_h) {
  if (is.null(cr_pt) || nrow(cr_pt) < 2) return(FALSE)
  crs <- cr_pt[cr_pt$offset_h >= 0 & cr_pt$offset_h <= t_h, ]
  if (nrow(crs) < 2) return(FALSE)
  baseline <- crs$labresult[1]
  if (is.na(baseline) || baseline <= 0) return(FALSE)
  any_aki <- any(crs$labresult[-1] - baseline >= 0.3, na.rm = TRUE) ||
             any(crs$labresult[-1] / baseline >= 1.5, na.rm = TRUE)
  return(any_aki)
}

safe_coeftest <- function(fit, varname = "treated") {
  ct <- tryCatch(suppressWarnings(coeftest(fit, vcov. = vcovCL(fit, cluster = fit$model$ctl_id))),
                 error = function(e) NULL)
  ok <- !is.null(ct) && is.matrix(ct) && varname %in% rownames(ct) &&
        ncol(ct) >= 4 && !any(is.nan(ct[varname, ]))
  if (ok) return(ct)
  ct <- tryCatch(suppressWarnings(coeftest(fit, vcov. = vcovHC(fit, type = "HC1"))),
                 error = function(e) NULL)
  ok <- !is.null(ct) && is.matrix(ct) && varname %in% rownames(ct) && ncol(ct) >= 4
  if (ok) return(ct)
  tryCatch(coeftest(fit), error = function(e) NULL)
}

# ============================================================================
# MAIN
# ============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) { cat("Usage: Rscript 01_did_riskset.R <db>\n"); quit(status = 1) }
db <- toupper(args[1]); tag <- tolower(db)

SEP <- paste(rep("=", 70), collapse = "")
cat(sprintf("\n%s\n01_did_riskset.R — Risk-Set PSM: %s\n%s\n", SEP, db, SEP))

# ── Load ──────────────────────────────────────────────────────────────────
all_pts <- read.csv(file.path(RESULTS, sprintf("did_all_%s.csv", tag)),
                    stringsAsFactors = FALSE)
cr_all  <- read.csv(file.path(RESULTS, sprintf("did_cr_all_%s.csv", tag)),
                    stringsAsFactors = FALSE)

# Normalize cr_all columns
cr_id <- if ("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
cr_all$pid <- cr_all[[cr_id]]
if (!"labresultoffset" %in% names(cr_all)) cr_all$labresultoffset <- cr_all$offset_min
cr_all$offset_h <- cr_all$labresultoffset / 60
cr_all <- cr_all[order(cr_all$pid, cr_all$offset_h), ]

N <- nrow(all_pts)
trt_idx  <- which(all_pts$treated == 1 & !is.na(all_pts$mg_offset_h))
n_trt    <- length(trt_idx)
cat(sprintf("  Loaded: %d patients (%d treated, %d control)\n",
            N, n_trt, N - n_trt))

# ── Pre-compute per-patient Cr indices ────────────────────────────────────
cat("  Pre-computing Cr indices...")
cr_list <- split(cr_all[, c("labresult", "offset_h")], cr_all$pid)

# earliest Cr time (for risk-set "has Cr before t_mg" check)
earliest_cr_h <- sapply(cr_list, function(x) min(x$offset_h, na.rm = TRUE))
all_pts$earliest_cr_h <- earliest_cr_h[as.character(all_pts$pid)]
all_pts$earliest_cr_h[is.na(all_pts$earliest_cr_h)] <- Inf

# first AKI time (for risk-set "no AKI at t_mg" check)
cat(" AKI onset...")
first_aki_h <- sapply(cr_list, function(cr_pt) {
  if (nrow(cr_pt) < 2) return(NA_real_)
  cr_pt <- cr_pt[order(cr_pt$offset_h), ]
  baseline <- cr_pt$labresult[1]
  if (is.na(baseline) || baseline <= 0) return(NA_real_)
  for (k in 2:nrow(cr_pt)) {
    delta <- cr_pt$labresult[k] - baseline
    ratio <- cr_pt$labresult[k] / baseline
    if (!is.na(delta) && (delta >= 0.3 || (!is.na(ratio) && ratio >= 1.5)))
      return(cr_pt$offset_h[k])
  }
  return(NA_real_)
})
all_pts$first_aki_h <- first_aki_h[as.character(all_pts$pid)]

# Exit time: when patient becomes ineligible as control
# = min(ICU discharge, first Mg, first AKI)
all_pts$exit_time <- pmin(
  all_pts$icu_discharge_h,
  ifelse(is.na(all_pts$mg_offset_h), Inf, all_pts$mg_offset_h),
  ifelse(is.na(all_pts$first_aki_h), Inf, all_pts$first_aki_h),
  na.rm = TRUE)

cat(sprintf(" done.\n"))
cat(sprintf("  Median exit_time: %.1fh (controls available until)\n",
            median(all_pts$exit_time[all_pts$treated == 0], na.rm = TRUE)))

# ── PS covariates ─────────────────────────────────────────────────────────
ps_vars <- intersect(PS_PRIMARY, names(all_pts))
ps_vars <- ps_vars[vapply(ps_vars, function(v) {
  x <- all_pts[[v]]; !all(is.na(x)) && var(x, na.rm = TRUE) > 1e-10
}, logical(1))]
cat(sprintf("  PS covariates: %d\n", length(ps_vars)))

# ── MICE ──────────────────────────────────────────────────────────────────
cat(sprintf("\n  ── MICE (m=%d) ──\n", M_IMP))
mice_vars <- intersect(ps_vars, names(all_pts))
to_impute <- mice_vars[vapply(mice_vars, function(v) any(is.na(all_pts[[v]])), logical(1))]
cat(sprintf("  Variables to impute: %s\n", paste(to_impute, collapse = ", ")))

if (length(to_impute) > 0) {
  mice_df <- all_pts[, c("treated", "mg_offset_h", mice_vars)]
  mice_df$mg_offset_h[is.na(mice_df$mg_offset_h)] <- 0  # controls = 0 for imputation
  meth <- rep("", ncol(mice_df))
  names(meth) <- names(mice_df)
  for (v in to_impute) meth[v] <- "pmm"
  imp <- mice(mice_df, m = M_IMP, method = meth, printFlag = FALSE, maxit = 10)
  cat(sprintf("  MICE done. Logged events: %d\n", nrow(imp$loggedEvents)))
} else {
  cat("  No missing data — skip MICE\n")
  imp <- NULL
}

# ============================================================================
# MATCHING + DiD — per imputation
# ============================================================================

run_one_imp <- function(m_idx) {
  cat(sprintf("\n  ── Imputation m=%d ──\n", m_idx))

  # Get imputed data
  if (!is.null(imp)) {
    imp_df <- complete(imp, m_idx)
    for (v in to_impute) all_pts[[v]] <- imp_df[[v]]
  }

  # Fit global PS
  ps_fml <- as.formula(paste("treated ~", paste(ps_vars, collapse = "+")))
  ps_fit <- glm(ps_fml, data = all_pts, family = binomial())
  all_pts$ps <- predict(ps_fit, type = "response")
  caliper <- CALIPER_SD * sd(all_pts$ps, na.rm = TRUE)

  # ── Risk-set matching ───────────────────────────────────────────
  # Sort treated by t_mg for reporting
  trt_sorted <- trt_idx[order(all_pts$mg_offset_h[trt_idx])]

  matches <- data.frame(
    trt_pid = integer(0), ctl_pid = integer(0), t_mg = numeric(0),
    ps_dist = numeric(0), risk_set_size = integer(0),
    stringsAsFactors = FALSE)

  n_matched <- 0; n_failed <- 0

  for (i in trt_sorted) {
    t_mg <- all_pts$mg_offset_h[i]
    ps_i <- all_pts$ps[i]

    # Risk set (vectorized)
    eligible <- (all_pts$exit_time > t_mg) &
                (all_pts$earliest_cr_h <= t_mg) &
                (seq_len(N) != i)
    rs_size <- sum(eligible, na.rm = TRUE)
    if (rs_size == 0) { n_failed <- n_failed + 1; next }

    # PS matching within risk set
    eligible_idx <- which(eligible)
    ps_dist <- abs(all_pts$ps[eligible_idx] - ps_i)
    within_cal <- ps_dist <= caliper
    if (sum(within_cal) == 0) { n_failed <- n_failed + 1; next }

    # Best match (nearest neighbor); tie-break: earliest Cr before t_mg
    cal_idx <- eligible_idx[within_cal]
    cal_dist <- ps_dist[within_cal]
    best_pos <- which.min(cal_dist)
    j <- cal_idx[best_pos]

    n_matched <- n_matched + 1
    matches[n_matched, ] <- list(
      trt_pid = all_pts$pid[i], ctl_pid = all_pts$pid[j],
      t_mg = t_mg, ps_dist = cal_dist[best_pos], risk_set_size = rs_size)
  }
  cat(sprintf("    Matched: %d/%d (%.1f%%), failed: %d\n",
              n_matched, n_trt, 100 * n_matched / n_trt, n_failed))

  if (n_matched < 50) {
    cat("    Too few matches — skipping this imputation\n")
    return(NULL)
  }

  # ── Compute outcomes ──────────────────────────────────────────────
  # For each matched pair, compute ΔCr at multiple time points
  tc_list <- list()
  for (target_h in TARGETS) {
    dcr_trt <- dcr_ctl <- cr_pre_diff <- numeric(n_matched)
    valid <- logical(n_matched)

    for (k in 1:n_matched) {
      tpid <- as.character(matches$trt_pid[k])
      cpid <- as.character(matches$ctl_pid[k])
      t_mg <- matches$t_mg[k]

      # Cr_pre: last Cr before t_mg
      pre_t <- find_cr_pre(cr_list[[tpid]], t_mg)
      pre_c <- find_cr_pre(cr_list[[cpid]], t_mg)

      # Cr_post: Cr closest to t_mg + target_h
      post_t <- find_cr(cr_list[[tpid]], t_mg + target_h)
      post_c <- find_cr(cr_list[[cpid]], t_mg + target_h)

      if (any(is.na(c(pre_t[1], pre_c[1], post_t[1], post_c[1])))) {
        valid[k] <- FALSE; next
      }

      dcr_trt[k] <- post_t[1] - pre_t[1]
      dcr_ctl[k] <- post_c[1] - pre_c[1]
      cr_pre_diff[k] <- abs(pre_t[2] - pre_c[2])  # time alignment quality
      valid[k] <- TRUE
    }

    n_valid <- sum(valid)
    if (n_valid < 30) {
      tc_list[[as.character(target_h)]] <- list(
        target_h = target_h, n = n_valid, did = NA, se = NA, p = NA)
      next
    }

    # DiD with cluster SE (control pid as cluster for with-replacement)
    pair_df <- data.frame(
      delta_cr = c(dcr_trt[valid], dcr_ctl[valid]),
      treated = rep(c(1, 0), each = n_valid),
      ctl_id = rep(matches$ctl_pid[valid], 2),
      stringsAsFactors = FALSE)

    fit <- lm(delta_cr ~ treated, data = pair_df)
    ct <- safe_coeftest(fit)
    if (!is.null(ct) && "treated" %in% rownames(ct)) {
      est <- ct["treated", 1]; se <- ct["treated", 2]; p <- ct["treated", ncol(ct)]
    } else {
      est <- coef(fit)["treated"]; se <- NA; p <- NA
    }

    tc_list[[as.character(target_h)]] <- list(
      target_h = target_h, n = n_valid,
      did = est, se = se, p = p,
      cr_pre_align_median = median(cr_pre_diff[valid]),
      mean_trt = mean(dcr_trt[valid]), mean_ctl = mean(dcr_ctl[valid]))
  }

  # ── Quality metrics ───────────────────────────────────────────────
  reuse <- table(matches$ctl_pid)
  unique_ctl <- length(reuse)
  max_reuse <- max(reuse)

  list(matches = matches, tc = tc_list,
       n_matched = n_matched, n_failed = n_failed,
       unique_ctl = unique_ctl, max_reuse = max_reuse)
}

# ============================================================================
# RUN ALL IMPUTATIONS + POOL
# ============================================================================

cat(sprintf("\n%s\n  RISK-SET MATCHING (m=%d imputations)\n%s\n", SEP, M_IMP, SEP))

all_results <- list()
for (m in 1:M_IMP) {
  all_results[[m]] <- run_one_imp(m)
}

# Pool with Rubin's rules
cat(sprintf("\n%s\n  POOLED RESULTS (Rubin's rules, %d imputations)\n%s\n", SEP, M_IMP, SEP))

cat("  target_h   DiD        SE       P        95%% CI                n    FMI\n")
cat("  ────────   ────────   ──────   ──────   ────────────────────   ───  ─────\n")

pooled_rows <- list()
for (target_h in TARGETS) {
  th_key <- as.character(target_h)
  dids <- ses <- numeric(0)

  for (m in 1:M_IMP) {
    r <- all_results[[m]]
    if (is.null(r) || is.null(r$tc[[th_key]])) next
    tc <- r$tc[[th_key]]
    if (!is.na(tc$did) && !is.na(tc$se)) {
      dids <- c(dids, tc$did)
      ses  <- c(ses, tc$se)
    }
  }

  if (length(dids) < 2) {
    cat(sprintf("  %4dh      (insufficient imputations)\n", target_h))
    next
  }

  # Rubin's rules
  m_valid <- length(dids)
  Q_bar <- mean(dids)
  U_bar <- mean(ses^2)
  B     <- var(dids)
  T_var <- U_bar + (1 + 1/m_valid) * B
  T_se  <- sqrt(T_var)

  # Barnard-Rubin df
  lambda <- ((1 + 1/m_valid) * B) / T_var
  df_old <- (m_valid - 1) / lambda^2
  n_avg  <- mean(sapply(all_results, function(r) {
    if (is.null(r$tc[[th_key]])) return(NA)
    r$tc[[th_key]]$n
  }), na.rm = TRUE)
  df_obs <- (n_avg - 1) * (1 - lambda) / (1 + lambda)  # simplified
  df     <- max(3, min(df_old, df_obs, na.rm = TRUE))

  p_val  <- 2 * pt(abs(Q_bar / T_se), df = df, lower.tail = FALSE)
  ci_lo  <- Q_bar - qt(0.975, df) * T_se
  ci_hi  <- Q_bar + qt(0.975, df) * T_se
  fmi    <- lambda

  is_primary <- (target_h == PRIMARY_H)
  sig <- if (!is.na(p_val) && p_val < 0.05) " *" else "  "
  tag_str <- if (is_primary) " << PRIMARY" else ""

  cat(sprintf("  %4dh      %+.4f   %.4f   %.4f%s [%+.4f,%+.4f]   %3.0f  %.3f%s\n",
              target_h, Q_bar, T_se, p_val, sig, ci_lo, ci_hi, n_avg, fmi, tag_str))

  pooled_rows[[length(pooled_rows) + 1]] <- data.frame(
    target_h = target_h, did = Q_bar, se = T_se, p = p_val,
    ci_lo = ci_lo, ci_hi = ci_hi, fmi = fmi,
    n_avg = n_avg, m_valid = m_valid,
    stringsAsFactors = FALSE)
}

if (length(pooled_rows) > 0) {
  pooled_df <- do.call(rbind, pooled_rows)
  write.csv(pooled_df, file.path(RESULTS, sprintf("did_riskset_%s.csv", tag)),
            row.names = FALSE)
}

# ── Quality report ──────────────────────────────────────────────────────
cat(sprintf("\n%s\n  QUALITY METRICS (from m=1)\n%s\n", SEP, SEP))

r1 <- all_results[[1]]
if (!is.null(r1)) {
  cat(sprintf("  Matched pairs: %d / %d treated (%.1f%%)\n",
              r1$n_matched, n_trt, 100 * r1$n_matched / n_trt))
  cat(sprintf("  Failed matches: %d\n", r1$n_failed))
  cat(sprintf("  Unique controls used: %d (of %d available)\n",
              r1$unique_ctl, N - n_trt))
  cat(sprintf("  Max control reuse: %dx\n", r1$max_reuse))

  reuse <- table(r1$matches$ctl_pid)
  cat(sprintf("  Control reuse: 1x=%d, 2x=%d, 3x=%d, 4+x=%d\n",
              sum(reuse == 1), sum(reuse == 2), sum(reuse == 3), sum(reuse >= 4)))

  cat(sprintf("\n  Risk set size: median=%.0f, IQR=[%.0f, %.0f], min=%d\n",
              median(r1$matches$risk_set_size),
              quantile(r1$matches$risk_set_size, 0.25),
              quantile(r1$matches$risk_set_size, 0.75),
              min(r1$matches$risk_set_size)))

  cat(sprintf("  PS distance: median=%.4f, max=%.4f\n",
              median(r1$matches$ps_dist), max(r1$matches$ps_dist)))

  cat(sprintf("  t_mg distribution (matched): median=%.1fh, IQR=[%.1f, %.1fh]\n",
              median(r1$matches$t_mg),
              quantile(r1$matches$t_mg, 0.25),
              quantile(r1$matches$t_mg, 0.75)))

  # Save matched pairs + quality metrics
  write.csv(r1$matches,
            file.path(RESULTS, sprintf("did_riskset_pairs_%s.csv", tag)),
            row.names = FALSE)

  # Cr timing alignment at primary target
  tc_primary <- r1$tc[[as.character(PRIMARY_H)]]
  if (!is.null(tc_primary) && !is.na(tc_primary$cr_pre_align_median))
    cat(sprintf("  Cr_pre alignment: median |t_diff| = %.1fh\n",
                tc_primary$cr_pre_align_median))
}

cat(sprintf("\n%s\n01_did_riskset.R — %s COMPLETE\n%s\n", SEP, db, SEP))

# Output manifest
cat("\n  Files written:\n")
for (f in c(sprintf("did_riskset_%s.csv", tag),
            sprintf("did_riskset_pairs_%s.csv", tag))) {
  fp <- file.path(RESULTS, f)
  if (file.exists(fp)) cat(sprintf("    %s\n", f))
}
