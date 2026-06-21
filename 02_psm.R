#!/usr/bin/env Rscript
# ============================================================================
# 02_psm.R — Risk-Set PSM with Dual Control Pool (v5 Design)
#
# Primary:     yet-untreated controls (sequential trial: "Mg now vs not now")
# Sensitivity: never-treated controls (parallel trial:   "Mg vs never Mg")
#
# Both share a single MICE imputation + global PS model.
# Each gets its own risk sets → matching → outcomes → balance.
#
# Lab timing: LAST (closest to T₀)
# PS: 21 covariates (no drug flags)
# MICE: m=20, averaged → single PS → single match
# Matching: 1:1 with replacement, caliper 0.2 SD, HC1 SE
#
# Usage: Rscript 02_psm.R eicu
#        Rscript 02_psm.R mimic
# ============================================================================

suppressPackageStartupMessages({ library(sandwich); library(lmtest); library(mice) })

RESULTS    <- path.expand("~/mg_aki/results")
PRIMARY_H  <- 24
CR_WINDOW  <- 12
CALIPER_SD <- 0.2
M_IMP      <- 20
TARGETS    <- c(6, 12, 18, 24, 30, 36)

PS_FIXED <- c("age","is_female","bmi",
              "surg_cabg","surg_valve","surg_combined",
              "heart_failure","hypertension","diabetes","ckd",
              "copd","pvd","stroke","liver_disease","egfr")
PS_LAB_BASES <- c("magnesium","potassium","calcium","lactate","heartrate")

# ═══════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════

find_cr <- function(cr_pt, target_h, window = CR_WINDOW) {
  if (is.null(cr_pt) || nrow(cr_pt) == 0) return(c(NA, NA))
  cand <- cr_pt[cr_pt$offset_h >= (target_h - window) &
                cr_pt$offset_h <= (target_h + window), ]
  if (nrow(cand) == 0) return(c(NA, NA))
  best <- cand[which.min(abs(cand$offset_h - target_h)), ]
  c(best$labresult, best$offset_h)
}

find_cr_pre <- function(cr_pt, t_h) {
  if (is.null(cr_pt) || nrow(cr_pt) == 0) return(c(NA, NA))
  cand <- cr_pt[cr_pt$offset_h >= 0 & cr_pt$offset_h < t_h, ]
  if (nrow(cand) == 0) return(c(NA, NA))
  best <- cand[which.max(cand$offset_h), ]
  c(best$labresult, best$offset_h)
}

safe_coeftest <- function(fit) {
  ct <- tryCatch(suppressWarnings(coeftest(fit, vcov. = vcovHC(fit, type = "HC1"))),
                 error = function(e) NULL)
  if (!is.null(ct) && is.matrix(ct) && "treated" %in% rownames(ct) &&
      ncol(ct) >= 4 && !any(is.nan(ct["treated", ]))) return(ct)
  tryCatch(coeftest(fit), error = function(e) NULL)
}

compute_balance <- function(all_pts, ps_vars, trt_m, ctl_m, label) {
  trt_all <- which(all_pts$treated == 1)
  ctl_all <- which(all_pts$treated == 0)
  n_viol <- 0; max_smd <- 0
  cat(sprintf("\n  ── COVARIATE BALANCE [%s] ──\n", label))
  cat("  Covariate                raw_SMD  matched_SMD  status\n")
  cat("  ───────────────────────  ───────  ───────────  ──────\n")
  for (v in ps_vars) {
    x1_raw <- all_pts[[v]][trt_all]; x0_raw <- all_pts[[v]][ctl_all]
    x1_mat <- all_pts[[v]][trt_m];   x0_mat <- all_pts[[v]][ctl_m]
    sp_raw <- sqrt((var(x1_raw, na.rm=T) + var(x0_raw, na.rm=T)) / 2)
    sp_mat <- sqrt((var(x1_mat, na.rm=T) + var(x0_mat, na.rm=T)) / 2)
    smd_raw <- if (!is.na(sp_raw) && sp_raw > 1e-10)
                 abs(mean(x1_raw,na.rm=T) - mean(x0_raw,na.rm=T)) / sp_raw else NA
    smd_mat <- if (!is.na(sp_mat) && sp_mat > 1e-10)
                 abs(mean(x1_mat,na.rm=T) - mean(x0_mat,na.rm=T)) / sp_mat else NA
    if (!is.na(smd_mat) && smd_mat > max_smd) max_smd <- smd_mat
    flag <- if (!is.na(smd_mat) && smd_mat > 0.1) { n_viol <- n_viol+1; "VIOL" } else "ok"
    cat(sprintf("  %-25s  %.3f    %.3f       %s\n", v,
                ifelse(is.na(smd_raw), NA, smd_raw),
                ifelse(is.na(smd_mat), NA, smd_mat), flag))
  }
  cat(sprintf("\n  Max matched SMD: %.3f | Violations (>0.1): %d/%d\n",
              max_smd, n_viol, length(ps_vars)))
  list(max_smd = max_smd, n_viol = n_viol)
}

# ═══════════════════════════════════════════════════════════════════
# run_pool — builds risk sets, matches, computes outcomes for one
#            control-pool definition
# ═══════════════════════════════════════════════════════════════════

run_pool <- function(pool_name, all_pts, trt_idx, cr_list, ps_vars, caliper) {

  n_trt   <- length(trt_idx)
  trt_pids <- all_pts$pid[trt_idx]
  trt_tmg  <- all_pts$mg_offset_h[trt_idx]

  sep <- paste(rep("=", 70), collapse = "")
  cat(sprintf("\n%s\n  POOL: %s  (n_trt = %d)\n%s\n", sep, toupper(pool_name), n_trt, sep))

  # ── Risk sets ─────────────────────────────────────────────────
  cat(sprintf("  Building %s risk sets...\n", pool_name))
  risk_sets <- vector("list", n_trt)

  for (k in seq_len(n_trt)) {
    t_mg <- trt_tmg[k]

    if (pool_name == "yet_untreated") {
      # Sequential trial: controls = anyone not yet treated at t_mg,
      # and Mg-free through t_mg + PRIMARY_H (no contamination at 24h)
      risk_sets[[k]] <- which(
        all_pts$icu_discharge_h > t_mg &
        (is.na(all_pts$first_aki_h) | all_pts$first_aki_h > t_mg) &
        all_pts$earliest_cr_h <= t_mg &
        (is.na(all_pts$mg_offset_h) | all_pts$mg_offset_h > t_mg + PRIMARY_H) &
        all_pts$pid != trt_pids[k])

    } else if (pool_name == "never_treated") {
      # Parallel trial: controls = never received IV Mg at all
      risk_sets[[k]] <- which(
        all_pts$treated == 0 &
        all_pts$icu_discharge_h > t_mg &
        (is.na(all_pts$first_aki_h) | all_pts$first_aki_h > t_mg) &
        all_pts$earliest_cr_h <= t_mg &
        all_pts$pid != trt_pids[k])
    }
  }

  rs_sizes <- vapply(risk_sets, length, integer(1))
  cat(sprintf("    Risk sets: median=%.0f, IQR=[%.0f,%.0f], empty=%d\n",
              median(rs_sizes), quantile(rs_sizes, 0.25),
              quantile(rs_sizes, 0.75), sum(rs_sizes == 0)))

  # ── Match ───────────────────────────────────────────────────────
  cat("  Matching...\n")
  match_trt <- integer(n_trt); match_ctl <- integer(n_trt)
  match_ps  <- numeric(n_trt); matched <- logical(n_trt)

  for (k in seq_len(n_trt)) {
    rs <- risk_sets[[k]]
    if (length(rs) == 0) next
    ps_i <- all_pts$ps[trt_idx[k]]
    ps_dist <- abs(all_pts$ps[rs] - ps_i)
    within_cal <- which(ps_dist <= caliper)
    if (length(within_cal) == 0) next
    best <- within_cal[which.min(ps_dist[within_cal])]
    match_trt[k] <- trt_idx[k]; match_ctl[k] <- rs[best]
    match_ps[k] <- ps_dist[best]; matched[k] <- TRUE
  }
  n_matched <- sum(matched)
  cat(sprintf("  Matched: %d/%d (%.1f%%)\n", n_matched, n_trt,
              100*n_matched/n_trt))

  if (n_matched < 50) {
    cat("  WARNING: <50 matches, skipping outcome computation\n")
    return(NULL)
  }

  # ── Temporal diagnostics ────────────────────────────────────────
  cat(sprintf("\n  ── TEMPORAL DIAGNOSTICS [%s] ──\n", pool_name))
  cat(sprintf("    Matched treated T₀: median=%.1fh, IQR=[%.1f,%.1f]\n",
              median(trt_tmg[matched], na.rm=TRUE),
              quantile(trt_tmg[matched], 0.25, na.rm=TRUE),
              quantile(trt_tmg[matched], 0.75, na.rm=TRUE)))

  ctl_indices <- match_ctl[matched]
  ctl_is_never <- all_pts$treated[ctl_indices] == 0
  ctl_is_future <- all_pts$treated[ctl_indices] == 1
  n_never <- sum(ctl_is_never, na.rm=TRUE)
  n_future <- sum(ctl_is_future, na.rm=TRUE)
  cat(sprintf("    Matched controls: %d never-treated (%.1f%%), %d yet-untreated (%.1f%%)\n",
              n_never, 100*n_never/n_matched,
              n_future, 100*n_future/n_matched))

  if (n_future > 0) {
    ctl_tmg <- all_pts$mg_offset_h[ctl_indices[ctl_is_future]]
    trt_tmg_m <- trt_tmg[matched][ctl_is_future]
    gap <- ctl_tmg - trt_tmg_m
    cat(sprintf("    Yet-untreated Mg gap: median=%.1fh later, IQR=[%.1f,%.1f]\n",
                median(gap, na.rm=TRUE), quantile(gap, 0.25, na.rm=TRUE),
                quantile(gap, 0.75, na.rm=TRUE)))
  }

  n_unique_ctl <- length(unique(ctl_indices))
  cat(sprintf("    Unique controls: %d (reuse: %.1fx)\n",
              n_unique_ctl, n_matched / n_unique_ctl))

  # Reuse distribution
  reuse_tab <- table(table(ctl_indices))
  cat("    Reuse distribution: ")
  for (nm in names(reuse_tab)) {
    cat(sprintf("%sx=%d ", nm, reuse_tab[nm]))
  }
  cat("\n")

  # ── Outcomes ────────────────────────────────────────────────────
  cat("  Computing outcomes...\n")
  results <- list()

  for (target_h in TARGETS) {
    dcr_trt <- dcr_ctl <- numeric(n_matched)
    valid <- logical(n_matched); idx <- 0

    for (kk in which(matched)) {
      idx <- idx + 1
      tpid <- as.character(all_pts$pid[match_trt[kk]])
      cpid <- as.character(all_pts$pid[match_ctl[kk]])
      t_mg <- trt_tmg[kk]

      # Contamination check: for yet-untreated pool, controls who
      # receive Mg before this target time are excluded.
      # For never-treated pool this is a no-op (controls never get Mg).
      if (pool_name == "yet_untreated") {
        ctl_mg_h <- all_pts$mg_offset_h[match_ctl[kk]]
        if (!is.na(ctl_mg_h) && ctl_mg_h < t_mg + target_h) {
          valid[idx] <- FALSE; next
        }
      }

      pre_t  <- find_cr_pre(cr_list[[tpid]], t_mg)
      pre_c  <- find_cr_pre(cr_list[[cpid]], t_mg)
      post_t <- find_cr(cr_list[[tpid]], t_mg + target_h)
      post_c <- find_cr(cr_list[[cpid]], t_mg + target_h)

      if (any(is.na(c(pre_t[1], pre_c[1], post_t[1], post_c[1])))) {
        valid[idx] <- FALSE; next
      }
      dcr_trt[idx] <- post_t[1] - pre_t[1]
      dcr_ctl[idx] <- post_c[1] - pre_c[1]
      valid[idx] <- TRUE
    }

    n_valid <- sum(valid)
    if (n_valid < 30) {
      results[[length(results)+1]] <- data.frame(
        pool=pool_name, target_h=target_h, n=n_valid,
        did=NA, se=NA, p=NA, ci_lo=NA, ci_hi=NA)
      next
    }

    pair_df <- data.frame(
      delta_cr = c(dcr_trt[valid], dcr_ctl[valid]),
      treated = rep(c(1,0), each = n_valid))
    fit <- lm(delta_cr ~ treated, data = pair_df)
    ct <- safe_coeftest(fit)
    est <- if (!is.null(ct) && "treated" %in% rownames(ct)) ct["treated",1] else coef(fit)["treated"]
    se  <- if (!is.null(ct) && "treated" %in% rownames(ct)) ct["treated",2] else NA
    pv  <- if (!is.null(ct) && "treated" %in% rownames(ct)) ct["treated",4] else NA
    ci_lo <- est - 1.96 * se; ci_hi <- est + 1.96 * se

    results[[length(results)+1]] <- data.frame(
      pool=pool_name, target_h=target_h, n=n_valid,
      did=est, se=se, p=pv, ci_lo=ci_lo, ci_hi=ci_hi)
  }

  res <- do.call(rbind, results)

  # ── Print results ───────────────────────────────────────────────
  cat(sprintf("\n  ── RESULTS [%s] ──\n", toupper(pool_name)))
  cat("  target   DiD        SE       P        95% CI                n\n")
  cat("  ──────   ────────   ──────   ──────   ────────────────────  ────\n")
  for (i in 1:nrow(res)) {
    r <- res[i, ]
    if (is.na(r$did)) {
      cat(sprintf("  %4dh    (insufficient data, n=%d)\n", r$target_h, r$n))
      next
    }
    sig <- if (!is.na(r$p) && r$p < 0.05) " *" else "  "
    pri <- if (r$target_h == PRIMARY_H) " << PRIMARY" else ""
    cat(sprintf("  %4dh    %+.4f   %.4f   %.4f%s [%+.4f,%+.4f]  %4d%s\n",
                r$target_h, r$did, r$se, r$p, sig,
                r$ci_lo, r$ci_hi, r$n, pri))
  }

  # ── Balance ─────────────────────────────────────────────────────
  trt_m <- match_trt[matched]; ctl_m <- match_ctl[matched]
  bal <- compute_balance(all_pts, ps_vars, trt_m, ctl_m, pool_name)

  # ── Save matched pairs for downstream HTE ───────────────────────
  pairs_df <- data.frame(
    match_pair_id = seq_len(n_matched),
    trt_row = match_trt[matched],
    ctl_row = match_ctl[matched],
    trt_pid = all_pts$pid[match_trt[matched]],
    ctl_pid = all_pts$pid[match_ctl[matched]],
    t_mg = trt_tmg[matched]
  )

  list(results = res, balance = bal, pairs = pairs_df, n_matched = n_matched)
}


# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) { cat("Usage: Rscript 02_psm.R <db>\n"); quit(status = 1) }
db <- toupper(args[1]); tag <- tolower(db)

SEP <- paste(rep("=", 70), collapse = "")
cat(sprintf("\n%s\n02_psm.R — Risk-Set PSM: %s\n", SEP, db))
cat(sprintf("  Primary: yet-untreated (sequential trial estimand)\n"))
cat(sprintf("  Sensitivity: never-treated (parallel trial estimand)\n"))
cat(sprintf("  PS: 21 covars | Labs: LAST | MICE m=%d averaged\n%s\n", M_IMP, SEP))

# ── Load data ─────────────────────────────────────────────────────
all_pts <- read.csv(file.path(RESULTS, sprintf("did_all_%s.csv", tag)),
                    stringsAsFactors = FALSE)
cr_all  <- read.csv(file.path(RESULTS, sprintf("did_cr_all_%s.csv", tag)),
                    stringsAsFactors = FALSE)

cr_id <- if ("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
cr_all$pid <- cr_all[[cr_id]]
if (!"offset_h" %in% names(cr_all)) cr_all$offset_h <- cr_all$labresultoffset / 60
cr_all <- cr_all[order(cr_all$pid, cr_all$offset_h), ]

cat("  Loading labs...\n")
labs_raw <- read.csv(file.path(RESULTS, sprintf("did_labs_all_%s.csv", tag)),
                     stringsAsFactors = FALSE)
pid_col_lab <- if ("patientunitstayid" %in% names(labs_raw)) "patientunitstayid" else "stay_id"
if (pid_col_lab %in% names(labs_raw)) labs_raw$pid <- labs_raw[[pid_col_lab]]
labs_elec <- labs_raw[labs_raw$lab_name %in% c("magnesium","potassium","calcium","lactate"), ]
labs_hr   <- labs_raw[labs_raw$lab_name == "heartrate", ]
if (nrow(labs_hr) > 500000) {
  labs_hr$hour_bin <- floor(labs_hr$offset_h)
  labs_hr <- labs_hr[order(labs_hr$pid, labs_hr$offset_h), ]
  labs_hr <- labs_hr[!duplicated(paste(labs_hr$pid, labs_hr$hour_bin)), ]
  labs_hr$hour_bin <- NULL
  cat(sprintf("    HR downsampled: %d → %d\n",
              sum(labs_raw$lab_name == "heartrate"), nrow(labs_hr)))
}
labs <- rbind(labs_elec, labs_hr)
rm(labs_raw, labs_elec, labs_hr); gc()

N <- nrow(all_pts)
cat(sprintf("  Patients: %d (%d treated, %d control)\n", N,
            sum(all_pts$treated == 1, na.rm=TRUE),
            sum(all_pts$treated == 0, na.rm=TRUE)))

# ── LAST lab values (closest to T₀) ──────────────────────────────
cat("  Computing LAST lab values...\n")
for (ln in PS_LAB_BASES) {
  sub <- labs[labs$lab_name == ln, ]
  if (nrow(sub) == 0) next
  sub$mg_offset_h <- all_pts$mg_offset_h[match(sub$pid, all_pts$pid)]
  sub <- sub[sub$offset_h >= 0 &
             (is.na(sub$mg_offset_h) | sub$offset_h < sub$mg_offset_h), ]
  if (nrow(sub) == 0) next
  s <- sub[order(-sub$offset_h), ]  # descending: LAST value first
  s <- s[!duplicated(s$pid), ]
  idx <- match(all_pts$pid, s$pid)
  all_pts[[paste0("first_", ln)]] <- s$value[idx]
  nf <- sum(!is.na(all_pts[[paste0("first_", ln)]]))
  cat(sprintf("    %s: %d (%.0f%%)\n", ln, nf, 100*nf/N))
}
all_pts$first_lactate_missing <- as.integer(is.na(all_pts$first_lactate))

# ── Treatment indices ─────────────────────────────────────────────
trt_idx <- which(all_pts$treated == 1 & !is.na(all_pts$mg_offset_h))
n_trt <- length(trt_idx)
cat(sprintf("  Treated eligible for matching: %d\n", n_trt))

# ── Cr list + AKI pre-computation (shared) ────────────────────────
cat("  Pre-computing Cr lists and AKI times...\n")
cr_list <- split(cr_all[, c("labresult","offset_h")], cr_all$pid)

earliest_cr <- sapply(cr_list, function(x) min(x$offset_h, na.rm = TRUE))
all_pts$earliest_cr_h <- earliest_cr[as.character(all_pts$pid)]
all_pts$earliest_cr_h[is.na(all_pts$earliest_cr_h)] <- Inf

first_aki <- sapply(cr_list, function(cr_pt) {
  if (nrow(cr_pt) < 2) return(NA_real_)
  cr_pt <- cr_pt[order(cr_pt$offset_h), ]
  bl <- cr_pt$labresult[1]
  if (is.na(bl) || bl <= 0) return(NA_real_)
  for (k in 2:nrow(cr_pt)) {
    d <- cr_pt$labresult[k] - bl
    r <- cr_pt$labresult[k] / bl
    if (!is.na(d) && (d >= 0.3 || (!is.na(r) && r >= 1.5)))
      return(cr_pt$offset_h[k])
  }
  NA_real_
})
all_pts$first_aki_h <- first_aki[as.character(all_pts$pid)]

# ── PS covariates ─────────────────────────────────────────────────
lab_cols <- paste0("first_", PS_LAB_BASES)
ps_vars <- c(intersect(PS_FIXED, names(all_pts)),
             intersect(c(lab_cols, "first_lactate_missing"), names(all_pts)))
ps_vars <- ps_vars[vapply(ps_vars, function(v) {
  x <- all_pts[[v]]; !all(is.na(x)) && var(x, na.rm = TRUE) > 1e-10
}, logical(1))]
cat(sprintf("  PS covariates (%d): %s\n", length(ps_vars),
            paste(ps_vars, collapse = ", ")))

# ── MICE m=20, average imputed values ─────────────────────────────
to_impute <- ps_vars[vapply(ps_vars, function(v)
  any(is.na(all_pts[[v]])), logical(1))]
cat(sprintf("  MICE m=%d on %d vars: %s\n", M_IMP, length(to_impute),
            paste(to_impute, collapse = ", ")))

if (length(to_impute) > 0) {
  mice_df <- all_pts[, c("treated", ps_vars)]
  meth <- rep("", ncol(mice_df)); names(meth) <- names(mice_df)
  for (v in to_impute) meth[v] <- "pmm"
  imp <- mice(mice_df, m = M_IMP, method = meth, printFlag = FALSE, maxit = 10)
  cat(sprintf("  MICE done. Logged events: %d\n", nrow(imp$loggedEvents)))
  for (v in to_impute) {
    vals <- sapply(1:M_IMP, function(m) complete(imp, m)[[v]])
    all_pts[[v]] <- rowMeans(vals, na.rm = TRUE)
  }
  cat("  Imputed values averaged across m=20\n")
} else {
  cat("  No imputation needed\n")
}

# ── Fit global PS (ONCE) ──────────────────────────────────────────
cat("  Fitting global PS...\n")
ps_fml <- as.formula(paste("treated ~", paste(ps_vars, collapse = "+")))
ps_fit <- glm(ps_fml, data = all_pts, family = binomial())
all_pts$ps <- predict(ps_fit, type = "response")
caliper <- CALIPER_SD * sd(all_pts$ps, na.rm = TRUE)
cat(sprintf("  PS: AUC-proxy range [%.3f, %.3f], caliper=%.4f\n",
            min(all_pts$ps, na.rm=T), max(all_pts$ps, na.rm=T), caliper))

# ═══════════════════════════════════════════════════════════════════
# RUN BOTH POOLS
# ═══════════════════════════════════════════════════════════════════

out_yt <- run_pool("yet_untreated", all_pts, trt_idx, cr_list, ps_vars, caliper)
out_nt <- run_pool("never_treated", all_pts, trt_idx, cr_list, ps_vars, caliper)

# ═══════════════════════════════════════════════════════════════════
# COMPARISON TABLE
# ═══════════════════════════════════════════════════════════════════

cat(sprintf("\n%s\n  HEAD-TO-HEAD COMPARISON: %s\n%s\n", SEP, db, SEP))
cat("  target   yet-untreated                     never-treated\n")
cat("           DiD       P      n                DiD       P      n\n")
cat("  ──────   ────────  ─────  ────             ────────  ─────  ────\n")

for (h in TARGETS) {
  yt <- if (!is.null(out_yt)) out_yt$results[out_yt$results$target_h == h, ] else NULL
  nt <- if (!is.null(out_nt)) out_nt$results[out_nt$results$target_h == h, ] else NULL

  yt_str <- if (!is.null(yt) && nrow(yt)==1 && !is.na(yt$did))
              sprintf("%+.4f  %.3f  %4d", yt$did, yt$p, yt$n)
            else "    —       —     —  "
  nt_str <- if (!is.null(nt) && nrow(nt)==1 && !is.na(nt$did))
              sprintf("%+.4f  %.3f  %4d", nt$did, nt$p, nt$n)
            else "    —       —     —  "

  pri <- if (h == PRIMARY_H) " << PRIMARY" else ""
  cat(sprintf("  %4dh     %s             %s%s\n", h, yt_str, nt_str, pri))
}

# Matched counts
cat(sprintf("\n  Matched pairs: yet-untreated=%s, never-treated=%s\n",
            if (!is.null(out_yt)) as.character(out_yt$n_matched) else "—",
            if (!is.null(out_nt)) as.character(out_nt$n_matched) else "—"))
cat(sprintf("  Max SMD: yet-untreated=%.3f (%d viol), never-treated=%.3f (%d viol)\n",
            if (!is.null(out_yt)) out_yt$balance$max_smd else NA,
            if (!is.null(out_yt)) out_yt$balance$n_viol else NA,
            if (!is.null(out_nt)) out_nt$balance$max_smd else NA,
            if (!is.null(out_nt)) out_nt$balance$n_viol else NA))

# ═══════════════════════════════════════════════════════════════════
# SAVE
# ═══════════════════════════════════════════════════════════════════

all_res <- rbind(
  if (!is.null(out_yt)) out_yt$results else NULL,
  if (!is.null(out_nt)) out_nt$results else NULL)
all_res$db <- db

write.csv(all_res,
          file.path(RESULTS, sprintf("did_riskset_%s.csv", tag)),
          row.names = FALSE)

# Save matched pairs for HTE
if (!is.null(out_yt)) {
  write.csv(out_yt$pairs,
            file.path(RESULTS, sprintf("did_pairs_yt_%s.csv", tag)),
            row.names = FALSE)
}
if (!is.null(out_nt)) {
  write.csv(out_nt$pairs,
            file.path(RESULTS, sprintf("did_pairs_nt_%s.csv", tag)),
            row.names = FALSE)
}

cat(sprintf("\n%s\n02_psm.R — %s COMPLETE\n", SEP, db))
cat(sprintf("  Output: did_riskset_%s.csv  (both pools)\n", tag))
cat(sprintf("          did_pairs_yt_%s.csv (yet-untreated matched pairs)\n", tag))
cat(sprintf("          did_pairs_nt_%s.csv (never-treated matched pairs)\n", tag))
cat(sprintf("  Next:   python 04_fig_timecourse.py\n"))
cat(sprintf("%s\n", SEP))
