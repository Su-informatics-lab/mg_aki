#!/usr/bin/env Rscript
# ============================================================================
# 02_psm.R — Risk-Set PSM with Temporal Alignment (v5 Design)
#
# Runs BOTH lab timing versions (last = primary, first = sensitivity)
# and compares results side-by-side.
#
# Design:
#   1. Pre-compute first+last lab values from did_labs_all
#   2. Pre-compute risk sets (shared, timing-independent)
#   3. For each timing {last, first}:
#      a. MICE m=20 on PS covariates
#      b. Per imputation: fit PS → risk-set match → compute ΔCr
#      c. Pool with Rubin's rules
#   4. Compare last vs first results
#
# Usage: Rscript 02_psm.R eicu
#        Rscript 02_psm.R mimic
# ============================================================================

suppressPackageStartupMessages({ library(sandwich); library(lmtest); library(mice) })

RESULTS    <- path.expand("~/mg_aki/results")
PRIMARY_H  <- 48       # primary: ΔCr at T₀ + 48h
CR_WINDOW  <- 12       # ±12h tolerance for Cr lookup
CALIPER_SD <- 0.2
M_IMP      <- 20
TARGETS    <- c(6, 12, 18, 24, 30, 36, 42, 48)

# PS: time-invariant covariates (always in model)
PS_FIXED <- c("age","is_female","bmi",
              "surg_cabg","surg_valve","surg_combined",
              "heart_failure","hypertension","diabetes","ckd",
              "copd","pvd","stroke","liver_disease","egfr",
              "ppi_chronic","loop_diuretic_chronic","acei_arb_chronic","nsaid_chronic")

# PS: time-varying labs (prefix with first_ or last_ depending on timing)
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

# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) { cat("Usage: Rscript 02_psm.R <db>\n"); quit(status = 1) }
db <- toupper(args[1]); tag <- tolower(db)

SEP <- paste(rep("=", 70), collapse = "")
cat(sprintf("\n%s\n02_psm.R — Risk-Set PSM: %s\n%s\n", SEP, db, SEP))

# ── Load data ─────────────────────────────────────────────────────────────
all_pts <- read.csv(file.path(RESULTS, sprintf("did_all_%s.csv", tag)), stringsAsFactors = FALSE)
cr_all  <- read.csv(file.path(RESULTS, sprintf("did_cr_all_%s.csv", tag)), stringsAsFactors = FALSE)

# Normalize cr_all
cr_id <- if ("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
cr_all$pid <- cr_all[[cr_id]]
if (!"labresultoffset" %in% names(cr_all)) cr_all$labresultoffset <- cr_all$offset_min
if (!"offset_h" %in% names(cr_all)) cr_all$offset_h <- cr_all$labresultoffset / 60
cr_all <- cr_all[order(cr_all$pid, cr_all$offset_h), ]

# Load labs (skip HR for now — 15M rows too large; load electrolytes only)
cat("  Loading labs (electrolytes only, skipping HR)...\n")
labs_raw <- read.csv(file.path(RESULTS, sprintf("did_labs_all_%s.csv", tag)), stringsAsFactors = FALSE)
pid_col_lab <- if ("patientunitstayid" %in% names(labs_raw)) "patientunitstayid" else "stay_id"
if (pid_col_lab %in% names(labs_raw)) labs_raw$pid <- labs_raw[[pid_col_lab]]
# Filter: electrolytes + HR (subsample HR to first per hour to reduce size)
labs_elec <- labs_raw[labs_raw$lab_name %in% c("magnesium","potassium","calcium","lactate"), ]
labs_hr   <- labs_raw[labs_raw$lab_name == "heartrate", ]
if (nrow(labs_hr) > 500000) {
  labs_hr$hour_bin <- floor(labs_hr$offset_h)
  labs_hr <- labs_hr[order(labs_hr$pid, labs_hr$offset_h), ]
  labs_hr <- labs_hr[!duplicated(paste(labs_hr$pid, labs_hr$hour_bin)), ]
  labs_hr$hour_bin <- NULL
  cat(sprintf("    HR downsampled: %d → %d\n", sum(labs_raw$lab_name=="heartrate"), nrow(labs_hr)))
}
labs <- rbind(labs_elec, labs_hr)
rm(labs_raw, labs_elec, labs_hr); gc()

N <- nrow(all_pts)
trt_idx <- which(all_pts$treated == 1 & !is.na(all_pts$mg_offset_h))
n_trt <- length(trt_idx)
cat(sprintf("  Patients: %d (%d treated, %d control)\n", N, n_trt, N - n_trt))

# ── Pre-compute per-patient lab values (BOTH first and last) ──────────────
cat("  Computing first/last lab values...\n")
for (ln in c("magnesium","potassium","calcium","lactate","heartrate")) {
  sub <- labs[labs$lab_name == ln, ]
  if (nrow(sub) == 0) next
  sub <- merge(sub, all_pts[, c("pid","mg_offset_h")], by = "pid")
  # For treated: [0, t_mg); for controls: [0, ∞)
  sub <- sub[sub$offset_h >= 0 & (is.na(sub$mg_offset_h) | sub$offset_h < sub$mg_offset_h), ]
  if (nrow(sub) == 0) next

  # First
  s1 <- sub[order(sub$offset_h), ]
  s1 <- s1[!duplicated(s1$pid), c("pid","value","offset_h")]
  names(s1) <- c("pid", paste0("first_",ln), paste0("first_",ln,"_time_h"))
  all_pts <- merge(all_pts, s1, by = "pid", all.x = TRUE)

  # Last
  s2 <- sub[order(-sub$offset_h), ]
  s2 <- s2[!duplicated(s2$pid), c("pid","value","offset_h")]
  names(s2) <- c("pid", paste0("last_",ln), paste0("last_",ln,"_time_h"))
  all_pts <- merge(all_pts, s2, by = "pid", all.x = TRUE)

  nf <- sum(!is.na(all_pts[[paste0("first_",ln)]]))
  nl <- sum(!is.na(all_pts[[paste0("last_",ln)]]))
  cat(sprintf("    %s: first=%d (%.0f%%), last=%d (%.0f%%)\n",
              ln, nf, 100*nf/N, nl, 100*nl/N))
}

# Lactate missing indicator
all_pts$first_lactate_missing <- as.integer(is.na(all_pts$first_lactate))
all_pts$last_lactate_missing  <- as.integer(is.na(all_pts$last_lactate))

# ── Pre-compute risk sets (timing-independent) ────────────────────────────
cat("  Pre-computing risk sets...\n")

# Per-patient Cr list (for AKI check + outcome computation)
cr_list <- split(cr_all[, c("labresult","offset_h")], cr_all$pid)

# Earliest Cr time
earliest_cr <- sapply(cr_list, function(x) min(x$offset_h, na.rm = TRUE))
all_pts$earliest_cr_h <- earliest_cr[as.character(all_pts$pid)]
all_pts$earliest_cr_h[is.na(all_pts$earliest_cr_h)] <- Inf

# First AKI time
cat("    Computing first AKI time...\n")
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

# Exit time: when patient becomes ineligible as control
all_pts$exit_time <- pmin(
  all_pts$icu_discharge_h,
  ifelse(is.na(all_pts$mg_offset_h), Inf, all_pts$mg_offset_h),
  ifelse(is.na(all_pts$first_aki_h), Inf, all_pts$first_aki_h),
  na.rm = TRUE)

cat(sprintf("    Median exit_time (controls): %.1fh\n",
            median(all_pts$exit_time[all_pts$treated == 0], na.rm = TRUE)))

# Pre-compute risk set INDICES for each treated patient
# (These don't depend on imputation or lab timing)
cat("    Building risk set indices...\n")
trt_pids <- all_pts$pid[trt_idx]
trt_tmg  <- all_pts$mg_offset_h[trt_idx]
risk_sets <- vector("list", n_trt)

for (k in seq_len(n_trt)) {
  t_mg <- trt_tmg[k]
  eligible <- which(
    all_pts$exit_time > t_mg &
    all_pts$earliest_cr_h <= t_mg &
    all_pts$pid != trt_pids[k])
  risk_sets[[k]] <- eligible
}

rs_sizes <- vapply(risk_sets, length, integer(1))
cat(sprintf("    Risk set sizes: median=%.0f, IQR=[%.0f,%.0f], min=%.0f, empty=%d\n",
            median(rs_sizes), quantile(rs_sizes, 0.25), quantile(rs_sizes, 0.75),
            min(rs_sizes), sum(rs_sizes == 0)))

# ── DIAGNOSTIC PROBE ──────────────────────────────────────────────
cat(sprintf("\n  ── DIAGNOSTIC: Why are risk sets empty? ──\n"))

has_rs  <- rs_sizes > 0
no_rs   <- rs_sizes == 0
cat(sprintf("    Matchable: %d (%.1f%%), Unmatchable: %d (%.1f%%)\n",
            sum(has_rs), 100*mean(has_rs), sum(no_rs), 100*mean(no_rs)))

# t_mg distribution
tmg_match   <- trt_tmg[has_rs]
tmg_nomatch <- trt_tmg[no_rs]
cat(sprintf("    t_mg (matchable):   median=%.1fh, IQR=[%.1f,%.1f]\n",
            median(tmg_match, na.rm=TRUE), quantile(tmg_match,0.25,na.rm=TRUE),
            quantile(tmg_match,0.75,na.rm=TRUE)))
cat(sprintf("    t_mg (unmatchable): median=%.1fh, IQR=[%.1f,%.1f]\n",
            median(tmg_nomatch, na.rm=TRUE), quantile(tmg_nomatch,0.25,na.rm=TRUE),
            quantile(tmg_nomatch,0.75,na.rm=TRUE)))

# Why unmatchable? Check each condition
cat("\n    Unmatchable breakdown (checking at t_mg):\n")
n_no_ctl_alive <- 0; n_no_ctl_cr <- 0
for (k in which(no_rs)) {
  t_mg <- trt_tmg[k]
  alive_icu <- sum(all_pts$exit_time > t_mg & all_pts$pid != trt_pids[k])
  has_cr    <- sum(all_pts$exit_time > t_mg & all_pts$earliest_cr_h <= t_mg &
                   all_pts$pid != trt_pids[k])
  if (alive_icu == 0) n_no_ctl_alive <- n_no_ctl_alive + 1
  else if (has_cr == 0) n_no_ctl_cr <- n_no_ctl_cr + 1
}
cat(sprintf("    No one alive/in-ICU at t_mg: %d\n", n_no_ctl_alive))
cat(sprintf("    Alive but no Cr before t_mg: %d\n", n_no_ctl_cr))
cat(sprintf("    (total unmatchable: %d)\n", sum(no_rs)))

# Covariate comparison: matchable vs unmatchable treated
cat("\n    Covariate means (matchable vs unmatchable treated):\n")
cat("    Variable              matchable  unmatch  diff\n")
probe_vars <- c("age","is_female","hypertension","diabetes","ckd",
                "heart_failure","surg_cabg","egfr","bmi")
probe_vars <- intersect(probe_vars, names(all_pts))
for (v in probe_vars) {
  m1 <- mean(all_pts[[v]][trt_idx[has_rs]], na.rm=TRUE)
  m2 <- mean(all_pts[[v]][trt_idx[no_rs]], na.rm=TRUE)
  cat(sprintf("    %-20s  %7.3f    %7.3f  %+.3f\n", v, m1, m2, m1-m2))
}

# ═══════════════════════════════════════════════════════════════════
# RUN ONE LAB TIMING VERSION
# ═══════════════════════════════════════════════════════════════════
run_timing <- function(timing) {
  cat(sprintf("\n%s\n  LAB TIMING: %s\n%s\n", SEP, toupper(timing), SEP))

  # Select PS covariates
  lab_cols <- paste0(timing, "_", PS_LAB_BASES)
  lac_miss <- paste0(timing, "_lactate_missing")
  ps_vars <- c(intersect(PS_FIXED, names(all_pts)),
               intersect(c(lab_cols, lac_miss), names(all_pts)))
  ps_vars <- ps_vars[vapply(ps_vars, function(v) {
    x <- all_pts[[v]]; !all(is.na(x)) && var(x, na.rm = TRUE) > 1e-10
  }, logical(1))]
  cat(sprintf("  PS covariates (%d): %s\n", length(ps_vars), paste(ps_vars, collapse=", ")))

  # ── MICE ──────────────────────────────────────────────────────
  to_impute <- ps_vars[vapply(ps_vars, function(v) any(is.na(all_pts[[v]])), logical(1))]
  cat(sprintf("  MICE m=%d on %d vars: %s\n", M_IMP, length(to_impute), paste(to_impute, collapse=", ")))

  if (length(to_impute) > 0) {
    mice_df <- all_pts[, c("treated", ps_vars)]
    meth <- rep("", ncol(mice_df)); names(meth) <- names(mice_df)
    for (v in to_impute) meth[v] <- "pmm"
    imp <- mice(mice_df, m = M_IMP, method = meth, printFlag = FALSE, maxit = 10)
    cat(sprintf("  MICE done. Logged events: %d\n", nrow(imp$loggedEvents)))
  } else {
    imp <- NULL
  }

  # ── Per-imputation matching + outcome ─────────────────────────
  all_dids <- list()
  m1_match_trt <- NULL; m1_match_ctl <- NULL; m1_matched <- NULL; m1_d <- NULL

  for (m_idx in 1:M_IMP) {
    if (m_idx %% 5 == 1) cat(sprintf("  m=%d...", m_idx))

    # Get imputed data
    d <- all_pts  # copy
    if (!is.null(imp)) {
      imp_df <- complete(imp, m_idx)
      for (v in to_impute) d[[v]] <- imp_df[[v]]
    }

    # Fit global PS
    ps_fml <- as.formula(paste("treated ~", paste(ps_vars, collapse = "+")))
    ps_fit <- tryCatch(glm(ps_fml, data = d, family = binomial()),
                       error = function(e) { cat("PS failed\n"); return(NULL) })
    if (is.null(ps_fit)) next
    d$ps <- predict(ps_fit, type = "response")
    caliper <- CALIPER_SD * sd(d$ps, na.rm = TRUE)

    # Risk-set matching (with replacement)
    match_trt <- integer(n_trt); match_ctl <- integer(n_trt)
    match_ps  <- numeric(n_trt); matched <- logical(n_trt)

    for (k in seq_len(n_trt)) {
      rs <- risk_sets[[k]]
      if (length(rs) == 0) next
      ps_i <- d$ps[trt_idx[k]]
      ps_dist <- abs(d$ps[rs] - ps_i)
      within_cal <- which(ps_dist <= caliper)
      if (length(within_cal) == 0) next
      best <- within_cal[which.min(ps_dist[within_cal])]
      match_trt[k] <- trt_idx[k]
      match_ctl[k] <- rs[best]
      match_ps[k]  <- ps_dist[best]
      matched[k]   <- TRUE
    }

    n_matched <- sum(matched)
    if (n_matched < 50) next

    # Save m=1 for balance reporting
    if (m_idx == 1) {
      m1_match_trt <- match_trt; m1_match_ctl <- match_ctl
      m1_matched <- matched; m1_d <- d
    }

    # Compute outcomes
    for (target_h in TARGETS) {
      dcr_trt <- dcr_ctl <- numeric(n_matched)
      valid <- logical(n_matched); idx <- 0
      matched_k <- which(matched)

      for (kk in matched_k) {
        idx <- idx + 1
        tpid <- as.character(all_pts$pid[match_trt[kk]])
        cpid <- as.character(all_pts$pid[match_ctl[kk]])
        t_mg <- trt_tmg[kk]

        pre_t <- find_cr_pre(cr_list[[tpid]], t_mg)
        pre_c <- find_cr_pre(cr_list[[cpid]], t_mg)
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
        all_dids[[length(all_dids)+1]] <- data.frame(
          timing=timing, m=m_idx, target_h=target_h, n=n_valid,
          did=NA, se=NA, stringsAsFactors=FALSE)
        next
      }

      pair_df <- data.frame(
        delta_cr = c(dcr_trt[valid], dcr_ctl[valid]),
        treated = rep(c(1,0), each = n_valid),
        stringsAsFactors = FALSE)
      fit <- lm(delta_cr ~ treated, data = pair_df)
      ct <- safe_coeftest(fit)
      est <- if (!is.null(ct) && "treated" %in% rownames(ct)) ct["treated",1] else coef(fit)["treated"]
      se  <- if (!is.null(ct) && "treated" %in% rownames(ct)) ct["treated",2] else NA

      all_dids[[length(all_dids)+1]] <- data.frame(
        timing=timing, m=m_idx, target_h=target_h, n=n_valid,
        did=est, se=se, n_matched=n_matched, stringsAsFactors=FALSE)
    }
  }
  cat("\n")

  res <- do.call(rbind, all_dids)

  # ── Rubin's rules pooling ─────────────────────────────────────
  cat(sprintf("\n  ── POOLED (%s) ──\n", toupper(timing)))
  cat("  target   DiD        SE       P        95%% CI                n    FMI\n")
  cat("  ──────   ────────   ──────   ──────   ────────────────────   ───  ─────\n")

  pooled <- list()
  for (th in TARGETS) {
    sub <- res[res$target_h == th & !is.na(res$did) & !is.na(res$se), ]
    if (nrow(sub) < 2) {
      cat(sprintf("  %4dh    (insufficient data)\n", th)); next
    }
    Q <- mean(sub$did); U <- mean(sub$se^2); B <- var(sub$did)
    m_v <- nrow(sub)
    T_var <- U + (1 + 1/m_v) * B; T_se <- sqrt(T_var)
    lam <- ((1 + 1/m_v) * B) / T_var
    df <- max(3, (m_v - 1) / lam^2)
    p <- 2 * pt(abs(Q/T_se), df = df, lower.tail = FALSE)
    ci_lo <- Q - qt(0.975, df) * T_se
    ci_hi <- Q + qt(0.975, df) * T_se
    sig <- if (!is.na(p) && p < 0.05) " *" else "  "
    pri <- if (th == PRIMARY_H) " << PRIMARY" else ""
    n_avg <- mean(sub$n)

    cat(sprintf("  %4dh    %+.4f   %.4f   %.4f%s [%+.4f,%+.4f]   %3.0f  %.3f%s\n",
                th, Q, T_se, p, sig, ci_lo, ci_hi, n_avg, lam, pri))
    pooled[[length(pooled)+1]] <- data.frame(
      timing=timing, target_h=th, did=Q, se=T_se, p=p,
      ci_lo=ci_lo, ci_hi=ci_hi, fmi=lam, n_avg=n_avg, stringsAsFactors=FALSE)
  }

  # Quality metrics (from m=1)
  m1 <- res[res$m == 1, ]
  if (nrow(m1) > 0) {
    cat(sprintf("\n  Matched (m=1): %d/%d (%.1f%%)\n",
                m1$n_matched[1], n_trt, 100*m1$n_matched[1]/n_trt))
  }

  # ── BALANCE (SMD before/after matching, m=1) ────────────────
  if (!is.null(m1_d) && !is.null(m1_matched) && sum(m1_matched) > 50) {
    cat("\n  ── COVARIATE BALANCE (m=1) ──\n")
    cat("  Covariate                raw_SMD  matched_SMD  status\n")
    cat("  ───────────────────────  ───────  ───────────  ──────\n")

    trt_all <- which(m1_d$treated == 1)
    ctl_all <- which(m1_d$treated == 0)
    trt_m <- m1_match_trt[m1_matched]
    ctl_m <- m1_match_ctl[m1_matched]
    n_viol <- 0

    balance_rows <- list()
    for (v in ps_vars) {
      x1_raw <- m1_d[[v]][trt_all]; x0_raw <- m1_d[[v]][ctl_all]
      x1_mat <- m1_d[[v]][trt_m];   x0_mat <- m1_d[[v]][ctl_m]
      sp_raw <- sqrt((var(x1_raw, na.rm=T) + var(x0_raw, na.rm=T)) / 2)
      sp_mat <- sqrt((var(x1_mat, na.rm=T) + var(x0_mat, na.rm=T)) / 2)
      smd_raw <- if (!is.na(sp_raw) && sp_raw > 1e-10)
                   abs(mean(x1_raw,na.rm=T) - mean(x0_raw,na.rm=T)) / sp_raw else NA
      smd_mat <- if (!is.na(sp_mat) && sp_mat > 1e-10)
                   abs(mean(x1_mat,na.rm=T) - mean(x0_mat,na.rm=T)) / sp_mat else NA
      flag <- if (!is.na(smd_mat) && smd_mat > 0.1) { n_viol <- n_viol+1; "VIOL" } else "ok"
      cat(sprintf("  %-25s  %.3f    %.3f       %s\n", v,
                  ifelse(is.na(smd_raw), NA, smd_raw),
                  ifelse(is.na(smd_mat), NA, smd_mat), flag))
      balance_rows[[length(balance_rows)+1]] <- data.frame(
        timing=timing, covariate=v, smd_raw=smd_raw, smd_mat=smd_mat, stringsAsFactors=FALSE)
    }
    cat(sprintf("\n  Max matched SMD: %.3f | Violations (>0.1): %d/%d\n",
                max(sapply(balance_rows, function(r) r$smd_mat), na.rm=TRUE),
                n_viol, length(ps_vars)))

    # Save balance
    bal_df <- do.call(rbind, balance_rows)
    write.csv(bal_df, file.path(RESULTS, sprintf("did_balance_%s_%s.csv", timing, tag)),
              row.names = FALSE)
  }

  do.call(rbind, pooled)
}

# ═══════════════════════════════════════════════════════════════════
# RUN BOTH TIMINGS
# ═══════════════════════════════════════════════════════════════════

cat(sprintf("\n%s\n  RUNNING BOTH LAB TIMINGS\n%s\n", SEP, SEP))

res_last  <- run_timing("last")
res_first <- run_timing("first")

# Save all results
all_res <- rbind(res_last, res_first)
write.csv(all_res, file.path(RESULTS, sprintf("did_riskset_%s.csv", tag)), row.names = FALSE)

# ── Side-by-side comparison ───────────────────────────────────────────────
cat(sprintf("\n%s\n  COMPARISON: LAST vs FIRST labs\n%s\n", SEP, SEP))
cat("  target   LAST_DiD    P        FIRST_DiD   P        |diff|\n")
cat("  ──────   ─────────   ──────   ─────────   ──────   ──────\n")

for (th in TARGETS) {
  rl <- res_last[res_last$target_h == th, ]
  rf <- res_first[res_first$target_h == th, ]
  if (nrow(rl) > 0 && nrow(rf) > 0) {
    d <- abs(rl$did - rf$did)
    pri <- if (th == PRIMARY_H) " << PRIMARY" else ""
    cat(sprintf("  %4dh    %+.4f   %.4f   %+.4f   %.4f   %.4f%s\n",
                th, rl$did, rl$p, rf$did, rf$p, d, pri))
  }
}

cat(sprintf("\n%s\n02_psm.R — %s COMPLETE\n%s\n", SEP, db, SEP))
cat(sprintf("  Output: did_riskset_%s.csv\n", tag))
