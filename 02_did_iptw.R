#!/usr/bin/env Rscript
# ============================================================================
# did_iptw.R — IPTW variant comparison (all ATE, ICU-time anchor)
#
# Variants:
#   sIPTW:        stabilized IPTW, no trimming
#   sIPTW_t99:    stabilized + trim at 1st/99th percentile
#   sIPTW_t95:    stabilized + trim at 5th/95th percentile
#   sIPTW_DR:     stabilized trimmed + covariate adjustment (doubly robust)
#   AIPW:         augmented IPW (outcome model backup)
#
# All use: same PS model, same ICU-time ΔCr definition, HC robust SEs
# Sweep: 6-36h every 3h from ICU admission
#
# Run:  Rscript did_iptw.R eicu
#       Rscript did_iptw.R mimic
# ============================================================================

suppressPackageStartupMessages({
  library(sandwich)
  library(lmtest)
})

RESULTS <- path.expand("~/mg_aki/results")
TARGETS <- seq(6, 36, by = 3)

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

median_impute <- function(d, vars) {
  for (v in vars)
    if (v %in% names(d) && any(is.na(d[[v]])))
      d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm = TRUE)
  d
}

# ============================================================================
run_iptw <- function(db) {
  tag <- tolower(db)
  SEP <- paste(rep("=", 70), collapse = "")
  cat(sprintf("\n%s\n%s: IPTW Variants (ATE, ICU-time anchor)\n%s\n", SEP, db, SEP))

  # ── Load ─────────────────────────────────────────────────────────────────
  trt <- read.csv(file.path(RESULTS, sprintf("did_treated_%s.csv", tag)), stringsAsFactors = F)
  ctl <- read.csv(file.path(RESULTS, sprintf("did_control_%s.csv", tag)), stringsAsFactors = F)
  cr_all <- read.csv(file.path(RESULTS, sprintf("did_cr_all_%s.csv", tag)), stringsAsFactors = F)

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
  n_trt <- sum(combined$treated == 1)
  n_ctl <- sum(combined$treated == 0)
  prev <- n_trt / nrow(combined)  # treatment prevalence for stabilization
  cat(sprintf("  %d treated + %d control (prevalence = %.1f%%)\n", n_trt, n_ctl, 100 * prev))

  # ── PS model ─────────────────────────────────────────────────────────────
  cat("  Fitting PS model...\n")
  avail <- ps_vars[sapply(ps_vars, function(v) v %in% names(combined) && var(combined[[v]], na.rm = T) > 1e-10)]
  ps_fml <- as.formula(paste("treated ~", paste(avail, collapse = "+")))
  ps_fit <- glm(ps_fml, data = combined, family = binomial())
  combined$ps <- predict(ps_fit, type = "response")

  # PS diagnostics
  ps_trt <- combined$ps[combined$treated == 1]
  ps_ctl <- combined$ps[combined$treated == 0]
  cat(sprintf("  PS distribution:\n"))
  cat(sprintf("    Treated:  median=%.3f, IQR=[%.3f,%.3f], range=[%.4f,%.4f]\n",
              median(ps_trt), quantile(ps_trt, 0.25), quantile(ps_trt, 0.75),
              min(ps_trt), max(ps_trt)))
  cat(sprintf("    Control:  median=%.3f, IQR=[%.3f,%.3f], range=[%.4f,%.4f]\n",
              median(ps_ctl), quantile(ps_ctl, 0.25), quantile(ps_ctl, 0.75),
              min(ps_ctl), max(ps_ctl)))

  # ── Compute ALL weight variants ──────────────────────────────────────────
  ps <- combined$ps; trt_flag <- combined$treated

  # Unstabilized IPTW
  w_unstab <- ifelse(trt_flag == 1, 1 / ps, 1 / (1 - ps))

  # Stabilized IPTW
  w_stab <- ifelse(trt_flag == 1, prev / ps, (1 - prev) / (1 - ps))

  # Trimmed at 1st/99th
  trim99 <- function(w) pmin(pmax(w, quantile(w, 0.01)), quantile(w, 0.99))
  w_t99 <- trim99(w_stab)

  # Trimmed at 5th/95th
  trim95 <- function(w) pmin(pmax(w, quantile(w, 0.05)), quantile(w, 0.95))
  w_t95 <- trim95(w_stab)

  cat(sprintf("\n  Weight diagnostics:\n"))
  for (nm in c("unstab", "stab", "t99", "t95")) {
    w <- get(paste0("w_", nm))
    cat(sprintf("    %-8s: median=%.2f, max=%.1f, CV=%.2f, ESS_trt=%.0f, ESS_ctl=%.0f\n",
                nm, median(w), max(w), sd(w) / mean(w),
                sum(w[trt_flag == 1])^2 / sum(w[trt_flag == 1]^2),
                sum(w[trt_flag == 0])^2 / sum(w[trt_flag == 0]^2)))
  }

  # ── Weighted balance check (sIPTW trimmed 99) ───────────────────────────
  cat("\n  Weighted balance (sIPTW_t99):\n")
  max_wsmd <- 0
  for (v in avail) {
    x1 <- combined[[v]][trt_flag == 1]; w1 <- w_t99[trt_flag == 1]
    x0 <- combined[[v]][trt_flag == 0]; w0 <- w_t99[trt_flag == 0]
    m1 <- weighted.mean(x1, w1, na.rm = T); m0 <- weighted.mean(x0, w0, na.rm = T)
    sp <- sqrt((var(x1, na.rm = T) + var(x0, na.rm = T)) / 2)
    if (!is.na(sp) && sp > 1e-10) {
      wsmd <- abs(m1 - m0) / sp
      if (wsmd > max_wsmd) max_wsmd <- wsmd
    }
  }
  cat(sprintf("    Max weighted SMD: %.4f\n", max_wsmd))

  # ── Cr_pre for everyone (first postop Cr) ───────────────────────────────
  # Treated: use their cr_pre (before IV Mg)
  # Controls: use first_postop_cr
  cr_pre_trt <- setNames(trt$cr_pre, trt$pid)
  cr_pre_off_trt <- setNames(trt$cr_pre_offset_min, trt$pid)

  cr_pre_ctl_col <- if ("first_postop_cr" %in% names(ctl)) "first_postop_cr" else "cr_pre"
  cr_pre_ctl <- setNames(ctl[[cr_pre_ctl_col]], ctl$pid)
  cr_pre_off_ctl_col <- if ("first_cr_offset_min" %in% names(ctl)) "first_cr_offset_min" else "cr_pre_offset_min"
  cr_pre_off_ctl <- if (cr_pre_off_ctl_col %in% names(ctl)) setNames(ctl[[cr_pre_off_ctl_col]], ctl$pid) else NULL

  combined$cr_pre <- ifelse(combined$treated == 1,
                            cr_pre_trt[as.character(combined$pid)],
                            cr_pre_ctl[as.character(combined$pid)])

  # ── Sweep timing targets ─────────────────────────────────────────────────
  cat(sprintf("\n  Sweeping %d time targets...\n", length(TARGETS)))
  all_results <- list(); ridx <- 0

  for (target_h in TARGETS) {
    target_min <- target_h * 60

    # Find Cr closest to target_h from ICU for each patient
    # Must be after their cr_pre
    cr_post_vals <- c()
    cr_post_pids <- c()

    for (pid_i in combined$pid) {
      cr_i <- cr_all[cr_all$pid == pid_i, ]
      cp <- combined$cr_pre[combined$pid == pid_i][1]
      if (is.na(cp) || nrow(cr_i) == 0) next

      # Get cr_pre offset for this patient
      if (as.character(pid_i) %in% names(cr_pre_off_trt) && combined$treated[combined$pid == pid_i][1] == 1) {
        pre_off <- cr_pre_off_trt[as.character(pid_i)]
      } else if (!is.null(cr_pre_off_ctl) && as.character(pid_i) %in% names(cr_pre_off_ctl)) {
        pre_off <- cr_pre_off_ctl[as.character(pid_i)]
      } else {
        pre_off <- 0
      }

      # Cr after cr_pre, closest to target
      cands <- cr_i[cr_i$labresultoffset > pre_off, ]
      if (nrow(cands) == 0) next
      cands$dist <- abs(cands$labresultoffset - target_min)
      best <- cands[which.min(cands$dist), ]
      cr_post_vals <- c(cr_post_vals, best$labresult)
      cr_post_pids <- c(cr_post_pids, pid_i)
    }

    cr_post_df <- data.frame(pid = cr_post_pids, cr_post = cr_post_vals, stringsAsFactors = F)
    df <- merge(combined, cr_post_df, by = "pid")
    df$delta_cr <- df$cr_post - df$cr_pre
    df <- df[!is.na(df$delta_cr), ]

    nt <- sum(df$treated == 1); nc <- sum(df$treated == 0)
    if (nt < 20 || nc < 20) next

    # Run each IPTW variant
    variants <- list(
      sIPTW = w_stab,
      sIPTW_t99 = w_t99,
      sIPTW_t95 = w_t95
    )

    for (vname in names(variants)) {
      w_full <- variants[[vname]]
      # Subset weights to df rows
      w_df <- w_full[match(df$pid, combined$pid)]

      # Weighted regression
      fit <- lm(delta_cr ~ treated, data = df, weights = w_df)
      cl <- vcovHC(fit, type = "HC1")
      ct <- coeftest(fit, vcov. = cl)

      ridx <- ridx + 1
      all_results[[ridx]] <- data.frame(
        method = vname, target_h = target_h,
        n_trt = nt, n_ctl = nc,
        did = round(ct["treated", "Estimate"], 4),
        se = round(ct["treated", "Std. Error"], 4),
        p = round(ct["treated", "Pr(>|t|)"], 4),
        ci_lo = round(ct["treated", "Estimate"] - 1.96 * ct["treated", "Std. Error"], 4),
        ci_hi = round(ct["treated", "Estimate"] + 1.96 * ct["treated", "Std. Error"], 4),
        stringsAsFactors = F)
    }

    # sIPTW_DR: trimmed + covariate adjustment
    w_df <- w_t99[match(df$pid, combined$pid)]
    adj_vars <- intersect(avail, names(df))
    adj_vars <- adj_vars[sapply(adj_vars, function(v) var(df[[v]], na.rm = T) > 1e-10)]
    fml_dr <- as.formula(paste("delta_cr ~ treated +", paste(adj_vars, collapse = "+")))
    fit_dr <- tryCatch(lm(fml_dr, data = df, weights = w_df), error = function(e) NULL)
    if (!is.null(fit_dr)) {
      cl_dr <- vcovHC(fit_dr, type = "HC1")
      ct_dr <- coeftest(fit_dr, vcov. = cl_dr)
      ridx <- ridx + 1
      all_results[[ridx]] <- data.frame(
        method = "sIPTW_DR", target_h = target_h,
        n_trt = nt, n_ctl = nc,
        did = round(ct_dr["treated", "Estimate"], 4),
        se = round(ct_dr["treated", "Std. Error"], 4),
        p = round(ct_dr["treated", "Pr(>|t|)"], 4),
        ci_lo = round(ct_dr["treated", "Estimate"] - 1.96 * ct_dr["treated", "Std. Error"], 4),
        ci_hi = round(ct_dr["treated", "Estimate"] + 1.96 * ct_dr["treated", "Std. Error"], 4),
        stringsAsFactors = F)
    }

    # AIPW: augmented IPW
    # Fit outcome models separately for treated and control
    df_t <- df[df$treated == 1, ]; df_c <- df[df$treated == 0, ]
    out_fml <- as.formula(paste("delta_cr ~", paste(adj_vars, collapse = "+")))
    mu1_fit <- tryCatch(lm(out_fml, data = df_t), error = function(e) NULL)
    mu0_fit <- tryCatch(lm(out_fml, data = df_c), error = function(e) NULL)

    if (!is.null(mu1_fit) && !is.null(mu0_fit)) {
      # Predict potential outcomes for everyone
      mu1_hat <- predict(mu1_fit, newdata = df)
      mu0_hat <- predict(mu0_fit, newdata = df)
      ps_df <- df$ps <- combined$ps[match(df$pid, combined$pid)]

      # AIPW estimate
      aipw_scores <- (mu1_hat - mu0_hat) +
        df$treated * (df$delta_cr - mu1_hat) / ps_df -
        (1 - df$treated) * (df$delta_cr - mu0_hat) / (1 - ps_df)

      # Trim extreme scores (cap at 1st/99th)
      q01 <- quantile(aipw_scores, 0.01, na.rm = T)
      q99 <- quantile(aipw_scores, 0.99, na.rm = T)
      aipw_scores <- pmin(pmax(aipw_scores, q01), q99)

      tau_aipw <- mean(aipw_scores, na.rm = T)
      se_aipw <- sd(aipw_scores, na.rm = T) / sqrt(sum(!is.na(aipw_scores)))
      p_aipw <- 2 * pnorm(-abs(tau_aipw / se_aipw))

      ridx <- ridx + 1
      all_results[[ridx]] <- data.frame(
        method = "AIPW", target_h = target_h,
        n_trt = nt, n_ctl = nc,
        did = round(tau_aipw, 4),
        se = round(se_aipw, 4),
        p = round(p_aipw, 4),
        ci_lo = round(tau_aipw - 1.96 * se_aipw, 4),
        ci_hi = round(tau_aipw + 1.96 * se_aipw, 4),
        stringsAsFactors = F)
    }

    # Progress
    sig_any <- any(sapply(all_results[max(1,ridx-4):ridx], function(r) r$p < 0.05))
    cat(sprintf("    %3dh: n=%d/%d %s\n", target_h, nt, nc,
                if (sig_any) "← has significant" else ""))
  }

  # ── Results ──────────────────────────────────────────────────────────────
  res <- do.call(rbind, all_results)
  write.csv(res, file.path(RESULTS, sprintf("did_iptw_%s.csv", tag)), row.names = F)

  # Print formatted tables by method
  cat(sprintf("\n%s\n%s: IPTW VARIANT RESULTS\n%s\n", SEP, db, SEP))

  methods <- unique(res$method)
  for (meth in methods) {
    sub <- res[res$method == meth, ]
    cat(sprintf("\n  ── %s ──\n", meth))
    cat("  target  n_trt  n_ctl     DiD      SE      P       95%% CI\n")
    cat("  ──────  ─────  ─────  ────────  ──────  ──────  ──────────────\n")
    for (i in seq_len(nrow(sub))) {
      r <- sub[i, ]
      sig <- if (r$p < 0.05) " *" else ""
      cat(sprintf("  %4dh   %5d  %5d  %+.4f  %.4f  %.4f  [%+.4f,%+.4f]%s\n",
                  r$target_h, r$n_trt, r$n_ctl, r$did, r$se, r$p,
                  r$ci_lo, r$ci_hi, sig))
    }
  }

  # Cross-method comparison at key time points
  cat(sprintf("\n%s\nCROSS-METHOD COMPARISON\n%s\n", SEP, SEP))
  cat("\n  target   sIPTW     sIPTW_t99  sIPTW_t95  sIPTW_DR   AIPW\n")
  cat("  ──────  ─────────  ─────────  ─────────  ─────────  ─────────\n")
  for (th in TARGETS) {
    vals <- character(5)
    for (j in seq_along(methods)) {
      r <- res[res$method == methods[j] & res$target_h == th, ]
      if (nrow(r) > 0) {
        sig <- if (r$p[1] < 0.05) "*" else " "
        vals[j] <- sprintf("%+.3f%s", r$did[1], sig)
      } else vals[j] <- "     NA"
    }
    cat(sprintf("  %4dh   %s  %s  %s  %s  %s\n",
                th, vals[1], vals[2], vals[3], vals[4], vals[5]))
  }

  cat(sprintf("\n  * P < 0.05\n"))
  cat(sprintf("  All methods estimate ATE using ICU-time anchor\n"))
  cat(sprintf("  sIPTW_DR = stabilized trimmed IPTW + covariate adjustment (doubly robust)\n"))
  cat(sprintf("  AIPW = augmented IPW with outcome model (theoretical gold standard)\n"))
  cat(sprintf("\n  Saved: did_iptw_%s.csv\n", tag))

  return(res)
}

# ============================================================================
cat("======================================================================\n")
cat("did_iptw.R — IPTW variant comparison (all ATE)\n")
cat("  Variants: sIPTW, sIPTW_t99, sIPTW_t95, sIPTW_DR, AIPW\n")
cat(sprintf("  Targets: %s h from ICU\n", paste(TARGETS, collapse = ",")))
cat("======================================================================\n")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) { cat("Usage: Rscript did_iptw.R eicu|mimic\n"); quit(status = 1) }
for (a in args) run_iptw(toupper(a))

cat("\n======================================================================\n")
cat("如果5个方法趋势一致 → 结果robust\n")
cat("如果AIPW和sIPTW_DR一致但sIPTW不一致 → extreme weights问题\n")
cat("如果都在同一个时间窗显著 → 信号真实\n")
cat("======================================================================\n")
