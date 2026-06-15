#!/usr/bin/env Rscript
# ============================================================================
# 08_mg_stratified.R — The Entanglement Test
#
# Question: Does the Mg supplementation → AKI protective effect survive
#           within narrow bands of first postoperative serum Mg?
#
# Why this matters:
#   Supplemented patients have lower serum Mg (1.88 vs 2.22 in eICU).
#   Lower serum Mg may mark less cardioplegia → simpler surgery → lower AKI.
#   If the effect persists WITHIN a narrow Mg band (where complexity is
#   approximately held constant), it can't be explained by this pathway.
#   If it vanishes, the reviewer is right.
#
# Design:
#   1. Stratify by clinical Mg bands: <1.8, 1.8–2.0, 2.0–2.3, >2.3
#   2. Within each stratum, re-estimate PS WITHOUT first_mg_value
#      (it's now held constant by stratification → remove from PS)
#   3. Compute OW and AC-OW treatment effects within each stratum
#   4. Formal interaction: treatment × continuous Mg in a single model
#   5. Cochran Q test for heterogeneity across strata
#
# Output: results/08_mg_stratified.csv
# Run:    Rscript 08_mg_stratified.R
# ============================================================================

suppressPackageStartupMessages({
  library(sandwich); library(lmtest)
})
RESULTS <- path.expand("~/mg_aki/results")

# ── Standardize column names ─────────────────────────────────────────────
stdz <- function(d) {
  rmap <- c(mg_supplementation="mg_supp", age_num="age",
    baseline_cr="baseline_creatinine", baseline_egfr="egfr",
    hx_chf="heart_failure", hx_hypertension="hypertension",
    hx_diabetes="diabetes", hx_ckd="ckd", hx_copd="copd",
    hx_pvd="pvd", hx_stroke="stroke", hx_liver="liver_disease",
    nephrotox_loop_diuretic="loop_diuretics", nephrotox_nsaid="nsaids",
    nephrotox_acei_arb="acei_arb", nephrotox_ppi="ppi",
    has_betablocker="beta_blockers", has_steroid="steroids",
    preop_antiarrhythmic="antiarrhythmics",
    first_k_value="first_potassium", first_ca_value="first_calcium",
    first_hr="first_heartrate", has_vasopressor="vasopressor_6h")
  for (old in names(rmap)) {
    new <- rmap[[old]]
    if (old %in% names(d) && !new %in% names(d)) names(d)[names(d)==old] <- new
  }
  if (is.character(d$age)) {
    d$age <- suppressWarnings(as.numeric(d$age))
    d$age[is.na(d$age)] <- 90
  }
  if ("surgery_type" %in% names(d)) {
    d$surg_cabg     <- as.integer(d$surgery_type == "cabg")
    d$surg_valve    <- as.integer(d$surgery_type == "valve")
    d$surg_combined <- as.integer(d$surgery_type == "combined")
  }
  # Lactate: median + indicator
  if ("first_lactate" %in% names(d)) {
    d$lactate_missing <- as.integer(is.na(d$first_lactate))
    d$first_lactate[is.na(d$first_lactate)] <- median(d$first_lactate, na.rm = TRUE)
  }
  # Median-impute MICE targets for simplicity in stratified analysis
  for (v in c("bmi", "first_heartrate", "first_calcium", "first_potassium")) {
    if (v %in% names(d) && any(is.na(d[[v]])))
      d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm = TRUE)
  }
  d
}

# ── PS covariates — WITHOUT first_mg_value ───────────────────────────────
# (removed because we're stratifying on it — it's held constant)
ps_covars_base <- c(
  "age", "is_female", "bmi",
  "surg_cabg", "surg_valve", "surg_combined",
  "heart_failure", "hypertension", "diabetes", "ckd",
  "copd", "pvd", "stroke", "liver_disease",
  "baseline_creatinine", "egfr",
  "loop_diuretics", "nsaids", "acei_arb", "ppi",
  "beta_blockers", "steroids", "antiarrhythmics",
  "first_potassium", "first_calcium", "first_heartrate",
  "vasopressor_6h"
  # first_mg_value deliberately EXCLUDED
)

# ── Helper: weighted GLM ─────────────────────────────────────────────────
wglm <- function(formula, data, w, cluster = NULL) {
  data$.w <- w
  fit <- glm(formula, data = data, weights = .w, family = quasibinomial())
  vc <- if (!is.null(cluster) && length(unique(cluster)) > 1)
    vcovCL(fit, cluster = cluster) else vcovHC(fit, type = "HC1")
  ct <- coeftest(fit, vcov. = vc)
  tr <- 2  # treatment is always second coefficient
  list(logOR = ct[tr, 1], se = ct[tr, 2],
       or = exp(ct[tr, 1]),
       lo = exp(ct[tr, 1] - 1.96 * ct[tr, 2]),
       hi = exp(ct[tr, 1] + 1.96 * ct[tr, 2]),
       p  = 2 * pnorm(-abs(ct[tr, 1] / ct[tr, 2])))
}

# ── Helper: SMD ──────────────────────────────────────────────────────────
smd_w <- function(x, trt, w) {
  m1 <- weighted.mean(x[trt == 1], w[trt == 1], na.rm = TRUE)
  m0 <- weighted.mean(x[trt == 0], w[trt == 0], na.rm = TRUE)
  sp <- sqrt((var(x[trt == 1], na.rm = TRUE) + var(x[trt == 0], na.rm = TRUE)) / 2)
  if (sp < 1e-10) return(0)
  abs(m1 - m0) / sp
}

# ── Mg strata ────────────────────────────────────────────────────────────
MG_CUTS   <- c(0, 1.8, 2.0, 2.3, Inf)
MG_LABELS <- c("<1.8", "1.8-2.0", "2.0-2.3", ">2.3")

# ============================================================================
# MAIN ANALYSIS FUNCTION
# ============================================================================
run_stratified <- function(dat, db_name, has_cluster = FALSE) {

  cat(sprintf("\n%s\n%s  (N=%d, trt=%d)\n%s\n",
              strrep("=", 65), db_name, nrow(dat), sum(dat$mg_supp), strrep("=", 65)))

  ps_covars <- intersect(ps_covars_base, names(dat))
  if ("first_lactate" %in% names(dat))
    ps_covars <- c(ps_covars, "first_lactate", "lactate_missing")

  cat(sprintf("  PS model: %d covariates (first_mg_value EXCLUDED)\n", length(ps_covars)))

  # Assign strata
  dat$mg_stratum <- cut(dat$first_mg_value, breaks = MG_CUTS, labels = MG_LABELS,
                        right = FALSE, include.lowest = TRUE)

  cat("\n  Mg stratum distribution:\n")
  for (s in MG_LABELS) {
    sub <- dat[dat$mg_stratum == s, ]
    cat(sprintf("    %-10s N=%5d  trt=%4d (%.1f%%)  AKI=%.1f%%  Mg: %.2f±%.2f\n",
                s, nrow(sub), sum(sub$mg_supp), 100*mean(sub$mg_supp),
                100*mean(sub$aki_kdigo1), mean(sub$first_mg_value), sd(sub$first_mg_value)))
  }

  results <- list()

  # ── Stratum-specific analyses ──────────────────────────────────────
  for (s in MG_LABELS) {
    d_s <- dat[dat$mg_stratum == s, ]
    d_s <- d_s[complete.cases(d_s[, ps_covars]), ]
    n_trt <- sum(d_s$mg_supp)
    n_ctrl <- sum(d_s$mg_supp == 0)

    cat(sprintf("\n  ── Stratum %s (N=%d, trt=%d, ctrl=%d) ──\n", s, nrow(d_s), n_trt, n_ctrl))

    if (n_trt < 15 || n_ctrl < 15) {
      cat("    SKIP: insufficient sample\n")
      next
    }

    # ── All-patient PS + OW ──────────────────────────────────────
    ps_fml <- as.formula(paste("mg_supp ~", paste(ps_covars, collapse = " + ")))
    tryCatch({
      ps_fit <- glm(ps_fml, data = d_s, family = binomial())
      d_s$ps <- pmax(pmin(fitted(ps_fit), 0.99), 0.01)
      d_s$ow <- ifelse(d_s$mg_supp == 1, 1 - d_s$ps, d_s$ps)

      # Balance check
      smds <- sapply(ps_covars, function(v)
        if (is.numeric(d_s[[v]])) smd_w(d_s[[v]], d_s$mg_supp, d_s$ow) else NA)
      max_smd <- max(smds, na.rm = TRUE)

      cl <- if (has_cluster && "hospitalid" %in% names(d_s)) d_s$hospitalid else NULL
      r <- wglm(aki_kdigo1 ~ mg_supp, d_s, d_s$ow, cl)
      sig <- ifelse(r$p < 0.05, " *", "")
      cat(sprintf("    All-pt OW: OR %.3f (%.3f–%.3f) P=%.4f  maxSMD=%.3f%s\n",
                  r$or, r$lo, r$hi, r$p, max_smd, sig))

      results[[length(results) + 1]] <- data.frame(
        db = db_name, stratum = s, analysis = "all_ow",
        n = nrow(d_s), n_trt = n_trt,
        mg_mean = mean(d_s$first_mg_value), mg_sd = sd(d_s$first_mg_value),
        or = round(r$or, 3), lo = round(r$lo, 3), hi = round(r$hi, 3),
        p = round(r$p, 4), logOR = round(r$logOR, 4), se = round(r$se, 4),
        max_smd = round(max_smd, 4))
    }, error = function(e) cat(sprintf("    All-pt OW FAILED: %s\n", e$message)))

    # ── AC: Mg+K vs K-only within stratum ────────────────────────
    if ("ac_group" %in% names(d_s)) {
      d_ac <- d_s[d_s$ac_group %in% c("mg_k", "k_only"), ]
      d_ac$ac_trt <- as.integer(d_ac$ac_group == "mg_k")
      n_ac_trt  <- sum(d_ac$ac_trt)
      n_ac_ctrl <- sum(d_ac$ac_trt == 0)

      if (n_ac_trt >= 10 && n_ac_ctrl >= 10) {
        tryCatch({
          ac_fml <- as.formula(paste("ac_trt ~", paste(ps_covars, collapse = " + ")))
          ac_fit <- glm(ac_fml, data = d_ac, family = binomial())
          d_ac$ac_ps <- pmax(pmin(fitted(ac_fit), 0.99), 0.01)
          d_ac$ac_ow <- ifelse(d_ac$ac_trt == 1, 1 - d_ac$ac_ps, d_ac$ac_ps)

          smds_ac <- sapply(ps_covars, function(v)
            if (is.numeric(d_ac[[v]])) smd_w(d_ac[[v]], d_ac$ac_trt, d_ac$ac_ow) else NA)
          max_smd_ac <- max(smds_ac, na.rm = TRUE)

          r_ac <- wglm(aki_kdigo1 ~ ac_trt, d_ac, d_ac$ac_ow, NULL)
          sig <- ifelse(r_ac$p < 0.05, " *", "")
          cat(sprintf("    AC OW:     OR %.3f (%.3f–%.3f) P=%.4f  maxSMD=%.3f%s  [n=%d, trt=%d]\n",
                      r_ac$or, r_ac$lo, r_ac$hi, r_ac$p, max_smd_ac, sig, nrow(d_ac), n_ac_trt))

          results[[length(results) + 1]] <- data.frame(
            db = db_name, stratum = s, analysis = "ac_ow",
            n = nrow(d_ac), n_trt = n_ac_trt,
            mg_mean = mean(d_ac$first_mg_value), mg_sd = sd(d_ac$first_mg_value),
            or = round(r_ac$or, 3), lo = round(r_ac$lo, 3), hi = round(r_ac$hi, 3),
            p = round(r_ac$p, 4), logOR = round(r_ac$logOR, 4), se = round(r_ac$se, 4),
            max_smd = round(max_smd_ac, 4))
        }, error = function(e) cat(sprintf("    AC OW FAILED: %s\n", e$message)))
      } else {
        cat(sprintf("    AC SKIP: trt=%d, ctrl=%d\n", n_ac_trt, n_ac_ctrl))
      }
    }
  }

  # ── Formal interaction test: treatment × continuous Mg ─────────────
  cat("\n  ── Interaction Test ──\n")
  d_int <- dat[complete.cases(dat[, ps_covars]), ]
  int_fml <- as.formula(paste("aki_kdigo1 ~ mg_supp * first_mg_value +",
                               paste(ps_covars, collapse = " + ")))
  tryCatch({
    int_fit <- glm(int_fml, data = d_int, family = binomial())
    ct <- coeftest(int_fit, vcov. = vcovHC(int_fit, type = "HC1"))
    int_row <- grep("mg_supp:first_mg_value", rownames(ct))
    if (length(int_row) > 0) {
      int_or <- exp(ct[int_row, 1])
      int_p  <- ct[int_row, 4]
      cat(sprintf("    mg_supp × first_mg_value: OR %.3f, P=%.4f\n", int_or, int_p))
      cat(sprintf("    Interpretation: %s\n",
                  ifelse(int_p < 0.05,
                         "SIGNIFICANT interaction — effect varies by Mg level",
                         "No significant interaction — effect is consistent across Mg levels")))
      results[[length(results) + 1]] <- data.frame(
        db = db_name, stratum = "interaction", analysis = "trt_x_mg",
        n = nrow(d_int), n_trt = sum(d_int$mg_supp),
        mg_mean = NA, mg_sd = NA,
        or = round(int_or, 3), lo = NA, hi = NA,
        p = round(int_p, 4), logOR = round(ct[int_row, 1], 4),
        se = round(ct[int_row, 2], 4), max_smd = NA)
    }
  }, error = function(e) cat(sprintf("    Interaction test FAILED: %s\n", e$message)))

  # ── Cochran Q for heterogeneity across strata ──────────────────────
  cat("\n  ── Heterogeneity (Cochran Q across strata) ──\n")
  ow_res <- do.call(rbind, results)
  ow_res <- ow_res[ow_res$db == db_name & ow_res$analysis == "all_ow" &
                    !is.na(ow_res$logOR), ]
  if (nrow(ow_res) >= 2) {
    w <- 1 / ow_res$se^2
    pool <- sum(w * ow_res$logOR) / sum(w)
    Q <- sum(w * (ow_res$logOR - pool)^2)
    df <- nrow(ow_res) - 1
    pQ <- pchisq(Q, df, lower.tail = FALSE)
    I2 <- max(0, (Q - df) / Q) * 100
    cat(sprintf("    Q=%.2f, df=%d, P=%.4f, I²=%.0f%%\n", Q, df, pQ, I2))
    cat(sprintf("    %s\n", ifelse(pQ < 0.10,
      "SIGNIFICANT heterogeneity — effect differs across Mg strata",
      "No significant heterogeneity — effect is consistent across Mg strata")))
  }

  do.call(rbind, results)
}

# ============================================================================
# RUN
# ============================================================================
cat(strrep("=", 65), "\n")
cat("08: SERUM-MG-STRATIFIED ANALYSIS (The Entanglement Test)\n")
cat(strrep("=", 65), "\n")
cat("If the protective effect persists within narrow Mg bands,\n")
cat("it cannot be explained by cardioplegia-driven complexity.\n")

dat_e <- stdz(read.csv(file.path(RESULTS, "01_analysis_a_cohort.csv"), stringsAsFactors = FALSE))
dat_m <- stdz(read.csv(file.path(RESULTS, "04_mimic_cohort.csv"), stringsAsFactors = FALSE))

res_e <- run_stratified(dat_e, "eICU", has_cluster = TRUE)
res_m <- run_stratified(dat_m, "MIMIC", has_cluster = FALSE)

all_res <- rbind(res_e, res_m)
outpath <- file.path(RESULTS, "08_mg_stratified.csv")
write.csv(all_res, outpath, row.names = FALSE)

# ============================================================================
# VERDICT
# ============================================================================
cat(sprintf("\n%s\nVERDICT\n%s\n", strrep("=", 65), strrep("=", 65)))

# Check if signal persists in the 1.8-2.0 or 2.0-2.3 bands
for (db in c("eICU", "MIMIC")) {
  cat(sprintf("\n  %s:\n", db))
  db_res <- all_res[all_res$db == db & all_res$analysis == "all_ow" &
                    all_res$stratum %in% MG_LABELS, ]
  if (nrow(db_res) == 0) { cat("    No results\n"); next }
  for (i in seq_len(nrow(db_res))) {
    r <- db_res[i, ]
    dir <- ifelse(r$or < 1, "protective", "harmful/null")
    sig <- ifelse(r$p < 0.05, "SIGNIFICANT", "ns")
    cat(sprintf("    Mg %s: OR=%.3f P=%.4f [%s, %s]\n",
                r$stratum, r$or, r$p, dir, sig))
  }
}

cat("\n  Decision framework:\n")
cat("    If OR<1 persists in >=2 strata → complexity confounding unlikely\n")
cat("    If OR<1 only in <1.8 stratum   → drug effect in depleted patients\n")
cat("    If OR~1 in all strata          → original signal was confounded\n")
cat("    If OR<1 only in >2.3 stratum   → reverse pattern, complexity artifact\n")

cat(sprintf("\n✓ Saved: %s (%d rows)\n", outpath, nrow(all_res)))
