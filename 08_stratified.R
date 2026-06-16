#!/usr/bin/env Rscript
# ============================================================================
# 08_stratified.R — Mg-stratified analysis + hospital random effects
#
# Consolidates: 08_mg_stratified.R, 08b_hospital_re.R
#
# Sections:
#   A. Mg-stratified treatment effects     → results/08_mg_stratified.csv
#   B. Hospital random effects (eICU)      → results/08b_hospital_re.csv
#
# Run: Rscript 08_stratified.R            # both sections
#      Rscript 08_stratified.R A          # stratified only
#      Rscript 08_stratified.R B          # hospital RE only
# ============================================================================

suppressPackageStartupMessages({
  library(sandwich); library(lmtest); library(lme4)
})
RESULTS <- path.expand("~/mg_aki/results")

# ── Shared: column standardization ────────────────────────────────────────
stdz <- function(d) {
  rmap <- c(
    mg_supplementation      = "mg_supp",
    age_num                  = "age",
    baseline_cr              = "baseline_creatinine",
    baseline_egfr            = "egfr",
    hx_chf                  = "heart_failure",
    hx_hypertension         = "hypertension",
    hx_diabetes             = "diabetes",
    hx_ckd                  = "ckd",
    hx_copd                 = "copd",
    hx_pvd                  = "pvd",
    hx_stroke               = "stroke",
    hx_liver                = "liver_disease",
    nephrotox_loop_diuretic  = "loop_diuretics",
    nephrotox_nsaid         = "nsaids",
    nephrotox_acei_arb      = "acei_arb",
    nephrotox_ppi           = "ppi",
    has_betablocker         = "beta_blockers",
    has_steroid             = "steroids",
    preop_antiarrhythmic    = "antiarrhythmics",
    first_k_value           = "first_potassium",
    first_ca_value          = "first_calcium",
    first_hr                = "first_heartrate",
    has_vasopressor         = "vasopressor_6h"
  )
  for (old in names(rmap)) {
    new <- rmap[[old]]
    if (old %in% names(d) && !new %in% names(d))
      names(d)[names(d) == old] <- new
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
  if ("first_lactate" %in% names(d)) {
    d$lactate_missing <- as.integer(is.na(d$first_lactate))
    d$first_lactate[is.na(d$first_lactate)] <- median(d$first_lactate,
                                                       na.rm = TRUE)
  }
  for (v in c("bmi", "first_heartrate", "first_calcium", "first_potassium"))
    if (v %in% names(d) && any(is.na(d[[v]])))
      d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm = TRUE)
  d
}

# ── Shared: PS covariates WITHOUT first_mg_value ──────────────────────────
# (removed because stratification holds it constant)
ps_covars_base <- c(
  "age", "is_female", "bmi",
  "surg_cabg", "surg_valve", "surg_combined",
  "heart_failure", "hypertension", "diabetes", "ckd",
  "copd", "pvd", "stroke", "liver_disease",
  "baseline_creatinine", "egfr",
  "loop_diuretics", "nsaids", "acei_arb", "ppi",
  "beta_blockers", "steroids", "antiarrhythmics",
  "first_potassium", "first_calcium", "first_heartrate",
  "vasopressor_6h")

# ── Shared: helpers ───────────────────────────────────────────────────────
wglm <- function(formula, data, w, cluster = NULL) {
  data$.w <- w
  fit <- glm(formula, data = data, weights = .w, family = quasibinomial())
  vc <- if (!is.null(cluster) && length(unique(cluster)) > 1)
    vcovCL(fit, cluster = cluster) else vcovHC(fit, type = "HC1")
  ct <- coeftest(fit, vcov. = vc)
  tr <- 2
  list(logOR = ct[tr, 1], se = ct[tr, 2],
       or = exp(ct[tr, 1]),
       lo = exp(ct[tr, 1] - 1.96 * ct[tr, 2]),
       hi = exp(ct[tr, 1] + 1.96 * ct[tr, 2]),
       p  = 2 * pnorm(-abs(ct[tr, 1] / ct[tr, 2])))
}

smd_w <- function(x, trt, w) {
  m1 <- weighted.mean(x[trt == 1], w[trt == 1], na.rm = TRUE)
  m0 <- weighted.mean(x[trt == 0], w[trt == 0], na.rm = TRUE)
  sp <- sqrt((var(x[trt == 1], na.rm = TRUE) +
              var(x[trt == 0], na.rm = TRUE)) / 2)
  if (sp < 1e-10) return(0)
  abs(m1 - m0) / sp
}

MG_CUTS   <- c(0, 1.8, 2.0, 2.3, Inf)
MG_LABELS <- c("<1.8", "1.8-2.0", "2.0-2.3", ">2.3")

# ── Load data ─────────────────────────────────────────────────────────────
cat("Loading cohorts...\n")
dat_e <- stdz(read.csv(file.path(RESULTS, "01_analysis_a_cohort.csv"),
                        stringsAsFactors = FALSE))
dat_m <- stdz(read.csv(file.path(RESULTS, "04_mimic_cohort.csv"),
                        stringsAsFactors = FALSE))
cat(sprintf("  eICU: N=%d  MIMIC: N=%d\n", nrow(dat_e), nrow(dat_m)))

args <- commandArgs(trailingOnly = TRUE)
run_all <- length(args) == 0
sections <- if (run_all) c("A", "B") else toupper(args)

# ============================================================================
# SECTION A: MG-STRATIFIED TREATMENT EFFECTS
# ============================================================================
if ("A" %in% sections) {
  cat(sprintf("\n%s\nA. SERUM-MG-STRATIFIED ANALYSIS\n%s\n",
              strrep("=", 65), strrep("=", 65)))

  run_stratified <- function(dat, db_name, has_cluster = FALSE) {
    cat(sprintf("\n  %s (N=%d, trt=%d)\n", db_name, nrow(dat),
                sum(dat$mg_supp)))

    ps_covars <- intersect(ps_covars_base, names(dat))
    if ("first_lactate" %in% names(dat))
      ps_covars <- c(ps_covars, "first_lactate", "lactate_missing")

    dat$mg_stratum <- cut(dat$first_mg_value, breaks = MG_CUTS,
                          labels = MG_LABELS, right = FALSE,
                          include.lowest = TRUE)

    for (s in MG_LABELS) {
      sub <- dat[dat$mg_stratum == s, ]
      cat(sprintf("    %-10s N=%5d trt=%4d (%.1f%%) AKI=%.1f%%\n",
                  s, nrow(sub), sum(sub$mg_supp),
                  100 * mean(sub$mg_supp), 100 * mean(sub$aki_kdigo1)))
    }

    results <- list()

    # ── Per-stratum PS + OW ─────────────────────────────────────
    for (s in MG_LABELS) {
      d_s <- dat[dat$mg_stratum == s, ]
      d_s <- d_s[complete.cases(d_s[, ps_covars]), ]
      n_trt <- sum(d_s$mg_supp)
      n_ctrl <- sum(d_s$mg_supp == 0)
      if (n_trt < 15 || n_ctrl < 15) {
        cat(sprintf("    %-10s SKIP\n", s)); next
      }

      # All-patient OW
      ps_fml <- as.formula(paste("mg_supp ~",
                                  paste(ps_covars, collapse = " + ")))
      tryCatch({
        ps_fit <- glm(ps_fml, data = d_s, family = binomial())
        d_s$ps <- pmax(pmin(fitted(ps_fit), 0.99), 0.01)
        d_s$ow <- ifelse(d_s$mg_supp == 1, 1 - d_s$ps, d_s$ps)
        smds <- sapply(ps_covars, function(v)
          if (is.numeric(d_s[[v]])) smd_w(d_s[[v]], d_s$mg_supp,
                                           d_s$ow) else NA)
        max_smd <- max(smds, na.rm = TRUE)
        cl <- if (has_cluster && "hospitalid" %in% names(d_s))
          d_s$hospitalid else NULL
        r <- wglm(aki_kdigo1 ~ mg_supp, d_s, d_s$ow, cl)
        sig <- ifelse(r$p < 0.05, " *", "")
        cat(sprintf("    %-10s OW:  OR %.3f (%.3f-%.3f) P=%.4f%s\n",
                    s, r$or, r$lo, r$hi, r$p, sig))
        results[[length(results) + 1]] <- data.frame(
          db = db_name, stratum = s, analysis = "all_ow",
          n = nrow(d_s), n_trt = n_trt,
          mg_mean = mean(d_s$first_mg_value),
          mg_sd = sd(d_s$first_mg_value),
          or = round(r$or, 3), lo = round(r$lo, 3), hi = round(r$hi, 3),
          p = round(r$p, 4), logOR = round(r$logOR, 4),
          se = round(r$se, 4), max_smd = round(max_smd, 4))
      }, error = function(e) cat(sprintf("    %-10s OW FAILED\n", s)))

      # AC OW within stratum
      if ("ac_group" %in% names(d_s)) {
        d_ac <- d_s[d_s$ac_group %in% c("mg_k", "k_only"), ]
        d_ac$ac_trt <- as.integer(d_ac$ac_group == "mg_k")
        if (sum(d_ac$ac_trt) >= 10 && sum(d_ac$ac_trt == 0) >= 10) {
          tryCatch({
            ac_fml <- as.formula(paste("ac_trt ~",
                                        paste(ps_covars, collapse = "+")))
            ac_fit <- glm(ac_fml, data = d_ac, family = binomial())
            d_ac$ac_ps <- pmax(pmin(fitted(ac_fit), 0.99), 0.01)
            d_ac$ac_ow <- ifelse(d_ac$ac_trt == 1,
                                  1 - d_ac$ac_ps, d_ac$ac_ps)
            smds_ac <- sapply(ps_covars, function(v)
              if (is.numeric(d_ac[[v]])) smd_w(d_ac[[v]], d_ac$ac_trt,
                                                d_ac$ac_ow) else NA)
            r_ac <- wglm(aki_kdigo1 ~ ac_trt, d_ac, d_ac$ac_ow, NULL)
            sig <- ifelse(r_ac$p < 0.05, " *", "")
            cat(sprintf("    %-10s AC:  OR %.3f (%.3f-%.3f) P=%.4f%s\n",
                        s, r_ac$or, r_ac$lo, r_ac$hi, r_ac$p, sig))
            results[[length(results) + 1]] <- data.frame(
              db = db_name, stratum = s, analysis = "ac_ow",
              n = nrow(d_ac), n_trt = sum(d_ac$ac_trt),
              mg_mean = mean(d_ac$first_mg_value),
              mg_sd = sd(d_ac$first_mg_value),
              or = round(r_ac$or, 3), lo = round(r_ac$lo, 3),
              hi = round(r_ac$hi, 3), p = round(r_ac$p, 4),
              logOR = round(r_ac$logOR, 4), se = round(r_ac$se, 4),
              max_smd = round(max(smds_ac, na.rm = TRUE), 4))
          }, error = function(e) cat(sprintf("    %-10s AC FAILED\n", s)))
        }
      }
    }

    # ── Interaction test ────────────────────────────────────────
    d_int <- dat[complete.cases(dat[, ps_covars]), ]
    int_fml <- as.formula(paste("aki_kdigo1 ~ mg_supp * first_mg_value +",
                                 paste(ps_covars, collapse = " + ")))
    tryCatch({
      int_fit <- glm(int_fml, data = d_int, family = binomial())
      ct <- coeftest(int_fit, vcov. = vcovHC(int_fit, type = "HC1"))
      int_row <- grep("mg_supp:first_mg_value", rownames(ct))
      if (length(int_row) > 0) {
        int_or <- exp(ct[int_row, 1]); int_p <- ct[int_row, 4]
        cat(sprintf("    Interaction: OR %.3f P=%.4f\n", int_or, int_p))
        results[[length(results) + 1]] <- data.frame(
          db = db_name, stratum = "interaction", analysis = "trt_x_mg",
          n = nrow(d_int), n_trt = sum(d_int$mg_supp),
          mg_mean = NA, mg_sd = NA,
          or = round(int_or, 3), lo = NA, hi = NA,
          p = round(int_p, 4), logOR = round(ct[int_row, 1], 4),
          se = round(ct[int_row, 2], 4), max_smd = NA)
      }
    }, error = function(e) cat("    Interaction FAILED\n"))

    # ── Cochran Q ───────────────────────────────────────────────
    ow_res <- do.call(rbind, results)
    ow_strata <- ow_res[ow_res$db == db_name & ow_res$analysis == "all_ow" &
                         !is.na(ow_res$logOR), ]
    if (nrow(ow_strata) >= 2) {
      w <- 1 / ow_strata$se^2
      pool <- sum(w * ow_strata$logOR) / sum(w)
      Q <- sum(w * (ow_strata$logOR - pool)^2)
      df <- nrow(ow_strata) - 1
      pQ <- pchisq(Q, df, lower.tail = FALSE)
      I2 <- max(0, (Q - df) / Q) * 100
      cat(sprintf("    Cochran Q=%.2f P=%.4f I²=%.0f%%\n", Q, pQ, I2))
    }

    do.call(rbind, results)
  }

  res_a_e <- run_stratified(dat_e, "eICU", has_cluster = TRUE)
  res_a_m <- run_stratified(dat_m, "MIMIC", has_cluster = FALSE)
  res_a <- rbind(res_a_e, res_a_m)
  write.csv(res_a, file.path(RESULTS, "08_mg_stratified.csv"),
            row.names = FALSE)
  cat(sprintf("  Saved: 08_mg_stratified.csv (%d rows)\n", nrow(res_a)))
}

# ============================================================================
# SECTION B: HOSPITAL RANDOM EFFECTS (eICU only)
# ============================================================================
if ("B" %in% sections) {
  cat(sprintf("\n%s\nB. HOSPITAL RANDOM EFFECTS (eICU)\n%s\n",
              strrep("=", 65), strrep("=", 65)))

  dat <- dat_e
  if (!"hospitalid" %in% names(dat)) {
    cat("  hospitalid not found — SKIP\n")
  } else {
    avail_covars <- intersect(ps_covars_base, names(dat))
    if ("first_lactate" %in% names(dat))
      avail_covars <- c(avail_covars, "first_lactate", "lactate_missing")
    dat <- dat[complete.cases(dat[, avail_covars]), ]
    dat$mg_stratum <- cut(dat$first_mg_value, breaks = MG_CUTS,
                          labels = MG_LABELS, right = FALSE,
                          include.lowest = TRUE)
    cov_str <- paste(avail_covars, collapse = " + ")
    cat(sprintf("  N=%d, hospitals=%d\n", nrow(dat),
                length(unique(dat$hospitalid))))

    results_b <- list()

    # ── ICC ─────────────────────────────────────────────────────
    icc_fit <- glmer(aki_kdigo1 ~ 1 + (1 | hospitalid), data = dat,
                     family = binomial, nAGQ = 1,
                     control = glmerControl(optimizer = "bobyqa"))
    vc <- as.data.frame(VarCorr(icc_fit))
    icc <- vc$vcov[1] / (vc$vcov[1] + pi^2 / 3)
    cat(sprintf("  ICC: %.3f (%.1f%% between hospitals)\n", icc, 100 * icc))

    # ── Hospital exposure variation ─────────────────────────────
    hosp_stats <- aggregate(mg_supp ~ hospitalid, data = dat,
      FUN = function(x) c(n = length(x), rate = mean(x)))
    hosp_stats <- do.call(data.frame, hosp_stats)
    names(hosp_stats) <- c("hospitalid", "n", "trt_rate")
    hosp_stats <- hosp_stats[hosp_stats$n >= 20, ]
    hosp_aki <- aggregate(aki_kdigo1 ~ hospitalid, data = dat, FUN = mean)
    names(hosp_aki) <- c("hospitalid", "aki_rate")
    hosp_merged <- merge(hosp_stats, hosp_aki)
    if (nrow(hosp_merged) >= 10) {
      cr <- cor.test(hosp_merged$trt_rate, hosp_merged$aki_rate)
      cat(sprintf("  Hospital correlation (supp rate vs AKI): r=%.3f P=%.4f\n",
                  cr$estimate, cr$p.value))
    }

    # ── Mixed-effects models per stratum ────────────────────────
    for (stratum in c("Overall", "<1.8", "2.0-2.3", ">2.3")) {
      d_s <- if (stratum == "Overall") dat else dat[dat$mg_stratum == stratum, ]
      n_trt <- sum(d_s$mg_supp)
      n_ctrl <- nrow(d_s) - n_trt
      n_hosp <- length(unique(d_s$hospitalid))
      if (n_trt < 15 || n_ctrl < 15 || n_hosp < 5) next

      hosp_n <- table(d_s$hospitalid)
      d_s <- d_s[d_s$hospitalid %in% names(hosp_n[hosp_n >= 2]), ]

      cat(sprintf("\n  ── %s (N=%d, trt=%d, hosp=%d) ──\n",
                  stratum, nrow(d_s), sum(d_s$mg_supp),
                  length(unique(d_s$hospitalid))))

      # Fixed (cluster-robust)
      fml_fixed <- as.formula(paste("aki_kdigo1 ~ mg_supp +", cov_str))
      tryCatch({
        fit_f <- glm(fml_fixed, data = d_s, family = binomial())
        ct_f <- coeftest(fit_f, vcov. = vcovCL(fit_f,
                                                cluster = d_s$hospitalid))
        or_f <- exp(ct_f[2, 1])
        lo_f <- exp(ct_f[2, 1] - 1.96 * ct_f[2, 2])
        hi_f <- exp(ct_f[2, 1] + 1.96 * ct_f[2, 2])
        p_f  <- 2 * pnorm(-abs(ct_f[2, 1] / ct_f[2, 2]))
        cat(sprintf("    Fixed:       OR %.3f (%.3f-%.3f) P=%.4f\n",
                    or_f, lo_f, hi_f, p_f))
        results_b[[length(results_b) + 1]] <- data.frame(
          stratum = stratum, model = "fixed_clusterSE",
          n = nrow(d_s), n_trt = sum(d_s$mg_supp),
          n_hosp = length(unique(d_s$hospitalid)),
          or = round(or_f, 3), lo = round(lo_f, 3), hi = round(hi_f, 3),
          p = round(p_f, 4))
      }, error = function(e) cat("    Fixed FAILED\n"))

      # Hospital random intercept
      fml_re <- as.formula(paste("aki_kdigo1 ~ mg_supp +", cov_str,
                                  "+ (1 | hospitalid)"))
      tryCatch({
        fit_re <- glmer(fml_re, data = d_s, family = binomial, nAGQ = 1,
          control = glmerControl(optimizer = "bobyqa",
                                 optCtrl = list(maxfun = 50000)))
        sf <- summary(fit_re)
        ct_re <- sf$coefficients
        tr <- which(rownames(ct_re) == "mg_supp")
        or_re <- exp(ct_re[tr, 1])
        lo_re <- exp(ct_re[tr, 1] - 1.96 * ct_re[tr, 2])
        hi_re <- exp(ct_re[tr, 1] + 1.96 * ct_re[tr, 2])
        p_re  <- 2 * pnorm(-abs(ct_re[tr, 1] / ct_re[tr, 2]))
        sig <- ifelse(p_re < 0.05, " *", "")
        vc_s <- as.data.frame(VarCorr(fit_re))
        icc_s <- vc_s$vcov[1] / (vc_s$vcov[1] + pi^2 / 3)
        cat(sprintf("    Hospital RE: OR %.3f (%.3f-%.3f) P=%.4f%s  ICC=%.3f\n",
                    or_re, lo_re, hi_re, p_re, sig, icc_s))
        if (exists("or_f")) {
          change <- (or_re - or_f) / or_f * 100
          cat(sprintf("    Delta: %+.1f%%\n", change))
        }
        results_b[[length(results_b) + 1]] <- data.frame(
          stratum = stratum, model = "hospital_RE",
          n = nrow(d_s), n_trt = sum(d_s$mg_supp),
          n_hosp = length(unique(d_s$hospitalid)),
          or = round(or_re, 3), lo = round(lo_re, 3),
          hi = round(hi_re, 3), p = round(p_re, 4))
      }, error = function(e) cat("    Hospital RE FAILED\n"))
    }

    # ── Narrow sub-bands within >2.3 ────────────────────────────
    cat(sprintf("\n  ── Narrow sub-bands within >2.3 ──\n"))
    d_hi <- dat[dat$mg_stratum == ">2.3", ]
    d_hi$mg_subband <- cut(d_hi$first_mg_value,
      breaks = c(2.3, 2.6, 3.0, Inf),
      labels = c("2.3-2.6", "2.6-3.0", ">3.0"),
      right = FALSE, include.lowest = TRUE)

    for (sb in c("2.3-2.6", "2.6-3.0", ">3.0")) {
      d_sb <- d_hi[d_hi$mg_subband == sb, ]
      n_trt <- sum(d_sb$mg_supp)
      n_ctrl <- nrow(d_sb) - n_trt
      if (n_trt < 10 || n_ctrl < 10) next
      fml <- as.formula(paste("aki_kdigo1 ~ mg_supp +", cov_str))
      tryCatch({
        fit <- glm(fml, data = d_sb, family = binomial())
        ct <- coeftest(fit, vcov. = vcovHC(fit, type = "HC1"))
        or <- exp(ct[2, 1])
        lo <- exp(ct[2, 1] - 1.96 * ct[2, 2])
        hi <- exp(ct[2, 1] + 1.96 * ct[2, 2])
        p  <- 2 * pnorm(-abs(ct[2, 1] / ct[2, 2]))
        sig <- ifelse(p < 0.05, " *", "")
        cat(sprintf("    %-10s N=%4d trt=%3d OR %.3f (%.3f-%.3f) P=%.4f  AKI: trt=%.0f%% ctrl=%.0f%%%s\n",
            sb, nrow(d_sb), n_trt, or, lo, hi, p,
            100 * mean(d_sb$aki_kdigo1[d_sb$mg_supp == 1]),
            100 * mean(d_sb$aki_kdigo1[d_sb$mg_supp == 0]), sig))
        results_b[[length(results_b) + 1]] <- data.frame(
          stratum = paste0(">2.3:", sb), model = "fixed_subband",
          n = nrow(d_sb), n_trt = n_trt,
          n_hosp = length(unique(d_sb$hospitalid)),
          or = round(or, 3), lo = round(lo, 3), hi = round(hi, 3),
          p = round(p, 4))
      }, error = function(e) cat(sprintf("    %-10s FAILED\n", sb)))
    }

    # ── Within-hospital meta-analysis (>2.3) ────────────────────
    cat(sprintf("\n  ── Within-hospital meta (>2.3) ──\n"))
    hosp_counts <- aggregate(mg_supp ~ hospitalid, data = d_hi,
      FUN = function(x) c(n = length(x), trt = sum(x), ctrl = sum(x == 0)))
    hosp_counts <- do.call(data.frame, hosp_counts)
    names(hosp_counts) <- c("hospitalid", "n", "n_trt", "n_ctrl")
    eligible_hosp <- hosp_counts$hospitalid[hosp_counts$n_trt >= 5 &
                                            hosp_counts$n_ctrl >= 5]
    cat(sprintf("    Hospitals with >=5 per arm: %d\n",
                length(eligible_hosp)))

    if (length(eligible_hosp) >= 3) {
      hosp_ors <- list()
      simple_covars <- intersect(c("age", "is_female",
        "baseline_creatinine", "surg_valve", "surg_combined",
        "vasopressor_6h"), names(d_hi))
      for (h in eligible_hosp) {
        d_h <- d_hi[d_hi$hospitalid == h, ]
        if (length(unique(d_h$aki_kdigo1)) < 2) next
        if (length(unique(d_h$mg_supp)) < 2) next
        fml_h <- as.formula(paste("aki_kdigo1 ~ mg_supp +",
                                   paste(simple_covars, collapse = "+")))
        tryCatch({
          fit_h <- glm(fml_h, data = d_h, family = binomial())
          ct_h <- coef(summary(fit_h))
          if ("mg_supp" %in% rownames(ct_h)) {
            hosp_ors[[length(hosp_ors) + 1]] <- data.frame(
              hospitalid = h, n = nrow(d_h), n_trt = sum(d_h$mg_supp),
              logOR = ct_h["mg_supp", 1], se = ct_h["mg_supp", 2])
          }
        }, error = function(e) {})
      }
      if (length(hosp_ors) >= 3) {
        ho <- do.call(rbind, hosp_ors)
        ho <- ho[abs(ho$logOR) < 5 & ho$se < 5, ]
        w <- 1 / ho$se^2
        pool_logOR <- sum(w * ho$logOR) / sum(w)
        pool_se <- sqrt(1 / sum(w))
        pool_or <- exp(pool_logOR)
        pool_lo <- exp(pool_logOR - 1.96 * pool_se)
        pool_hi <- exp(pool_logOR + 1.96 * pool_se)
        pool_p  <- 2 * pnorm(-abs(pool_logOR / pool_se))
        Q <- sum(w * (ho$logOR - pool_logOR)^2)
        I2 <- max(0, (Q - (nrow(ho) - 1)) / Q) * 100
        cat(sprintf("    Pooled: OR %.3f (%.3f-%.3f) P=%.4f I²=%.0f%%\n",
                    pool_or, pool_lo, pool_hi, pool_p, I2))
        results_b[[length(results_b) + 1]] <- data.frame(
          stratum = ">2.3", model = "within_hospital_meta",
          n = sum(ho$n), n_trt = sum(ho$n_trt), n_hosp = nrow(ho),
          or = round(pool_or, 3), lo = round(pool_lo, 3),
          hi = round(pool_hi, 3), p = round(pool_p, 4))
      }
    }

    all_b <- do.call(rbind, results_b)
    write.csv(all_b, file.path(RESULTS, "08b_hospital_re.csv"),
              row.names = FALSE)
    cat(sprintf("  Saved: 08b_hospital_re.csv (%d rows)\n", nrow(all_b)))
  }
}

cat(sprintf("\n%s\n08_stratified.R COMPLETE\n%s\n",
            strrep("=", 65), strrep("=", 65)))
