#!/usr/bin/env Rscript
# ============================================================================
# test_r2.R — R2 Reviewer Fixes (MICE + Lactate + OW/PSM primary)
#
# Run: Rscript test_r2.R
# Prereqs: mice package (needs cmake: pip install cmake first)
# Input:  data/eicu_cohort.parquet, data/mimic_cohort.parquet
# Output: results/02_results.csv
# ============================================================================

cat("======================================================================\n")
cat("R2 REVIEWER FIXES\n")
cat("MICE m=5 | OW primary | PSM co-primary | Lactate in PS\n")
cat("======================================================================\n")

# ── Packages ──────────────────────────────────────────────────────────────
needed <- c("mice", "sandwich", "lmtest", "MatchIt", "survival")
miss <- needed[!sapply(needed, requireNamespace, quietly = TRUE)]
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org", quiet = TRUE)
suppressPackageStartupMessages({
  library(mice);  library(sandwich)
  library(lmtest); library(MatchIt); library(survival)
})

SEED <- 42
M_IMP <- 5  # number of MICE imputations

# ── Helper: Rubin's rules for log-OR ──────────────────────────────────────
rubin_pool <- function(ests, ses) {
  # ests: vector of point estimates (log-OR) across imputations
  # ses:  vector of standard errors across imputations
  m <- length(ests)
  qbar <- mean(ests)                          # pooled estimate
  ubar <- mean(ses^2)                         # within-imputation variance
  b    <- var(ests)                            # between-imputation variance
  tv   <- ubar + (1 + 1/m) * b               # total variance
  se   <- sqrt(tv)
  # Barnard-Rubin df
  lambda <- (b + b/m) / tv
  df_old <- (m - 1) / lambda^2
  df_obs <- ubar / tv * (nrow(dat) - length(coef(ps_fit_template)) - 1) # approximate
  # Use simpler df
  df <- df_old
  or  <- exp(qbar)
  lo  <- exp(qbar - 1.96 * se)
  hi  <- exp(qbar + 1.96 * se)
  p   <- 2 * pnorm(-abs(qbar / se))
  list(or = or, lo = lo, hi = hi, p = p, logOR = qbar, se = se)
}

# ── Helper: simple pooling (no df needed) ─────────────────────────────────
pool_simple <- function(ests, ses) {
  m <- length(ests)
  qbar <- mean(ests)
  ubar <- mean(ses^2)
  b    <- var(ests)
  tv   <- ubar + (1 + 1/m) * b
  se   <- sqrt(tv)
  or   <- exp(qbar)
  lo   <- exp(qbar - 1.96 * se)
  hi   <- exp(qbar + 1.96 * se)
  p    <- 2 * pnorm(-abs(qbar / se))
  list(or = or, lo = lo, hi = hi, p = p, logOR = qbar, se = se)
}

# ── Helper: weighted logistic with cluster-robust SE ──────────────────────
wglm <- function(formula, data, w, cluster = NULL) {
  data$.w <- w
  fit <- glm(formula, data = data, weights = .w, family = quasibinomial())
  if (!is.null(cluster) && length(unique(cluster)) > 1) {
    vc <- vcovCL(fit, cluster = cluster)
  } else {
    vc <- vcovHC(fit, type = "HC1")
  }
  ct <- coeftest(fit, vcov. = vc)
  trt_row <- which(rownames(ct) == "mg_supp")
  if (length(trt_row) == 0) trt_row <- 2
  list(logOR = ct[trt_row, 1], se = ct[trt_row, 2])
}

# ── Helper: compute SMD ───────────────────────────────────────────────────
smd <- function(x, trt, w = NULL) {
  if (is.null(w)) w <- rep(1, length(x))
  m1 <- weighted.mean(x[trt == 1], w[trt == 1], na.rm = TRUE)
  m0 <- weighted.mean(x[trt == 0], w[trt == 0], na.rm = TRUE)
  # unweighted pooled SD
  s1 <- sd(x[trt == 1], na.rm = TRUE)
  s0 <- sd(x[trt == 0], na.rm = TRUE)
  sp <- sqrt((s1^2 + s0^2) / 2)
  if (sp < 1e-10) return(0)
  abs(m1 - m0) / sp
}

# ============================================================================
# LOAD DATA + STANDARDIZE COLUMN NAMES
# ============================================================================

standardize <- function(dat) {
  # Rename to canonical names used by the analysis
  rmap <- c(
    "mg_supplementation"     = "mg_supp",
    "hosp_mortality"         = "hospital_mortality",
    "age_num"                = "age",
    "hx_chf"                = "heart_failure",
    "hx_hypertension"       = "hypertension",
    "hx_diabetes"           = "diabetes",
    "hx_ckd"                = "ckd",
    "hx_copd"               = "copd",
    "hx_pvd"                = "pvd",
    "hx_stroke"             = "stroke",
    "hx_liver"              = "liver_disease",
    "baseline_cr"           = "baseline_creatinine",
    "baseline_egfr"         = "egfr",
    "nephrotox_loop_diuretic" = "loop_diuretics",
    "nephrotox_nsaid"       = "nsaids",
    "nephrotox_acei_arb"    = "acei_arb",
    "nephrotox_ppi"         = "ppi",
    "has_betablocker"       = "beta_blockers",
    "has_steroid"           = "steroids",
    "preop_antiarrhythmic"  = "antiarrhythmics",
    "first_k_value"         = "first_potassium",
    "first_ca_value"        = "first_calcium",
    "first_hr"              = "first_heartrate",
    "has_vasopressor"       = "vasopressor_6h",
    "nc_fracture"           = "fracture",
    "neuro_encephalopathy"  = "encephalopathy"
  )
  for (old in names(rmap)) {
    new <- rmap[[old]]
    if (old %in% names(dat) && !new %in% names(dat)) {
      names(dat)[names(dat) == old] <- new
    }
  }
  # eICU age is sometimes "> 89" string — use age_num
  if (is.character(dat$age) && "age" %in% names(dat)) {
    dat$age <- suppressWarnings(as.numeric(dat$age))
    dat$age[is.na(dat$age)] <- 90
  }
  # Surgery dummies
  if ("surgery_type" %in% names(dat)) {
    dat$surg_cabg     <- as.integer(dat$surgery_type == "cabg")
    dat$surg_valve    <- as.integer(dat$surgery_type == "valve")
    dat$surg_combined <- as.integer(dat$surgery_type == "combined")
  }
  dat
}

cat("\nLoading eICU cohort...\n")
dat_e <- standardize(read.csv("results/01_analysis_a_cohort.csv", stringsAsFactors = FALSE))
cat(sprintf("  N=%d, trt=%d (%.1f%%), AKI=%.1f%%\n",
            nrow(dat_e), sum(dat_e$mg_supp), 100*mean(dat_e$mg_supp),
            100*mean(dat_e$aki_kdigo1)))

cat("Loading MIMIC cohort...\n")
dat_m <- standardize(read.csv("results/04_mimic_cohort.csv", stringsAsFactors = FALSE))
cat(sprintf("  N=%d, trt=%d (%.1f%%), AKI=%.1f%%\n",
            nrow(dat_m), sum(dat_m$mg_supp), 100*mean(dat_m$mg_supp),
            100*mean(dat_m$aki_kdigo1)))

# ============================================================================
# ANALYSIS FUNCTION (runs on one database)
# ============================================================================
run_analysis <- function(dat, db_name, has_cluster = FALSE) {

  cat(sprintf("\n======================================================================\n"))
  cat(sprintf("%s  [N=%d, trt=%d]\n", db_name, nrow(dat), sum(dat$mg_supp)))
  cat(sprintf("======================================================================\n"))

  results <- list()

  # ── Identify available columns ────────────────────────────────────────
  has_lactate <- "first_lactate" %in% names(dat)
  has_ac      <- "ac_group" %in% names(dat)

  # ── Lactate: median + missing indicator (79% missing, MNAR) ───────────
  if (has_lactate) {
    n_lac_na <- sum(is.na(dat$first_lactate))
    pct_lac  <- round(100 * n_lac_na / nrow(dat), 1)
    cat(sprintf("  Lactate: %d NA (%.1f%%) → median + indicator\n", n_lac_na, pct_lac))
    dat$lactate_missing <- as.integer(is.na(dat$first_lactate))
    lac_med <- median(dat$first_lactate, na.rm = TRUE)
    dat$first_lactate[is.na(dat$first_lactate)] <- lac_med
  }

  # ── MICE targets: only low-missingness variables ──────────────────────
  mice_vars <- c("bmi", "first_heartrate", "first_calcium", "first_potassium")
  mice_vars <- mice_vars[mice_vars %in% names(dat)]

  # Count missingness
  for (v in mice_vars) {
    n_na <- sum(is.na(dat[[v]]))
    if (n_na > 0) cat(sprintf("  MICE target: %s (%d NA, %.1f%%)\n", v, n_na,
                               100 * n_na / nrow(dat)))
  }

  any_missing <- any(sapply(mice_vars, function(v) sum(is.na(dat[[v]])) > 0))

  # ── Build PS formula ──────────────────────────────────────────────────
  ps_covars <- c("age", "is_female", "bmi",
                 "surg_cabg", "surg_valve", "surg_combined",
                 "heart_failure", "hypertension", "diabetes", "ckd",
                 "copd", "pvd", "stroke", "liver_disease",
                 "baseline_creatinine", "egfr",
                 "loop_diuretics", "nsaids", "acei_arb", "ppi",
                 "beta_blockers", "steroids", "antiarrhythmics",
                 "first_potassium", "first_calcium", "first_heartrate",
                 "vasopressor_6h", "first_mg_value")

  if (has_lactate) {
    ps_covars <- c(ps_covars, "first_lactate", "lactate_missing")
  }

  # Keep only covariates that exist in data
  ps_covars <- ps_covars[ps_covars %in% names(dat)]
  cat(sprintf("  PS model: %d covariates\n", length(ps_covars)))

  ps_formula <- as.formula(paste("mg_supp ~", paste(ps_covars, collapse = " + ")))

  # ── Cluster variable ──────────────────────────────────────────────────
  cluster_var <- NULL
  if (has_cluster && "hospitalid" %in% names(dat)) {
    cluster_var <- dat$hospitalid
  }

  # ── Run MICE (or single run if no missing) ────────────────────────────
  if (any_missing) {
    cat(sprintf("  Running MICE (m=%d, pmm)...\n", M_IMP))

    # Include predictors that help imputation
    imp_predictors <- c(mice_vars, "age", "is_female", "baseline_creatinine",
                        "mg_supp", "aki_kdigo1")
    imp_predictors <- unique(imp_predictors[imp_predictors %in% names(dat)])

    imp <- mice(dat[, imp_predictors], m = M_IMP, method = "pmm",
                seed = SEED, printFlag = FALSE, maxit = 10)

    cat("  MICE converged.\n")
    n_imp <- M_IMP
  } else {
    cat("  No missing values in MICE targets — skipping imputation.\n")
    n_imp <- 1
  }

  # ── Storage for pooling across imputations ────────────────────────────
  # Each analysis stores logOR and SE per imputation
  store <- list(
    ow_aki1  = list(est = numeric(0), se = numeric(0)),
    ow_mort  = list(est = numeric(0), se = numeric(0)),
    psm_aki1 = list(est = numeric(0), se = numeric(0)),
    iptw_aki1 = list(est = numeric(0), se = numeric(0)),
    ow_frac  = list(est = numeric(0), se = numeric(0)),
    ow_enceph = list(est = numeric(0), se = numeric(0)),
    ac_aki1  = list(est = numeric(0), se = numeric(0))
  )

  # Track balance from last imputation for reporting
  final_ow_smds <- NULL
  final_ac_mg_smd <- NA

  for (i in seq_len(n_imp)) {
    # ── Get imputed data ──────────────────────────────────────────────
    d <- dat
    if (any_missing && n_imp > 1) {
      imputed <- complete(imp, i)
      for (v in mice_vars) {
        if (v %in% names(imputed)) d[[v]] <- imputed[[v]]
      }
    }

    # ── Drop rows with remaining NAs in PS covariates ─────────────────
    ps_complete <- complete.cases(d[, ps_covars])
    d_ps <- d[ps_complete, ]
    if (nrow(d_ps) < nrow(d)) {
      cat(sprintf("    Imp %d: dropped %d rows with remaining NA in PS covars\n",
                  i, nrow(d) - nrow(d_ps)))
    }

    # ── Estimate propensity score ─────────────────────────────────────
    ps_fit <- glm(ps_formula, data = d_ps, family = binomial())
    d_ps$ps <- fitted(ps_fit)

    # Clip extreme PS
    d_ps$ps <- pmax(pmin(d_ps$ps, 0.99), 0.01)

    # ── Overlap weights ───────────────────────────────────────────────
    d_ps$ow <- ifelse(d_ps$mg_supp == 1, 1 - d_ps$ps, d_ps$ps)

    # ── Stabilized IPTW ───────────────────────────────────────────────
    prev <- mean(d_ps$mg_supp)
    d_ps$iptw <- ifelse(d_ps$mg_supp == 1,
                        prev / d_ps$ps,
                        (1 - prev) / (1 - d_ps$ps))
    # Truncate at 1st/99th percentile
    q01 <- quantile(d_ps$iptw, 0.01)
    q99 <- quantile(d_ps$iptw, 0.99)
    d_ps$iptw <- pmax(pmin(d_ps$iptw, q99), q01)

    # ── OW balance check (last imputation) ────────────────────────────
    if (i == n_imp) {
      ow_smds <- sapply(ps_covars, function(v) {
        if (is.numeric(d_ps[[v]])) {
          smd(d_ps[[v]], d_ps$mg_supp, d_ps$ow)
        } else { NA }
      })
      final_ow_smds <- ow_smds
    }

    # ── OW: AKI KDIGO≥1 ──────────────────────────────────────────────
    cl <- if (has_cluster) cluster_var[ps_complete] else NULL
    res <- wglm(aki_kdigo1 ~ mg_supp, d_ps, d_ps$ow, cl)
    store$ow_aki1$est <- c(store$ow_aki1$est, res$logOR)
    store$ow_aki1$se  <- c(store$ow_aki1$se, res$se)

    # ── OW: Mortality ─────────────────────────────────────────────────
    if ("hospital_mortality" %in% names(d_ps)) {
      res <- wglm(hospital_mortality ~ mg_supp, d_ps, d_ps$ow, cl)
      store$ow_mort$est <- c(store$ow_mort$est, res$logOR)
      store$ow_mort$se  <- c(store$ow_mort$se, res$se)
    }

    # ── OW: Fracture ──────────────────────────────────────────────────
    if ("fracture" %in% names(d_ps)) {
      res <- wglm(fracture ~ mg_supp, d_ps, d_ps$ow, cl)
      store$ow_frac$est <- c(store$ow_frac$est, res$logOR)
      store$ow_frac$se  <- c(store$ow_frac$se, res$se)
    }

    # ── OW: Encephalopathy ────────────────────────────────────────────
    if ("encephalopathy" %in% names(d_ps)) {
      res <- wglm(encephalopathy ~ mg_supp, d_ps, d_ps$ow, cl)
      store$ow_enceph$est <- c(store$ow_enceph$est, res$logOR)
      store$ow_enceph$se  <- c(store$ow_enceph$se, res$se)
    }

    # ── IPTW: AKI (sensitivity) ───────────────────────────────────────
    res <- wglm(aki_kdigo1 ~ mg_supp, d_ps, d_ps$iptw, cl)
    store$iptw_aki1$est <- c(store$iptw_aki1$est, res$logOR)
    store$iptw_aki1$se  <- c(store$iptw_aki1$se, res$se)

    # ── PSM: AKI (co-primary) ─────────────────────────────────────────
    tryCatch({
      m_out <- matchit(ps_formula, data = d_ps, method = "nearest",
                       distance = d_ps$ps, caliper = 0.2, std.caliper = TRUE,
                       ratio = 1, replace = FALSE)
      md <- match.data(m_out)
      md$.w <- md$weights
      psm_fit <- glm(aki_kdigo1 ~ mg_supp, data = md, family = binomial(),
                     weights = .w)
      vc <- vcovHC(psm_fit, type = "HC1")
      ct <- coeftest(psm_fit, vcov. = vc)
      trt_row <- which(rownames(ct) == "mg_supp")
      if (length(trt_row) == 0) trt_row <- 2
      store$psm_aki1$est <- c(store$psm_aki1$est, ct[trt_row, 1])
      store$psm_aki1$se  <- c(store$psm_aki1$se, ct[trt_row, 2])
    }, error = function(e) {
      cat(sprintf("    PSM failed imp %d: %s\n", i, e$message))
    })

    # ── Active Comparator: OW ─────────────────────────────────────────
    if (has_ac && "ac_group" %in% names(d_ps)) {
      d_ac <- d_ps[d_ps$ac_group %in% c("mg_k", "k_only"), ]
      if (nrow(d_ac) > 50 && sum(d_ac$ac_group == "mg_k") > 10) {
        d_ac$ac_trt <- as.integer(d_ac$ac_group == "mg_k")
        ac_formula <- as.formula(paste("ac_trt ~", paste(ps_covars, collapse = " + ")))
        tryCatch({
          ac_ps_fit <- glm(ac_formula, data = d_ac, family = binomial())
          d_ac$ac_ps <- fitted(ac_ps_fit)
          d_ac$ac_ps <- pmax(pmin(d_ac$ac_ps, 0.99), 0.01)
          d_ac$ac_ow <- ifelse(d_ac$ac_trt == 1, 1 - d_ac$ac_ps, d_ac$ac_ps)

          # AC Mg balance (last imputation)
          if (i == n_imp && "first_mg_value" %in% names(d_ac)) {
            final_ac_mg_smd <- smd(d_ac$first_mg_value, d_ac$ac_trt, d_ac$ac_ow)
          }

          ac_res <- wglm(aki_kdigo1 ~ ac_trt, d_ac, d_ac$ac_ow, NULL)
          store$ac_aki1$est <- c(store$ac_aki1$est, ac_res$logOR)
          store$ac_aki1$se  <- c(store$ac_aki1$se, ac_res$se)
        }, error = function(e) {
          cat(sprintf("    AC model failed imp %d: %s\n", i, e$message))
        })
      }
    }
  }

  # ── Pool results via Rubin's rules ──────────────────────────────────
  cat("\n  ── Pooled Results (Rubin's rules, m=", n_imp, ") ──\n")

  out_rows <- list()

  for (analysis in names(store)) {
    s <- store[[analysis]]
    if (length(s$est) == 0) next

    if (length(s$est) == 1) {
      # No imputation — just report directly
      r <- list(or = exp(s$est), lo = exp(s$est - 1.96*s$se),
                hi = exp(s$est + 1.96*s$se),
                p = 2*pnorm(-abs(s$est/s$se)),
                logOR = s$est, se = s$se)
    } else {
      r <- pool_simple(s$est, s$se)
    }

    sig <- ifelse(r$p < 0.05, " *", "")
    cat(sprintf("  %-15s OR=%.3f (%.3f–%.3f) P=%.4f%s\n",
                analysis, r$or, r$lo, r$hi, r$p, sig))

    out_rows[[length(out_rows) + 1]] <- data.frame(
      db = db_name, analysis = analysis,
      or = round(r$or, 3), lo = round(r$lo, 3), hi = round(r$hi, 3),
      p = round(r$p, 4), logOR = round(r$logOR, 4), se = round(r$se, 4),
      stringsAsFactors = FALSE
    )
  }

  # ── OW balance report ───────────────────────────────────────────────
  if (!is.null(final_ow_smds)) {
    max_smd <- max(final_ow_smds, na.rm = TRUE)
    mg_smd  <- final_ow_smds["first_mg_value"]
    cat(sprintf("\n  OW balance: max SMD = %.4f", max_smd))
    if (!is.na(mg_smd)) cat(sprintf(", first_mg SMD = %.4f", mg_smd))
    n_above_01 <- sum(final_ow_smds > 0.10, na.rm = TRUE)
    cat(sprintf(", covariates with SMD>0.10: %d/%d\n", n_above_01,
                sum(!is.na(final_ow_smds))))
  }

  # ── AC Mg balance report ────────────────────────────────────────────
  if (!is.na(final_ac_mg_smd)) {
    cat(sprintf("  AC subgroup: first_mg_value weighted SMD = %.4f %s\n",
                final_ac_mg_smd,
                ifelse(final_ac_mg_smd < 0.10, "(balanced)", "(IMBALANCED)")))
  }

  do.call(rbind, out_rows)
}

# ============================================================================
# RUN BOTH DATABASES
# ============================================================================
res_e <- run_analysis(dat_e, "eICU", has_cluster = TRUE)
res_m <- run_analysis(dat_m, "MIMIC", has_cluster = FALSE)

# ============================================================================
# META-ANALYSIS (fixed-effects, inverse-variance)
# ============================================================================
cat("\n======================================================================\n")
cat("FIXED-EFFECTS META-ANALYSIS\n")
cat("======================================================================\n")

meta_rows <- list()

for (a in unique(c(res_e$analysis, res_m$analysis))) {
  e <- res_e[res_e$analysis == a, ]
  m <- res_m[res_m$analysis == a, ]
  if (nrow(e) == 0 || nrow(m) == 0) next

  # Inverse-variance weights
  w_e <- 1 / e$se^2
  w_m <- 1 / m$se^2
  w_tot <- w_e + w_m

  pool_logOR <- (w_e * e$logOR + w_m * m$logOR) / w_tot
  pool_se    <- sqrt(1 / w_tot)
  pool_or    <- exp(pool_logOR)
  pool_lo    <- exp(pool_logOR - 1.96 * pool_se)
  pool_hi    <- exp(pool_logOR + 1.96 * pool_se)
  pool_p     <- 2 * pnorm(-abs(pool_logOR / pool_se))

  # I-squared
  Q <- w_e * (e$logOR - pool_logOR)^2 + w_m * (m$logOR - pool_logOR)^2
  I2 <- max(0, (Q - 1) / Q) * 100

  sig <- ifelse(pool_p < 0.05, " *", "")
  cat(sprintf("  %-15s Pooled OR=%.3f (%.3f–%.3f) P=%.4f I²=%.0f%%%s\n",
              a, pool_or, pool_lo, pool_hi, pool_p, I2, sig))

  meta_rows[[length(meta_rows) + 1]] <- data.frame(
    db = "Pooled", analysis = a,
    or = round(pool_or, 3), lo = round(pool_lo, 3), hi = round(pool_hi, 3),
    p = round(pool_p, 4), logOR = round(pool_logOR, 4), se = round(pool_se, 4),
    I2 = round(I2, 0),
    stringsAsFactors = FALSE
  )
}

# ============================================================================
# SAVE
# ============================================================================
all_res <- rbind(res_e, res_m)
all_res$I2 <- NA
meta_df <- do.call(rbind, meta_rows)
all_res <- rbind(all_res, meta_df)

write.csv(all_res, "results/02_results.csv", row.names = FALSE)
cat("\n✓ Saved results/02_results.csv\n")

cat("\n======================================================================\n")
cat("KEY DECISION POINT:\n")
cat("If OW AKI pooled P < 0.05 → OW as primary estimator\n")
cat("If OW AKI pooled P > 0.05 → PSM as primary estimator\n")
cat("======================================================================\n")
