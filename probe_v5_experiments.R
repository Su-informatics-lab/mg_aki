#!/usr/bin/env Rscript
# ============================================================================
# probe_v5_experiments.R — Pre-submission analysis probes (standalone)
#
#   A. BMI sensitivity: compare primary results with raw vs capped BMI
#   B. APACHE IV: add to PS model, compare with/without
#   C. Weighted SMD verification for AC population
#   D. Summary and pipeline fix recommendations
#
# Prereqs: run probe_bmi_apache.py first (generates probe_cohort_bmifixed.csv
#          and probe_apache_scores.csv)
#
# Outputs:
#   results/probe_bmi_sensitivity.csv
#   results/probe_apache_sensitivity.csv
#   results/probe_weighted_smds.csv
#
# Run on Tempest: Rscript probe_v5_experiments.R
# Does NOT modify any pipeline files.
# ============================================================================

suppressPackageStartupMessages({
  library(sandwich); library(lmtest); library(mice)
})

RESULTS <- path.expand("~/mg_aki/results")
SEED <- 42

# ── Shared: column standardization ────────────────────────────────────────
stdz <- function(d) {
  rmap <- c(
    mg_supplementation      = "mg_supp",
    hosp_mortality           = "hospital_mortality",
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
    has_vasopressor         = "vasopressor_6h",
    nc_fracture             = "fracture",
    neuro_encephalopathy    = "encephalopathy"
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
    d$first_lactate[is.na(d$first_lactate)] <- median(d$first_lactate, na.rm = TRUE)
  }
  d
}

# ── Shared: PS covariates ─────────────────────────────────────────────────
ps_covars_base <- c(
  "age", "is_female", "bmi",
  "surg_cabg", "surg_valve", "surg_combined",
  "heart_failure", "hypertension", "diabetes", "ckd",
  "copd", "pvd", "stroke", "liver_disease",
  "baseline_creatinine", "egfr",
  "loop_diuretics", "nsaids", "acei_arb", "ppi",
  "beta_blockers", "steroids", "antiarrhythmics",
  "first_potassium", "first_calcium", "first_heartrate",
  "vasopressor_6h", "first_mg_value"
)

# ── Shared: weighted GLM ──────────────────────────────────────────────────
wglm <- function(formula, data, w, cluster = NULL) {
  data$.w <- w
  fit <- glm(formula, data = data, weights = .w, family = quasibinomial())
  vc <- if (!is.null(cluster) && length(unique(cluster)) > 1)
    vcovCL(fit, cluster = cluster) else vcovHC(fit, type = "HC1")
  ct <- coeftest(fit, vcov. = vc)
  tr <- which(rownames(ct) == names(coef(fit))[2])
  if (length(tr) == 0) tr <- 2
  list(or = exp(ct[tr, 1]),
       lo = exp(ct[tr, 1] - 1.96 * ct[tr, 2]),
       hi = exp(ct[tr, 1] + 1.96 * ct[tr, 2]),
       p  = 2 * pnorm(-abs(ct[tr, 1] / ct[tr, 2])),
       logOR = ct[tr, 1], se = ct[tr, 2])
}

smd_w <- function(x, trt, w) {
  m1 <- weighted.mean(x[trt == 1], w[trt == 1], na.rm = TRUE)
  m0 <- weighted.mean(x[trt == 0], w[trt == 0], na.rm = TRUE)
  sp <- sqrt((var(x[trt == 1], na.rm = TRUE) +
              var(x[trt == 0], na.rm = TRUE)) / 2)
  if (sp < 1e-10) return(0)
  abs(m1 - m0) / sp
}

# ── Shared: run primary analyses on a given dataset ───────────────────────
# Returns data.frame with one row per analysis
run_primary <- function(dat, db_name, ps_covars, has_cluster = FALSE) {
  # Median-impute MICE targets
  for (v in c("bmi", "first_heartrate", "first_calcium", "first_potassium"))
    if (v %in% names(dat) && any(is.na(dat[[v]])))
      dat[[v]][is.na(dat[[v]])] <- median(dat[[v]], na.rm = TRUE)

  ps_covars <- intersect(ps_covars, names(dat))
  if ("first_lactate" %in% names(dat))
    ps_covars <- c(ps_covars, "first_lactate", "lactate_missing")

  d <- dat[complete.cases(dat[, ps_covars]), ]
  ps_fml <- as.formula(paste("mg_supp ~", paste(ps_covars, collapse = " + ")))
  ps_fit <- glm(ps_fml, data = d, family = binomial())
  d$ps <- pmax(pmin(fitted(ps_fit), 0.99), 0.01)

  cl <- if (has_cluster && "hospitalid" %in% names(d)) d$hospitalid else NULL
  results <- list()

  # All-patient OW
  d$ow <- ifelse(d$mg_supp == 1, 1 - d$ps, d$ps)
  r <- wglm(aki_kdigo1 ~ mg_supp, d, d$ow, cl)
  results[[1]] <- data.frame(db = db_name, analysis = "ow_aki1", n = nrow(d),
    or = round(r$or, 4), lo = round(r$lo, 4), hi = round(r$hi, 4),
    p = round(r$p, 5))

  # IPTW
  prev <- mean(d$mg_supp)
  d$iptw <- ifelse(d$mg_supp == 1, prev / d$ps, (1 - prev) / (1 - d$ps))
  q01 <- quantile(d$iptw, 0.01); q99 <- quantile(d$iptw, 0.99)
  d$iptw <- pmax(pmin(d$iptw, q99), q01)
  r <- wglm(aki_kdigo1 ~ mg_supp, d, d$iptw, cl)
  results[[2]] <- data.frame(db = db_name, analysis = "iptw_aki1", n = nrow(d),
    or = round(r$or, 4), lo = round(r$lo, 4), hi = round(r$hi, 4),
    p = round(r$p, 5))

  # AC OW
  if ("ac_group" %in% names(d)) {
    d_ac <- d[d$ac_group %in% c("mg_k", "k_only"), ]
    d_ac$ac_trt <- as.integer(d_ac$ac_group == "mg_k")
    if (sum(d_ac$ac_trt) >= 10 && sum(d_ac$ac_trt == 0) >= 10) {
      tryCatch({
        ac_fml <- as.formula(paste("ac_trt ~", paste(ps_covars, collapse = "+")))
        ac_fit <- glm(ac_fml, data = d_ac, family = binomial())
        d_ac$ac_ps <- pmax(pmin(fitted(ac_fit), 0.99), 0.01)
        d_ac$ac_ow <- ifelse(d_ac$ac_trt == 1, 1 - d_ac$ac_ps, d_ac$ac_ps)
        r <- wglm(aki_kdigo1 ~ ac_trt, d_ac, d_ac$ac_ow, NULL)
        results[[3]] <- data.frame(db = db_name, analysis = "ac_aki1",
          n = nrow(d_ac), or = round(r$or, 4), lo = round(r$lo, 4),
          hi = round(r$hi, 4), p = round(r$p, 5))
      }, error = function(e) {})
    }
  }

  do.call(rbind, results)
}

# ============================================================================
# SECTION A: BMI SENSITIVITY
# ============================================================================
cat(sprintf("%s\nA. BMI SENSITIVITY ANALYSIS\n%s\n", strrep("=", 65), strrep("=", 65)))

bmi_results <- list()

# ── Load original cohort ──────────────────────────────────────────────────
cat("Loading original eICU cohort...\n")
dat_orig <- stdz(read.csv(file.path(RESULTS, "01_analysis_a_cohort.csv"),
                           stringsAsFactors = FALSE))

# Audit original BMI
bmi_orig <- dat_orig$bmi[!is.na(dat_orig$bmi)]
cat(sprintf("  Original BMI: n=%d, range=%.1f–%.1f, mean=%.1f, median=%.1f\n",
            length(bmi_orig), min(bmi_orig), max(bmi_orig),
            mean(bmi_orig), median(bmi_orig)))
n_out <- sum(bmi_orig < 10 | bmi_orig > 80)
cat(sprintf("  Outside [10,80]: %d patients\n", n_out))

# ── Run with original BMI ────────────────────────────────────────────────
cat("\n  Running primary analyses with ORIGINAL BMI...\n")
res_orig <- run_primary(dat_orig, "eICU_bmi_original", ps_covars_base,
                         has_cluster = TRUE)
for (i in seq_len(nrow(res_orig))) {
  r <- res_orig[i, ]
  cat(sprintf("    %-12s OR=%.4f (%.4f–%.4f) P=%.5f\n",
              r$analysis, r$or, r$lo, r$hi, r$p))
}
bmi_results[[1]] <- res_orig

# ── Cap BMI at [10, 80] and re-run ───────────────────────────────────────
cat("\n  Running primary analyses with CAPPED BMI [10,80]...\n")
dat_cap <- dat_orig
dat_cap$bmi[!is.na(dat_cap$bmi) & (dat_cap$bmi < 10 | dat_cap$bmi > 80)] <- NA

bmi_cap <- dat_cap$bmi[!is.na(dat_cap$bmi)]
cat(sprintf("  Capped BMI: n=%d, range=%.1f–%.1f, mean=%.1f\n",
            length(bmi_cap), min(bmi_cap), max(bmi_cap), mean(bmi_cap)))

res_cap <- run_primary(dat_cap, "eICU_bmi_capped", ps_covars_base,
                        has_cluster = TRUE)
for (i in seq_len(nrow(res_cap))) {
  r <- res_cap[i, ]
  cat(sprintf("    %-12s OR=%.4f (%.4f–%.4f) P=%.5f\n",
              r$analysis, r$or, r$lo, r$hi, r$p))
}
bmi_results[[2]] <- res_cap

# ── Run WITHOUT BMI entirely (worst-case comparison) ─────────────────────
cat("\n  Running primary analyses WITHOUT BMI in PS model...\n")
ps_no_bmi <- setdiff(ps_covars_base, "bmi")
res_nobmi <- run_primary(dat_orig, "eICU_no_bmi", ps_no_bmi,
                          has_cluster = TRUE)
for (i in seq_len(nrow(res_nobmi))) {
  r <- res_nobmi[i, ]
  cat(sprintf("    %-12s OR=%.4f (%.4f–%.4f) P=%.5f\n",
              r$analysis, r$or, r$lo, r$hi, r$p))
}
bmi_results[[3]] <- res_nobmi

# ── Compare ──────────────────────────────────────────────────────────────
cat(sprintf("\n  ── BMI SENSITIVITY COMPARISON ──\n"))
cat(sprintf("  %-15s %-12s %-12s %-12s  |ΔOR|\n",
            "Analysis", "Original", "Capped", "No BMI"))
for (a in c("ow_aki1", "iptw_aki1", "ac_aki1")) {
  r_o <- res_orig[res_orig$analysis == a, ]
  r_c <- res_cap[res_cap$analysis == a, ]
  r_n <- res_nobmi[res_nobmi$analysis == a, ]
  if (nrow(r_o) == 0) next
  delta_cap  <- if (nrow(r_c) > 0) abs(r_o$or - r_c$or) else NA
  delta_nobmi <- if (nrow(r_n) > 0) abs(r_o$or - r_n$or) else NA
  cat(sprintf("  %-15s %-12.4f %-12s %-12s  %.4f / %.4f\n",
              a, r_o$or,
              if (nrow(r_c) > 0) sprintf("%.4f", r_c$or) else "---",
              if (nrow(r_n) > 0) sprintf("%.4f", r_n$or) else "---",
              ifelse(is.na(delta_cap), NA, delta_cap),
              ifelse(is.na(delta_nobmi), NA, delta_nobmi)))
}

bmi_all <- do.call(rbind, bmi_results)
write.csv(bmi_all, file.path(RESULTS, "probe_bmi_sensitivity.csv"),
          row.names = FALSE)
cat("  Saved: probe_bmi_sensitivity.csv\n")


# ============================================================================
# SECTION B: APACHE IV IN PS MODEL
# ============================================================================
cat(sprintf("\n%s\nB. APACHE IV IN PS MODEL\n%s\n", strrep("=", 65), strrep("=", 65)))

apache_path <- file.path(RESULTS, "probe_apache_scores.csv")
if (!file.exists(apache_path)) {
  cat("  probe_apache_scores.csv not found — run probe_bmi_apache.py first\n")
  cat("  SKIPPING section B\n")
} else {
  apache <- read.csv(apache_path, stringsAsFactors = FALSE)
  cat(sprintf("  APACHE scores loaded: %d patients\n", nrow(apache)))

  # Detect score column name
  score_col <- intersect(c("apachescore", "acutephysiologyscore"), names(apache))[1]
  if (is.na(score_col)) {
    cat("  ERROR: no score column found\n")
  } else {
    cat(sprintf("  Score column: %s\n", score_col))

    # Use capped-BMI cohort as base (the right version)
    dat_base <- dat_cap

    # Merge APACHE
    dat_apache <- merge(dat_base, apache[, c("patientunitstayid", score_col)],
                        by = "patientunitstayid", all.x = TRUE)
    names(dat_apache)[names(dat_apache) == score_col] <- "apache_score"

    n_avail <- sum(!is.na(dat_apache$apache_score))
    pct <- 100 * n_avail / nrow(dat_apache)
    cat(sprintf("  APACHE available: %d/%d (%.1f%%)\n", n_avail, nrow(dat_apache), pct))

    if (pct < 30) {
      cat("  Too sparse (<30%) — skipping APACHE probe\n")
    } else {
      # Median-impute APACHE for missing
      apache_median <- median(dat_apache$apache_score, na.rm = TRUE)
      dat_apache$apache_missing <- as.integer(is.na(dat_apache$apache_score))
      dat_apache$apache_score[is.na(dat_apache$apache_score)] <- apache_median
      cat(sprintf("  APACHE median: %.0f (used for imputation of %d missing)\n",
                  apache_median, sum(dat_apache$apache_missing)))

      # ── Without APACHE (baseline, using capped BMI) ────────────
      cat("\n  Primary analyses WITHOUT APACHE (capped BMI):\n")
      res_no_apache <- run_primary(dat_apache, "eICU_no_apache", ps_covars_base,
                                    has_cluster = TRUE)
      for (i in seq_len(nrow(res_no_apache))) {
        r <- res_no_apache[i, ]
        cat(sprintf("    %-12s OR=%.4f P=%.5f\n", r$analysis, r$or, r$p))
      }

      # ── With APACHE ────────────────────────────────────────────
      cat("\n  Primary analyses WITH APACHE (capped BMI):\n")
      ps_with_apache <- c(ps_covars_base, "apache_score", "apache_missing")
      res_with_apache <- run_primary(dat_apache, "eICU_with_apache",
                                      ps_with_apache, has_cluster = TRUE)
      for (i in seq_len(nrow(res_with_apache))) {
        r <- res_with_apache[i, ]
        cat(sprintf("    %-12s OR=%.4f P=%.5f\n", r$analysis, r$or, r$p))
      }

      # ── Compare ────────────────────────────────────────────────
      cat(sprintf("\n  ── APACHE SENSITIVITY COMPARISON ──\n"))
      cat(sprintf("  %-15s %-12s %-12s  |ΔOR|\n",
                  "Analysis", "Without", "With APACHE"))
      for (a in c("ow_aki1", "iptw_aki1", "ac_aki1")) {
        r_no  <- res_no_apache[res_no_apache$analysis == a, ]
        r_yes <- res_with_apache[res_with_apache$analysis == a, ]
        if (nrow(r_no) == 0 || nrow(r_yes) == 0) next
        delta <- abs(r_no$or - r_yes$or)
        pct_change <- 100 * delta / r_no$or
        cat(sprintf("  %-15s %-12.4f %-12.4f  %.4f (%.1f%%)\n",
                    a, r_no$or, r_yes$or, delta, pct_change))
      }

      apache_all <- rbind(res_no_apache, res_with_apache)
      write.csv(apache_all, file.path(RESULTS, "probe_apache_sensitivity.csv"),
                row.names = FALSE)
      cat("  Saved: probe_apache_sensitivity.csv\n")
    }
  }
}


# ============================================================================
# SECTION C: WEIGHTED SMD VERIFICATION (AC POPULATION)
# ============================================================================
cat(sprintf("\n%s\nC. WEIGHTED SMD VERIFICATION (AC POPULATION)\n%s\n",
            strrep("=", 65), strrep("=", 65)))

smd_results <- list()

for (db_info in list(
  list(file = "01_analysis_a_cohort.csv", name = "eICU", cluster = TRUE),
  list(file = "04_mimic_cohort.csv", name = "MIMIC", cluster = FALSE)
)) {
  path <- file.path(RESULTS, db_info$file)
  if (!file.exists(path)) next

  d <- stdz(read.csv(path, stringsAsFactors = FALSE))

  # Cap BMI for eICU
  if (db_info$name == "eICU")
    d$bmi[!is.na(d$bmi) & (d$bmi < 10 | d$bmi > 80)] <- NA

  ps_covars <- intersect(ps_covars_base, names(d))
  if ("first_lactate" %in% names(d))
    ps_covars <- c(ps_covars, "first_lactate", "lactate_missing")

  # Median-impute
  for (v in c("bmi", "first_heartrate", "first_calcium", "first_potassium"))
    if (v %in% names(d) && any(is.na(d[[v]])))
      d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm = TRUE)

  if (!"ac_group" %in% names(d)) next

  d_ac <- d[d$ac_group %in% c("mg_k", "k_only"), ]
  d_ac$ac_trt <- as.integer(d_ac$ac_group == "mg_k")
  d_ac <- d_ac[complete.cases(d_ac[, ps_covars]), ]

  if (sum(d_ac$ac_trt) < 10) next

  cat(sprintf("\n  %s AC: N=%d (trt=%d, ctrl=%d)\n",
              db_info$name, nrow(d_ac), sum(d_ac$ac_trt),
              sum(d_ac$ac_trt == 0)))

  # Fit AC PS
  ac_fml <- as.formula(paste("ac_trt ~", paste(ps_covars, collapse = "+")))
  tryCatch({
    ac_fit <- glm(ac_fml, data = d_ac, family = binomial())
    d_ac$ac_ps <- pmax(pmin(fitted(ac_fit), 0.99), 0.01)
    d_ac$ac_ow <- ifelse(d_ac$ac_trt == 1, 1 - d_ac$ac_ps, d_ac$ac_ps)

    # Compute unweighted and weighted SMDs
    cat(sprintf("  %-25s %8s %8s  %s\n", "Covariate", "Unwt SMD", "Wt SMD", "Status"))
    cat(sprintf("  %s\n", strrep("-", 60)))

    max_unwt <- 0; max_wt <- 0
    n_wt_above_01 <- 0; n_wt_above_001 <- 0

    for (v in ps_covars) {
      if (!is.numeric(d_ac[[v]])) next
      smd_unwt <- smd_w(d_ac[[v]], d_ac$ac_trt, rep(1, nrow(d_ac)))
      smd_wt   <- smd_w(d_ac[[v]], d_ac$ac_trt, d_ac$ac_ow)
      max_unwt <- max(max_unwt, smd_unwt)
      max_wt   <- max(max_wt, smd_wt)
      if (smd_wt > 0.10) n_wt_above_01 <- n_wt_above_01 + 1
      if (smd_wt > 0.01) n_wt_above_001 <- n_wt_above_001 + 1

      flag <- if (smd_wt > 0.10) "✗ IMBALANCED" else
              if (smd_wt > 0.01) "~ marginal" else "✓"
      # Only print if notable
      if (smd_unwt > 0.05 || smd_wt > 0.005) {
        cat(sprintf("  %-25s %8.4f %8.4f  %s\n", v, smd_unwt, smd_wt, flag))
      }

      smd_results[[length(smd_results) + 1]] <- data.frame(
        db = db_info$name, covariate = v,
        smd_unweighted = round(smd_unwt, 5),
        smd_weighted = round(smd_wt, 5))
    }

    cat(sprintf("\n  Summary:\n"))
    cat(sprintf("    Max unweighted SMD: %.4f\n", max_unwt))
    cat(sprintf("    Max weighted SMD:   %.4f\n", max_wt))
    cat(sprintf("    Covariates with weighted SMD > 0.10: %d\n", n_wt_above_01))
    cat(sprintf("    Covariates with weighted SMD > 0.01: %d\n", n_wt_above_001))

    # ── VERDICT ──────────────────────────────────────────────────
    if (max_wt < 0.01) {
      cat(sprintf("    ✓ VERIFIED: all weighted SMDs < 0.01\n"))
      cat(sprintf("    Manuscript claim 'all SMDs < 0.01' is SUPPORTED\n"))
    } else if (max_wt < 0.10) {
      cat(sprintf("    ⚠ Some weighted SMDs between 0.01 and 0.10\n"))
      cat(sprintf("    Manuscript should say 'all SMDs < 0.10' or identify exceptions\n"))
    } else {
      cat(sprintf("    ✗ WARNING: weighted SMD > 0.10 detected\n"))
      cat(sprintf("    Manuscript claim needs revision\n"))
    }
  }, error = function(e) {
    cat(sprintf("  AC PS fitting failed: %s\n", e$message))
  })
}

if (length(smd_results) > 0) {
  smd_df <- do.call(rbind, smd_results)
  write.csv(smd_df, file.path(RESULTS, "probe_weighted_smds.csv"),
            row.names = FALSE)
  cat("\n  Saved: probe_weighted_smds.csv\n")
}


# ============================================================================
# SECTION D: SUMMARY AND PIPELINE FIX RECOMMENDATIONS
# ============================================================================
cat(sprintf("\n%s\nD. SUMMARY AND RECOMMENDATIONS\n%s\n",
            strrep("=", 65), strrep("=", 65)))

cat("
  1. BMI FIX
     If |ΔOR| between original and capped BMI is < 0.02:
       → BMI outliers had minimal impact on results
       → Still fix the ETL (add bmi.between(10,80) filter) for correctness
       → No need to re-run full pipeline — just fix ETL + re-run steps 1-7
     If |ΔOR| > 0.02:
       → BMI outliers materially affected PS estimates
       → Must fix ETL and re-run ENTIRE pipeline
       → Report capped-BMI results as primary

  2. APACHE
     If |ΔOR| between with/without APACHE is < 0.03:
       → APACHE adds minimal information beyond existing covariates
       → Consider adding as sensitivity analysis in supplement
       → Report: 'Adding APACHE IV score did not change estimates (ΔOR < X)'
     If |ΔOR| > 0.05:
       → APACHE captures important confounding not in current PS
       → Must add to primary PS model
       → Re-run entire pipeline with APACHE in ps_covars

  3. WEIGHTED SMDs
     If all < 0.01: manuscript claim verified
     If some > 0.01: adjust manuscript language

  4. PIPELINE FIX CHECKLIST (if probes pass):
     a. 01_etl.py: add BMI cap in eICU section (one line)
     b. Decide on APACHE based on probe results
     c. Re-run: bash run.sh 1 2 3 4 5 7 8 9 fig
     d. Verify 07d CSV shows clean BMI
     e. Update eTable 3 with new 07d output
")

cat(sprintf("\n%s\nprobe_v5_experiments.R COMPLETE\n%s\n",
            strrep("=", 65), strrep("=", 65)))
