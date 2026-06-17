#!/usr/bin/env Rscript
# ============================================================================
# 07_sensitivity.R — Sensitivity & supplementary analyses
#
# Sections:
#   A. E-value sensitivity analysis        → results/07_evalues.csv
#   B. Prognostic Mg→AKI by severity       → results/07b_prognostic.csv
#   C. MICE stability (m=10, 20)           → results/07c_mice_stability.csv
#   D. AC baseline characteristics          → results/07d_ac_table1.csv
#
# Run: Rscript 07_sensitivity.R            # all sections
#      Rscript 07_sensitivity.R A B        # specific sections
# ============================================================================

suppressPackageStartupMessages({
  library(sandwich); library(lmtest); library(tableone); library(mice)
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
  # Lactate: median + missing indicator
  if ("first_lactate" %in% names(d)) {
    d$lactate_missing <- as.integer(is.na(d$first_lactate))
    d$first_lactate[is.na(d$first_lactate)] <- median(d$first_lactate, na.rm = TRUE)
  }
  # Median-impute MICE targets
  for (v in c("bmi", "first_heartrate", "first_calcium", "first_potassium"))
    if (v %in% names(d) && any(is.na(d[[v]])))
      d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm = TRUE)
  d
}

# ── Shared: PS covariate list (◆ 31 covariates) ──────────────────────────
ps_covars_full <- c(
  "age", "is_female", "bmi",
  "surg_cabg", "surg_valve", "surg_combined",
  "heart_failure", "hypertension", "diabetes", "ckd",
  "copd", "pvd", "stroke", "liver_disease",
  "baseline_creatinine", "egfr",
  "loop_diuretics", "nsaids", "acei_arb", "ppi",
  "beta_blockers", "steroids", "antiarrhythmics",
  "first_potassium", "first_calcium", "first_heartrate",
  "vasopressor_6h", "transfusion_6h", "first_mg_value"
)

# ── Shared: weighted GLM ──────────────────────────────────────────────────
wglm <- function(formula, data, w, cluster = NULL) {
  data$.w <- w
  fit <- glm(formula, data = data, weights = .w, family = quasibinomial())
  vc <- if (!is.null(cluster) && length(unique(cluster)) > 1)
    vcovCL(fit, cluster = cluster) else vcovHC(fit, type = "HC1")
  ct <- coeftest(fit, vcov. = vc)
  tr <- which(rownames(ct) == "mg_supp")
  if (length(tr) == 0) tr <- 2
  list(logOR = ct[tr, 1], se = ct[tr, 2])
}

pool_simple <- function(ests, ses) {
  m <- length(ests); qbar <- mean(ests); ubar <- mean(ses^2)
  b <- var(ests); tv <- ubar + (1 + 1/m) * b; se <- sqrt(tv)
  list(or = exp(qbar), lo = exp(qbar - 1.96 * se),
       hi = exp(qbar + 1.96 * se), p = 2 * pnorm(-abs(qbar / se)))
}

# ── Shared: load both cohorts ─────────────────────────────────────────────
cat("Loading cohorts...\n")
dat_e <- stdz(read.csv(file.path(RESULTS, "01_analysis_a_cohort.csv"),
                        stringsAsFactors = FALSE))
dat_m <- stdz(read.csv(file.path(RESULTS, "04_mimic_cohort.csv"),
                        stringsAsFactors = FALSE))
cat(sprintf("  eICU: N=%d  MIMIC: N=%d\n", nrow(dat_e), nrow(dat_m)))

# ── Parse section arguments ───────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
run_all <- length(args) == 0
sections <- if (run_all) c("A", "B", "C", "D") else toupper(args)

# ============================================================================
# SECTION A: E-VALUE SENSITIVITY ANALYSIS
# ============================================================================
if ("A" %in% sections) {
  cat(sprintf("\n%s\nA. E-VALUE SENSITIVITY ANALYSIS\n%s\n",
              strrep("=", 65), strrep("=", 65)))

  e_value <- function(rr) {
    if (is.na(rr) || rr <= 1.0) return(1.0)
    rr + sqrt(rr * (rr - 1))
  }
  or_to_rr <- function(or_val, p0) {
    or_val / (1 - p0 + p0 * or_val)
  }

  res <- read.csv(file.path(RESULTS, "02_results.csv"), stringsAsFactors = FALSE)
  or_col <- if ("or" %in% names(res)) "or" else "or_"

  p0_e <- mean(dat_e$aki_kdigo1[dat_e$mg_supp == 0], na.rm = TRUE)
  p0_m <- mean(dat_m$aki_kdigo1[dat_m$mg_supp == 0], na.rm = TRUE)
  p0_ac_e <- if ("ac_group" %in% names(dat_e))
    mean(dat_e$aki_kdigo1[dat_e$ac_group == "k_only"], na.rm = TRUE) else NA
  p0_ac_m <- if ("ac_group" %in% names(dat_m))
    mean(dat_m$aki_kdigo1[dat_m$ac_group == "k_only"], na.rm = TRUE) else NA

  cat(sprintf("  Prevalence: eICU all=%.3f, MIMIC all=%.3f\n", p0_e, p0_m))
  cat(sprintf("  Prevalence: eICU AC=%.3f, MIMIC AC=%.3f\n",
              ifelse(is.na(p0_ac_e), NA, p0_ac_e),
              ifelse(is.na(p0_ac_m), NA, p0_ac_m)))

  pooled <- res[res$db == "Pooled", ]
  p0_all <- mean(c(p0_e, p0_m), na.rm = TRUE)
  p0_ac  <- mean(c(p0_ac_e, p0_ac_m), na.rm = TRUE)

  analyses <- list(
    list(key = "ac_aki1",   label = "Primary: AC AKI",       p0 = p0_ac),
    list(key = "iptw_aki1", label = "Sensitivity: IPTW AKI", p0 = p0_all),
    list(key = "ow_aki1",   label = "Sensitivity: OW AKI",   p0 = p0_all),
    list(key = "ow_enceph", label = "Exploratory: Enceph",   p0 = NA),
    list(key = "ow_frac",   label = "Control: Fracture",     p0 = NA)
  )

  mg_strat_path <- file.path(RESULTS, "08_mg_stratified.csv")
  if (file.exists(mg_strat_path)) {
    mg_strat <- read.csv(mg_strat_path, stringsAsFactors = FALSE)
    mg_or_col <- if ("or" %in% names(mg_strat)) "or" else "or_"
    gt23 <- mg_strat[mg_strat$db == "eICU" & mg_strat$stratum == ">2.3" &
                     mg_strat$analysis == "all_ow", ]
    gt23_ac <- mg_strat[mg_strat$db == "eICU" & mg_strat$stratum == ">2.3" &
                        mg_strat$analysis == "ac_ow", ]
    if (nrow(gt23) > 0) {
      analyses[[length(analyses) + 1]] <- list(
        key = "mg_gt23_ow", label = "Mg >2.3 OW", p0 = p0_all,
        custom_or = gt23[[mg_or_col]][1], custom_lo = gt23$lo[1],
        custom_hi = gt23$hi[1])
    }
    if (nrow(gt23_ac) > 0) {
      analyses[[length(analyses) + 1]] <- list(
        key = "mg_gt23_ac", label = "Mg >2.3 AC", p0 = p0_ac,
        custom_or = gt23_ac[[mg_or_col]][1], custom_lo = gt23_ac$lo[1],
        custom_hi = gt23_ac$hi[1])
    }
  }

  ev_rows <- list()
  for (a in analyses) {
    if (!is.null(a$custom_or)) {
      or_est <- a$custom_or; or_lo <- a$custom_lo; or_hi <- a$custom_hi
    } else {
      r <- pooled[pooled$analysis == a$key, ]
      if (nrow(r) == 0) next
      or_est <- r[[or_col]][1]; or_lo <- r$lo[1]; or_hi <- r$hi[1]
    }

    ci_null <- if (or_est < 1) or_hi else or_lo
    rr_pt <- if (or_est < 1) 1 / or_est else or_est
    rr_ci <- if (ci_null < 1) 1 / ci_null else ci_null
    e_pt <- e_value(rr_pt)
    e_ci <- e_value(rr_ci)

    cat(sprintf("  %-25s OR=%.3f (%.3f-%.3f)  E=%.2f (CI limit %.2f)\n",
                a$label, or_est, or_lo, or_hi, e_pt, e_ci))

    ev_rows[[length(ev_rows) + 1]] <- data.frame(
      analysis = a$key, label = a$label,
      or = round(or_est, 3), lo = round(or_lo, 3), hi = round(or_hi, 3),
      e_point = round(e_pt, 2), e_ci_limit = round(e_ci, 2),
      stringsAsFactors = FALSE)
  }

  ev_df <- do.call(rbind, ev_rows)
  write.csv(ev_df, file.path(RESULTS, "07_evalues.csv"), row.names = FALSE)
  cat("  Saved: 07_evalues.csv\n")
}

# ============================================================================
# SECTION B: PROGNOSTIC Mg→AKI BY SEVERITY
# ============================================================================
if ("B" %in% sections) {
  cat(sprintf("\n%s\nB. PROGNOSTIC Mg→AKI ASSOCIATION\n%s\n",
              strrep("=", 65), strrep("=", 65)))

  adj_vars <- c("age", "is_female", "bmi",
    "surg_cabg", "surg_valve", "surg_combined",
    "heart_failure", "hypertension", "diabetes", "ckd",
    "copd", "pvd", "stroke", "liver_disease",
    "baseline_creatinine", "egfr",
    "loop_diuretics", "nsaids", "acei_arb", "ppi",
    "beta_blockers", "steroids", "antiarrhythmics",
    "first_potassium", "first_calcium", "first_heartrate",
    "vasopressor_6h", "transfusion_6h")

  aki_defs <- c(
    "KDIGO stage >=1"    = "aki_kdigo1",
    "Cr ratio >=1.5x"    = "aki_primary",
    "Delta Cr >=0.3/48h" = "aki_delta03",
    "KDIGO stage >=2"    = "aki_stage2",
    "KDIGO stage >=3"    = "aki_stage3")

  run_prognostic <- function(d, db_name) {
    cat(sprintf("\n  %s (N=%d)\n", db_name, nrow(d)))
    avail <- intersect(adj_vars, names(d))
    if ("first_lactate" %in% names(d))
      avail <- c(avail, "first_lactate", "lactate_missing")

    rows <- list()
    for (lbl in names(aki_defs)) {
      outcome <- aki_defs[[lbl]]
      if (!outcome %in% names(d)) next
      fml <- as.formula(paste(outcome, "~ first_mg_value +",
                              paste(avail, collapse = "+")))
      d_cc <- d[complete.cases(d[, c("first_mg_value", avail, outcome)]), ]
      tryCatch({
        fit <- glm(fml, data = d_cc, family = binomial())
        ct <- coeftest(fit, vcov. = vcovHC(fit, type = "HC1"))
        mg_row <- which(rownames(ct) == "first_mg_value")
        or <- exp(ct[mg_row, 1])
        lo <- exp(ct[mg_row, 1] - 1.96 * ct[mg_row, 2])
        hi <- exp(ct[mg_row, 1] + 1.96 * ct[mg_row, 2])
        p  <- 2 * pnorm(-abs(ct[mg_row, 1] / ct[mg_row, 2]))
        sig <- ifelse(p < 0.05, " *", "")
        cat(sprintf("    %-25s OR %.3f (%.3f-%.3f) P=%.4f%s\n",
                    lbl, or, lo, hi, p, sig))
        rows[[length(rows) + 1]] <- data.frame(
          db = db_name, outcome = lbl, surgery = "All",
          or = round(or, 3), lo = round(lo, 3), hi = round(hi, 3),
          p = round(p, 4))
      }, error = function(e) cat(sprintf("    %-25s FAILED: %s\n",
                                          lbl, e$message)))
    }

    if ("aki_kdigo1" %in% names(d) && "surgery_type" %in% names(d)) {
      d$complex <- as.integer(d$surgery_type %in% c("valve", "combined"))
      for (stype in c("Simple", "Complex")) {
        sub <- if (stype == "Simple") d[d$complex == 0, ] else d[d$complex == 1, ]
        if (nrow(sub) < 50) next
        fml <- as.formula(paste("aki_kdigo1 ~ first_mg_value +",
                                paste(avail, collapse = "+")))
        sub_cc <- sub[complete.cases(sub[, c("first_mg_value", avail,
                                             "aki_kdigo1")]), ]
        tryCatch({
          fit <- glm(fml, data = sub_cc, family = binomial())
          ct <- coeftest(fit, vcov. = vcovHC(fit, type = "HC1"))
          mg_row <- which(rownames(ct) == "first_mg_value")
          or <- exp(ct[mg_row, 1])
          lo <- exp(ct[mg_row, 1] - 1.96 * ct[mg_row, 2])
          hi <- exp(ct[mg_row, 1] + 1.96 * ct[mg_row, 2])
          p  <- 2 * pnorm(-abs(ct[mg_row, 1] / ct[mg_row, 2]))
          cat(sprintf("    %-25s OR %.3f (%.3f-%.3f) P=%.4f\n",
                      stype, or, lo, hi, p))
          rows[[length(rows) + 1]] <- data.frame(
            db = db_name, outcome = "KDIGO stage >=1", surgery = stype,
            or = round(or, 3), lo = round(lo, 3), hi = round(hi, 3),
            p = round(p, 4))
        }, error = function(e) cat(sprintf("    %-25s FAILED\n", stype)))
      }
    }
    do.call(rbind, rows)
  }

  res_b <- rbind(run_prognostic(dat_e, "eICU"),
                 run_prognostic(dat_m, "MIMIC"))
  write.csv(res_b, file.path(RESULTS, "07b_prognostic.csv"), row.names = FALSE)
  cat(sprintf("  Saved: 07b_prognostic.csv (%d rows)\n", nrow(res_b)))
}

# ============================================================================
# SECTION C: MICE STABILITY (m=10, 20)
# ============================================================================
if ("C" %in% sections) {
  cat(sprintf("\n%s\nC. MICE STABILITY (m=10, 20)\n%s\n",
              strrep("=", 65), strrep("=", 65)))

  run_m <- function(dat, m_val, db_name, has_cluster) {
    cat(sprintf("  %s m=%d...\n", db_name, m_val))
    ps_covars <- intersect(ps_covars_full, names(dat))
    if ("first_lactate" %in% names(dat))
      ps_covars <- c(ps_covars, "first_lactate", "lactate_missing")

    mice_vars <- intersect(c("bmi", "first_heartrate", "first_calcium",
                              "first_potassium"), names(dat))
    any_missing <- any(sapply(mice_vars,
                              function(v) sum(is.na(dat[[v]])) > 0))

    dat_raw <- read.csv(file.path(RESULTS,
      ifelse(db_name == "eICU", "01_analysis_a_cohort.csv",
             "04_mimic_cohort.csv")), stringsAsFactors = FALSE)
    dat_raw <- stdz(dat_raw)
    for (v in mice_vars)
      if (v %in% names(dat_raw)) dat[[v]] <- dat_raw[[v]][seq_len(nrow(dat))]
    any_missing <- any(sapply(mice_vars,
                              function(v) sum(is.na(dat[[v]])) > 0))

    imp_preds <- unique(c(mice_vars, "age", "is_female",
                           "baseline_creatinine", "mg_supp", "aki_kdigo1"))
    imp_preds <- imp_preds[imp_preds %in% names(dat)]

    if (any_missing) {
      imp <- mice(dat[, imp_preds], m = m_val, method = "pmm",
                  seed = SEED, printFlag = FALSE, maxit = 10)
    }

    ps_formula <- as.formula(paste("mg_supp ~",
                                    paste(ps_covars, collapse = "+")))

    results <- list()
    for (aname in c("ow_aki1", "iptw_aki1", "ac_aki1")) {
      ests <- ses <- numeric(0)
      for (i in seq_len(m_val)) {
        d <- dat
        if (any_missing) {
          imputed <- complete(imp, i)
          for (v in mice_vars)
            if (v %in% names(imputed)) d[[v]] <- imputed[[v]]
        }
        for (v in ps_covars)
          if (v %in% names(d) && any(is.na(d[[v]])))
            d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm = TRUE)
        d_ps <- d[complete.cases(d[, ps_covars]), ]
        ps_fit <- glm(ps_formula, data = d_ps, family = binomial())
        d_ps$ps <- pmax(pmin(fitted(ps_fit), 0.99), 0.01)

        if (aname == "ac_aki1" && "ac_group" %in% names(d_ps)) {
          d_ac <- d_ps[d_ps$ac_group %in% c("mg_k", "k_only"), ]
          d_ac$ac_trt <- as.integer(d_ac$ac_group == "mg_k")
          ac_fml <- as.formula(paste("ac_trt ~",
                                      paste(ps_covars, collapse = "+")))
          tryCatch({
            ac_fit <- glm(ac_fml, data = d_ac, family = binomial())
            d_ac$ac_ps <- pmax(pmin(fitted(ac_fit), 0.99), 0.01)
            d_ac$ac_ow <- ifelse(d_ac$ac_trt == 1,
                                  1 - d_ac$ac_ps, d_ac$ac_ps)
            d_ac$.w <- d_ac$ac_ow
            fit <- glm(aki_kdigo1 ~ ac_trt, data = d_ac, weights = .w,
                       family = quasibinomial())
            ct <- coeftest(fit, vcov. = vcovHC(fit, type = "HC1"))
            ests <- c(ests, ct[2, 1]); ses <- c(ses, ct[2, 2])
          }, error = function(e) {})
        } else {
          d_ps$ow <- ifelse(d_ps$mg_supp == 1, 1 - d_ps$ps, d_ps$ps)
          prev <- mean(d_ps$mg_supp)
          d_ps$iptw <- ifelse(d_ps$mg_supp == 1,
                               prev / d_ps$ps, (1 - prev) / (1 - d_ps$ps))
          q01 <- quantile(d_ps$iptw, 0.01)
          q99 <- quantile(d_ps$iptw, 0.99)
          d_ps$iptw <- pmax(pmin(d_ps$iptw, q99), q01)
          w <- if (aname == "iptw_aki1") d_ps$iptw else d_ps$ow
          cl <- if (has_cluster && "hospitalid" %in% names(d_ps))
            d_ps$hospitalid else NULL
          res <- wglm(aki_kdigo1 ~ mg_supp, d_ps, w, cl)
          ests <- c(ests, res$logOR); ses <- c(ses, res$se)
        }
      }
      if (length(ests) > 0) {
        r <- pool_simple(ests, ses)
        results[[length(results) + 1]] <- data.frame(
          db = db_name, analysis = aname, m = m_val,
          or = round(r$or, 3), lo = round(r$lo, 3),
          hi = round(r$hi, 3), p = round(r$p, 4))
        cat(sprintf("    %-12s OR=%.3f (%.3f-%.3f) P=%.4f\n",
                    aname, r$or, r$lo, r$hi, r$p))
      }
    }
    do.call(rbind, results)
  }

  all_c <- list()
  for (m_val in c(10, 20)) {
    all_c[[length(all_c) + 1]] <- run_m(dat_e, m_val, "eICU", TRUE)
    all_c[[length(all_c) + 1]] <- run_m(dat_m, m_val, "MIMIC", FALSE)
  }
  res_c <- do.call(rbind, all_c)

  res5 <- read.csv(file.path(RESULTS, "02_results.csv"), stringsAsFactors = FALSE)
  or5_col <- if ("or" %in% names(res5)) "or" else "or_"
  if (or5_col == "or_") names(res5)[names(res5) == "or_"] <- "or"
  m5 <- res5[res5$analysis %in% c("ow_aki1", "iptw_aki1", "ac_aki1") &
             res5$db != "Pooled", ]
  if (nrow(m5) > 0) {
    m5$m <- 5
    m5 <- m5[, intersect(names(res_c), names(m5))]
    res_c <- rbind(m5, res_c)
  }

  write.csv(res_c, file.path(RESULTS, "07c_mice_stability.csv"),
            row.names = FALSE)
  cat(sprintf("  Saved: 07c_mice_stability.csv (%d rows)\n", nrow(res_c)))
}

# ============================================================================
# SECTION D: AC BASELINE CHARACTERISTICS (eTable 3)
# ============================================================================
if ("D" %in% sections) {
  cat(sprintf("\n%s\nD. AC BASELINE CHARACTERISTICS (eTable 3)\n%s\n",
              strrep("=", 65), strrep("=", 65)))

  make_ac_table <- function(d, db_label) {
    if (!"ac_group" %in% names(d) && !"mg_supp" %in% names(d)) {
      cat(sprintf("  %s: no ac_group\n", db_label)); return(NULL)
    }
    if ("ac_group" %in% names(d)) {
      d_ac <- d[d$ac_group %in% c("mg_k", "k_only"), ]
      d_ac$trt_label <- ifelse(d_ac$ac_group == "mg_k", "Mg+K", "K-only")
    } else {
      return(NULL)
    }

    vars <- intersect(c("age", "is_female", "bmi", "surgery_type",
      "heart_failure", "hypertension", "diabetes", "ckd",
      "copd", "pvd", "stroke", "liver_disease",
      "baseline_creatinine", "egfr",
      "loop_diuretics", "nsaids", "acei_arb", "ppi",
      "beta_blockers", "steroids", "antiarrhythmics", "vasopressor_6h",
      "transfusion_6h",
      "first_mg_value", "first_potassium", "first_calcium",
      "first_heartrate",
      "aki_kdigo1", "hospital_mortality"), names(d_ac))

    cat_vars <- intersect(c("surgery_type", "is_female",
      "heart_failure", "hypertension", "diabetes", "ckd",
      "copd", "pvd", "stroke", "liver_disease",
      "loop_diuretics", "nsaids", "acei_arb", "ppi",
      "beta_blockers", "steroids", "antiarrhythmics", "vasopressor_6h",
      "transfusion_6h",
      "aki_kdigo1", "hospital_mortality"), vars)

    cat(sprintf("  %s AC: N=%d (Mg+K=%d, K-only=%d)\n",
        db_label, nrow(d_ac), sum(d_ac$trt_label == "Mg+K"),
        sum(d_ac$trt_label == "K-only")))

    t1 <- CreateTableOne(vars = vars, strata = "trt_label",
                          factorVars = cat_vars, data = d_ac, test = FALSE)
    out <- print(t1, smd = TRUE, printToggle = FALSE)
    cat("\n"); print(print(t1, smd = TRUE, printToggle = TRUE))
    df <- as.data.frame(out)
    df$Variable <- rownames(df)
    df$Database <- db_label
    df
  }

  t1_e <- make_ac_table(dat_e, "eICU")
  t1_m <- make_ac_table(dat_m, "MIMIC-IV")
  combined <- rbind(t1_e, t1_m)
  write.csv(combined, file.path(RESULTS, "07d_ac_table1.csv"),
            row.names = FALSE)
  cat("  Saved: 07d_ac_table1.csv\n")
}

cat(sprintf("\n%s\n07_sensitivity.R COMPLETE\n%s\n",
            strrep("=", 65), strrep("=", 65)))
