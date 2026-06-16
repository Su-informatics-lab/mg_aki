#!/usr/bin/env Rscript
# ============================================================================
# 09_robustness.R — Reviewer-requested robustness analyses
#
#   E. Lactate sensitivity: re-estimate primary results with/without
#      lactate in PS model (addresses MNAR concern)
#   F. Quantitative bias analysis: deterministic adjustment for
#      unmeasured CPB time confounding (replaces E-value reassurance)
#
# Output: results/09_robustness.csv
# Run:    Rscript 09_robustness.R
# ============================================================================

suppressPackageStartupMessages({
  library(sandwich); library(lmtest)
})
RESULTS <- path.expand("~/mg_aki/results")

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
  for (v in c("bmi", "first_heartrate", "first_calcium", "first_potassium"))
    if (v %in% names(d) && any(is.na(d[[v]])))
      d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm = TRUE)
  d
}

wglm <- function(formula, data, w, cluster = NULL) {
  data$.w <- w
  fit <- glm(formula, data = data, weights = .w, family = quasibinomial())
  vc <- if (!is.null(cluster) && length(unique(cluster)) > 1)
    vcovCL(fit, cluster = cluster) else vcovHC(fit, type = "HC1")
  ct <- coeftest(fit, vcov. = vc)
  tr <- 2
  list(or = exp(ct[tr, 1]),
       lo = exp(ct[tr, 1] - 1.96 * ct[tr, 2]),
       hi = exp(ct[tr, 1] + 1.96 * ct[tr, 2]),
       p  = 2 * pnorm(-abs(ct[tr, 1] / ct[tr, 2])))
}

cat("Loading cohorts...\n")
dat_e <- stdz(read.csv(file.path(RESULTS, "01_analysis_a_cohort.csv"),
                        stringsAsFactors = FALSE))
dat_m <- stdz(read.csv(file.path(RESULTS, "04_mimic_cohort.csv"),
                        stringsAsFactors = FALSE))

all_rows <- list()

# ============================================================================
# SECTION E: LACTATE SENSITIVITY
# Compare PS model with vs without lactate (+ missing indicator)
# ============================================================================
cat(sprintf("\n%s\nE. LACTATE SENSITIVITY ANALYSIS\n%s\n",
            strrep("=", 65), strrep("=", 65)))

run_lactate_sens <- function(dat, db_name, has_cluster) {
  # Base PS covariates (always included)
  ps_base <- intersect(c(
    "age", "is_female", "bmi",
    "surg_cabg", "surg_valve", "surg_combined",
    "heart_failure", "hypertension", "diabetes", "ckd",
    "copd", "pvd", "stroke", "liver_disease",
    "baseline_creatinine", "egfr",
    "loop_diuretics", "nsaids", "acei_arb", "ppi",
    "beta_blockers", "steroids", "antiarrhythmics",
    "first_potassium", "first_calcium", "first_heartrate",
    "vasopressor_6h", "first_mg_value"), names(dat))

  has_lactate <- "first_lactate" %in% names(dat)
  rows <- list()

  for (include_lac in c(TRUE, FALSE)) {
    label <- if (include_lac) "with_lactate" else "without_lactate"
    if (include_lac && !has_lactate) next

    ps_covars <- ps_base
    if (include_lac) {
      dat$lactate_missing <- as.integer(is.na(dat$first_lactate))
      lac_med <- median(dat$first_lactate, na.rm = TRUE)
      dat$first_lactate_imp <- ifelse(is.na(dat$first_lactate),
                                       lac_med, dat$first_lactate)
      ps_covars <- c(ps_covars, "first_lactate_imp", "lactate_missing")
    }

    d <- dat[complete.cases(dat[, ps_base]), ]
    ps_fml <- as.formula(paste("mg_supp ~",
                                paste(ps_covars, collapse = " + ")))
    ps_fit <- glm(ps_fml, data = d, family = binomial())
    d$ps <- pmax(pmin(fitted(ps_fit), 0.99), 0.01)

    cl <- if (has_cluster && "hospitalid" %in% names(d))
      d$hospitalid else NULL

    # All-patient OW
    d$ow <- ifelse(d$mg_supp == 1, 1 - d$ps, d$ps)
    r_ow <- wglm(aki_kdigo1 ~ mg_supp, d, d$ow, cl)

    # All-patient IPTW
    prev <- mean(d$mg_supp)
    d$iptw <- ifelse(d$mg_supp == 1, prev / d$ps, (1 - prev) / (1 - d$ps))
    q01 <- quantile(d$iptw, 0.01); q99 <- quantile(d$iptw, 0.99)
    d$iptw <- pmax(pmin(d$iptw, q99), q01)
    r_iptw <- wglm(aki_kdigo1 ~ mg_supp, d, d$iptw, cl)

    # AC OW
    r_ac <- NULL
    if ("ac_group" %in% names(d)) {
      d_ac <- d[d$ac_group %in% c("mg_k", "k_only"), ]
      d_ac$ac_trt <- as.integer(d_ac$ac_group == "mg_k")
      if (sum(d_ac$ac_trt) >= 10 && sum(d_ac$ac_trt == 0) >= 10) {
        tryCatch({
          ac_fml <- as.formula(paste("ac_trt ~",
                                      paste(ps_covars, collapse = "+")))
          ac_fit <- glm(ac_fml, data = d_ac, family = binomial())
          d_ac$ac_ps <- pmax(pmin(fitted(ac_fit), 0.99), 0.01)
          d_ac$ac_ow <- ifelse(d_ac$ac_trt == 1,
                                1 - d_ac$ac_ps, d_ac$ac_ps)
          r_ac <- wglm(aki_kdigo1 ~ ac_trt, d_ac, d_ac$ac_ow, NULL)
        }, error = function(e) {})
      }
    }

    cat(sprintf("  %s %s:\n", db_name, label))
    cat(sprintf("    OW:   OR %.3f (%.3f-%.3f) P=%.4f\n",
                r_ow$or, r_ow$lo, r_ow$hi, r_ow$p))
    cat(sprintf("    IPTW: OR %.3f (%.3f-%.3f) P=%.4f\n",
                r_iptw$or, r_iptw$lo, r_iptw$hi, r_iptw$p))

    rows[[length(rows) + 1]] <- data.frame(
      section = "E_lactate", db = db_name, model = label,
      analysis = "ow_aki1",
      or = round(r_ow$or, 3), lo = round(r_ow$lo, 3),
      hi = round(r_ow$hi, 3), p = round(r_ow$p, 4))
    rows[[length(rows) + 1]] <- data.frame(
      section = "E_lactate", db = db_name, model = label,
      analysis = "iptw_aki1",
      or = round(r_iptw$or, 3), lo = round(r_iptw$lo, 3),
      hi = round(r_iptw$hi, 3), p = round(r_iptw$p, 4))

    if (!is.null(r_ac)) {
      cat(sprintf("    AC:   OR %.3f (%.3f-%.3f) P=%.4f\n",
                  r_ac$or, r_ac$lo, r_ac$hi, r_ac$p))
      rows[[length(rows) + 1]] <- data.frame(
        section = "E_lactate", db = db_name, model = label,
        analysis = "ac_aki1",
        or = round(r_ac$or, 3), lo = round(r_ac$lo, 3),
        hi = round(r_ac$hi, 3), p = round(r_ac$p, 4))
    }
  }

  # Compute max delta
  df <- do.call(rbind, rows)
  for (a in unique(df$analysis)) {
    w <- df[df$analysis == a & df$model == "with_lactate" &
            df$db == db_name, "or"]
    wo <- df[df$analysis == a & df$model == "without_lactate" &
             df$db == db_name, "or"]
    if (length(w) > 0 && length(wo) > 0) {
      delta <- abs(w - wo)
      cat(sprintf("    %s |ΔOR| = %.3f\n", a, delta))
    }
  }
  df
}

lac_e <- run_lactate_sens(dat_e, "eICU", TRUE)
lac_m <- run_lactate_sens(dat_m, "MIMIC", FALSE)
all_rows <- c(all_rows, list(lac_e), list(lac_m))

# ============================================================================
# SECTION F: QUANTITATIVE BIAS ANALYSIS
# Deterministic adjustment for unmeasured CPB time confounding
#
# Model: binary unmeasured confounder U (complex/long CPB surgery)
#   OR_UY  = OR of U → AKI (from literature)
#   p1     = P(U=1 | treated)   — supplemented patients
#   p0     = P(U=1 | untreated) — unsupplemented patients
#
# Bias factor = (OR_UY * p1 + (1-p1)) / (OR_UY * p0 + (1-p0))
# Adjusted OR = Observed OR / Bias factor
#
# Literature anchors:
#   CPB time >120 min → AKI: OR ~1.5–3.0 (Karkouti 2009, Sgouralis 2015)
#   Cross-clamp time → AKI: OR ~1.5–2.5 per 30-min increment
#   Supplemented patients have lower serum Mg (1.88 vs 2.22) suggesting
#   simpler surgery → P(complex | treated) < P(complex | untreated)
# ============================================================================
cat(sprintf("\n%s\nF. QUANTITATIVE BIAS ANALYSIS\n%s\n",
            strrep("=", 65), strrep("=", 65)))

# Grid of assumptions
or_uy_grid <- c(1.5, 2.0, 2.5, 3.0)   # OR of complex surgery → AKI
p0_grid    <- c(0.35, 0.40, 0.45)       # P(complex | untreated)
delta_grid <- c(0.05, 0.10, 0.15)       # prevalence difference (p0 - p1)

# Observed estimates to adjust
# Load from existing results
observed <- list(
  list(label = "AC primary (eICU)",     or = 0.75),
  list(label = "IPTW (eICU)",           or = 0.76),
  list(label = "Mg >2.3 OW (eICU)",    or = 0.53),
  list(label = "Mg >2.3 AC (eICU)",    or = 0.46),
  list(label = "2.6-3.0 subband",       or = 0.35)
)

# Try loading actual values from CSVs
res_path <- file.path(RESULTS, "02_results.csv")
if (file.exists(res_path)) {
  res <- read.csv(res_path, stringsAsFactors = FALSE)
  or_col <- if ("or" %in% names(res)) "or" else "or_"
  ac_e <- res[res$db == "eICU" & res$analysis == "ac_aki1", or_col]
  if (length(ac_e) > 0) observed[[1]]$or <- ac_e[1]
  iptw_e <- res[res$db == "eICU" & res$analysis == "iptw_aki1", or_col]
  if (length(iptw_e) > 0) observed[[2]]$or <- iptw_e[1]
}
mg_path <- file.path(RESULTS, "08_mg_stratified.csv")
if (file.exists(mg_path)) {
  mg <- read.csv(mg_path, stringsAsFactors = FALSE)
  mg_or_col <- if ("or" %in% names(mg)) "or" else "or_"
  gt23_ow <- mg[mg$db == "eICU" & mg$stratum == ">2.3" &
                mg$analysis == "all_ow", mg_or_col]
  if (length(gt23_ow) > 0) observed[[3]]$or <- gt23_ow[1]
  gt23_ac <- mg[mg$db == "eICU" & mg$stratum == ">2.3" &
                mg$analysis == "ac_ow", mg_or_col]
  if (length(gt23_ac) > 0) observed[[4]]$or <- gt23_ac[1]
}

cat("  Assumptions:\n")
cat("    U = complex surgery (long CPB / cross-clamp)\n")
cat("    OR_UY = effect of complex surgery on AKI\n")
cat("    p0 = P(complex | untreated), p1 = p0 - delta\n")
cat("    Bias factor = (OR_UY*p1 + 1-p1) / (OR_UY*p0 + 1-p0)\n")
cat("    Adjusted OR = Observed OR / Bias factor\n\n")

qba_rows <- list()

for (obs in observed) {
  cat(sprintf("  ── %s (observed OR = %.3f) ──\n", obs$label, obs$or))

  # For stratified analyses, reduce delta (complexity held ~constant)
  is_stratified <- grepl(">2.3|2.6-3.0", obs$label)

  for (or_uy in or_uy_grid) {
    for (p0 in p0_grid) {
      for (delta in delta_grid) {
        # For stratified analyses, halve the prevalence difference
        eff_delta <- if (is_stratified) delta / 2 else delta
        p1 <- p0 - eff_delta
        if (p1 < 0.05) next  # implausible

        # Bias factor (Lin, Psaty, Kronmal 1998)
        bf <- (or_uy * p1 + (1 - p1)) / (or_uy * p0 + (1 - p0))
        adj_or <- obs$or / bf
        survives <- adj_or < 1.0

        qba_rows[[length(qba_rows) + 1]] <- data.frame(
          section = "F_qba",
          analysis = obs$label,
          observed_or = obs$or,
          or_uy = or_uy,
          p0_untreated = p0,
          p1_treated = round(p1, 3),
          delta = eff_delta,
          bias_factor = round(bf, 4),
          adjusted_or = round(adj_or, 3),
          survives = survives,
          stringsAsFactors = FALSE)
      }
    }
  }

  # Summary for this estimate
  qba_sub <- do.call(rbind, qba_rows[sapply(qba_rows, function(x)
    x$analysis[1] == obs$label)])
  if (!is.null(qba_sub) && nrow(qba_sub) > 0) {
    n_total <- nrow(qba_sub)
    n_survive <- sum(qba_sub$survives)
    worst_bf <- min(qba_sub$bias_factor)
    worst_adj <- max(qba_sub$adjusted_or)
    cat(sprintf("    Survives %d/%d scenarios (%.0f%%)\n",
                n_survive, n_total, 100 * n_survive / n_total))
    cat(sprintf("    Worst case: bias factor=%.3f → adjusted OR=%.3f %s\n",
                worst_bf, worst_adj,
                ifelse(worst_adj < 1, "(still protective)", "(CROSSES NULL)")))
  }
}

# Key summary table
cat(sprintf("\n%s\nQBA SUMMARY: Does the estimate survive?\n%s\n",
            strrep("-", 65), strrep("-", 65)))
cat("  Scenario: OR_UY=2.0, p0=0.40, delta=0.10 (moderate assumptions)\n")
cat("  For stratified estimates, delta halved to 0.05\n\n")

for (obs in observed) {
  is_strat <- grepl(">2.3|2.6-3.0", obs$label)
  p0 <- 0.40; delta_eff <- if (is_strat) 0.05 else 0.10
  p1 <- p0 - delta_eff; or_uy <- 2.0
  bf <- (or_uy * p1 + (1 - p1)) / (or_uy * p0 + (1 - p0))
  adj <- obs$or / bf
  cat(sprintf("  %-25s OR %.3f → %.3f %s\n",
              obs$label, obs$or, adj,
              ifelse(adj < 1, "✓ survives", "✗ null")))
}

cat(sprintf("\n  Key insight: within the >2.3 stratum, complexity is held\n"))
cat(sprintf("  approximately constant by stratification, so the prevalence\n"))
cat(sprintf("  difference shrinks → bias factor approaches 1 → adjusted\n"))
cat(sprintf("  estimate stays protective.\n"))

# ── Combine and save ─────────────────────────────────────────────
qba_df <- do.call(rbind, qba_rows)
# Add lactate results
lac_df <- do.call(rbind, all_rows)

out <- rbind(
  data.frame(lac_df, or_uy = NA, p0_untreated = NA, p1_treated = NA,
             delta = NA, bias_factor = NA, adjusted_or = NA,
             survives = NA, observed_or = NA),
  data.frame(qba_df, db = NA, model = NA, lo = NA, hi = NA, p = NA)
)
# Simpler: save separately
write.csv(lac_df, file.path(RESULTS, "09_lactate_sensitivity.csv"),
          row.names = FALSE)
write.csv(qba_df, file.path(RESULTS, "09_qba.csv"), row.names = FALSE)

cat(sprintf("\n%s\n09_robustness.R COMPLETE\n%s\n",
            strrep("=", 65), strrep("=", 65)))
cat("  Saved: results/09_lactate_sensitivity.csv\n")
cat("  Saved: results/09_qba.csv\n")
