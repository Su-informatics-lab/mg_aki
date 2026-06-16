#!/usr/bin/env Rscript
# ============================================================================
# 08d_poaf.R — POAF as Downstream Exploratory Outcome
#
# Mg supplementation is guideline-recommended for POAF prevention.
# Question: does the POAF protective effect show the same Mg threshold
# pattern as AKI, or does it appear across all Mg levels?
#
# Run: Rscript 08d_poaf.R
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
    if (old %in% names(d) && !new %in% names(d)) names(d)[names(d)==old] <- new
  }
  if (is.character(d$age)) {
    d$age <- suppressWarnings(as.numeric(d$age))
    d$age[is.na(d$age)] <- 90
  }
  if ("surgery_type" %in% names(d)) {
    d$surg_cabg <- as.integer(d$surgery_type == "cabg")
    d$surg_valve <- as.integer(d$surgery_type == "valve")
    d$surg_combined <- as.integer(d$surgery_type == "combined")
  }
  if ("first_lactate" %in% names(d)) {
    d$lactate_missing <- as.integer(is.na(d$first_lactate))
    d$first_lactate[is.na(d$first_lactate)] <- median(d$first_lactate, na.rm = TRUE)
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
  list(or = exp(ct[tr, 1]), lo = exp(ct[tr, 1] - 1.96 * ct[tr, 2]),
       hi = exp(ct[tr, 1] + 1.96 * ct[tr, 2]),
       p = 2 * pnorm(-abs(ct[tr, 1] / ct[tr, 2])),
       logOR = ct[tr, 1], se = ct[tr, 2])
}

# ── PS covariates (WITH first_mg_value for overall; WITHOUT for stratified) ──
ps_covars_full <- c("age", "is_female", "bmi",
  "surg_cabg", "surg_valve", "surg_combined",
  "heart_failure", "hypertension", "diabetes", "ckd",
  "copd", "pvd", "stroke", "liver_disease",
  "baseline_creatinine", "egfr",
  "loop_diuretics", "nsaids", "acei_arb", "ppi",
  "beta_blockers", "steroids", "antiarrhythmics",
  "first_potassium", "first_calcium", "first_heartrate",
  "vasopressor_6h", "first_mg_value")

ps_covars_no_mg <- setdiff(ps_covars_full, "first_mg_value")

MG_CUTS   <- c(0, 1.8, 2.0, 2.3, Inf)
MG_LABELS <- c("<1.8", "1.8-2.0", "2.0-2.3", ">2.3")

# ============================================================================
run_poaf <- function(dat, db_name, has_cluster = FALSE) {

  cat(sprintf("\n%s\n%s\n%s\n", strrep("=", 65), db_name, strrep("=", 65)))

  # ── Identify POAF columns ──────────────────────────────────────
  poaf_cols <- intersect(c("poaf", "poaf_composite", "poaf_cardioversion"), names(dat))
  cat(sprintf("  POAF columns available: %s\n", paste(poaf_cols, collapse = ", ")))

  if (length(poaf_cols) == 0) {
    cat("  No POAF data — SKIP\n")
    return(NULL)
  }

  # ── Exclude pre-existing AF ────────────────────────────────────
  if ("preexisting_af" %in% names(dat)) {
    n_pre <- sum(dat$preexisting_af == 1, na.rm = TRUE)
    dat_poaf <- dat[dat$preexisting_af == 0, ]
    cat(sprintf("  Excluded pre-existing AF: %d → POAF-eligible: %d\n", n_pre, nrow(dat_poaf)))
  } else {
    dat_poaf <- dat
    cat("  No preexisting_af column — using full cohort\n")
  }

  ps_covars <- intersect(ps_covars_full, names(dat_poaf))
  if ("first_lactate" %in% names(dat_poaf))
    ps_covars <- c(ps_covars, "first_lactate", "lactate_missing")
  dat_poaf <- dat_poaf[complete.cases(dat_poaf[, ps_covars]), ]

  # ── Fit PS + OW on POAF-eligible cohort ────────────────────────
  ps_fml <- as.formula(paste("mg_supp ~", paste(ps_covars, collapse = " + ")))
  ps_fit <- glm(ps_fml, data = dat_poaf, family = binomial())
  dat_poaf$ps <- pmax(pmin(fitted(ps_fit), 0.99), 0.01)
  dat_poaf$ow <- ifelse(dat_poaf$mg_supp == 1, 1 - dat_poaf$ps, dat_poaf$ps)

  # Mg strata
  dat_poaf$mg_stratum <- cut(dat_poaf$first_mg_value,
    breaks = MG_CUTS, labels = MG_LABELS,
    right = FALSE, include.lowest = TRUE)

  results <- list()

  # ── Overall POAF rates ─────────────────────────────────────────
  cat(sprintf("\n  ── POAF Prevalence ──\n"))
  for (pc in poaf_cols) {
    n_event <- sum(dat_poaf[[pc]], na.rm = TRUE)
    n_elig <- sum(!is.na(dat_poaf[[pc]]))
    cat(sprintf("    %-25s %d / %d (%.1f%%)\n", pc, n_event, n_elig,
        100 * n_event / n_elig))
    # By treatment
    trt_rate <- mean(dat_poaf[[pc]][dat_poaf$mg_supp == 1], na.rm = TRUE)
    ctrl_rate <- mean(dat_poaf[[pc]][dat_poaf$mg_supp == 0], na.rm = TRUE)
    cat(sprintf("      Supplemented: %.1f%%  Not supplemented: %.1f%%\n",
        100 * trt_rate, 100 * ctrl_rate))
  }

  # ── Overall OW estimates ───────────────────────────────────────
  cat(sprintf("\n  ── Overall OW Treatment Effect ──\n"))
  cl <- if (has_cluster && "hospitalid" %in% names(dat_poaf)) dat_poaf$hospitalid else NULL

  for (pc in poaf_cols) {
    if (sum(dat_poaf[[pc]], na.rm = TRUE) < 10) {
      cat(sprintf("    %-25s SKIP (<%d events)\n", pc, 10))
      next
    }
    d_pc <- dat_poaf[!is.na(dat_poaf[[pc]]), ]
    cl_pc <- if (!is.null(cl)) cl[!is.na(dat_poaf[[pc]])] else NULL
    tryCatch({
      r <- wglm(as.formula(paste(pc, "~ mg_supp")), d_pc, d_pc$ow, cl_pc)
      sig <- ifelse(r$p < 0.05, " *", "")
      cat(sprintf("    %-25s OR %.3f (%.3f–%.3f) P=%.4f%s\n",
          pc, r$or, r$lo, r$hi, r$p, sig))
      results[[length(results) + 1]] <- data.frame(
        db = db_name, outcome = pc, stratum = "Overall", analysis = "all_ow",
        n = nrow(d_pc), n_trt = sum(d_pc$mg_supp), n_events = sum(d_pc[[pc]], na.rm = TRUE),
        or = round(r$or, 3), lo = round(r$lo, 3), hi = round(r$hi, 3),
        p = round(r$p, 4), logOR = round(r$logOR, 4), se = round(r$se, 4))
    }, error = function(e) cat(sprintf("    %-25s FAILED: %s\n", pc, e$message)))
  }

  # ── AC: Mg+K vs K-only ────────────────────────────────────────
  if ("ac_group" %in% names(dat_poaf)) {
    d_ac <- dat_poaf[dat_poaf$ac_group %in% c("mg_k", "k_only"), ]
    d_ac$ac_trt <- as.integer(d_ac$ac_group == "mg_k")
    if (sum(d_ac$ac_trt) >= 10 && sum(d_ac$ac_trt == 0) >= 10) {
      ac_fml <- as.formula(paste("ac_trt ~", paste(ps_covars, collapse = " + ")))
      tryCatch({
        ac_fit <- glm(ac_fml, data = d_ac, family = binomial())
        d_ac$ac_ps <- pmax(pmin(fitted(ac_fit), 0.99), 0.01)
        d_ac$ac_ow <- ifelse(d_ac$ac_trt == 1, 1 - d_ac$ac_ps, d_ac$ac_ps)

        cat(sprintf("\n  ── AC (Mg+K vs K-only, N=%d, trt=%d) ──\n",
            nrow(d_ac), sum(d_ac$ac_trt)))
        for (pc in poaf_cols) {
          if (sum(d_ac[[pc]], na.rm = TRUE) < 5) next
          d_pc <- d_ac[!is.na(d_ac[[pc]]), ]
          tryCatch({
            r <- wglm(as.formula(paste(pc, "~ ac_trt")), d_pc, d_pc$ac_ow, NULL)
            sig <- ifelse(r$p < 0.05, " *", "")
            cat(sprintf("    %-25s OR %.3f (%.3f–%.3f) P=%.4f%s\n",
                pc, r$or, r$lo, r$hi, r$p, sig))
            results[[length(results) + 1]] <- data.frame(
              db = db_name, outcome = pc, stratum = "Overall", analysis = "ac_ow",
              n = nrow(d_pc), n_trt = sum(d_pc$ac_trt),
              n_events = sum(d_pc[[pc]], na.rm = TRUE),
              or = round(r$or, 3), lo = round(r$lo, 3), hi = round(r$hi, 3),
              p = round(r$p, 4), logOR = round(r$logOR, 4), se = round(r$se, 4))
          }, error = function(e) cat(sprintf("    %-25s FAILED: %s\n", pc, e$message)))
        }
      }, error = function(e) cat(sprintf("  AC PS failed: %s\n", e$message)))
    }
  }

  # ── Mg-stratified POAF (primary definition only) ───────────────
  primary_poaf <- poaf_cols[1]  # "poaf" (dx-only)
  cat(sprintf("\n  ── Mg-Stratified POAF (%s) ──\n", primary_poaf))

  ps_covars_strat <- intersect(ps_covars_no_mg, names(dat_poaf))
  if ("first_lactate" %in% names(dat_poaf))
    ps_covars_strat <- c(ps_covars_strat, "first_lactate", "lactate_missing")

  for (s in MG_LABELS) {
    d_s <- dat_poaf[dat_poaf$mg_stratum == s, ]
    d_s <- d_s[!is.na(d_s[[primary_poaf]]), ]
    n_trt <- sum(d_s$mg_supp)
    n_ctrl <- nrow(d_s) - n_trt
    n_events <- sum(d_s[[primary_poaf]], na.rm = TRUE)

    if (n_trt < 10 || n_ctrl < 10 || n_events < 5) {
      cat(sprintf("    %-10s SKIP (trt=%d, ctrl=%d, events=%d)\n", s, n_trt, n_ctrl, n_events))
      next
    }

    d_s <- d_s[complete.cases(d_s[, ps_covars_strat]), ]
    tryCatch({
      ps_fml_s <- as.formula(paste("mg_supp ~", paste(ps_covars_strat, collapse = " + ")))
      ps_fit_s <- glm(ps_fml_s, data = d_s, family = binomial())
      d_s$ps_s <- pmax(pmin(fitted(ps_fit_s), 0.99), 0.01)
      d_s$ow_s <- ifelse(d_s$mg_supp == 1, 1 - d_s$ps_s, d_s$ps_s)

      cl_s <- if (has_cluster && "hospitalid" %in% names(d_s)) d_s$hospitalid else NULL
      r <- wglm(as.formula(paste(primary_poaf, "~ mg_supp")), d_s, d_s$ow_s, cl_s)

      trt_rate <- 100 * mean(d_s[[primary_poaf]][d_s$mg_supp == 1], na.rm = TRUE)
      ctrl_rate <- 100 * mean(d_s[[primary_poaf]][d_s$mg_supp == 0], na.rm = TRUE)
      sig <- ifelse(r$p < 0.05, " *", "")
      cat(sprintf("    %-10s N=%4d trt=%3d events=%3d  OR %.3f (%.3f–%.3f) P=%.4f  [trt %.1f%% ctrl %.1f%%]%s\n",
          s, nrow(d_s), n_trt, n_events, r$or, r$lo, r$hi, r$p, trt_rate, ctrl_rate, sig))

      results[[length(results) + 1]] <- data.frame(
        db = db_name, outcome = primary_poaf, stratum = s, analysis = "all_ow_stratified",
        n = nrow(d_s), n_trt = n_trt, n_events = n_events,
        or = round(r$or, 3), lo = round(r$lo, 3), hi = round(r$hi, 3),
        p = round(r$p, 4), logOR = round(r$logOR, 4), se = round(r$se, 4))
    }, error = function(e) cat(sprintf("    %-10s FAILED: %s\n", s, e$message)))
  }

  if (length(results) > 0) do.call(rbind, results) else NULL
}

# ============================================================================
# RUN
# ============================================================================
cat(strrep("=", 65), "\n")
cat("08d: POAF EXPLORATORY ANALYSIS\n")
cat(strrep("=", 65), "\n")

dat_e <- stdz(read.csv(file.path(RESULTS, "01_analysis_a_cohort.csv"), stringsAsFactors = FALSE))
dat_m <- stdz(read.csv(file.path(RESULTS, "04_mimic_cohort.csv"), stringsAsFactors = FALSE))

res_e <- run_poaf(dat_e, "eICU", has_cluster = TRUE)
res_m <- run_poaf(dat_m, "MIMIC", has_cluster = FALSE)

all_res <- rbind(res_e, res_m)
if (!is.null(all_res) && nrow(all_res) > 0) {
  outpath <- file.path(RESULTS, "08d_poaf.csv")
  write.csv(all_res, outpath, row.names = FALSE)
  cat(sprintf("\n✓ Saved: %s (%d rows)\n", outpath, nrow(all_res)))
} else {
  cat("\nNo POAF results to save.\n")
}

cat(sprintf("\n%s\nINTERPRETATION KEY\n%s\n", strrep("=", 65), strrep("=", 65)))
cat("  If POAF protection across ALL Mg strata → AF mechanism has lower threshold than AKI\n")
cat("  If POAF protection only in >2.3         → same threshold for both outcomes\n")
cat("  If POAF null everywhere                 → underpowered or our phenotype is noisy\n")
