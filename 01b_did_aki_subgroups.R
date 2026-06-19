#!/usr/bin/env Rscript
# ============================================================================
# 01b_did_aki_subgroups.R — AKI binary endpoint sweep across subgroups
#
# Uses the primary matched dataset (did_matched_{db}_24h.csv)
# Defines KDIGO AKI from cr_pre / cr_post already in the matched set
# Sweeps subgroups: age, sex, eGFR, Cr, diabetes, CKD, HF, surgery, Mg
#
# Run:  Rscript 01b_did_aki_subgroups.R eicu
#       Rscript 01b_did_aki_subgroups.R mimic
# ============================================================================

suppressPackageStartupMessages({
  library(sandwich)
  library(lmtest)
})

RESULTS <- path.expand("~/mg_aki/results")

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

# ── Binary outcome helper ────────────────────────────────────────────────
run_binary <- function(df, outcome_col, ps_vars) {
  df$y <- as.numeric(df[[outcome_col]])
  df <- df[!is.na(df$y), ]
  nt <- sum(df$treated == 1); nc <- sum(df$treated == 0)
  if (nt < 15 || nc < 15) return(NULL)
  r1 <- mean(df$y[df$treated == 1]); r0 <- mean(df$y[df$treated == 0])
  events_t <- sum(df$y[df$treated == 1]); events_c <- sum(df$y[df$treated == 0])
  if (events_t + events_c < 5) return(NULL)

  # Adjusted GLM
  avail <- intersect(ps_vars, names(df))
  avail <- avail[sapply(avail, function(v) var(df[[v]], na.rm = T) > 1e-10)]
  smds <- sapply(avail, function(v) {
    x1 <- df[[v]][df$treated==1]; x0 <- df[[v]][df$treated==0]
    sp <- sqrt((var(x1,na.rm=T)+var(x0,na.rm=T))/2)
    if(is.na(sp)||sp<1e-10) 0 else abs(mean(x1,na.rm=T)-mean(x0,na.rm=T))/sp
  })
  adj <- names(smds[smds > 0.05])

  if (length(adj) > 0) {
    fml <- as.formula(paste("y ~ treated +", paste(adj, collapse = "+")))
  } else {
    fml <- y ~ treated
  }

  fit <- tryCatch(glm(fml, data = df, family = quasibinomial()), error = function(e) NULL)
  if (is.null(fit) || !"treated" %in% names(coef(fit))) return(NULL)

  vc <- tryCatch({
    if ("match_pair_id" %in% names(df) && length(unique(df$match_pair_id)) > 1)
      vcovCL(fit, cluster = df$match_pair_id)
    else vcovHC(fit, type = "HC1")
  }, error = function(e) tryCatch(vcovHC(fit, type = "HC1"), error = function(e2) vcov(fit)))

  ct <- tryCatch(coeftest(fit, vcov. = vc), error = function(e) NULL)
  if (is.null(ct)) return(NULL)
  trt_row <- which(rownames(ct) == "treated")
  if (length(trt_row) == 0) return(NULL)
  p_col <- grep("^Pr", colnames(ct))
  if (length(p_col) == 0) return(NULL)

  or <- exp(ct[trt_row, "Estimate"])
  or_lo <- exp(ct[trt_row, "Estimate"] - 1.96 * ct[trt_row, "Std. Error"])
  or_hi <- exp(ct[trt_row, "Estimate"] + 1.96 * ct[trt_row, "Std. Error"])
  p <- ct[trt_row, p_col]
  ard <- r1 - r0
  nnt <- if (ard < 0 && abs(ard) > 0.001) round(1 / abs(ard)) else NA

  list(n_trt = nt, n_ctl = nc, events_trt = events_t, events_ctl = events_c,
       rate_trt = r1, rate_ctl = r0, ard = ard,
       or = or, or_lo = or_lo, or_hi = or_hi, p = p, nnt = nnt)
}

# ============================================================================
run_sweep <- function(db) {
  tag <- tolower(db)
  SEP <- paste(rep("=", 70), collapse = "")
  cat(sprintf("\n%s\n%s: AKI Subgroup Sweep\n%s\n", SEP, db, SEP))

  path <- file.path(RESULTS, sprintf("did_matched_%s_24h.csv", tag))
  if (!file.exists(path)) { cat("  File not found:", path, "\n"); return(NULL) }

  df <- read.csv(path, stringsAsFactors = FALSE)
  cat(sprintf("  Loaded: %d rows (%d treated, %d control)\n",
              nrow(df), sum(df$treated == 1), sum(df$treated == 0)))

  # Median impute PS covariates
  for (v in PS_COVARS)
    if (v %in% names(df) && any(is.na(df[[v]])))
      df[[v]][is.na(df[[v]])] <- median(df[[v]], na.rm = TRUE)

  # ── Define AKI endpoints ─────────────────────────────────────────────
  df$aki_kdigo1 <- as.integer(df$delta_cr >= 0.3)
  df$aki_kdigo2 <- as.integer(!is.na(df$cr_pre) & df$cr_pre > 0 &
                               df$cr_post >= 2.0 * df$cr_pre)
  df$aki_kdigo3 <- as.integer((!is.na(df$cr_pre) & df$cr_pre > 0 &
                                df$cr_post >= 3.0 * df$cr_pre) |
                               df$cr_post >= 4.0)
  df$aki_any <- df$aki_kdigo1  # KDIGO ≥ Stage 1

  cat(sprintf("  AKI rates (overall):\n"))
  for (oc in c("aki_kdigo1", "aki_kdigo2", "aki_kdigo3")) {
    r1 <- mean(df[[oc]][df$treated == 1], na.rm = T)
    r0 <- mean(df[[oc]][df$treated == 0], na.rm = T)
    cat(sprintf("    %s: treated %.1f%%, control %.1f%%\n", oc, 100*r1, 100*r0))
  }

  # ── Define subgroups ─────────────────────────────────────────────────
  pair_ids <- df$match_pair_id[df$treated == 1]

  subgroups <- list(
    # Overall
    list(name = "Overall", var = NULL, val = NULL),

    # Age
    list(name = "Age < 65", var = "age", op = "<", val = 65),
    list(name = "Age >= 65", var = "age", op = ">=", val = 65),
    list(name = "Age >= 75", var = "age", op = ">=", val = 75),

    # Sex
    list(name = "Female", var = "is_female", op = "==", val = 1),
    list(name = "Male", var = "is_female", op = "==", val = 0),

    # eGFR
    list(name = "eGFR < 45", var = "egfr", op = "<", val = 45),
    list(name = "eGFR 45-60", var = "egfr", op = "range", val = c(45, 60)),
    list(name = "eGFR 60-90", var = "egfr", op = "range", val = c(60, 90)),
    list(name = "eGFR >= 90", var = "egfr", op = ">=", val = 90),
    list(name = "eGFR < 60", var = "egfr", op = "<", val = 60),

    # Baseline Cr
    list(name = "Cr_pre <= 1.0", var = "cr_pre", op = "<=", val = 1.0),
    list(name = "Cr_pre > 1.0", var = "cr_pre", op = ">", val = 1.0),
    list(name = "Cr_pre > 1.2", var = "cr_pre", op = ">", val = 1.2),
    list(name = "Cr_pre > 1.5", var = "cr_pre", op = ">", val = 1.5),

    # Comorbidities
    list(name = "Diabetes", var = "diabetes", op = "==", val = 1),
    list(name = "No diabetes", var = "diabetes", op = "==", val = 0),
    list(name = "CKD", var = "ckd", op = "==", val = 1),
    list(name = "No CKD", var = "ckd", op = "==", val = 0),
    list(name = "Heart failure", var = "heart_failure", op = "==", val = 1),
    list(name = "No HF", var = "heart_failure", op = "==", val = 0),
    list(name = "Hypertension", var = "hypertension", op = "==", val = 1),
    list(name = "No HTN", var = "hypertension", op = "==", val = 0),

    # Surgery
    list(name = "CABG", var = "surg_cabg", op = "==", val = 1),
    list(name = "Valve", var = "surg_valve", op = "==", val = 1),

    # Mg strata
    list(name = "Mg < 1.8", var = "first_mg_value", op = "<", val = 1.8),
    list(name = "Mg >= 2.0", var = "first_mg_value", op = ">=", val = 2.0),

    # BMI
    list(name = "BMI < 25", var = "bmi", op = "<", val = 25),
    list(name = "BMI 25-30", var = "bmi", op = "range", val = c(25, 30)),
    list(name = "BMI >= 30", var = "bmi", op = ">=", val = 30)
  )

  # ── Sweep ────────────────────────────────────────────────────────────
  aki_endpoints <- c("aki_kdigo1", "aki_kdigo2", "aki_kdigo3")
  aki_labels <- c(aki_kdigo1 = "KDIGO>=1 (dCr>=0.3)",
                   aki_kdigo2 = "KDIGO>=2 (Cr>=2x)",
                   aki_kdigo3 = "KDIGO>=3 (Cr>=3x|>=4)")

  all_results <- list(); ridx <- 0

  for (sg in subgroups) {
    # Subset by treated patient's value
    if (is.null(sg$var)) {
      # Overall
      sub <- df
    } else if (!sg$var %in% names(df)) {
      next
    } else {
      trt_vals <- df[[sg$var]][df$treated == 1]
      if (sg$op == "<") in_sg <- !is.na(trt_vals) & trt_vals < sg$val
      else if (sg$op == "<=") in_sg <- !is.na(trt_vals) & trt_vals <= sg$val
      else if (sg$op == ">") in_sg <- !is.na(trt_vals) & trt_vals > sg$val
      else if (sg$op == ">=") in_sg <- !is.na(trt_vals) & trt_vals >= sg$val
      else if (sg$op == "==") in_sg <- !is.na(trt_vals) & trt_vals == sg$val
      else if (sg$op == "range") in_sg <- !is.na(trt_vals) & trt_vals >= sg$val[1] & trt_vals < sg$val[2]
      else next

      pairs_in <- pair_ids[in_sg]
      sub <- df[df$match_pair_id %in% pairs_in, ]
    }

    for (aki in aki_endpoints) {
      res <- run_binary(sub, aki, PS_COVARS)
      ridx <- ridx + 1
      if (!is.null(res)) {
        all_results[[ridx]] <- data.frame(
          subgroup = sg$name, endpoint = aki,
          n_trt = res$n_trt, n_ctl = res$n_ctl,
          events_trt = res$events_trt, events_ctl = res$events_ctl,
          rate_trt = round(100 * res$rate_trt, 1),
          rate_ctl = round(100 * res$rate_ctl, 1),
          ard_pct = round(100 * res$ard, 1),
          or = round(res$or, 2), or_lo = round(res$or_lo, 2),
          or_hi = round(res$or_hi, 2),
          p = round(res$p, 4),
          nnt = ifelse(is.na(res$nnt), NA, res$nnt),
          stringsAsFactors = FALSE)
      } else {
        all_results[[ridx]] <- data.frame(
          subgroup = sg$name, endpoint = aki,
          n_trt = sum(sub$treated == 1), n_ctl = sum(sub$treated == 0),
          events_trt = NA, events_ctl = NA,
          rate_trt = NA, rate_ctl = NA, ard_pct = NA,
          or = NA, or_lo = NA, or_hi = NA, p = NA, nnt = NA,
          stringsAsFactors = FALSE)
      }
    }
  }

  results <- do.call(rbind, all_results)
  write.csv(results, file.path(RESULTS, sprintf("did_aki_subgroups_%s.csv", tag)), row.names = FALSE)

  # ── Print KDIGO ≥1 results (most clinically relevant) ────────────────
  cat(sprintf("\n%s\n%s: KDIGO ≥ Stage 1 (ΔCr ≥ 0.3 mg/dL)\n%s\n",
              SEP, db, SEP))
  cat("\n  subgroup              n_trt  n_ctl  AKI_trt  AKI_ctl   ARD    OR (95%CI)          P      NNT\n")
  cat("  ───────────────────  ─────  ─────  ───────  ───────  ─────  ──────────────────  ──────  ────\n")

  k1 <- results[results$endpoint == "aki_kdigo1", ]
  for (i in seq_len(nrow(k1))) {
    r <- k1[i, ]
    if (is.na(r$or)) {
      cat(sprintf("  %-21s  %5d  %5d     —        —       —         —               —     —\n",
                  r$subgroup, r$n_trt, r$n_ctl))
    } else {
      sig <- if (r$p < 0.05) " *" else ""
      nnt_str <- if (!is.na(r$nnt)) sprintf("%4d", r$nnt) else "   —"
      cat(sprintf("  %-21s  %5d  %5d   %5.1f%%   %5.1f%%  %+.1f%%  %5.2f (%5.2f–%5.2f)  %.4f%s %s\n",
                  r$subgroup, r$n_trt, r$n_ctl,
                  r$rate_trt, r$rate_ctl, r$ard_pct,
                  r$or, r$or_lo, r$or_hi, r$p, sig, nnt_str))
    }
  }

  # ── Print significant KDIGO ≥2 and ≥3 ────────────────────────────────
  for (ep in c("aki_kdigo2", "aki_kdigo3")) {
    sub <- results[results$endpoint == ep & !is.na(results$p) & results$p < 0.1, ]
    if (nrow(sub) > 0) {
      cat(sprintf("\n  %s — notable results (P < 0.1):\n", aki_labels[ep]))
      for (i in seq_len(nrow(sub))) {
        r <- sub[i, ]
        sig <- if (r$p < 0.05) " *" else ""
        cat(sprintf("    %s: %5.1f%% vs %5.1f%%, OR=%.2f (%.2f–%.2f), P=%.4f%s\n",
                    r$subgroup, r$rate_trt, r$rate_ctl,
                    r$or, r$or_lo, r$or_hi, r$p, sig))
      }
    }
  }

  cat(sprintf("\n  Saved: did_aki_subgroups_%s.csv (%d rows)\n", tag, nrow(results)))
  return(results)
}

# ============================================================================
cat("======================================================================\n")
cat("01b_did_aki_subgroups.R — AKI binary endpoint subgroup sweep\n")
cat("  Endpoints: KDIGO ≥1, ≥2, ≥3\n")
cat("  Subgroups: age, sex, eGFR, Cr, diabetes, CKD, HF, surgery, Mg, BMI\n")
cat("======================================================================\n")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) { cat("Usage: Rscript 01b_did_aki_subgroups.R eicu|mimic\n"); quit(status = 1) }
for (a in args) run_sweep(toupper(a))
