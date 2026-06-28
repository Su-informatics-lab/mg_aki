#!/usr/bin/env Rscript
# ============================================================================
# 03_hte.R — Heterogeneous Treatment Effects (PAIR-PRESERVING)
#
# CRITICAL FIX: Subgroup analyses subset by the TREATED patient's covariate
# and compare against their MATCHED CONTROL, regardless of the control's
# covariate value. This preserves the propensity-score-matched pair structure.
#
# WRONG:  df[df$ckd == 1, ]  → breaks pairs (filters both patients by own CKD)
# RIGHT:  idx <- which(ckd_trt == 1); run_or(aki_trt[idx], aki_ctl[idx])
#
# Outputs:
#   results/did_hte_{db}.csv         — single subgroup ORs
#   results/did_hte_crossed_{db}.csv — pairwise crossed phenotype ORs
#   results/did_hte_interact_{db}.csv — interaction test P-values
#
# Usage: Rscript 03_hte.R mimic
#        Rscript 03_hte.R eicu
# ============================================================================

suppressPackageStartupMessages({ library(sandwich); library(lmtest) })

RESULTS <- path.expand("~/mg_aki/results")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) { cat("Usage: Rscript 03_hte.R <db>\n"); quit(status = 1) }
tag <- tolower(args[1]); db <- toupper(tag)

SEP <- paste(rep("=", 70), collapse = "")
cat(sprintf("\n%s\n03_hte.R — %s (pair-preserving)\n%s\n", SEP, db, SEP))

# ══════════════════════════════════════════════════════════════════
# LOAD DATA
# ══════════════════════════════════════════════════════════════════
all_pts <- read.csv(file.path(RESULTS, sprintf("did_all_%s.csv", tag)), stringsAsFactors = FALSE)
pairs   <- read.csv(file.path(RESULTS, sprintf("did_pairs_primary_yet_untreated_%s.csv", tag)),
                     stringsAsFactors = FALSE)
cr_all  <- read.csv(file.path(RESULTS, sprintf("did_cr_all_%s.csv", tag)), stringsAsFactors = FALSE)

cr_id <- if ("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
cr_all$pid <- cr_all[[cr_id]]
if (!"offset_h" %in% names(cr_all)) cr_all$offset_h <- cr_all$labresultoffset / 60
cr_all <- cr_all[order(cr_all$pid, cr_all$offset_h), ]
cr_list <- split(cr_all[, c("labresult", "offset_h")], cr_all$pid)

trt_rows <- match(pairs$trt_pid, all_pts$pid)
ctl_rows <- match(pairs$ctl_pid, all_pts$pid)
n_pairs  <- nrow(pairs)

cat(sprintf("  Pairs: %d | Patients: %d\n", n_pairs, nrow(all_pts)))

# ══════════════════════════════════════════════════════════════════
# COMPUTE ALL OUTCOMES
# ══════════════════════════════════════════════════════════════════

# ── AKI stages (from Cr) ──────────────────────────────────────────
compute_aki <- function(pid, t_mg) {
  cr <- cr_list[[as.character(pid)]]
  if (is.null(cr) || nrow(cr) < 1) return(c(NA, NA, NA, NA))
  pre <- cr[cr$offset_h >= 0 & cr$offset_h < t_mg, ]
  if (nrow(pre) == 0) return(c(NA, NA, NA, NA))
  bl <- pre$labresult[which.max(pre$offset_h)]
  if (is.na(bl) || bl <= 0) return(c(NA, NA, NA, NA))

  # 48h AKI, 7d AKI (Stage 1+), Stage 2+, Stage 3+
  aki_48h <- 0; aki_7d <- 0; stage2 <- 0; stage3 <- 0
  post <- cr[cr$offset_h >= t_mg & cr$offset_h <= (t_mg + 168), ]
  if (nrow(post) == 0) return(c(0, 0, 0, 0))

  for (i in seq_len(nrow(post))) {
    h <- post$offset_h[i] - t_mg
    val <- post$labresult[i]
    delta <- val - bl
    ratio <- val / bl

    # KDIGO Stage 1+ consolidated
    if (h <= 48 && (delta >= 0.3 || ratio >= 1.5)) { aki_7d <- 1; aki_48h <- 1 }
    if (h > 48 && ratio >= 1.5) aki_7d <- 1

    # 48h AKI only (absolute + ratio, within 48h)
    if (h <= 48 && (delta >= 0.3 || ratio >= 1.5)) aki_48h <- 1

    # Stage 2+
    if (ratio >= 2.0) stage2 <- 1
    # Stage 3+
    if (ratio >= 3.0 || val >= 4.0) stage3 <- 1
  }
  c(aki_48h, aki_7d, stage2, stage3)
}

cat("  Computing AKI outcomes...\n")
aki_trt <- aki_ctl <- matrix(NA, n_pairs, 4)
colnames(aki_trt) <- colnames(aki_ctl) <- c("aki_48h", "aki_7d", "aki_stage2", "aki_stage3")
for (i in seq_len(n_pairs)) {
  aki_trt[i, ] <- compute_aki(pairs$trt_pid[i], pairs$t_mg[i])
  aki_ctl[i, ] <- compute_aki(pairs$ctl_pid[i], pairs$t_mg[i])
}

# ── Pre-computed binary outcomes (from did_all) ───────────────────
outcome_names <- c("aki_48h", "aki_7d", "aki_stage2", "aki_stage3",
                    "hosp_mortality", "poaf", "encephalopathy_delirium",
                    "transfusion", "reintubation", "vent_arrhythmia",
                    "poaf_icd", "encephalopathy_icd")

# Build outcome matrices: trt and ctl vectors, indexed by pair
out_trt <- list(); out_ctl <- list()
for (j in 1:4) {
  out_trt[[colnames(aki_trt)[j]]] <- aki_trt[, j]
  out_ctl[[colnames(aki_ctl)[j]]] <- aki_ctl[, j]
}
for (oc in c("hosp_mortality", "poaf", "encephalopathy_delirium",
             "transfusion", "reintubation", "vent_arrhythmia",
             "poaf_icd", "encephalopathy_icd")) {
  if (oc %in% names(all_pts)) {
    out_trt[[oc]] <- all_pts[[oc]][trt_rows]
    out_ctl[[oc]] <- all_pts[[oc]][ctl_rows]
  }
}
cat(sprintf("  Outcomes: %s\n", paste(names(out_trt), collapse = ", ")))

# ══════════════════════════════════════════════════════════════════
# TREATED PATIENT COVARIATES (for subgroup definition)
# ══════════════════════════════════════════════════════════════════

# All covariates indexed by pair (from TREATED patient only)
egfr_trt   <- all_pts$egfr[trt_rows]
age_trt    <- all_pts$age[trt_rows]
female_trt <- all_pts$is_female[trt_rows]
dm_trt     <- all_pts$diabetes[trt_rows]
ckd_trt    <- all_pts$ckd[trt_rows]
hf_trt     <- all_pts$heart_failure[trt_rows]
cabg_trt   <- if ("surg_cabg" %in% names(all_pts)) all_pts$surg_cabg[trt_rows] else rep(0, n_pairs)
bmi_trt    <- all_pts$bmi[trt_rows]
htn_trt    <- all_pts$hypertension[trt_rows]
copd_trt   <- all_pts$copd[trt_rows]

# Baseline Mg (try multiple column names)
mg_trt <- rep(NA, n_pairs)
for (mc in c("last_magnesium", "first_mg_value", "first_magnesium")) {
  if (mc %in% names(all_pts)) { mg_trt <- all_pts[[mc]][trt_rows]; break }
}

# ── Define subgroup flags ─────────────────────────────────────────
subgroups <- list(
  # Overall
  list(name = "Overall",           idx = seq_len(n_pairs)),
  # Age
  list(name = "Age < 65",          idx = which(!is.na(age_trt) & age_trt < 65)),
  list(name = "Age >= 65",         idx = which(!is.na(age_trt) & age_trt >= 65)),
  # eGFR (5-bin stratification — THE FIX)
  list(name = "eGFR >= 90",        idx = which(!is.na(egfr_trt) & egfr_trt >= 90)),
  list(name = "eGFR 60-89",        idx = which(!is.na(egfr_trt) & egfr_trt >= 60 & egfr_trt < 90)),
  list(name = "eGFR 45-59",        idx = which(!is.na(egfr_trt) & egfr_trt >= 45 & egfr_trt < 60)),
  list(name = "eGFR 30-44",        idx = which(!is.na(egfr_trt) & egfr_trt >= 30 & egfr_trt < 45)),
  list(name = "eGFR < 30",         idx = which(!is.na(egfr_trt) & egfr_trt < 30)),
  # Also keep the old binary for backward compat, but now correct
  list(name = "eGFR < 60",         idx = which(!is.na(egfr_trt) & egfr_trt < 60)),
  list(name = "eGFR >= 60",        idx = which(!is.na(egfr_trt) & egfr_trt >= 60)),
  # Mg
  list(name = "Mg < 1.8",          idx = which(!is.na(mg_trt) & mg_trt < 1.8)),
  list(name = "Mg >= 1.8",         idx = which(!is.na(mg_trt) & mg_trt >= 1.8)),
  # Surgery
  list(name = "CABG",              idx = which(cabg_trt == 1)),
  list(name = "Non-CABG",          idx = which(cabg_trt == 0)),
  # Comorbidities
  list(name = "Diabetes",          idx = which(dm_trt == 1)),
  list(name = "No diabetes",       idx = which(dm_trt == 0)),
  list(name = "CKD",               idx = which(ckd_trt == 1)),
  list(name = "No CKD",            idx = which(ckd_trt == 0)),
  list(name = "Heart failure",     idx = which(hf_trt == 1)),
  list(name = "No HF",             idx = which(hf_trt == 0)),
  list(name = "BMI >= 30",         idx = which(!is.na(bmi_trt) & bmi_trt >= 30)),
  list(name = "BMI < 30",          idx = which(!is.na(bmi_trt) & bmi_trt < 30)),
  list(name = "Female",            idx = which(female_trt == 1)),
  list(name = "Male",              idx = which(female_trt == 0)),
  # Crossed phenotypes
  list(name = "DM + CKD",          idx = which(dm_trt == 1 & ckd_trt == 1)),
  list(name = "HF + CABG",         idx = which(hf_trt == 1 & cabg_trt == 1)),
  list(name = "Mg<1.8 + CKD",      idx = which(!is.na(mg_trt) & mg_trt < 1.8 & ckd_trt == 1))
)

# ══════════════════════════════════════════════════════════════════
# OR HELPER (pair-preserving)
# ══════════════════════════════════════════════════════════════════
run_or <- function(ot, oc) {
  valid <- !is.na(ot) & !is.na(oc)
  ot <- ot[valid]; oc <- oc[valid]
  n <- sum(valid); et <- sum(ot); ec <- sum(oc)
  if (n < 30 || (et + ec) == 0)
    return(data.frame(or = NA, or_lo = NA, or_hi = NA, p = NA,
                      rate_trt = NA, rate_ctl = NA, n = n, events_trt = et, events_ctl = ec))
  r1 <- mean(ot); r0 <- mean(oc)
  df <- data.frame(outcome = c(ot, oc), treated = rep(c(1, 0), each = sum(valid)))
  fit <- tryCatch(glm(outcome ~ treated, data = df, family = quasibinomial()),
                  error = function(e) NULL)
  if (is.null(fit))
    return(data.frame(or = NA, or_lo = NA, or_hi = NA, p = NA,
                      rate_trt = r1, rate_ctl = r0, n = n, events_trt = et, events_ctl = ec))
  ct <- tryCatch(coeftest(fit, vcov. = vcovHC(fit, type = "HC1")),
                 error = function(e) tryCatch(coeftest(fit), error = function(e2) NULL))
  if (is.null(ct))
    return(data.frame(or = NA, or_lo = NA, or_hi = NA, p = NA,
                      rate_trt = r1, rate_ctl = r0, n = n, events_trt = et, events_ctl = ec))
  or <- exp(ct["treated", "Estimate"])
  lo <- exp(ct["treated", "Estimate"] - 1.96 * ct["treated", "Std. Error"])
  hi <- exp(ct["treated", "Estimate"] + 1.96 * ct["treated", "Std. Error"])
  p  <- ct["treated", ncol(ct)]
  data.frame(or = round(or, 4), or_lo = round(lo, 4), or_hi = round(hi, 4), p = round(p, 6),
             rate_trt = round(r1, 4), rate_ctl = round(r0, 4), n = n,
             events_trt = et, events_ctl = ec)
}

# ══════════════════════════════════════════════════════════════════
# SECTION 1: SINGLE SUBGROUP ANALYSIS
# ══════════════════════════════════════════════════════════════════
cat("\n── Section 1: Single subgroup ORs ──\n")

hte_results <- list(); ridx <- 0

for (sg in subgroups) {
  idx <- sg$idx
  if (length(idx) < 30) next
  for (oc in names(out_trt)) {
    res <- run_or(out_trt[[oc]][idx], out_ctl[[oc]][idx])
    res$subgroup <- sg$name
    res$outcome <- oc
    res$db <- db
    ridx <- ridx + 1
    hte_results[[ridx]] <- res
  }
}

hte_df <- do.call(rbind, hte_results)
hte_df$sig <- !is.na(hte_df$p) & hte_df$p < 0.05

# Print summary for key outcomes
for (oc in c("aki_7d", "hosp_mortality")) {
  cat(sprintf("\n  [%s]\n", oc))
  sub <- hte_df[hte_df$outcome == oc, ]
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    if (is.na(r$or)) { cat(sprintf("    %-20s  n=%d (skip)\n", r$subgroup, r$n)); next }
    sig <- if (r$sig) " *" else "  "
    cat(sprintf("    %-20s  OR=%.3f [%.3f,%.3f]  P=%.4f%s  %.1f%% vs %.1f%%  n=%d\n",
                r$subgroup, r$or, r$or_lo, r$or_hi, r$p, sig,
                100 * r$rate_trt, 100 * r$rate_ctl, r$n))
  }
}

# Save
outpath <- file.path(RESULTS, sprintf("did_hte_%s.csv", tag))
write.csv(hte_df, outpath, row.names = FALSE)
cat(sprintf("\n  Saved: %s (%d rows)\n", outpath, nrow(hte_df)))

# ══════════════════════════════════════════════════════════════════
# SECTION 2: PAIRWISE CROSSED PHENOTYPE SWEEP
# ══════════════════════════════════════════════════════════════════
cat("\n── Section 2: Crossed phenotype sweep (n_trt >= 30) ──\n")

subvar_list <- list(
  age_ge65     = which(!is.na(age_trt) & age_trt >= 65),
  is_female    = which(female_trt == 1),
  egfr_lt60    = which(!is.na(egfr_trt) & egfr_trt < 60),
  mg_lt18      = which(!is.na(mg_trt) & mg_trt < 1.8),
  surg_cabg    = which(cabg_trt == 1),
  diabetes     = which(dm_trt == 1),
  ckd          = which(ckd_trt == 1),
  heart_failure = which(hf_trt == 1),
  bmi_ge30     = which(!is.na(bmi_trt) & bmi_trt >= 30)
)

combos <- combn(names(subvar_list), 2, simplify = FALSE)
crossed_results <- list(); cidx <- 0

for (cb in combos) {
  idx <- intersect(subvar_list[[cb[1]]], subvar_list[[cb[2]]])
  if (length(idx) < 30) next
  label <- paste(cb, collapse = " + ")
  # Primary outcome: 7d AKI
  res <- run_or(out_trt[["aki_7d"]][idx], out_ctl[["aki_7d"]][idx])
  res$phenotype <- label
  res$outcome <- "aki_7d"
  res$db <- db
  cidx <- cidx + 1
  crossed_results[[cidx]] <- res
}

if (cidx > 0) {
  crossed_df <- do.call(rbind, crossed_results)
  crossed_df <- crossed_df[order(crossed_df$or), ]
  crossed_df$bonferroni_sig <- !is.na(crossed_df$p) & crossed_df$p < (0.05 / 36)

  cat(sprintf("\n  %-30s  %6s  %12s  %8s\n", "Phenotype", "OR", "95% CI", "P"))
  cat(paste(rep("-", 70), collapse = ""), "\n")
  for (i in seq_len(nrow(crossed_df))) {
    r <- crossed_df[i, ]
    if (is.na(r$or)) next
    sig <- if (!is.na(r$p) && r$p < 0.05) " *" else "  "
    bonf <- if (!is.na(r$bonferroni_sig) && r$bonferroni_sig) " **" else ""
    cat(sprintf("  %-30s  %6.3f  [%.3f,%.3f]  %8.4f%s%s  n=%d\n",
                r$phenotype, r$or, r$or_lo, r$or_hi, r$p, sig, bonf, r$n))
  }

  outpath2 <- file.path(RESULTS, sprintf("did_hte_crossed_%s.csv", tag))
  write.csv(crossed_df, outpath2, row.names = FALSE)
  cat(sprintf("\n  Saved: %s (%d rows)\n", outpath2, nrow(crossed_df)))
}

# ══════════════════════════════════════════════════════════════════
# SECTION 3: INTERACTION TESTS
# ══════════════════════════════════════════════════════════════════
cat("\n── Section 3: Interaction tests ──\n")

# Build stacked data for interaction model (pair-preserving)
# Each pair contributes two rows: treated + their matched control
# The subgroup variable is the TREATED patient's value for BOTH rows
# (because we're asking: does the treatment effect differ BY THIS CHARACTERISTIC?)

interact_vars <- list(
  "Age >= 65"  = as.integer(!is.na(age_trt) & age_trt >= 65),
  "eGFR < 60"  = as.integer(!is.na(egfr_trt) & egfr_trt < 60),
  "Mg < 1.8"   = as.integer(!is.na(mg_trt) & mg_trt < 1.8),
  "CABG"        = as.integer(cabg_trt == 1),
  "Diabetes"    = as.integer(dm_trt == 1),
  "CKD"         = as.integer(ckd_trt == 1),
  "Heart failure" = as.integer(hf_trt == 1),
  "BMI >= 30"   = as.integer(!is.na(bmi_trt) & bmi_trt >= 30)
)

interact_results <- list(); iidx <- 0

for (iv_name in names(interact_vars)) {
  sg_val <- interact_vars[[iv_name]]  # length = n_pairs, from TREATED patient

  for (oc in c("aki_7d", "aki_48h", "hosp_mortality")) {
    if (!(oc %in% names(out_trt))) next
    ot <- out_trt[[oc]]; oc_val <- out_ctl[[oc]]

    # Stack: treated rows then control rows
    valid <- !is.na(ot) & !is.na(oc_val) & !is.na(sg_val)
    n_valid <- sum(valid)
    if (n_valid < 50) next

    idf <- data.frame(
      outcome  = c(ot[valid], oc_val[valid]),
      treated  = rep(c(1, 0), each = n_valid),
      sg       = rep(sg_val[valid], 2)  # TREATED patient's value for BOTH rows
    )

    fit <- tryCatch(
      glm(outcome ~ treated * sg, data = idf, family = quasibinomial()),
      error = function(e) NULL)
    if (is.null(fit)) next

    ct <- tryCatch(coeftest(fit, vcov. = vcovHC(fit, type = "HC1")),
                   error = function(e) tryCatch(coeftest(fit), error = function(e2) NULL))
    if (is.null(ct) || !("treated:sg" %in% rownames(ct))) next

    p_interact <- ct["treated:sg", ncol(ct)]
    iidx <- iidx + 1
    interact_results[[iidx]] <- data.frame(
      variable = iv_name, outcome = oc,
      p_interaction = round(p_interact, 4),
      sig = p_interact < 0.05,
      db = db)
  }
}

if (iidx > 0) {
  interact_df <- do.call(rbind, interact_results)

  cat(sprintf("\n  %-18s  %-15s  %10s\n", "Variable", "Outcome", "P_interact"))
  cat(paste(rep("-", 50), collapse = ""), "\n")
  for (i in seq_len(nrow(interact_df))) {
    r <- interact_df[i, ]
    sig <- if (r$sig) " *" else "  "
    cat(sprintf("  %-18s  %-15s  %10.4f%s\n", r$variable, r$outcome, r$p_interaction, sig))
  }

  outpath3 <- file.path(RESULTS, sprintf("did_hte_interact_%s.csv", tag))
  write.csv(interact_df, outpath3, row.names = FALSE)
  cat(sprintf("\n  Saved: %s (%d rows)\n", outpath3, nrow(interact_df)))
}

# ══════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════
cat(sprintf("\n%s\n03_hte.R — %s DONE\n", SEP, db))
cat(sprintf("  did_hte_%s.csv:         %d rows (single subgroups × outcomes)\n", tag, nrow(hte_df)))
if (cidx > 0) cat(sprintf("  did_hte_crossed_%s.csv:  %d rows (crossed phenotypes)\n", tag, nrow(crossed_df)))
if (iidx > 0) cat(sprintf("  did_hte_interact_%s.csv: %d rows (interaction tests)\n", tag, nrow(interact_df)))
cat(sprintf("%s\n", SEP))
