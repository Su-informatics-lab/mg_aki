#!/usr/bin/env Rscript
# ============================================================================
# probe_complexity_nc.R — Complexity-specific negative control probe
#
# Tests: if Mg supplementation is NOT associated with outcomes driven
# by surgical complexity (but not by Mg), then complexity confounding
# is unlikely to explain the AKI finding.
#
# Outcomes (from eICU treatment/intakeOutput/diagnosis tables):
#   1. Reoperation / surgical re-exploration
#   2. Perioperative blood transfusion (PRBC/FFP/platelets)
#
# These are caused by complex/prolonged surgery, NOT by postoperative
# Mg supplementation. A null here = strong evidence against the
# "supplemented patients just had simpler surgery" critique.
#
# Run: Rscript probe_complexity_nc.R
# ============================================================================

suppressPackageStartupMessages({
  library(sandwich); library(lmtest)
})

RESULTS <- path.expand("~/mg_aki/results")
EICU    <- path.expand("~/mg_aki/eicu-crd-2.0")

# ── Find file helper ─────────────────────────────────────────────────
find_csv <- function(name) {
  for (n in c(name, tolower(name)))
    for (ext in c(".csv.gz", ".csv")) {
      p <- file.path(EICU, paste0(n, ext))
      if (file.exists(p)) return(p)
    }
  stop(paste("Not found:", name))
}

matches_any <- function(x, patterns) {
  grepl(paste(patterns, collapse = "|"), x, ignore.case = TRUE)
}

# ── Load cohort ──────────────────────────────────────────────────────
cat("Loading eICU cohort...\n")
dat <- read.csv(file.path(RESULTS, "01_analysis_a_cohort.csv"),
                stringsAsFactors = FALSE)

# Standardize
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
  if (old %in% names(dat) && !new %in% names(dat))
    names(dat)[names(dat) == old] <- new
}
if (is.character(dat$age)) {
  dat$age <- suppressWarnings(as.numeric(dat$age))
  dat$age[is.na(dat$age)] <- 90
}
if ("surgery_type" %in% names(dat)) {
  dat$surg_cabg     <- as.integer(dat$surgery_type == "cabg")
  dat$surg_valve    <- as.integer(dat$surgery_type == "valve")
  dat$surg_combined <- as.integer(dat$surgery_type == "combined")
}
if ("first_lactate" %in% names(dat)) {
  dat$lactate_missing <- as.integer(is.na(dat$first_lactate))
  dat$first_lactate[is.na(dat$first_lactate)] <- median(dat$first_lactate, na.rm = TRUE)
}
for (v in c("bmi", "first_heartrate", "first_calcium", "first_potassium"))
  if (v %in% names(dat) && any(is.na(dat[[v]])))
    dat[[v]][is.na(dat[[v]])] <- median(dat[[v]], na.rm = TRUE)

pids <- dat$patientunitstayid
cat(sprintf("  N = %d\n", nrow(dat)))

# ── Extract negative control outcomes from raw tables ────────────────
cat("\nExtracting complexity-driven outcomes...\n")

# 1. Reoperation / re-exploration (from treatment table)
cat("  1. Reoperation...\n")
tx <- read.csv(find_csv("treatment"), stringsAsFactors = FALSE)
names(tx) <- tolower(names(tx))
tx_elig <- tx[tx$patientunitstayid %in% pids, ]

reop_patterns <- c("reoperation", "re-exploration", "reexploration",
                    "return to or", "return to operating",
                    "mediastinal exploration", "surgical exploration",
                    "take back", "re-sternotomy", "resternotomy")
reop_pids <- unique(tx_elig$patientunitstayid[
  matches_any(tx_elig$treatmentstring, reop_patterns) &
  tx_elig$treatmentoffset > 0])
dat$nc_reoperation <- as.integer(dat$patientunitstayid %in% reop_pids)
cat(sprintf("    Events: %d (%.1f%%)\n", sum(dat$nc_reoperation),
            100 * mean(dat$nc_reoperation)))

# 2. Blood transfusion (from treatment + intakeOutput)
cat("  2. Transfusion...\n")
transfusion_patterns <- c("packed red blood", "prbc", "red blood cell",
                           "fresh frozen plasma", "ffp",
                           "platelet", "cryoprecipitate",
                           "blood product", "transfus")

# Treatment table
tx_transfusion_pids <- unique(tx_elig$patientunitstayid[
  matches_any(tx_elig$treatmentstring, transfusion_patterns) &
  tx_elig$treatmentoffset > 0])

# intakeOutput table (blood products administered)
io_transfusion_pids <- c()
tryCatch({
  io <- read.csv(find_csv("intakeOutput"), stringsAsFactors = FALSE)
  names(io) <- tolower(names(io))
  io_elig <- io[io$patientunitstayid %in% pids, ]
  io_transfusion_pids <- unique(io_elig$patientunitstayid[
    matches_any(io_elig$celllabel, transfusion_patterns) &
    io_elig$intakeoutputoffset > 0])
  cat(sprintf("    intakeOutput source: %d patients\n",
              length(io_transfusion_pids)))
}, error = function(e) cat(sprintf("    intakeOutput: %s\n", e$message)))

all_transfusion <- unique(c(tx_transfusion_pids, io_transfusion_pids))
dat$nc_transfusion <- as.integer(dat$patientunitstayid %in% all_transfusion)
cat(sprintf("    Events: %d (%.1f%%)\n", sum(dat$nc_transfusion),
            100 * mean(dat$nc_transfusion)))

# ── Estimate PS + OW (same model as primary analysis) ────────────────
cat("\nFitting PS + OW...\n")
ps_covars <- intersect(c("age", "is_female", "bmi",
  "surg_cabg", "surg_valve", "surg_combined",
  "heart_failure", "hypertension", "diabetes", "ckd",
  "copd", "pvd", "stroke", "liver_disease",
  "baseline_creatinine", "egfr",
  "loop_diuretics", "nsaids", "acei_arb", "ppi",
  "beta_blockers", "steroids", "antiarrhythmics",
  "first_potassium", "first_calcium", "first_heartrate",
  "vasopressor_6h", "first_mg_value",
  "first_lactate", "lactate_missing"), names(dat))

d <- dat[complete.cases(dat[, ps_covars]), ]
ps_fml <- as.formula(paste("mg_supp ~", paste(ps_covars, collapse = " + ")))
ps_fit <- glm(ps_fml, data = d, family = binomial())
d$ps <- pmax(pmin(fitted(ps_fit), 0.99), 0.01)
d$ow <- ifelse(d$mg_supp == 1, 1 - d$ps, d$ps)

# AC weights
if ("ac_group" %in% names(d)) {
  d$ac_trt <- NA_integer_
  d$ac_trt[d$ac_group == "mg_k"] <- 1L
  d$ac_trt[d$ac_group == "k_only"] <- 0L
  d_ac <- d[!is.na(d$ac_trt), ]
  ac_fml <- as.formula(paste("ac_trt ~", paste(ps_covars, collapse = "+")))
  ac_fit <- glm(ac_fml, data = d_ac, family = binomial())
  d_ac$ac_ps <- pmax(pmin(fitted(ac_fit), 0.99), 0.01)
  d_ac$ac_ow <- ifelse(d_ac$ac_trt == 1, 1 - d_ac$ac_ps, d_ac$ac_ps)
}

cat(sprintf("  N = %d (PS complete cases)\n", nrow(d)))

# ── Test each negative control ───────────────────────────────────────
cat(sprintf("\n%s\nCOMPLEXITY-SPECIFIC NEGATIVE CONTROLS\n%s\n",
            strrep("=", 60), strrep("=", 60)))

for (nc in c("nc_reoperation", "nc_transfusion")) {
  nice <- gsub("nc_", "", nc)
  n_events <- sum(d[[nc]], na.rm = TRUE)

  if (n_events < 5) {
    cat(sprintf("  %-20s SKIP (%d events)\n", nice, n_events))
    next
  }

  cat(sprintf("\n  ── %s (%d events, %.1f%%) ──\n",
              nice, n_events, 100 * mean(d[[nc]])))

  # Crude rates
  rate_trt  <- 100 * mean(d[[nc]][d$mg_supp == 1], na.rm = TRUE)
  rate_ctrl <- 100 * mean(d[[nc]][d$mg_supp == 0], na.rm = TRUE)
  cat(sprintf("    Crude: supplemented %.1f%% vs not %.1f%%\n",
              rate_trt, rate_ctrl))

  # All-patient OW
  tryCatch({
    d$.w <- d$ow
    fit <- glm(as.formula(paste(nc, "~ mg_supp")), data = d,
               weights = .w, family = quasibinomial())
    vc <- vcovCL(fit, cluster = d$hospitalid)
    ct <- coeftest(fit, vcov. = vc)
    or <- exp(ct[2, 1])
    lo <- exp(ct[2, 1] - 1.96 * ct[2, 2])
    hi <- exp(ct[2, 1] + 1.96 * ct[2, 2])
    p  <- 2 * pnorm(-abs(ct[2, 1] / ct[2, 2]))
    null_flag <- ifelse(p > 0.05, "NULL ✓", "SIGNIFICANT ✗")
    cat(sprintf("    All-patient OW: OR %.3f (%.3f-%.3f) P=%.4f  [%s]\n",
                or, lo, hi, p, null_flag))
  }, error = function(e) cat(sprintf("    OW failed: %s\n", e$message)))

  # AC OW
  if (exists("d_ac") && nc %in% names(d_ac) && sum(d_ac[[nc]], na.rm = TRUE) >= 3) {
    tryCatch({
      d_ac$.w <- d_ac$ac_ow
      fit_ac <- glm(as.formula(paste(nc, "~ ac_trt")), data = d_ac,
                     weights = .w, family = quasibinomial())
      vc_ac <- vcovHC(fit_ac, type = "HC1")
      ct_ac <- coeftest(fit_ac, vcov. = vc_ac)
      or_ac <- exp(ct_ac[2, 1])
      lo_ac <- exp(ct_ac[2, 1] - 1.96 * ct_ac[2, 2])
      hi_ac <- exp(ct_ac[2, 1] + 1.96 * ct_ac[2, 2])
      p_ac  <- 2 * pnorm(-abs(ct_ac[2, 1] / ct_ac[2, 2]))
      null_ac <- ifelse(p_ac > 0.05, "NULL ✓", "SIGNIFICANT ✗")
      cat(sprintf("    AC OW:          OR %.3f (%.3f-%.3f) P=%.4f  [%s]\n",
                  or_ac, lo_ac, hi_ac, p_ac, null_ac))
    }, error = function(e) cat(sprintf("    AC failed: %s\n", e$message)))
  }

  # Compare with AKI effect as reference
  cat(sprintf("    (Reference: AKI AC OR = 0.75, P=.02)\n"))
}

# ── Verdict ──────────────────────────────────────────────────────────
cat(sprintf("\n%s\nVERDICT\n%s\n", strrep("=", 60), strrep("=", 60)))
cat("  If BOTH null → complexity confounding unlikely to explain AKI finding\n")
cat("  If EITHER significant → residual complexity confounding plausible\n")
cat("  (Supplement as eTable alongside fracture negative control)\n")
