# ============================================================================
# did_covars.R — PS covariate specifications (multi-model)
#
# Models:
#   v3:       24 covars — chronic drug flags, keep steroids/K+
#   yan:      23 covars — original 28 minus {nsaids,ppi,acei_arb,beta_blockers,steroids}
#   original: 28 covars — all pre-t0 drugs + lactate
# ============================================================================

# ── v3: chronic drug revision (24 covariates) ────────────────────────────
PS_COVARS_V3 <- c(
  "age", "is_female", "bmi",
  "surg_cabg", "surg_valve", "surg_combined",
  "heart_failure", "hypertension", "diabetes", "ckd",
  "copd", "pvd", "stroke", "liver_disease",
  "egfr",
  "steroids",
  "ppi_chronic", "loop_diuretic_chronic", "acei_arb_chronic", "nsaid_chronic",
  "first_potassium", "first_calcium", "first_heartrate", "first_mg_value"
)

# ── yan: drop 5 drugs from original (23 covariates) ─────────────────────
# Dropped: nsaids, ppi, acei_arb, beta_blockers, steroids
# Kept: loop_diuretics, antiarrhythmics (+ all non-drug covariates)
PS_COVARS_YAN <- c(
  "age", "is_female", "bmi",
  "surg_cabg", "surg_valve", "surg_combined",
  "heart_failure", "hypertension", "diabetes", "ckd",
  "copd", "pvd", "stroke", "liver_disease",
  "egfr",
  "loop_diuretics", "antiarrhythmics",
  "first_potassium", "first_calcium", "first_heartrate",
  "first_mg_value", "first_lactate", "lactate_missing"
)

# ── original 28 ──────────────────────────────────────────────────────────
PS_COVARS_ORIGINAL <- c(
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

# ── Model selector ───────────────────────────────────────────────────────
MODEL_REGISTRY <- list(
  v3       = PS_COVARS_V3,
  yan      = PS_COVARS_YAN,
  original = PS_COVARS_ORIGINAL
)

select_model <- function(name) {
  name <- tolower(name)
  if (!name %in% names(MODEL_REGISTRY))
    stop(sprintf("Unknown model '%s'. Available: %s", name,
                 paste(names(MODEL_REGISTRY), collapse=", ")))
  covars <- MODEL_REGISTRY[[name]]
  cat(sprintf("  Model '%s': %d covariates\n", name, length(covars)))
  covars
}

COVAR_NICE_NAMES <- c(
  age="Age", is_female="Female sex", bmi="BMI",
  surg_cabg="CABG", surg_valve="Valve surgery", surg_combined="Combined surgery",
  heart_failure="Heart failure", hypertension="Hypertension", diabetes="Diabetes",
  ckd="CKD", copd="COPD", pvd="PVD", stroke="Stroke", liver_disease="Liver disease",
  egfr="eGFR", steroids="Steroids (pre-t0)",
  ppi_chronic="PPI (chronic)", loop_diuretic_chronic="Loop diuretic (chronic)",
  acei_arb_chronic="ACEi/ARB (chronic)", nsaid_chronic="NSAID (chronic)",
  loop_diuretics="Loop diuretics", nsaids="NSAIDs", acei_arb="ACEi/ARB",
  ppi="PPI", beta_blockers="Beta blockers", antiarrhythmics="Antiarrhythmics",
  first_potassium="Potassium", first_calcium="Calcium",
  first_heartrate="Heart rate", first_mg_value="Serum Mg",
  first_lactate="Lactate", lactate_missing="Lactate missing"
)

cat(sprintf("  did_covars.R loaded: v3=%d, yan=%d, original=%d\n",
            length(PS_COVARS_V3), length(PS_COVARS_YAN), length(PS_COVARS_ORIGINAL)))
