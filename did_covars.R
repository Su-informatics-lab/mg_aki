# ============================================================================
# did_covars.R — Final PS covariate specifications
#
# Primary: spec241 (21 covars) — labs only, no drug covariates
#   Root cause: ICU drug records = protocol co-prescriptions, not patient
#   variation. Labs (K+, Ca, Mg, lactate) capture patient physiology directly.
#   256-spec AIPW sweep: ALL specs concordant negative, drugs irrelevant.
#
# Method: AIPW primary, sIPTW_DR sensitivity, ICU-time anchor
# ============================================================================

# ── Primary: 21 covariates (spec 241, labs only) ─────────────────────────
# BASE(16) + K+ + Ca + lactate + lactate_missing + Mg
PS_PRIMARY <- c(
  "age", "is_female", "bmi",
  "surg_cabg", "surg_valve", "surg_combined",
  "heart_failure", "hypertension", "diabetes", "ckd",
  "copd", "pvd", "stroke", "liver_disease",
  "egfr", "first_heartrate",
  "first_potassium", "first_calcium",
  "first_lactate", "lactate_missing",
  "first_mg_value"
)

# ── Sensitivity A: 18 covars (spec 49, K+ + Ca only) ────────────────────
PS_SENS_A <- c(
  "age", "is_female", "bmi",
  "surg_cabg", "surg_valve", "surg_combined",
  "heart_failure", "hypertension", "diabetes", "ckd",
  "copd", "pvd", "stroke", "liver_disease",
  "egfr", "first_heartrate",
  "first_potassium", "first_calcium"
)

# ── Sensitivity B: 20 covars (spec 113, K+ + Ca + lactate, no Mg) ───────
PS_SENS_B <- c(
  "age", "is_female", "bmi",
  "surg_cabg", "surg_valve", "surg_combined",
  "heart_failure", "hypertension", "diabetes", "ckd",
  "copd", "pvd", "stroke", "liver_disease",
  "egfr", "first_heartrate",
  "first_potassium", "first_calcium",
  "first_lactate", "lactate_missing"
)

# ── Sensitivity C: 26 covars (spec 244, + chronic drugs + steroids) ─────
PS_SENS_C <- c(
  "age", "is_female", "bmi",
  "surg_cabg", "surg_valve", "surg_combined",
  "heart_failure", "hypertension", "diabetes", "ckd",
  "copd", "pvd", "stroke", "liver_disease",
  "egfr", "first_heartrate",
  "ppi_chronic", "loop_diuretic_chronic", "acei_arb_chronic", "nsaid_chronic",
  "steroids",
  "first_potassium", "first_calcium",
  "first_lactate", "lactate_missing",
  "first_mg_value"
)

# ── Sensitivity D: 16 covars (spec 1, base only) ────────────────────────
PS_SENS_D <- c(
  "age", "is_female", "bmi",
  "surg_cabg", "surg_valve", "surg_combined",
  "heart_failure", "hypertension", "diabetes", "ckd",
  "copd", "pvd", "stroke", "liver_disease",
  "egfr", "first_heartrate"
)

# ── Original 28 (for legacy comparison) ──────────────────────────────────
PS_ORIGINAL <- c(
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

MODEL_REGISTRY <- list(
  primary  = PS_PRIMARY,
  sens_a   = PS_SENS_A,
  sens_b   = PS_SENS_B,
  sens_c   = PS_SENS_C,
  sens_d   = PS_SENS_D,
  original = PS_ORIGINAL
)

select_model <- function(name) {
  name <- tolower(name)
  if (!name %in% names(MODEL_REGISTRY))
    stop(sprintf("Unknown model '%s'. Available: %s",
                 name, paste(names(MODEL_REGISTRY), collapse=", ")))
  covars <- MODEL_REGISTRY[[name]]
  cat(sprintf("  Model '%s': %d covariates\n", name, length(covars)))
  covars
}

COVAR_NICE_NAMES <- c(
  age="Age", is_female="Female sex", bmi="BMI",
  surg_cabg="CABG", surg_valve="Valve surgery", surg_combined="Combined surgery",
  heart_failure="Heart failure", hypertension="Hypertension", diabetes="Diabetes",
  ckd="CKD", copd="COPD", pvd="PVD", stroke="Stroke", liver_disease="Liver disease",
  egfr="eGFR", first_heartrate="Heart rate",
  first_potassium="Potassium", first_calcium="Calcium",
  first_lactate="Lactate", lactate_missing="Lactate missing",
  first_mg_value="Serum Mg",
  steroids="Steroids", ppi_chronic="PPI (chronic)",
  loop_diuretic_chronic="Loop diuretic (chronic)",
  acei_arb_chronic="ACEi/ARB (chronic)", nsaid_chronic="NSAID (chronic)",
  loop_diuretics="Loop diuretics", nsaids="NSAIDs", acei_arb="ACEi/ARB",
  ppi="PPI", beta_blockers="Beta blockers", antiarrhythmics="Antiarrhythmics"
)

cat(sprintf("  did_covars.R: primary=%d, sensA=%d, sensB=%d, sensC=%d, sensD=%d\n",
            length(PS_PRIMARY), length(PS_SENS_A), length(PS_SENS_B),
            length(PS_SENS_C), length(PS_SENS_D)))
