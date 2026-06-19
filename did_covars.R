# ============================================================================
# did_covars.R — Shared PS covariate specification (v4)
#
# v4: Dropped first_potassium (K+ repletion protocol confounds Mg treatment)
#     Dropped steroids (intraop dexamethasone protocol-driven, worst SMD eICU)
#     Primary model: 22 covariates
#
# v3: Nephrotoxic drugs -> chronic (home med) flags; dropped BB/AA/lactate
# ============================================================================

PS_COVARS_PRIMARY <- c(
  "age", "is_female", "bmi",
  "surg_cabg", "surg_valve", "surg_combined",
  "heart_failure", "hypertension", "diabetes", "ckd",
  "copd", "pvd", "stroke", "liver_disease",
  "egfr",
  "ppi_chronic", "loop_diuretic_chronic", "acei_arb_chronic", "nsaid_chronic",
  "first_calcium", "first_heartrate", "first_mg_value"
)

# Sensitivity A: +steroids +K+ (= v3 primary, 24 covars)
PS_COVARS_SENS_A <- c(PS_COVARS_PRIMARY, "steroids", "first_potassium")

# Sensitivity B: +lactate (26 covars)
PS_COVARS_SENS_B <- c(PS_COVARS_SENS_A, "first_lactate", "lactate_missing")

# Sensitivity C: original 28
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

cat(sprintf("  PS covariate config loaded: primary=%d, sensA=%d, sensB=%d, original=%d\n",
            length(PS_COVARS_PRIMARY), length(PS_COVARS_SENS_A),
            length(PS_COVARS_SENS_B), length(PS_COVARS_ORIGINAL)))
