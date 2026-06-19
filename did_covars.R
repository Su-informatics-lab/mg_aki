# ============================================================================
# did_covars.R — Shared PS covariate specification (v2: chronic drug revision)
#
# Source this from all DiD analysis scripts:
#   source("did_covars.R")
#
# Revision per Dr. Yan (2026-06-19):
#   - Nephrotoxic drugs split into chronic (home med) vs ICU protocol
#   - 4 chronic flags REPLACE combined flags in primary model
#   - Steroids kept as pre-t0 (intraoperative + home)
#   - Beta blockers, antiarrhythmics: DROPPED (positivity / collider)
#   - Lactate + lactate_missing: DROPPED (missingness = protocol noise)
#
# Chronic flag source:
#   eICU:  admissionDrug table (medication reconciliation)
#   MIMIC: prescriptions pre-admission + oral route
#   Created by: 00b_patch_chronic_drugs.py
# ============================================================================

# ── Primary model: 24 covariates ──────────────────────────────────────────
# Clinically parsimonious per Dr. Yan review
PS_COVARS_PRIMARY <- c(
  # Demographics (3)
  "age", "is_female", "bmi",
  # Surgery type (3)
  "surg_cabg", "surg_valve", "surg_combined",
  # Comorbidities (8)
  "heart_failure", "hypertension", "diabetes", "ckd",
  "copd", "pvd", "stroke", "liver_disease",
  # Renal function (1)
  "egfr",
  # Drugs — chronic/home use only (5)
  "steroids",              # pre-t0: intraop + home (already balanced, SMD ~0.1)
  "ppi_chronic",           # home PPI (interstitial nephritis, hypoMg)
  "loop_diuretic_chronic", # home diuretic (HF severity marker)
  "acei_arb_chronic",      # home RAAS blockade (renal hemodynamics)
  "nsaid_chronic",         # home nephrotoxin
  # Labs (4)
  "first_potassium", "first_calcium", "first_heartrate",
  "first_mg_value"
)

# ── Sensitivity A: + combined (pre-t0) ACEi/ARB and NSAIDs ───────────────
# Tests whether the timing distinction matters for nephrotoxic drugs
PS_COVARS_SENS_A <- c(
  PS_COVARS_PRIMARY,
  "acei_arb",   # combined pre-t0 (chronic + ICU)
  "nsaids"      # combined pre-t0 (chronic + ICU)
)

# ── Sensitivity B: + lactate ──────────────────────────────────────────────
# Tests whether lactate missingness was distorting the PS
PS_COVARS_SENS_B <- c(
  PS_COVARS_PRIMARY,
  "first_lactate", "lactate_missing"
)

# ── Sensitivity C: original 28-covariate model (full robustness) ─────────
# All pre-t0 drugs + lactate (same as v1 for comparison)
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

# ── Nice names for figures / tables ───────────────────────────────────────
COVAR_NICE_NAMES <- c(
  age="Age", is_female="Female sex", bmi="BMI",
  surg_cabg="CABG", surg_valve="Valve surgery", surg_combined="Combined surgery",
  heart_failure="Heart failure", hypertension="Hypertension", diabetes="Diabetes",
  ckd="CKD", copd="COPD", pvd="PVD", stroke="Stroke", liver_disease="Liver disease",
  egfr="eGFR",
  steroids="Steroids (pre-t0)",
  ppi_chronic="PPI (chronic)", loop_diuretic_chronic="Loop diuretic (chronic)",
  acei_arb_chronic="ACEi/ARB (chronic)", nsaid_chronic="NSAID (chronic)",
  # Original combined flags (for sensitivity)
  loop_diuretics="Loop diuretics", nsaids="NSAIDs", acei_arb="ACEi/ARB",
  ppi="PPI", beta_blockers="Beta blockers", antiarrhythmics="Antiarrhythmics",
  # Labs
  first_potassium="Potassium", first_calcium="Calcium",
  first_heartrate="Heart rate", first_mg_value="Serum Mg",
  first_lactate="Lactate", lactate_missing="Lactate missing"
)

# ── Helper: select covariates available in data ───────────────────────────
select_available <- function(covars, data_names) {
  avail <- intersect(covars, data_names)
  missing <- setdiff(covars, data_names)
  if (length(missing) > 0) {
    cat(sprintf("  WARNING: %d covariates not in data: %s\n",
                length(missing), paste(missing, collapse=", ")))
  }
  avail
}

cat(sprintf("  PS covariate config loaded: primary=%d, sensA=%d, sensB=%d, original=%d\n",
            length(PS_COVARS_PRIMARY), length(PS_COVARS_SENS_A),
            length(PS_COVARS_SENS_B), length(PS_COVARS_ORIGINAL)))
