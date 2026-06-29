"""
00_config.py — Shared constants for Mg → AKI study (v5 design)

Used by: 01_etl.py, 02_psm.R (via did_covars.R), 03_hte.R
"""

import os

# ═══════════════════════════════════════════════════════════════════
# PATHS
# ═══════════════════════════════════════════════════════════════════
RESULTS = os.path.expanduser("~/mg_aki/results")
os.makedirs(RESULTS, exist_ok=True)

_FULL = os.path.expanduser("~/mg_aki/eicu-crd-2.0")
_DEMO = os.path.expanduser("~/mg_aki/eicu-collaborative-research-database-demo-2.0.1")
EICU_ROOT = _FULL if os.path.isdir(_FULL) else _DEMO

MIMIC_ROOT = os.path.expanduser("~/mg_aki/mimic-iv-3.1")
MIMIC_HOSP = os.path.join(MIMIC_ROOT, "hosp")
MIMIC_ICU = os.path.join(MIMIC_ROOT, "icu")

# ═══════════════════════════════════════════════════════════════════
# CLINICAL CONSTANTS
# ═══════════════════════════════════════════════════════════════════
MIN_AGE = 18
CR_MIN, CR_MAX = 0.1, 25.0
CR_POST_PLAUSIBLE_MAX = 15.0
BASELINE_CR_MAX = 4.0

# ═══════════════════════════════════════════════════════════════════
# MIMIC ITEM IDS
# ═══════════════════════════════════════════════════════════════════
LAB_CR_MIMIC = [50912, 52546]
LAB_MG_MIMIC = [50960]
LAB_K_MIMIC = [50971]
LAB_CA_MIMIC = [50893]
LAB_LAC_MIMIC = [50813]
VITAL_HR_MIMIC = [220045]

MG_SUPP_ITEMS_MIMIC = [222011, 227523]
K_SUPP_ITEMS_MIMIC = [225166, 225168, 222139, 227521, 227522]

# ── RRT detection (for KDIGO ≥3) ─────────────────────────────────
# MIMIC procedureevents: dialysis procedure item IDs
# VERIFY against d_items on Tempest before first run:
#   SELECT itemid, label FROM mimiciv_icu.d_items
#   WHERE LOWER(label) LIKE '%dialysis%' OR LOWER(label) LIKE '%crrt%'
RRT_PROCEDURE_ITEMS_MIMIC = [
    225441,  # Hemodialysis
    225802,  # Dialysis - CRRT
    225803,  # Dialysis - CVVHD
    225805,  # Peritoneal Dialysis
    225809,  # Dialysis - CVVHDF
    225955,  # Dialysis - SCUF
]
# MIMIC inputevents: CRRT replacement/dialysate items
# (presence = CRRT running)
# VERIFIED via probe_schema.py 2026-06-29:
#   227711 (Drain Removed) and 225183 (Current Goal) were WRONG — removed
#   230044 (Heparin Sodium CRRT-Prefilter) added as genuine CRRT marker
RRT_INPUT_ITEMS_MIMIC = [
    227536,  # KCl (CRRT) — 72,419 rows
    227525,  # Calcium Gluconate (CRRT) — 47,378 rows
    230044,  # Heparin Sodium (CRRT-Prefilter)
]

# MIMIC chartevents: Hemodialysis Output (presence = HD running)
# VERIFIED: adds 988 patients missed by procedureevents+inputevents
RRT_CHART_ITEMS_MIMIC = [
    226499,  # Hemodialysis Output
]

# All lab items needed for DuckDB filtering
ALL_LAB_ITEMS_MIMIC = set(
    LAB_CR_MIMIC + LAB_MG_MIMIC + LAB_K_MIMIC + LAB_CA_MIMIC + LAB_LAC_MIMIC
)

# ── Table 1 descriptive labs (not in PS model) ────────────────────
LAB_HGB_MIMIC = [51222, 50811]  # Hemoglobin (serum, blood gas)
LAB_WBC_MIMIC = [51301, 51300]  # WBC (auto, manual)
LAB_PLT_MIMIC = [51265]  # Platelet count
LAB_ALB_MIMIC = [50862]  # Albumin
TABLE1_LAB_ITEMS_MIMIC = set(
    LAB_HGB_MIMIC + LAB_WBC_MIMIC + LAB_PLT_MIMIC + LAB_ALB_MIMIC
)

# ═══════════════════════════════════════════════════════════════════
# TEXT PATTERNS
# ═══════════════════════════════════════════════════════════════════
MG_SUPP_PATTERNS = [
    "magnesium",
    "mag sulfate",
    "mgso4",
    "mag oxide",
    "mag gluconate",
    "mag hydroxide",
    "mag chloride",
]

CARDIAC_DX_PATTERNS = [
    "cabg",
    "valve",
    "cardiac surgery",
    "open heart",
    "coronary artery bypass",
    "aortic valve",
    "mitral valve",
    "cardiothoracic",
    "aortic surgery",
    "tricuspid",
    "pulmonic valve",
]
CARDIAC_UNIT_TYPES = {"CSICU", "CTICU", "CCU-CTICU"}

ESKD_PATTERNS = [
    "dialysis",
    "esrd",
    "end stage renal",
    "end-stage renal",
    "renal transplant",
    "kidney transplant",
]

# ═══════════════════════════════════════════════════════════════════
# ICD CODES
# ═══════════════════════════════════════════════════════════════════
CABG_ICD9 = ["3610", "3611", "3612", "3613", "3614", "3615", "3616", "3617", "3619"]
VALVE_ICD9 = [
    "3521",
    "3522",
    "3523",
    "3524",
    "3525",
    "3526",
    "3527",
    "3528",
    "3511",
    "3512",
    "3513",
    "3514",
]
CABG_ICD10 = ["0210", "0211", "0212", "0213"]
VALVE_ICD10 = ["02RF", "02RG", "02RH", "02RJ", "02QF", "02QG", "02QH", "02QJ"]
CVICU = "Cardiac Vascular Intensive Care Unit (CVICU)"

# ── Comorbidities ─────────────────────────────────────────────────
EICU_COMORB = {
    "heart_failure": ["heart failure", "chf", "cardiomyopathy"],
    "hypertension": ["hypertension"],
    "diabetes": ["diabetes"],
    "ckd": ["chronic kidney", "chronic renal", "ckd"],
    "copd": ["copd", "chronic obstructive", "emphysema"],
    "pvd": ["peripheral vascular", "pvd", "claudication"],
    "stroke": ["stroke", "cva", "cerebrovascular"],
    "liver_disease": ["cirrhosis", "hepatitis", "liver disease", "liver failure"],
}

MIMIC_COMORB_ICD = {
    "heart_failure": {9: ["4280", "4281", "4289", "428"], 10: ["I50"]},
    "hypertension": {
        9: ["401", "402", "403", "404", "405"],
        10: ["I10", "I11", "I12", "I13", "I15"],
    },
    "diabetes": {9: ["250"], 10: ["E08", "E09", "E10", "E11", "E12", "E13"]},
    "ckd": {9: ["585", "586"], 10: ["N18", "N19"]},
    "copd": {
        9: ["490", "491", "492", "493", "494", "496"],
        10: ["J40", "J41", "J42", "J43", "J44", "J45", "J47"],
    },
    "pvd": {9: ["4431", "4439", "4471"], 10: ["I73"]},
    "stroke": {
        9: ["430", "431", "432", "433", "434", "435", "436"],
        10: ["I60", "I61", "I62", "I63", "I64", "I65", "I66", "G45"],
    },
    "liver_disease": {
        9: ["571"],
        10: ["K70", "K71", "K72", "K73", "K74", "K75", "K76"],
    },
}

# ── Secondary outcome codes ───────────────────────────────────────
AF_PATTERNS_EICU = [
    "atrial fibrillation",
    "atrial flutter",
    "a-fib",
    "afib",
    "a fib",
    "new onset af",
]
AF_PRIOR_PATTERNS_EICU = [
    "atrial fibrillation",
    "atrial flutter",
    "a-fib",
    "afib",
    "chronic af",
]
AF_ICD9 = ["42731", "42732"]
AF_ICD10_PREFIX = "I48"

ENCEPH_PATTERNS_EICU = [
    "encephalopathy",
    "delirium",
    "altered mental",
    "acute confusional",
    "metabolic encephalopathy",
]
ENCEPH_ICD9 = ["3481", "3489", "2930", "2931"]
ENCEPH_ICD10_PREFIX = ["G93", "F05"]

VARR_PATTERNS_EICU = [
    "ventricular tachycardia",
    "ventricular fibrillation",
    "v-tach",
    "v-fib",
    "vtach",
    "vfib",
    "cardiac arrest",
]
VARR_ICD9 = ["4271", "42741", "42742"]
VARR_ICD10 = ["I472", "I490"]

# ── Chronic drug classes (home medications) ───────────────────────
CHRONIC_DRUG_CLASSES = {
    "ppi_chronic": [
        "omeprazole",
        "pantoprazole",
        "lansoprazole",
        "esomeprazole",
        "rabeprazole",
    ],
    "loop_diuretic_chronic": ["furosemide", "bumetanide", "torsemide", "lasix"],
    "acei_arb_chronic": [
        "lisinopril",
        "enalapril",
        "ramipril",
        "captopril",
        "losartan",
        "valsartan",
        "irbesartan",
        "candesartan",
        "olmesartan",
        "telmisartan",
    ],
    "nsaid_chronic": [
        "ibuprofen",
        "ketorolac",
        "naproxen",
        "diclofenac",
        "celecoxib",
        "indomethacin",
        "meloxicam",
    ],
}
ORAL_ROUTE_RE = r"oral|po\b|tablet|capsule|cap\b|tab\b"

ESKD_ICD = {
    9: ["5856", "V4511", "V560", "V561", "V562"],
    10: ["N186", "Z491", "Z492", "Z9911", "Z940"],
}

# ═══════════════════════════════════════════════════════════════════
# PS COVARIATES (for did_covars.R and R scripts)
# ═══════════════════════════════════════════════════════════════════

# Time-invariant: extracted in ETL, fixed per patient
PS_TIME_INVARIANT = [
    "age",
    "is_female",
    "bmi",
    "surg_cabg",
    "surg_valve",
    "surg_combined",
    "heart_failure",
    "hypertension",
    "diabetes",
    "ckd",
    "copd",
    "pvd",
    "stroke",
    "liver_disease",
    "egfr",
    "ppi_chronic",
    "loop_diuretic_chronic",
    "acei_arb_chronic",
    "nsaid_chronic",
]

# Time-varying: computed at match time in R from did_labs_all
PS_TIME_VARYING_LABS = [
    "mg_value",
    "potassium",
    "calcium",
    "lactate",
    "lactate_missing",
    "heartrate",
]

# eICU lab name patterns (for extracting from lab table)
EICU_LAB_PATTERNS = {
    "magnesium": ["magnesium"],
    "potassium": ["potassium"],
    "calcium": ["calcium"],
    "lactate": ["lactate"],
}

# eICU Table 1 descriptive lab patterns (not in PS model)
EICU_TABLE1_LAB_PATTERNS = {
    "hemoglobin": ["hgb", "hemoglobin"],
    "wbc": ["wbc"],
    "platelets": ["platelet"],
    "albumin": ["albumin"],
}

# eICU RRT detection patterns (treatment table + intakeOutput)
EICU_RRT_TREATMENT_PATTERNS = [
    "dialysis",
    "crrt",
    "cvvh",
    "cvvhd",
    "cvvhdf",
    "hemodialysis",
    "scuf",
    "ultrafiltration",
]

# Output column names for did_all_{db}.csv
ALL_PATIENTS_COLS = [
    "pid",
    "hadm_id",
    "treated",
    "mg_offset_h",
    "mg_offset_min",
    "icu_discharge_h",
    "icu_outcome",
    "age",
    "is_female",
    "bmi",
    "surgery_type",
    "surg_cabg",
    "surg_valve",
    "surg_combined",
    *EICU_COMORB.keys(),
    "first_cr",
    "egfr",
    *CHRONIC_DRUG_CLASSES.keys(),
    "hosp_mortality",
    "poaf",
    "encephalopathy_delirium",
    "transfusion",
    "reintubation",
    "poaf_icd",
    "encephalopathy_icd",
    "vent_arrhythmia",
    "rrt_offset_h",
    "has_rrt",
    "death_offset_h",
]
