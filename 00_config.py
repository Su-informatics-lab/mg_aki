#!/usr/bin/env python3
"""
Mg Reserve → Cardiac Surgery AKI  (eICU-CRD)
00_config.py — shared constants, paths, phenotype parameters

Usage: imported by 01_etl.py, 02_analysis.py
"""

import os

# =====================================================================
# PATHS — Tempest (g91p721)
# =====================================================================
_FULL = os.path.expanduser("~/mg_aki/eicu-crd-2.0")
_DEMO = os.path.expanduser("~/mg_aki/eicu-collaborative-research-database-demo-2.0.1")
DATA_ROOT = _FULL if os.path.isdir(_FULL) else _DEMO
RESULTS = os.path.expanduser("~/mg_aki/results")
os.makedirs(RESULTS, exist_ok=True)

# =====================================================================
# COHORT PARAMETERS
# =====================================================================
MIN_AGE = 18
AGE_CAP = 90  # "> 89" coded as 90

# Cardiac surgery identification patterns
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

# =====================================================================
# EXPOSURE: MAGNESIUM
# =====================================================================
MG_LABNAMES = ["magnesium"]  # case-insensitive match via str.contains
MG_PLAUSIBLE_MIN = 0.5  # mg/dL
MG_PLAUSIBLE_MAX = 10.0
MG_WINDOW_HOURS = 6  # first Mg within this window post-ICU-admit
MG_WINDOW_MIN = MG_WINDOW_HOURS * 60  # 360 minutes

# Clinical cutpoints (mg/dL)
MG_HYPO_THRESHOLD = 2.0  # for TTE eligibility (Analysis B)
MG_NORMAL_RANGE = (1.8, 2.3)

# =====================================================================
# OUTCOME: AKI
# =====================================================================
CR_LABNAMES = ["creatinine"]  # case-insensitive
CR_PLAUSIBLE_MIN = 0.1  # mg/dL
CR_PLAUSIBLE_MAX = 25.0

# -- Baseline Cr strategies --
# Primary (Approach 2): pre-ICU Cr if available, else first post-ICU
BASELINE_PRE_ICU_WINDOW_MIN = -720  # 12h before ICU admit
BASELINE_PRE_ICU_WINDOW_MAX = 0  # ICU admission

# Sensitivity (Approach 1 / Eadon): lowest Cr within 48h
BASELINE_EADON_WINDOW_MIN = -60  # 1h before ICU (generous)
BASELINE_EADON_WINDOW_MAX = 2880  # 48h post-admit

# -- AKI thresholds (KDIGO Cr-based) --
AKI_RATIO_STAGE1 = 1.5  # primary: ≥1.5x baseline
AKI_DELTA_48H = 0.3  # secondary: ≥0.3 mg/dL within 48h
AKI_RATIO_STAGE2 = 2.0
AKI_RATIO_STAGE3 = 3.0
AKI_CR_ABSOLUTE = 4.0  # Stage 3 absolute threshold

# Follow-up window for AKI detection (minutes from ICU admit)
AKI_WINDOW_7D_MIN = 7 * 24 * 60  # 10080
AKI_WINDOW_48H_MIN = 48 * 60  # 2880

# =====================================================================
# EXCLUSIONS
# =====================================================================
ESKD_DX_PATTERNS = [
    "dialysis",
    "esrd",
    "end stage renal",
    "end-stage renal",
    "renal transplant",
    "kidney transplant",
]
BASELINE_CR_MAX = 4.0  # exclude baseline Cr ≥ 4.0 mg/dL

# =====================================================================
# TREATMENT (TTE Analysis B): Mg supplementation
# =====================================================================
MG_SUPP_DRUG_PATTERNS = [
    "magnesium",
    "mag sulfate",
    "mgso4",
    "mag oxide",
    "mag gluconate",
    "mag hydroxide",
    "mag chloride",
]
MG_SUPP_GRACE_HOURS = 6  # supplementation within 6h of time zero
MG_SUPP_GRACE_MIN = MG_SUPP_GRACE_HOURS * 60

# =====================================================================
# COVARIATES
# =====================================================================
# Comorbidity keywords for pastHistory table
COMORBIDITY_KEYWORDS = {
    "chf": ["chf", "heart failure", "cardiomyopathy", "hf "],
    "hypertension": ["hypertension", "htn", "high blood pressure"],
    "diabetes": ["diabetes", "dm ", "insulin dependent", "iddm", "niddm"],
    "ckd": [
        "renal insufficiency",
        "ckd",
        "chronic kidney",
        "renal failure",
        "kidney disease",
    ],
    "copd": ["copd", "chronic obstructive", "emphysema"],
    "pvd": ["peripheral vascular", "pvd", "claudication"],
    "stroke": ["stroke", "cva", "cerebrovascular", "tia"],
    "liver": ["cirrhosis", "hepatic failure", "liver disease", "hepatitis"],
    "afib": ["atrial fibrillation", "afib", "a-fib", "af "],
}

# Nephrotoxin patterns in medication/admissionDrug
NEPHROTOXIN_CLASSES = {
    "ppi": [
        "omeprazole",
        "pantoprazole",
        "lansoprazole",
        "esomeprazole",
        "rabeprazole",
    ],
    "nsaid": [
        "ibuprofen",
        "naproxen",
        "diclofenac",
        "meloxicam",
        "celecoxib",
        "indomethacin",
        "ketorolac",
    ],
    "acei_arb": [
        "lisinopril",
        "enalapril",
        "ramipril",
        "captopril",
        "losartan",
        "valsartan",
        "irbesartan",
        "olmesartan",
    ],
    "loop_diuretic": [
        "furosemide",
        "bumetanide",
        "torsemide",
        "lasix",
    ],
}

# Surgery type classification from apacheAdmissionDx
SURGERY_TYPE_MAP = {
    "cabg": ["cabg", "coronary artery bypass"],
    "valve": ["valve", "aortic valve", "mitral valve", "tricuspid", "pulmonic"],
    "combined": [],  # assigned if both CABG + valve keywords match
}

# =====================================================================
# POSITIVE CONTROL: POAF (postoperative atrial fibrillation)
# =====================================================================
# Diagnosis table patterns for NEW-ONSET POAF
POAF_DX_PATTERNS = [
    "atrial fibrillation",
    "atrial flutter",
    "a fib",
    "afib",
    "a-fib",
    "a. fib",
    "rapid afib",
]
# Pre-existing AF patterns in pastHistory (for exclusion)
PREEXISTING_AF_PATTERNS = [
    "atrial fibrillation",
    "atrial flutter",
    "afib",
    "a-fib",
    "chronic atrial",
    "paroxysmal atrial",
]
POAF_WINDOW_MIN = 7 * 24 * 60  # 7 days in minutes

# =====================================================================
# EFFECT MODIFIERS: β-blockers, steroids
# =====================================================================
BETA_BLOCKER_PATTERNS = [
    "metoprolol",
    "atenolol",
    "propranolol",
    "carvedilol",
    "bisoprolol",
    "labetalol",
    "nadolol",
    "sotalol",
    "esmolol",
    "timolol",
    "nebivolol",
    "betaxolol",
]

STEROID_PATTERNS = [
    "methylprednisolone",
    "dexamethasone",
    "hydrocortisone",
    "prednisone",
    "prednisolone",
    "solumedrol",
]

# =====================================================================
# NEGATIVE CONTROLS (bias calibration)
# =====================================================================
NEGATIVE_CONTROL_DX = {
    "fracture": ["fracture"],
    "skin_infection": ["cellulitis", "skin infection", "abscess"],
    "uti": ["urinary tract infection", "uti", "cystitis"],
}

# =====================================================================
# EXPLORATORY: Neurological outcomes (descriptive)
# =====================================================================
NEURO_DX_PATTERNS = {
    "delirium": ["delirium", "confusion", "altered mental", "agitation"],
    "seizure": ["seizure", "convulsion", "epilep"],
    "stroke_postop": [
        "stroke",
        "cerebrovascular accident",
        "cva",
        "intracranial hemorrhage",
        "cerebral infarct",
    ],
    "encephalopathy": ["encephalopathy", "metabolic encephalop"],
}
RANDOM_SEED = 42
PS_CALIPER = 0.2  # SD of logit(PS)
IPTW_TRIM_PERCENTILE = (1, 99)  # truncate at 1st/99th percentile
SMD_THRESHOLD = 0.1  # acceptable standardized mean difference

# Sensitivity analysis: Mg window variants (hours)
MG_WINDOW_SENSITIVITY = [1, 2, 3, 6, 12]

# Negative control outcomes (not plausibly affected by Mg in 7d)
NEGATIVE_CONTROL_DX = {
    "fracture": ["fracture"],
    "skin_infection": ["cellulitis", "skin infection", "abscess"],
    "uti": ["urinary tract infection", "uti", "cystitis"],
}

# =====================================================================
# COMPOSITE POAF PHENOTYPE (v2 — multi-source)
# =====================================================================
# Treatment table patterns for AF management
POAF_TREATMENT_PATTERNS = [
    "cardioversion",
    "rhythm control",
    "defibrillat",
    "atrial fibrillation",
    "atrial flutter",
]

# NEW postop antiarrhythmic medications (if NOT in admissionDrug → new AF)
POAF_MED_PATTERNS = [
    "amiodarone",
    "sotalol",
    "flecainide",
    "propafenone",
    "ibutilide",
    "dofetilide",
    "vernakalant",
    "dronedarone",
]

# InfusionDrug patterns for AF treatment
POAF_INFUSION_PATTERNS = [
    "amiodarone",
    "cardizem",
    "diltiazem",
]

# Pre-operative antiarrhythmic flag (from admissionDrug)
PREOP_ANTIARRHYTHMIC_PATTERNS = [
    "amiodarone",
    "sotalol",
    "flecainide",
    "propafenone",
    "dofetilide",
    "dronedarone",
    "digoxin",
]

# Potassium lab for electrolyte covariate
K_LABNAMES = ["potassium"]
K_PLAUSIBLE_MIN = 1.5
K_PLAUSIBLE_MAX = 8.0

# =====================================================================
# POSITIVE CONTROL: Ventricular arrhythmia (RR 0.52 from 10 RCTs)
# =====================================================================
VENT_ARRHYTHMIA_DX = [
    "ventricular tachycardia",
    "ventricular fibrillation",
    "v-tach",
    "v-fib",
    "vtach",
    "vfib",
    "v tach",
    "v fib",
    "torsade",
    "torsades",
    "ventricular arrhythmia",
]
