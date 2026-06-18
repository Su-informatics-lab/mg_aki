#!/usr/bin/env python3
"""
20_did_etl.py — DiD cohort construction for Mg → AKI study (v1)

  Design change from TTE (01_etl.py):
    - Time zero = first postop IV Mg administration (patient-specific)
    - Cr_pre = latest Cr between ICU admission and IV Mg (no hosp fallback)
    - Cr_post = first Cr after IV Mg at multiple windows (6-12h, 6-24h, 6-48h)
    - Outcome = ΔCr (continuous) not KDIGO (binary)
    - Controls = patients who never received postop IV Mg
    - Matching (1:4 PSM + temporal alignment) done in a separate R script

  Outputs:
    results/20_did_treated_eicu.csv     — treated arm with Cr_pre, Cr_post, all Cr times
    results/20_did_control_eicu.csv     — control arm with all Cr measurements inventoried
    results/20_did_treated_mimic.csv    — same for MIMIC
    results/20_did_control_mimic.csv    — same for MIMIC
    results/20_did_consort.csv          — CONSORT numbers for both databases

  Run: python 20_did_etl.py
"""

import os
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ═══════════════════════════════════════════════════════════════════
# PATHS & CONSTANTS
# ═══════════════════════════════════════════════════════════════════
RESULTS = os.path.expanduser("~/mg_aki/results")
os.makedirs(RESULTS, exist_ok=True)

# eICU
_FULL = os.path.expanduser("~/mg_aki/eicu-crd-2.0")
_DEMO = os.path.expanduser("~/mg_aki/eicu-collaborative-research-database-demo-2.0.1")
EICU_ROOT = _FULL if os.path.isdir(_FULL) else _DEMO

# MIMIC
MIMIC_ROOT = os.path.expanduser("~/mg_aki/mimic-iv-3.1")
MIMIC_HOSP = os.path.join(MIMIC_ROOT, "hosp")
MIMIC_ICU = os.path.join(MIMIC_ROOT, "icu")

# Shared constants
MIN_AGE = 18
CR_MIN, CR_MAX = 0.1, 25.0
MG_LAB_MIN, MG_LAB_MAX = 0.5, 10.0
BASELINE_CR_MAX = 4.0

# Cardiac surgery identification
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

# IV Mg drug patterns (eICU)
MG_SUPP_PATTERNS = [
    "magnesium",
    "mag sulfate",
    "mgso4",
    "mag oxide",
    "mag gluconate",
    "mag hydroxide",
    "mag chloride",
]

# IV Mg item IDs (MIMIC)
MG_SUPP_ITEMS_MIMIC = [222011, 227523]

# Cr outcome windows (hours after IV Mg)
CR_POST_WINDOWS = {
    "6_12h": (6, 12),
    "6_24h": (6, 24),  # primary
    "6_48h": (6, 48),
    "6_72h": (6, 72),
    "0_24h": (0, 24),  # sensitivity: no floor
    "0_48h": (0, 48),
}

# Comorbidity ICD codes (MIMIC)
COMORB_ICD = {
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

# eICU comorbidity patterns (pastHistory)
EICU_COMORB_PATTERNS = {
    "heart_failure": ["heart failure", "chf", "cardiomyopathy"],
    "hypertension": ["hypertension"],
    "diabetes": ["diabetes"],
    "ckd": ["chronic kidney", "chronic renal", "ckd"],
    "copd": ["copd", "chronic obstructive", "emphysema"],
    "pvd": ["peripheral vascular", "pvd", "claudication"],
    "stroke": ["stroke", "cva", "cerebrovascular"],
    "liver_disease": ["cirrhosis", "hepatitis", "liver disease", "liver failure"],
}

# Nephrotoxin / perioperative drug patterns (eICU)
DRUG_CLASSES_EICU = {
    "loop_diuretics": ["furosemide", "bumetanide", "torsemide", "lasix"],
    "nsaids": ["ibuprofen", "ketorolac", "naproxen", "diclofenac", "celecoxib"],
    "acei_arb": [
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
    "ppi": [
        "omeprazole",
        "pantoprazole",
        "lansoprazole",
        "esomeprazole",
        "rabeprazole",
    ],
    "beta_blockers": [
        "metoprolol",
        "atenolol",
        "propranolol",
        "carvedilol",
        "labetalol",
        "bisoprolol",
        "esmolol",
        "nadolol",
    ],
    "steroids": [
        "methylprednisolone",
        "hydrocortisone",
        "dexamethasone",
        "prednisone",
        "prednisolone",
        "solumedrol",
    ],
    "antiarrhythmics": [
        "amiodarone",
        "lidocaine",
        "procainamide",
        "flecainide",
        "sotalol",
        "dronedarone",
        "digoxin",
    ],
}

# ESKD patterns
ESKD_PATTERNS = [
    "dialysis",
    "esrd",
    "end stage renal",
    "end-stage renal",
    "renal transplant",
    "kidney transplant",
]


# ═══════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════
def gz(p):
    """Return .csv.gz path if exists, else .csv."""
    return p if os.path.exists(p) else p.replace(".csv.gz", ".csv")


def matches_any(series, patterns):
    """Case-insensitive multi-pattern match on a string Series."""
    s = series.astype(str).str.lower()
    mask = pd.Series(False, index=series.index)
    for p in patterns:
        mask |= s.str.contains(p.lower(), na=False)
    return mask


def compute_egfr(cr, age, is_female):
    """CKD-EPI 2021 (race-free)."""
    cr = np.asarray(cr, dtype=np.float64)
    age = np.asarray(age, dtype=np.float64)
    fem = np.asarray(is_female, dtype=bool)
    kappa = np.where(fem, 0.7, 0.9)
    alpha = np.where(fem, -0.241, -0.302)
    ratio = cr / kappa
    return (
        142
        * np.power(np.minimum(ratio, 1.0), alpha)
        * np.power(np.maximum(ratio, 1.0), -1.200)
        * np.power(0.9938, age)
        * np.where(fem, 1.012, 1.0)
    )


def pct(n, total):
    return f"{n:,} ({100*n/max(total,1):.1f}%)"


def describe_col(s, name, indent="    "):
    """Print summary stats for a numeric Series."""
    valid = s.dropna()
    if len(valid) == 0:
        print(f"{indent}{name}: all missing")
        return
    print(
        f"{indent}{name}: n={len(valid)}, "
        f"mean={valid.mean():.2f}, median={valid.median():.2f}, "
        f"IQR=[{valid.quantile(0.25):.2f}–{valid.quantile(0.75):.2f}], "
        f"range=[{valid.min():.2f}–{valid.max():.2f}]"
    )


# ═══════════════════════════════════════════════════════════════════
# eICU ETL
# ═══════════════════════════════════════════════════════════════════
def run_eicu():
    SEP = "=" * 70
    print(f"\n{SEP}")
    print("eICU-CRD: DiD Cohort Construction")
    print(f"  Data: {EICU_ROOT}")
    print(SEP)

    consort = {}

    # ── Load tables ──────────────────────────────────────────────
    def load(name, **kw):
        for ext in [".csv.gz", ".csv"]:
            p = os.path.join(EICU_ROOT, name + ext)
            if os.path.exists(p):
                df = pd.read_csv(p, low_memory=False, **kw)
                df.columns = df.columns.str.lower()
                print(f"  Loaded {name}: {len(df):,} rows")
                return df
        print(f"  WARNING: {name} not found")
        return pd.DataFrame()

    patient = load("patient")
    lab = load("lab")
    med = load(
        "medication",
        usecols=[
            "patientunitstayid",
            "drugstartoffset",
            "drugstopoffset",
            "drugname",
            "drugordercancelled",
            "routeadmin",
            "dosage",
        ],
    )
    inf = load("infusionDrug")
    diag = load("diagnosis")
    pasthx = load("pastHistory")
    treatment = load("treatment")

    consort["total_icu"] = len(patient)
    print(f"\n  Total ICU admissions: {len(patient):,}")

    # ── Cardiac surgery identification ───────────────────────────
    cardiac_mask = matches_any(
        patient.apacheadmissiondx, CARDIAC_DX_PATTERNS
    ) | patient.unittype.isin(CARDIAC_UNIT_TYPES)
    cardiac = patient[cardiac_mask].copy()

    # Age filter + first stay
    cardiac["age_num"] = pd.to_numeric(
        cardiac["age"].astype(str).str.replace(">", ""), errors="coerce"
    ).fillna(90)
    cardiac = cardiac[cardiac.age_num >= MIN_AGE]
    cardiac = (
        cardiac.sort_values("hospitaladmitoffset")
        .groupby("uniquepid")
        .first()
        .reset_index()
    )
    pids = set(cardiac.patientunitstayid)
    consort["cardiac_adult_first"] = len(cardiac)
    print(f"  Cardiac surgery, adult, first stay: {len(cardiac):,}")

    # ── ESKD exclusion ───────────────────────────────────────────
    eskd_pids = set()
    if len(pasthx) > 0 and "pasthistorypath" in pasthx.columns:
        eskd_pids = set(
            pasthx[
                pasthx.patientunitstayid.isin(pids)
                & matches_any(pasthx.pasthistorypath, ESKD_PATTERNS)
            ].patientunitstayid
        )
    if len(diag) > 0:
        eskd_pids |= set(
            diag[
                diag.patientunitstayid.isin(pids)
                & matches_any(diag.diagnosisstring, ESKD_PATTERNS)
            ].patientunitstayid
        )
    cardiac = cardiac[~cardiac.patientunitstayid.isin(eskd_pids)]
    pids = set(cardiac.patientunitstayid)
    consort["excl_eskd"] = len(eskd_pids & set(cardiac.patientunitstayid)) + len(
        eskd_pids
    )
    consort["post_eskd"] = len(cardiac)
    print(f"  ESKD excluded: {len(eskd_pids):,} → remaining: {len(cardiac):,}")

    # ── Classify surgery type ────────────────────────────────────
    def classify_surgery(dx_str):
        s = str(dx_str).lower()
        has_cabg = any(p in s for p in ["cabg", "coronary artery bypass"])
        has_valve = any(
            p in s
            for p in [
                "valve",
                "aortic valve",
                "mitral valve",
                "tricuspid",
                "pulmonic valve",
            ]
        )
        if has_cabg and has_valve:
            return "combined"
        if has_cabg:
            return "cabg"
        if has_valve:
            return "valve"
        return "other_cardiac"

    cardiac["surgery_type"] = cardiac.apacheadmissiondx.apply(classify_surgery)
    cardiac["is_female"] = (cardiac.gender.str.lower() == "female").astype(int)
    cardiac["age"] = cardiac["age_num"]

    # ══════════════════════════════════════════════════════════════
    # IV Mg IDENTIFICATION — No 6h window; any postop IV Mg
    # ══════════════════════════════════════════════════════════════
    print("\n── IV Magnesium identification (no 6h window) ──")

    # Filter medication table
    med_elig = (
        med[med.patientunitstayid.isin(pids) & (med.drugordercancelled != "Yes")].copy()
        if len(med) > 0
        else pd.DataFrame()
    )

    # Find IV Mg from medication table
    mg_med = pd.DataFrame()
    if len(med_elig) > 0:
        mg_mask = matches_any(med_elig.drugname, MG_SUPP_PATTERNS)
        # IV route filter (if available)
        if "routeadmin" in med_elig.columns:
            iv_mask = (
                med_elig.routeadmin.str.lower().str.contains(
                    "iv|intravenous|inject", na=False
                )
                | med_elig.routeadmin.isna()
            )
        else:
            iv_mask = pd.Series(True, index=med_elig.index)
        mg_med = med_elig[mg_mask & iv_mask & (med_elig.drugstartoffset >= 0)].copy()
        mg_med = mg_med.rename(columns={"drugstartoffset": "mg_offset_min"})

    # Also check infusionDrug table
    mg_inf = pd.DataFrame()
    if len(inf) > 0 and "drugname" in inf.columns:
        mg_inf = inf[
            inf.patientunitstayid.isin(pids)
            & matches_any(inf.drugname, MG_SUPP_PATTERNS)
            & (inf.infusionoffset >= 0)
        ].copy()
        if len(mg_inf) > 0:
            mg_inf = mg_inf.rename(columns={"infusionoffset": "mg_offset_min"})

    # Combine and find FIRST IV Mg per patient
    mg_all_cols = ["patientunitstayid", "mg_offset_min"]
    frames = []
    if len(mg_med) > 0 and all(c in mg_med.columns for c in mg_all_cols):
        frames.append(mg_med[mg_all_cols])
    if len(mg_inf) > 0 and all(c in mg_inf.columns for c in mg_all_cols):
        frames.append(mg_inf[mg_all_cols])

    if len(frames) == 0:
        print("  ERROR: No IV Mg records found. Check data paths.")
        return None

    mg_all = pd.concat(frames, ignore_index=True)
    first_iv_mg = (
        mg_all.sort_values("mg_offset_min")
        .groupby("patientunitstayid")
        .first()
        .reset_index()
    )
    first_iv_mg["mg_offset_h"] = first_iv_mg.mg_offset_min / 60.0

    treated_pids = set(first_iv_mg.patientunitstayid)
    control_pids = pids - treated_pids

    print(
        f"  IV Mg recipients (any postop): {len(treated_pids):,} "
        f"({100*len(treated_pids)/len(pids):.1f}%)"
    )
    print(f"  No IV Mg (controls):           {len(control_pids):,}")
    describe_col(first_iv_mg.mg_offset_h, "IV Mg timing (h from ICU)")
    consort["treated_any_ivmg"] = len(treated_pids)
    consort["control_no_ivmg"] = len(control_pids)

    # Timing distribution
    for cutoff in [6, 12, 24, 48]:
        n = (first_iv_mg.mg_offset_h <= cutoff).sum()
        print(f"    IV Mg within {cutoff:2d}h: {pct(n, len(first_iv_mg))}")

    # ══════════════════════════════════════════════════════════════
    # Cr EXTRACTION — All measurements for all patients
    # ══════════════════════════════════════════════════════════════
    print("\n── Creatinine extraction ──")
    cr = lab[
        lab.patientunitstayid.isin(pids)
        & lab.labname.str.lower().str.contains("creatinine", na=False)
        & lab.labresult.between(CR_MIN, CR_MAX)
    ].copy()
    cr["cr_offset_h"] = cr.labresultoffset / 60.0
    print(
        f"  Total Cr measurements: {len(cr):,} across {cr.patientunitstayid.nunique():,} patients"
    )

    # ── Cr_pre for TREATED: latest Cr between ICU admit and IV Mg ─
    print("\n── Cr_pre for treated (latest ICU Cr before IV Mg) ──")
    cr_treated = cr[cr.patientunitstayid.isin(treated_pids)].merge(
        first_iv_mg[["patientunitstayid", "mg_offset_min"]], on="patientunitstayid"
    )
    # Cr between ICU admission (offset >= 0) and before IV Mg
    cr_pre_cand = cr_treated[
        (cr_treated.labresultoffset >= 0)
        & (cr_treated.labresultoffset < cr_treated.mg_offset_min)
    ].copy()

    # Latest (closest to IV Mg)
    cr_pre = (
        cr_pre_cand.sort_values("labresultoffset", ascending=False)
        .groupby("patientunitstayid")
        .first()
        .reset_index()
    )
    cr_pre = cr_pre.rename(
        columns={
            "labresult": "cr_pre",
            "labresultoffset": "cr_pre_offset_min",
        }
    )
    cr_pre["cr_pre_offset_h"] = cr_pre.cr_pre_offset_min / 60.0
    cr_pre["gap_to_ivmg_h"] = (cr_pre.mg_offset_min - cr_pre.cr_pre_offset_min) / 60.0

    n_has_cr_pre = len(cr_pre)
    n_no_cr_pre = len(treated_pids) - n_has_cr_pre
    print(f"  Treated with ICU Cr before IV Mg: {pct(n_has_cr_pre, len(treated_pids))}")
    print(
        f"  Treated WITHOUT ICU Cr before IV Mg: {pct(n_no_cr_pre, len(treated_pids))}"
    )
    describe_col(cr_pre.cr_pre, "Cr_pre (mg/dL)")
    describe_col(cr_pre.gap_to_ivmg_h, "Gap Cr_pre → IV Mg (h)")
    consort["treated_has_cr_pre"] = n_has_cr_pre
    consort["treated_no_cr_pre"] = n_no_cr_pre

    # ── Compare with hospitalization Cr (for information) ────────
    print("\n── Comparison: ICU Cr_pre vs hospitalization Cr ──")
    hosp_off = cardiac.set_index("patientunitstayid")["hospitaladmitoffset"].to_dict()
    cr_treated["hosp_off"] = cr_treated.patientunitstayid.map(hosp_off)

    # Hospitalization Cr: closest to hospital admission (±6h), before ICU
    cr_hosp_cand = cr_treated[
        (cr_treated.labresultoffset >= cr_treated.hosp_off - 360)
        & (cr_treated.labresultoffset <= cr_treated.hosp_off + 360)
        & (cr_treated.labresultoffset < 0)  # before ICU
    ].copy()
    cr_hosp_cand["dist_to_admit"] = (
        cr_hosp_cand.labresultoffset - cr_hosp_cand.hosp_off
    ).abs()
    cr_hosp = (
        cr_hosp_cand.sort_values("dist_to_admit")
        .groupby("patientunitstayid")
        .first()
        .reset_index()
    )

    # Merge for comparison
    both = cr_pre[["patientunitstayid", "cr_pre"]].merge(
        cr_hosp[["patientunitstayid", "labresult"]].rename(
            columns={"labresult": "cr_hosp"}
        ),
        on="patientunitstayid",
        how="left",
    )
    n_both = both.cr_hosp.notna().sum()
    if n_both > 0:
        both_valid = both.dropna(subset=["cr_hosp"])
        delta = both_valid.cr_pre - both_valid.cr_hosp
        print(f"  Patients with both ICU Cr_pre and hosp Cr: {n_both}")
        describe_col(both_valid.cr_hosp, "Hospitalization Cr")
        describe_col(both_valid.cr_pre, "ICU Cr_pre (latest before IV Mg)")
        describe_col(delta, "Δ (ICU Cr_pre − hosp Cr)")
        print(f"    Correlation: r = {both_valid.cr_pre.corr(both_valid.cr_hosp):.3f}")
        print(f"    |Δ| > 0.3 mg/dL: {pct((delta.abs() > 0.3).sum(), n_both)}")
    else:
        print(f"  No patients with both measures available for comparison")

    n_hosp_only = len(treated_pids) - n_has_cr_pre
    n_hosp_has = cr_hosp.patientunitstayid.isin(
        treated_pids - set(cr_pre.patientunitstayid)
    ).sum()
    print(f"\n  ◆ Patients lost (no ICU Cr_pre): {n_no_cr_pre}")
    print(f"    Of those, {n_hosp_has} have hospitalization Cr (sensitivity recovery)")

    # ── Cr_post for TREATED: multiple windows ────────────────────
    print("\n── Cr_post for treated (multiple windows after IV Mg) ──")
    cr_post_all = cr_treated[
        cr_treated.labresultoffset > cr_treated.mg_offset_min
    ].copy()
    cr_post_all["post_ivmg_h"] = (
        cr_post_all.labresultoffset - cr_post_all.mg_offset_min
    ) / 60.0

    for wname, (lo_h, hi_h) in CR_POST_WINDOWS.items():
        cand = cr_post_all[
            (cr_post_all.post_ivmg_h >= lo_h) & (cr_post_all.post_ivmg_h <= hi_h)
        ]
        first_post = (
            cand.sort_values("labresultoffset")
            .groupby("patientunitstayid")
            .first()
            .reset_index()
        )
        n = len(first_post)
        primary = " ◀ PRIMARY" if wname == "6_24h" else ""
        print(f"  {wname:>6s}: {pct(n, n_has_cr_pre)} of Cr_pre patients{primary}")
        if n > 0:
            describe_col(first_post.labresult, f"Cr_post ({wname})")

    # ── Baseline Cr exclusion ────────────────────────────────────
    cr_pre_elig = cr_pre[cr_pre.cr_pre < BASELINE_CR_MAX].copy()
    n_excl_cr = n_has_cr_pre - len(cr_pre_elig)
    print(f"\n  Excluded: Cr_pre ≥ {BASELINE_CR_MAX}: {n_excl_cr}")
    consort["excl_cr_high"] = n_excl_cr

    # ── Prevalent AKI check ──────────────────────────────────────
    # Any Cr before IV Mg that already shows ≥1.5x the earliest postop Cr
    print("\n── Prevalent AKI check (pre-IV-Mg AKI) ──")
    # Use the earliest postop Cr as reference baseline
    cr_earliest = (
        cr_pre_cand.sort_values("labresultoffset")
        .groupby("patientunitstayid")
        .first()
        .reset_index()
        .rename(columns={"labresult": "cr_earliest"})
    )
    prev_check = cr_pre_elig.merge(
        cr_earliest[["patientunitstayid", "cr_earliest"]],
        on="patientunitstayid",
        how="left",
    )
    # If cr_pre (latest) is ≥1.5x cr_earliest → AKI already happening
    prev_check["prevalent_aki"] = (
        (prev_check.cr_pre / prev_check.cr_earliest.clip(lower=0.1)) >= 1.5
    ).astype(int)
    n_prevalent = prev_check.prevalent_aki.sum()
    print(
        f"  Prevalent AKI (Cr_pre ≥ 1.5× earliest postop Cr): {pct(n_prevalent, len(prev_check))}"
    )
    consort["excl_prevalent_aki"] = n_prevalent

    # ══════════════════════════════════════════════════════════════
    # BUILD TREATED COHORT
    # ══════════════════════════════════════════════════════════════
    print(f"\n{'─'*50}")
    print("Building treated cohort...")
    treated = cardiac[
        cardiac.patientunitstayid.isin(
            set(cr_pre_elig.patientunitstayid)
            - set(prev_check[prev_check.prevalent_aki == 1].patientunitstayid)
        )
    ].copy()

    # Merge IV Mg timing
    treated = treated.merge(
        first_iv_mg[["patientunitstayid", "mg_offset_min", "mg_offset_h"]],
        on="patientunitstayid",
    )
    # Merge Cr_pre
    treated = treated.merge(
        cr_pre[
            [
                "patientunitstayid",
                "cr_pre",
                "cr_pre_offset_min",
                "cr_pre_offset_h",
                "gap_to_ivmg_h",
            ]
        ],
        on="patientunitstayid",
    )

    # Compute eGFR from Cr_pre
    treated["egfr"] = compute_egfr(treated.cr_pre, treated.age, treated.is_female)

    # Add Cr_post at multiple windows
    for wname, (lo_h, hi_h) in CR_POST_WINDOWS.items():
        cand = cr_post_all[
            cr_post_all.patientunitstayid.isin(set(treated.patientunitstayid))
            & (cr_post_all.post_ivmg_h >= lo_h)
            & (cr_post_all.post_ivmg_h <= hi_h)
        ]
        first_post = (
            cand.sort_values("labresultoffset")
            .groupby("patientunitstayid")
            .first()
            .reset_index()[["patientunitstayid", "labresult", "labresultoffset"]]
            .rename(
                columns={
                    "labresult": f"cr_post_{wname}",
                    "labresultoffset": f"cr_post_offset_{wname}",
                }
            )
        )
        treated = treated.merge(first_post, on="patientunitstayid", how="left")
        # ΔCr
        treated[f"delta_cr_{wname}"] = treated[f"cr_post_{wname}"] - treated.cr_pre

    consort["treated_final"] = len(treated)

    # ── Covariates for treated ───────────────────────────────────
    print(f"\n── Covariates (pre-t0 for treated) ──")
    # Surgery type dummies
    treated["surg_cabg"] = (treated.surgery_type == "cabg").astype(int)
    treated["surg_valve"] = (treated.surgery_type == "valve").astype(int)
    treated["surg_combined"] = (treated.surgery_type == "combined").astype(int)

    # Comorbidities from pastHistory
    for como, patterns in EICU_COMORB_PATTERNS.items():
        if len(pasthx) > 0 and "pasthistorypath" in pasthx.columns:
            como_pids = set(
                pasthx[
                    pasthx.patientunitstayid.isin(set(treated.patientunitstayid))
                    & matches_any(pasthx.pasthistorypath, patterns)
                ].patientunitstayid
            )
        else:
            como_pids = set()
        treated[como] = treated.patientunitstayid.isin(como_pids).astype(int)

    # Perioperative drugs: defined as any before IV Mg (patient-specific t0)
    mg_offsets_dict = dict(zip(treated.patientunitstayid, treated.mg_offset_min))
    trt_pids = set(treated.patientunitstayid)

    for drug_class, patterns in DRUG_CLASSES_EICU.items():
        flagged = set()
        if len(med_elig) > 0:
            for pid_i in trt_pids:
                mg_off = mg_offsets_dict.get(pid_i)
                if mg_off is None:
                    continue
                pt_meds = med_elig[
                    (med_elig.patientunitstayid == pid_i)
                    & (med_elig.drugstartoffset <= mg_off)
                ]
                if len(pt_meds) > 0 and matches_any(pt_meds.drugname, patterns).any():
                    flagged.add(pid_i)
        treated[drug_class] = treated.patientunitstayid.isin(flagged).astype(int)

    # First postop labs BEFORE IV Mg (patient-specific)
    print("  Extracting pre-t0 lab values...")
    for lab_name, col_name, lab_patterns in [
        ("first_mg_value", "first_mg_value", ["magnesium"]),
        ("first_potassium", "first_potassium", ["potassium"]),
        ("first_calcium", "first_calcium", ["calcium"]),
        ("first_heartrate", "first_heartrate", ["heart rate"]),
        ("first_lactate", "first_lactate", ["lactate"]),
    ]:
        if lab_name == "first_heartrate":
            # Heart rate from vitalPeriodic if available, else skip
            # For now, extract from lab table if available
            pass

        lab_sub = lab[
            lab.patientunitstayid.isin(trt_pids)
            & matches_any(lab.labname, lab_patterns)
        ].copy()

        vals = {}
        for pid_i in trt_pids:
            mg_off = mg_offsets_dict.get(pid_i, 0)
            pt_lab = lab_sub[
                (lab_sub.patientunitstayid == pid_i)
                & (lab_sub.labresultoffset >= 0)
                & (lab_sub.labresultoffset < mg_off)
            ]
            if len(pt_lab) > 0:
                # Use first (earliest postop) value
                row = pt_lab.sort_values("labresultoffset").iloc[0]
                vals[pid_i] = row.labresult
        treated[col_name] = treated.patientunitstayid.map(vals)
        n_avail = treated[col_name].notna().sum()
        print(f"    {col_name}: {pct(n_avail, len(treated))}")

    # Lactate missing indicator
    treated["lactate_missing"] = treated.first_lactate.isna().astype(int)

    # Vasopressor before IV Mg (patient-specific)
    if len(med_elig) > 0:
        vaso_patterns = [
            "norepinephrine",
            "vasopressin",
            "epinephrine",
            "phenylephrine",
            "dopamine",
            "dobutamine",
            "milrinone",
        ]
        vaso_pids = set()
        for pid_i in trt_pids:
            mg_off = mg_offsets_dict.get(pid_i, 0)
            pt = med_elig[
                (med_elig.patientunitstayid == pid_i)
                & (med_elig.drugstartoffset >= 0)
                & (med_elig.drugstartoffset <= mg_off)
            ]
            if len(pt) > 0 and matches_any(pt.drugname, vaso_patterns).any():
                vaso_pids.add(pid_i)
        treated["vasopressor_pre_t0"] = treated.patientunitstayid.isin(
            vaso_pids
        ).astype(int)
        print(
            f"    vasopressor_pre_t0: {pct(treated.vasopressor_pre_t0.sum(), len(treated))}"
        )

    # Transfusion before IV Mg (patient-specific)
    transfusion_patterns = [
        "transfusion",
        "blood product",
        "packed red",
        "prbc",
        "red blood cell",
        "ffp",
        "fresh frozen",
        "platelet",
        "cryoprecipitate",
    ]
    if len(treatment) > 0:
        tx_pids = set()
        tx_match = treatment[
            treatment.patientunitstayid.isin(trt_pids)
            & matches_any(treatment.treatmentstring, transfusion_patterns)
        ]
        for pid_i in trt_pids:
            mg_off = mg_offsets_dict.get(pid_i, 0)
            pt_tx = tx_match[
                (tx_match.patientunitstayid == pid_i)
                & (tx_match.treatmentoffset >= 0)
                & (tx_match.treatmentoffset <= mg_off)
            ]
            if len(pt_tx) > 0:
                tx_pids.add(pid_i)
        treated["transfusion_pre_t0"] = treated.patientunitstayid.isin(tx_pids).astype(
            int
        )
        print(
            f"    transfusion_pre_t0: {pct(treated.transfusion_pre_t0.sum(), len(treated))}"
        )

    # BMI
    if "admissionweight" in cardiac.columns and "admissionheight" in cardiac.columns:
        wt = cardiac.set_index("patientunitstayid")["admissionweight"].to_dict()
        ht = cardiac.set_index("patientunitstayid")["admissionheight"].to_dict()
        treated["bmi"] = treated.patientunitstayid.apply(
            lambda x: (
                wt.get(x, np.nan) / ((ht.get(x, np.nan) / 100) ** 2)
                if pd.notna(wt.get(x)) and pd.notna(ht.get(x)) and ht.get(x, 0) > 0
                else np.nan
            )
        )
        treated.loc[~treated.bmi.between(10, 80), "bmi"] = np.nan

    # Outcomes (for secondary analyses)
    treated["icu_mortality"] = (
        (
            cardiac.set_index("patientunitstayid")["unitdischargestatus"]
            .str.lower()
            .eq("expired")
        )
        .reindex(treated.patientunitstayid)
        .fillna(False)
        .astype(int)
        .values
    )
    treated["hosp_mortality"] = (
        (
            cardiac.set_index("patientunitstayid")["hospitaldischargestatus"]
            .str.lower()
            .eq("expired")
        )
        .reindex(treated.patientunitstayid)
        .fillna(False)
        .astype(int)
        .values
    )

    # KDIGO staging (for comparison with DiD)
    treated["aki_kdigo1"] = (
        (treated.cr_post_6_24h / treated.cr_pre >= 1.5)
        | ((treated.cr_post_6_24h - treated.cr_pre) >= 0.3)
    ).astype(int)
    treated.loc[treated.cr_post_6_24h.isna(), "aki_kdigo1"] = np.nan

    # ══════════════════════════════════════════════════════════════
    # BUILD CONTROL COHORT
    # ══════════════════════════════════════════════════════════════
    print(f"\n{'─'*50}")
    print("Building control cohort...")

    control = cardiac[cardiac.patientunitstayid.isin(control_pids)].copy()

    # Require at least 2 Cr measurements postop (for DiD)
    cr_ctrl = cr[cr.patientunitstayid.isin(control_pids) & (cr.labresultoffset >= 0)]
    cr_counts = cr_ctrl.groupby("patientunitstayid").size()
    ctrl_with_2cr = set(cr_counts[cr_counts >= 2].index)
    control = control[control.patientunitstayid.isin(ctrl_with_2cr)].copy()
    n_no_2cr = len(control_pids) - len(ctrl_with_2cr)
    print(f"  Controls with ≥2 postop Cr: {pct(len(control), len(control_pids))}")
    print(f"  Excluded (<2 Cr): {n_no_2cr}")
    consort["control_has_2cr"] = len(control)

    # Baseline Cr exclusion
    ctrl_first_cr = (
        cr_ctrl[cr_ctrl.patientunitstayid.isin(set(control.patientunitstayid))]
        .sort_values("labresultoffset")
        .groupby("patientunitstayid")
        .first()
        .reset_index()
    )
    ctrl_excl_cr = set(
        ctrl_first_cr[ctrl_first_cr.labresult >= BASELINE_CR_MAX].patientunitstayid
    )
    control = control[~control.patientunitstayid.isin(ctrl_excl_cr)].copy()
    print(f"  Excluded (first Cr ≥ {BASELINE_CR_MAX}): {len(ctrl_excl_cr)}")

    # ESKD already excluded above

    # eGFR for controls (from first postop Cr)
    ctrl_first = ctrl_first_cr[
        ctrl_first_cr.patientunitstayid.isin(set(control.patientunitstayid))
    ].rename(
        columns={
            "labresult": "first_postop_cr",
            "labresultoffset": "first_cr_offset_min",
        }
    )
    control = control.merge(
        ctrl_first[["patientunitstayid", "first_postop_cr", "first_cr_offset_min"]],
        on="patientunitstayid",
        how="left",
    )
    control["egfr"] = compute_egfr(
        control.first_postop_cr, control.age_num, control.is_female
    )
    control["is_female"] = (control.gender.str.lower() == "female").astype(int)
    control["age"] = control["age_num"]

    # Surgery type
    control["surgery_type"] = control.apacheadmissiondx.apply(classify_surgery)
    control["surg_cabg"] = (control.surgery_type == "cabg").astype(int)
    control["surg_valve"] = (control.surgery_type == "valve").astype(int)
    control["surg_combined"] = (control.surgery_type == "combined").astype(int)

    # Comorbidities
    for como, patterns in EICU_COMORB_PATTERNS.items():
        if len(pasthx) > 0 and "pasthistorypath" in pasthx.columns:
            como_pids = set(
                pasthx[
                    pasthx.patientunitstayid.isin(set(control.patientunitstayid))
                    & matches_any(pasthx.pasthistorypath, patterns)
                ].patientunitstayid
            )
        else:
            como_pids = set()
        control[como] = control.patientunitstayid.isin(como_pids).astype(int)

    # Drugs for controls: use first 6h of ICU (common anchor)
    ctrl_pids_set = set(control.patientunitstayid)
    for drug_class, patterns in DRUG_CLASSES_EICU.items():
        flagged = set()
        if len(med_elig) > 0:
            early = med_elig[
                med_elig.patientunitstayid.isin(ctrl_pids_set)
                & (med_elig.drugstartoffset >= 0)
                & (med_elig.drugstartoffset <= 360)  # first 6h
            ]
            if len(early) > 0 and matches_any(early.drugname, patterns).any():
                flagged = set(
                    early[matches_any(early.drugname, patterns)].patientunitstayid
                )
        control[drug_class] = control.patientunitstayid.isin(flagged).astype(int)

    # First postop labs for controls (first 6h, same anchor)
    for lab_name, col_name, lab_patterns in [
        ("first_mg_value", "first_mg_value", ["magnesium"]),
        ("first_potassium", "first_potassium", ["potassium"]),
        ("first_calcium", "first_calcium", ["calcium"]),
        ("first_lactate", "first_lactate", ["lactate"]),
    ]:
        lab_sub = lab[
            lab.patientunitstayid.isin(ctrl_pids_set)
            & matches_any(lab.labname, lab_patterns)
            & (lab.labresultoffset >= 0)
            & (lab.labresultoffset <= 360)
        ]
        first_val = (
            lab_sub.sort_values("labresultoffset")
            .groupby("patientunitstayid")
            .first()
            .reset_index()[["patientunitstayid", "labresult"]]
            .rename(columns={"labresult": col_name})
        )
        control = control.merge(first_val, on="patientunitstayid", how="left")

    control["lactate_missing"] = control.first_lactate.isna().astype(int)

    # Vasopressor / transfusion for controls (first 6h)
    if len(med_elig) > 0:
        vaso_early = med_elig[
            med_elig.patientunitstayid.isin(ctrl_pids_set)
            & (med_elig.drugstartoffset >= 0)
            & (med_elig.drugstartoffset <= 360)
        ]
        vaso_patterns = [
            "norepinephrine",
            "vasopressin",
            "epinephrine",
            "phenylephrine",
            "dopamine",
            "dobutamine",
            "milrinone",
        ]
        vp = set(
            vaso_early[
                matches_any(vaso_early.drugname, vaso_patterns)
            ].patientunitstayid
        )
        control["vasopressor_pre_t0"] = control.patientunitstayid.isin(vp).astype(int)

    if len(treatment) > 0:
        tx_early = treatment[
            treatment.patientunitstayid.isin(ctrl_pids_set)
            & matches_any(treatment.treatmentstring, transfusion_patterns)
            & (treatment.treatmentoffset >= 0)
            & (treatment.treatmentoffset <= 360)
        ]
        tp = set(tx_early.patientunitstayid)
        control["transfusion_pre_t0"] = control.patientunitstayid.isin(tp).astype(int)

    # BMI for controls
    if "admissionweight" in cardiac.columns:
        control["bmi"] = control.patientunitstayid.apply(
            lambda x: (
                wt.get(x, np.nan) / ((ht.get(x, np.nan) / 100) ** 2)
                if pd.notna(wt.get(x)) and pd.notna(ht.get(x)) and ht.get(x, 0) > 0
                else np.nan
            )
        )
        control.loc[~control.bmi.between(10, 80), "bmi"] = np.nan

    consort["control_final"] = len(control)

    # ══════════════════════════════════════════════════════════════
    # ALL Cr MEASUREMENTS (for temporal matching in R)
    # ══════════════════════════════════════════════════════════════
    print(f"\n── Exporting all Cr measurements for temporal matching ──")
    all_cohort_pids = set(treated.patientunitstayid) | set(control.patientunitstayid)
    cr_all = cr[cr.patientunitstayid.isin(all_cohort_pids) & (cr.labresultoffset >= 0)][
        ["patientunitstayid", "labresult", "labresultoffset"]
    ].copy()
    cr_all["cr_offset_h"] = cr_all.labresultoffset / 60.0
    cr_all.to_csv(os.path.join(RESULTS, "20_did_cr_all_eicu.csv"), index=False)
    print(
        f"  Exported {len(cr_all):,} Cr measurements for {cr_all.patientunitstayid.nunique()} patients"
    )

    # ══════════════════════════════════════════════════════════════
    # COVARIATE TIMING REPORT
    # ══════════════════════════════════════════════════════════════
    print(f"\n{'─'*50}")
    print("COVARIATE TIMING ASSESSMENT")
    print("─" * 50)
    print("  For TREATED: covariates measured before patient-specific t0 (IV Mg)")
    print("  For CONTROLS: covariates measured in first 6h of ICU (common anchor)")
    print()

    # Compare vasopressor/transfusion availability
    for var_name in ["vasopressor_pre_t0", "transfusion_pre_t0"]:
        if var_name in treated.columns and var_name in control.columns:
            t_rate = treated[var_name].mean()
            c_rate = control[var_name].mean()
            smd = abs(t_rate - c_rate) / np.sqrt(
                (t_rate * (1 - t_rate) + c_rate * (1 - c_rate)) / 2 + 1e-10
            )
            print(f"  {var_name}:")
            print(f"    Treated:  {100*t_rate:.1f}%")
            print(f"    Control:  {100*c_rate:.1f}%")
            print(f"    Raw SMD:  {smd:.3f}")
            if smd > 0.25:
                print(f"    ⚠ High imbalance — consider dropping if PS can't fix")
            print()

    # ══════════════════════════════════════════════════════════════
    # SAVE
    # ══════════════════════════════════════════════════════════════
    treated["treated"] = 1
    control["treated"] = 0

    treated.to_csv(os.path.join(RESULTS, "20_did_treated_eicu.csv"), index=False)
    control.to_csv(os.path.join(RESULTS, "20_did_control_eicu.csv"), index=False)
    print(
        f"\n  Saved: 20_did_treated_eicu.csv ({len(treated)} patients, {treated.shape[1]} cols)"
    )
    print(
        f"  Saved: 20_did_control_eicu.csv ({len(control)} patients, {control.shape[1]} cols)"
    )

    # ══════════════════════════════════════════════════════════════
    # CONSORT SUMMARY
    # ══════════════════════════════════════════════════════════════
    print(f"\n{'='*70}")
    print("eICU CONSORT SUMMARY")
    print("=" * 70)
    steps = [
        ("Total ICU admissions", consort.get("total_icu", 0)),
        ("Cardiac surgery, adult, 1st stay", consort.get("cardiac_adult_first", 0)),
        ("Post-ESKD exclusion", consort.get("post_eskd", 0)),
        ("Received any postop IV Mg", consort.get("treated_any_ivmg", 0)),
        ("  Has ICU Cr before IV Mg", consort.get("treated_has_cr_pre", 0)),
        ("  No ICU Cr before IV Mg (LOST)", consort.get("treated_no_cr_pre", 0)),
        ("  Excl: Cr_pre ≥ 4.0", consort.get("excl_cr_high", 0)),
        ("  Excl: prevalent AKI", consort.get("excl_prevalent_aki", 0)),
        ("  TREATED FINAL", consort.get("treated_final", 0)),
        ("Never received IV Mg", consort.get("control_no_ivmg", 0)),
        ("  Has ≥2 postop Cr", consort.get("control_has_2cr", 0)),
        ("  CONTROL FINAL", consort.get("control_final", 0)),
    ]
    for label, n in steps:
        print(f"  {label:<40s} {n:>8,}")

    # Cr_post availability by window
    print(f"\n  Cr_post availability (treated, within final cohort):")
    for wname in CR_POST_WINDOWS:
        col = f"cr_post_{wname}"
        if col in treated.columns:
            n = treated[col].notna().sum()
            print(f"    {wname:>6s}: {pct(n, len(treated))}")

    # DiD-ready sample (Cr_pre + Cr_post both available)
    for wname in ["6_24h", "6_48h", "0_24h"]:
        col = f"cr_post_{wname}"
        if col in treated.columns:
            n_ready = treated[col].notna().sum()
            print(f"\n  DiD-ready (Cr_pre + {wname}): {pct(n_ready, len(treated))}")
            if n_ready > 0:
                sub = treated[treated[col].notna()]
                describe_col(sub[f"delta_cr_{wname}"], f"ΔCr ({wname})")

    # Subgroup sizes for power
    print(f"\n  Surgery type distribution (treated):")
    for stype in ["cabg", "valve", "combined", "other_cardiac"]:
        n = (treated.surgery_type == stype).sum()
        print(f"    {stype:<15s}: {pct(n, len(treated))}")

    print(f"\n  First serum Mg strata (treated):")
    if "first_mg_value" in treated.columns:
        for lo, hi, label in [
            (0, 1.8, "<1.8"),
            (1.8, 2.0, "1.8-2.0"),
            (2.0, 2.3, "2.0-2.3"),
            (2.3, 99, ">2.3"),
        ]:
            n = treated.first_mg_value.between(lo, hi, inclusive="left").sum()
            print(f"    {label:<10s}: {pct(n, len(treated))}")

    return consort


# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("=" * 70)
    print("20_did_etl.py — DiD Cohort Construction for Mg → AKI")
    print("  Design: patient-specific t0, continuous Cr outcome")
    print("=" * 70)

    args = [a.lower() for a in sys.argv[1:]]
    run_all = len(args) == 0

    consort_rows = []

    if run_all or "eicu" in args:
        c = run_eicu()
        if c:
            c["db"] = "eICU"
            consort_rows.append(c)

    # MIMIC ETL would follow the same structure but using
    # inputevents for IV Mg timing and labevents for Cr.
    # Skeleton included in comments; expand after eICU validation.

    if consort_rows:
        consort_df = pd.DataFrame(consort_rows)
        consort_df.to_csv(os.path.join(RESULTS, "20_did_consort.csv"), index=False)
        print(f"\n  Saved: 20_did_consort.csv")

    print(f"\n{'='*70}")
    print("NEXT STEPS")
    print("=" * 70)
    print("  1. Review CONSORT numbers above — check attrition at each step")
    print("  2. Decide on primary Cr_post window (recommend 6-24h)")
    print("  3. Assess vasopressor/transfusion SMDs — drop if >0.25 raw")
    print("  4. Run 21_did_matching.R for 1:4 PSM + temporal alignment")
    print("  5. Run 22_did_analysis.R for DiD estimation")
    print(f"{'='*70}")
