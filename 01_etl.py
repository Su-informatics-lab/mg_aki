#!/usr/bin/env python3
"""
01_etl.py — Cohort construction for eICU-CRD and MIMIC-IV
  ◆ v2: Three critical design fixes (gtgs/ziyue review)
    1. Baseline Cr = pre-operative (admission), not first ICU Cr
    2. AKI outcome Cr only counted after 6 h landmark (washout)
    3. 6 h landmark exclusion in eligibility (not sensitivity)

  Section A: eICU-CRD → results/01_analysis_a_cohort.csv
  Section B: MIMIC-IV  → results/04_mimic_cohort.csv

Run:  python 01_etl.py            # both databases
      python 01_etl.py eicu       # eICU only
      python 01_etl.py mimic      # MIMIC only
"""

import os
import sys
import warnings

warnings.filterwarnings("ignore")
from importlib.util import module_from_spec, spec_from_file_location

import numpy as np
import pandas as pd

# ── Load 00_config.py ─────────────────────────────────────────────────
_cfg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "00_config.py")
_spec = spec_from_file_location("config", _cfg_path)
cfg = module_from_spec(_spec)
_spec.loader.exec_module(cfg)

# ◆ CHANGED: 6 h landmark is now a first-class constant
LANDMARK_MIN = 360  # 6 hours in minutes — time zero for all analyses
LANDMARK_HOURS = 6


# =====================================================================
# SHARED HELPERS (unchanged)
# =====================================================================
def save(df, name):
    path = os.path.join(cfg.RESULTS, name)
    df.to_csv(path, index=False)
    print(f"  → {path}  ({len(df):,} × {df.shape[1]})")


def load_eicu_csv(table_name, usecols=None):
    for pattern in [
        f"{table_name}.csv",
        f"{table_name}.csv.gz",
        f"{table_name.lower()}.csv",
        f"{table_name.lower()}.csv.gz",
    ]:
        path = os.path.join(cfg.DATA_ROOT, pattern)
        if os.path.exists(path):
            uc = [c.lower() for c in usecols] if usecols else None
            df = pd.read_csv(path, low_memory=False, usecols=uc)
            df.columns = df.columns.str.lower()
            print(f"  {table_name}: {len(df):,} rows")
            return df
    print(f"  WARNING: {table_name} not found in {cfg.DATA_ROOT}")
    return pd.DataFrame()


def gz(path):
    return path if os.path.exists(path) else path.replace(".csv.gz", ".csv")


def age_numeric(age_str):
    if pd.isna(age_str):
        return np.nan
    s = str(age_str).strip()
    if s.startswith(">"):
        return cfg.AGE_CAP
    try:
        return float(s)
    except ValueError:
        return np.nan


def matches_any(series, patterns):
    return series.str.lower().str.contains("|".join(patterns), na=False)


def drug_match(series, patterns):
    s = series.str.lower().fillna("")
    mask = pd.Series(False, index=series.index)
    for p in patterns:
        mask |= s.str.contains(p, na=False)
    return mask


def compute_egfr(cr, age, is_female):
    """CKD-EPI 2021 race-free (vectorized)."""
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


def compute_egfr_scalar(cr, age, is_female):
    if pd.isna(cr) or pd.isna(age) or cr <= 0:
        return np.nan
    k = 0.7 if is_female else 0.9
    a = -0.241 if is_female else -0.302
    r = cr / k
    val = 142 * (min(r, 1) ** a) * (max(r, 1) ** (-1.200)) * (0.9938**age)
    return val * 1.012 if is_female else val


# =====================================================================
# MIMIC-IV CONSTANTS (unchanged)
# =====================================================================
MIMIC_ROOT = os.path.expanduser("~/mg_aki/mimic-iv-3.1")
MIMIC_HOSP = os.path.join(MIMIC_ROOT, "hosp")
MIMIC_ICU = os.path.join(MIMIC_ROOT, "icu")

LAB_MG = [50960]
LAB_CR = [50912, 52546]
LAB_K = [50971]
LAB_CA = [50893]
LAB_LAC = [50813]
VITAL_HR = [220045]
VITAL_WEIGHT = [226512]
VITAL_HEIGHT = [226730]
VASO_ITEMS = [221906, 221289, 222315, 221749, 221662, 221653, 221986]
MG_SUPP_ITEMS = [222011, 227523]
K_SUPP_ITEMS = [225166, 225168, 222139, 227521, 227522]
METO_ITEMS = [225974]
AMIO_ITEMS = [221347, 228339, 229654, 230034]

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

COMORB_ICD = {
    "hx_chf": {
        9: ["428", "4254", "4255", "4257", "4258", "4259"],
        10: ["I50", "I110", "I130", "I132"],
    },
    "hx_hypertension": {
        9: ["401", "402", "403", "404", "405"],
        10: ["I10", "I11", "I12", "I13", "I15"],
    },
    "hx_diabetes": {9: ["250"], 10: ["E10", "E11", "E13"]},
    "hx_ckd": {9: ["585", "586"], 10: ["N18", "N19"]},
    "hx_copd": {9: ["491", "492", "496"], 10: ["J43", "J44"]},
    "hx_pvd": {9: ["4431", "4432", "4438", "4439"], 10: ["I73", "I771", "I790"]},
    "hx_stroke": {9: ["430", "431", "434", "436"], 10: ["I60", "I61", "I63", "I64"]},
    "hx_liver": {
        9: ["5712", "5714", "5715", "5716", "5718", "5719"],
        10: ["K70", "K713", "K714", "K715", "K721", "K729", "K73", "K74"],
    },
    "hx_afib": {9: ["42731"], 10: ["I48"]},
}
NEPHROTOX_MIMIC = {
    "nephrotox_loop_diuretic": ["furosemide", "bumetanide", "torsemide"],
    "nephrotox_nsaid": [
        "ibuprofen",
        "ketorolac",
        "naproxen",
        "indomethacin",
        "diclofenac",
        "meloxicam",
        "celecoxib",
    ],
    "nephrotox_acei_arb": [
        "lisinopril",
        "enalapril",
        "captopril",
        "ramipril",
        "benazepril",
        "losartan",
        "valsartan",
        "irbesartan",
        "olmesartan",
        "candesartan",
    ],
    "nephrotox_ppi": [
        "pantoprazole",
        "omeprazole",
        "esomeprazole",
        "lansoprazole",
        "rabeprazole",
    ],
}
BB_DRUGS = [
    "metoprolol",
    "atenolol",
    "propranolol",
    "carvedilol",
    "bisoprolol",
    "labetalol",
    "nadolol",
    "sotalol",
    "esmolol",
    "nebivolol",
]
STEROID_DRUGS = [
    "methylprednisolone",
    "dexamethasone",
    "hydrocortisone",
    "prednisone",
    "prednisolone",
    "solumedrol",
]
ANTIARR_DRUGS = [
    "amiodarone",
    "sotalol",
    "flecainide",
    "propafenone",
    "dofetilide",
    "dronedarone",
    "digoxin",
]
AF_ICD9 = ["42731"]
AF_ICD10 = ["I48"]
NC_ICD = {
    "nc_fracture": {
        9: [str(i) for i in range(800, 830)],
        10: ["S12", "S22", "S32", "S42", "S52", "S62", "S72", "S82", "S92"],
    },
    "nc_uti": {9: ["5990"], 10: ["N390"]},
}
NEURO_ICD = {
    "neuro_delirium": {9: ["2930", "2931"], 10: ["F05"]},
    "neuro_seizure": {9: ["345", "780.3"], 10: ["G40", "R56"]},
    "neuro_stroke_postop": {9: ["434", "436"], 10: ["I63", "I64"]},
    "neuro_encephalopathy": {9: ["3481", "3489"], 10: ["G93"]},
}
VT_ICD = {9: ["4271", "42741", "42742"], 10: ["I472", "I490", "I4901", "I4902"]}


def matches_icd(dx_df, hadm_ids, code_map):
    sub = dx_df[dx_df.hadm_id.isin(hadm_ids)]
    pids = set()
    for ver, prefixes in code_map.items():
        v = sub[sub.icd_version == ver]
        for p in prefixes:
            pids |= set(v[v.icd_code.str.startswith(p)].hadm_id)
    return pids


# =====================================================================
# SECTION A: eICU-CRD ETL
# =====================================================================
def run_eicu():
    print("=" * 70)
    print("SECTION A: eICU-CRD Cohort Construction (v2: pre-op Cr + 6h landmark)")
    print("=" * 70)
    print(f"  DATA_ROOT: {cfg.DATA_ROOT}")

    # ── Load tables ──────────────────────────────────────────────
    t = {}
    t["patient"] = load_eicu_csv("patient")
    t["lab"] = load_eicu_csv("lab")
    t["medication"] = load_eicu_csv(
        "medication",
        usecols=[
            "patientunitstayid",
            "medicationid",
            "drugstartoffset",
            "drugstopoffset",
            "drugname",
            "drugordercancelled",
            "routeadmin",
            "dosage",
        ],
    )
    t["infusionDrug"] = load_eicu_csv("infusionDrug")
    t["diagnosis"] = load_eicu_csv("diagnosis")
    t["pastHistory"] = load_eicu_csv("pastHistory")
    t["treatment"] = load_eicu_csv("treatment")
    t["admissionDrug"] = load_eicu_csv("admissionDrug")
    t["hospital"] = load_eicu_csv("hospital")
    t["apachePatientResult"] = load_eicu_csv("apachePatientResult")
    t["apacheApsVar"] = load_eicu_csv("apacheApsVar")
    t["apachePredVar"] = load_eicu_csv("apachePredVar")
    t["intakeOutput"] = load_eicu_csv("intakeOutput")

    # ── Cardiac surgery cohort (unchanged) ───────────────────────
    consort = {}
    pt = t["patient"].copy()
    consort["total_icu_stays"] = len(pt)
    pt["age_num"] = pt["age"].apply(age_numeric)
    adults = pt[pt.age_num >= cfg.MIN_AGE].copy()
    consort["adults"] = len(adults)

    by_dx = matches_any(adults.apacheadmissiondx, cfg.CARDIAC_DX_PATTERNS)
    by_unit = adults.unittype.isin(cfg.CARDIAC_UNIT_TYPES)
    cardiac = adults[by_dx | by_unit].copy()
    consort["cardiac_surgery"] = len(cardiac)
    print(f"  Cardiac surgery: {len(cardiac):,}")

    cardiac = (
        cardiac.sort_values("hospitaladmitoffset", ascending=False)
        .groupby("uniquepid")
        .first()
        .reset_index()
    )
    consort["first_stay"] = len(cardiac)

    def _surgery_type(dx):
        if pd.isna(dx):
            return "unknown"
        dx = str(dx).lower()
        has_cabg = any(k in dx for k in cfg.SURGERY_TYPE_MAP["cabg"])
        has_valve = any(k in dx for k in cfg.SURGERY_TYPE_MAP["valve"])
        if has_cabg and has_valve:
            return "combined"
        if has_cabg:
            return "cabg"
        if has_valve:
            return "valve"
        return "other_cardiac"

    cardiac["surgery_type"] = cardiac.apacheadmissiondx.apply(_surgery_type)

    # ── Mg exposure (unchanged) ──────────────────────────────────
    pids = set(cardiac.patientunitstayid)
    lab = t["lab"]
    mg = lab[
        lab.patientunitstayid.isin(pids)
        & matches_any(lab.labname, cfg.MG_LABNAMES)
        & lab.labresult.between(cfg.MG_PLAUSIBLE_MIN, cfg.MG_PLAUSIBLE_MAX)
    ]
    mg_window = mg[
        (mg.labresultoffset >= 0) & (mg.labresultoffset <= cfg.MG_WINDOW_MIN)
    ]
    first_mg = (
        mg_window.sort_values("labresultoffset")
        .groupby("patientunitstayid")
        .first()
        .reset_index()[["patientunitstayid", "labresult", "labresultoffset"]]
        .rename(columns={"labresult": "first_mg_value", "labresultoffset": "mg_offset"})
    )
    first_mg["mg_hours"] = first_mg.mg_offset / 60.0
    first_mg["mg_quartile"] = pd.qcut(
        first_mg.first_mg_value, 4, labels=["Q1", "Q2", "Q3", "Q4"]
    )
    first_mg["mg_category"] = pd.cut(
        first_mg.first_mg_value,
        bins=[0, 1.8, 2.3, 999],
        labels=["hypo", "normal", "hyper"],
    )
    cohort = cardiac.merge(first_mg, on="patientunitstayid")
    consort["has_mg"] = len(cohort)
    print(f"  Mg within {cfg.MG_WINDOW_HOURS}h: {len(cohort):,}")

    # ══════════════════════════════════════════════════════════════
    # ◆ CHANGED: Baseline creatinine = PRE-OPERATIVE (admission)
    #   Priority 1: Cr near hospital admission (within ±6 h of hosp admit)
    #   Priority 2: Any pre-ICU Cr (offset < 0), closest to hosp admit
    #   Priority 3: First ICU Cr (offset 0–60 min) — flagged, not ideal
    #   Primary: most recent pre-operative Cr (Yan/ADQI standard)
    #   Also saves admission Cr as sensitivity column
    #   Reference: Nadim 2018 ADQI, KDIGO CSA-AKI consensus
    # ══════════════════════════════════════════════════════════════
    cr = lab[
        lab.patientunitstayid.isin(pids)
        & matches_any(lab.labname, cfg.CR_LABNAMES)
        & lab.labresult.between(cfg.CR_PLAUSIBLE_MIN, cfg.CR_PLAUSIBLE_MAX)
    ].copy()

    hosp_off = cardiac.set_index("patientunitstayid")["hospitaladmitoffset"].to_dict()
    cr_cohort = cr[cr.patientunitstayid.isin(set(cohort.patientunitstayid))].copy()
    cr_cohort["hosp_off"] = cr_cohort.patientunitstayid.map(hosp_off)

    # ── PRIMARY: most recent pre-ICU Cr (术前最后一个) ────────────
    pre_icu = cr_cohort[cr_cohort.labresultoffset < 0]
    bl_primary = (
        pre_icu.sort_values("labresultoffset", ascending=False)
        .groupby("patientunitstayid")
        .first()
        .reset_index()[["patientunitstayid", "labresult", "labresultoffset"]]
        .rename(
            columns={
                "labresult": "baseline_cr",
                "labresultoffset": "baseline_cr_offset",
            }
        )
    )
    bl_primary["baseline_source"] = "last_preop"
    n_primary = len(bl_primary)

    cohort = cohort.merge(bl_primary, on="patientunitstayid")
    consort["has_baseline_cr"] = len(cohort)

    # ── SENSITIVITY: admission Cr (closest to hospital admission) ─
    admit_window = cr_cohort[
        cr_cohort.patientunitstayid.isin(set(cohort.patientunitstayid))
        & (cr_cohort.labresultoffset >= cr_cohort.hosp_off - LANDMARK_MIN)
        & (cr_cohort.labresultoffset <= cr_cohort.hosp_off + LANDMARK_MIN)
    ].copy()
    admit_window["dist_to_admit"] = abs(
        admit_window.labresultoffset - admit_window.hosp_off
    )
    bl_admit_sens = (
        admit_window.sort_values("dist_to_admit")
        .groupby("patientunitstayid")
        .first()
        .reset_index()[["patientunitstayid", "labresult"]]
        .rename(columns={"labresult": "baseline_cr_admit"})
    )
    cohort = cohort.merge(bl_admit_sens, on="patientunitstayid", how="left")
    n_admit = bl_admit_sens.patientunitstayid.nunique()

    print(f"  ◆ Primary baseline: last pre-op Cr (n={n_primary})")
    print(
        f"  ◆ Sensitivity: admission Cr available for {n_admit} "
        f"({100*n_admit/len(cohort):.0f}%)"
    )
    for src, n in cohort.baseline_source.value_counts().items():
        print(f"    {src}: {n} ({100*n/len(cohort):.1f}%)")

    # Eadon sensitivity baseline (unchanged)
    eadon = cr[
        (cr.labresultoffset >= cfg.BASELINE_EADON_WINDOW_MIN)
        & (cr.labresultoffset <= cfg.BASELINE_EADON_WINDOW_MAX)
    ]
    bl_eadon = (
        eadon.sort_values("labresult")
        .groupby("patientunitstayid")
        .first()
        .reset_index()[["patientunitstayid", "labresult"]]
        .rename(columns={"labresult": "baseline_cr_eadon"})
    )
    cohort = cohort.merge(bl_eadon, on="patientunitstayid", how="left")

    # Exclude high Cr + ESKD (unchanged)
    cohort = cohort[cohort.baseline_cr < cfg.BASELINE_CR_MAX].copy()
    ph = t["pastHistory"]
    eskd_pids = set()
    if len(ph) > 0:
        eskd_pids = set(
            ph[
                ph.patientunitstayid.isin(cohort.patientunitstayid)
                & (
                    matches_any(ph.pasthistorypath, cfg.ESKD_DX_PATTERNS)
                    | matches_any(ph.pasthistoryvalue, cfg.ESKD_DX_PATTERNS)
                )
            ].patientunitstayid
        )
    cohort = cohort[~cohort.patientunitstayid.isin(eskd_pids)].copy()

    apv = t["apachePredVar"]
    if len(apv) > 0:
        apv_elig = apv[apv.patientunitstayid.isin(set(cohort.patientunitstayid))]
        if "dialysis" in apv_elig.columns:
            dial_pids = set(apv_elig[apv_elig.dialysis == 1].patientunitstayid)
            cohort = cohort[~cohort.patientunitstayid.isin(dial_pids)].copy()

    consort["eligible_pre_landmark"] = len(cohort)

    # ══════════════════════════════════════════════════════════════
    # ◆ CHANGED: 6 h landmark exclusion (eligibility, not sensitivity)
    #   Exclude patients discharged, dead, or with AKI before 6 h.
    #   This is now part of eligibility, not a sensitivity analysis.
    # ══════════════════════════════════════════════════════════════
    n_before = len(cohort)

    # Discharged before landmark
    if "unitdischargeoffset" in cohort.columns:
        discharged_pre = cohort.unitdischargeoffset <= LANDMARK_MIN
        cohort = cohort[~discharged_pre].copy()

    # Died before landmark
    cohort["death_offset_min_raw"] = np.where(
        cohort.hospitaldischargestatus.str.lower() == "expired",
        cohort.hospitaldischargeoffset,
        np.nan,
    )
    died_pre = cohort.death_offset_min_raw.notna() & (
        cohort.death_offset_min_raw <= LANDMARK_MIN
    )
    cohort = cohort[~died_pre].copy()

    consort["landmark_excluded"] = n_before - len(cohort)
    consort["eligible_post_landmark"] = len(cohort)
    print(
        f"  ◆ 6h landmark: excluded {n_before - len(cohort)} "
        f"(discharge/death before 6h), {len(cohort)} remain"
    )

    # ══════════════════════════════════════════════════════════════
    # ◆ CHANGED: AKI phenotype — only Cr AFTER 6 h landmark
    #   Washout: Cr in first 6 h may reflect CPB hemodilution or
    #   surgical muscle damage, not true renal injury.
    #   Reference: BMC Anesthesiol 2021 (pCr at H6 post-CPB)
    # ══════════════════════════════════════════════════════════════
    cr_all = cr
    results = []
    for _, row in cohort.iterrows():
        pid = row.patientunitstayid
        bl_cr = row.baseline_cr
        pt_cr = cr_all[cr_all.patientunitstayid == pid]

        # ◆ Follow-up: only Cr AFTER 6h landmark, up to 7 days
        fu = pt_cr[
            (pt_cr.labresultoffset > LANDMARK_MIN)
            & (pt_cr.labresultoffset <= cfg.AKI_WINDOW_7D_MIN)
        ].sort_values("labresultoffset")

        # ◆ Prevalent AKI: Cr between ICU admission and 6h that
        #   already meets KDIGO vs pre-op baseline
        pre_lm = pt_cr[
            (pt_cr.labresultoffset >= 0) & (pt_cr.labresultoffset <= LANDMARK_MIN)
        ]
        prevalent = int(
            len(pre_lm) > 0
            and bl_cr > 0
            and (pre_lm.labresult / bl_cr).max() >= cfg.AKI_RATIO_STAGE1
        )

        aki_15x = 0
        aki_time = np.nan
        aki_delta03 = 0
        if len(fu) > 0 and bl_cr > 0:
            # Ratio criterion: any fu Cr >= 1.5× baseline
            hits = fu[fu.labresult / bl_cr >= cfg.AKI_RATIO_STAGE1]
            if len(hits) > 0:
                aki_15x = 1
                aki_time = hits.labresultoffset.iloc[0]

            # Delta criterion: any fu Cr >= baseline + 0.3 within 48h
            fu_48 = fu[fu.labresultoffset <= LANDMARK_MIN + cfg.AKI_WINDOW_48H_MIN]
            if len(fu_48) > 0 and (fu_48.labresult - bl_cr).max() >= cfg.AKI_DELTA_48H:
                aki_delta03 = 1

        max_cr = fu.labresult.max() if len(fu) > 0 else np.nan
        max_ratio = max_cr / bl_cr if bl_cr > 0 and not np.isnan(max_cr) else 0

        results.append(
            {
                "patientunitstayid": pid,
                "aki_primary": aki_15x,
                "aki_time_offset": aki_time,
                "aki_delta03": aki_delta03,
                "aki_stage2": int(max_ratio >= cfg.AKI_RATIO_STAGE2),
                "aki_stage3": int(
                    max_ratio >= cfg.AKI_RATIO_STAGE3
                    or (
                        max_cr >= cfg.AKI_CR_ABSOLUTE if not np.isnan(max_cr) else False
                    )
                ),
                "max_followup_cr": max_cr,
                "max_cr_ratio": max_ratio,
                "n_followup_cr": len(fu),
                "prevalent_aki": prevalent,
            }
        )

    aki_df = pd.DataFrame(results)
    cohort = cohort.merge(aki_df, on="patientunitstayid")

    # ◆ CHANGED: exclude prevalent AKI (vs pre-op baseline, within 0–6h)
    n_prevalent = cohort.prevalent_aki.sum()
    cohort = cohort[cohort.prevalent_aki == 0].copy()
    consort["excluded_prevalent_aki"] = n_prevalent

    cohort["aki_kdigo1"] = (
        (cohort.aki_primary == 1) | (cohort.aki_delta03 == 1)
    ).astype(int)
    cohort["time_to_aki_hours"] = cohort.aki_time_offset / 60.0

    if "unitdischargeoffset" in cohort.columns:
        cohort["censor_offset"] = cohort.unitdischargeoffset.clip(
            upper=cfg.AKI_WINDOW_7D_MIN
        )
        cohort["time_to_event_hours"] = np.where(
            cohort.aki_primary == 1,
            cohort.time_to_aki_hours,
            (cohort.censor_offset - LANDMARK_MIN) / 60.0,
        )  # ◆ from landmark

    consort["eligible_final"] = len(cohort)
    print(f"  ◆ Prevalent AKI (pre-op→6h): {n_prevalent} excluded")
    print(
        f"  AKI KDIGO≥1: {cohort.aki_kdigo1.sum()} "
        f"({cohort.aki_kdigo1.mean()*100:.1f}%)"
    )

    # ── Covariates (unchanged except BMI cap) ────────────────────
    pids = set(cohort.patientunitstayid)
    cohort["sex"] = cohort.gender.map({"Male": "Male", "Female": "Female"}).fillna(
        "Other"
    )
    cohort["is_female"] = (cohort.sex == "Female").astype(int)
    cohort["baseline_egfr"] = compute_egfr(
        cohort.baseline_cr, cohort.age_num, cohort.is_female == 1
    )
    cohort["bmi"] = np.where(
        (cohort.admissionheight > 0) & (cohort.admissionweight > 0),
        cohort.admissionweight / (cohort.admissionheight / 100) ** 2,
        np.nan,
    )
    n_bmi_out = int(cohort.bmi.notna().sum() - cohort.bmi.between(10, 80).sum())
    cohort.loc[cohort.bmi.notna() & ~cohort.bmi.between(10, 80), "bmi"] = np.nan
    print(f"  BMI: capped {n_bmi_out} outliers outside [10, 80] to NaN")

    # Comorbidities
    ph_elig = ph[ph.patientunitstayid.isin(pids)] if len(ph) > 0 else pd.DataFrame()
    for cmorb, keywords in cfg.COMORBIDITY_KEYWORDS.items():
        flagged = set()
        if len(ph_elig) > 0:
            flagged = set(
                ph_elig[
                    matches_any(ph_elig.pasthistorypath, keywords)
                    | matches_any(ph_elig.pasthistoryvalue, keywords)
                ].patientunitstayid
            )
        cohort[f"hx_{cmorb}"] = cohort.patientunitstayid.isin(flagged).astype(int)

    # Nephrotoxins (unchanged)
    ad = t["admissionDrug"]
    med = t["medication"]
    ad_elig = ad[ad.patientunitstayid.isin(pids)] if len(ad) > 0 else pd.DataFrame()
    med_elig = (
        med[med.patientunitstayid.isin(pids) & (med.drugordercancelled != "Yes")]
        if len(med) > 0
        else pd.DataFrame()
    )
    mg_offsets = dict(zip(cohort.patientunitstayid, cohort.mg_offset))

    for drug_class, patterns in cfg.NEPHROTOXIN_CLASSES.items():
        flagged = set()
        if len(ad_elig) > 0 and "drugname" in ad_elig.columns:
            flagged |= set(
                ad_elig[matches_any(ad_elig.drugname, patterns)].patientunitstayid
            )
        if len(med_elig) > 0:
            for pid_i in pids:
                mg_off = mg_offsets.get(pid_i)
                if mg_off is None:
                    continue
                pt_meds = med_elig[
                    (med_elig.patientunitstayid == pid_i)
                    & (med_elig.drugstartoffset <= mg_off)
                ]
                if len(pt_meds) > 0 and matches_any(pt_meds.drugname, patterns).any():
                    flagged.add(pid_i)
        cohort[f"nephrotox_{drug_class}"] = cohort.patientunitstayid.isin(
            flagged
        ).astype(int)

    # Mg supplementation (unchanged)
    mg_supp_pids = set()
    if len(med_elig) > 0:
        mg_meds = med_elig[matches_any(med_elig.drugname, cfg.MG_SUPP_DRUG_PATTERNS)]
        for pid_i in pids:
            mo = mg_offsets.get(pid_i)
            if mo is None:
                continue
            pt_mg = mg_meds[
                (mg_meds.patientunitstayid == pid_i)
                & (mg_meds.drugstartoffset >= mo)
                & (mg_meds.drugstartoffset <= mo + cfg.MG_SUPP_GRACE_MIN)
            ]
            if len(pt_mg) > 0:
                mg_supp_pids.add(pid_i)
    cohort["mg_supplementation"] = cohort.patientunitstayid.isin(mg_supp_pids).astype(
        int
    )

    # ◆ Flag: pre-lab Mg supplementation (for strict sensitivity exclusion)
    pre_lab_pids = set()
    if len(med_elig) > 0:
        for pid_i in pids:
            mo = mg_offsets.get(pid_i)
            if mo is None:
                continue
            pt_pre = mg_meds[
                (mg_meds.patientunitstayid == pid_i)
                & (mg_meds.drugstartoffset < mo)
                & (mg_meds.drugstartoffset >= 0)
            ]
            if len(pt_pre) > 0:
                pre_lab_pids.add(pid_i)
    # Also check infusionDrug
    _inf = t["infusionDrug"]
    if len(_inf) > 0 and "drugname" in _inf.columns:
        mg_inf_all = _inf[
            _inf.patientunitstayid.isin(pids)
            & matches_any(_inf.drugname, cfg.MG_SUPP_DRUG_PATTERNS)
        ]
        for pid_i in pids:
            mo = mg_offsets.get(pid_i)
            if mo is None:
                continue
            pt_inf = mg_inf_all[
                (mg_inf_all.patientunitstayid == pid_i)
                & (mg_inf_all.infusionoffset < mo)
                & (mg_inf_all.infusionoffset >= 0)
            ]
            if len(pt_inf) > 0:
                pre_lab_pids.add(pid_i)
    cohort["pre_lab_mg_supp"] = cohort.patientunitstayid.isin(pre_lab_pids).astype(int)
    n_pre = cohort.pre_lab_mg_supp.sum()
    print(f"  Mg supplementation: {cohort.mg_supplementation.sum()}")
    print(
        f"  ◆ Pre-lab Mg supp (strict exclusion flag): {n_pre} ({100*n_pre/len(cohort):.1f}%)"
    )

    # K+ supplementation (unchanged)
    k_supp_pids = set()
    if len(med_elig) > 0:
        k_meds = med_elig[matches_any(med_elig.drugname, cfg.K_SUPP_DRUG_PATTERNS)]
        k_early = k_meds[
            (k_meds.drugstartoffset >= 0)
            & (k_meds.drugstartoffset <= cfg.MG_SUPP_GRACE_MIN)
        ]
        k_supp_pids |= set(k_early.patientunitstayid)
    inf = t["infusionDrug"]
    if len(inf) > 0 and "drugname" in inf.columns:
        k_inf = inf[
            inf.patientunitstayid.isin(pids)
            & matches_any(inf.drugname, cfg.K_SUPP_DRUG_PATTERNS)
            & (inf.infusionoffset >= 0)
            & (inf.infusionoffset <= cfg.MG_SUPP_GRACE_MIN)
        ]
        k_supp_pids |= set(k_inf.patientunitstayid)
    cohort["k_supp"] = cohort.patientunitstayid.isin(k_supp_pids).astype(int)
    cohort["ac_group"] = "neither"
    cohort.loc[(cohort.mg_supplementation == 1) & (cohort.k_supp == 1), "ac_group"] = (
        "mg_k"
    )
    cohort.loc[(cohort.mg_supplementation == 1) & (cohort.k_supp == 0), "ac_group"] = (
        "mg_only"
    )
    cohort.loc[(cohort.mg_supplementation == 0) & (cohort.k_supp == 1), "ac_group"] = (
        "k_only"
    )

    # β-blocker, steroid, vasopressor, antiarrhythmic (unchanged)
    for drug_name, patterns, col_name in [
        ("betablocker", cfg.BETA_BLOCKER_PATTERNS, "has_betablocker"),
        ("steroid", cfg.STEROID_PATTERNS, "has_steroid"),
    ]:
        dpids = set()
        if len(ad_elig) > 0 and "drugname" in ad_elig.columns:
            dpids |= set(
                ad_elig[matches_any(ad_elig.drugname, patterns)].patientunitstayid
            )
        if len(med_elig) > 0:
            early = med_elig[
                matches_any(med_elig.drugname, patterns)
                & (med_elig.drugstartoffset <= cfg.MG_WINDOW_MIN)
            ]
            dpids |= set(early.patientunitstayid)
        cohort[col_name] = cohort.patientunitstayid.isin(dpids).astype(int)

    vaso_pids = set()
    if len(inf) > 0 and "drugname" in inf.columns:
        vaso_inf = inf[
            inf.patientunitstayid.isin(pids)
            & matches_any(inf.drugname, cfg.VASOPRESSOR_PATTERNS)
            & (inf.infusionoffset >= 0)
            & (inf.infusionoffset <= cfg.MG_WINDOW_MIN)
        ]
        vaso_pids |= set(vaso_inf.patientunitstayid)
    if len(med_elig) > 0:
        vaso_med = med_elig[
            matches_any(med_elig.drugname, cfg.VASOPRESSOR_PATTERNS)
            & (med_elig.drugstartoffset >= 0)
            & (med_elig.drugstartoffset <= cfg.MG_WINDOW_MIN)
        ]
        vaso_pids |= set(vaso_med.patientunitstayid)
    cohort["has_vasopressor"] = cohort.patientunitstayid.isin(vaso_pids).astype(int)

    preop_aa_pids = set()
    if len(ad_elig) > 0 and "drugname" in ad_elig.columns:
        preop_aa_pids = set(
            ad_elig[
                matches_any(ad_elig.drugname, cfg.POAF_MED_PATTERNS)
            ].patientunitstayid
        )
    cohort["preop_antiarrhythmic"] = cohort.patientunitstayid.isin(
        preop_aa_pids
    ).astype(int)

    # First HR, Ca, K, lactate (unchanged)
    try:
        vp = pd.read_csv(
            os.path.join(cfg.DATA_ROOT, "vitalPeriodic.csv.gz"),
            usecols=["patientunitstayid", "observationoffset", "heartrate"],
            dtype={"patientunitstayid": int},
        )
        vp.columns = vp.columns.str.lower()
        vp_first = (
            vp[
                vp.patientunitstayid.isin(pids)
                & (vp.observationoffset >= 0)
                & (vp.observationoffset <= 60)
                & vp.heartrate.between(20, 250)
            ]
            .sort_values("observationoffset")
            .groupby("patientunitstayid")
            .first()
            .reset_index()
        )
        cohort = cohort.merge(
            vp_first[["patientunitstayid", "heartrate"]].rename(
                columns={"heartrate": "first_hr"}
            ),
            on="patientunitstayid",
            how="left",
        )
    except Exception:
        pass

    for lab_patterns, col, lo, hi in [
        (cfg.CA_LABNAMES, "first_ca_value", cfg.CA_PLAUSIBLE_MIN, cfg.CA_PLAUSIBLE_MAX),
        (cfg.K_LABNAMES, "first_k_value", cfg.K_PLAUSIBLE_MIN, cfg.K_PLAUSIBLE_MAX),
        (
            cfg.LACTATE_LABNAMES,
            "first_lactate",
            cfg.LACTATE_PLAUSIBLE_MIN,
            cfg.LACTATE_PLAUSIBLE_MAX,
        ),
    ]:
        elec = lab[
            lab.patientunitstayid.isin(pids)
            & matches_any(lab.labname, lab_patterns)
            & lab.labresult.between(lo, hi)
            & (lab.labresultoffset >= 0)
            & (lab.labresultoffset <= cfg.MG_WINDOW_MIN)
        ]
        first_e = (
            elec.sort_values("labresultoffset")
            .groupby("patientunitstayid")
            .first()
            .reset_index()[["patientunitstayid", "labresult"]]
            .rename(columns={"labresult": col})
        )
        cohort = cohort.merge(first_e, on="patientunitstayid", how="left")

    # Mortality (unchanged)
    cohort["icu_mortality"] = (
        cohort.unitdischargestatus.str.lower() == "expired"
    ).astype(int)
    cohort["hosp_mortality"] = (
        cohort.hospitaldischargestatus.str.lower() == "expired"
    ).astype(int)
    cohort["death_offset_min"] = np.where(
        cohort.hospitaldischargestatus.str.lower() == "expired",
        cohort.hospitaldischargeoffset,
        np.nan,
    )

    # RRT (unchanged)
    tx = t["treatment"]
    rrt_pids = set()
    if len(tx) > 0:
        rrt_pids = set(
            tx[
                tx.patientunitstayid.isin(pids)
                & matches_any(
                    tx.treatmentstring, ["dialysis", "crrt", "hemodialysis", "cvvh"]
                )
            ].patientunitstayid
        )
    cohort["rrt_7d"] = cohort.patientunitstayid.isin(rrt_pids).astype(int)

    # ◆ Perioperative transfusion (complexity-specific negative control)
    transfusion_pids = set()
    if len(tx) > 0:
        transfusion_pids = set(
            tx[
                tx.patientunitstayid.isin(pids)
                & matches_any(
                    tx.treatmentstring,
                    [
                        "transfusion",
                        "blood product",
                        "packed red",
                        "prbc",
                        "red blood cell",
                        "ffp",
                        "fresh frozen",
                        "platelet",
                        "cryoprecipitate",
                        "blood administration",
                    ],
                )
            ].patientunitstayid
        )
    cohort["nc_transfusion"] = cohort.patientunitstayid.isin(transfusion_pids).astype(
        int
    )
    print(
        f"  Perioperative transfusion: {cohort.nc_transfusion.sum()} "
        f"({100*cohort.nc_transfusion.mean():.1f}%)"
    )

    # POAF (unchanged)
    dx = t["diagnosis"]
    dx_elig = dx[dx.patientunitstayid.isin(pids)] if len(dx) > 0 else pd.DataFrame()
    preexist_af_pids = set()
    if len(ph_elig) > 0:
        preexist_af_pids = set(
            ph_elig[
                matches_any(ph_elig.pasthistorypath, cfg.PREEXISTING_AF_PATTERNS)
                | matches_any(ph_elig.pasthistoryvalue, cfg.PREEXISTING_AF_PATTERNS)
            ].patientunitstayid
        )
    cohort["preexisting_af"] = cohort.patientunitstayid.isin(preexist_af_pids).astype(
        int
    )

    poaf_dx_pids = set()
    if len(dx_elig) > 0 and "diagnosisstring" in dx_elig.columns:
        poaf_dx = dx_elig[
            matches_any(dx_elig.diagnosisstring, cfg.POAF_DX_PATTERNS)
            & (dx_elig.diagnosisoffset >= 0)
            & (dx_elig.diagnosisoffset <= cfg.POAF_WINDOW_MIN)
        ]
        poaf_dx_pids = set(poaf_dx.patientunitstayid) - preexist_af_pids
    cohort["poaf"] = cohort.patientunitstayid.isin(poaf_dx_pids).astype(int)
    cohort.loc[cohort.preexisting_af == 1, "poaf"] = np.nan

    # Negative controls + neuro (unchanged)
    for nc_name, nc_patterns in cfg.NEGATIVE_CONTROL_DX.items():
        nc_pids = set()
        if len(dx_elig) > 0 and "diagnosisstring" in dx_elig.columns:
            nc_pids = set(
                dx_elig[
                    matches_any(dx_elig.diagnosisstring, nc_patterns)
                    & (dx_elig.diagnosisoffset >= 0)
                ].patientunitstayid
            )
        cohort[f"nc_{nc_name}"] = cohort.patientunitstayid.isin(nc_pids).astype(int)
    for neuro_name, neuro_patterns in cfg.NEURO_DX_PATTERNS.items():
        neuro_pids = set()
        if len(dx_elig) > 0 and "diagnosisstring" in dx_elig.columns:
            neuro_pids = set(
                dx_elig[
                    matches_any(dx_elig.diagnosisstring, neuro_patterns)
                    & (dx_elig.diagnosisoffset >= 0)
                ].patientunitstayid
            )
        cohort[f"neuro_{neuro_name}"] = cohort.patientunitstayid.isin(
            neuro_pids
        ).astype(int)

    # Follow-up Mg (unchanged)
    mg_fu = lab[
        lab.patientunitstayid.isin(pids)
        & matches_any(lab.labname, cfg.MG_LABNAMES)
        & lab.labresult.between(cfg.MG_PLAUSIBLE_MIN, cfg.MG_PLAUSIBLE_MAX)
        & (lab.labresultoffset >= cfg.MG_WINDOW_MIN)
        & (lab.labresultoffset <= 2880)
    ]
    fu_first = (
        mg_fu.sort_values("labresultoffset")
        .groupby("patientunitstayid")
        .first()
        .reset_index()[["patientunitstayid", "labresult"]]
        .rename(columns={"labresult": "followup_mg_value"})
    )
    cohort = cohort.merge(fu_first, on="patientunitstayid", how="left")
    cohort["delta_mg"] = cohort.followup_mg_value - cohort.first_mg_value

    # ◆ Composite outcome: AKI or death (competing risk sensitivity)
    cohort["aki_or_death"] = (
        (cohort.aki_kdigo1 == 1) | (cohort.hosp_mortality == 1)
    ).astype(int)

    # ── Save ─────────────────────────────────────────────────────
    consort_df = pd.DataFrame([{"step": k, "n": int(v)} for k, v in consort.items()])
    save(consort_df, "00_consort.csv")
    save(cohort, "01_analysis_a_cohort.csv")
    print(f"\n  eICU COMPLETE: {len(cohort):,} patients, {cohort.shape[1]} cols")
    print(f"  AKI: {cohort.aki_kdigo1.sum()} ({cohort.aki_kdigo1.mean()*100:.1f}%)")
    print(
        f"  Mg supp: {cohort.mg_supplementation.sum()} "
        f"({cohort.mg_supplementation.mean()*100:.1f}%)"
    )
    print(f"  Baseline source distribution:")
    for src, n in cohort.baseline_source.value_counts().items():
        print(f"    {src}: {n} ({100*n/len(cohort):.1f}%)")


# =====================================================================
# SECTION B: MIMIC-IV ETL
# =====================================================================
def run_mimic():
    print("\n" + "=" * 70)
    print("SECTION B: MIMIC-IV Cohort Construction (v2: pre-op Cr + 6h landmark)")
    print("=" * 70)

    # ── Load tables (unchanged) ──────────────────────────────────
    patients = pd.read_csv(gz(f"{MIMIC_HOSP}/patients.csv.gz"))
    admissions = pd.read_csv(gz(f"{MIMIC_HOSP}/admissions.csv.gz"))
    icustays = pd.read_csv(gz(f"{MIMIC_ICU}/icustays.csv.gz"))
    dx = pd.read_csv(gz(f"{MIMIC_HOSP}/diagnoses_icd.csv.gz"))
    px = pd.read_csv(gz(f"{MIMIC_HOSP}/procedures_icd.csv.gz"))

    needed_labs = set(LAB_MG + LAB_CR + LAB_K + LAB_CA + LAB_LAC)
    print("  Loading labevents...")
    lab_chunks = []
    for chunk in pd.read_csv(
        gz(f"{MIMIC_HOSP}/labevents.csv.gz"),
        usecols=["subject_id", "hadm_id", "itemid", "charttime", "valuenum"],
        dtype={"subject_id": int, "hadm_id": "Int64", "itemid": int},
        chunksize=5_000_000,
    ):
        lab_chunks.append(chunk[chunk.itemid.isin(needed_labs)])
    labs = pd.concat(lab_chunks, ignore_index=True)

    inputevents = pd.read_csv(
        gz(f"{MIMIC_ICU}/inputevents.csv.gz"),
        usecols=[
            "subject_id",
            "hadm_id",
            "stay_id",
            "itemid",
            "starttime",
            "endtime",
            "amount",
            "amountuom",
        ],
    )
    presc = pd.read_csv(
        gz(f"{MIMIC_HOSP}/prescriptions.csv.gz"),
        usecols=[
            "subject_id",
            "hadm_id",
            "starttime",
            "stoptime",
            "drug",
            "route",
            "dose_val_rx",
            "dose_unit_rx",
        ],
    )

    needed_vitals = set(VITAL_HR + VITAL_WEIGHT + VITAL_HEIGHT)
    ce_path = gz(f"{MIMIC_ICU}/chartevents.csv.gz")
    chartevents = pd.DataFrame()
    if os.path.exists(ce_path):
        print("  Loading chartevents...")
        ce_chunks = []
        for chunk in pd.read_csv(
            ce_path,
            usecols=["subject_id", "stay_id", "itemid", "charttime", "valuenum"],
            dtype={"subject_id": int, "stay_id": int, "itemid": int},
            chunksize=10_000_000,
        ):
            ce_chunks.append(chunk[chunk.itemid.isin(needed_vitals)])
        chartevents = pd.concat(ce_chunks, ignore_index=True)

    # ── Cardiac surgery cohort (unchanged) ───────────────────────
    px["icd_code"] = px["icd_code"].astype(str).str.strip()
    all_cardiac = CABG_ICD9 + VALVE_ICD9 + CABG_ICD10 + VALVE_ICD10
    cardiac_hadm_px = set(
        px[px.icd_code.str[:4].isin(all_cardiac)].hadm_id.dropna().astype(int)
    )
    cvicu_hadm = set(
        icustays[icustays.first_careunit == CVICU].hadm_id.dropna().astype(int)
    )
    cardiac_hadm = cardiac_hadm_px | cvicu_hadm
    print(f"  Cardiac (ICD + CVICU): {len(cardiac_hadm)}")

    def classify_surgery(hadm_id):
        pt_codes = set(px[px.hadm_id == hadm_id].icd_code.str[:4])
        has_cabg = bool(pt_codes & set(CABG_ICD9 + CABG_ICD10))
        has_valve = bool(pt_codes & set(VALVE_ICD9 + VALVE_ICD10))
        if has_cabg and has_valve:
            return "combined"
        if has_cabg:
            return "cabg"
        if has_valve:
            return "valve"
        return "other_cardiac"

    icustays["hadm_id"] = icustays["hadm_id"].astype("Int64")
    cardiac_icu = icustays[icustays.hadm_id.isin(cardiac_hadm)].copy()
    cardiac_icu = (
        cardiac_icu.sort_values("intime").groupby("subject_id").first().reset_index()
    )
    cardiac_icu["intime"] = pd.to_datetime(cardiac_icu["intime"])
    cardiac_icu["outtime"] = pd.to_datetime(cardiac_icu["outtime"])
    cardiac_icu = cardiac_icu.merge(
        patients[["subject_id", "gender", "anchor_age"]], on="subject_id"
    )
    cardiac_icu = cardiac_icu[cardiac_icu.anchor_age >= 18]
    cardiac_icu["age"] = cardiac_icu["anchor_age"]
    cardiac_icu["is_female"] = (cardiac_icu.gender == "F").astype(int)
    cardiac_icu["surgery_type"] = cardiac_icu.hadm_id.apply(classify_surgery)
    cohort = cardiac_icu.copy()
    hadms = set(cohort.hadm_id.dropna().astype(int))
    stays = set(cohort.stay_id)
    print(f"  Adults, first stay: {len(cohort)}")

    # ── Mg exposure (unchanged) ──────────────────────────────────
    mg_labs = labs[
        labs.itemid.isin(LAB_MG)
        & labs.hadm_id.isin(hadms)
        & labs.valuenum.between(0.5, 5.0)
    ].copy()
    mg_labs["charttime"] = pd.to_datetime(mg_labs["charttime"])
    mg_labs = mg_labs.merge(cohort[["stay_id", "hadm_id", "intime"]], on="hadm_id")
    mg_labs["offset_h"] = (mg_labs.charttime - mg_labs.intime).dt.total_seconds() / 3600
    mg_early = (
        mg_labs[(mg_labs.offset_h >= -1) & (mg_labs.offset_h <= LANDMARK_HOURS)]
        .sort_values("offset_h")
        .groupby("stay_id")
        .first()
        .reset_index()
        .rename(columns={"valuenum": "first_mg_value", "charttime": "mg_charttime"})
    )
    cohort = cohort.merge(
        mg_early[["stay_id", "first_mg_value", "mg_charttime"]],
        on="stay_id",
        how="inner",
    )
    print(f"  Mg within {LANDMARK_HOURS}h: {len(cohort)}")

    # ══════════════════════════════════════════════════════════════
    # ◆ CHANGED: Baseline Cr = PRE-OPERATIVE (hospital admission)
    #   MIMIC labevents spans entire hospitalization, so we can get
    #   pre-surgery labs by looking near admittime (before ICU).
    # ══════════════════════════════════════════════════════════════
    hadms = set(cohort.hadm_id.dropna().astype(int))
    cr_labs = labs[
        labs.itemid.isin(LAB_CR)
        & labs.hadm_id.isin(hadms)
        & labs.valuenum.between(0.1, 25.0)
    ].copy()
    cr_labs["charttime"] = pd.to_datetime(cr_labs["charttime"])

    # Merge both intime (ICU) and admittime (hospital)
    adm_times = admissions[["hadm_id", "admittime"]].copy()
    adm_times["admittime"] = pd.to_datetime(adm_times["admittime"])
    cr_labs = cr_labs.merge(cohort[["stay_id", "hadm_id", "intime"]], on="hadm_id")
    cr_labs = cr_labs.merge(adm_times, on="hadm_id", how="left")

    cr_labs["offset_h_icu"] = (
        cr_labs.charttime - cr_labs.intime
    ).dt.total_seconds() / 3600
    cr_labs["offset_h_admit"] = (
        cr_labs.charttime - cr_labs.admittime
    ).dt.total_seconds() / 3600

    # ── PRIMARY: most recent pre-ICU Cr (术前最后一个) ────────────
    pre_icu = cr_labs[cr_labs.offset_h_icu < 0]
    bl_primary = (
        pre_icu.sort_values("offset_h_icu", ascending=False)
        .groupby("stay_id")
        .first()
        .reset_index()[["stay_id", "valuenum"]]
        .rename(columns={"valuenum": "baseline_cr"})
    )
    bl_primary["baseline_source"] = "last_preop"
    n_primary = len(bl_primary)

    cohort = cohort.merge(bl_primary, on="stay_id", how="inner")
    cohort = cohort[cohort.baseline_cr < 4.0]
    cohort["baseline_egfr"] = cohort.apply(
        lambda r: compute_egfr_scalar(r.baseline_cr, r.age, r.is_female), axis=1
    )

    # ── SENSITIVITY: admission Cr (closest to admittime) ──────────
    admit_cr = cr_labs[
        cr_labs.stay_id.isin(set(cohort.stay_id))
        & (cr_labs.offset_h_admit >= -24)
        & (cr_labs.offset_h_admit <= 24)
        & (cr_labs.offset_h_icu <= 0)
    ].copy()
    admit_cr["dist_admit"] = abs(admit_cr.offset_h_admit)
    bl_admit_sens = (
        admit_cr.sort_values("dist_admit")
        .groupby("stay_id")
        .first()
        .reset_index()[["stay_id", "valuenum"]]
        .rename(columns={"valuenum": "baseline_cr_admit"})
    )
    cohort = cohort.merge(bl_admit_sens, on="stay_id", how="left")
    n_admit = bl_admit_sens.stay_id.nunique()

    stays = set(cohort.stay_id)
    hadms = set(cohort.hadm_id.dropna().astype(int))

    print(f"  ◆ Primary baseline: last pre-op Cr (n={n_primary})")
    print(
        f"  ◆ Sensitivity: admission Cr available for {n_admit} "
        f"({100*n_admit/len(cohort):.0f}%)"
    )
    print(f"  With baseline Cr: {len(cohort)}")

    # ══════════════════════════════════════════════════════════════
    # ◆ CHANGED: 6h landmark exclusion
    # ══════════════════════════════════════════════════════════════
    n_before = len(cohort)
    cohort["los_min"] = (cohort.outtime - cohort.intime).dt.total_seconds() / 60

    # Death offset
    adm = admissions[["hadm_id", "hospital_expire_flag", "deathtime", "dischtime"]]
    cohort = cohort.merge(adm, on="hadm_id", how="left")
    cohort["hosp_mortality"] = cohort.hospital_expire_flag.fillna(0).astype(int)
    cohort["death_offset_min"] = (
        pd.to_datetime(cohort.deathtime, errors="coerce") - cohort.intime
    ).dt.total_seconds() / 60

    # Exclude discharged before landmark
    cohort = cohort[cohort.los_min > LANDMARK_MIN].copy()
    # Exclude died before landmark
    died_pre = cohort.death_offset_min.notna() & (
        cohort.death_offset_min <= LANDMARK_MIN
    )
    cohort = cohort[~died_pre].copy()

    stays = set(cohort.stay_id)
    hadms = set(cohort.hadm_id.dropna().astype(int))
    print(f"  ◆ 6h landmark: excluded {n_before - len(cohort)}, {len(cohort)} remain")

    # ══════════════════════════════════════════════════════════════
    # ◆ CHANGED: AKI phenotype — only Cr after 6h landmark
    #   + saves aki_time_offset for KM curves
    # ══════════════════════════════════════════════════════════════
    landmark_dt = cohort[["stay_id", "intime"]].copy()
    landmark_dt["landmark_time"] = landmark_dt.intime + pd.Timedelta(
        minutes=LANDMARK_MIN
    )

    cr_post = cr_labs[cr_labs.stay_id.isin(stays)].merge(
        cohort[["stay_id", "baseline_cr"]], on="stay_id"
    )
    cr_post = cr_post.merge(landmark_dt[["stay_id", "landmark_time"]], on="stay_id")

    # ◆ Only Cr AFTER 6h landmark
    cr_post = cr_post[cr_post.charttime > cr_post.landmark_time]
    cr_post["cr_ratio"] = cr_post.valuenum / cr_post.baseline_cr
    cr_post["cr_delta"] = cr_post.valuenum - cr_post.baseline_cr
    cr_post["hours_post_lm"] = (
        cr_post.charttime - cr_post.landmark_time
    ).dt.total_seconds() / 3600
    cr_post["offset_min_from_icu"] = (
        cr_post.charttime - cr_post.intime
    ).dt.total_seconds() / 60

    # ◆ Prevalent AKI: Cr in 0–6h that meets KDIGO vs pre-op baseline
    #   cr_labs already has offset_h_icu from earlier merge — reuse it
    cr_prelm = cr_labs[cr_labs.stay_id.isin(stays)].merge(
        cohort[["stay_id", "baseline_cr"]], on="stay_id"
    )
    cr_prelm = cr_prelm[
        (cr_prelm.offset_h_icu >= 0) & (cr_prelm.offset_h_icu <= LANDMARK_HOURS)
    ]
    prevalent_stays = set()
    for sid, g in cr_prelm.groupby("stay_id"):
        bl = g.baseline_cr.iloc[0]
        if bl > 0 and (g.valuenum / bl).max() >= 1.5:
            prevalent_stays.add(sid)

    aki = (
        cr_post.groupby("stay_id")
        .apply(
            lambda g: pd.Series(
                {
                    "aki_primary": int((g.cr_ratio >= 1.5).any()),
                    "aki_kdigo1": int(
                        (g.cr_ratio >= 1.5).any()
                        | ((g.cr_delta >= 0.3) & (g.hours_post_lm <= 48)).any()
                    ),
                    "aki_delta03": int(
                        ((g.cr_delta >= 0.3) & (g.hours_post_lm <= 48)).any()
                    ),
                    "aki_stage2": int((g.cr_ratio >= 2.0).any()),
                    "aki_stage3": int((g.cr_ratio >= 3.0).any()),
                    "peak_cr": g.valuenum.max(),
                    # ◆ NEW: save AKI onset time for KM
                    "aki_time_offset": (
                        g[g.cr_ratio >= 1.5].offset_min_from_icu.min()
                        if (g.cr_ratio >= 1.5).any()
                        else np.nan
                    ),
                }
            )
        )
        .reset_index()
    )

    cohort = cohort.merge(aki, on="stay_id", how="left")
    for c in ["aki_primary", "aki_kdigo1", "aki_delta03", "aki_stage2", "aki_stage3"]:
        cohort[c] = cohort[c].fillna(0).astype(int)
    cohort["peak_cr_ratio"] = cohort["peak_cr"] / cohort["baseline_cr"]

    # ◆ Exclude prevalent AKI
    cohort["prevalent_aki"] = cohort.stay_id.isin(prevalent_stays).astype(int)
    n_prevalent = cohort.prevalent_aki.sum()
    cohort = cohort[cohort.prevalent_aki == 0].copy()
    stays = set(cohort.stay_id)
    hadms = set(cohort.hadm_id.dropna().astype(int))

    print(f"  ◆ Prevalent AKI (pre-op→6h): {n_prevalent} excluded")
    print(
        f"  AKI KDIGO≥1: {cohort.aki_kdigo1.sum()} ({cohort.aki_kdigo1.mean()*100:.1f}%)"
    )

    # ── Covariates (unchanged) ───────────────────────────────────
    dx["icd_code"] = dx["icd_code"].astype(str).str.strip()
    for name, code_map in COMORB_ICD.items():
        cohort[name] = cohort.hadm_id.isin(matches_icd(dx, hadms, code_map)).astype(int)

    presc["starttime"] = pd.to_datetime(presc["starttime"], errors="coerce")
    presc_with = presc[presc.hadm_id.isin(hadms)].merge(
        cohort[["hadm_id", "stay_id", "intime"]], on="hadm_id"
    )
    presc_with["offset_h"] = (
        presc_with.starttime - presc_with.intime
    ).dt.total_seconds() / 3600
    presc_early = presc_with[presc_with.offset_h.between(-24, LANDMARK_HOURS)]

    for name, patterns in NEPHROTOX_MIMIC.items():
        cohort[name] = cohort.stay_id.isin(
            set(presc_early[drug_match(presc_early.drug, patterns)].stay_id)
        ).astype(int)

    ie_with = inputevents[inputevents.stay_id.isin(stays)].copy()
    ie_with["starttime"] = pd.to_datetime(ie_with["starttime"], errors="coerce")
    ie_with = ie_with.merge(cohort[["stay_id", "intime"]], on="stay_id")
    ie_with["offset_h"] = (ie_with.starttime - ie_with.intime).dt.total_seconds() / 3600
    ie_early = ie_with[ie_with.offset_h.between(0, LANDMARK_HOURS)]

    bb_stays = set(presc_early[drug_match(presc_early.drug, BB_DRUGS)].stay_id)
    bb_stays |= set(ie_early[ie_early.itemid.isin(METO_ITEMS)].stay_id)
    cohort["has_betablocker"] = cohort.stay_id.isin(bb_stays).astype(int)
    cohort["has_steroid"] = cohort.stay_id.isin(
        set(presc_early[drug_match(presc_early.drug, STEROID_DRUGS)].stay_id)
    ).astype(int)
    cohort["preop_antiarrhythmic"] = cohort.stay_id.isin(
        set(presc_early[drug_match(presc_early.drug, ANTIARR_DRUGS)].stay_id)
    ).astype(int)
    cohort["has_vasopressor"] = cohort.stay_id.isin(
        set(ie_early[ie_early.itemid.isin(VASO_ITEMS)].stay_id)
    ).astype(int)

    # Electrolytes (unchanged)
    for lab_ids, col, lo, hi in [
        (LAB_K, "first_k_value", 1.5, 8.0),
        (LAB_CA, "first_ca_value", 4.0, 15.0),
        (LAB_LAC, "first_lactate", 0.1, 30.0),
    ]:
        elec = labs[
            labs.itemid.isin(lab_ids)
            & labs.hadm_id.isin(hadms)
            & labs.valuenum.between(lo, hi)
        ].copy()
        elec["charttime"] = pd.to_datetime(elec["charttime"])
        elec = elec.merge(cohort[["stay_id", "hadm_id", "intime"]], on="hadm_id")
        elec["offset_h"] = (elec.charttime - elec.intime).dt.total_seconds() / 3600
        first_e = (
            elec[elec.offset_h.between(-1, LANDMARK_HOURS)]
            .sort_values("offset_h")
            .groupby("stay_id")
            .first()
            .reset_index()
        )
        cohort = cohort.merge(
            first_e[["stay_id", "valuenum"]].rename(columns={"valuenum": col}),
            on="stay_id",
            how="left",
        )

    # HR + BMI from chartevents (unchanged)
    if len(chartevents) > 0:
        hr = chartevents[
            chartevents.itemid.isin(VITAL_HR) & chartevents.stay_id.isin(stays)
        ].copy()
        hr["charttime"] = pd.to_datetime(hr["charttime"])
        hr = hr.merge(cohort[["stay_id", "intime"]], on="stay_id")
        hr["offset_h"] = (hr.charttime - hr.intime).dt.total_seconds() / 3600
        first_hr = (
            hr[hr.offset_h.between(0, 1) & hr.valuenum.between(20, 250)]
            .sort_values("offset_h")
            .groupby("stay_id")
            .first()
            .reset_index()
        )
        cohort = cohort.merge(
            first_hr[["stay_id", "valuenum"]].rename(columns={"valuenum": "first_hr"}),
            on="stay_id",
            how="left",
        )
        for vitals, col in [(VITAL_WEIGHT, "weight_kg"), (VITAL_HEIGHT, "height_cm")]:
            vt = chartevents[
                chartevents.itemid.isin(vitals) & chartevents.stay_id.isin(stays)
            ]
            if len(vt) > 0:
                cohort = cohort.merge(
                    vt.groupby("stay_id")["valuenum"]
                    .first()
                    .reset_index()
                    .rename(columns={"valuenum": col}),
                    on="stay_id",
                    how="left",
                )
        if "weight_kg" in cohort.columns and "height_cm" in cohort.columns:
            cohort["bmi"] = cohort.weight_kg / ((cohort.height_cm / 100) ** 2)
            cohort.loc[~cohort.bmi.between(10, 80), "bmi"] = np.nan

    # Mg supplementation + dose (unchanged)
    mg_supp = ie_early[ie_early.itemid.isin(MG_SUPP_ITEMS)]
    cohort["mg_supplementation"] = cohort.stay_id.isin(set(mg_supp.stay_id)).astype(int)
    if "amount" in mg_supp.columns:
        dose = (
            mg_supp.groupby("stay_id")["amount"]
            .sum()
            .reset_index()
            .rename(columns={"amount": "mg_total_dose"})
        )
        cohort = cohort.merge(dose, on="stay_id", how="left")
        cohort["mg_total_dose"] = cohort.mg_total_dose.fillna(0)

    # ◆ Flag: pre-lab Mg supplementation (for strict sensitivity)
    mg_all_ie = ie_with[ie_with.itemid.isin(MG_SUPP_ITEMS)]
    pre_lab_stays = set()
    for _, row in cohort.iterrows():
        sid = row.stay_id
        mg_time = row.mg_charttime
        if pd.isna(mg_time):
            continue
        pt_ie = mg_all_ie[
            (mg_all_ie.stay_id == sid) & (mg_all_ie.starttime < pd.to_datetime(mg_time))
        ]
        if len(pt_ie) > 0:
            pre_lab_stays.add(sid)
    cohort["pre_lab_mg_supp"] = cohort.stay_id.isin(pre_lab_stays).astype(int)
    n_pre = cohort.pre_lab_mg_supp.sum()
    print(
        f"  ◆ Pre-lab Mg supp (strict exclusion flag): {n_pre} ({100*n_pre/len(cohort):.1f}%)"
    )

    # K+ supplementation (unchanged)
    cohort["k_supp"] = cohort.stay_id.isin(
        set(ie_early[ie_early.itemid.isin(K_SUPP_ITEMS)].stay_id)
    ).astype(int)
    cohort["ac_group"] = "neither"
    cohort.loc[(cohort.mg_supplementation == 1) & (cohort.k_supp == 1), "ac_group"] = (
        "mg_k"
    )
    cohort.loc[(cohort.mg_supplementation == 1) & (cohort.k_supp == 0), "ac_group"] = (
        "mg_only"
    )
    cohort.loc[(cohort.mg_supplementation == 0) & (cohort.k_supp == 1), "ac_group"] = (
        "k_only"
    )

    # Negative controls + neuro + VT (unchanged)
    for name, code_map in {**NC_ICD, **NEURO_ICD}.items():
        cohort[name] = cohort.hadm_id.isin(matches_icd(dx, hadms, code_map)).astype(int)
    cohort["vent_arrhythmia"] = cohort.hadm_id.isin(
        matches_icd(dx, hadms, VT_ICD)
    ).astype(int)

    # ◆ Perioperative transfusion (complexity-specific NC)
    BLOOD_ITEMS = [
        225168,
        220970,
        225170,
        225171,  # PRBC, FFP, PLT, Cryo
        226368,
        226369,
        226370,
        226371,
    ]  # OR blood products
    blood_stays = set(
        ie_with[ie_with.itemid.isin(BLOOD_ITEMS) & ie_with.stay_id.isin(stays)].stay_id
    )
    cohort["nc_transfusion"] = cohort.stay_id.isin(blood_stays).astype(int)
    print(
        f"  Perioperative transfusion: {cohort.nc_transfusion.sum()} "
        f"({100*cohort.nc_transfusion.mean():.1f}%)"
    )

    # Follow-up Mg (unchanged)
    mg_fu = mg_labs[
        mg_labs.stay_id.isin(stays)
        & (mg_labs.offset_h >= LANDMARK_HOURS)
        & (mg_labs.offset_h <= 48)
    ]
    mg_fu_first = mg_fu.sort_values("offset_h").groupby("stay_id").first().reset_index()
    cohort = cohort.merge(
        mg_fu_first[["stay_id", "valuenum"]].rename(
            columns={"valuenum": "followup_mg_value"}
        ),
        on="stay_id",
        how="left",
    )
    cohort["delta_mg"] = cohort.followup_mg_value - cohort.first_mg_value

    # ── POAF with prior-admission fix (unchanged) ────────────────
    print("\n  POAF phenotype (prior-admission fix)...")
    admissions_full = admissions.copy()
    admissions_full["admittime"] = pd.to_datetime(admissions_full["admittime"])
    subjects = set(cohort.subject_id)
    all_subj_adm = admissions_full[admissions_full.subject_id.isin(subjects)]
    all_subj_dx = dx[dx.hadm_id.isin(set(all_subj_adm.hadm_id))]

    all_af_hadms = set()
    for ver, prefixes in [(9, AF_ICD9), (10, AF_ICD10)]:
        v = all_subj_dx[all_subj_dx.icd_version == ver]
        for p in prefixes:
            all_af_hadms |= set(v[v.icd_code.str.startswith(p)].hadm_id)

    cohort_admit = cohort[["subject_id", "hadm_id"]].merge(
        admissions_full[["hadm_id", "admittime"]], on="hadm_id"
    )
    preexist_af_subjects = set()
    current_af_hadms = set()
    for _, row in cohort_admit.iterrows():
        sid = row.subject_id
        current_hadm = int(row.hadm_id)
        prior = all_subj_adm[
            (all_subj_adm.subject_id == sid)
            & (all_subj_adm.admittime < row.admittime)
            & (all_subj_adm.hadm_id != current_hadm)
        ]
        if set(prior.hadm_id) & all_af_hadms:
            preexist_af_subjects.add(sid)
        if current_hadm in all_af_hadms:
            current_af_hadms.add(current_hadm)

    poaf_hadms = current_af_hadms - {
        int(r.hadm_id)
        for _, r in cohort[cohort.subject_id.isin(preexist_af_subjects)].iterrows()
    }
    cohort["preexisting_af"] = cohort.subject_id.isin(preexist_af_subjects).astype(int)
    cohort["poaf"] = cohort.hadm_id.astype(int).isin(poaf_hadms).astype(int)
    cohort.loc[cohort.preexisting_af == 1, "poaf"] = np.nan

    n_elig = (cohort.preexisting_af == 0).sum()
    n_poaf = cohort.poaf.sum()
    print(f"    Pre-existing AF: {cohort.preexisting_af.sum()}")
    print(f"    New-onset POAF: {int(n_poaf)} / {n_elig} ({100*n_poaf/n_elig:.1f}%)")
    print(
        f"  Mg supp: {cohort.mg_supplementation.sum()} "
        f"({cohort.mg_supplementation.mean()*100:.1f}%)"
    )

    # ◆ Composite outcome: AKI or death (competing risk sensitivity)
    cohort["aki_or_death"] = (
        (cohort.aki_kdigo1 == 1) | (cohort.hosp_mortality == 1)
    ).astype(int)

    # ── Save ─────────────────────────────────────────────────────
    out = os.path.join(cfg.RESULTS, "04_mimic_cohort.csv")
    cohort.to_csv(out, index=False)
    print(f"\n  MIMIC COMPLETE: {out}  ({len(cohort)} × {len(cohort.columns)})")
    print(f"  AKI: {cohort.aki_kdigo1.sum()} ({cohort.aki_kdigo1.mean()*100:.1f}%)")
    print(f"  Baseline source distribution:")
    for src, n in cohort.baseline_source.value_counts().items():
        print(f"    {src}: {n} ({100*n/len(cohort):.1f}%)")


# =====================================================================
# MAIN
# =====================================================================
if __name__ == "__main__":
    args = [a.lower() for a in sys.argv[1:]]
    run_both = len(args) == 0
    if run_both or "eicu" in args:
        run_eicu()
    if run_both or "mimic" in args:
        run_mimic()
    print("\n✓ ETL complete (v2: pre-op Cr + 6h landmark). Next: Rscript 02_analysis.R")
