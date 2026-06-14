#!/usr/bin/env python3
"""
MIMIC-IV External Validation ETL — Full PS Model
=================================================
Mirrors eICU 01_etl.py with 30-covariate PS model:
  Demographics (4): age, sex, BMI, surgery_type
  Comorbidities (8): CHF, HTN, DM, CKD, COPD, PVD, stroke, liver
  Renal (2): baseline Cr, eGFR
  Nephrotoxins (4): loop diuretics, NSAIDs, ACEi/ARB, PPIs
  Medications (3): β-blockers, steroids, pre-op antiarrhythmics
  Electrolytes (3): K+, Ca2+, Mg (TTE-B)
  Hemodynamics (2): first HR, vasopressor use
"""

import os
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

MIMIC = os.path.expanduser("~/mg_aki/mimic-iv-3.1")
RESULTS = os.path.expanduser("~/mg_aki/results")
os.makedirs(RESULTS, exist_ok=True)

MG_WINDOW_H = 6
HOSP = os.path.join(MIMIC, "hosp")
ICU = os.path.join(MIMIC, "icu")


def gz(path):
    if os.path.exists(path):
        return path
    return path.replace(".csv.gz", ".csv")


# ─── Lab / Med ItemIDs (from skill) ────────────────────────────────
LAB_MG = [50960]  # serum Mg only
LAB_CR = [50912, 52546]  # serum Cr only
LAB_K = [50971]  # serum K
LAB_CA = [50893]  # total Ca
LAB_LAC = [50813]  # lactate

VITAL_HR = [220045]
VITAL_WEIGHT = [226512]  # admission weight kg
VITAL_HEIGHT = [226730]  # height cm

VASO_ITEMS = [221906, 221289, 222315, 221749, 221662, 221653, 221986]
MG_SUPP_ITEMS = [222011, 227523]  # MgSO4 + bolus (exclude 227524 OB-GYN)
AMIO_ITEMS = [221347, 228339, 229654, 230034]
METO_ITEMS = [225974]

# ICD codes
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

NEPHROTOX = {
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

# Negative controls / neuro
NC_ICD = {
    "nc_fracture": {
        9: [
            "800",
            "801",
            "802",
            "803",
            "804",
            "805",
            "806",
            "807",
            "808",
            "809",
            "810",
            "811",
            "812",
            "813",
            "814",
            "815",
            "816",
            "817",
            "818",
            "819",
            "820",
            "821",
            "822",
            "823",
            "824",
            "825",
            "826",
            "827",
            "828",
            "829",
        ],
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
    """Return set of hadm_ids matching ICD codes."""
    sub = dx_df[dx_df.hadm_id.isin(hadm_ids)]
    pids = set()
    for ver, prefixes in code_map.items():
        v = sub[sub.icd_version == ver]
        for p in prefixes:
            pids |= set(v[v.icd_code.str.startswith(p)].hadm_id)
    return pids


def drug_match(series, patterns):
    """Case-insensitive substring match on drug name series."""
    s = series.str.lower().fillna("")
    mask = pd.Series(False, index=series.index)
    for p in patterns:
        mask |= s.str.contains(p, na=False)
    return mask


def ckd_epi_2021(cr, age, is_female):
    """CKD-EPI 2021 (race-free)."""
    if pd.isna(cr) or pd.isna(age) or cr <= 0:
        return np.nan
    k = 0.7 if is_female else 0.9
    a = -0.241 if is_female else -0.302
    cr_k = cr / k
    eGFR = 142 * (min(cr_k, 1) ** a) * (max(cr_k, 1) ** (-1.200)) * (0.9938**age)
    if is_female:
        eGFR *= 1.012
    return eGFR


def main():
    print("=" * 70)
    print("MIMIC-IV FULL ETL — Matching eICU 30-Variable PS Model")
    print("=" * 70)

    # ================================================================
    # STEP 0: Load tables
    # ================================================================
    print("\nSTEP 0: Loading tables...")
    patients = pd.read_csv(gz(f"{HOSP}/patients.csv.gz"))
    admissions = pd.read_csv(gz(f"{HOSP}/admissions.csv.gz"))
    icustays = pd.read_csv(gz(f"{ICU}/icustays.csv.gz"))
    dx = pd.read_csv(gz(f"{HOSP}/diagnoses_icd.csv.gz"))
    px = pd.read_csv(gz(f"{HOSP}/procedures_icd.csv.gz"))
    d_lab = pd.read_csv(gz(f"{HOSP}/d_labitems.csv.gz"))
    d_items = pd.read_csv(gz(f"{ICU}/d_items.csv.gz"))

    # Labs — load only needed itemids
    needed_labs = set(LAB_MG + LAB_CR + LAB_K + LAB_CA + LAB_LAC)
    print("  Loading labevents (filtered by itemid)...")
    lab_chunks = []
    for chunk in pd.read_csv(
        gz(f"{HOSP}/labevents.csv.gz"),
        usecols=["subject_id", "hadm_id", "itemid", "charttime", "valuenum"],
        dtype={"subject_id": int, "hadm_id": "Int64", "itemid": int},
        chunksize=5_000_000,
    ):
        lab_chunks.append(chunk[chunk.itemid.isin(needed_labs)])
    labs = pd.concat(lab_chunks, ignore_index=True)
    print(f"  labevents (filtered): {len(labs):,}")

    # InputEvents
    inputevents = pd.read_csv(
        gz(f"{ICU}/inputevents.csv.gz"),
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
    print(f"  inputevents: {len(inputevents):,}")

    # Prescriptions
    presc = pd.read_csv(
        gz(f"{HOSP}/prescriptions.csv.gz"),
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
    print(f"  prescriptions: {len(presc):,}")

    # Chartevents — load only HR, weight, height via chunks
    needed_vitals = set(VITAL_HR + VITAL_WEIGHT + VITAL_HEIGHT)
    print("  Loading chartevents (HR + weight + height only)...")
    ce_chunks = []
    ce_path = gz(f"{ICU}/chartevents.csv.gz")
    if os.path.exists(ce_path):
        for chunk in pd.read_csv(
            ce_path,
            usecols=["subject_id", "stay_id", "itemid", "charttime", "valuenum"],
            dtype={"subject_id": int, "stay_id": int, "itemid": int},
            chunksize=10_000_000,
        ):
            ce_chunks.append(chunk[chunk.itemid.isin(needed_vitals)])
        chartevents = pd.concat(ce_chunks, ignore_index=True)
        print(f"  chartevents (filtered): {len(chartevents):,}")
    else:
        chartevents = pd.DataFrame()
        print("  chartevents: NOT FOUND (vitals will be unavailable)")

    for name, df in [
        ("patients", patients),
        ("admissions", admissions),
        ("icustays", icustays),
        ("dx", dx),
        ("px", px),
    ]:
        print(f"  {name}: {len(df):,}")

    # ================================================================
    # STEP 1: Cardiac surgery cohort
    # ================================================================
    print("\n" + "=" * 70)
    print("STEP 1: Cardiac Surgery Cohort")
    print("=" * 70)

    px["icd_code"] = px["icd_code"].astype(str).str.strip()
    all_cardiac_codes = CABG_ICD9 + VALVE_ICD9 + CABG_ICD10 + VALVE_ICD10
    cardiac_px = px[px.icd_code.str[:4].isin(all_cardiac_codes)]
    cardiac_hadm_px = set(cardiac_px.hadm_id.dropna().astype(int))
    cvicu_hadm = set(
        icustays[icustays.first_careunit == CVICU].hadm_id.dropna().astype(int)
    )
    cardiac_hadm = cardiac_hadm_px | cvicu_hadm
    print(f"  By ICD procedure: {len(cardiac_hadm_px)}")
    print(f"  By CVICU unit: {len(cvicu_hadm)}")
    print(f"  Combined: {len(cardiac_hadm)}")

    # Surgery type classification
    def classify_surgery(hadm_id):
        pt_codes = set(cardiac_px[cardiac_px.hadm_id == hadm_id].icd_code.str[:4])
        has_cabg = bool(pt_codes & set(CABG_ICD9 + CABG_ICD10))
        has_valve = bool(pt_codes & set(VALVE_ICD9 + VALVE_ICD10))
        if has_cabg and has_valve:
            return "combined"
        if has_cabg:
            return "cabg"
        if has_valve:
            return "valve"
        return "other_cardiac"

    # Build cohort from ICU stays
    icustays["hadm_id"] = icustays["hadm_id"].astype("Int64")
    cardiac_icu = icustays[icustays.hadm_id.isin(cardiac_hadm)].copy()
    cardiac_icu = (
        cardiac_icu.sort_values("intime").groupby("subject_id").first().reset_index()
    )
    cardiac_icu["intime"] = pd.to_datetime(cardiac_icu["intime"])
    cardiac_icu["outtime"] = pd.to_datetime(cardiac_icu["outtime"])

    # Demographics
    cardiac_icu = cardiac_icu.merge(
        patients[["subject_id", "gender", "anchor_age"]], on="subject_id"
    )
    cardiac_icu = cardiac_icu[cardiac_icu.anchor_age >= 18]
    cardiac_icu["age"] = cardiac_icu["anchor_age"]
    cardiac_icu["is_female"] = (cardiac_icu.gender == "F").astype(int)
    cardiac_icu["surgery_type"] = cardiac_icu.hadm_id.apply(classify_surgery)
    print(f"  Adults, first ICU stay: {len(cardiac_icu)}")
    print(f"  Surgery types: {cardiac_icu.surgery_type.value_counts().to_dict()}")

    cohort = cardiac_icu.copy()
    hadms = set(cohort.hadm_id.dropna().astype(int))
    stays = set(cohort.stay_id)

    # ================================================================
    # STEP 2: Magnesium Exposure
    # ================================================================
    print("\n" + "=" * 70)
    print("STEP 2: Magnesium Exposure")
    print("=" * 70)

    mg_labs = labs[
        labs.itemid.isin(LAB_MG)
        & labs.hadm_id.isin(hadms)
        & labs.valuenum.between(0.5, 5.0)
    ].copy()
    mg_labs["charttime"] = pd.to_datetime(mg_labs["charttime"])
    mg_labs = mg_labs.merge(cohort[["stay_id", "hadm_id", "intime"]], on="hadm_id")
    mg_labs["offset_h"] = (mg_labs.charttime - mg_labs.intime).dt.total_seconds() / 3600

    mg_early = (
        mg_labs[(mg_labs.offset_h >= -1) & (mg_labs.offset_h <= MG_WINDOW_H)]
        .sort_values("offset_h")
        .groupby("stay_id")
        .first()
        .reset_index()
    )
    mg_early = mg_early.rename(
        columns={"valuenum": "first_mg_value", "charttime": "mg_charttime"}
    )
    cohort = cohort.merge(
        mg_early[["stay_id", "first_mg_value", "mg_charttime"]],
        on="stay_id",
        how="inner",
    )
    print(f"  Mg within 6h: {len(cohort)}")
    print(
        f"  Mg: median={cohort.first_mg_value.median():.1f}, "
        f"IQR=[{cohort.first_mg_value.quantile(.25):.1f}, {cohort.first_mg_value.quantile(.75):.1f}]"
    )

    # ================================================================
    # STEP 3: Baseline Creatinine
    # ================================================================
    print("\n" + "=" * 70)
    print("STEP 3: Baseline Creatinine")
    print("=" * 70)

    cr_labs = labs[
        labs.itemid.isin(LAB_CR)
        & labs.hadm_id.isin(hadms)
        & labs.valuenum.between(0.1, 25.0)
    ].copy()
    cr_labs["charttime"] = pd.to_datetime(cr_labs["charttime"])
    cr_labs = cr_labs.merge(cohort[["stay_id", "hadm_id", "intime"]], on="hadm_id")
    cr_labs["offset_h"] = (cr_labs.charttime - cr_labs.intime).dt.total_seconds() / 3600

    pre_cr = (
        cr_labs[cr_labs.offset_h <= 0]
        .sort_values("valuenum")
        .groupby("stay_id")
        .first()
    )
    fb_cr = (
        cr_labs[cr_labs.offset_h.between(-1, 12)]
        .sort_values("offset_h")
        .groupby("stay_id")
        .first()
    )
    baseline = pre_cr[["valuenum"]].rename(columns={"valuenum": "baseline_cr"})
    baseline = baseline.combine_first(
        fb_cr[["valuenum"]].rename(columns={"valuenum": "baseline_cr"})
    )
    baseline = baseline.reset_index()

    cohort = cohort.merge(baseline, on="stay_id", how="inner")
    cohort = cohort[cohort.baseline_cr < 4.0]
    print(f"  With baseline Cr: {len(cohort)}")

    # eGFR (CKD-EPI 2021)
    cohort["baseline_egfr"] = cohort.apply(
        lambda r: ckd_epi_2021(r.baseline_cr, r.age, r.is_female), axis=1
    )
    print(f"  eGFR: median={cohort.baseline_egfr.median():.1f}")

    stays = set(cohort.stay_id)
    hadms = set(cohort.hadm_id.dropna().astype(int))

    # ================================================================
    # STEP 4: AKI Phenotyping
    # ================================================================
    print("\n" + "=" * 70)
    print("STEP 4: AKI Phenotyping")
    print("=" * 70)

    cr_post = cr_labs[cr_labs.stay_id.isin(stays)].merge(
        cohort[["stay_id", "mg_charttime", "baseline_cr"]], on="stay_id"
    )
    cr_post = cr_post[cr_post.charttime > cr_post.mg_charttime]
    cr_post["cr_ratio"] = cr_post.valuenum / cr_post.baseline_cr
    cr_post["cr_delta"] = cr_post.valuenum - cr_post.baseline_cr
    cr_post["hours_post_mg"] = (
        cr_post.charttime - cr_post.mg_charttime
    ).dt.total_seconds() / 3600

    aki = (
        cr_post.groupby("stay_id")
        .apply(
            lambda g: pd.Series(
                {
                    "aki_primary": int((g.cr_ratio >= 1.5).any()),
                    "aki_kdigo1": int(
                        (g.cr_ratio >= 1.5).any()
                        | ((g.cr_delta >= 0.3) & (g.hours_post_mg <= 48)).any()
                    ),
                    "aki_delta03": int(
                        ((g.cr_delta >= 0.3) & (g.hours_post_mg <= 48)).any()
                    ),
                    "aki_stage2": int((g.cr_ratio >= 2.0).any()),
                    "aki_stage3": int((g.cr_ratio >= 3.0).any()),
                    "aki_primary_24h": (
                        int((g[g.hours_post_mg <= 24].cr_ratio >= 1.5).any())
                        if len(g[g.hours_post_mg <= 24]) > 0
                        else 0
                    ),
                    "aki_primary_48h": (
                        int((g[g.hours_post_mg <= 48].cr_ratio >= 1.5).any())
                        if len(g[g.hours_post_mg <= 48]) > 0
                        else 0
                    ),
                    "aki_primary_72h": (
                        int((g[g.hours_post_mg <= 72].cr_ratio >= 1.5).any())
                        if len(g[g.hours_post_mg <= 72]) > 0
                        else 0
                    ),
                    "peak_cr": g.valuenum.max(),
                }
            )
        )
        .reset_index()
    )
    cohort = cohort.merge(aki, on="stay_id", how="left")
    for c in [
        "aki_primary",
        "aki_kdigo1",
        "aki_delta03",
        "aki_stage2",
        "aki_stage3",
        "aki_primary_24h",
        "aki_primary_48h",
        "aki_primary_72h",
    ]:
        cohort[c] = cohort[c].fillna(0).astype(int)
    cohort["peak_cr_ratio"] = cohort["peak_cr"] / cohort["baseline_cr"]

    print(
        f"  KDIGO ≥1: {cohort.aki_kdigo1.sum()} ({cohort.aki_kdigo1.mean()*100:.1f}%)"
    )
    print(
        f"  Ratio ≥1.5×: {cohort.aki_primary.sum()} ({cohort.aki_primary.mean()*100:.1f}%)"
    )
    print(f"  48h: {cohort.aki_primary_48h.sum()}, 72h: {cohort.aki_primary_72h.sum()}")

    # ================================================================
    # STEP 5: Covariates (30-variable PS model)
    # ================================================================
    print("\n" + "=" * 70)
    print("STEP 5: Covariates")
    print("=" * 70)

    # ── 5a. Comorbidities from ICD ────────────────────────────────
    print("  Comorbidities (from ICD):")
    dx["icd_code"] = dx["icd_code"].astype(str).str.strip()
    for name, code_map in COMORB_ICD.items():
        matched = matches_icd(dx, hadms, code_map)
        cohort[name] = cohort.hadm_id.isin(matched).astype(int)
        print(f"    {name}: {cohort[name].sum()}")

    # ── 5b. Nephrotoxins from prescriptions ───────────────────────
    print("  Nephrotoxins (from prescriptions):")
    presc["starttime"] = pd.to_datetime(presc["starttime"], errors="coerce")
    presc_hadm = presc[presc.hadm_id.isin(hadms)]
    presc_with_time = presc_hadm.merge(
        cohort[["hadm_id", "stay_id", "intime"]], on="hadm_id"
    )
    presc_with_time["offset_h"] = (
        presc_with_time.starttime - presc_with_time.intime
    ).dt.total_seconds() / 3600
    presc_early = presc_with_time[presc_with_time.offset_h.between(-24, MG_WINDOW_H)]

    for name, patterns in NEPHROTOX.items():
        matched_stays = set(presc_early[drug_match(presc_early.drug, patterns)].stay_id)
        cohort[name] = cohort.stay_id.isin(matched_stays).astype(int)
        print(f"    {name}: {cohort[name].sum()}")

    # ── 5c. β-blockers ────────────────────────────────────────────
    bb_stays = set(presc_early[drug_match(presc_early.drug, BB_DRUGS)].stay_id)
    # Also from inputevents
    ie_with_time = inputevents[inputevents.stay_id.isin(stays)].copy()
    ie_with_time["starttime"] = pd.to_datetime(
        ie_with_time["starttime"], errors="coerce"
    )
    ie_with_time = ie_with_time.merge(cohort[["stay_id", "intime"]], on="stay_id")
    ie_with_time["offset_h"] = (
        ie_with_time.starttime - ie_with_time.intime
    ).dt.total_seconds() / 3600
    ie_early = ie_with_time[ie_with_time.offset_h.between(0, MG_WINDOW_H)]
    bb_stays |= set(ie_early[ie_early.itemid.isin(METO_ITEMS)].stay_id)
    cohort["has_betablocker"] = cohort.stay_id.isin(bb_stays).astype(int)
    print(
        f"  β-blocker: {cohort.has_betablocker.sum()} ({cohort.has_betablocker.mean()*100:.1f}%)"
    )

    # ── 5d. Steroids ──────────────────────────────────────────────
    steroid_stays = set(
        presc_early[drug_match(presc_early.drug, STEROID_DRUGS)].stay_id
    )
    cohort["has_steroid"] = cohort.stay_id.isin(steroid_stays).astype(int)
    print(
        f"  Steroid: {cohort.has_steroid.sum()} ({cohort.has_steroid.mean()*100:.1f}%)"
    )

    # ── 5e. Pre-op antiarrhythmics ────────────────────────────────
    aa_stays = set(presc_early[drug_match(presc_early.drug, ANTIARR_DRUGS)].stay_id)
    cohort["preop_antiarrhythmic"] = cohort.stay_id.isin(aa_stays).astype(int)
    print(f"  Pre-op antiarrhythmic: {cohort.preop_antiarrhythmic.sum()}")

    # ── 5f. Vasopressors from inputevents ─────────────────────────
    vaso_stays = set(ie_early[ie_early.itemid.isin(VASO_ITEMS)].stay_id)
    cohort["has_vasopressor"] = cohort.stay_id.isin(vaso_stays).astype(int)
    print(
        f"  Vasopressor: {cohort.has_vasopressor.sum()} ({cohort.has_vasopressor.mean()*100:.1f}%)"
    )

    # ── 5g. Electrolytes (K+, Ca2+) ──────────────────────────────
    for lab_ids, col, lo, hi in [
        (LAB_K, "first_k_value", 1.5, 8.0),
        (LAB_CA, "first_ca_value", 4.0, 15.0),
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
            elec[elec.offset_h.between(-1, MG_WINDOW_H)]
            .sort_values("offset_h")
            .groupby("stay_id")
            .first()
            .reset_index()
        )
        first_e = first_e[["stay_id", "valuenum"]].rename(columns={"valuenum": col})
        cohort = cohort.merge(first_e, on="stay_id", how="left")
        print(f"  {col}: {cohort[col].notna().sum()} available")

    # ── 5h. First HR from chartevents ─────────────────────────────
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
        first_hr = first_hr[["stay_id", "valuenum"]].rename(
            columns={"valuenum": "first_hr"}
        )
        cohort = cohort.merge(first_hr, on="stay_id", how="left")
        print(f"  First HR: {cohort.first_hr.notna().sum()} available")

        # BMI from weight + height
        wt = chartevents[
            chartevents.itemid.isin(VITAL_WEIGHT) & chartevents.stay_id.isin(stays)
        ]
        ht = chartevents[
            chartevents.itemid.isin(VITAL_HEIGHT) & chartevents.stay_id.isin(stays)
        ]
        if len(wt) > 0:
            wt_first = (
                wt.groupby("stay_id")["valuenum"]
                .first()
                .reset_index()
                .rename(columns={"valuenum": "weight_kg"})
            )
            cohort = cohort.merge(wt_first, on="stay_id", how="left")
        if len(ht) > 0:
            ht_first = (
                ht.groupby("stay_id")["valuenum"]
                .first()
                .reset_index()
                .rename(columns={"valuenum": "height_cm"})
            )
            cohort = cohort.merge(ht_first, on="stay_id", how="left")
        if "weight_kg" in cohort.columns and "height_cm" in cohort.columns:
            cohort["bmi"] = cohort.weight_kg / ((cohort.height_cm / 100) ** 2)
            cohort.loc[~cohort.bmi.between(10, 80), "bmi"] = np.nan
            print(
                f"  BMI: {cohort.bmi.notna().sum()} available, median={cohort.bmi.median():.1f}"
                if cohort.bmi.notna().sum() > 0
                else "  BMI: unavailable"
            )
    else:
        print("  Chartevents unavailable — HR and BMI will be missing")

    # ── 5i. Mg supplementation ────────────────────────────────────
    mg_supp = ie_early[ie_early.itemid.isin(MG_SUPP_ITEMS)]
    mg_supp_stays = set(mg_supp.stay_id)
    cohort["mg_supplementation"] = cohort.stay_id.isin(mg_supp_stays).astype(int)
    print(
        f"\n  Mg supplementation: {cohort.mg_supplementation.sum()} "
        f"({cohort.mg_supplementation.mean()*100:.1f}%)"
    )

    # Dose (MIMIC advantage)
    if "amount" in mg_supp.columns:
        dose = (
            mg_supp.groupby("stay_id")["amount"]
            .sum()
            .reset_index()
            .rename(columns={"amount": "mg_total_dose"})
        )
        cohort = cohort.merge(dose, on="stay_id", how="left")
        cohort["mg_total_dose"] = cohort.mg_total_dose.fillna(0)

    # ── 5j. Mortality ─────────────────────────────────────────────
    adm = admissions[["hadm_id", "hospital_expire_flag"]]
    cohort = cohort.merge(adm, on="hadm_id", how="left")
    cohort["hosp_mortality"] = cohort.hospital_expire_flag.fillna(0).astype(int)
    print(
        f"  Hospital mortality: {cohort.hosp_mortality.sum()} ({cohort.hosp_mortality.mean()*100:.1f}%)"
    )

    # ── 5k. POAF (dx-only) ───────────────────────────────────────
    preexist_af = matches_icd(dx, hadms, COMORB_ICD["hx_afib"])
    cohort["preexisting_af"] = cohort.hadm_id.isin(preexist_af).astype(int)
    af_current = matches_icd(dx, hadms, {9: ["42731"], 10: ["I48"]})
    poaf_hadm = af_current - preexist_af  # This is imperfect — same-admission AF codes
    cohort["poaf"] = cohort.hadm_id.isin(poaf_hadm).astype(int)
    cohort.loc[cohort.preexisting_af == 1, "poaf"] = np.nan
    print(f"  Pre-existing AF: {cohort.preexisting_af.sum()}")
    n_poaf = cohort.poaf.sum()
    print(f"  POAF (new-onset): {int(n_poaf)}")

    # ── 5l. Ventricular arrhythmia ────────────────────────────────
    vt_hadm = matches_icd(dx, hadms, VT_ICD)
    cohort["vent_arrhythmia"] = cohort.hadm_id.isin(vt_hadm).astype(int)
    print(
        f"  VT/VF: {cohort.vent_arrhythmia.sum()} ({cohort.vent_arrhythmia.mean()*100:.1f}%)"
    )

    # ── 5m. Negative controls + neuro ─────────────────────────────
    for name, code_map in {**NC_ICD, **NEURO_ICD}.items():
        matched = matches_icd(dx, hadms, code_map)
        cohort[name] = cohort.hadm_id.isin(matched).astype(int)
    print(
        f"  Negative controls: fracture={cohort.nc_fracture.sum()}, UTI={cohort.nc_uti.sum()}"
    )
    print(
        f"  Neuro: delirium={cohort.neuro_delirium.sum()}, enceph={cohort.neuro_encephalopathy.sum()}"
    )

    # ── 5n. Follow-up Mg (positive control) ───────────────────────
    mg_followup = mg_labs[mg_labs.stay_id.isin(stays)].copy()
    mg_fu = (
        mg_followup[
            (mg_followup.offset_h >= MG_WINDOW_H) & (mg_followup.offset_h <= 48)
        ]
        .sort_values("offset_h")
        .groupby("stay_id")
        .first()
        .reset_index()
    )
    mg_fu = mg_fu[["stay_id", "valuenum"]].rename(
        columns={"valuenum": "followup_mg_value"}
    )
    cohort = cohort.merge(mg_fu, on="stay_id", how="left")
    cohort["delta_mg"] = cohort.followup_mg_value - cohort.first_mg_value

    # ================================================================
    # SAVE
    # ================================================================
    out = os.path.join(RESULTS, "04_mimic_cohort.csv")
    cohort.to_csv(out, index=False)
    print(f"\n{'='*70}")
    print(f"SAVED: {out}  ({len(cohort)} × {len(cohort.columns)})")
    print(f"{'='*70}")

    # Summary
    print(f"\n  N = {len(cohort)}")
    print(f"  Surgery: {cohort.surgery_type.value_counts().to_dict()}")
    print(
        f"  Mg supp: {cohort.mg_supplementation.sum()} ({cohort.mg_supplementation.mean()*100:.1f}%)"
    )
    print(
        f"  AKI KDIGO≥1: {cohort.aki_kdigo1.sum()} ({cohort.aki_kdigo1.mean()*100:.1f}%)"
    )
    print(
        f"  AKI 1.5×: {cohort.aki_primary.sum()} ({cohort.aki_primary.mean()*100:.1f}%)"
    )
    print(
        f"  Mortality: {cohort.hosp_mortality.sum()} ({cohort.hosp_mortality.mean()*100:.1f}%)"
    )
    print(f"  Columns available for PS: {len(cohort.columns)}")

    # Covariate availability
    ps_vars = [
        "age",
        "is_female",
        "bmi",
        "surgery_type",
        "hx_chf",
        "hx_hypertension",
        "hx_diabetes",
        "hx_ckd",
        "hx_copd",
        "hx_pvd",
        "hx_stroke",
        "hx_liver",
        "baseline_cr",
        "baseline_egfr",
        "nephrotox_loop_diuretic",
        "nephrotox_nsaid",
        "nephrotox_acei_arb",
        "nephrotox_ppi",
        "has_betablocker",
        "has_steroid",
        "preop_antiarrhythmic",
        "first_k_value",
        "first_ca_value",
        "first_hr",
        "has_vasopressor",
        "first_mg_value",
    ]
    avail = sum(1 for v in ps_vars if v in cohort.columns)
    print(f"  PS covariates matched: {avail}/{len(ps_vars)}")
    for v in ps_vars:
        if v in cohort.columns:
            na = cohort[v].isna().sum()
            print(f"    {v}: {len(cohort)-na}/{len(cohort)} available ({na} NA)")
        else:
            print(f"    {v}: MISSING")

    print(f"\nNext: Rscript 05_mimic_tte.R")


if __name__ == "__main__":
    main()
