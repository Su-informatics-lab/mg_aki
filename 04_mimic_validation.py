#!/usr/bin/env python3
"""
MIMIC-IV External Validation: Mg Supplementation → AKI After Cardiac Surgery
=============================================================================
Replicates the eICU TTE-B analysis in MIMIC-IV (BIDMC, 2008-2019).

Usage:
    python 04_mimic_validation.py

Expects MIMIC-IV v3.1 at ~/mg_aki/mimic-iv-3.1/
"""

import os
import warnings

import pandas as pd

warnings.filterwarnings("ignore")

# ─── Configuration ──────────────────────────────────────────────────
MIMIC_ROOT = os.path.expanduser("~/mg_aki/mimic-iv-3.1")
RESULTS = os.path.expanduser("~/mg_aki/results")
os.makedirs(RESULTS, exist_ok=True)

MG_WINDOW_H = 6  # Mg lab/supplementation within 6h of ICU admission
FOLLOWUP_WINDOW_H = 168  # 7 days

# ICD-9/10 codes for cardiac surgery
CARDIAC_SURGERY_ICD9_PROC = [
    "3610",
    "3611",
    "3612",
    "3613",
    "3614",
    "3615",
    "3616",
    "3617",
    "3619",  # CABG
    "3521",
    "3522",
    "3523",
    "3524",
    "3525",
    "3526",
    "3527",
    "3528",  # Valve repair
    "3511",
    "3512",
    "3513",
    "3514",  # Valve replacement
]
CARDIAC_SURGERY_ICD10_PROC = [
    "0210",
    "0211",
    "0212",
    "0213",  # CABG (bypass)
    "02RF",
    "02RG",
    "02RH",
    "02RJ",  # Valve replacement
    "02QF",
    "02QG",
    "02QH",
    "02QJ",  # Valve repair
]

# Lab item IDs (from d_labitems — will be verified at runtime)
# Common MIMIC-IV itemids:
#   Mg: 50960 (labevents)
#   Cr: 50912 (labevents)
#   K:  50971 (labevents)
#   Ca: 50893 (labevents)

# Mg supplementation itemids in inputevents (from d_items)
# Will search for "magnesium" in d_items at runtime


def load_gz(path, **kwargs):
    """Load a gzipped CSV, handling both .csv.gz and .csv"""
    if os.path.exists(path):
        return pd.read_csv(path, **kwargs)
    elif os.path.exists(path.replace(".csv.gz", ".csv")):
        return pd.read_csv(path.replace(".csv.gz", ".csv"), **kwargs)
    else:
        print(f"  WARNING: {path} not found")
        return pd.DataFrame()


def main():
    print("=" * 70)
    print("MIMIC-IV External Validation: Mg Supplementation → AKI")
    print("=" * 70)

    # ── Step 0: Load tables ─────────────────────────────────────────
    print("\nStep 0: Loading MIMIC-IV tables...")

    hosp = os.path.join(MIMIC_ROOT, "hosp")
    icu = os.path.join(MIMIC_ROOT, "icu")

    patients = load_gz(os.path.join(hosp, "patients.csv.gz"))
    admissions = load_gz(os.path.join(hosp, "admissions.csv.gz"))
    dx = load_gz(os.path.join(hosp, "diagnoses_icd.csv.gz"))
    px = load_gz(os.path.join(hosp, "procedures_icd.csv.gz"))
    labevents = load_gz(
        os.path.join(hosp, "labevents.csv.gz"),
        usecols=[
            "subject_id",
            "hadm_id",
            "itemid",
            "charttime",
            "valuenum",
            "valueuom",
        ],
        dtype={"subject_id": int, "hadm_id": "Int64", "itemid": int},
    )
    d_labitems = load_gz(os.path.join(hosp, "d_labitems.csv.gz"))

    icustays = load_gz(os.path.join(icu, "icustays.csv.gz"))
    inputevents = load_gz(os.path.join(icu, "inputevents.csv.gz"))
    d_items = load_gz(os.path.join(icu, "d_items.csv.gz"))

    for name, df in [
        ("patients", patients),
        ("admissions", admissions),
        ("diagnoses_icd", dx),
        ("procedures_icd", px),
        ("labevents", labevents),
        ("icustays", icustays),
        ("inputevents", inputevents),
    ]:
        print(f"  {name}: {len(df):,} rows")

    # ── Step 1: Identify cardiac surgery patients ───────────────────
    print("\nStep 1: Cardiac surgery cohort...")

    # By ICD procedure codes
    px["icd_code"] = px["icd_code"].astype(str).str.strip()
    cardiac_px_hadm = set()
    for code in CARDIAC_SURGERY_ICD9_PROC + CARDIAC_SURGERY_ICD10_PROC:
        matches = px[px.icd_code.str.startswith(code)]
        cardiac_px_hadm |= set(matches.hadm_id.dropna().astype(int))
    print(f"  Cardiac surgery by ICD procedure: {len(cardiac_px_hadm)} admissions")

    # By CSRU (cardiac surgery recovery unit) careunit
    csru_stays = icustays[
        icustays.first_careunit.str.contains(
            "CSRU|Cardiac|CVICU|SICU", case=False, na=False
        )
    ]
    cardiac_icu_hadm = set(csru_stays.hadm_id.dropna().astype(int))
    print(f"  Cardiac surgery by ICU unit type: {len(cardiac_icu_hadm)} admissions")

    cardiac_hadm = cardiac_px_hadm | cardiac_icu_hadm
    print(f"  Combined: {len(cardiac_hadm)} unique admissions")

    # Link to ICU stays
    icustays["hadm_id"] = icustays["hadm_id"].astype("Int64")
    cardiac_icu = icustays[icustays.hadm_id.isin(cardiac_hadm)].copy()

    # First ICU stay per patient
    cardiac_icu = (
        cardiac_icu.sort_values("intime").groupby("subject_id").first().reset_index()
    )
    print(f"  First ICU stay per patient: {len(cardiac_icu)}")

    # Parse times
    cardiac_icu["intime"] = pd.to_datetime(cardiac_icu["intime"])
    cardiac_icu["outtime"] = pd.to_datetime(cardiac_icu["outtime"])

    # Adults only
    cardiac_icu = cardiac_icu.merge(
        patients[["subject_id", "anchor_age", "gender"]], on="subject_id", how="left"
    )
    cardiac_icu = cardiac_icu[cardiac_icu.anchor_age >= 18]
    print(f"  Adults: {len(cardiac_icu)}")

    pids = set(cardiac_icu.subject_id)
    hadms = set(cardiac_icu.hadm_id.dropna().astype(int))

    # ── Step 2: Mg labs ─────────────────────────────────────────────
    print("\nStep 2: Magnesium exposure...")

    # Find Mg lab itemid
    mg_items = d_labitems[
        d_labitems.label.str.contains("magnesium", case=False, na=False)
    ]
    print(f"  Mg lab items found: {len(mg_items)}")
    if len(mg_items) > 0:
        print(f"    {mg_items[['itemid', 'label']].to_string(index=False)}")
    mg_itemids = set(mg_items.itemid)

    # Get Mg labs for cardiac surgery patients
    mg_labs = labevents[
        labevents.hadm_id.isin(hadms)
        & labevents.itemid.isin(mg_itemids)
        & labevents.valuenum.between(0.5, 5.0)
    ].copy()
    mg_labs["charttime"] = pd.to_datetime(mg_labs["charttime"])
    print(f"  Mg measurements (cardiac patients): {len(mg_labs):,}")

    # Merge with ICU admission time to compute offset
    mg_labs = mg_labs.merge(
        cardiac_icu[["subject_id", "hadm_id", "stay_id", "intime"]],
        on=["subject_id", "hadm_id"],
        how="inner",
    )
    mg_labs["offset_h"] = (mg_labs.charttime - mg_labs.intime).dt.total_seconds() / 3600

    # First Mg within 6h
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
    print(f"  Patients with Mg within 6h: {len(mg_early)}")

    cohort = cardiac_icu.merge(
        mg_early[["stay_id", "first_mg_value", "mg_charttime"]],
        on="stay_id",
        how="inner",
    )
    print(f"  Cohort with Mg: {len(cohort)}")

    # ── Step 3: Baseline creatinine ─────────────────────────────────
    print("\nStep 3: Baseline creatinine...")

    cr_items = d_labitems[
        d_labitems.label.str.contains("creatinine", case=False, na=False)
        & ~d_labitems.label.str.contains("urine|ratio|clearance", case=False, na=False)
    ]
    print(f"  Cr lab items: {cr_items[['itemid', 'label']].to_string(index=False)}")
    cr_itemids = set(cr_items.itemid)

    cr_labs = labevents[
        labevents.hadm_id.isin(hadms)
        & labevents.itemid.isin(cr_itemids)
        & labevents.valuenum.between(0.1, 25.0)
    ].copy()
    cr_labs["charttime"] = pd.to_datetime(cr_labs["charttime"])
    cr_labs = cr_labs.merge(
        cohort[["stay_id", "subject_id", "hadm_id", "intime"]],
        on=["subject_id", "hadm_id"],
        how="inner",
    )
    cr_labs["offset_h"] = (cr_labs.charttime - cr_labs.intime).dt.total_seconds() / 3600

    # Baseline: lowest pre-ICU Cr (offset ≤ 0), fallback first 12h
    pre_cr = (
        cr_labs[cr_labs.offset_h <= 0]
        .sort_values("valuenum")
        .groupby("stay_id")
        .first()
    )
    fallback_cr = (
        cr_labs[(cr_labs.offset_h >= -1) & (cr_labs.offset_h <= 12)]
        .sort_values("offset_h")
        .groupby("stay_id")
        .first()
    )

    baseline = pre_cr[["valuenum"]].rename(columns={"valuenum": "baseline_cr"})
    fb = fallback_cr[["valuenum"]].rename(columns={"valuenum": "baseline_cr"})
    baseline = baseline.combine_first(fb).reset_index()

    cohort = cohort.merge(baseline, on="stay_id", how="inner")
    cohort = cohort[cohort.baseline_cr < 4.0]
    print(f"  With baseline Cr (< 4.0): {len(cohort)}")

    # ── Step 4: AKI phenotyping ─────────────────────────────────────
    print("\nStep 4: AKI phenotyping...")

    # Post-Mg Cr measurements
    cr_post = cr_labs.merge(
        cohort[["stay_id", "mg_charttime"]], on="stay_id", how="inner"
    )
    cr_post = cr_post[cr_post.charttime > cr_post.mg_charttime]
    cr_post = cr_post.merge(
        cohort[["stay_id", "baseline_cr"]], on="stay_id", how="inner"
    )
    cr_post["cr_ratio"] = cr_post.valuenum / cr_post.baseline_cr
    cr_post["cr_delta"] = cr_post.valuenum - cr_post.baseline_cr
    cr_post["hours_post_mg"] = (
        cr_post.charttime - cr_post.mg_charttime
    ).dt.total_seconds() / 3600

    # AKI definitions
    for stay, grp in cr_post.groupby("stay_id"):
        mask = cohort.stay_id == stay
        cohort.loc[mask, "aki_primary"] = int((grp.cr_ratio >= 1.5).any())
        cohort.loc[mask, "aki_kdigo1"] = int(
            (grp.cr_ratio >= 1.5).any()
            | ((grp.cr_delta >= 0.3) & (grp.hours_post_mg <= 48)).any()
        )
        cohort.loc[mask, "aki_stage2"] = int((grp.cr_ratio >= 2.0).any())
        cohort.loc[mask, "aki_stage3"] = int((grp.cr_ratio >= 3.0).any())
        # Time-windowed
        g48 = grp[grp.hours_post_mg <= 48]
        cohort.loc[mask, "aki_primary_48h"] = (
            int((g48.cr_ratio >= 1.5).any()) if len(g48) > 0 else 0
        )

    for col in [
        "aki_primary",
        "aki_kdigo1",
        "aki_stage2",
        "aki_stage3",
        "aki_primary_48h",
    ]:
        cohort[col] = cohort[col].fillna(0).astype(int)

    print(
        f"  AKI KDIGO ≥1: {cohort.aki_kdigo1.sum()} ({cohort.aki_kdigo1.mean()*100:.1f}%)"
    )
    print(
        f"  AKI ratio ≥1.5×: {cohort.aki_primary.sum()} ({cohort.aki_primary.mean()*100:.1f}%)"
    )
    print(f"  AKI ≤48h: {cohort.aki_primary_48h.sum()}")

    # ── Step 5: Mg supplementation ──────────────────────────────────
    print("\nStep 5: Mg supplementation...")

    # Find Mg supplementation itemids in d_items (inputevents)
    mg_input_items = d_items[
        d_items.label.str.contains("magnesium", case=False, na=False)
        & ~d_items.label.str.contains("lab|level|test", case=False, na=False)
    ]
    print(f"  Mg input items found:")
    if len(mg_input_items) > 0:
        print(f"    {mg_input_items[['itemid', 'label']].to_string(index=False)}")
    mg_input_ids = set(mg_input_items.itemid)

    # Mg supplementation within 6h of ICU admission
    mg_inputs = inputevents[
        inputevents.stay_id.isin(set(cohort.stay_id))
        & inputevents.itemid.isin(mg_input_ids)
    ].copy()
    mg_inputs["starttime"] = pd.to_datetime(mg_inputs["starttime"])
    mg_inputs = mg_inputs.merge(
        cohort[["stay_id", "intime"]], on="stay_id", how="inner"
    )
    mg_inputs["offset_h"] = (
        mg_inputs.starttime - mg_inputs.intime
    ).dt.total_seconds() / 3600

    mg_supp_early = mg_inputs[
        (mg_inputs.offset_h >= 0) & (mg_inputs.offset_h <= MG_WINDOW_H)
    ]
    mg_supp_stays = set(mg_supp_early.stay_id)
    cohort["mg_supplementation"] = cohort.stay_id.isin(mg_supp_stays).astype(int)
    print(
        f"  Mg supplementation within 6h: {cohort.mg_supplementation.sum()} "
        f"({cohort.mg_supplementation.mean()*100:.1f}%)"
    )

    # Dose information (MIMIC advantage over eICU!)
    if "amount" in mg_supp_early.columns:
        dose_per_stay = mg_supp_early.groupby("stay_id")["amount"].sum().reset_index()
        dose_per_stay = dose_per_stay.rename(columns={"amount": "mg_total_dose"})
        cohort = cohort.merge(dose_per_stay, on="stay_id", how="left")
        cohort["mg_total_dose"] = cohort["mg_total_dose"].fillna(0)
        treated = cohort[cohort.mg_supplementation == 1]
        print(
            f"  Dose (treated): median={treated.mg_total_dose.median():.0f}, "
            f"IQR=[{treated.mg_total_dose.quantile(0.25):.0f}, "
            f"{treated.mg_total_dose.quantile(0.75):.0f}]"
        )

    # ── Step 6: Covariates ──────────────────────────────────────────
    print("\nStep 6: Covariates...")

    # Demographics
    cohort["age"] = cohort["anchor_age"]
    cohort["is_female"] = (cohort["gender"] == "F").astype(int)

    # Hospital mortality
    adm_info = admissions[["hadm_id", "hospital_expire_flag", "deathtime"]]
    cohort = cohort.merge(adm_info, on="hadm_id", how="left")
    cohort["hosp_mortality"] = cohort["hospital_expire_flag"].fillna(0).astype(int)
    print(f"  Hospital mortality: {cohort.hosp_mortality.sum()}")

    # Electrolytes (K+, Ca2+)
    for lab_name, col_name, lo, hi in [
        ("potassium", "first_k_value", 1.5, 8.0),
        ("calcium", "first_ca_value", 4.0, 15.0),
    ]:
        items = d_labitems[
            d_labitems.label.str.lower().str.contains(lab_name, na=False)
            & ~d_labitems.label.str.contains("urine|ratio", case=False, na=False)
        ]
        if len(items) > 0:
            elec = labevents[
                labevents.hadm_id.isin(hadms)
                & labevents.itemid.isin(set(items.itemid))
                & labevents.valuenum.between(lo, hi)
            ].copy()
            elec["charttime"] = pd.to_datetime(elec["charttime"])
            elec = elec.merge(
                cohort[["stay_id", "hadm_id", "intime"]], on="hadm_id", how="inner"
            )
            elec["offset_h"] = (elec.charttime - elec.intime).dt.total_seconds() / 3600
            first_elec = (
                elec[(elec.offset_h >= -1) & (elec.offset_h <= MG_WINDOW_H)]
                .sort_values("offset_h")
                .groupby("stay_id")
                .first()
                .reset_index()
            )
            first_elec = first_elec[["stay_id", "valuenum"]].rename(
                columns={"valuenum": col_name}
            )
            cohort = cohort.merge(first_elec, on="stay_id", how="left")
            print(f"  {col_name}: {cohort[col_name].notna().sum()} available")

    # ── Step 7: Save ────────────────────────────────────────────────
    out_path = os.path.join(RESULTS, "04_mimic_cohort.csv")
    cohort.to_csv(out_path, index=False)
    print(f"\nSaved: {out_path}  ({len(cohort)} × {len(cohort.columns)})")

    # ── Summary ─────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("MIMIC-IV VALIDATION COHORT SUMMARY")
    print("=" * 70)
    print(f"  N = {len(cohort)}")
    print(
        f"  Treated (Mg supp): {cohort.mg_supplementation.sum()} "
        f"({cohort.mg_supplementation.mean()*100:.1f}%)"
    )
    print(
        f"  AKI KDIGO ≥1: {cohort.aki_kdigo1.sum()} ({cohort.aki_kdigo1.mean()*100:.1f}%)"
    )
    print(
        f"  AKI ratio ≥1.5×: {cohort.aki_primary.sum()} ({cohort.aki_primary.mean()*100:.1f}%)"
    )
    print(f"  AKI ≤48h: {cohort.aki_primary_48h.sum()}")
    print(
        f"  Hospital mortality: {cohort.hosp_mortality.sum()} "
        f"({cohort.hosp_mortality.mean()*100:.1f}%)"
    )
    print(f"\n  Next: Run IPTW analysis using 05_mimic_tte.R")


if __name__ == "__main__":
    main()
