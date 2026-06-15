#!/usr/bin/env python3
"""
subgroup_extract.py — Tier 2 variable extraction

Extracts vent duration, IABP/ECMO, post-treatment HR, max post-Mg,
alcohol history from raw eICU/MIMIC tables. Merges into cohort CSVs.

Run on Tempest: python subgroup_extract.py
"""

import os
import warnings

import pandas as pd

warnings.filterwarnings("ignore")

RESULTS = os.path.expanduser("~/mg_aki/results")
MIMIC_ROOT = os.path.expanduser("~/mg_aki/mimic-iv-3.1")

# Resolve eICU path
EICU_ROOT = None
for p in [
    "~/mg_aki/eicu-crd-2.0",
    "~/mg_aki/eicu-collaborative-research-database-demo-2.0.1",
]:
    if os.path.isdir(os.path.expanduser(p)):
        EICU_ROOT = os.path.expanduser(p)
        break


def find_file(root, name):
    """Try name.csv.gz, name.csv, lowercase variants."""
    for n in [name, name.lower()]:
        for ext in [".csv.gz", ".csv"]:
            p = os.path.join(root, n + ext)
            if os.path.exists(p):
                return p
    raise FileNotFoundError(f"{name} not found in {root}")


def gz(path):
    return path if os.path.exists(path) else path.replace(".csv.gz", ".csv")


def matches_any(series, patterns):
    return series.str.lower().str.contains("|".join(patterns), na=False)


# =====================================================================
# eICU
# =====================================================================
def extract_eicu():
    print("=" * 60)
    print("eICU: Extracting Tier 2 variables")
    print("=" * 60)

    cohort = pd.read_csv(os.path.join(RESULTS, "01_analysis_a_cohort.csv"))
    pids = set(cohort.patientunitstayid)
    print(f"  Cohort: {len(cohort)} patients")

    # ── Load treatment table (vent + IABP/ECMO) ──────────────────
    tx_elig = None
    try:
        tx = pd.read_csv(find_file(EICU_ROOT, "treatment"), low_memory=False)
        tx.columns = tx.columns.str.lower()
        tx_elig = tx[tx.patientunitstayid.isin(pids)]
        print(f"  Treatment table: {len(tx_elig)} rows")
    except Exception as e:
        print(f"  Treatment table: {e}")

    # Vent duration
    print("\n  Vent duration...")
    if tx_elig is not None:
        try:
            extub = tx_elig[
                matches_any(tx_elig.treatmentstring, ["extubat", "self-extub"])
                & (tx_elig.treatmentoffset > 0)
            ]
            first_extub = (
                extub.sort_values("treatmentoffset")
                .groupby("patientunitstayid")
                .first()
                .reset_index()[["patientunitstayid", "treatmentoffset"]]
                .rename(columns={"treatmentoffset": "extub_offset_min"})
            )
            first_extub["vent_duration_h"] = first_extub.extub_offset_min / 60.0
            cohort = cohort.merge(first_extub, on="patientunitstayid", how="left")
            print(
                f"    {cohort.vent_duration_h.notna().sum()} patients, median {cohort.vent_duration_h.median():.1f}h"
            )
        except Exception as e:
            print(f"    Failed: {e}")

    # IABP / ECMO
    print("  IABP / ECMO...")
    if tx_elig is not None:
        try:
            iabp = set(
                tx_elig[
                    matches_any(
                        tx_elig.treatmentstring,
                        ["iabp", "intra-aortic", "balloon pump", "intraaortic"],
                    )
                ].patientunitstayid
            )
            ecmo = set(
                tx_elig[
                    matches_any(
                        tx_elig.treatmentstring,
                        ["ecmo", "extracorporeal membrane", "ecls"],
                    )
                ].patientunitstayid
            )
            cohort["has_iabp"] = cohort.patientunitstayid.isin(iabp).astype(int)
            cohort["has_ecmo"] = cohort.patientunitstayid.isin(ecmo).astype(int)
            print(f"    IABP: {cohort.has_iabp.sum()}, ECMO: {cohort.has_ecmo.sum()}")
        except Exception as e:
            print(f"    Failed: {e}")

    # Post-treatment HR (>6h)
    print("  Post-treatment HR...")
    try:
        vp = pd.read_csv(
            find_file(EICU_ROOT, "vitalPeriodic"),
            usecols=["patientunitstayid", "observationoffset", "heartrate"],
            dtype={"patientunitstayid": int},
        )
        vp.columns = vp.columns.str.lower()
        post_hr = (
            vp[
                vp.patientunitstayid.isin(pids)
                & (vp.observationoffset > 360)
                & (vp.observationoffset <= 1440)
                & vp.heartrate.between(20, 250)
            ]
            .sort_values("observationoffset")
            .groupby("patientunitstayid")
            .first()
            .reset_index()[["patientunitstayid", "heartrate"]]
            .rename(columns={"heartrate": "post_hr_6h"})
        )
        cohort = cohort.merge(post_hr, on="patientunitstayid", how="left")
        print(
            f"    {cohort.post_hr_6h.notna().sum()} available, median {cohort.post_hr_6h.median():.0f}"
        )
    except Exception as e:
        print(f"    Failed: {e}")

    # Max post-treatment Mg
    print("  Max post-treatment Mg...")
    try:
        lab = pd.read_csv(find_file(EICU_ROOT, "lab"), low_memory=False)
        lab.columns = lab.columns.str.lower()
        mg = lab[
            lab.patientunitstayid.isin(pids)
            & lab.labname.str.lower().str.contains("magnesium", na=False)
            & lab.labresult.between(0.5, 10.0)
            & (lab.labresultoffset > 360)
        ]
        max_mg = (
            mg.groupby("patientunitstayid")["labresult"]
            .max()
            .reset_index()
            .rename(columns={"labresult": "max_posttreat_mg"})
        )
        cohort = cohort.merge(max_mg, on="patientunitstayid", how="left")
        n = cohort.max_posttreat_mg.notna().sum()
        print(f"    {n} available, median {cohort.max_posttreat_mg.median():.1f}")
        print(f"    >4.8 mg/dL: {(cohort.max_posttreat_mg > 4.8).sum()}")
    except Exception as e:
        print(f"    Failed: {e}")

    # Alcohol history
    print("  Alcohol history...")
    try:
        ph = pd.read_csv(find_file(EICU_ROOT, "pastHistory"), low_memory=False)
        ph.columns = ph.columns.str.lower()
        ph_elig = ph[ph.patientunitstayid.isin(pids)]
        alc = set(
            ph_elig[
                matches_any(
                    ph_elig.pasthistorypath,
                    ["alcohol", "etoh", "alcoholism", "drinking"],
                )
                | matches_any(
                    ph_elig.pasthistoryvalue,
                    ["alcohol", "etoh", "alcoholism", "heavy drink"],
                )
            ].patientunitstayid
        )
        cohort["alcohol_history"] = cohort.patientunitstayid.isin(alc).astype(int)
        print(f"    {cohort.alcohol_history.sum()} patients")
    except Exception as e:
        print(f"    Failed: {e}")

    # ICU LOS
    if "unitdischargeoffset" in cohort.columns:
        cohort["icu_los_h"] = cohort.unitdischargeoffset / 60.0
        print(f"  ICU LOS: median {cohort.icu_los_h.median():.1f}h")

    out = os.path.join(RESULTS, "01_analysis_a_cohort_enriched.csv")
    cohort.to_csv(out, index=False)
    print(f"\n  Saved: {out} ({len(cohort)} x {cohort.shape[1]})")


# =====================================================================
# MIMIC (unchanged — already worked)
# =====================================================================
def extract_mimic():
    print("\n" + "=" * 60)
    print("MIMIC-IV: Extracting Tier 2 variables")
    print("=" * 60)

    cohort = pd.read_csv(os.path.join(RESULTS, "04_mimic_cohort.csv"))
    stays = set(cohort.stay_id)
    hadms = set(cohort.hadm_id.dropna().astype(int))
    cohort["intime"] = pd.to_datetime(cohort["intime"])
    print(f"  Cohort: {len(cohort)} patients")
    ICU = os.path.join(MIMIC_ROOT, "icu")
    HOSP = os.path.join(MIMIC_ROOT, "hosp")

    # Vent duration
    print("\n  Vent duration...")
    try:
        pe = pd.read_csv(
            gz(f"{ICU}/procedureevents.csv.gz"),
            usecols=["stay_id", "itemid", "starttime", "endtime", "value"],
        )
        pe_vent = pe[pe.stay_id.isin(stays) & pe.itemid.isin([225792, 225794])].copy()
        pe_vent["starttime"] = pd.to_datetime(pe_vent["starttime"])
        pe_vent["endtime"] = pd.to_datetime(pe_vent["endtime"])
        pe_vent["vent_h"] = (
            pe_vent.endtime - pe_vent.starttime
        ).dt.total_seconds() / 3600
        vent = (
            pe_vent.groupby("stay_id")["vent_h"]
            .sum()
            .reset_index()
            .rename(columns={"vent_h": "vent_duration_h"})
        )
        cohort = cohort.merge(vent, on="stay_id", how="left")
        print(
            f"    {cohort.vent_duration_h.notna().sum()} patients, median {cohort.vent_duration_h.median():.1f}h"
        )
    except Exception as e:
        print(f"    Failed: {e}")

    # IABP / ECMO
    print("  IABP / ECMO...")
    try:
        iabp_stays = set(pe[pe.stay_id.isin(stays) & pe.itemid.isin([225977])].stay_id)
        cohort["has_iabp"] = cohort.stay_id.isin(iabp_stays).astype(int)
        px = pd.read_csv(gz(f"{HOSP}/procedures_icd.csv.gz"))
        px["icd_code"] = px["icd_code"].astype(str).str.strip()
        ecmo_hadm = set()
        for ver, codes in [
            (9, ["3965", "3966"]),
            (10, ["5A15223", "5A1522F", "5A1522G", "5A1522H"]),
        ]:
            sub = px[(px.hadm_id.isin(hadms)) & (px.icd_version == ver)]
            for c in codes:
                ecmo_hadm |= set(sub[sub.icd_code.str.startswith(c)].hadm_id)
        cohort["has_ecmo"] = cohort.hadm_id.isin(ecmo_hadm).astype(int)
        print(f"    IABP: {cohort.has_iabp.sum()}, ECMO: {cohort.has_ecmo.sum()}")
    except Exception as e:
        print(f"    Failed: {e}")

    # Post-treatment HR
    print("  Post-treatment HR...")
    try:
        hr_chunks = []
        for chunk in pd.read_csv(
            gz(f"{ICU}/chartevents.csv.gz"),
            usecols=["stay_id", "itemid", "charttime", "valuenum"],
            dtype={"stay_id": int, "itemid": int},
            chunksize=10_000_000,
        ):
            hr_chunks.append(
                chunk[chunk.itemid.isin([220045]) & chunk.stay_id.isin(stays)]
            )
        hr = pd.concat(hr_chunks, ignore_index=True)
        hr["charttime"] = pd.to_datetime(hr["charttime"])
        hr = hr.merge(cohort[["stay_id", "intime"]], on="stay_id")
        hr["offset_h"] = (hr.charttime - hr.intime).dt.total_seconds() / 3600
        post_hr = (
            hr[(hr.offset_h > 6) & (hr.offset_h <= 24) & hr.valuenum.between(20, 250)]
            .sort_values("offset_h")
            .groupby("stay_id")
            .first()
            .reset_index()[["stay_id", "valuenum"]]
            .rename(columns={"valuenum": "post_hr_6h"})
        )
        cohort = cohort.merge(post_hr, on="stay_id", how="left")
        print(
            f"    {cohort.post_hr_6h.notna().sum()} available, median {cohort.post_hr_6h.median():.0f}"
        )
    except Exception as e:
        print(f"    Failed: {e}")

    # Max post-treatment Mg
    print("  Max post-treatment Mg...")
    try:
        mg_chunks = []
        for chunk in pd.read_csv(
            gz(f"{HOSP}/labevents.csv.gz"),
            usecols=["subject_id", "hadm_id", "itemid", "charttime", "valuenum"],
            dtype={"hadm_id": "Int64", "itemid": int},
            chunksize=5_000_000,
        ):
            mg_chunks.append(
                chunk[chunk.itemid.isin([50960]) & chunk.hadm_id.isin(hadms)]
            )
        mg = pd.concat(mg_chunks, ignore_index=True)
        mg["charttime"] = pd.to_datetime(mg["charttime"])
        mg = mg.merge(cohort[["stay_id", "hadm_id", "intime"]], on="hadm_id")
        mg["offset_h"] = (mg.charttime - mg.intime).dt.total_seconds() / 3600
        max_mg = (
            mg[(mg.offset_h > 6) & mg.valuenum.between(0.5, 10.0)]
            .groupby("stay_id")["valuenum"]
            .max()
            .reset_index()
            .rename(columns={"valuenum": "max_posttreat_mg"})
        )
        cohort = cohort.merge(max_mg, on="stay_id", how="left")
        print(
            f"    {cohort.max_posttreat_mg.notna().sum()}, median {cohort.max_posttreat_mg.median():.1f}"
        )
        print(f"    >4.8 mg/dL: {(cohort.max_posttreat_mg > 4.8).sum()}")
    except Exception as e:
        print(f"    Failed: {e}")

    # Alcohol history
    print("  Alcohol history...")
    try:
        dx = pd.read_csv(gz(f"{HOSP}/diagnoses_icd.csv.gz"))
        dx["icd_code"] = dx["icd_code"].astype(str).str.strip()
        alc_hadm = set()
        for ver, prefixes in [(9, ["303", "3050", "291"]), (10, ["F10"])]:
            sub = dx[(dx.hadm_id.isin(hadms)) & (dx.icd_version == ver)]
            for p in prefixes:
                alc_hadm |= set(sub[sub.icd_code.str.startswith(p)].hadm_id)
        cohort["alcohol_history"] = cohort.hadm_id.isin(alc_hadm).astype(int)
        print(f"    {cohort.alcohol_history.sum()} patients")
    except Exception as e:
        print(f"    Failed: {e}")

    # ICU LOS
    if "outtime" in cohort.columns:
        cohort["outtime"] = pd.to_datetime(cohort["outtime"])
        cohort["icu_los_h"] = (cohort.outtime - cohort.intime).dt.total_seconds() / 3600
        print(f"  ICU LOS: median {cohort.icu_los_h.median():.1f}h")

    out = os.path.join(RESULTS, "04_mimic_cohort_enriched.csv")
    cohort.to_csv(out, index=False)
    print(f"\n  Saved: {out} ({len(cohort)} x {cohort.shape[1]})")


if __name__ == "__main__":
    if EICU_ROOT:
        extract_eicu()
    extract_mimic()
    print("\n✓ COMPLETE — Next: Rscript subgroup_analysis.R")
