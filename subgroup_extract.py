#!/usr/bin/env python3
"""
subgroup_extract.py — Tier 2 variable extraction for subgroup/safety analyses

Extracts from raw eICU/MIMIC tables:
  - Vent duration (extubation time proxy)
  - IABP / ECMO use
  - Post-treatment HR (first HR after 6h)
  - Max post-treatment serum Mg (safety)
  - Alcohol history (for MDS score)

Merges into existing cohort CSVs → saves enriched versions.

Run on Tempest: python subgroup_extract.py
"""

import os
import warnings

import pandas as pd

warnings.filterwarnings("ignore")

RESULTS = os.path.expanduser("~/mg_aki/results")
EICU_ROOT = None
MIMIC_ROOT = os.path.expanduser("~/mg_aki/mimic-iv-3.1")

# Resolve eICU path (same logic as 00_config.py)
for p in [
    "~/mg_aki/eicu-crd-2.0",
    "~/mg_aki/eicu-collaborative-research-database-demo-2.0.1",
]:
    if os.path.isdir(os.path.expanduser(p)):
        EICU_ROOT = os.path.expanduser(p)
        break
if EICU_ROOT is None:
    print("WARNING: eICU data not found, skipping eICU extraction")


def gz(path):
    return path if os.path.exists(path) else path.replace(".csv.gz", ".csv")


def matches_any(series, patterns):
    pat = "|".join(patterns)
    return series.str.lower().str.contains(pat, na=False)


# =====================================================================
# eICU EXTRACTION
# =====================================================================
def extract_eicu():
    print("=" * 60)
    print("eICU: Extracting Tier 2 variables")
    print("=" * 60)

    cohort = pd.read_csv(os.path.join(RESULTS, "01_analysis_a_cohort.csv"))
    pids = set(cohort.patientunitstayid)
    print(f"  Cohort: {len(cohort)} patients")

    # ── Vent duration: first extubation from treatment table ──────
    print("\n  Vent duration (from treatment table)...")
    try:
        tx = pd.read_csv(os.path.join(EICU_ROOT, "treatment.csv"), low_memory=False)
        tx.columns = tx.columns.str.lower()
        tx_elig = tx[tx.patientunitstayid.isin(pids)]
        extub = tx_elig[
            matches_any(
                tx_elig.treatmentstring, ["extubat", "self-extub", "unplanned extub"]
            )
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
        n = cohort.vent_duration_h.notna().sum()
        print(
            f"    Extubation found: {n} patients, median {cohort.vent_duration_h.median():.1f}h"
        )
    except Exception as e:
        print(f"    Failed: {e}")

    # ── IABP / ECMO from treatment table ──────────────────────────
    print("\n  IABP / ECMO...")
    try:
        iabp_pids = set(
            tx_elig[
                matches_any(
                    tx_elig.treatmentstring,
                    ["iabp", "intra-aortic", "balloon pump", "intraaortic"],
                )
            ].patientunitstayid
        )
        ecmo_pids = set(
            tx_elig[
                matches_any(
                    tx_elig.treatmentstring, ["ecmo", "extracorporeal membrane", "ecls"]
                )
            ].patientunitstayid
        )
        cohort["has_iabp"] = cohort.patientunitstayid.isin(iabp_pids).astype(int)
        cohort["has_ecmo"] = cohort.patientunitstayid.isin(ecmo_pids).astype(int)
        print(f"    IABP: {cohort.has_iabp.sum()}, ECMO: {cohort.has_ecmo.sum()}")
    except Exception as e:
        print(f"    Failed: {e}")

    # ── Post-treatment HR (first HR after 6h) ─────────────────────
    print("\n  Post-treatment HR (>6h)...")
    try:
        vp = pd.read_csv(
            os.path.join(EICU_ROOT, "vitalPeriodic.csv.gz"),
            usecols=["patientunitstayid", "observationoffset", "heartrate"],
            dtype={"patientunitstayid": int},
        )
        vp.columns = vp.columns.str.lower()
        vp_elig = vp[
            vp.patientunitstayid.isin(pids)
            & (vp.observationoffset > 360)  # > 6h
            & (vp.observationoffset <= 1440)  # ≤ 24h
            & vp.heartrate.between(20, 250)
        ]
        post_hr = (
            vp_elig.sort_values("observationoffset")
            .groupby("patientunitstayid")
            .first()
            .reset_index()[["patientunitstayid", "heartrate"]]
            .rename(columns={"heartrate": "post_hr_6h"})
        )
        cohort = cohort.merge(post_hr, on="patientunitstayid", how="left")
        print(
            f"    Post-6h HR: {cohort.post_hr_6h.notna().sum()} available, "
            f"median {cohort.post_hr_6h.median():.0f}"
        )
    except Exception as e:
        print(f"    Failed: {e}")

    # ── Max post-treatment Mg ─────────────────────────────────────
    print("\n  Max post-treatment Mg...")
    try:
        lab = pd.read_csv(os.path.join(EICU_ROOT, "lab.csv"), low_memory=False)
        lab.columns = lab.columns.str.lower()
        mg_labs = lab[
            lab.patientunitstayid.isin(pids)
            & lab.labname.str.lower().str.contains("magnesium", na=False)
            & lab.labresult.between(0.5, 10.0)
            & (lab.labresultoffset > 360)  # after 6h (post-treatment)
        ]
        max_mg = (
            mg_labs.groupby("patientunitstayid")["labresult"]
            .max()
            .reset_index()
            .rename(columns={"labresult": "max_posttreat_mg"})
        )
        cohort = cohort.merge(max_mg, on="patientunitstayid", how="left")
        n = cohort.max_posttreat_mg.notna().sum()
        print(
            f"    Max post-Mg: {n} available, median {cohort.max_posttreat_mg.median():.1f}"
        )
        # Safety flag
        n_high = (cohort.max_posttreat_mg > 4.8).sum()
        print(f"    >4.8 mg/dL (symptomatic threshold): {n_high}")
    except Exception as e:
        print(f"    Failed: {e}")

    # ── Alcohol history (for MDS) ─────────────────────────────────
    print("\n  Alcohol history...")
    try:
        ph = pd.read_csv(os.path.join(EICU_ROOT, "pastHistory.csv"), low_memory=False)
        ph.columns = ph.columns.str.lower()
        ph_elig = ph[ph.patientunitstayid.isin(pids)]
        alc_pids = set(
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
        cohort["alcohol_history"] = cohort.patientunitstayid.isin(alc_pids).astype(int)
        print(f"    Alcohol history: {cohort.alcohol_history.sum()}")
    except Exception as e:
        print(f"    Failed: {e}")

    # ── ICU LOS (from existing data, just make explicit) ──────────
    if "unitdischargeoffset" in cohort.columns:
        cohort["icu_los_h"] = cohort.unitdischargeoffset / 60.0
        print(f"\n  ICU LOS: median {cohort.icu_los_h.median():.1f}h")

    # ── Save ──────────────────────────────────────────────────────
    out = os.path.join(RESULTS, "01_analysis_a_cohort_enriched.csv")
    cohort.to_csv(out, index=False)
    print(f"\n  Saved: {out} ({len(cohort)} × {cohort.shape[1]})")
    return cohort


# =====================================================================
# MIMIC EXTRACTION
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

    # ── Vent duration from procedureevents ─────────────────────────
    print("\n  Vent duration (from procedureevents)...")
    try:
        pe = pd.read_csv(
            gz(f"{ICU}/procedureevents.csv.gz"),
            usecols=["stay_id", "itemid", "starttime", "endtime", "value"],
        )
        pe_vent = pe[
            pe.stay_id.isin(stays)
            & pe.itemid.isin([225792, 225794])  # InvasiveVent, NonInvasiveVent
        ].copy()
        pe_vent["starttime"] = pd.to_datetime(pe_vent["starttime"])
        pe_vent["endtime"] = pd.to_datetime(pe_vent["endtime"])
        pe_vent = pe_vent.merge(cohort[["stay_id", "intime"]], on="stay_id")

        # Total vent hours = sum of all vent episodes
        pe_vent["vent_h"] = (
            pe_vent.endtime - pe_vent.starttime
        ).dt.total_seconds() / 3600
        vent_total = (
            pe_vent.groupby("stay_id")["vent_h"]
            .sum()
            .reset_index()
            .rename(columns={"vent_h": "vent_duration_h"})
        )
        cohort = cohort.merge(vent_total, on="stay_id", how="left")
        n = cohort.vent_duration_h.notna().sum()
        print(
            f"    Vent records: {n} patients, median {cohort.vent_duration_h.median():.1f}h"
        )
    except Exception as e:
        print(f"    Failed: {e}")

    # ── IABP / ECMO from procedureevents ──────────────────────────
    print("\n  IABP / ECMO...")
    try:
        # IABP: itemid 225977 or ICD procedure
        iabp_stays = set(pe[pe.stay_id.isin(stays) & pe.itemid.isin([225977])].stay_id)
        cohort["has_iabp"] = cohort.stay_id.isin(iabp_stays).astype(int)

        # ECMO: from ICD procedure codes
        px = pd.read_csv(gz(f"{HOSP}/procedures_icd.csv.gz"))
        px["icd_code"] = px["icd_code"].astype(str).str.strip()
        ecmo_icd9 = ["3965", "3966"]
        ecmo_icd10 = ["5A15223", "5A1522F", "5A1522G", "5A1522H"]
        ecmo_hadm = set()
        for ver, codes in [(9, ecmo_icd9), (10, ecmo_icd10)]:
            sub = px[(px.hadm_id.isin(hadms)) & (px.icd_version == ver)]
            for c in codes:
                ecmo_hadm |= set(sub[sub.icd_code.str.startswith(c)].hadm_id)
        cohort["has_ecmo"] = cohort.hadm_id.isin(ecmo_hadm).astype(int)

        print(f"    IABP: {cohort.has_iabp.sum()}, ECMO: {cohort.has_ecmo.sum()}")
    except Exception as e:
        print(f"    Failed: {e}")

    # ── Post-treatment HR (>6h) ───────────────────────────────────
    print("\n  Post-treatment HR (>6h)...")
    try:
        HR_ITEMID = [220045]
        ce_path = gz(f"{ICU}/chartevents.csv.gz")
        hr_chunks = []
        for chunk in pd.read_csv(
            ce_path,
            usecols=["stay_id", "itemid", "charttime", "valuenum"],
            dtype={"stay_id": int, "itemid": int},
            chunksize=10_000_000,
        ):
            hr_chunks.append(
                chunk[chunk.itemid.isin(HR_ITEMID) & chunk.stay_id.isin(stays)]
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
            f"    Post-6h HR: {cohort.post_hr_6h.notna().sum()} available, "
            f"median {cohort.post_hr_6h.median():.0f}"
        )
    except Exception as e:
        print(f"    Failed: {e}")

    # ── Max post-treatment Mg ─────────────────────────────────────
    print("\n  Max post-treatment Mg...")
    try:
        MG_ITEMID = [50960]
        labs = pd.read_csv(
            gz(f"{HOSP}/labevents.csv.gz"),
            usecols=["subject_id", "hadm_id", "itemid", "charttime", "valuenum"],
            dtype={"hadm_id": "Int64", "itemid": int},
            chunksize=5_000_000,
        )
        mg_chunks = []
        for chunk in labs:
            mg_chunks.append(
                chunk[chunk.itemid.isin(MG_ITEMID) & chunk.hadm_id.isin(hadms)]
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
        n = cohort.max_posttreat_mg.notna().sum()
        print(f"    Max post-Mg: {n}, median {cohort.max_posttreat_mg.median():.1f}")
        print(f"    >4.8 mg/dL: {(cohort.max_posttreat_mg > 4.8).sum()}")
    except Exception as e:
        print(f"    Failed: {e}")

    # ── Alcohol history from ICD ──────────────────────────────────
    print("\n  Alcohol history...")
    try:
        dx = pd.read_csv(gz(f"{HOSP}/diagnoses_icd.csv.gz"))
        dx["icd_code"] = dx["icd_code"].astype(str).str.strip()
        alc_hadm = set()
        # ICD-9: 303.x (dependence), 305.0x (abuse), 291.x (withdrawal)
        # ICD-10: F10.x
        for ver, prefixes in [(9, ["303", "3050", "291"]), (10, ["F10"])]:
            sub = dx[(dx.hadm_id.isin(hadms)) & (dx.icd_version == ver)]
            for p in prefixes:
                alc_hadm |= set(sub[sub.icd_code.str.startswith(p)].hadm_id)
        cohort["alcohol_history"] = cohort.hadm_id.isin(alc_hadm).astype(int)
        print(f"    Alcohol history: {cohort.alcohol_history.sum()}")
    except Exception as e:
        print(f"    Failed: {e}")

    # ── ICU LOS ───────────────────────────────────────────────────
    if "outtime" in cohort.columns:
        cohort["outtime"] = pd.to_datetime(cohort["outtime"])
        cohort["icu_los_h"] = (cohort.outtime - cohort.intime).dt.total_seconds() / 3600
        print(f"\n  ICU LOS: median {cohort.icu_los_h.median():.1f}h")

    # ── Save ──────────────────────────────────────────────────────
    out = os.path.join(RESULTS, "04_mimic_cohort_enriched.csv")
    cohort.to_csv(out, index=False)
    print(f"\n  Saved: {out} ({len(cohort)} × {cohort.shape[1]})")
    return cohort


# =====================================================================
if __name__ == "__main__":
    if EICU_ROOT:
        extract_eicu()
    extract_mimic()
    print("\n✓ subgroup_extract.py COMPLETE")
    print("  Next: Rscript subgroup_analysis.R")
