#!/usr/bin/env python3
"""
01c_dose_extract.py — Extract administered Mg dose and add to cohorts

eICU:  medication + infusionDrug tables → parse dose text → grams
MIMIC: inputevents → filter amountuom='grams' → sum per stay

Also counts n_mg_administrations as a unit-free dose proxy.

New columns added to cohorts:
  mg_dose_grams      — total administered dose in grams (0 for untreated)
  n_mg_doses         — number of IV Mg administrations (0 for untreated)

Run: python 01c_dose_extract.py
Then: Rscript probe_three_analyses.R   (dose-response will now work)
"""

import os
import re
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

RESULTS = os.path.expanduser("~/mg_aki/results")
MIMIC_ICU = os.path.expanduser("~/mg_aki/mimic-iv-3.1/icu")
EICU_ROOT = None
for p in [
    "~/mg_aki/eicu-crd-2.0",
    "~/mg_aki/eicu-collaborative-research-database-demo-2.0.1",
]:
    if os.path.isdir(os.path.expanduser(p)):
        EICU_ROOT = os.path.expanduser(p)
        break


def gz(path):
    if os.path.exists(path):
        return path
    alt = path.replace(".csv.gz", ".csv") if path.endswith(".gz") else path + ".gz"
    return alt if os.path.exists(alt) else path


def parse_dose_grams(text):
    """Parse eICU dosage text → grams. Returns None if unparseable."""
    if pd.isna(text) or not isinstance(text, str):
        return None
    text = text.strip().lower()

    # Try to extract number + unit
    m = re.match(r"([\d.]+)\s*(g|gm|gram|grams|mg|milligram|meq|ml)\b", text)
    if m:
        val = float(m.group(1))
        unit = m.group(2)
        if unit in ("g", "gm", "gram", "grams"):
            return val if 0.1 <= val <= 20 else None  # plausible range
        elif unit in ("mg", "milligram"):
            g = val / 1000
            return g if 0.1 <= g <= 20 else None
        elif unit == "meq":
            # MgSO4: 1g ≈ 8.12 mEq. 16 mEq ≈ 2g
            g = val / 8.12
            return g if 0.1 <= g <= 20 else None
        elif unit == "ml":
            # Standard 40mg/mL (4%) → 50mL = 2g
            g = val * 0.04
            return g if 0.1 <= g <= 20 else None

    # Bare number — guess based on magnitude
    m2 = re.match(r"^([\d.]+)$", text)
    if m2:
        val = float(m2.group(1))
        if 0.5 <= val <= 10:
            return val  # likely grams
        elif 500 <= val <= 10000:
            return val / 1000  # likely mg

    return None


# =====================================================================
# eICU DOSE EXTRACTION
# =====================================================================
def extract_eicu():
    print("=" * 60)
    print("eICU: DOSE EXTRACTION")
    print("=" * 60)

    cohort = pd.read_csv(os.path.join(RESULTS, "01_analysis_a_cohort.csv"))
    pids = set(cohort.patientunitstayid)
    print(
        f"  Cohort: {len(cohort)} patients, {cohort.mg_supplementation.sum()} treated"
    )

    dose_data = {}  # pid → {grams: float, n_doses: int}

    # ── medication table ─────────────────────────────────────────
    for fname in ["medication.csv.gz", "medication.csv"]:
        mpath = os.path.join(EICU_ROOT, fname)
        if os.path.exists(mpath):
            print(f"  Loading {fname}...")
            med = pd.read_csv(mpath, low_memory=False)
            med.columns = med.columns.str.lower()
            med = med[med.patientunitstayid.isin(pids)]

            # Filter to IV magnesium
            mg_mask = med.drugname.str.lower().str.contains("magnesium", na=False)
            iv_mask = pd.Series(True, index=med.index)
            if "routeadmin" in med.columns:
                iv_mask = (
                    med.routeadmin.str.lower().str.contains(
                        "iv|intravenous|inject", na=False
                    )
                    | med.routeadmin.isna()
                )
            mg_med = med[mg_mask & iv_mask]

            # Time filter: within 6h of admission (360 min)
            if "drugstartoffset" in mg_med.columns:
                mg_med = mg_med[mg_med.drugstartoffset.between(0, 360)]

            print(f"  IV Mg medication entries: {len(mg_med)}")

            # Parse dose
            dose_col = None
            for c in ["dosage", "drugdosage", "dose"]:
                if c in mg_med.columns:
                    dose_col = c
                    break

            for _, row in mg_med.iterrows():
                pid = row.patientunitstayid
                if pid not in dose_data:
                    dose_data[pid] = {"grams": 0.0, "n_doses": 0, "parsed": 0}

                dose_data[pid]["n_doses"] += 1

                if dose_col and pd.notna(row.get(dose_col)):
                    g = parse_dose_grams(str(row[dose_col]))
                    if g is not None:
                        dose_data[pid]["grams"] += g
                        dose_data[pid]["parsed"] += 1
            break

    # ── infusionDrug table ───────────────────────────────────────
    for fname in [
        "infusionDrug.csv.gz",
        "infusionDrug.csv",
        "infusiondrug.csv.gz",
        "infusiondrug.csv",
    ]:
        ipath = os.path.join(EICU_ROOT, fname)
        if os.path.exists(ipath):
            print(f"  Loading {fname}...")
            inf = pd.read_csv(ipath, low_memory=False)
            inf.columns = inf.columns.str.lower()
            inf = inf[inf.patientunitstayid.isin(pids)]

            mg_inf = inf[inf.drugname.str.lower().str.contains("magnesium", na=False)]

            if "infusionoffset" in mg_inf.columns:
                mg_inf = mg_inf[mg_inf.infusionoffset.between(0, 360)]

            print(f"  IV Mg infusion entries: {len(mg_inf)}")

            for _, row in mg_inf.iterrows():
                pid = row.patientunitstayid
                if pid not in dose_data:
                    dose_data[pid] = {"grams": 0.0, "n_doses": 0, "parsed": 0}

                dose_data[pid]["n_doses"] += 1

                # Try drugamount
                for c in ["drugamount", "drugrate"]:
                    if c in row.index and pd.notna(row[c]):
                        val = (
                            float(row[c]) if isinstance(row[c], (int, float)) else None
                        )
                        if val and 0.5 <= val <= 10:
                            dose_data[pid]["grams"] += val
                            dose_data[pid]["parsed"] += 1
                            break
            break

    # ── Build dose DataFrame ─────────────────────────────────────
    dose_df = pd.DataFrame(
        [
            {
                "patientunitstayid": pid,
                "mg_dose_grams": d["grams"],
                "n_mg_doses": d["n_doses"],
                "n_parsed": d["parsed"],
            }
            for pid, d in dose_data.items()
        ]
    )

    # For patients with parsed=0 but n_doses>0, estimate: 2g per dose
    if len(dose_df) > 0:
        unparsed = dose_df[(dose_df.n_mg_doses > 0) & (dose_df.n_parsed == 0)]
        if len(unparsed) > 0:
            dose_df.loc[unparsed.index, "mg_dose_grams"] = unparsed.n_mg_doses * 2.0
            print(
                f"  {len(unparsed)} patients had unparseable dose → estimated at 2g/dose"
            )

    # ── Merge into cohort ────────────────────────────────────────
    for c in ["mg_dose_grams", "n_mg_doses", "mg_total_dose"]:
        if c in cohort.columns:
            cohort = cohort.drop(columns=[c])

    if len(dose_df) > 0:
        cohort = cohort.merge(
            dose_df[["patientunitstayid", "mg_dose_grams", "n_mg_doses"]],
            on="patientunitstayid",
            how="left",
        )
    else:
        cohort["mg_dose_grams"] = 0.0
        cohort["n_mg_doses"] = 0

    cohort.mg_dose_grams = cohort.mg_dose_grams.fillna(0)
    cohort.n_mg_doses = cohort.n_mg_doses.fillna(0).astype(int)

    # Set untreated to 0
    cohort.loc[cohort.mg_supplementation == 0, "mg_dose_grams"] = 0
    cohort.loc[cohort.mg_supplementation == 0, "n_mg_doses"] = 0

    # Also save as mg_total_dose for R script compatibility
    cohort["mg_total_dose"] = cohort["mg_dose_grams"]

    # ── Diagnostics ──────────────────────────────────────────────
    trt = cohort[cohort.mg_supplementation == 1]
    print(f"\n  ── eICU Dose Summary (treated, n={len(trt)}) ──")
    print(
        f"  mg_dose_grams: median={trt.mg_dose_grams.median():.1f}, "
        f"mean={trt.mg_dose_grams.mean():.1f}, "
        f"IQR={trt.mg_dose_grams.quantile(0.25):.1f}-{trt.mg_dose_grams.quantile(0.75):.1f}"
    )
    print(
        f"  n_mg_doses:    median={trt.n_mg_doses.median():.0f}, "
        f"range={trt.n_mg_doses.min()}-{trt.n_mg_doses.max()}"
    )

    # Dose distribution
    for lo, hi, label in [
        (0.01, 2, "0.1-2g"),
        (2.01, 4, "2.1-4g"),
        (4.01, 8, "4.1-8g"),
        (8.01, 100, ">8g"),
    ]:
        n = ((trt.mg_dose_grams >= lo) & (trt.mg_dose_grams <= hi)).sum()
        if n > 0:
            aki = trt[
                (trt.mg_dose_grams >= lo) & (trt.mg_dose_grams <= hi)
            ].aki_kdigo1.mean()
            print(f"    {label:>8s}: n={n:4d}, AKI={100*aki:.1f}%")

    # Save
    cohort.to_csv(os.path.join(RESULTS, "01_analysis_a_cohort.csv"), index=False)
    print(f"\n  Updated cohort with dose columns")


# =====================================================================
# MIMIC DOSE EXTRACTION
# =====================================================================
def extract_mimic():
    print(f"\n{'=' * 60}")
    print("MIMIC-IV: DOSE EXTRACTION (grams only)")
    print("=" * 60)

    cohort = pd.read_csv(os.path.join(RESULTS, "04_mimic_cohort.csv"))
    stays = set(cohort.stay_id)
    print(
        f"  Cohort: {len(cohort)} patients, {cohort.mg_supplementation.sum()} treated"
    )

    # Read inputevents
    ie_path = gz(os.path.join(MIMIC_ICU, "inputevents.csv.gz"))
    print(f"  Loading inputevents...")
    ie = pd.read_csv(
        ie_path,
        usecols=[
            "stay_id",
            "itemid",
            "amount",
            "amountuom",
            "starttime",
            "cancelreason",
        ],
        low_memory=False,
    )

    # Filter to Mg items, cohort stays, not cancelled
    mg_items = [222011, 227523]
    ie = ie[ie.itemid.isin(mg_items) & ie.stay_id.isin(stays)]
    if "cancelreason" in ie.columns:
        ie = ie[ie.cancelreason.isna() | (ie.cancelreason == 0)]
    ie = ie[ie.amount.notna() & (ie.amount > 0)]

    print(f"  Mg entries in cohort: {len(ie)}")
    print(f"  Unit distribution:")
    for unit, n in ie.amountuom.value_counts().items():
        sub = ie[ie.amountuom == unit]
        print(f"    {unit}: n={n}, median={sub.amount.median():.1f}")

    # Convert everything to grams
    ie["dose_g"] = np.nan
    ie.loc[ie.amountuom.str.lower() == "grams", "dose_g"] = ie.loc[
        ie.amountuom.str.lower() == "grams", "amount"
    ]
    # mL → grams: standard 40mg/mL = 0.04g/mL
    ie.loc[ie.amountuom.str.lower() == "ml", "dose_g"] = (
        ie.loc[ie.amountuom.str.lower() == "ml", "amount"] * 0.04
    )
    # mg → grams
    ie.loc[ie.amountuom.str.lower() == "mg", "dose_g"] = (
        ie.loc[ie.amountuom.str.lower() == "mg", "amount"] / 1000
    )
    # dose (assume grams if ~2)
    ie.loc[ie.amountuom.str.lower() == "dose", "dose_g"] = ie.loc[
        ie.amountuom.str.lower() == "dose", "amount"
    ]

    # Filter plausible (0.1-20g per administration)
    ie = ie[ie.dose_g.between(0.1, 20)]

    # Sum per stay
    dose_per_stay = (
        ie.groupby("stay_id")
        .agg(mg_dose_grams=("dose_g", "sum"), n_mg_doses=("dose_g", "count"))
        .reset_index()
    )

    print(f"\n  Stays with valid dose: {len(dose_per_stay)}")

    # Merge
    for c in ["mg_dose_grams", "n_mg_doses", "mg_total_dose"]:
        if c in cohort.columns:
            cohort = cohort.drop(columns=[c])

    cohort = cohort.merge(dose_per_stay, on="stay_id", how="left")
    cohort.mg_dose_grams = cohort.mg_dose_grams.fillna(0)
    cohort.n_mg_doses = cohort.n_mg_doses.fillna(0).astype(int)
    cohort.loc[cohort.mg_supplementation == 0, "mg_dose_grams"] = 0
    cohort.loc[cohort.mg_supplementation == 0, "n_mg_doses"] = 0
    cohort["mg_total_dose"] = cohort["mg_dose_grams"]

    # Diagnostics
    trt = cohort[cohort.mg_supplementation == 1]
    print(f"\n  ── MIMIC Dose Summary (treated, n={len(trt)}) ──")
    print(
        f"  mg_dose_grams: median={trt.mg_dose_grams.median():.1f}, "
        f"mean={trt.mg_dose_grams.mean():.1f}, "
        f"IQR={trt.mg_dose_grams.quantile(0.25):.1f}-{trt.mg_dose_grams.quantile(0.75):.1f}"
    )
    print(
        f"  n_mg_doses:    median={trt.n_mg_doses.median():.0f}, "
        f"range={trt.n_mg_doses.min()}-{trt.n_mg_doses.max()}"
    )

    for lo, hi, label in [
        (0.01, 2, "0.1-2g"),
        (2.01, 4, "2.1-4g"),
        (4.01, 8, "4.1-8g"),
        (8.01, 100, ">8g"),
    ]:
        n = ((trt.mg_dose_grams >= lo) & (trt.mg_dose_grams <= hi)).sum()
        if n > 0:
            aki = trt[
                (trt.mg_dose_grams >= lo) & (trt.mg_dose_grams <= hi)
            ].aki_kdigo1.mean()
            print(f"    {label:>8s}: n={n:4d}, AKI={100*aki:.1f}%")

    cohort.to_csv(os.path.join(RESULTS, "04_mimic_cohort.csv"), index=False)
    print(f"\n  Updated cohort with dose columns")


# =====================================================================
if __name__ == "__main__":
    extract_eicu()
    extract_mimic()
    print(f"\n{'=' * 60}")
    print("DONE. Now re-run:")
    print("  Rscript probe_three_analyses.R")
    print("Dose-response section will now find mg_total_dose column.")
    print("=" * 60)
