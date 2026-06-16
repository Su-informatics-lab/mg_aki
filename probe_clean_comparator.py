#!/usr/bin/env python3
"""
probe_clean_comparator.py
Flag 'untreated' patients who actually received Mg outside the 6h window.
Then call R to re-run analyses: treated vs true-never-treated only.

Run: python probe_clean_comparator.py
"""

import os

import pandas as pd

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


# ── MIMIC: find ALL stays that ever received Mg ─────────────────
print("=" * 60)
print("MIMIC-IV: IDENTIFYING COMPARATOR CONTAMINATION")
print("=" * 60)

cohort_m = pd.read_csv(os.path.join(RESULTS, "04_mimic_cohort.csv"))
stays = set(cohort_m.stay_id)

ie_path = gz(os.path.join(MIMIC_ICU, "inputevents.csv.gz"))
print("Loading inputevents...")
ie = pd.read_csv(
    ie_path, usecols=lambda c: c in {"stay_id", "itemid", "amount"}, low_memory=False
)
mg_items = [222011, 227523]
ie_mg = ie[
    (ie.itemid.isin(mg_items))
    & (ie.stay_id.isin(stays))
    & (ie.amount.notna())
    & (ie.amount > 0)
]

ever_mg_stays = set(ie_mg.stay_id.unique())
print(
    f"  Stays with ANY Mg administration: {len(ever_mg_stays)}/{len(stays)} "
    f"({100*len(ever_mg_stays)/len(stays):.1f}%)"
)

# Classify
cohort_m["comparator_status"] = "clean_untreated"
cohort_m.loc[cohort_m.mg_supplementation == 1, "comparator_status"] = "treated"
cohort_m.loc[
    (cohort_m.mg_supplementation == 0) & (cohort_m.stay_id.isin(ever_mg_stays)),
    "comparator_status",
] = "contaminated"

for s in ["treated", "contaminated", "clean_untreated"]:
    sub = cohort_m[cohort_m.comparator_status == s]
    print(f"  {s:>20s}: n={len(sub):5d}  AKI={100*sub.aki_kdigo1.mean():.1f}%")

cohort_m.to_csv(os.path.join(RESULTS, "04_mimic_cohort.csv"), index=False)

# ── eICU: same check ────────────────────────────────────────────
print(f"\n{'=' * 60}")
print("eICU: IDENTIFYING COMPARATOR CONTAMINATION")
print("=" * 60)

cohort_e = pd.read_csv(os.path.join(RESULTS, "01_analysis_a_cohort.csv"))
pids = set(cohort_e.patientunitstayid)

# Check medication table (NO time filter)
ever_mg_pids_e = set()
for fname in ["medication.csv.gz", "medication.csv"]:
    mpath = os.path.join(EICU_ROOT, fname)
    if os.path.exists(mpath):
        print(f"Loading {fname} (all timepoints)...")
        med = pd.read_csv(mpath, low_memory=False)
        med.columns = med.columns.str.lower()
        med = med[med.patientunitstayid.isin(pids)]
        mg_mask = med.drugname.str.lower().str.contains("magnesium", na=False)
        iv_mask = True
        if "routeadmin" in med.columns:
            iv_mask = (
                med.routeadmin.str.lower().str.contains(
                    "iv|intravenous|inject", na=False
                )
                | med.routeadmin.isna()
            )
        mg_med = med[mg_mask & iv_mask]
        ever_mg_pids_e.update(mg_med.patientunitstayid.unique())
        break

# Also check infusionDrug
for fname in [
    "infusionDrug.csv.gz",
    "infusionDrug.csv",
    "infusiondrug.csv.gz",
    "infusiondrug.csv",
]:
    ipath = os.path.join(EICU_ROOT, fname)
    if os.path.exists(ipath):
        inf = pd.read_csv(ipath, low_memory=False)
        inf.columns = inf.columns.str.lower()
        inf = inf[inf.patientunitstayid.isin(pids)]
        mg_inf = inf[inf.drugname.str.lower().str.contains("magnesium", na=False)]
        ever_mg_pids_e.update(mg_inf.patientunitstayid.unique())
        break

print(
    f"  Patients with ANY Mg: {len(ever_mg_pids_e)}/{len(pids)} "
    f"({100*len(ever_mg_pids_e)/len(pids):.1f}%)"
)

cohort_e["comparator_status"] = "clean_untreated"
cohort_e.loc[cohort_e.mg_supplementation == 1, "comparator_status"] = "treated"
cohort_e.loc[
    (cohort_e.mg_supplementation == 0)
    & (cohort_e.patientunitstayid.isin(ever_mg_pids_e)),
    "comparator_status",
] = "contaminated"

for s in ["treated", "contaminated", "clean_untreated"]:
    sub = cohort_e[cohort_e.comparator_status == s]
    print(f"  {s:>20s}: n={len(sub):5d}  AKI={100*sub.aki_kdigo1.mean():.1f}%")

cohort_e.to_csv(os.path.join(RESULTS, "01_analysis_a_cohort.csv"), index=False)

# ── Summary ─────────────────────────────────────────────────────
print(f"\n{'=' * 60}")
print("SUMMARY — COMPARATOR CONTAMINATION")
print("=" * 60)
for db, coh, id_col in [
    ("eICU", cohort_e, "patientunitstayid"),
    ("MIMIC", cohort_m, "stay_id"),
]:
    trt = len(coh[coh.comparator_status == "treated"])
    contam = len(coh[coh.comparator_status == "contaminated"])
    clean = len(coh[coh.comparator_status == "clean_untreated"])
    total_untrt = contam + clean
    print(
        f"  {db}: treated={trt}, untreated={total_untrt} "
        f"(contaminated={contam} [{100*contam/total_untrt:.0f}%], "
        f"clean={clean} [{100*clean/total_untrt:.0f}%])"
    )

print(f"\nNow run: Rscript probe_clean_comparator.R")
