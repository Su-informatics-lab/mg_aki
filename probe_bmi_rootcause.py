#!/usr/bin/env python3
"""
probe_bmi_rootcause.py — Root cause analysis of eICU BMI outliers

Investigates WHY extreme BMIs occur by examining raw admissionheight
and admissionweight patterns across the full cardiac surgery cohort.

Root cause categories:
  1. Height entered in inches (not cm) — multi-site unit confusion
  2. Height data-entry error (decimal shift, typo)
  3. Weight entered in lbs (not kg)
  4. Weight data-entry error
  5. Both height and weight anomalous

For each category, tests whether heuristic correction (e.g. inches→cm)
recovers plausible BMI, or whether NaN is the only safe option.

Run: python probe_bmi_rootcause.py
"""

import os
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

RESULTS = os.path.expanduser("~/mg_aki/results")
EICU_ROOT = None
for p in [
    "~/mg_aki/eicu-crd-2.0",
    "~/mg_aki/eicu-collaborative-research-database-demo-2.0.1",
]:
    if os.path.isdir(os.path.expanduser(p)):
        EICU_ROOT = os.path.expanduser(p)
        break


def find_csv(root, name):
    for n in [name, name.lower()]:
        for ext in [".csv.gz", ".csv"]:
            p = os.path.join(root, n + ext)
            if os.path.exists(p):
                return p
    return None


print("=" * 65)
print("BMI OUTLIER ROOT CAUSE ANALYSIS")
print("=" * 65)

# ── Load cohort + raw patient table ──────────────────────────────────
cohort = pd.read_csv(os.path.join(RESULTS, "01_analysis_a_cohort.csv"))
pids = set(cohort.patientunitstayid)

pt = pd.read_csv(
    find_csv(EICU_ROOT, "patient"),
    usecols=["patientunitstayid", "admissionheight", "admissionweight"],
    low_memory=False,
)
pt.columns = pt.columns.str.lower()
pt_cohort = pt[pt.patientunitstayid.isin(pids)].copy()
print(f"Cohort: {len(cohort)} patients")
print(f"Raw patient records: {len(pt_cohort)}")


# =====================================================================
# SECTION 1: FULL HEIGHT/WEIGHT DISTRIBUTION
# =====================================================================
print(f"\n{'=' * 65}")
print("1. FULL HEIGHT/WEIGHT DISTRIBUTION IN COHORT")
print("=" * 65)

h = pt_cohort.admissionheight
w = pt_cohort.admissionweight

print(f"\n  Height (admissionheight):")
print(f"    Total records: {len(h)}")
print(f"    Missing/zero:  {(h.isna() | (h <= 0)).sum()}")
h_pos = h[h > 0]
print(f"    Valid (>0):     {len(h_pos)}")
print(f"    Range:          {h_pos.min():.1f} – {h_pos.max():.1f}")

# Height buckets
buckets_h = [
    (0, 50, "Implausible (<50cm / data error)"),
    (50, 85, "Likely INCHES (50–85 → 127–216cm)"),
    (85, 120, "Ambiguous (85–120cm: short adult or wrong unit)"),
    (120, 240, "Plausible CM (120–240)"),
    (240, 999, "Implausible (>240cm / data error)"),
]
print(f"\n    Height buckets:")
for lo, hi, label in buckets_h:
    n = ((h_pos >= lo) & (h_pos < hi)).sum()
    pct = 100 * n / len(h_pos) if len(h_pos) > 0 else 0
    if n > 0:
        examples = h_pos[(h_pos >= lo) & (h_pos < hi)].head(5).tolist()
        print(f"      {label}: {n} ({pct:.2f}%)  e.g. {examples}")

print(f"\n  Weight (admissionweight):")
print(f"    Total records: {len(w)}")
print(f"    Missing/zero:  {(w.isna() | (w <= 0)).sum()}")
w_pos = w[w > 0]
print(f"    Valid (>0):     {len(w_pos)}")
print(f"    Range:          {w_pos.min():.1f} – {w_pos.max():.1f}")

buckets_w = [
    (0, 30, "Implausible (<30kg / data error)"),
    (30, 200, "Plausible KG (30–200)"),
    (200, 350, "Likely LBS or extreme obesity (200–350)"),
    (350, 9999, "Implausible (>350kg / data error)"),
]
print(f"\n    Weight buckets:")
for lo, hi, label in buckets_w:
    n = ((w_pos >= lo) & (w_pos < hi)).sum()
    pct = 100 * n / len(w_pos) if len(w_pos) > 0 else 0
    if n > 0:
        examples = w_pos[(w_pos >= lo) & (w_pos < hi)].head(5).tolist()
        print(f"      {label}: {n} ({pct:.2f}%)  e.g. {examples}")


# =====================================================================
# SECTION 2: PER-OUTLIER ROOT CAUSE CLASSIFICATION
# =====================================================================
print(f"\n{'=' * 65}")
print("2. PER-OUTLIER ROOT CAUSE CLASSIFICATION")
print("=" * 65)

# Merge height/weight with cohort BMI
merged = cohort[["patientunitstayid", "bmi", "mg_supplementation"]].merge(
    pt_cohort, on="patientunitstayid", how="left"
)
# Identify outliers (BMI outside [10, 80])
outliers = merged[
    merged.bmi.notna() & ((merged.bmi < 10) | (merged.bmi > 80))
].copy()

print(f"\n  Total outliers (BMI outside [10,80]): {len(outliers)}")


def classify_outlier(row):
    """Classify root cause and attempt correction."""
    h = row.admissionheight
    w = row.admissionweight
    bmi = row.bmi

    causes = []
    corrections = {}

    # Height analysis
    if h < 50:
        causes.append("height_data_error")
        # Cannot correct — no plausible interpretation
    elif 50 <= h < 85:
        causes.append("height_in_inches")
        h_corrected = h * 2.54
        bmi_corrected = w / (h_corrected / 100) ** 2
        corrections["inches_to_cm"] = {
            "h_corrected": round(h_corrected, 1),
            "bmi_corrected": round(bmi_corrected, 1),
            "plausible": 10 <= bmi_corrected <= 80,
        }
    elif 85 <= h < 120:
        causes.append("height_ambiguous_short")
        # Try inches interpretation
        h_as_inches = h * 2.54
        bmi_if_inches = w / (h_as_inches / 100) ** 2
        # Try as-is (just short)
        corrections["as_inches"] = {
            "h_corrected": round(h_as_inches, 1),
            "bmi_corrected": round(bmi_if_inches, 1),
            "plausible": 10 <= bmi_if_inches <= 80,
        }
    elif h > 240:
        causes.append("height_data_error")

    # Weight analysis
    if 0 < w < 30:
        causes.append("weight_implausibly_low")
    elif w > 200:
        causes.append("weight_possibly_lbs")
        w_corrected = w * 0.4536
        h_for_bmi = h if (120 <= h <= 240) else None
        if h_for_bmi:
            bmi_corrected = w_corrected / (h_for_bmi / 100) ** 2
            corrections["lbs_to_kg"] = {
                "w_corrected": round(w_corrected, 1),
                "bmi_corrected": round(bmi_corrected, 1),
                "plausible": 10 <= bmi_corrected <= 80,
            }
    elif w > 350:
        causes.append("weight_data_error")

    if not causes:
        # Height and weight both in normal ranges but BMI still extreme
        # → must be a combination issue
        causes.append("combination_marginal")

    return causes, corrections


print(f"\n  {'PID':>12s} {'h':>6s} {'w':>7s} {'BMI':>8s} {'trt':>3s}  Root cause → Correction")
print(f"  {'-'*80}")

cause_counts = {}
correctable = 0
uncorrectable = 0

for _, row in outliers.iterrows():
    causes, corrections = classify_outlier(row)
    for c in causes:
        cause_counts[c] = cause_counts.get(c, 0) + 1

    # Best correction
    best = None
    for method, result in corrections.items():
        if result["plausible"]:
            best = (method, result)
            break

    if best:
        correctable += 1
        corr_str = f"{best[0]} → BMI={best[1]['bmi_corrected']:.1f}"
    else:
        uncorrectable += 1
        corr_str = "UNCORRECTABLE → set NaN"

    cause_str = " + ".join(causes)
    print(
        f"  {int(row.patientunitstayid):>12d} "
        f"{row.admissionheight:>6.1f} {row.admissionweight:>7.1f} "
        f"{row.bmi:>8.1f} {int(row.mg_supplementation):>3d}  "
        f"{cause_str} → {corr_str}"
    )

print(f"\n  Root cause summary:")
for cause, n in sorted(cause_counts.items(), key=lambda x: -x[1]):
    print(f"    {cause}: {n}")
print(f"\n  Correctable by unit conversion: {correctable}")
print(f"  Uncorrectable (must NaN): {uncorrectable}")


# =====================================================================
# SECTION 3: HOSPITAL-LEVEL PATTERN CHECK
# =====================================================================
print(f"\n{'=' * 65}")
print("3. HOSPITAL-LEVEL PATTERN CHECK")
print("=" * 65)
print("  (Do outliers cluster in specific hospitals → systematic unit issue?)")

if "hospitalid" in cohort.columns:
    outlier_pids = set(outliers.patientunitstayid)
    outlier_hosp = cohort[cohort.patientunitstayid.isin(outlier_pids)][
        ["patientunitstayid", "hospitalid"]
    ]
    hosp_counts = outlier_hosp.hospitalid.value_counts()
    n_hosp_with_outliers = len(hosp_counts)
    print(f"  Outliers spread across {n_hosp_with_outliers} hospitals")
    if hosp_counts.max() > 1:
        print(f"  Hospitals with >1 outlier:")
        for hid, n in hosp_counts[hosp_counts > 1].items():
            # Check if this hospital has systematic unit issues
            hosp_heights = pt_cohort.merge(
                cohort[cohort.hospitalid == hid][["patientunitstayid"]],
                on="patientunitstayid",
            ).admissionheight
            hosp_h_valid = hosp_heights[hosp_heights > 0]
            median_h = hosp_h_valid.median() if len(hosp_h_valid) > 0 else np.nan
            pct_inches = (
                (hosp_h_valid.between(50, 85)).mean() * 100
                if len(hosp_h_valid) > 0
                else 0
            )
            print(
                f"    Hospital {hid}: {n} outliers, "
                f"all heights median={median_h:.1f}cm, "
                f"in inches range: {pct_inches:.0f}%"
            )
    else:
        print(f"  No hospital has >1 outlier → sporadic data entry errors, not systematic")


# =====================================================================
# SECTION 4: WHAT HAPPENS IN THE ETL — FULL TRACE
# =====================================================================
print(f"\n{'=' * 65}")
print("4. ETL TRACE: HOW OUTLIERS PROPAGATE")
print("=" * 65)

print("""
  Current ETL (01_etl.py):
  
    1. Raw height/weight from eICU patient table
       → No plausibility filter on height or weight
    
    2. BMI = weight / (height/100)^2
       → Computed for ALL patients with height>0 and weight>0
       → No range check on computed BMI
    
    3. BMI enters PS model via MICE imputation
       → MICE uses PMM (predictive mean matching)
       → Extreme observed BMIs become potential donor values
       → But PMM matches on predicted value, not raw BMI
       → So extreme donors rarely match to typical patients
    
    4. PS model: glm(mg_supp ~ ... + bmi + ...)
       → Extreme BMI → extreme linear predictor → PS near 0 or 1
       → Warning: "fitted probabilities numerically 0 or 1 occurred"
    
    5. Overlap weights: OW = ps if untreated, (1-ps) if treated
       → PS near 0 or 1 → OW near 0
       → Outlier patients effectively removed from weighted analysis
    
    6. Net effect: |ΔOR| = 0.0009
       → OW self-corrected, but the warning is a red flag
       → A reviewer seeing that warning could question PS model validity
""")


# =====================================================================
# SECTION 5: RECOMMENDED FIX (MULTI-LAYER)
# =====================================================================
print("=" * 65)
print("5. RECOMMENDED ETL FIX")
print("=" * 65)

print("""
  OPTION A: Simple cap (minimal, sufficient for results)
  ─────────────────────────────────────────────────────────
    cohort.loc[cohort.bmi.notna() & ~cohort.bmi.between(10, 80), "bmi"] = np.nan
    
    Pro: One line, catches all 22 outliers
    Con: Discards potentially recoverable data
    Impact: |ΔOR| = 0.0009 — zero practical effect
  
  OPTION B: Height/weight plausibility + heuristic correction
  ─────────────────────────────────────────────────────────
    # Step 1: Height plausibility
    # Heights 50-85 with plausible weight → likely inches, convert
    inches_mask = (
        cohort.admissionheight.between(50, 85) &
        cohort.admissionweight.between(30, 200)
    )
    cohort.loc[inches_mask, "admissionheight"] *= 2.54
    
    # Step 2: Height range filter
    cohort.loc[
        ~cohort.admissionheight.between(120, 240), "admissionheight"
    ] = np.nan
    
    # Step 3: Weight range filter  
    cohort.loc[
        ~cohort.admissionweight.between(30, 300), "admissionweight"
    ] = np.nan
    
    # Step 4: Recompute BMI from cleaned height/weight
    cohort["bmi"] = np.where(
        cohort.admissionheight.notna() & cohort.admissionweight.notna(),
        cohort.admissionweight / (cohort.admissionheight / 100) ** 2,
        np.nan,
    )
    
    # Step 5: Safety cap
    cohort.loc[cohort.bmi.notna() & ~cohort.bmi.between(10, 80), "bmi"] = np.nan
    
    Pro: Recovers ~6 patients with inches-as-cm error
    Con: Heuristic correction adds complexity, marginal benefit
    Impact: Same — |ΔOR| unchanged either way
  
  RECOMMENDATION: OPTION A
  ─────────────────────────────────────────────────────────
  The 22 outliers are 0.27% of the cohort. Even the "correctable" ones
  (height in inches) would only recover ~6 patients' BMI. With |ΔOR|
  = 0.0009, neither option changes any result. Option A is auditable,
  defensible, and mirrors the MIMIC ETL's existing approach.
  
  The one-line fix aligns both databases:
    eICU:  cohort.loc[~cohort.bmi.between(10, 80), "bmi"] = np.nan  (NEW)
    MIMIC: cohort.loc[~cohort.bmi.between(10, 80), "bmi"] = np.nan  (EXISTING)
""")

# ── Verify: how many patients would each option recover? ─────────────
print("\n  Recovery comparison:")
# Option A: all 22 → NaN
# Option B: attempt corrections
n_recovered = 0
for _, row in outliers.iterrows():
    causes, corrections = classify_outlier(row)
    for method, result in corrections.items():
        if result["plausible"]:
            n_recovered += 1
            break

print(f"    Option A (cap only): 0 recovered, 22 → NaN")
print(f"    Option B (heuristic): {n_recovered} potentially recovered, {len(outliers) - n_recovered} → NaN")
print(f"    Net difference: {n_recovered} patients out of 8109 (0.{n_recovered:02d}%)")

print(f"\n{'=' * 65}")
print("DONE")
print("=" * 65)
