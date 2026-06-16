#!/usr/bin/env python3
"""
01b_uo_kdigo.py — Add KDIGO urine-output AKI staging to cohort CSVs

  Reads existing cohort CSVs + raw UO data from both databases.
  Computes KDIGO UO staging following MIT-LCP/mimic-code canonical logic
  (kdigo_uo.sql → kdigo_stages.sql).

  KDIGO UO criteria:
    Stage 1: <0.5 mL/kg/h for ≥6 consecutive hours
    Stage 2: <0.5 mL/kg/h for ≥12 consecutive hours
    Stage 3: <0.3 mL/kg/h for ≥24h OR anuria (0 mL) for ≥12h
  Combined: max(creatinine_stage, uo_stage)

  New columns added to cohort CSVs:
    n_uo_measures     — number of UO measurements in follow-up window
    uo_total_ml       — total UO in follow-up (mL)
    weight_kg         — patient weight used for mL/kg/h
    aki_uo            — 1 if KDIGO UO stage ≥1 during follow-up
    aki_uo_stage      — max KDIGO UO stage (0-3)
    aki_combined      — 1 if aki_kdigo1 (Cr) OR aki_uo

Prereqs: 01_etl.py must have run (cohort CSVs exist)
Outputs:
    results/01_analysis_a_cohort.csv  (updated in place)
    results/04_mimic_cohort.csv       (updated in place)
    results/probe_uo_diagnostics.csv  (availability/rate summary)

Run: python 01b_uo_kdigo.py
Time: ~5 minutes (large outputevents/intakeOutput reads)
"""

import os
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ── Paths ────────────────────────────────────────────────────────────
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

# ── Constants ────────────────────────────────────────────────────────
# MIMIC-IV UO item IDs (from mimic-code/measurement/urine_output.sql)
MIMIC_UO_ITEMS = [
    226559,  # Foley
    226560,  # Void
    226561,  # Condom Cath
    226563,  # Suprapubic
    226564,  # R Nephrostomy
    226565,  # Straight Cath
    226567,  # L Nephrostomy
    226584,  # Ileoconduit
    227489,  # GU Irrigant/Urine Volume Out
]

# MIMIC weight item in chartevents
WEIGHT_ITEM = 226512

# KDIGO thresholds
KDIGO_STAGE1_RATE = 0.5  # mL/kg/h
KDIGO_STAGE2_RATE = 0.5  # same rate, longer duration
KDIGO_STAGE3_RATE = 0.3
KDIGO_ANURIA_RATE = 0.0  # literal zero for anuria

# Follow-up window (hours), matching Cr-based AKI
FOLLOWUP_HOURS = 168  # 7 days


def gz(path):
    """Check .csv.gz then .csv."""
    if os.path.exists(path):
        return path
    alt = path.replace(".csv.gz", ".csv") if path.endswith(".gz") else path + ".gz"
    return alt if os.path.exists(alt) else path


# =====================================================================
# KDIGO UO STAGING ENGINE (shared by both databases)
# =====================================================================
def compute_kdigo_uo(uo_df, weight_kg, mg_offset_h, max_h=FOLLOWUP_HOURS):
    """
    Compute KDIGO UO staging for a single patient.

    Args:
        uo_df: DataFrame with columns [offset_h, value_ml]
               offset_h = hours from ICU admission
        weight_kg: patient weight in kg
        mg_offset_h: offset of Mg measurement (hours from admission)
                     — UO events before this are excluded (same as Cr)
        max_h: maximum follow-up hours

    Returns:
        dict with n_uo, uo_total, aki_uo_stage
    """
    if weight_kg is None or np.isnan(weight_kg) or weight_kg <= 10:
        return {"n_uo": 0, "uo_total": 0, "aki_uo_stage": np.nan}

    # Filter to follow-up window (after Mg measurement, up to max_h)
    fu = uo_df[
        (uo_df.offset_h > mg_offset_h) & (uo_df.offset_h <= mg_offset_h + max_h)
    ].copy()

    if len(fu) == 0:
        return {"n_uo": 0, "uo_total": 0, "aki_uo_stage": np.nan}

    n_uo = len(fu)
    uo_total = fu.value_ml.sum()

    # ── Hourly binning ───────────────────────────────────────────
    # Bin relative to Mg measurement
    fu["rel_h"] = fu.offset_h - mg_offset_h
    fu["hour_bin"] = np.floor(fu.rel_h).astype(int)

    # Sum UO per hour
    hourly = fu.groupby("hour_bin").value_ml.sum()

    # Create complete hourly grid (missing hours = 0 UO if patient
    # is in ICU; this is conservative — assumes nurse would have
    # charted if patient was producing urine)
    max_hour = min(int(hourly.index.max()) + 1, int(max_h))
    min_hour = max(0, int(hourly.index.min()))
    full_idx = range(min_hour, max_hour)
    hourly = hourly.reindex(full_idx, fill_value=0.0)

    if len(hourly) < 6:
        # Need at least 6 hours for Stage 1
        return {"n_uo": n_uo, "uo_total": uo_total, "aki_uo_stage": 0}

    # ── Rolling rates ────────────────────────────────────────────
    max_stage = 0

    # Stage 1: <0.5 mL/kg/h for ≥6h
    if len(hourly) >= 6:
        roll_6 = hourly.rolling(6, min_periods=6).sum()
        rate_6 = roll_6 / (weight_kg * 6)
        if (rate_6.dropna() < KDIGO_STAGE1_RATE).any():
            max_stage = max(max_stage, 1)

    # Stage 2: <0.5 mL/kg/h for ≥12h
    if len(hourly) >= 12:
        roll_12 = hourly.rolling(12, min_periods=12).sum()
        rate_12 = roll_12 / (weight_kg * 12)
        if (rate_12.dropna() < KDIGO_STAGE2_RATE).any():
            max_stage = max(max_stage, 2)

    # Stage 3: <0.3 mL/kg/h for ≥24h
    if len(hourly) >= 24:
        roll_24 = hourly.rolling(24, min_periods=24).sum()
        rate_24 = roll_24 / (weight_kg * 24)
        if (rate_24.dropna() < KDIGO_STAGE3_RATE).any():
            max_stage = max(max_stage, 3)

    # Stage 3 (alt): anuria for ≥12h
    if len(hourly) >= 12:
        roll_12_total = hourly.rolling(12, min_periods=12).sum()
        if (roll_12_total.dropna() == 0).any():
            max_stage = max(max_stage, 3)

    return {"n_uo": n_uo, "uo_total": round(uo_total, 1), "aki_uo_stage": max_stage}


# =====================================================================
# eICU UO EXTRACTION
# =====================================================================
def run_eicu():
    print(f"\n{'=' * 65}")
    print("eICU: URINE OUTPUT EXTRACTION + KDIGO STAGING")
    print("=" * 65)

    # ── Load cohort ──────────────────────────────────────────────
    cohort_path = os.path.join(RESULTS, "01_analysis_a_cohort.csv")
    if not os.path.exists(cohort_path):
        print("  ERROR: eICU cohort not found")
        return
    cohort = pd.read_csv(cohort_path)
    pids = set(cohort.patientunitstayid)
    print(f"  Cohort: {len(cohort)} patients")

    # ── Get weight ───────────────────────────────────────────────
    # admissionweight should be in the cohort from the patient table
    if "admissionweight" in cohort.columns:
        print(f"  Weight from cohort.admissionweight")
        wt = cohort.set_index("patientunitstayid")["admissionweight"]
    else:
        print("  Loading weight from patient table...")
        pt_path = None
        for n in ["patient.csv.gz", "patient.csv"]:
            p = os.path.join(EICU_ROOT, n)
            if os.path.exists(p):
                pt_path = p
                break
        if pt_path is None:
            print("  ERROR: patient table not found")
            return
        pt = pd.read_csv(
            pt_path, usecols=["patientunitstayid", "admissionweight"], low_memory=False
        )
        pt.columns = pt.columns.str.lower()
        pt = pt[pt.patientunitstayid.isin(pids)]
        wt = pt.set_index("patientunitstayid")["admissionweight"]

    n_wt = (wt > 10).sum()
    print(f"  Valid weight (>10 kg): {n_wt}/{len(pids)} ({100*n_wt/len(pids):.1f}%)")

    # ── Load intakeOutput ────────────────────────────────────────
    if not EICU_ROOT:
        print("  ERROR: eICU data root not found")
        return

    io_path = None
    for n in [
        "intakeOutput.csv.gz",
        "intakeOutput.csv",
        "intakeoutput.csv.gz",
        "intakeoutput.csv",
    ]:
        p = os.path.join(EICU_ROOT, n)
        if os.path.exists(p):
            io_path = p
            break
    if io_path is None:
        print("  ERROR: intakeOutput not found")
        return

    print(f"  Loading intakeOutput from {os.path.basename(io_path)}...")
    io = pd.read_csv(io_path, low_memory=False)
    io.columns = io.columns.str.lower()
    print(f"  intakeOutput: {len(io):,} total rows")

    # Filter to cohort
    io = io[io.patientunitstayid.isin(pids)]
    print(f"  In cohort: {len(io):,} rows")

    # ── Filter to urine output ───────────────────────────────────
    # Use celllabel or cellpath containing urine-related terms
    # (following MIT-LCP/eicu-code/concepts/pivoted/pivoted-uo.sql)
    uo_mask = pd.Series(False, index=io.index)
    for col in ["celllabel", "cellpath"]:
        if col in io.columns:
            uo_mask |= (
                io[col]
                .str.lower()
                .str.contains(
                    "urine|foley|void|condom cath|suprapubic|nephrostomy|straight cath",
                    na=False,
                )
            )

    uo = io[uo_mask].copy()
    print(f"  UO rows: {len(uo):,}")

    # Get the value column
    val_col = None
    for c in ["cellvaluenumeric", "cellvaluetext"]:
        if c in uo.columns:
            val_col = c
            break
    if val_col is None:
        print("  ERROR: no value column found")
        return

    uo["value_ml"] = pd.to_numeric(uo[val_col], errors="coerce")
    uo = uo[uo.value_ml.notna() & (uo.value_ml >= 0)]
    uo["offset_h"] = uo["intakeoutputoffset"] / 60.0  # minutes → hours

    n_pts_uo = uo.patientunitstayid.nunique()
    print(
        f"  Patients with any UO data: {n_pts_uo}/{len(pids)} ({100*n_pts_uo/len(pids):.1f}%)"
    )

    # ── Compute KDIGO UO per patient ─────────────────────────────
    print("  Computing KDIGO UO staging...")
    mg_off_col = "mg_offset" if "mg_offset" in cohort.columns else None
    results = []

    for _, row in cohort.iterrows():
        pid = row.patientunitstayid
        mg_off_h = (
            (row[mg_off_col] / 60.0)
            if mg_off_col and not np.isnan(row[mg_off_col])
            else 0
        )

        pt_uo = uo[uo.patientunitstayid == pid][["offset_h", "value_ml"]]
        pt_wt = wt.get(pid, np.nan)
        if np.isnan(pt_wt) or pt_wt <= 10:
            pt_wt = np.nan

        if len(pt_uo) == 0 or np.isnan(pt_wt):
            results.append(
                {
                    "patientunitstayid": pid,
                    "n_uo_measures": len(pt_uo) if len(pt_uo) > 0 else 0,
                    "uo_total_ml": 0,
                    "weight_kg": pt_wt if not np.isnan(pt_wt) else np.nan,
                    "aki_uo_stage": np.nan,
                }
            )
        else:
            r = compute_kdigo_uo(pt_uo, pt_wt, mg_off_h)
            r["patientunitstayid"] = pid
            r["weight_kg"] = pt_wt
            results.append(r)

    res_df = pd.DataFrame(results)
    res_df["aki_uo"] = (res_df.aki_uo_stage >= 1).astype(int)
    res_df.loc[res_df.aki_uo_stage.isna(), "aki_uo"] = np.nan

    # ── Merge into cohort ────────────────────────────────────────
    # Drop old UO columns if they exist
    for c in [
        "n_uo_measures",
        "uo_total_ml",
        "weight_kg",
        "aki_uo_stage",
        "aki_uo",
        "aki_combined",
    ]:
        if c in cohort.columns:
            cohort = cohort.drop(columns=[c])

    cohort = cohort.merge(
        res_df[
            [
                "patientunitstayid",
                "n_uo_measures",
                "uo_total_ml",
                "weight_kg",
                "aki_uo_stage",
                "aki_uo",
            ]
        ],
        on="patientunitstayid",
        how="left",
    )
    cohort["aki_combined"] = ((cohort.aki_kdigo1 == 1) | (cohort.aki_uo == 1)).astype(
        int
    )
    # If UO data missing, combined = Cr only
    cohort.loc[cohort.aki_uo.isna(), "aki_combined"] = cohort.loc[
        cohort.aki_uo.isna(), "aki_kdigo1"
    ]

    # ── Diagnostics ──────────────────────────────────────────────
    n_uo_avail = (cohort.aki_uo.notna()).sum()
    n_aki_cr = cohort.aki_kdigo1.sum()
    n_aki_uo = cohort.aki_uo.sum() if cohort.aki_uo.notna().any() else 0
    n_aki_comb = cohort.aki_combined.sum()

    print(f"\n  ── eICU RESULTS ──")
    print(
        f"  UO data available for KDIGO staging: {n_uo_avail}/{len(cohort)} ({100*n_uo_avail/len(cohort):.1f}%)"
    )
    print(f"  AKI by Cr only:    {n_aki_cr} ({100*n_aki_cr/len(cohort):.1f}%)")
    print(
        f"  AKI by UO only:    {n_aki_uo} ({100*n_aki_uo/max(n_uo_avail,1):.1f}% of those with UO data)"
    )
    print(f"  AKI combined:      {n_aki_comb} ({100*n_aki_comb/len(cohort):.1f}%)")

    # UO staging distribution
    if n_uo_avail > 0:
        staged = cohort[cohort.aki_uo.notna()]
        for s in [0, 1, 2, 3]:
            n = (staged.aki_uo_stage == s).sum()
            print(f"    UO Stage {s}: {n} ({100*n/len(staged):.1f}%)")

    # By treatment
    trt_col = "mg_supplementation"
    for g in [1, 0]:
        sub = cohort[(cohort[trt_col] == g) & cohort.aki_uo.notna()]
        if len(sub) > 0:
            label = "Trt" if g == 1 else "Ctrl"
            print(
                f"  {label}: AKI-Cr={sub.aki_kdigo1.mean()*100:.1f}%, "
                f"AKI-UO={sub.aki_uo.mean()*100:.1f}%, "
                f"AKI-combined={sub.aki_combined.mean()*100:.1f}% "
                f"(n={len(sub)})"
            )

    # ── Save ─────────────────────────────────────────────────────
    cohort.to_csv(cohort_path, index=False)
    print(f"\n  Updated: {cohort_path}")

    return {
        "db": "eICU",
        "n": len(cohort),
        "n_uo_avail": n_uo_avail,
        "aki_cr": n_aki_cr,
        "aki_uo": n_aki_uo,
        "aki_combined": n_aki_comb,
    }


# =====================================================================
# MIMIC-IV UO EXTRACTION
# =====================================================================
def run_mimic():
    print(f"\n{'=' * 65}")
    print("MIMIC-IV: URINE OUTPUT EXTRACTION + KDIGO STAGING")
    print("=" * 65)

    # ── Load cohort ──────────────────────────────────────────────
    cohort_path = os.path.join(RESULTS, "04_mimic_cohort.csv")
    if not os.path.exists(cohort_path):
        print("  ERROR: MIMIC cohort not found")
        return
    cohort = pd.read_csv(cohort_path)
    stays = set(cohort.stay_id)
    print(f"  Cohort: {len(cohort)} patients")

    # Parse datetimes
    cohort["intime"] = pd.to_datetime(cohort["intime"])
    if "mg_charttime" in cohort.columns:
        cohort["mg_charttime"] = pd.to_datetime(cohort["mg_charttime"])

    # ── Get weight from chartevents ──────────────────────────────
    print("  Loading weight from chartevents...")
    ce_path = gz(os.path.join(MIMIC_ICU, "chartevents.csv.gz"))
    # chartevents is huge; read in chunks filtered to weight + our stays
    wt_dict = {}
    chunk_iter = pd.read_csv(
        ce_path,
        usecols=["stay_id", "itemid", "charttime", "valuenum"],
        chunksize=2_000_000,
        low_memory=False,
    )
    for chunk in chunk_iter:
        chunk = chunk[
            (chunk.itemid == WEIGHT_ITEM)
            & chunk.stay_id.isin(stays)
            & chunk.valuenum.between(20, 300)  # plausible kg
        ]
        for _, r in chunk.iterrows():
            sid = r.stay_id
            if sid not in wt_dict:
                wt_dict[sid] = r.valuenum  # take first valid weight

    n_wt = len(wt_dict)
    print(f"  Weight available: {n_wt}/{len(stays)} ({100*n_wt/len(stays):.1f}%)")

    # Fallback: compute weight from BMI + height if available
    if "bmi" in cohort.columns:
        for _, row in cohort.iterrows():
            sid = row.stay_id
            if sid not in wt_dict and not np.isnan(row.get("bmi", np.nan)):
                # Can't reliably back-compute weight from BMI without height
                pass

    # ── Load outputevents ────────────────────────────────────────
    oe_path = gz(os.path.join(MIMIC_ICU, "outputevents.csv.gz"))
    if not os.path.exists(oe_path):
        print(f"  ERROR: outputevents not found at {oe_path}")
        return

    print(f"  Loading outputevents...")
    oe = pd.read_csv(
        oe_path,
        usecols=["stay_id", "charttime", "itemid", "value"],
        low_memory=False,
    )
    print(f"  outputevents: {len(oe):,} total rows")

    # Filter to cohort + UO items
    oe = oe[oe.stay_id.isin(stays) & oe.itemid.isin(MIMIC_UO_ITEMS)]
    oe = oe[oe.value.notna() & (oe.value >= 0)]
    oe["charttime"] = pd.to_datetime(oe["charttime"])
    print(f"  UO rows in cohort: {len(oe):,}")

    n_pts_uo = oe.stay_id.nunique()
    print(
        f"  Patients with any UO data: {n_pts_uo}/{len(stays)} ({100*n_pts_uo/len(stays):.1f}%)"
    )

    # Compute offset_h relative to ICU admission
    oe = oe.merge(cohort[["stay_id", "intime"]], on="stay_id")
    oe["offset_h"] = (oe.charttime - oe.intime).dt.total_seconds() / 3600
    oe = oe.rename(columns={"value": "value_ml"})

    # ── Compute KDIGO UO per patient ─────────────────────────────
    print("  Computing KDIGO UO staging...")
    results = []

    for _, row in cohort.iterrows():
        sid = row.stay_id
        # Mg measurement offset (hours from admission)
        if "mg_charttime" in row.index and pd.notna(row.mg_charttime):
            mg_off_h = (row.mg_charttime - row.intime).total_seconds() / 3600
        else:
            mg_off_h = 0

        pt_uo = oe[oe.stay_id == sid][["offset_h", "value_ml"]]
        pt_wt = wt_dict.get(sid, np.nan)

        if len(pt_uo) == 0 or np.isnan(pt_wt):
            results.append(
                {
                    "stay_id": sid,
                    "n_uo_measures": len(pt_uo),
                    "uo_total_ml": 0,
                    "weight_kg": pt_wt if not np.isnan(pt_wt) else np.nan,
                    "aki_uo_stage": np.nan,
                }
            )
        else:
            r = compute_kdigo_uo(pt_uo, pt_wt, mg_off_h)
            r["stay_id"] = sid
            r["weight_kg"] = pt_wt
            results.append(r)

    res_df = pd.DataFrame(results)
    res_df["aki_uo"] = (res_df.aki_uo_stage >= 1).astype(int)
    res_df.loc[res_df.aki_uo_stage.isna(), "aki_uo"] = np.nan

    # ── Merge into cohort ────────────────────────────────────────
    for c in [
        "n_uo_measures",
        "uo_total_ml",
        "weight_kg",
        "aki_uo_stage",
        "aki_uo",
        "aki_combined",
    ]:
        if c in cohort.columns:
            cohort = cohort.drop(columns=[c])

    cohort = cohort.merge(
        res_df[
            [
                "stay_id",
                "n_uo_measures",
                "uo_total_ml",
                "weight_kg",
                "aki_uo_stage",
                "aki_uo",
            ]
        ],
        on="stay_id",
        how="left",
    )
    cohort["aki_combined"] = ((cohort.aki_kdigo1 == 1) | (cohort.aki_uo == 1)).astype(
        int
    )
    cohort.loc[cohort.aki_uo.isna(), "aki_combined"] = cohort.loc[
        cohort.aki_uo.isna(), "aki_kdigo1"
    ]

    # ── Diagnostics ──────────────────────────────────────────────
    n_uo_avail = (cohort.aki_uo.notna()).sum()
    n_aki_cr = cohort.aki_kdigo1.sum()
    n_aki_uo = cohort.aki_uo.sum() if cohort.aki_uo.notna().any() else 0
    n_aki_comb = cohort.aki_combined.sum()

    print(f"\n  ── MIMIC RESULTS ──")
    print(
        f"  UO data available for KDIGO staging: {n_uo_avail}/{len(cohort)} ({100*n_uo_avail/len(cohort):.1f}%)"
    )
    print(f"  AKI by Cr only:    {n_aki_cr} ({100*n_aki_cr/len(cohort):.1f}%)")
    print(
        f"  AKI by UO only:    {n_aki_uo} ({100*n_aki_uo/max(n_uo_avail,1):.1f}% of those with UO data)"
    )
    print(f"  AKI combined:      {n_aki_comb} ({100*n_aki_comb/len(cohort):.1f}%)")

    if n_uo_avail > 0:
        staged = cohort[cohort.aki_uo.notna()]
        for s in [0, 1, 2, 3]:
            n = (staged.aki_uo_stage == s).sum()
            print(f"    UO Stage {s}: {n} ({100*n/len(staged):.1f}%)")

    trt_col = "mg_supplementation"
    for g in [1, 0]:
        sub = cohort[(cohort[trt_col] == g) & cohort.aki_uo.notna()]
        if len(sub) > 0:
            label = "Trt" if g == 1 else "Ctrl"
            print(
                f"  {label}: AKI-Cr={sub.aki_kdigo1.mean()*100:.1f}%, "
                f"AKI-UO={sub.aki_uo.mean()*100:.1f}%, "
                f"AKI-combined={sub.aki_combined.mean()*100:.1f}% "
                f"(n={len(sub)})"
            )

    # Also check MIMIC dose unit while we're here
    print("\n  ── MIMIC Mg DOSE UNIT CHECK ──")
    if "mg_total_dose" in cohort.columns:
        ie_path = gz(os.path.join(MIMIC_ICU, "inputevents.csv.gz"))
        try:
            ie = pd.read_csv(
                ie_path,
                usecols=["stay_id", "itemid", "amount", "amountuom"],
                low_memory=False,
            )
            mg_ie = ie[ie.itemid.isin([222011, 227523])]
            print(f"  amountuom distribution for Mg items:")
            for unit, n in mg_ie.amountuom.value_counts().items():
                sub_amt = mg_ie[mg_ie.amountuom == unit].amount
                print(
                    f"    {unit}: n={n}, median={sub_amt.median():.1f}, "
                    f"mean={sub_amt.mean():.1f}, range={sub_amt.min():.1f}-{sub_amt.max():.1f}"
                )
        except Exception as e:
            print(f"  Could not check dose unit: {e}")

    # ── Save ─────────────────────────────────────────────────────
    cohort.to_csv(cohort_path, index=False)
    print(f"\n  Updated: {cohort_path}")

    return {
        "db": "MIMIC",
        "n": len(cohort),
        "n_uo_avail": n_uo_avail,
        "aki_cr": n_aki_cr,
        "aki_uo": n_aki_uo,
        "aki_combined": n_aki_comb,
    }


# =====================================================================
# MAIN
# =====================================================================
if __name__ == "__main__":
    print("=" * 65)
    print("01b_uo_kdigo.py — KDIGO Urine Output AKI Staging")
    print("  Reference: MIT-LCP/mimic-code kdigo_stages.sql")
    print("=" * 65)

    diag = []
    args = [a.lower() for a in sys.argv[1:]]
    run_both = len(args) == 0

    if run_both or "eicu" in args:
        r = run_eicu()
        if r:
            diag.append(r)

    if run_both or "mimic" in args:
        r = run_mimic()
        if r:
            diag.append(r)

    # ── Save diagnostics ─────────────────────────────────────────
    if diag:
        diag_df = pd.DataFrame(diag)
        diag_path = os.path.join(RESULTS, "probe_uo_diagnostics.csv")
        diag_df.to_csv(diag_path, index=False)
        print(f"\n  Diagnostics: {diag_path}")

    # ── Summary ──────────────────────────────────────────────────
    print(f"\n{'=' * 65}")
    print("SUMMARY")
    print("=" * 65)
    for d in diag:
        print(
            f"  {d['db']}: UO available {d['n_uo_avail']}/{d['n']} "
            f"({100*d['n_uo_avail']/d['n']:.0f}%)"
        )
        print(
            f"    AKI Cr-only: {d['aki_cr']} → Combined: {d['aki_combined']} "
            f"(+{d['aki_combined']-d['aki_cr']} from UO)"
        )
    print(f"\n  Next steps:")
    print(f"  1. If UO availability >50%: run primary analyses with aki_combined")
    print(
        f"  2. Rscript 02_analysis.R  (it reads aki_kdigo1; add aki_combined sensitivity)"
    )
    print(f"  3. Compare Cr-only vs Combined results in supplement eTable")
    print(f"\n{'=' * 65}")
    print("✓ UO-KDIGO staging complete.")
