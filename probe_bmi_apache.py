#!/usr/bin/env python3
"""
probe_bmi_apache.py — Pre-submission QC probes (standalone, reads only)

  Section A: BMI distribution audit + outlier quantification
  Section B: APACHE IV score extraction + availability check
  Section C: Achieved-Mg distribution (baseline ≠ achieved quantification)

Outputs:
  results/probe_bmi_audit.csv         — per-patient BMI audit
  results/probe_cohort_bmifixed.csv   — eICU cohort with capped BMI
  results/probe_apache_scores.csv     — APACHE IV merged with cohort
  results/probe_achieved_mg.csv       — follow-up Mg distributions

Run on Tempest: python probe_bmi_apache.py
Does NOT modify any pipeline files.
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
MIMIC_ROOT = os.path.expanduser("~/mg_aki/mimic-iv-3.1")


def find_csv(root, name):
    for n in [name, name.lower()]:
        for ext in [".csv.gz", ".csv"]:
            p = os.path.join(root, n + ext)
            if os.path.exists(p):
                return p
    return None


# =====================================================================
# SECTION A: BMI DISTRIBUTION AUDIT
# =====================================================================
def section_a_bmi_audit():
    print("=" * 65)
    print("A. BMI DISTRIBUTION AUDIT")
    print("=" * 65)

    # ── Load cohort ──────────────────────────────────────────────
    cohort_path = os.path.join(RESULTS, "01_analysis_a_cohort.csv")
    if not os.path.exists(cohort_path):
        print("  ERROR: cohort not found")
        return
    cohort = pd.read_csv(cohort_path)
    print(f"  eICU cohort: {len(cohort)} patients")

    # ── BMI distribution ─────────────────────────────────────────
    has_bmi = cohort.bmi.notna()
    n_bmi = has_bmi.sum()
    pct_bmi = 100 * n_bmi / len(cohort)
    print(f"\n  BMI available: {n_bmi} ({pct_bmi:.1f}%)")

    bmi = cohort.bmi[has_bmi]
    print(f"  BMI distribution:")
    print(f"    Mean:   {bmi.mean():.1f}")
    print(f"    Median: {bmi.median():.1f}")
    print(f"    SD:     {bmi.std():.1f}")
    print(f"    Min:    {bmi.min():.1f}")
    print(f"    Max:    {bmi.max():.1f}")
    for q in [0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99]:
        print(f"    P{int(q*100):02d}:    {bmi.quantile(q):.1f}")

    # ── Outlier counts at thresholds ─────────────────────────────
    print(f"\n  Outlier counts:")
    thresholds = [
        (10, "below", lambda x: x < 10),
        (15, "below", lambda x: x < 15),
        (60, "above", lambda x: x > 60),
        (80, "above", lambda x: x > 80),
        (100, "above", lambda x: x > 100),
        (200, "above", lambda x: x > 200),
        (500, "above", lambda x: x > 500),
    ]
    for thresh, direction, fn in thresholds:
        n = fn(bmi).sum()
        if n > 0:
            print(f"    BMI {direction} {thresh}: {n} patients")
    n_outside_10_80 = ((bmi < 10) | (bmi > 80)).sum()
    print(f"\n  Total outside [10, 80]: {n_outside_10_80} patients")

    # ── Investigate outlier patients ─────────────────────────────
    outliers = cohort[has_bmi & ((cohort.bmi < 10) | (cohort.bmi > 80))].copy()
    if len(outliers) > 0:
        print(f"\n  ── Outlier patients (n={len(outliers)}) ──")
        trt_col = (
            "mg_supplementation"
            if "mg_supplementation" in outliers.columns
            else "mg_supp"
        )
        print(f"    Treatment distribution:")
        print(f"      Supplemented:   {(outliers[trt_col] == 1).sum()}")
        print(f"      Unsupplemented: {(outliers[trt_col] == 0).sum()}")
        if "ac_group" in outliers.columns:
            print(f"    AC group distribution:")
            for g in ["mg_k", "k_only", "mg_only", "neither"]:
                n = (outliers.ac_group == g).sum()
                if n > 0:
                    print(f"      {g}: {n}")
        # Show raw height/weight
        hw_cols = ["admissionheight", "admissionweight"]
        avail_hw = [c for c in hw_cols if c in outliers.columns]
        if avail_hw:
            print(f"\n    Raw height/weight for outliers:")
            for _, row in outliers[avail_hw + ["bmi", trt_col]].head(20).iterrows():
                h = row.get("admissionheight", np.nan)
                w = row.get("admissionweight", np.nan)
                print(
                    f"      h={h:.1f}cm  w={w:.1f}kg  → BMI={row['bmi']:.1f}  trt={int(row[trt_col])}"
                )

    # ── Investigate raw data source ──────────────────────────────
    print(f"\n  ── Raw eICU patient table audit ──")
    if EICU_ROOT:
        pt_path = find_csv(EICU_ROOT, "patient")
        if pt_path:
            pt = pd.read_csv(
                pt_path,
                usecols=["patientunitstayid", "admissionheight", "admissionweight"],
                low_memory=False,
            )
            pt.columns = pt.columns.str.lower()
            pt_cohort = pt[pt.patientunitstayid.isin(set(cohort.patientunitstayid))]
            h = pt_cohort.admissionheight[pt_cohort.admissionheight > 0]
            w = pt_cohort.admissionweight[pt_cohort.admissionweight > 0]
            print(f"    Height (>0): n={len(h)}")
            print(f"      Range: {h.min():.1f} – {h.max():.1f}")
            print(f"      Likely inches (40–85): {((h >= 40) & (h <= 85)).sum()}")
            print(f"      Likely cm (100–220): {((h >= 100) & (h <= 220)).sum()}")
            print(f"      Suspicious (<100 or >220): {((h < 100) | (h > 220)).sum()}")
            if ((h < 100) | (h > 220)).sum() > 0:
                susp = pt_cohort[
                    (pt_cohort.admissionheight > 0)
                    & (
                        (pt_cohort.admissionheight < 100)
                        | (pt_cohort.admissionheight > 220)
                    )
                ]
                print(f"      Examples:")
                for _, row in susp.head(10).iterrows():
                    bmi_calc = (
                        row.admissionweight / (row.admissionheight / 100) ** 2
                        if row.admissionheight > 0 and row.admissionweight > 0
                        else np.nan
                    )
                    print(
                        f"        h={row.admissionheight:.1f} w={row.admissionweight:.1f} → BMI={bmi_calc:.1f}"
                    )
            print(f"\n    Weight (>0): n={len(w)}")
            print(f"      Range: {w.min():.1f} – {w.max():.1f}")
            print(f"      Likely lbs (>200): {(w > 200).sum()}")
        else:
            print("    patient table not found")
    else:
        print("    eICU data root not found")

    # ── Create fixed cohort ──────────────────────────────────────
    print(f"\n  ── Creating BMI-fixed cohort ──")
    cohort_fix = cohort.copy()
    n_capped = ((cohort_fix.bmi < 10) | (cohort_fix.bmi > 80)).sum()
    cohort_fix.loc[
        cohort_fix.bmi.notna() & ((cohort_fix.bmi < 10) | (cohort_fix.bmi > 80)), "bmi"
    ] = np.nan
    print(f"    Set {n_capped} BMI values outside [10,80] to NaN")
    print(
        f"    BMI availability: {cohort_fix.bmi.notna().sum()} ({100*cohort_fix.bmi.notna().mean():.1f}%)"
    )
    bmi_fix = cohort_fix.bmi[cohort_fix.bmi.notna()]
    print(f"    New range: {bmi_fix.min():.1f} – {bmi_fix.max():.1f}")
    print(f"    New mean (SD): {bmi_fix.mean():.1f} ({bmi_fix.std():.1f})")

    fix_path = os.path.join(RESULTS, "probe_cohort_bmifixed.csv")
    cohort_fix.to_csv(fix_path, index=False)
    print(f"    Saved: {fix_path}")

    # ── BMI audit CSV ────────────────────────────────────────────
    audit = pd.DataFrame(
        {
            "metric": [
                "n_total",
                "n_bmi_available",
                "pct_available",
                "mean_raw",
                "median_raw",
                "sd_raw",
                "max_raw",
                "n_outside_10_80",
                "n_above_80",
                "n_above_100",
                "mean_capped",
                "median_capped",
                "sd_capped",
            ],
            "value": [
                len(cohort),
                n_bmi,
                f"{pct_bmi:.1f}",
                f"{bmi.mean():.1f}",
                f"{bmi.median():.1f}",
                f"{bmi.std():.1f}",
                f"{bmi.max():.1f}",
                n_outside_10_80,
                (bmi > 80).sum(),
                (bmi > 100).sum(),
                f"{bmi_fix.mean():.1f}",
                f"{bmi_fix.median():.1f}",
                f"{bmi_fix.std():.1f}",
            ],
        }
    )
    audit_path = os.path.join(RESULTS, "probe_bmi_audit.csv")
    audit.to_csv(audit_path, index=False)
    print(f"    Saved: {audit_path}")

    # ── MIMIC BMI check (should be clean) ────────────────────────
    print(f"\n  ── MIMIC BMI check ──")
    mimic_path = os.path.join(RESULTS, "04_mimic_cohort.csv")
    if os.path.exists(mimic_path):
        m = pd.read_csv(mimic_path)
        bmi_m = m.bmi[m.bmi.notna()]
        print(f"    n={len(bmi_m)}, range: {bmi_m.min():.1f} – {bmi_m.max():.1f}")
        print(f"    Outside [10,80]: {((bmi_m < 10) | (bmi_m > 80)).sum()}")
        print(f"    MIMIC ETL already caps at [10,80] → clean")


# =====================================================================
# SECTION B: APACHE IV SCORE EXTRACTION
# =====================================================================
def section_b_apache():
    print(f"\n{'=' * 65}")
    print("B. APACHE IV SCORE EXTRACTION")
    print("=" * 65)

    cohort_path = os.path.join(RESULTS, "01_analysis_a_cohort.csv")
    cohort = pd.read_csv(cohort_path)
    pids = set(cohort.patientunitstayid)

    # ── Check if APACHE already in cohort ────────────────────────
    apache_cols = [c for c in cohort.columns if "apache" in c.lower()]
    if apache_cols:
        print(f"  APACHE-related columns already in cohort: {apache_cols}")
        for c in apache_cols:
            n_avail = cohort[c].notna().sum()
            print(f"    {c}: {n_avail} available ({100*n_avail/len(cohort):.1f}%)")

    # ── Load apachePatientResult ─────────────────────────────────
    if not EICU_ROOT:
        print("  ERROR: eICU data root not found")
        return

    apr_path = find_csv(EICU_ROOT, "apachePatientResult")
    if not apr_path:
        print("  ERROR: apachePatientResult not found")
        return

    apr = pd.read_csv(apr_path, low_memory=False)
    apr.columns = apr.columns.str.lower()
    print(f"  apachePatientResult: {len(apr)} rows")

    # Filter to cohort
    apr_cohort = apr[apr.patientunitstayid.isin(pids)].copy()
    print(
        f"  In cohort: {len(apr_cohort)} rows, {apr_cohort.patientunitstayid.nunique()} unique patients"
    )

    # Available columns
    score_cols = [
        c for c in apr_cohort.columns if "score" in c.lower() or "apache" in c.lower()
    ]
    print(f"  Score-related columns: {score_cols}")

    # ── APACHE version distribution ──────────────────────────────
    if "apacheversion" in apr_cohort.columns:
        print(f"\n  APACHE version distribution:")
        for v, n in apr_cohort.apacheversion.value_counts().items():
            print(f"    {v}: {n}")

    # ── Extract APACHE IV score ──────────────────────────────────
    # Try different column names
    score_col = None
    for candidate in ["apachescore", "acutephysiologyscore", "apache_score"]:
        if candidate in apr_cohort.columns:
            score_col = candidate
            break

    if score_col is None:
        print("  WARNING: no APACHE score column found. Available columns:")
        print(f"    {list(apr_cohort.columns)}")
        # Try to find any numeric column that looks like a score
        for c in apr_cohort.columns:
            if apr_cohort[c].dtype in [np.float64, np.int64, "Int64"]:
                vals = apr_cohort[c].dropna()
                if len(vals) > 0 and vals.median() > 10 and vals.median() < 200:
                    print(
                        f"    Candidate: {c} — median={vals.median():.0f}, range={vals.min():.0f}–{vals.max():.0f}"
                    )
        return

    print(f"\n  Using score column: {score_col}")
    scores = apr_cohort[apr_cohort[score_col].notna()].copy()

    # If multiple APACHE versions per patient, prefer IV
    if "apacheversion" in scores.columns:
        scores_iv = scores[
            scores.apacheversion.astype(str).str.contains("IV|4", na=False)
        ]
        if len(scores_iv) > 0:
            print(f"  APACHE IV rows: {len(scores_iv)}")
            scores = scores_iv

    # One score per patient (take first / highest version)
    scores_uniq = (
        scores.sort_values(score_col, ascending=False)
        .groupby("patientunitstayid")
        .first()
        .reset_index()
    )

    n_available = len(scores_uniq)
    pct = 100 * n_available / len(cohort)
    print(f"\n  APACHE score available for: {n_available}/{len(cohort)} ({pct:.1f}%)")

    if n_available > 0:
        sc = scores_uniq[score_col]
        print(f"  Distribution:")
        print(f"    Mean (SD): {sc.mean():.1f} ({sc.std():.1f})")
        print(
            f"    Median [IQR]: {sc.median():.0f} [{sc.quantile(0.25):.0f}–{sc.quantile(0.75):.0f}]"
        )
        print(f"    Range: {sc.min():.0f} – {sc.max():.0f}")

        # Check by treatment group
        merged = cohort[["patientunitstayid", "mg_supplementation"]].merge(
            scores_uniq[["patientunitstayid", score_col]],
            on="patientunitstayid",
            how="left",
        )
        for trt in [0, 1]:
            vals = merged[score_col][merged.mg_supplementation == trt].dropna()
            label = "Supplemented" if trt == 1 else "Unsupplemented"
            if len(vals) > 0:
                print(
                    f"    {label}: mean={vals.mean():.1f}, median={vals.median():.0f} (n={len(vals)})"
                )

        # Also extract predicted mortality for potential use
        mort_cols = [
            c
            for c in scores_uniq.columns
            if "mortality" in c.lower() or "predicted" in c.lower()
        ]

    # ── Save for R probe ─────────────────────────────────────────
    out_cols = ["patientunitstayid", score_col]
    # Add predicted mortality if available
    for mc in ["predictedhospitalmortality", "predictedicumortality"]:
        if mc in scores_uniq.columns:
            out_cols.append(mc)
    # Add APACHE version
    if "apacheversion" in scores_uniq.columns:
        out_cols.append("apacheversion")

    out_path = os.path.join(RESULTS, "probe_apache_scores.csv")
    scores_uniq[out_cols].to_csv(out_path, index=False)
    print(f"\n  Saved: {out_path} ({len(scores_uniq)} rows)")

    # ── Feasibility assessment ───────────────────────────────────
    print(f"\n  ── FEASIBILITY FOR PS MODEL ──")
    if pct >= 80:
        print(f"    ✓ {pct:.0f}% availability — safe to add to PS model")
        print(f"    Strategy: add as covariate with median imputation for missing")
    elif pct >= 50:
        print(f"    ⚠ {pct:.0f}% availability — usable with imputation caveat")
        print(f"    Strategy: add with missing indicator, report sensitivity")
    else:
        print(f"    ✗ {pct:.0f}% availability — too sparse for PS model")
        print(f"    Strategy: skip, note in limitations")


# =====================================================================
# SECTION C: ACHIEVED MG DISTRIBUTION
# (quantifies the baseline ≠ achieved gap)
# =====================================================================
def section_c_achieved_mg():
    print(f"\n{'=' * 65}")
    print("C. ACHIEVED MG DISTRIBUTION (baseline ≠ achieved)")
    print("=" * 65)

    all_rows = []

    for db_name, cohort_file in [
        ("eICU", "01_analysis_a_cohort.csv"),
        ("MIMIC", "04_mimic_cohort.csv"),
    ]:
        path = os.path.join(RESULTS, cohort_file)
        if not os.path.exists(path):
            print(f"  {db_name}: cohort not found, skipping")
            continue

        d = pd.read_csv(path)
        trt_col = "mg_supplementation"
        fu_col = "followup_mg_value"
        bl_col = "first_mg_value"

        if fu_col not in d.columns:
            print(f"  {db_name}: {fu_col} not found, skipping")
            continue

        print(f"\n  ── {db_name} (N={len(d)}) ──")

        for trt_val, trt_label in [(1, "Supplemented"), (0, "Unsupplemented")]:
            sub = d[(d[trt_col] == trt_val) & d[fu_col].notna()]
            n = len(sub)
            if n == 0:
                continue

            fu = sub[fu_col]
            bl = sub[bl_col]
            delta = fu - bl

            print(f"\n    {trt_label} (n={n} with follow-up Mg):")
            print(f"      Baseline Mg:  mean={bl.mean():.2f}, median={bl.median():.2f}")
            print(f"      Follow-up Mg: mean={fu.mean():.2f}, median={fu.median():.2f}")
            print(
                f"      ΔMg:          mean={delta.mean():.2f}, median={delta.median():.2f}"
            )
            print(f"      Follow-up Mg quantiles:")
            for q in [0.10, 0.25, 0.50, 0.75, 0.90]:
                print(f"        P{int(q*100):02d}: {fu.quantile(q):.2f}")

            # Key thresholds
            print(f"      Threshold crossing (follow-up Mg):")
            for thresh in [2.0, 2.3, 2.6, 3.0]:
                n_above = (fu > thresh).sum()
                pct = 100 * n_above / n
                print(f"        >{thresh}: {n_above}/{n} ({pct:.1f}%)")

            all_rows.append(
                {
                    "db": db_name,
                    "group": trt_label,
                    "n": n,
                    "bl_mean": round(bl.mean(), 3),
                    "bl_median": round(bl.median(), 3),
                    "fu_mean": round(fu.mean(), 3),
                    "fu_median": round(fu.median(), 3),
                    "delta_mean": round(delta.mean(), 3),
                    "delta_median": round(delta.median(), 3),
                    "pct_fu_gt_2.0": round(100 * (fu > 2.0).mean(), 1),
                    "pct_fu_gt_2.3": round(100 * (fu > 2.3).mean(), 1),
                    "pct_fu_gt_2.6": round(100 * (fu > 2.6).mean(), 1),
                    "pct_fu_gt_3.0": round(100 * (fu > 3.0).mean(), 1),
                }
            )

        # ── By baseline Mg stratum ───────────────────────────────
        print(f"\n    By baseline Mg stratum (supplemented only):")
        trt = d[(d[trt_col] == 1) & d[fu_col].notna()]
        if len(trt) > 0:
            cuts = [0, 1.8, 2.0, 2.3, float("inf")]
            labels = ["<1.8", "1.8-2.0", "2.0-2.3", ">2.3"]
            trt = trt.copy()
            trt["mg_stratum"] = pd.cut(
                trt[bl_col], bins=cuts, labels=labels, right=False, include_lowest=True
            )
            for s in labels:
                sub_s = trt[trt.mg_stratum == s]
                if len(sub_s) == 0:
                    continue
                fu_s = sub_s[fu_col]
                delta_s = fu_s - sub_s[bl_col]
                pct_gt23 = 100 * (fu_s > 2.3).mean()
                print(
                    f"      Baseline {s}: n={len(sub_s)}, "
                    f"follow-up mean={fu_s.mean():.2f}, "
                    f"Δ={delta_s.mean():+.2f}, "
                    f">{2.3}: {pct_gt23:.0f}%"
                )

    # ── Save ─────────────────────────────────────────────────────
    if all_rows:
        out = pd.DataFrame(all_rows)
        out_path = os.path.join(RESULTS, "probe_achieved_mg.csv")
        out.to_csv(out_path, index=False)
        print(f"\n  Saved: {out_path}")

    # ── Key verdict ──────────────────────────────────────────────
    print(f"\n  ── VERDICT FOR DISCUSSION ──")
    print(f"  If supplemented patients rarely cross >2.3 on follow-up:")
    print(f"    → baseline≠achieved argument is strong")
    print(f"    → threshold observation is clearly hypothesis-generating")
    print(f"  If supplemented patients commonly cross >2.3:")
    print(f"    → baseline≠achieved argument needs softening")
    print(f"    → threshold may be more directly interpretable")


# =====================================================================
# MAIN
# =====================================================================
if __name__ == "__main__":
    print("=" * 65)
    print("probe_bmi_apache.py — Pre-submission QC probes")
    print("=" * 65)

    section_a_bmi_audit()
    section_b_apache()
    section_c_achieved_mg()

    print(f"\n{'=' * 65}")
    print("DONE. Next: Rscript probe_v5_experiments.R")
    print("=" * 65)
