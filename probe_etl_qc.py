#!/usr/bin/env python3
"""
probe_etl_qc.py — QC the ETL output after bash run.sh 1

Checks the three new column families:
  1. rrt_offset_h / has_rrt  — RRT detection
  2. death_offset_h          — Death timing
  3. Table 1 labs in did_labs_all (hemoglobin, wbc, platelets, albumin)

Also sanity-checks existing columns haven't broken.

Run: python probe_etl_qc.py
"""

import os

import pandas as pd

RESULTS = os.path.expanduser("~/mg_aki/results")
SEP = "=" * 70


def pct(n, total):
    return f"{n:,} ({100*n/max(total,1):.1f}%)" if total > 0 else "0"


for db in ["mimic", "eicu"]:
    path_all = os.path.join(RESULTS, f"did_all_{db}.csv")
    path_cr = os.path.join(RESULTS, f"did_cr_all_{db}.csv")
    path_lab = os.path.join(RESULTS, f"did_labs_all_{db}.csv")

    if not os.path.exists(path_all):
        print(f"\n  {db}: did_all not found, skipping")
        continue

    print(f"\n{SEP}")
    print(f"  QC: {db.upper()}")
    print(SEP)

    df = pd.read_csv(path_all)
    n = len(df)
    n_trt = (df.treated == 1).sum()
    n_ctl = (df.treated == 0).sum()
    print(f"\n  did_all: {n:,} patients ({n_trt:,} treated, {n_ctl:,} control)")

    # ── 1. Basic columns (sanity check) ──
    print(f"\n  1. Basic columns:")
    for col in [
        "pid",
        "treated",
        "age",
        "is_female",
        "first_cr",
        "egfr",
        "mg_offset_h",
        "icu_discharge_h",
        "hosp_mortality",
    ]:
        if col in df.columns:
            nn = df[col].notna().sum()
            print(f"    {col:<25s} non-null: {pct(nn, n)}")
        else:
            print(f"    {col:<25s} *** MISSING ***")

    # ── 2. RRT columns ──
    print(f"\n  2. RRT detection:")
    if "has_rrt" in df.columns:
        n_rrt = df.has_rrt.sum()
        rrt_trt = df[df.treated == 1].has_rrt.sum()
        rrt_ctl = df[df.treated == 0].has_rrt.sum()
        print(f"    has_rrt total:   {pct(n_rrt, n)}")
        print(f"    has_rrt treated: {pct(rrt_trt, n_trt)}")
        print(f"    has_rrt control: {pct(rrt_ctl, n_ctl)}")
        if "rrt_offset_h" in df.columns:
            rrt_vals = df.rrt_offset_h.dropna()
            if len(rrt_vals) > 0:
                print(
                    f"    rrt_offset_h:    median={rrt_vals.median():.1f}h, "
                    f"IQR=[{rrt_vals.quantile(.25):.1f}-{rrt_vals.quantile(.75):.1f}]"
                )
            # Sanity: RRT should be rare in cardiac surgery (expect 1-5%)
            rate = 100 * n_rrt / n
            if rate > 15:
                print(
                    f"    ⚠ WARNING: RRT rate {rate:.1f}% seems high for cardiac surgery"
                )
            elif rate < 0.1:
                print(f"    ⚠ WARNING: RRT rate {rate:.1f}% seems very low")
            else:
                print(f"    ✓ RRT rate {rate:.1f}% plausible for cardiac surgery")
    else:
        print(f"    *** has_rrt MISSING — ETL did not produce it ***")

    # ── 3. Death offset ──
    print(f"\n  3. Death timing:")
    if "death_offset_h" in df.columns:
        n_death = df.death_offset_h.notna().sum()
        n_mort = df.hosp_mortality.sum() if "hosp_mortality" in df.columns else "?"
        print(f"    death_offset_h non-null: {pct(n_death, n)}")
        print(f"    hosp_mortality == 1:     {n_mort}")
        # These should be close but not necessarily equal
        # (hosp_mortality can be 1 even if death time is missing in eICU)
        if isinstance(n_mort, int) and n_death > 0:
            d_vals = df.death_offset_h.dropna()
            print(
                f"    death_offset_h:  median={d_vals.median():.1f}h, "
                f"IQR=[{d_vals.quantile(.25):.1f}-{d_vals.quantile(.75):.1f}]"
            )
            # Negative offsets would be a bug (died before ICU admission?)
            n_neg = (d_vals < 0).sum()
            if n_neg > 0:
                print(f"    ⚠ WARNING: {n_neg} negative death_offset_h values!")
            else:
                print(f"    ✓ No negative offsets")
            # Check consistency: death_offset_h present → hosp_mortality should be 1
            has_death = df[df.death_offset_h.notna()]
            inconsistent = (has_death.hosp_mortality != 1).sum()
            if inconsistent > 0:
                print(
                    f"    ⚠ WARNING: {inconsistent} patients with death_offset but hosp_mortality != 1"
                )
            else:
                print(f"    ✓ Consistent with hosp_mortality")
    else:
        print(f"    *** death_offset_h MISSING — ETL did not produce it ***")

    # ── 4. Existing secondary outcomes (not broken?) ──
    print(f"\n  4. Secondary outcomes (existing):")
    for col in [
        "poaf",
        "encephalopathy_delirium",
        "transfusion",
        "reintubation",
        "poaf_icd",
        "encephalopathy_icd",
        "vent_arrhythmia",
    ]:
        if col in df.columns:
            nn = df[col].sum()
            print(f"    {col:<30s} events: {pct(nn, n)}")
        else:
            print(f"    {col:<30s} not present")

    # ── 5. Table 1 labs in did_labs_all ──
    print(f"\n  5. Table 1 labs in did_labs_all:")
    if os.path.exists(path_lab):
        labs = pd.read_csv(path_lab)
        pid_col = (
            "patientunitstayid" if "patientunitstayid" in labs.columns else "stay_id"
        )
        total_pts = labs[pid_col].nunique()
        print(f"    Total measurements: {len(labs):,} across {total_pts:,} patients")
        for ln in sorted(labs.lab_name.unique()):
            sub = labs[labs.lab_name == ln]
            n_pts = sub[pid_col].nunique()
            vals = sub["value"].dropna()
            if len(vals) > 0:
                print(
                    f"    {ln:<20s} {len(sub):>10,} rows  {n_pts:>7,} pts  "
                    f"median={vals.median():>8.1f}  IQR=[{vals.quantile(.25):.1f}-{vals.quantile(.75):.1f}]"
                )
            else:
                print(
                    f"    {ln:<20s} {len(sub):>10,} rows  {n_pts:>7,} pts  (all null)"
                )

        # Check new labs specifically
        new_labs = ["hemoglobin", "wbc", "platelets", "albumin"]
        missing = [l for l in new_labs if l not in labs.lab_name.values]
        if missing:
            print(f"\n    ⚠ MISSING new labs: {missing}")
        else:
            print(f"\n    ✓ All 4 Table 1 labs present")
    else:
        print(f"    *** did_labs_all not found ***")

    # ── 6. Creatinine (not broken?) ──
    print(f"\n  6. Creatinine data:")
    if os.path.exists(path_cr):
        cr = pd.read_csv(path_cr)
        cr_pid = "patientunitstayid" if "patientunitstayid" in cr.columns else "stay_id"
        print(
            f"    did_cr_all: {len(cr):,} measurements, {cr[cr_pid].nunique():,} patients"
        )
        print(
            f"    Cr median: {cr.labresult.median():.2f}, "
            f"IQR=[{cr.labresult.quantile(.25):.2f}-{cr.labresult.quantile(.75):.2f}]"
        )
    else:
        print(f"    *** did_cr_all not found ***")

print(f"\n{SEP}")
print("  ETL QC COMPLETE")
print(f"  If all checks pass → run: bash run.sh 2 3 4 5 6 7 8")
print(SEP)
