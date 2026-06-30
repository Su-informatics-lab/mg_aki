#!/usr/bin/env python3
"""
gen_consort.py -- CONSORT flow numbers for draw.io annotation

Reads did_consort.csv (authoritative ETL-stage counts from 01_etl.py)
+ did_all / did_pairs CSVs for matching-stage verification.

Prints ALL exclusion steps with exact numbers.
Does NOT generate a figure; populate draw.io manually.

Usage: python gen_consort.py
"""

import os

import pandas as pd

RESULTS = os.path.expanduser("~/mg_aki/results")


def load(name):
    p = os.path.join(RESULTS, name)
    if not os.path.exists(p):
        print(f"  WARN: {name} not found")
        return None
    return pd.read_csv(p)


def run(tag):
    db_label = "MIMIC-IV" if tag == "MIMIC" else "eICU-CRD"
    sep = "=" * 60

    # -- ETL-stage numbers from did_consort.csv --
    consort_df = load("did_consort.csv")
    if consort_df is None:
        print("  did_consort.csv missing -- run 01_etl.py first")
        return
    row = consort_df[consort_df.db == tag]
    if len(row) == 0:
        print(f"  No {tag} row in did_consort.csv")
        return
    c = row.iloc[0].to_dict()

    # -- Verification data --
    tag_lc = tag.lower()
    did_all = load(f"did_all_{tag_lc}.csv")
    pairs = load(f"did_pairs_primary_yet_untreated_{tag_lc}.csv")

    # Derived
    total_icu = int(c["total_icu"])
    cardiac = int(c["cardiac_adult_first"])
    excl_noncardiac = total_icu - cardiac
    eskd = cardiac - int(c["post_eskd"])
    post_eskd = int(c["post_eskd"])
    ivmg = int(c["treated_any_ivmg"])
    no_ivmg = int(c["control_no_ivmg"])
    cr_pre_icu = int(c["treated_has_cr_pre"])
    no_cr_pre = int(c["treated_no_cr_pre"])
    fallback = int(c["treated_hosp_fallback"])
    cr_pre_total = int(c["treated_has_cr_pre_with_fallback"])
    excl_cr_high_trt = int(c["excl_cr_high"])
    treated_final = int(c["treated_final"])
    ctrl_has_2cr = int(c["control_has_2cr"])
    control_final = int(c["control_final"])
    eligible_total = treated_final + control_final

    # Control Cr>=4.0 exclusion (not tracked separately in consort dict)
    excl_cr_high_ctrl = ctrl_has_2cr - control_final
    no_2cr = no_ivmg - ctrl_has_2cr

    # Matching
    n_pairs = len(pairs) if pairs is not None else 0
    n_unique_ctl = pairs.ctl_pid.nunique() if pairs is not None else 0
    n_unique_trt = pairs.trt_pid.nunique() if pairs is not None else 0
    match_rate = 100 * n_unique_trt / treated_final if treated_final > 0 else 0
    reuse = n_pairs / n_unique_ctl if n_unique_ctl > 0 else 0

    # -- Verify consistency --
    checks = []
    if did_all is not None:
        da_trt = int((did_all.treated == 1).sum())
        da_ctl = int((did_all.treated == 0).sum())
        checks.append(("treated_final vs did_all treated", treated_final, da_trt))
        checks.append(("control_final vs did_all control", control_final, da_ctl))
        checks.append(("eligible_total vs did_all total", eligible_total, len(did_all)))

    # -- Print --
    print(f"\n{sep}")
    print(f"  CONSORT: {db_label}")
    print(sep)
    print(f"\n  [ETL STAGE -- 01_etl.py]\n")
    print(f"  1. Total ICU admissions:                 {total_icu:>10,}")
    print(f"     (-) Non-cardiac / age<18 / repeat:    {excl_noncardiac:>10,}")
    print(f"  2. Cardiac surgery, adult, 1st stay:     {cardiac:>10,}")
    print(f"     (-) ESKD (pre-surgery dialysis):      {eskd:>10,}")
    print(f"  3. Post-ESKD:                            {post_eskd:>10,}")
    print()
    print(f"  4. TREATMENT SPLIT (during entire ICU stay):")
    print(f"     Received IV Mg at any time:           {ivmg:>10,}")
    print(f"     Never received IV Mg:                 {no_ivmg:>10,}")
    print()
    print(f"  5. CR FILTERS -- treated arm:")
    print(f"     Has ICU pre-T0 Cr:                    {cr_pre_icu:>10,}")
    print(f"     No ICU pre-T0 Cr:                     {no_cr_pre:>10,}")
    print(f"     + Hospital Cr fallback:               {fallback:>10,}")
    print(f"     Total with any pre-T0 Cr:             {cr_pre_total:>10,}")
    print(f"     (-) Baseline Cr >= 4.0:               {excl_cr_high_trt:>10,}")
    print(f"     = Treated eligible:                   {treated_final:>10,}")
    print()
    print(f"  6. CR FILTERS -- control arm (never-treated):")
    print(f"     Has >=2 postop Cr:                    {ctrl_has_2cr:>10,}")
    print(f"     <2 postop Cr (excluded):              {no_2cr:>10,}")
    print(f"     (-) Baseline Cr >= 4.0:               {excl_cr_high_ctrl:>10,}")
    print(f"     = Control eligible (never-treated):   {control_final:>10,}")
    print()
    print(f"  7. ELIGIBLE COHORT:                      {eligible_total:>10,}")
    print(
        f"     (treated {treated_final:,} + never-treated controls {control_final:,})"
    )
    print()
    print(f"  [MATCHING STAGE -- 02_psm.R]\n")
    print(f"  8. Risk-set criteria at each T0 (time-varying):")
    print(f"     - Still in ICU at T0")
    print(f"     - No prevalent AKI at T0")
    print(f"     - Has Cr measurement by T0")
    print(f"     - Not yet received IV Mg by T0+48h")
    print(f"     - 1:1 PSM, replacement, caliper 0.2 SD, 19 covariates")
    print(f"     NOTE: Not-yet-treated patients from the 'received IV Mg'")
    print(f"           group ARE eligible as controls at earlier T0 values.")
    print()
    print(f"  9. MATCHED PAIRS:                        {n_pairs:>10,}")
    print(f"     Matched treated (unique):             {n_unique_trt:>10,}")
    print(f"     Match rate:                           {match_rate:>9.1f}%")
    print(f"     Unique controls used:                 {n_unique_ctl:>10,}")
    print(f"     Control reuse ratio:                  {reuse:>10.2f}")

    # Consistency
    print(f"\n  [CONSISTENCY CHECKS]")
    all_ok = True
    for label, a, b in checks:
        ok = "OK" if a == b else f"MISMATCH ({a} vs {b})"
        if a != b:
            all_ok = False
        print(f"     {label}: {ok}")
    if all_ok:
        print(f"     All checks passed.")
    print(sep)
    return c


def main():
    print("gen_consort.py -- CONSORT numbers for draw.io\n")
    for tag in ["MIMIC", "eICU"]:
        run(tag)


if __name__ == "__main__":
    main()
