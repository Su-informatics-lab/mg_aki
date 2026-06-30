#!/usr/bin/env python3
"""
gen_consort.py -- CONSORT flow numbers for draw.io annotation

Computes all cohort-construction counts from the patient-level CSVs
and prints them in a structured format.  Does NOT generate a figure;
use the output to populate a draw.io template.

Upstream (ETL-stage) counts that cannot be derived from results CSVs
are hardcoded from the 01_etl.py run log and verified where possible.

Reads:
  results/did_all_{db}.csv
  results/did_pairs_primary_yet_untreated_{db}.csv

Outputs:
  results/consort_numbers.csv   (machine-readable)
  stdout                        (human-readable, for draw.io)

Usage:
  python gen_consort.py
"""

import os

import pandas as pd

RESULTS = os.path.expanduser("~/mg_aki/results")

# ================================================================
# ETL-STAGE COUNTS (from 01_etl.py run log, not in results CSVs)
# These are verified against the manuscript Table 1 / abstract.
# ================================================================
ETL = {
    "mimic": {
        "total_icu": 94458,
        "cardiac_surgery": 13404,
        "eskd_excluded": 429,
    },
    "eicu": {
        "total_icu": 200859,
        "cardiac_surgery": 26725,
        "eskd_excluded": 903,
    },
}


def compute_consort(tag):
    """Compute CONSORT numbers for one database from results CSVs."""
    etl = ETL[tag]

    all_path = os.path.join(RESULTS, f"did_all_{tag}.csv")
    pairs_path = os.path.join(RESULTS, f"did_pairs_primary_yet_untreated_{tag}.csv")

    if not os.path.exists(all_path):
        print(f"  {tag}: did_all not found, skipping")
        return None

    df = pd.read_csv(all_path)
    pairs = pd.read_csv(pairs_path) if os.path.exists(pairs_path) else None

    # -- Derived from did_all --
    n_total = len(df)
    trt = df[df.treated == 1]
    ctl = df[df.treated == 0]
    n_treated = len(trt)
    n_control = len(ctl)
    n_eligible = n_total

    # ETL-stage derived
    post_eskd = etl["cardiac_surgery"] - etl["eskd_excluded"]
    n_excl_noncardiac = etl["total_icu"] - etl["cardiac_surgery"]

    # -- From pairs --
    n_matched_pairs = len(pairs) if pairs is not None else 0
    n_matched_trt = pairs.trt_pid.nunique() if pairs is not None else 0
    match_rate = 100 * n_matched_trt / n_treated if n_treated > 0 else 0
    n_unique_ctl = pairs.ctl_pid.nunique() if pairs is not None else 0
    reuse_ratio = n_matched_pairs / n_unique_ctl if n_unique_ctl > 0 else 0

    nums = {
        "db": tag,
        # ETL stage
        "total_icu": etl["total_icu"],
        "excl_noncardiac": n_excl_noncardiac,
        "cardiac_surgery": etl["cardiac_surgery"],
        "eskd_excluded": etl["eskd_excluded"],
        "post_eskd": post_eskd,
        # Treatment split
        "treated_eligible": n_treated,
        "control_eligible": n_control,
        "eligible_total": n_eligible,
        # Matching
        "matched_pairs": n_matched_pairs,
        "matched_treated": n_matched_trt,
        "match_rate_pct": round(match_rate, 1),
        "unique_controls_used": n_unique_ctl,
        "control_reuse_ratio": round(reuse_ratio, 2),
    }
    return nums


def print_consort(nums):
    """Print CONSORT numbers in a structured human-readable format."""
    tag = nums["db"]
    db_label = "MIMIC-IV" if tag == "mimic" else "eICU-CRD"
    sep = "=" * 55

    print(f"\n{sep}")
    print(f"  CONSORT NUMBERS: {db_label}")
    print(sep)
    print(f"  Total ICU admissions:           {nums['total_icu']:>10,}")
    print(f"  (-) Non-cardiac / <18 / repeat: {nums['excl_noncardiac']:>10,}")
    print(f"  = Cardiac surgery, first stay:  {nums['cardiac_surgery']:>10,}")
    print(f"  (-) ESKD (chronic dialysis):    {nums['eskd_excluded']:>10,}")
    print(f"  = Post-ESKD:                    {nums['post_eskd']:>10,}")
    print(f"  {'~' * 50}")
    print(f"  Treated (IV Mg, eligible):      {nums['treated_eligible']:>10,}")
    print(f"  Control (no IV Mg, eligible):   {nums['control_eligible']:>10,}")
    print(f"  Eligible total:                 {nums['eligible_total']:>10,}")
    print(f"  {'~' * 50}")
    print(f"  Matched pairs:                  {nums['matched_pairs']:>10,}")
    print(f"  Matched treated (unique):       {nums['matched_treated']:>10,}")
    print(f"  Match rate:                     {nums['match_rate_pct']:>9.1f}%")
    print(f"  Unique controls in matching:    {nums['unique_controls_used']:>10,}")
    print(f"  Control reuse ratio:            {nums['control_reuse_ratio']:>10.2f}")
    print(sep)


def main():
    print("gen_consort.py -- CONSORT numbers for draw.io\n")

    all_nums = []
    for tag in ["mimic", "eicu"]:
        nums = compute_consort(tag)
        if nums is not None:
            print_consort(nums)
            all_nums.append(nums)

    # Save machine-readable CSV
    if all_nums:
        out = pd.DataFrame(all_nums)
        out_path = os.path.join(RESULTS, "consort_numbers.csv")
        out.to_csv(out_path, index=False)
        print(f"\nSaved: {out_path}")

    # Combined summary for quick copy-paste
    if len(all_nums) == 2:
        m, e = all_nums[0], all_nums[1]
        print(f"\n{'=' * 55}")
        print("  COMBINED SUMMARY (for manuscript text)")
        print(f"{'=' * 55}")
        print(f"  Total:  {m['eligible_total'] + e['eligible_total']:,} eligible")
        print(f"          {m['matched_pairs'] + e['matched_pairs']:,} matched pairs")
        print(
            f"  MIMIC:  {m['treated_eligible']:,} treated, {m['control_eligible']:,} control"
        )
        print(
            f"          {m['matched_pairs']:,} pairs ({m['match_rate_pct']}% matched)"
        )
        print(
            f"  eICU:   {e['treated_eligible']:,} treated, {e['control_eligible']:,} control"
        )
        print(
            f"          {e['matched_pairs']:,} pairs ({e['match_rate_pct']}% matched)"
        )


if __name__ == "__main__":
    main()
