#!/usr/bin/env python3
"""
gen_etables.py — Generate eTables 2 and 3 for supplement

  eTable 2: Covariate balance (raw vs matched SMD, all PS vars)
  eTable 3: Complete HTE matrix (all subgroups × all outcomes)

Reads: did_all_{db}.csv, did_hte_data_{db}.csv, did_hte_{db}.csv
Output: etable2_balance.csv, etable3_hte_matrix.csv

Usage: python gen_etables.py
"""

import os

import numpy as np
import pandas as pd

RESULTS = os.path.expanduser("~/mg_aki/results")
DBS = ["mimic", "eicu"]

PS_VARS = [
    "age",
    "is_female",
    "bmi",
    "surg_cabg",
    "surg_valve",
    "surg_combined",
    "heart_failure",
    "hypertension",
    "diabetes",
    "ckd",
    "copd",
    "pvd",
    "stroke",
    "liver_disease",
    "egfr",
]

NICE = {
    "age": "Age",
    "is_female": "Female sex",
    "bmi": "BMI",
    "surg_cabg": "CABG",
    "surg_valve": "Valve",
    "surg_combined": "Combined",
    "heart_failure": "Heart failure",
    "hypertension": "Hypertension",
    "diabetes": "Diabetes",
    "ckd": "CKD",
    "copd": "COPD",
    "pvd": "PVD",
    "stroke": "Stroke/TIA",
    "liver_disease": "Liver disease",
    "egfr": "eGFR",
    "last_calcium": "Calcium",
    "last_lactate": "Lactate",
    "last_lactate_missing": "Lactate missing",
    "last_heartrate": "Heart rate",
    "last_magnesium": "Magnesium",
    "last_potassium": "Potassium",
}


def abs_smd(df, var):
    t1 = df.loc[df.treated == 1, var].dropna()
    t0 = df.loc[df.treated == 0, var].dropna()
    if len(t1) < 5 or len(t0) < 5:
        return np.nan
    sp = np.sqrt((t1.var() + t0.var()) / 2)
    return abs(t1.mean() - t0.mean()) / sp if sp > 1e-10 else 0


# ═══════════════════════════════════════════════════════════════════
# eTable 2: Covariate balance
# ═══════════════════════════════════════════════════════════════════
def etable2():
    print("\n── eTable 2: Covariate Balance ──")
    rows = []
    for tag in DBS:
        db = "MIMIC-IV" if tag == "mimic" else "eICU-CRD"
        raw = pd.read_csv(os.path.join(RESULTS, f"did_all_{tag}.csv"))
        matched = pd.read_csv(os.path.join(RESULTS, f"did_hte_data_{tag}.csv"))

        # Find available PS vars
        avail = [v for v in PS_VARS if v in raw.columns]
        # Add lab vars from matched if available
        for v in matched.columns:
            if v.startswith("last_") and v in NICE and v not in avail:
                avail.append(v)

        for v in avail:
            raw_smd = abs_smd(raw, v) if v in raw.columns else np.nan
            mat_smd = abs_smd(matched, v) if v in matched.columns else np.nan
            rows.append(
                {
                    "Database": db,
                    "Variable": NICE.get(v, v),
                    "Raw SMD": f"{raw_smd:.3f}" if not np.isnan(raw_smd) else "—",
                    "Matched SMD": f"{mat_smd:.3f}" if not np.isnan(mat_smd) else "—",
                }
            )

    df = pd.DataFrame(rows)
    out = os.path.join(RESULTS, "etable2_balance.csv")
    df.to_csv(out, index=False)
    print(f"  Saved: {out}")

    # Print
    for db in ["MIMIC-IV", "eICU-CRD"]:
        sub = df[df.Database == db]
        print(f"\n  {db}:")
        print(f"  {'Variable':<25s} {'Raw':>8s} {'Matched':>8s}")
        print("  " + "─" * 43)
        for _, r in sub.iterrows():
            print(f"  {r['Variable']:<25s} {r['Raw SMD']:>8s} {r['Matched SMD']:>8s}")


# ═══════════════════════════════════════════════════════════════════
# eTable 3: Complete HTE matrix
# ═══════════════════════════════════════════════════════════════════
def etable3():
    print("\n── eTable 3: Complete HTE Matrix ──")
    oc_order = [
        "dcr_48h",
        "aki_48h",
        "aki_7d",
        "hosp_mortality",
        "poaf",
        "vent_arrhythmia",
    ]
    oc_labels = {
        "dcr_48h": "ΔCr 48h",
        "aki_48h": "48h AKI",
        "aki_7d": "7d AKI",
        "hosp_mortality": "Mortality",
        "poaf": "POAF",
        "vent_arrhythmia": "Vent arrhythmia",
    }

    all_rows = []
    for tag in DBS:
        db = "MIMIC-IV" if tag == "mimic" else "eICU-CRD"
        hte_path = os.path.join(RESULTS, f"did_hte_{tag}.csv")
        if not os.path.exists(hte_path):
            print(f"  {hte_path} not found")
            continue
        hte = pd.read_csv(hte_path)

        for _, r in hte.iterrows():
            if r.outcome not in oc_order:
                continue
            if r.type == "continuous":
                est_str = f"{r.est:+.4f}" if not pd.isna(r.est) else "—"
            else:
                est_str = f"{r.est:.2f}" if not pd.isna(r.est) else "—"
            p_str = f"{r.p:.4f}" if not pd.isna(r.p) else "—"
            sig = "*" if not pd.isna(r.p) and r.p < 0.05 else ""

            all_rows.append(
                {
                    "Database": db,
                    "Subgroup": r.subgroup,
                    "Outcome": oc_labels.get(r.outcome, r.outcome),
                    "n_trt": int(r.n_trt) if not pd.isna(r.n_trt) else 0,
                    "n_ctl": int(r.n_ctl) if not pd.isna(r.n_ctl) else 0,
                    "Estimate": est_str,
                    "P": p_str + sig,
                    "Rate_trt": (
                        f"{100*r.rate_trt:.1f}%" if not pd.isna(r.rate_trt) else ""
                    ),
                    "Rate_ctl": (
                        f"{100*r.rate_ctl:.1f}%" if not pd.isna(r.rate_ctl) else ""
                    ),
                    "RD": f"{100*r.rd:+.1f}%" if not pd.isna(r.rd) else "",
                }
            )

    df = pd.DataFrame(all_rows)
    out = os.path.join(RESULTS, "etable3_hte_matrix.csv")
    df.to_csv(out, index=False)
    print(f"  Saved: {out} ({len(df)} rows)")


if __name__ == "__main__":
    print("=" * 60)
    print("gen_etables.py — eTables for supplement")
    print("=" * 60)
    etable2()
    etable3()
    print(f"\n{'='*60}\nDONE\n{'='*60}")
