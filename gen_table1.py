#!/usr/bin/env python3
"""
gen_table1.py — Generate Table 1: Baseline Characteristics of Matched Pairs

Reads:
  did_pairs_primary_yet_untreated_{db}.csv  — matched pair IDs
  did_all_{db}.csv                          — full covariate data
  did_hte_data_{db}.csv                     — matched-pair-level data (fallback)

Outputs:
  table1_{db}.csv          — machine-readable
  table1_combined.csv      — MIMIC + eICU side-by-side (for manuscript)

Usage:
  python gen_table1.py              # both databases
  python gen_table1.py mimic        # one database
"""

import os
import sys

import numpy as np
import pandas as pd

RESULTS = os.path.expanduser("~/mg_aki/results")

# ── Variables to include in Table 1 ──
# (varname, display_name, type)
# type: "continuous" → mean (SD)
#         "binary"     → n (%)
#         "categorical"→ handled specially

DEMOG_VARS = [
    ("age", "Age, mean (SD), y", "continuous"),
    ("is_female", "Female sex, n (%)", "binary"),
    ("bmi", "BMI, mean (SD), kg/m²", "continuous"),
]

SURGERY_VARS = [
    ("surg_cabg", "CABG", "binary"),
    ("surg_valve", "Valve surgery", "binary"),
    ("surg_combined", "Combined CABG + valve", "binary"),
]

COMORBIDITY_VARS = [
    ("heart_failure", "Heart failure", "binary"),
    ("hypertension", "Hypertension", "binary"),
    ("diabetes", "Diabetes mellitus", "binary"),
    ("ckd", "Chronic kidney disease", "binary"),
    ("copd", "COPD", "binary"),
    ("pvd", "Peripheral vascular disease", "binary"),
    ("stroke", "Stroke/TIA", "binary"),
    ("liver_disease", "Liver disease", "binary"),
]

RENAL_VARS = [
    ("egfr", "eGFR, mean (SD), mL/min/1.73 m²", "continuous"),
    ("cr_pre", "Baseline creatinine, mean (SD), mg/dL", "continuous"),
]

LAB_VARS = [
    ("last_calcium", "Serum calcium, mean (SD), mg/dL", "continuous"),
    ("last_lactate", "Lactate, mean (SD), mmol/L", "continuous"),
    ("last_lactate_missing", "Lactate missing, n (%)", "binary"),
    ("last_heartrate", "Heart rate, mean (SD), bpm", "continuous"),
]

# Supplementary: Mg and K+ (not in PS but reported)
SUPP_VARS = [
    ("last_magnesium", "Serum magnesium, mean (SD), mg/dL", "continuous"),
    ("last_potassium", "Serum potassium, mean (SD), mEq/L", "continuous"),
]

TABLE1_LAB_VARS = [
    ("hemoglobin", "Hemoglobin, mean (SD), g/dL", "continuous"),
    ("wbc", "White blood cell count, mean (SD), ×10⁹/L", "continuous"),
    ("platelets", "Platelet count, mean (SD), ×10⁹/L", "continuous"),
    ("albumin", "Albumin, mean (SD), g/dL", "continuous"),
]

OUTCOME_VARS = [
    ("aki1_48h", "48-hour AKI (KDIGO ≥1), n (%)", "binary"),
    ("aki1_7d", "7-day AKI (KDIGO ≥1), n (%)", "binary"),
    ("aki2_7d", "7-day AKI (KDIGO ≥2), n (%)", "binary"),
    ("aki3_7d", "7-day AKI (KDIGO ≥3), n (%)", "binary"),
    ("hosp_mortality", "Hospital mortality, n (%)", "binary"),
    ("death_7d", "7-day mortality, n (%)", "binary"),
    ("vent_arrhythmia", "Ventricular arrhythmia, n (%)", "binary"),
]

ALL_SECTIONS = [
    ("Demographics", DEMOG_VARS),
    ("Surgery type, n (%)", SURGERY_VARS),
    ("Comorbidities, n (%)", COMORBIDITY_VARS),
    ("Renal function", RENAL_VARS),
    ("Laboratory values (closest to T₀)", LAB_VARS),
    ("Supplementary labs (not in PS model)", SUPP_VARS),
    ("Descriptive labs (Table 1 only)", TABLE1_LAB_VARS),
    ("Outcomes", OUTCOME_VARS),
]


def compute_smd(x1, x0, vtype="continuous"):
    """Compute absolute standardized mean difference."""
    x1 = x1.dropna()
    x0 = x0.dropna()
    if len(x1) < 5 or len(x0) < 5:
        return np.nan
    m1, m0 = x1.mean(), x0.mean()
    if vtype == "binary":
        p1, p0 = m1, m0
        sp = np.sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
    else:
        sp = np.sqrt((x1.var() + x0.var()) / 2)
    if sp < 1e-10:
        return 0.0
    return abs(m1 - m0) / sp


def format_continuous(series):
    """Format as 'mean (SD)' with appropriate precision."""
    m = series.mean()
    s = series.std()
    n_obs = series.notna().sum()
    if pd.isna(m):
        return "—"
    # Choose precision based on magnitude
    if abs(m) >= 10:
        return f"{m:.1f} ({s:.1f})"
    else:
        return f"{m:.2f} ({s:.2f})"


def format_binary(series, n_total):
    """Format as 'n (%)' ."""
    n_events = int(series.sum())
    pct = 100 * n_events / n_total if n_total > 0 else 0
    return f"{n_events:,} ({pct:.1f})"


def build_table1(db_tag):
    """Build Table 1 for one database."""
    tag = db_tag.lower()
    db_label = "MIMIC-IV" if tag == "mimic" else "eICU-CRD"
    print(f"\n{'='*60}")
    print(f"Table 1: {db_label}")
    print(f"{'='*60}")

    # Load matched pair IDs
    pairs_path = os.path.join(RESULTS, f"did_pairs_primary_yet_untreated_{tag}.csv")
    all_path = os.path.join(RESULTS, f"did_all_{tag}.csv")
    hte_path = os.path.join(RESULTS, f"did_hte_data_{tag}.csv")
    trt_path = os.path.join(RESULTS, f"did_treated_{tag}.csv")
    ctl_path = os.path.join(RESULTS, f"did_control_{tag}.csv")

    # Strategy: try hte_data first (has matched pairs with covariates),
    # then pairs + all_pts, then pairs + treated/control files
    trt = None
    ctl = None

    if not os.path.exists(pairs_path):
        print(f"  WARN: {pairs_path} not found, trying hte_data")
        if not os.path.exists(hte_path):
            print(f"  ERROR: no data files found for {tag}")
            return None
        # Fallback: use hte_data directly
        df = pd.read_csv(hte_path)
        df["pid"] = df["pid"].astype(str)
        trt = df[df.treated == 1]
        ctl = df[df.treated == 0]
        print(f"  Using hte_data: {len(trt)} treated, {len(ctl)} control")
    else:
        pairs = pd.read_csv(pairs_path)
        # Detect pair ID columns
        trt_pid_col = [
            c for c in pairs.columns if "trt" in c.lower() and "pid" in c.lower()
        ]
        ctl_pid_col = [
            c for c in pairs.columns if "ctl" in c.lower() and "pid" in c.lower()
        ]
        if not trt_pid_col or not ctl_pid_col:
            print(f"  Pairs columns: {list(pairs.columns)}")
            trt_pid_col = [pairs.columns[0]]  # assume first two are trt, ctl
            ctl_pid_col = [pairs.columns[1]]
        trt_pid_col = trt_pid_col[0]
        ctl_pid_col = ctl_pid_col[0]
        print(
            f"  Loaded {len(pairs)} matched pairs (cols: {trt_pid_col}, {ctl_pid_col})"
        )

        # Load full patient data — try did_all first, then combine treated+control
        all_pts = None
        if os.path.exists(all_path):
            all_pts = pd.read_csv(all_path)
        elif os.path.exists(trt_path) and os.path.exists(ctl_path):
            print(f"  did_all not found; combining did_treated + did_control")
            _t = pd.read_csv(trt_path)
            _c = pd.read_csv(ctl_path)
            shared_cols = [c for c in _t.columns if c in _c.columns]
            all_pts = pd.concat([_t[shared_cols], _c[shared_cols]], ignore_index=True)
        else:
            # Last resort: use hte_data
            print(f"  No patient data files found, falling back to hte_data")
            if os.path.exists(hte_path):
                df = pd.read_csv(hte_path)
                df["pid"] = df["pid"].astype(str)
                trt = df[df.treated == 1]
                ctl = df[df.treated == 0]
                print(f"  Using hte_data: {len(trt)} treated, {len(ctl)} control")
                # Skip the all_pts merge path
                all_pts = None

        if all_pts is not None:
            # Detect ID column: canonical 01_etl.py uses 'pid'
            id_candidates = ["pid", "stay_id", "patientunitstayid"]
            id_col = None
            for c in id_candidates:
                if c in all_pts.columns:
                    id_col = c
                    break
            if id_col is None:
                print(
                    f"  ERROR: no ID column found. Columns: {list(all_pts.columns[:15])}"
                )
                return None
            if id_col != "pid":
                all_pts["pid"] = all_pts[id_col].astype(str)
            else:
                all_pts["pid"] = all_pts["pid"].astype(str)
            print(f"  ID column: {id_col} ({len(all_pts)} patients)")

            # Also load hte_data for outcomes (aki_48h, aki_7d computed there)
            hte_data = None
            if os.path.exists(hte_path):
                hte_data = pd.read_csv(hte_path)
                hte_data["pid"] = hte_data["pid"].astype(str)

            # Get matched patient IDs
            trt_pids = set(pairs[trt_pid_col].astype(str))
            ctl_pids = set(pairs[ctl_pid_col].astype(str))

            trt = all_pts[all_pts.pid.isin(trt_pids)].copy()
            ctl = all_pts[all_pts.pid.isin(ctl_pids)].copy()

            # Merge outcome columns from hte_data if available
            if hte_data is not None:
                for oc in ["aki_48h", "aki_7d", "dcr_48h", "cr_pre"]:
                    if oc in hte_data.columns and oc not in trt.columns:
                        ht = hte_data[hte_data.treated == 1][
                            ["pid", oc]
                        ].drop_duplicates("pid")
                        hc = hte_data[hte_data.treated == 0][
                            ["pid", oc]
                        ].drop_duplicates("pid")
                        trt = trt.merge(ht, on="pid", how="left")
                        ctl = ctl.merge(hc, on="pid", how="left")

            print(f"  Matched: {len(trt)} treated, {len(ctl)} control")
        else:
            # all_pts was None — we already loaded from hte_data above
            if trt is None or ctl is None:
                print(f"  ERROR: could not load any data for {tag}")
                return None

    # ── Merge Table 1 descriptive labs from did_labs_all ──
    labs_path = os.path.join(RESULTS, f"did_labs_all_{tag}.csv")
    if os.path.exists(labs_path):
        labs_all = pd.read_csv(labs_path)
        pid_col_labs = (
            "patientunitstayid"
            if "patientunitstayid" in labs_all.columns
            else "stay_id"
        )
        labs_all["pid"] = labs_all[pid_col_labs].astype(str)
        for lab_name in ["hemoglobin", "wbc", "platelets", "albumin"]:
            sub = labs_all[labs_all.lab_name == lab_name].copy()
            if len(sub) == 0:
                continue
            # First postop value per patient
            first_val = (
                sub.sort_values("offset_h")
                .groupby("pid")["value"]
                .first()
                .reset_index()
                .rename(columns={"value": lab_name})
            )
            trt = trt.merge(first_val, on="pid", how="left")
            ctl = ctl.merge(first_val, on="pid", how="left")
        print(f"  Merged Table 1 labs from {labs_path}")

    # ── Merge binary outcomes from did_binary_pairs ──
    bin_path = os.path.join(RESULTS, f"did_binary_pairs_{tag}.csv")
    if os.path.exists(bin_path):
        bo = pd.read_csv(bin_path)
        bo["trt_pid"] = bo["trt_pid"].astype(str)
        bo["ctl_pid"] = bo["ctl_pid"].astype(str)
        # Extract trt outcomes
        trt_cols = [c for c in bo.columns if c.endswith("_trt") and c != "trt_pid"]
        for c in trt_cols:
            oname = c.replace("_trt", "")
            mapping = bo.set_index("trt_pid")[c].to_dict()
            trt[oname] = trt.pid.map(mapping)
        # Extract ctl outcomes
        ctl_cols = [c for c in bo.columns if c.endswith("_ctl") and c != "ctl_pid"]
        for c in ctl_cols:
            oname = c.replace("_ctl", "")
            mapping = bo.set_index("ctl_pid")[c].to_dict()
            ctl[oname] = ctl.pid.map(mapping)
        print(f"  Merged binary outcomes from {bin_path}")

    # Build table rows
    rows = []
    for section_name, var_list in ALL_SECTIONS:
        rows.append(
            {
                "variable": f"**{section_name}**",
                "trt_value": "",
                "ctl_value": "",
                "smd": "",
            }
        )
        for varname, display, vtype in var_list:
            if varname not in trt.columns and varname not in ctl.columns:
                # Try alternative names
                alt_names = {
                    "last_calcium": ["first_calcium", "calcium"],
                    "last_lactate": ["first_lactate", "lactate"],
                    "last_lactate_missing": [
                        "lactate_missing",
                        "first_lactate_missing",
                    ],
                    "last_heartrate": ["first_heartrate", "heartrate"],
                    "last_magnesium": [
                        "first_mg_value",
                        "magnesium",
                        "first_magnesium",
                    ],
                    "last_potassium": ["first_potassium", "potassium"],
                    "cr_pre": ["cr_pre", "baseline_cr"],
                }
                found = False
                for alt in alt_names.get(varname, []):
                    if alt in trt.columns:
                        varname = alt
                        found = True
                        break
                if not found:
                    rows.append(
                        {
                            "variable": f"  {display}",
                            "trt_value": "—",
                            "ctl_value": "—",
                            "smd": "—",
                        }
                    )
                    continue

            x1 = (
                trt[varname].astype(float)
                if varname in trt.columns
                else pd.Series(dtype=float)
            )
            x0 = (
                ctl[varname].astype(float)
                if varname in ctl.columns
                else pd.Series(dtype=float)
            )

            if vtype == "continuous":
                trt_val = format_continuous(x1)
                ctl_val = format_continuous(x0)
            else:  # binary
                trt_val = format_binary(x1.dropna(), len(trt))
                ctl_val = format_binary(x0.dropna(), len(ctl))

            smd = compute_smd(x1, x0, vtype)
            smd_str = f"{smd:.3f}" if not np.isnan(smd) else "—"
            flag = " *" if not np.isnan(smd) and smd > 0.10 else ""

            rows.append(
                {
                    "variable": f"  {display}",
                    "trt_value": trt_val,
                    "ctl_value": ctl_val,
                    "smd": smd_str + flag,
                }
            )

    result = pd.DataFrame(rows)
    result.columns = [
        "Variable",
        f"IV Mg (n={len(trt):,})",
        f"Control (n={len(ctl):,})",
        "SMD",
    ]

    # Print
    print(f"\n  {'Variable':<50s} {'IV Mg':>18s} {'Control':>18s} {'SMD':>8s}")
    print("  " + "─" * 96)
    for _, r in result.iterrows():
        print(f"  {r.iloc[0]:<50s} {r.iloc[1]:>18s} {r.iloc[2]:>18s} {r.iloc[3]:>8s}")

    # Save
    out_path = os.path.join(RESULTS, f"table1_{tag}.csv")
    result.to_csv(out_path, index=False)
    print(f"\n  Saved: {out_path}")

    return result


def combine_tables(tables):
    """Combine MIMIC + eICU tables side-by-side."""
    if len(tables) < 2:
        print("  Only one database available, skipping combined table")
        return

    mimic = tables.get("mimic")
    eicu = tables.get("eicu")
    if mimic is None or eicu is None:
        return

    # Rename columns for combined view
    combined = mimic.copy()
    combined.columns = ["Variable", "MIMIC IV Mg", "MIMIC Control", "MIMIC SMD"]

    # Add eICU columns
    combined["eICU IV Mg"] = eicu.iloc[:, 1]
    combined["eICU Control"] = eicu.iloc[:, 2]
    combined["eICU SMD"] = eicu.iloc[:, 3]

    out = os.path.join(RESULTS, "table1_combined.csv")
    combined.to_csv(out, index=False)
    print(f"\n  Combined table saved: {out}")


if __name__ == "__main__":
    print("=" * 60)
    print("gen_table1.py — Table 1: Baseline Characteristics")
    print("=" * 60)

    args = [a.lower() for a in sys.argv[1:]]
    dbs = args if args else ["mimic", "eicu"]

    tables = {}
    for db in dbs:
        t = build_table1(db)
        if t is not None:
            tables[db] = t

    if len(tables) == 2:
        combine_tables(tables)

    print(f"\n{'='*60}")
    print("DONE")
    print(f"{'='*60}")
