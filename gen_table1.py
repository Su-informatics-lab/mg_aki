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
    ("last_magnesium", "Serum magnesium, median (IQR), mg/dL", "median_iqr"),
    ("last_potassium", "Serum potassium, mean (SD), mEq/L", "continuous"),
]
# Mg and K moved to LAB_VARS above; kept for backward compatibility
SUPP_VARS = []
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

# ── Plausible lab ranges (same logic as ETL HR filter) ──
LAB_RANGES = {
    "magnesium": (0.5, 10.0),
    "potassium": (1.0, 12.0),
    "calcium": (4.0, 16.0),
    "lactate": (0.1, 30.0),
    "heartrate": (20, 250),
    "hemoglobin": (3.0, 25.0),
    "wbc": (0.1, 100.0),
    "platelets": (5.0, 1500.0),
    "albumin": (0.5, 6.0),
}


def compute_smd(x1, x0, vtype="continuous"):
    x1, x0 = x1.dropna(), x0.dropna()
    if len(x1) < 5 or len(x0) < 5:
        return np.nan
    m1, m0 = x1.mean(), x0.mean()
    if vtype == "binary":
        sp = np.sqrt((m1 * (1 - m1) + m0 * (1 - m0)) / 2)
    else:
        sp = np.sqrt((x1.var() + x0.var()) / 2)
    return abs(m1 - m0) / sp if sp > 1e-10 else 0.0


def format_continuous(series):
    m, s = series.mean(), series.std()
    if pd.isna(m):
        return "—"
    return f"{m:.1f} ({s:.1f})" if abs(m) >= 10 else f"{m:.2f} ({s:.2f})"


def format_median_iqr(series):
    s = series.dropna()
    if len(s) < 5:
        return "—"
    q25, med, q75 = s.quantile([0.25, 0.50, 0.75])
    return f"{med:.2f} ({q25:.2f}–{q75:.2f})"


def format_binary(series, n_total):
    n_events = int(series.sum())
    pct = 100 * n_events / n_total if n_total > 0 else 0
    return f"{n_events:,} ({pct:.1f})"


def merge_labs(trt, ctl, tag):
    """
    Merge labs from did_labs_all into trt/ctl DataFrames.
    Two kinds:
      1. PS-model labs (Ca, lactate, HR, Mg, K): T0-relative
         Treated: last value where offset_h < mg_offset_h
         Control: mg_offset_h is NaN → all values qualify → take last
         (Same logic as 02_psm.R extract_labs())
      2. Descriptive labs (Hb, WBC, plt, albumin): first postop value
    Both apply plausible range filtering from LAB_RANGES.
    """
    labs_path = os.path.join(RESULTS, f"did_labs_all_{tag}.csv")
    if not os.path.exists(labs_path):
        print(f"  WARN: {labs_path} not found, skipping lab merge")
        return trt, ctl

    labs_all = pd.read_csv(labs_path)
    pid_col_labs = (
        "patientunitstayid" if "patientunitstayid" in labs_all.columns else "stay_id"
    )
    labs_all["pid"] = labs_all[pid_col_labs].astype(str)

    # Get mg_offset_h for T0-relative computation
    all_path = os.path.join(RESULTS, f"did_all_{tag}.csv")
    mg_off_map = {}
    if os.path.exists(all_path):
        _da = pd.read_csv(all_path)
        _da["pid"] = _da["pid"].astype(str)
        mg_off_map = _da.set_index("pid")["mg_offset_h"].to_dict()

    # 1. PS-model labs: T0-relative (last before T0)
    for lab_name in ["magnesium", "potassium", "calcium", "lactate", "heartrate"]:
        col_name = f"last_{lab_name}"
        sub = labs_all[labs_all.lab_name == lab_name].copy()
        if len(sub) == 0:
            continue
        lo, hi = LAB_RANGES.get(lab_name, (None, None))
        if lo is not None:
            sub = sub[sub.value.between(lo, hi)]
        sub["mg_off"] = sub.pid.map(mg_off_map)
        # Treated: offset_h < mg_offset_h; Control: mg_off is NaN → keep all
        sub = sub[
            (sub.offset_h >= 0) & (sub.mg_off.isna() | (sub.offset_h < sub.mg_off))
        ]
        # Take LAST value (closest to T0)
        last_val = (
            sub.sort_values("offset_h", ascending=False)
            .groupby("pid")["value"]
            .first()
            .reset_index()
            .rename(columns={"value": col_name})
        )
        trt = trt.merge(last_val, on="pid", how="left")
        ctl = ctl.merge(last_val, on="pid", how="left")

    # Lactate missing indicator (matches 02_psm.R extract_labs logic)
    trt["last_lactate_missing"] = trt["last_lactate"].isna().astype(int)
    ctl["last_lactate_missing"] = ctl["last_lactate"].isna().astype(int)

    # 2. Descriptive labs: first postop value with range filter
    for lab_name in ["hemoglobin", "wbc", "platelets", "albumin"]:
        sub = labs_all[labs_all.lab_name == lab_name].copy()
        if len(sub) == 0:
            continue
        lo, hi = LAB_RANGES.get(lab_name, (None, None))
        if lo is not None:
            n_before = len(sub)
            sub = sub[sub.value.between(lo, hi)]
            n_dropped = n_before - len(sub)
            if n_dropped > 0:
                print(
                    f"    {lab_name}: dropped {n_dropped:,} outliers "
                    f"outside [{lo}-{hi}]"
                )
        first_val = (
            sub.sort_values("offset_h")
            .groupby("pid")["value"]
            .first()
            .reset_index()
            .rename(columns={"value": lab_name})
        )
        trt = trt.merge(first_val, on="pid", how="left")
        ctl = ctl.merge(first_val, on="pid", how="left")

    print(f"  Merged labs (PS-model + descriptive) from {labs_path}")
    return trt, ctl


def build_table1(db_tag):
    tag = db_tag.lower()
    db_label = "MIMIC-IV" if tag == "mimic" else "eICU-CRD"
    print(f"\n{'='*60}")
    print(f"Table 1: {db_label}")
    print(f"{'='*60}")

    pairs_path = os.path.join(RESULTS, f"did_pairs_primary_yet_untreated_{tag}.csv")
    all_path = os.path.join(RESULTS, f"did_all_{tag}.csv")
    hte_path = os.path.join(RESULTS, f"did_hte_data_{tag}.csv")

    trt = ctl = None

    if not os.path.exists(pairs_path):
        if not os.path.exists(hte_path):
            print(f"  ERROR: no data files found for {tag}")
            return None
        df = pd.read_csv(hte_path)
        df["pid"] = df["pid"].astype(str)
        trt = df[df.treated == 1]
        ctl = df[df.treated == 0]
    else:
        pairs = pd.read_csv(pairs_path)
        trt_pid_col = [
            c for c in pairs.columns if "trt" in c.lower() and "pid" in c.lower()
        ]
        ctl_pid_col = [
            c for c in pairs.columns if "ctl" in c.lower() and "pid" in c.lower()
        ]
        trt_pid_col = trt_pid_col[0] if trt_pid_col else pairs.columns[0]
        ctl_pid_col = ctl_pid_col[0] if ctl_pid_col else pairs.columns[1]
        print(
            f"  Loaded {len(pairs)} matched pairs (cols: {trt_pid_col}, {ctl_pid_col})"
        )

        all_pts = None
        if os.path.exists(all_path):
            all_pts = pd.read_csv(all_path)

        if all_pts is not None:
            id_col = next(
                (
                    c
                    for c in ["pid", "stay_id", "patientunitstayid"]
                    if c in all_pts.columns
                ),
                None,
            )
            if id_col is None:
                print(f"  ERROR: no ID column found")
                return None
            all_pts["pid"] = all_pts[id_col].astype(str)
            print(f"  ID column: {id_col} ({len(all_pts)} patients)")

            trt_pids = set(pairs[trt_pid_col].astype(str))
            ctl_pids = set(pairs[ctl_pid_col].astype(str))
            trt = all_pts[all_pts.pid.isin(trt_pids)].copy()
            ctl = all_pts[all_pts.pid.isin(ctl_pids)].copy()

            # Merge outcomes from hte_data if available
            if os.path.exists(hte_path):
                hte_data = pd.read_csv(hte_path)
                hte_data["pid"] = hte_data["pid"].astype(str)
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

    if trt is None or ctl is None:
        print(f"  ERROR: could not load data for {tag}")
        return None

    # ── Merge labs ──
    trt, ctl = merge_labs(trt, ctl, tag)

    # ── Merge binary outcomes from did_binary_pairs ──
    bin_path = os.path.join(RESULTS, f"did_binary_pairs_{tag}.csv")
    if os.path.exists(bin_path):
        bo = pd.read_csv(bin_path)
        bo["trt_pid"] = bo["trt_pid"].astype(str)
        bo["ctl_pid"] = bo["ctl_pid"].astype(str)
        for c in [c for c in bo.columns if c.endswith("_trt") and c != "trt_pid"]:
            trt[c.replace("_trt", "")] = trt.pid.map(
                bo.set_index("trt_pid")[c].to_dict()
            )
        for c in [c for c in bo.columns if c.endswith("_ctl") and c != "ctl_pid"]:
            ctl[c.replace("_ctl", "")] = ctl.pid.map(
                bo.set_index("ctl_pid")[c].to_dict()
            )
        print(f"  Merged binary outcomes from {bin_path}")

    # ── Build table rows ──
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
                    "cr_pre": ["cr_pre", "baseline_cr", "first_cr"],
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
            if vtype == "median_iqr":
                trt_val = format_median_iqr(x1)
                ctl_val = format_median_iqr(x0)
            elif vtype == "continuous":
                trt_val = format_continuous(x1)
                ctl_val = format_continuous(x0)
            else:
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

    print(f"\n  {'Variable':<50s} {'IV Mg':>18s} {'Control':>18s} {'SMD':>8s}")
    print("  " + "─" * 96)
    for _, r in result.iterrows():
        print(f"  {r.iloc[0]:<50s} {r.iloc[1]:>18s} {r.iloc[2]:>18s} {r.iloc[3]:>8s}")

    out_path = os.path.join(RESULTS, f"table1_{tag}.csv")
    result.to_csv(out_path, index=False)
    print(f"\n  Saved: {out_path}")
    return result


def combine_tables(tables):
    mimic, eicu = tables.get("mimic"), tables.get("eicu")
    if mimic is None or eicu is None:
        return
    combined = mimic.copy()
    combined.columns = ["Variable", "MIMIC IV Mg", "MIMIC Control", "MIMIC SMD"]
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
    print(f"\n{'='*60}\nDONE\n{'='*60}")
