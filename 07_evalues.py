#!/usr/bin/env python3
"""
07_evalues.py — E-value sensitivity analysis (VanderWeele & Ding 2017)

Reads:  results/02_results.csv
        results/01_analysis_a_cohort.csv   (for outcome prevalence)
        results/04_mimic_cohort.csv        (for outcome prevalence)
Writes: results/07_evalues.csv

Run: python 07_evalues.py
"""

import math
import os

import numpy as np
import pandas as pd

RESULTS = os.path.expanduser("~/mg_aki/results")


# =====================================================================
# E-value formula (VanderWeele & Ding, Ann Intern Med 2017)
# =====================================================================
def or_to_rr(or_val, p0):
    """Convert OR → RR using baseline risk p0 (Zhang & Yu 1998)."""
    return or_val / (1 - p0 + p0 * or_val)


def e_value(rr):
    """E-value for RR > 1. For protective (RR<1), caller inverts first."""
    if rr <= 1.0:
        return 1.0
    return rr + math.sqrt(rr * (rr - 1))


def compute_evalues(or_est, or_lo, or_hi, p0=None):
    """
    Compute E-values for point estimate and CI limit closest to null.

    For protective effects (OR < 1):
      - CI limit closest to null = upper bound (or_hi)
      - Invert to RR > 1 before computing E-value

    Returns dict with both conservative (OR-as-RR) and adjusted (OR→RR).
    """
    # Determine which CI limit is closest to null
    ci_null = or_hi if or_est < 1 else or_lo

    results = {}

    # Conservative: treat OR as RR
    rr_point = 1.0 / or_est if or_est < 1 else or_est
    rr_ci = 1.0 / ci_null if ci_null < 1 else ci_null
    results["e_point_conservative"] = e_value(rr_point)
    results["e_ci_conservative"] = e_value(rr_ci)

    # Adjusted: convert OR → RR using outcome prevalence
    if p0 is not None and p0 > 0:
        rr_adj_est = or_to_rr(or_est, p0)
        rr_adj_ci = or_to_rr(ci_null, p0)
        # Invert if protective
        rr_adj_point = 1.0 / rr_adj_est if rr_adj_est < 1 else rr_adj_est
        rr_adj_ci_f = 1.0 / rr_adj_ci if rr_adj_ci < 1 else rr_adj_ci
        results["e_point_adjusted"] = e_value(rr_adj_point)
        results["e_ci_adjusted"] = e_value(rr_adj_ci_f)
        results["p0"] = p0
    else:
        results["e_point_adjusted"] = np.nan
        results["e_ci_adjusted"] = np.nan
        results["p0"] = np.nan

    return results


# =====================================================================
# Compute outcome prevalence from cohort CSVs
# =====================================================================
def get_prevalence(cohort_path, outcome_col, subset_col=None, subset_val=None):
    """Get outcome prevalence in control arm, optionally within a subset."""
    try:
        df = pd.read_csv(cohort_path)
        # Standardize treatment column
        if "mg_supplementation" in df.columns:
            df["mg_supp"] = df["mg_supplementation"]
        if subset_col and subset_val:
            df = df[df[subset_col] == subset_val]
        ctrl = df[df["mg_supp"] == 0]
        if outcome_col in ctrl.columns:
            return ctrl[outcome_col].mean()
    except Exception:
        pass
    return None


def main():
    print("=" * 65)
    print("E-VALUE SENSITIVITY ANALYSIS")
    print("VanderWeele & Ding, Ann Intern Med 2017;167:268-274")
    print("=" * 65)

    # ── Load results ─────────────────────────────────────────────
    res = pd.read_csv(os.path.join(RESULTS, "02_results.csv"))
    if "m" in res.columns:
        m_val = res["m"].dropna().min()
        res = res[(res["m"].isna()) | (res["m"] == m_val)]

    pooled = res[res["db"] == "Pooled"].copy()
    or_col = "or_" if "or_" in pooled.columns else "or"

    # ── Outcome prevalence from cohorts ──────────────────────────
    eicu_path = os.path.join(RESULTS, "01_analysis_a_cohort.csv")
    mimic_path = os.path.join(RESULTS, "04_mimic_cohort.csv")

    # All-patient AKI prevalence in controls
    p0_eicu = get_prevalence(eicu_path, "aki_kdigo1")
    p0_mimic = get_prevalence(mimic_path, "aki_kdigo1")

    # AC population AKI prevalence in K-only controls
    p0_ac_eicu = get_prevalence(eicu_path, "aki_kdigo1", "ac_group", "k_only")
    p0_ac_mimic = get_prevalence(mimic_path, "aki_kdigo1", "ac_group", "k_only")

    # Mortality prevalence in controls
    p0_mort_eicu = get_prevalence(eicu_path, "hosp_mortality")
    p0_mort_mimic = get_prevalence(mimic_path, "hosp_mortality")
    if p0_mort_eicu is None:
        p0_mort_eicu = get_prevalence(eicu_path, "hospital_mortality")
    if p0_mort_mimic is None:
        p0_mort_mimic = get_prevalence(mimic_path, "hospital_mortality")

    print("\n── Outcome Prevalence (control arm) ──")
    for label, p in [
        ("eICU all-patient AKI", p0_eicu),
        ("MIMIC all-patient AKI", p0_mimic),
        ("eICU AC (K-only) AKI", p0_ac_eicu),
        ("MIMIC AC (K-only) AKI", p0_ac_mimic),
        ("eICU mortality", p0_mort_eicu),
        ("MIMIC mortality", p0_mort_mimic),
    ]:
        print(f"  {label}: {p:.3f}" if p else f"  {label}: unavailable")

    # ── Pooled prevalence for E-value ────────────────────────────
    def pool_prev(*ps):
        vals = [p for p in ps if p is not None]
        return sum(vals) / len(vals) if vals else None

    p0_ac = pool_prev(p0_ac_eicu, p0_ac_mimic)
    p0_all = pool_prev(p0_eicu, p0_mimic)
    p0_mort = pool_prev(p0_mort_eicu, p0_mort_mimic)

    # ── Define analyses for E-value ──────────────────────────────
    analyses = [
        ("ac_aki1", "Primary: AC AKI (Mg+K vs K-only)", p0_ac),
        ("iptw_aki1", "Sensitivity: IPTW AKI", p0_all),
        ("ow_aki1", "Sensitivity: OW AKI", p0_all),
        ("psm_aki1", "Sensitivity: PSM AKI", p0_all),
        ("ow_mort", "Exploratory: Hospital mortality", p0_mort),
        ("ow_enceph", "Exploratory: Encephalopathy", None),
        ("ow_frac", "Control: Fracture", None),
    ]

    # ── Compute E-values ─────────────────────────────────────────
    print("\n" + "=" * 65)
    print("E-VALUES")
    print("=" * 65)

    rows = []
    for analysis_key, label, p0 in analyses:
        r = pooled[pooled["analysis"] == analysis_key]
        if len(r) == 0:
            continue
        r = r.iloc[0]
        or_est = r[or_col]
        or_lo = r["lo"]
        or_hi = r["hi"]
        p_val = r["p"]

        ev = compute_evalues(or_est, or_lo, or_hi, p0)

        print(f"\n  {label}")
        print(f"    OR = {or_est:.3f} ({or_lo:.3f}–{or_hi:.3f}), P = {p_val:.4f}")
        print(
            f"    E-value (conservative): point = {ev['e_point_conservative']:.2f}, "
            f"CI limit = {ev['e_ci_conservative']:.2f}"
        )
        if not np.isnan(ev["e_point_adjusted"]):
            print(
                f"    E-value (adjusted, p0={p0:.3f}): point = {ev['e_point_adjusted']:.2f}, "
                f"CI limit = {ev['e_ci_adjusted']:.2f}"
            )

        rows.append(
            {
                "analysis": analysis_key,
                "label": label,
                "or": round(or_est, 3),
                "lo": round(or_lo, 3),
                "hi": round(or_hi, 3),
                "p": round(p_val, 4),
                "p0_control": round(p0, 3) if p0 else np.nan,
                "e_point_conservative": round(ev["e_point_conservative"], 2),
                "e_ci_conservative": round(ev["e_ci_conservative"], 2),
                "e_point_adjusted": (
                    round(ev["e_point_adjusted"], 2)
                    if not np.isnan(ev["e_point_adjusted"])
                    else np.nan
                ),
                "e_ci_adjusted": (
                    round(ev["e_ci_adjusted"], 2)
                    if not np.isnan(ev["e_ci_adjusted"])
                    else np.nan
                ),
            }
        )

    # ── Save ─────────────────────────────────────────────────────
    out = pd.DataFrame(rows)
    outpath = os.path.join(RESULTS, "07_evalues.csv")
    out.to_csv(outpath, index=False)
    print(f"\n{'='*65}")
    print(f"Saved: {outpath}")
    print(f"{'='*65}")

    # ── Manuscript sentence ──────────────────────────────────────
    ac = out[out.analysis == "ac_aki1"]
    iptw = out[out.analysis == "iptw_aki1"]
    if len(ac) > 0 and len(iptw) > 0:
        ac = ac.iloc[0]
        iptw = iptw.iloc[0]
        # Use adjusted if available, else conservative
        ac_pt = (
            ac.e_point_adjusted
            if not np.isnan(ac.e_point_adjusted)
            else ac.e_point_conservative
        )
        ac_ci = (
            ac.e_ci_adjusted if not np.isnan(ac.e_ci_adjusted) else ac.e_ci_conservative
        )
        ip_pt = (
            iptw.e_point_adjusted
            if not np.isnan(iptw.e_point_adjusted)
            else iptw.e_point_conservative
        )
        ip_ci = (
            iptw.e_ci_adjusted
            if not np.isnan(iptw.e_ci_adjusted)
            else iptw.e_ci_conservative
        )
        print("\n── Suggested manuscript sentence (Discussion, limitations) ──")
        print(
            f"  The E-value for the primary active-comparator estimate was "
            f"{ac_pt:.2f} (CI limit, {ac_ci:.2f}), and for the all-patient "
            f"IPTW estimate, {ip_pt:.2f} (CI limit, {ip_ci:.2f}), indicating "
            f"that an unmeasured confounder associated with both supplementation "
            f"and AKI by a risk ratio of {ac_pt:.2f}-fold or greater could "
            f"explain away the observed association, but weaker confounding "
            f"could not."
        )


if __name__ == "__main__":
    main()
