#!/usr/bin/env python3
"""
05_replicate_xiong_koh.py — Replicate prior Mg-AKI associations

Purpose:
  Show that Xiong 2023 (higher post-op Mg → more AKI) and Koh 2022
  (lower pre-op Mg → more AKI) are both reproducible in our data,
  AND that eGFR stratification resolves the apparent contradiction.

Reads (patient-level, never committed):
  results/did_all_{db}.csv         — full cohort with covariates
  results/did_labs_all_{db}.csv    — individual lab measurements (optional)

Outputs (aggregate only, safe to commit):
  results/replication_mg_quartile_{db}.csv  — AKI rate by Mg quartile × eGFR
  results/replication_regression_{db}.csv   — logistic regression ORs

Usage on Tempest:
  module purge && module load Python/3.10.8-GCCcore-12.2.0
  source ~/alcrx/.venv/bin/activate
  cd ~/mg_aki
  python 05_replicate_xiong_koh.py
"""

import os
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

RESULTS = os.path.expanduser("~/mg_aki/results")
DBS = ["mimic", "eicu"]


def load_cohort(db):
    """Load full cohort with covariates and outcomes."""
    path = os.path.join(RESULTS, f"did_all_{db}.csv")
    if not os.path.exists(path):
        print(f"  ✗ {path} not found — skipping {db}")
        return None
    df = pd.read_csv(path)
    print(
        f"  Loaded {db}: {len(df):,} patients, "
        f"{df.columns.tolist()[:5]}... ({len(df.columns)} cols)"
    )
    return df


def compute_first_mg(db, cohort):
    """
    Compute first post-admission Mg from did_labs_all if available.
    Falls back to last_magnesium from cohort if labs file missing.
    """
    labs_path = os.path.join(RESULTS, f"did_labs_all_{db}.csv")
    if not os.path.exists(labs_path):
        print(f"  did_labs_all not found — using last_magnesium only")
        return None

    labs = pd.read_csv(labs_path)
    pid_col = "patientunitstayid" if "patientunitstayid" in labs.columns else "stay_id"
    mg_labs = labs[labs.lab_name == "magnesium"].copy()
    if len(mg_labs) == 0:
        print(f"  No magnesium measurements in labs — using last_magnesium")
        return None

    # First Mg measurement after ICU admission (offset_h ≥ 0)
    mg_labs = mg_labs[mg_labs.offset_h >= 0].sort_values("offset_h")
    first_mg = mg_labs.groupby(pid_col).first().reset_index()
    first_mg = first_mg[[pid_col, "value"]].rename(
        columns={pid_col: "pid", "value": "first_magnesium"}
    )
    print(f"  Computed first_magnesium for {len(first_mg):,} patients")
    return first_mg


def assign_quartiles(series, prefix="Q"):
    """Assign quartile labels. Returns (labels, cut_points)."""
    try:
        labels, bins = pd.qcut(series, 4, labels=False, retbins=True)
        q_labels = []
        for i in range(4):
            lo = f"{bins[i]:.2f}"
            hi = f"{bins[i+1]:.2f}"
            q_labels.append(f"{prefix}{i+1} [{lo}-{hi}]")
        return labels.map(lambda x: q_labels[int(x)] if pd.notna(x) else np.nan), bins
    except Exception as e:
        print(f"    Quartile computation failed: {e}")
        return pd.Series([np.nan] * len(series)), None


def crosstab_rates(df, mg_col, outcome, egfr_strata):
    """
    Compute AKI rates by Mg quartile × eGFR stratum.
    Returns a DataFrame with one row per (quartile, eGFR_stratum).
    """
    rows = []
    df = df.dropna(subset=[mg_col, outcome])
    q_labels, bins = assign_quartiles(df[mg_col])
    df = df.copy()
    df["mg_quartile"] = q_labels

    for stratum_name, mask in egfr_strata.items():
        sub = df[mask(df)] if callable(mask) else df
        for q in sorted(df["mg_quartile"].dropna().unique()):
            qsub = sub[sub.mg_quartile == q]
            n = len(qsub)
            events = qsub[outcome].sum() if n > 0 else 0
            rate = events / n if n > 0 else np.nan
            rows.append(
                {
                    "mg_quartile": q,
                    "egfr_stratum": stratum_name,
                    "n": n,
                    "events": int(events),
                    "rate": round(rate, 4) if pd.notna(rate) else np.nan,
                }
            )

    return pd.DataFrame(rows)


def logistic_regression(df, mg_col, outcome, covariates, label=""):
    """
    Multivariable logistic regression: outcome ~ mg_col + covariates.
    Returns a dict with OR, CI, P for the Mg coefficient.
    """
    try:
        import statsmodels.api as sm
    except ImportError:
        print("  statsmodels not available — skipping regression")
        return None

    sub = df.dropna(subset=[mg_col, outcome] + covariates).copy()
    if len(sub) < 100:
        return {
            "label": label,
            "n": len(sub),
            "or": np.nan,
            "ci_lo": np.nan,
            "ci_hi": np.nan,
            "p": np.nan,
            "note": "insufficient data",
        }

    X = sub[[mg_col] + covariates].astype(float)
    X = sm.add_constant(X)
    y = sub[outcome].astype(float)

    try:
        model = sm.Logit(y, X).fit(disp=0, maxiter=100)
        coef = model.params[mg_col]
        se = model.bse[mg_col]
        p = model.pvalues[mg_col]
        or_val = np.exp(coef)
        ci_lo = np.exp(coef - 1.96 * se)
        ci_hi = np.exp(coef + 1.96 * se)
        return {
            "label": label,
            "n": len(sub),
            "or": round(or_val, 4),
            "ci_lo": round(ci_lo, 4),
            "ci_hi": round(ci_hi, 4),
            "p": round(p, 6),
        }
    except Exception as e:
        return {
            "label": label,
            "n": len(sub),
            "or": np.nan,
            "ci_lo": np.nan,
            "ci_hi": np.nan,
            "p": np.nan,
            "note": str(e),
        }


def run_db(db):
    """Run full replication for one database."""
    print(f"\n{'=' * 60}")
    print(f"  {db.upper()}")
    print(f"{'=' * 60}")

    df = load_cohort(db)
    if df is None:
        return

    # Identify pid column
    pid_col = (
        "pid"
        if "pid" in df.columns
        else ("patientunitstayid" if "patientunitstayid" in df.columns else "stay_id")
    )

    # ── Mg variables ──
    # last_magnesium: closest to T0 (Xiong-like: post-op, close to treatment)
    mg_cols = {}
    if "last_magnesium" in df.columns:
        mg_cols["last_magnesium"] = "Post-op Mg (closest to T0)"
        n_avail = df["last_magnesium"].notna().sum()
        print(f"  last_magnesium available: {n_avail:,}/{len(df):,}")

    # Try to compute first_magnesium from labs
    first_mg = compute_first_mg(db, df)
    if first_mg is not None:
        df["_pid_merge"] = df[pid_col].astype(str)
        first_mg["pid"] = first_mg["pid"].astype(str)
        df = df.merge(
            first_mg,
            left_on="_pid_merge",
            right_on="pid",
            how="left",
            suffixes=("", "_fm"),
        )
        mg_cols["first_magnesium"] = "First post-admission Mg (Koh-like)"
        n_avail = df["first_magnesium"].notna().sum()
        print(f"  first_magnesium available: {n_avail:,}/{len(df):,}")

    if not mg_cols:
        print("  ✗ No Mg variable available — skipping")
        return

    # ── AKI outcome ──
    outcome = "aki1_7d"
    if outcome not in df.columns:
        # Try alternative names
        for alt in ["aki_7d", "aki1_7d", "kdigo1_7d"]:
            if alt in df.columns:
                outcome = alt
                break
    print(f"  Outcome: {outcome} (rate = {df[outcome].mean():.3f})")

    # ── Covariates for regression (Xiong-style) ──
    covar_candidates = [
        "age",
        "is_female",
        "bmi",
        "egfr",
        "heart_failure",
        "hypertension",
        "diabetes",
        "ckd",
        "copd",
        "pvd",
        "stroke",
        "liver_disease",
        "surg_cabg",
        "surg_valve",
        "surg_combined",
    ]
    covariates = [c for c in covar_candidates if c in df.columns]
    print(f"  Covariates for regression: {len(covariates)}")

    # ── eGFR strata for cross-tabulation ──
    egfr_strata = {
        "Overall": lambda d: pd.Series([True] * len(d), index=d.index),
        "eGFR >= 60": lambda d: d["egfr"] >= 60,
        "eGFR 45-59": lambda d: (d["egfr"] >= 45) & (d["egfr"] < 60),
        "eGFR < 45": lambda d: d["egfr"] < 45,
    }

    # ═══════════════════════════════════════════════════════════
    # 1. CROSS-TABULATION: AKI rate by Mg quartile × eGFR
    # ═══════════════════════════════════════════════════════════
    all_xtab = []
    for mg_col, mg_label in mg_cols.items():
        print(f"\n  Cross-tab: {mg_label}")
        xtab = crosstab_rates(df, mg_col, outcome, egfr_strata)
        xtab["mg_variable"] = mg_label
        xtab["db"] = db.upper()
        all_xtab.append(xtab)
        # Print summary
        overall = xtab[xtab.egfr_stratum == "Overall"]
        for _, row in overall.iterrows():
            print(
                f"    {row.mg_quartile}: n={row.n:,}, "
                f"AKI={row.events}, rate={row.rate:.3f}"
            )

    xtab_df = pd.concat(all_xtab, ignore_index=True)
    out_path = os.path.join(RESULTS, f"replication_mg_quartile_{db}.csv")
    xtab_df.to_csv(out_path, index=False)
    print(f"\n  ✓ Saved {out_path}")

    # ═══════════════════════════════════════════════════════════
    # 2. LOGISTIC REGRESSIONS
    # ═══════════════════════════════════════════════════════════
    reg_rows = []
    for mg_col, mg_label in mg_cols.items():
        print(f"\n  Regression: {mg_label}")

        # 2a. Overall (Xiong-style): AKI ~ Mg + covariates
        res = logistic_regression(
            df, mg_col, outcome, covariates, label=f"{mg_label} — Overall"
        )
        if res:
            res["db"] = db.upper()
            res["mg_variable"] = mg_label
            res["stratum"] = "Overall"
            reg_rows.append(res)
            print(
                f"    Overall: OR {res['or']} "
                f"({res['ci_lo']}-{res['ci_hi']}), P={res['p']}"
            )

        # 2b. Stratified by eGFR
        for sname, smask in egfr_strata.items():
            if sname == "Overall":
                continue
            sub = df[smask(df)]
            if len(sub) < 50:
                continue
            res = logistic_regression(
                sub, mg_col, outcome, covariates, label=f"{mg_label} — {sname}"
            )
            if res:
                res["db"] = db.upper()
                res["mg_variable"] = mg_label
                res["stratum"] = sname
                reg_rows.append(res)
                print(
                    f"    {sname}: OR {res['or']} "
                    f"({res['ci_lo']}-{res['ci_hi']}), P={res['p']}"
                )

        # 2c. Interaction model: AKI ~ Mg + Mg×eGFR_low + covariates
        sub = df.dropna(subset=[mg_col, outcome, "egfr"] + covariates).copy()
        sub["egfr_below_60"] = (sub["egfr"] < 60).astype(float)
        sub["mg_x_egfr_low"] = sub[mg_col] * sub["egfr_below_60"]
        interaction_covars = covariates + ["egfr_below_60", "mg_x_egfr_low"]
        # Remove egfr from covariates to avoid collinearity with binary
        interaction_covars_clean = [c for c in interaction_covars if c != "egfr"]
        res = logistic_regression(
            sub,
            mg_col,
            outcome,
            interaction_covars_clean,
            label=f"{mg_label} — Mg main effect " f"(interaction model)",
        )
        if res:
            res["db"] = db.upper()
            res["mg_variable"] = mg_label
            res["stratum"] = "Interaction (Mg main)"
            reg_rows.append(res)

        # Get the interaction term OR
        res_int = logistic_regression(
            sub,
            "mg_x_egfr_low",
            outcome,
            [mg_col] + [c for c in interaction_covars_clean if c != "mg_x_egfr_low"],
            label=f"{mg_label} — Mg×eGFR<60 interaction",
        )
        if res_int:
            res_int["db"] = db.upper()
            res_int["mg_variable"] = mg_label
            res_int["stratum"] = "Interaction (Mg × eGFR<60)"
            reg_rows.append(res_int)
            print(
                f"    Interaction Mg×eGFR<60: OR {res_int['or']} "
                f"({res_int['ci_lo']}-{res_int['ci_hi']}), P={res_int['p']}"
            )

    reg_df = pd.DataFrame(reg_rows)
    out_path = os.path.join(RESULTS, f"replication_regression_{db}.csv")
    reg_df.to_csv(out_path, index=False)
    print(f"  ✓ Saved {out_path}")


# ═══════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("=" * 60)
    print("05_replicate_xiong_koh.py")
    print("Replicating Xiong 2023 (Ren Fail) and Koh 2022 (AJKD)")
    print("=" * 60)

    targets = sys.argv[1:] if len(sys.argv) > 1 else DBS
    for db in targets:
        if db.lower() in DBS:
            run_db(db.lower())
        else:
            print(f"Unknown db: {db}. Use 'mimic' or 'eicu'.")

    print(f"\n{'=' * 60}")
    print("DONE — outputs are aggregate only, safe to commit.")
    print(f"{'=' * 60}")
