#!/usr/bin/env python3
"""
06_replicate_xiong_koh.py — Replicate prior Mg-AKI associations

Purpose:
  Show that Xiong 2023 (higher post-op Mg → more AKI) and Koh 2022
  (lower pre-op Mg → more AKI) are both reproducible in our data,
  AND that eGFR stratification resolves the apparent contradiction.

Reads (patient-level, never committed):
  results/did_all_{db}.csv       — demographics, comorbidities, first_cr, egfr
  results/did_cr_all_{db}.csv    — creatinine measurements (for AKI computation)
  results/did_labs_all_{db}.csv  — Mg measurements

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


# ═══════════════════════════════════════════════════════════════════
#  COMPUTE AKI FROM CREATININE
# ═══════════════════════════════════════════════════════════════════
def compute_aki(did_all, did_cr, db):
    """
    Compute 7-day AKI (KDIGO ≥1) for each patient from creatinine
    trajectories. Uses first_cr from did_all as baseline.

    KDIGO criteria:
      - Absolute: any Cr ≥ baseline + 0.3 within 48h of ICU admission
      - Relative: any Cr ≥ 1.5 × baseline within 168h (7d) of admission

    Returns: Series indexed by pid with 0/1 AKI flag.
    """
    pid_col = "stay_id" if "stay_id" in did_cr.columns else "pid"
    if pid_col == "stay_id":
        did_cr = did_cr.rename(columns={"stay_id": "pid"})
    if "patientunitstayid" in did_cr.columns:
        did_cr = did_cr.rename(columns={"patientunitstayid": "pid"})

    # Map baseline creatinine
    baseline = did_all.set_index("pid")["first_cr"].to_dict()

    # Filter: only first 7 days (168h)
    cr = did_cr[did_cr.offset_h.between(0, 168)].copy()
    cr["baseline"] = cr.pid.map(baseline)
    cr = cr.dropna(subset=["baseline", "labresult"])

    # Absolute criterion: ΔCr ≥ 0.3 within 48h
    cr_48 = cr[cr.offset_h <= 48]
    abs_hit = cr_48.groupby("pid").apply(
        lambda g: int((g.labresult >= g.baseline + 0.3).any())
    )

    # Relative criterion: Cr ≥ 1.5× baseline within 168h
    rel_hit = cr.groupby("pid").apply(
        lambda g: int((g.labresult >= g.baseline * 1.5).any())
    )

    # Combine
    aki = (
        abs_hit.reindex(did_all.pid, fill_value=0)
        | rel_hit.reindex(did_all.pid, fill_value=0)
    ).astype(int)
    aki.name = "aki_7d"

    n_aki = aki.sum()
    rate = n_aki / len(aki)
    print(f"  AKI computed: {n_aki:,}/{len(aki):,} = {rate:.3f}")
    return aki


# ═══════════════════════════════════════════════════════════════════
#  EXTRACT Mg VALUES FROM LABS
# ═══════════════════════════════════════════════════════════════════
def extract_mg(did_all, labs_path, db):
    """
    Extract first and last post-admission Mg measurement per patient.
    Returns DataFrame with pid, first_mg, last_mg.
    """
    if not os.path.exists(labs_path):
        print(f"  ✗ Labs file not found: {labs_path}")
        return None

    labs = pd.read_csv(labs_path)
    pid_col = (
        "patientunitstayid"
        if "patientunitstayid" in labs.columns
        else "stay_id" if "stay_id" in labs.columns else "pid"
    )

    mg = labs[labs.lab_name == "magnesium"].copy()
    if len(mg) == 0:
        print("  ✗ No magnesium measurements in labs")
        return None

    mg = mg.rename(columns={pid_col: "pid"})
    mg = mg[mg.offset_h >= 0].sort_values("offset_h")

    # First Mg (closest to Koh's "pre-op" — earliest ICU measurement)
    first_mg = (
        mg.groupby("pid").first()[["value"]].rename(columns={"value": "first_mg"})
    )

    # Last Mg (closest to treatment decision — Xiong-like)
    last_mg = mg.groupby("pid").last()[["value"]].rename(columns={"value": "last_mg"})

    result = first_mg.join(last_mg, how="outer")
    print(
        f"  Mg extracted: first={result.first_mg.notna().sum():,}, "
        f"last={result.last_mg.notna().sum():,}"
    )
    return result


# ═══════════════════════════════════════════════════════════════════
#  CROSS-TAB: AKI BY MG QUARTILE × eGFR
# ═══════════════════════════════════════════════════════════════════
def crosstab(df, mg_col, outcome="aki_7d"):
    """AKI rate by Mg quartile × eGFR stratum."""
    sub = df.dropna(subset=[mg_col, outcome, "egfr"]).copy()
    if len(sub) < 100:
        print(f"    Skipping {mg_col}: only {len(sub)} complete cases")
        return pd.DataFrame()

    try:
        sub["mg_q"], bins = pd.qcut(
            sub[mg_col], 4, labels=["Q1", "Q2", "Q3", "Q4"], retbins=True
        )
    except ValueError:
        sub["mg_q"], bins = pd.qcut(
            sub[mg_col].rank(method="first"),
            4,
            labels=["Q1", "Q2", "Q3", "Q4"],
            retbins=True,
        )

    egfr_cuts = {
        "Overall": sub.index,
        "eGFR>=60": sub[sub.egfr >= 60].index,
        "eGFR_45-59": sub[sub.egfr.between(45, 59.999)].index,
        "eGFR<45": sub[sub.egfr < 45].index,
    }

    rows = []
    for stratum, idx in egfr_cuts.items():
        s = sub.loc[idx]
        for q in ["Q1", "Q2", "Q3", "Q4"]:
            qs = s[s.mg_q == q]
            n = len(qs)
            events = int(qs[outcome].sum()) if n > 0 else 0
            rate = events / n if n > 0 else np.nan
            # Get Mg range for this quartile
            qvals = sub[sub.mg_q == q][mg_col]
            mg_range = f"{qvals.min():.2f}-{qvals.max():.2f}" if len(qvals) > 0 else ""
            rows.append(
                {
                    "mg_quartile": q,
                    "mg_range": mg_range,
                    "egfr_stratum": stratum,
                    "n": n,
                    "events": events,
                    "rate": round(rate, 4),
                }
            )

    return pd.DataFrame(rows)


# ═══════════════════════════════════════════════════════════════════
#  LOGISTIC REGRESSION
# ═══════════════════════════════════════════════════════════════════
def run_logistic(df, mg_col, outcome, covariates, label):
    """Multivariable logistic: outcome ~ mg + covariates. Returns dict."""
    try:
        import statsmodels.api as sm
    except ImportError:
        print("  ✗ statsmodels not available")
        return None

    cols = [mg_col, outcome] + covariates
    sub = df.dropna(subset=cols).copy()
    if len(sub) < 100:
        return {
            "label": label,
            "n": len(sub),
            "or": np.nan,
            "p": np.nan,
            "note": f"n={len(sub)} too small",
        }

    X = sm.add_constant(sub[[mg_col] + covariates].astype(float))
    y = sub[outcome].astype(float)

    try:
        fit = sm.Logit(y, X).fit(disp=0, maxiter=100)
        b = fit.params[mg_col]
        se = fit.bse[mg_col]
        return {
            "label": label,
            "n": len(sub),
            "or": round(np.exp(b), 4),
            "ci_lo": round(np.exp(b - 1.96 * se), 4),
            "ci_hi": round(np.exp(b + 1.96 * se), 4),
            "p": round(fit.pvalues[mg_col], 6),
        }
    except Exception as e:
        return {
            "label": label,
            "n": len(sub),
            "or": np.nan,
            "p": np.nan,
            "note": str(e)[:80],
        }


# ═══════════════════════════════════════════════════════════════════
#  MAIN PER-DB
# ═══════════════════════════════════════════════════════════════════
def run_db(db):
    print(f"\n{'=' * 60}\n  {db.upper()}\n{'=' * 60}")

    # ── Load data ──
    all_path = os.path.join(RESULTS, f"did_all_{db}.csv")
    cr_path = os.path.join(RESULTS, f"did_cr_all_{db}.csv")
    labs_path = os.path.join(RESULTS, f"did_labs_all_{db}.csv")

    for p, name in [(all_path, "did_all"), (cr_path, "did_cr_all")]:
        if not os.path.exists(p):
            print(f"  ✗ {name} not found at {p}")
            return

    did_all = pd.read_csv(all_path)
    did_cr = pd.read_csv(cr_path)
    print(f"  did_all: {len(did_all):,} patients")
    print(f"  did_cr:  {len(did_cr):,} measurements")

    # Normalise cr column name
    if "labresult" not in did_cr.columns:
        for alt in ["value", "creatinine", "cr"]:
            if alt in did_cr.columns:
                did_cr = did_cr.rename(columns={alt: "labresult"})
                break

    # ── Compute AKI ──
    aki = compute_aki(did_all, did_cr, db)
    did_all = did_all.set_index("pid")
    did_all["aki_7d"] = aki
    did_all = did_all.reset_index()

    # ── Extract Mg ──
    mg_df = extract_mg(did_all, labs_path, db)
    if mg_df is not None:
        did_all = did_all.set_index("pid").join(mg_df, how="left").reset_index()

    # Determine which Mg columns are available
    mg_cols = []
    for c, label in [
        ("first_mg", "First post-admission Mg"),
        ("last_mg", "Last Mg before discharge/treatment"),
    ]:
        if c in did_all.columns and did_all[c].notna().sum() > 100:
            mg_cols.append((c, label))
            desc = did_all[c].describe()
            print(
                f"  {label}: median={desc['50%']:.2f}, "
                f"IQR=[{desc['25%']:.2f}-{desc['75%']:.2f}], "
                f"n={int(desc['count']):,}"
            )

    if not mg_cols:
        print("  ✗ No Mg variable available")
        return

    # ── Covariates ──
    covar_pool = [
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
    covariates = [c for c in covar_pool if c in did_all.columns]

    # ═══════════════════════════════════════════════════════════
    #  1. CROSS-TABS
    # ═══════════════════════════════════════════════════════════
    all_xtab = []
    for mg_col, mg_label in mg_cols:
        print(f"\n  --- Cross-tab: {mg_label} ---")
        xt = crosstab(did_all, mg_col)
        if len(xt) > 0:
            xt["mg_variable"] = mg_label
            xt["db"] = db.upper()
            all_xtab.append(xt)
            # Print overall quartile rates
            ov = xt[xt.egfr_stratum == "Overall"]
            for _, r in ov.iterrows():
                print(
                    f"    {r.mg_quartile} ({r.mg_range}): "
                    f"n={r.n:,}, AKI={r.events}, rate={r.rate:.3f}"
                )

    if all_xtab:
        out = pd.concat(all_xtab, ignore_index=True)
        path = os.path.join(RESULTS, f"replication_mg_quartile_{db}.csv")
        out.to_csv(path, index=False)
        print(f"\n  ✓ {path}")

    # ═══════════════════════════════════════════════════════════
    #  2. REGRESSIONS
    # ═══════════════════════════════════════════════════════════
    reg_rows = []
    for mg_col, mg_label in mg_cols:
        print(f"\n  --- Regression: {mg_label} ---")

        # Overall (Xiong-style)
        r = run_logistic(did_all, mg_col, "aki_7d", covariates, f"{mg_label} | Overall")
        if r:
            r.update({"db": db.upper(), "mg_variable": mg_label, "stratum": "Overall"})
            reg_rows.append(r)
            print(
                f"    Overall: OR={r['or']} ({r.get('ci_lo', '')}-"
                f"{r.get('ci_hi', '')}), P={r['p']}"
            )

        # By eGFR stratum
        for sname, mask in [
            ("eGFR>=60", did_all.egfr >= 60),
            ("eGFR_45-59", did_all.egfr.between(45, 59.999)),
            ("eGFR<45", did_all.egfr < 45),
        ]:
            sub = did_all[mask]
            r = run_logistic(sub, mg_col, "aki_7d", covariates, f"{mg_label} | {sname}")
            if r:
                r.update({"db": db.upper(), "mg_variable": mg_label, "stratum": sname})
                reg_rows.append(r)
                print(
                    f"    {sname}: OR={r['or']} ({r.get('ci_lo', '')}-"
                    f"{r.get('ci_hi', '')}), P={r['p']}"
                )

        # Interaction: Mg × eGFR<60
        sub = did_all.dropna(subset=[mg_col, "aki_7d", "egfr"] + covariates).copy()
        sub["egfr_low"] = (sub.egfr < 60).astype(float)
        sub["mg_x_egfr_low"] = sub[mg_col] * sub["egfr_low"]
        int_covars = [c for c in covariates if c != "egfr"] + ["egfr_low"]
        r = run_logistic(
            sub,
            "mg_x_egfr_low",
            "aki_7d",
            [mg_col] + int_covars,
            f"{mg_label} | Interaction Mg×eGFR<60",
        )
        if r:
            r.update(
                {
                    "db": db.upper(),
                    "mg_variable": mg_label,
                    "stratum": "Interaction Mg×eGFR<60",
                }
            )
            reg_rows.append(r)
            print(f"    Interaction Mg×eGFR<60: OR={r['or']}, P={r['p']}")

    if reg_rows:
        out = pd.DataFrame(reg_rows)
        path = os.path.join(RESULTS, f"replication_regression_{db}.csv")
        out.to_csv(path, index=False)
        print(f"\n  ✓ {path}")


# ═══════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("=" * 60)
    print("05_replicate_xiong_koh.py")
    print("Replicating Xiong 2023 (Ren Fail) & Koh 2022 (AJKD)")
    print("=" * 60)

    targets = [a.lower() for a in sys.argv[1:]] or DBS
    for db in targets:
        if db in DBS:
            run_db(db)

    print(f"\n{'=' * 60}")
    print("DONE — outputs are aggregate only, safe to commit.")
    print("=" * 60)
