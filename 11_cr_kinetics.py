#!/usr/bin/env python3
"""
11_cr_kinetics.py — Post-operative creatinine kinetics (v3)

  Saves 6 individual panel PDFs:
    fig_cr_trajectory_{db}.pdf    — Cr over time (median + IQR)
    fig_cr_ratio_{db}.pdf         — Cr/baseline P75/P90 fan
    fig_cumulative_aki_{db}.pdf   — Cumulative KDIGO AKI (from precomputed)

Run: python 11_cr_kinetics.py
"""

import os
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

mpl.rcParams.update(
    {
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 8,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "axes.linewidth": 0.5,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "legend.frameon": False,
        "legend.fontsize": 6,
        "lines.linewidth": 1.0,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
)

C_BLUE = "#0072B2"
C_VERM = "#D55E00"
C_GRAY = "#999999"
W_SINGLE = 3.504

RESULTS = os.path.expanduser("~/mg_aki/results")
FIGS = os.path.expanduser("~/mg_aki/figs")
os.makedirs(FIGS, exist_ok=True)

T0_HOURS = 6  # time zero


# ── Data loaders ─────────────────────────────────────────────────
def gz(p):
    return p if os.path.exists(p) else p.replace(".csv.gz", ".csv")


def load_eicu():
    from importlib.util import module_from_spec, spec_from_file_location

    _spec = spec_from_file_location(
        "config",
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "00_config.py"),
    )
    cfg = module_from_spec(_spec)
    _spec.loader.exec_module(cfg)

    cohort = pd.read_csv(os.path.join(RESULTS, "01_analysis_a_cohort.csv"))
    pids = set(cohort.patientunitstayid)
    bl = dict(zip(cohort.patientunitstayid, cohort.baseline_cr))
    trt_col = (
        "mg_supplementation" if "mg_supplementation" in cohort.columns else "mg_supp"
    )
    trt = dict(zip(cohort.patientunitstayid, cohort[trt_col]))

    lab = pd.read_csv(gz(os.path.join(cfg.DATA_ROOT, "lab.csv.gz")), low_memory=False)
    lab.columns = lab.columns.str.lower()
    cr = lab[
        lab.patientunitstayid.isin(pids)
        & lab.labname.str.lower().str.contains("creatinine", na=False)
        & lab.labresult.between(0.1, 25.0)
    ].copy()
    cr["hours"] = cr.labresultoffset / 60.0
    cr["bl_cr"] = cr.patientunitstayid.map(bl)
    cr["cr_ratio"] = cr.labresult / cr.bl_cr
    cr["trt"] = cr.patientunitstayid.map(trt)
    cr = cr[(cr.hours >= -6) & (cr.hours <= 168)].copy()
    print(
        f"  eICU: {len(cr):,} Cr measurements, {cr.patientunitstayid.nunique()} patients"
    )
    return cr, cohort, "eICU"


def load_mimic():
    MIMIC_HOSP = os.path.expanduser("~/mg_aki/mimic-iv-3.1/hosp")
    cohort = pd.read_csv(os.path.join(RESULTS, "04_mimic_cohort.csv"))
    cohort["intime"] = pd.to_datetime(cohort["intime"])
    hadms = set(cohort.hadm_id.dropna().astype(int))
    bl = dict(zip(cohort.stay_id, cohort.baseline_cr))
    trt_col = (
        "mg_supplementation" if "mg_supplementation" in cohort.columns else "mg_supp"
    )
    trt = dict(zip(cohort.stay_id, cohort[trt_col]))
    intime_d = dict(zip(cohort.hadm_id.astype(int), cohort.intime))
    stay_d = dict(zip(cohort.hadm_id.astype(int), cohort.stay_id))

    LAB_CR = [50912, 52546]
    chunks = []
    for chunk in pd.read_csv(
        gz(f"{MIMIC_HOSP}/labevents.csv.gz"),
        usecols=["hadm_id", "itemid", "charttime", "valuenum"],
        dtype={"hadm_id": "Int64", "itemid": int},
        chunksize=5_000_000,
    ):
        chunks.append(chunk[chunk.itemid.isin(LAB_CR) & chunk.hadm_id.isin(hadms)])
    cr = pd.concat(chunks, ignore_index=True)
    cr["charttime"] = pd.to_datetime(cr["charttime"])
    cr = cr[cr.valuenum.between(0.1, 25.0)]
    cr["stay_id"] = cr.hadm_id.astype(int).map(stay_d)
    cr["intime"] = cr.hadm_id.astype(int).map(intime_d)
    cr["hours"] = (cr.charttime - cr.intime).dt.total_seconds() / 3600
    cr["bl_cr"] = cr.stay_id.map(bl)
    cr["cr_ratio"] = cr.valuenum / cr.bl_cr
    cr["trt"] = cr.stay_id.map(trt)
    cr = cr.rename(columns={"valuenum": "labresult"})
    cr = cr[(cr.hours >= -6) & (cr.hours <= 168)].copy()
    print(f"  MIMIC: {len(cr):,} Cr measurements, {cr.stay_id.nunique()} patients")
    return cr, cohort, "MIMIC"


# ── Panel functions (each saves its own PDF) ────────────────────


def _t0_annotation(ax, ypos=None):
    """Add T₀ vertical line."""
    ax.axvline(T0_HOURS, color=C_GRAY, linestyle="--", linewidth=0.5, alpha=0.7)
    if ypos is not None:
        ax.text(T0_HOURS + 1, ypos, "T₀", fontsize=6, color=C_GRAY, va="top")


def plot_trajectory(cr, cohort, db):
    """Panel a/d: Cr trajectory (median + IQR) by treatment group."""
    trt_col = (
        "mg_supplementation" if "mg_supplementation" in cohort.columns else "mg_supp"
    )
    n_trt = int(cohort[trt_col].sum())
    n_ctrl = len(cohort) - n_trt

    bins = np.arange(-6, 174, 6)
    cr["tbin"] = pd.cut(cr.hours, bins=bins, labels=(bins[:-1] + bins[1:]) / 2).astype(
        float
    )

    fig, ax = plt.subplots(figsize=(W_SINGLE, 2.6))
    for g, color, label, n in [
        (1, C_VERM, f"Mg supp (n={n_trt})", n_trt),
        (0, C_BLUE, f"No supp (n={n_ctrl})", n_ctrl),
    ]:
        sub = cr[cr.trt == g]
        s = (
            sub.groupby("tbin")["labresult"]
            .agg(
                [
                    "median",
                    "count",
                    lambda x: x.quantile(0.25),
                    lambda x: x.quantile(0.75),
                ]
            )
            .reset_index()
        )
        s.columns = ["t", "med", "n", "q25", "q75"]
        s = s[s.n >= 20]
        ax.fill_between(s.t, s.q25, s.q75, alpha=0.15, color=color)
        ax.plot(s.t, s.med, color=color, linewidth=1.2, label=label)
    _t0_annotation(ax, ax.get_ylim()[1] * 0.97)
    ax.set_xlabel("Hours from ICU admission")
    ax.set_ylabel("Serum creatinine (mg/dL)")
    ax.set_title(
        f"{db}: Creatinine trajectory (median, IQR)", fontsize=7, fontweight="bold"
    )
    ax.legend(loc="upper right", fontsize=5)
    path = os.path.join(FIGS, f"fig_cr_trajectory_{db.lower()}.pdf")
    fig.savefig(path)
    fig.savefig(path.replace(".pdf", ".png"))
    plt.close(fig)
    print(f"  Saved: {path}")


def plot_ratio(cr, db):
    """Panel b/e: Cr/baseline ratio percentile fan (P50 dotted, P75 dashed, P90 solid)."""
    bins = np.arange(0, 174, 6)
    post = cr[(cr.hours > 0) & cr.cr_ratio.notna() & (cr.cr_ratio < 10)].copy()
    post["tbin"] = pd.cut(
        post.hours, bins=bins, labels=(bins[:-1] + bins[1:]) / 2
    ).astype(float)

    fig, ax = plt.subplots(figsize=(W_SINGLE, 2.6))
    for g, color, glabel in [(1, C_VERM, "Mg supp"), (0, C_BLUE, "No supp")]:
        sub = post[post.trt == g]
        s = (
            sub.groupby("tbin")
            .agg(
                p50=("cr_ratio", "median"),
                p75=("cr_ratio", lambda x: x.quantile(0.75)),
                p90=("cr_ratio", lambda x: x.quantile(0.90)),
                n=("cr_ratio", "count"),
            )
            .reset_index()
        )
        s = s[s.n >= 20]
        t = s.tbin
        ax.fill_between(t, s.p75, s.p90, alpha=0.12, color=color)
        ax.plot(
            t, s.p90, color=color, linewidth=1.2, linestyle="-", label=f"{glabel} P90"
        )
        ax.plot(
            t, s.p75, color=color, linewidth=0.8, linestyle="--", label=f"{glabel} P75"
        )
        ax.plot(
            t,
            s.p50,
            color=color,
            linewidth=0.5,
            linestyle=":",
            label=f"{glabel} P50",
            alpha=0.6,
        )

    ax.axhline(1.5, color="red", linestyle=":", linewidth=0.7, alpha=0.6)
    ax.text(155, 1.52, "KDIGO 1.5×", fontsize=5, color="red", alpha=0.7)
    ax.axhline(1.0, color=C_GRAY, linestyle="-", linewidth=0.3, alpha=0.4)
    _t0_annotation(ax, 2.9)
    ax.set_xlabel("Hours from ICU admission")
    ax.set_ylabel("Cr / baseline ratio")
    ax.set_ylim(0.7, 3.0)
    ax.set_title(f"{db}: Cr ratio percentiles", fontsize=7, fontweight="bold")
    ax.legend(loc="upper left", fontsize=5, ncol=2)
    path = os.path.join(FIGS, f"fig_cr_ratio_{db.lower()}.pdf")
    fig.savefig(path)
    fig.savefig(path.replace(".pdf", ".png"))
    plt.close(fig)
    print(f"  Saved: {path}")


def plot_cumulative_aki(cohort, db):
    """Panel c/f: Cumulative KDIGO AKI using precomputed aki_kdigo1 + time_to_aki_hours."""
    trt_col = (
        "mg_supplementation" if "mg_supplementation" in cohort.columns else "mg_supp"
    )

    # Use precomputed KDIGO (includes both ratio AND delta criteria)
    has_timing = "time_to_aki_hours" in cohort.columns
    if not has_timing:
        # Fallback: aki_time_offset in minutes
        if "aki_time_offset" in cohort.columns:
            cohort = cohort.copy()
            cohort["time_to_aki_hours"] = cohort.aki_time_offset / 60.0
            has_timing = True

    fig, ax = plt.subplots(figsize=(W_SINGLE, 2.6))
    time_grid = np.arange(T0_HOURS, 169, 1)

    for g, color, glabel in [(1, C_VERM, "Mg supp"), (0, C_BLUE, "No supp")]:
        sub = cohort[cohort[trt_col] == g]
        n_total = len(sub)
        n_aki = int(sub.aki_kdigo1.sum())

        if has_timing:
            aki_pts = sub[sub.aki_kdigo1 == 1].copy()
            # For patients with AKI but no timing, assume they hit criteria
            # by the end of the observation window
            aki_times = aki_pts["time_to_aki_hours"].fillna(168).values
            cum = [np.sum(aki_times <= t) / n_total * 100 for t in time_grid]
        else:
            # No timing info: flat line at final AKI rate (crude fallback)
            final_pct = n_aki / n_total * 100
            cum = [0 if t <= T0_HOURS else final_pct for t in time_grid]

        final_pct = n_aki / n_total * 100
        ax.plot(
            time_grid,
            cum,
            color=color,
            linewidth=1.2,
            label=f"{glabel} (n={n_total}, AKI={final_pct:.1f}%)",
        )

    _t0_annotation(ax, ax.get_ylim()[1] * 0.95 if ax.get_ylim()[1] > 0 else 1)
    ax.set_xlabel("Hours from ICU admission")
    ax.set_ylabel("Cumulative AKI (%)")
    ax.set_title(f"{db}: Cumulative KDIGO Stage ≥1 AKI", fontsize=7, fontweight="bold")
    ax.legend(loc="upper left", fontsize=5)
    ax.set_xlim(0, 168)
    path = os.path.join(FIGS, f"fig_cumulative_aki_{db.lower()}.pdf")
    fig.savefig(path)
    fig.savefig(path.replace(".pdf", ".png"))
    plt.close(fig)
    print(f"  Saved: {path}")


# ── Main ────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("Cr kinetics (v3): separate panels, precomputed AKI timing\n")

    cr_e, coh_e, _ = load_eicu()
    plot_trajectory(cr_e, coh_e, "eICU")
    plot_ratio(cr_e, "eICU")
    plot_cumulative_aki(coh_e, "eICU")

    cr_m, coh_m, _ = load_mimic()
    plot_trajectory(cr_m, coh_m, "MIMIC")
    plot_ratio(cr_m, "MIMIC")
    plot_cumulative_aki(coh_m, "MIMIC")

    print("\nDone. 6 panels saved individually.")
