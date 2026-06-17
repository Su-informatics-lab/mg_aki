#!/usr/bin/env python3
"""
11_cr_kinetics.py — Post-operative creatinine kinetics visualization

  Panel a: Cr trajectory (median + IQR) by treatment group over time
  Panel b: Cr/baseline ratio over time (shows when KDIGO 1.5× threshold is crossed)
  Panel c: Cumulative proportion meeting AKI criteria over time

Reads: results/01_analysis_a_cohort.csv + raw lab tables
Output: figs/fig_cr_kinetics.pdf

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
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.linewidth": 0.5,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "legend.frameon": False,
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
C_GREEN = "#009E73"
RESULTS = os.path.expanduser("~/mg_aki/results")
FIGS = os.path.expanduser("~/mg_aki/figs")


def gz(p):
    return p if os.path.exists(p) else p.replace(".csv.gz", ".csv")


def run_eicu():
    print("eICU: Loading Cr kinetics data...")
    from importlib.util import module_from_spec, spec_from_file_location

    _cfg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "00_config.py")
    _spec = spec_from_file_location("config", _cfg_path)
    cfg = module_from_spec(_spec)
    _spec.loader.exec_module(cfg)

    cohort = pd.read_csv(os.path.join(RESULTS, "01_analysis_a_cohort.csv"))
    pids = set(cohort.patientunitstayid)

    # Map pid → baseline_cr and treatment
    bl = dict(zip(cohort.patientunitstayid, cohort.baseline_cr))
    trt = dict(
        zip(
            cohort.patientunitstayid,
            (
                cohort.mg_supplementation
                if "mg_supplementation" in cohort.columns
                else cohort.mg_supp
            ),
        )
    )
    aki = dict(zip(cohort.patientunitstayid, cohort.aki_kdigo1))

    # Load all Cr labs for cohort
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
    cr["aki_flag"] = cr.patientunitstayid.map(aki)

    # Filter: 0 to 168h (7 days)
    cr = cr[(cr.hours >= -6) & (cr.hours <= 168)].copy()
    print(f"  {len(cr):,} Cr measurements, {cr.patientunitstayid.nunique()} patients")
    return cr, cohort, "eICU"


def run_mimic():
    print("MIMIC: Loading Cr kinetics data...")
    MIMIC_HOSP = os.path.expanduser("~/mg_aki/mimic-iv-3.1/hosp")

    cohort = pd.read_csv(os.path.join(RESULTS, "04_mimic_cohort.csv"))
    cohort["intime"] = pd.to_datetime(cohort["intime"])
    stays = set(cohort.stay_id)
    hadms = set(cohort.hadm_id.dropna().astype(int))

    bl = dict(zip(cohort.stay_id, cohort.baseline_cr))
    trt_col = (
        "mg_supplementation" if "mg_supplementation" in cohort.columns else "mg_supp"
    )
    trt = dict(zip(cohort.stay_id, cohort[trt_col]))
    aki_d = dict(zip(cohort.stay_id, cohort.aki_kdigo1))
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
    cr["aki_flag"] = cr.stay_id.map(aki_d)
    cr = cr.rename(columns={"valuenum": "labresult"})
    cr = cr[(cr.hours >= -6) & (cr.hours <= 168)].copy()
    print(f"  {len(cr):,} Cr measurements, {cr.stay_id.nunique()} patients")
    return cr, cohort, "MIMIC"


def plot_kinetics(cr, cohort, db_name, axes_row):
    """Plot 3-panel Cr kinetics on a row of axes."""
    # Time bins: every 6h from -6 to 168
    bins = np.arange(-6, 174, 6)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    cr["time_bin"] = pd.cut(cr.hours, bins=bins, labels=bin_centers)
    cr["time_bin"] = cr.time_bin.astype(float)

    # ── Panel a: Cr trajectory by treatment group ────────────────
    ax = axes_row[0]
    for g, color, label in [(1, C_VERM, "Mg supp"), (0, C_BLUE, "No supp")]:
        sub = cr[cr.trt == g]
        stats = (
            sub.groupby("time_bin")["labresult"]
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
        stats.columns = ["t", "median", "n", "q25", "q75"]
        stats = stats[stats.n >= 20]  # require ≥20 per bin
        ax.fill_between(stats.t, stats.q25, stats.q75, alpha=0.15, color=color)
        ax.plot(stats.t, stats["median"], color=color, linewidth=1.2, label=label)

    ax.axvline(6, color=C_GRAY, linestyle="--", linewidth=0.5, alpha=0.7)
    ax.text(7, ax.get_ylim()[1] * 0.95, "6h\nlandmark", fontsize=5, color=C_GRAY)
    ax.set_xlabel("Hours from ICU admission")
    ax.set_ylabel("Serum creatinine (mg/dL)")
    ax.set_title(f"{db_name}: Cr trajectory", fontsize=7, fontweight="bold")
    ax.legend(fontsize=5, loc="upper right")

    # ── Panel b: Cr/baseline ratio by treatment ──────────────────
    ax = axes_row[1]
    for g, color, label in [(1, C_VERM, "Mg supp"), (0, C_BLUE, "No supp")]:
        sub = cr[(cr.trt == g) & cr.cr_ratio.notna() & (cr.cr_ratio < 5)]
        stats = (
            sub.groupby("time_bin")["cr_ratio"]
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
        stats.columns = ["t", "median", "n", "q25", "q75"]
        stats = stats[stats.n >= 20]
        ax.fill_between(stats.t, stats.q25, stats.q75, alpha=0.15, color=color)
        ax.plot(stats.t, stats["median"], color=color, linewidth=1.2, label=label)

    ax.axhline(1.5, color="red", linestyle=":", linewidth=0.7, alpha=0.6)
    ax.text(150, 1.52, "KDIGO 1.5×", fontsize=5, color="red", alpha=0.7)
    ax.axhline(1.0, color=C_GRAY, linestyle="-", linewidth=0.3, alpha=0.5)
    ax.axvline(6, color=C_GRAY, linestyle="--", linewidth=0.5, alpha=0.7)
    ax.set_xlabel("Hours from ICU admission")
    ax.set_ylabel("Cr / baseline ratio")
    ax.set_ylim(0.4, 2.5)
    ax.set_title(f"{db_name}: Cr/baseline ratio", fontsize=7, fontweight="bold")
    ax.legend(fontsize=5, loc="upper right")

    # ── Panel c: Cumulative AKI proportion over time ─────────────
    ax = axes_row[2]
    # For each time point, what fraction of patients have met KDIGO by then
    post_cr = cr[cr.hours > 6].copy()  # only after landmark

    for g, color, label in [(1, C_VERM, "Mg supp"), (0, C_BLUE, "No supp")]:
        sub = post_cr[post_cr.trt == g]
        patient_ids = set(
            sub.patientunitstayid if "patientunitstayid" in sub.columns else sub.stay_id
        )
        n_total = len(patient_ids)
        if n_total == 0:
            continue

        cum_times = []
        cum_props = []
        aki_pids = set()
        for t in sorted(bins[bins > 6]):
            window = sub[sub.hours <= t]
            # Patients meeting KDIGO by time t
            hits = window[window.cr_ratio >= 1.5]
            id_col = (
                "patientunitstayid"
                if "patientunitstayid" in hits.columns
                else "stay_id"
            )
            aki_pids |= set(hits[id_col])
            cum_times.append(t)
            cum_props.append(len(aki_pids) / n_total)

        ax.plot(
            cum_times,
            [p * 100 for p in cum_props],
            color=color,
            linewidth=1.2,
            label=f"{label} (n={n_total})",
        )

    ax.axvline(6, color=C_GRAY, linestyle="--", linewidth=0.5, alpha=0.7)
    ax.set_xlabel("Hours from ICU admission")
    ax.set_ylabel("Cumulative AKI (%)")
    ax.set_title(f"{db_name}: Cumulative AKI incidence", fontsize=7, fontweight="bold")
    ax.legend(fontsize=5, loc="upper left")


if __name__ == "__main__":
    fig, axes = plt.subplots(2, 3, figsize=(7.205, 5.0))
    fig.subplots_adjust(hspace=0.45, wspace=0.35)

    cr_e, coh_e, _ = run_eicu()
    plot_kinetics(cr_e, coh_e, "eICU", axes[0])

    cr_m, coh_m, _ = run_mimic()
    plot_kinetics(cr_m, coh_m, "MIMIC", axes[1])

    for i, ax in enumerate(axes.flat):
        ax.text(
            -0.12,
            1.06,
            chr(ord("a") + i),
            transform=ax.transAxes,
            fontsize=8,
            fontweight="bold",
            va="top",
            ha="right",
        )

    path = os.path.join(FIGS, "fig_cr_kinetics.pdf")
    fig.savefig(path)
    fig.savefig(path.replace(".pdf", ".png"))
    plt.close(fig)
    print(f"\n  Saved: {path}")
    print("Done. Show this to Drs. Su and Meng.")
