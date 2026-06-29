#!/usr/bin/env python3
"""
04c_fig_km_egfr.py — KM-style cumulative AKI by eGFR stratum

The key figure: shows the REVERSAL from protection to harm as eGFR decreases.

Panel a: eGFR ≥ 90  → solid (treated) BELOW dashed (control) = protection
Panel b: eGFR 60–89 → lines overlap = neutral
Panel c: eGFR < 45  → solid ABOVE dashed = HARM

Reads:
  results/did_all_{db}.csv
  results/did_pairs_primary_yet_untreated_{db}.csv
  results/did_cr_all_{db}.csv

Outputs:
  results/fig_km_egfr.{pdf,png}
  results/fig_km_egfr_mimic.{pdf,png}  (MIMIC-only, cleaner for slides)

Usage: python 04c_fig_km_egfr.py
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

# ── Nature / JAMA style (consistent with 04_figures.py, 04b_fig_egfr.py) ──
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
        "legend.fontsize": 6,
        "legend.frameon": False,
        "axes.linewidth": 0.5,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    }
)

WONG = {
    "blue": "#0072B2",
    "vermil": "#D55E00",
    "orange": "#E69F00",
    "skyblue": "#56B4E9",
    "green": "#009E73",
    "purple": "#CC79A7",
}

RESULTS = os.path.expanduser("~/mg_aki/results")
FIG_DIR = os.path.join(RESULTS, "figures")
os.makedirs(FIG_DIR, exist_ok=True)
TIME_GRID = np.arange(0, 169, 2)  # every 2h from 0 to 168h (7d)


# ══════════════════════════════════════════════════════════════════
# DATA LOADING + AKI TIME COMPUTATION
# ══════════════════════════════════════════════════════════════════
def load_db(tag):
    """Load matched pairs and compute time-to-first-AKI for each patient."""
    all_pts = pd.read_csv(os.path.join(RESULTS, f"did_all_{tag}.csv"))
    pairs = pd.read_csv(
        os.path.join(RESULTS, f"did_pairs_primary_yet_untreated_{tag}.csv")
    )
    cr_all = pd.read_csv(os.path.join(RESULTS, f"did_cr_all_{tag}.csv"))

    pid_col = (
        "patientunitstayid" if "patientunitstayid" in cr_all.columns else "stay_id"
    )
    cr_all["pid"] = cr_all[pid_col]
    if "offset_h" not in cr_all.columns:
        cr_all["offset_h"] = cr_all["labresultoffset"] / 60.0
    cr_all = cr_all.sort_values(["pid", "offset_h"])

    # Build pid → creatinine lookup
    cr_dict = {
        pid: grp[["labresult", "offset_h"]].values for pid, grp in cr_all.groupby("pid")
    }

    # Map pairs to eGFR
    pid_to_row = {pid: i for i, pid in enumerate(all_pts["pid"])}
    egfr_vals = all_pts["egfr"].values

    records = []
    for _, pair in pairs.iterrows():
        trt_pid = pair["trt_pid"]
        ctl_pid = pair["ctl_pid"]
        t_mg = pair["t_mg"]

        trt_row = pid_to_row.get(trt_pid)
        if trt_row is None:
            continue
        egfr = egfr_vals[trt_row]

        # Assign eGFR stratum
        if pd.isna(egfr):
            continue
        if egfr >= 90:
            stratum = "eGFR ≥ 90"
        elif egfr >= 60:
            stratum = "eGFR 60–89"
        elif egfr >= 45:
            stratum = "eGFR 45–59"
        else:
            stratum = "eGFR < 45"

        for pid, arm in [(trt_pid, "Treated"), (ctl_pid, "Control")]:
            cr = cr_dict.get(pid)
            if cr is None or len(cr) == 0:
                continue

            # Baseline: last Cr before t_mg
            pre_mask = (cr[:, 1] >= 0) & (cr[:, 1] < t_mg)
            if not np.any(pre_mask):
                continue
            pre = cr[pre_mask]
            bl = pre[np.argmax(pre[:, 1]), 0]
            if bl <= 0 or np.isnan(bl):
                continue

            # Post-T0 creatinine
            post_mask = (cr[:, 1] >= t_mg) & (cr[:, 1] <= t_mg + 168)
            post = cr[post_mask]
            if len(post) == 0:
                # Censored at T0 with no follow-up
                records.append(
                    {
                        "stratum": stratum,
                        "arm": arm,
                        "aki_time": np.nan,
                        "censor_time": 0,
                        "event": 0,
                    }
                )
                continue

            # Find time-to-first-AKI
            aki_time = np.nan
            for j in range(len(post)):
                val = post[j, 0]
                h = post[j, 1] - t_mg  # hours from T0
                delta = val - bl
                ratio = val / bl if bl > 0 else 0

                # KDIGO: absolute within 48h, ratio within 7d
                if h <= 48 and (delta >= 0.3 or ratio >= 1.5):
                    aki_time = h
                    break
                if h > 48 and ratio >= 1.5:
                    aki_time = h
                    break

            last_cr_time = post[-1, 1] - t_mg
            records.append(
                {
                    "stratum": stratum,
                    "arm": arm,
                    "aki_time": aki_time,
                    "censor_time": min(last_cr_time, 168),
                    "event": 0 if np.isnan(aki_time) else 1,
                }
            )

    df = pd.DataFrame(records)
    print(f"  {tag}: {len(df)} patient-records across {df.stratum.nunique()} strata")
    return df


def compute_cumulative_incidence(df, time_grid):
    """Compute 1 - KM estimator (cumulative incidence) on a time grid."""
    events = df[df.event == 1]["aki_time"].values
    censors = df[df.event == 0]["censor_time"].values
    n_total = len(df)
    if n_total == 0:
        return np.zeros_like(time_grid, dtype=float)

    # Sort all event/censor times
    times = np.concatenate([events, censors])
    event_flags = np.concatenate([np.ones(len(events)), np.zeros(len(censors))])
    order = np.argsort(times)
    times = times[order]
    event_flags = event_flags[order]

    # KM estimator
    n_at_risk = n_total
    surv = 1.0
    surv_curve = np.ones_like(time_grid, dtype=float)

    t_idx = 0
    for i in range(len(times)):
        # Advance time grid
        while t_idx < len(time_grid) and time_grid[t_idx] <= times[i]:
            surv_curve[t_idx] = surv
            t_idx += 1
        if event_flags[i] == 1:
            surv *= 1 - 1 / max(n_at_risk, 1)
        n_at_risk -= 1

    # Fill remaining
    while t_idx < len(time_grid):
        surv_curve[t_idx] = surv
        t_idx += 1

    return 1 - surv_curve  # cumulative incidence


def compute_ci_bootstrap(df, time_grid, n_boot=200, seed=2026):
    """Bootstrap 95% CI for cumulative incidence."""
    rng = np.random.RandomState(seed)
    n = len(df)
    boot_curves = np.zeros((n_boot, len(time_grid)))
    for b in range(n_boot):
        idx = rng.choice(n, size=n, replace=True)
        boot_curves[b] = compute_cumulative_incidence(df.iloc[idx], time_grid)
    lo = np.percentile(boot_curves, 2.5, axis=0)
    hi = np.percentile(boot_curves, 97.5, axis=0)
    return lo, hi


# ══════════════════════════════════════════════════════════════════
# PLOTTING
# ══════════════════════════════════════════════════════════════════
def plot_km_panel(
    ax, data, stratum, db_color, db_label, panel_label, show_ci=True, n_boot=200
):
    """Plot one panel: treated vs control cumulative AKI."""

    sub = data[data.stratum == stratum]
    trt = sub[sub.arm == "Treated"]
    ctl = sub[sub.arm == "Control"]

    ci_trt = compute_cumulative_incidence(trt, TIME_GRID) * 100
    ci_ctl = compute_cumulative_incidence(ctl, TIME_GRID) * 100

    # Plot control (dashed) first so treated is on top
    ax.step(
        TIME_GRID,
        ci_ctl,
        where="post",
        color=db_color,
        ls="--",
        lw=0.9,
        alpha=0.8,
        label=f"Control",
    )
    ax.step(
        TIME_GRID, ci_trt, where="post", color=db_color, ls="-", lw=1.2, label=f"IV Mg"
    )

    # CI bands
    if show_ci and len(trt) > 50:
        lo_t, hi_t = compute_ci_bootstrap(trt, TIME_GRID, n_boot)
        lo_c, hi_c = compute_ci_bootstrap(ctl, TIME_GRID, n_boot)
        ax.fill_between(
            TIME_GRID, lo_t * 100, hi_t * 100, step="post", alpha=0.10, color=db_color
        )
        ax.fill_between(
            TIME_GRID, lo_c * 100, hi_c * 100, step="post", alpha=0.06, color=db_color
        )

    # 48h marker
    ax.axvspan(46, 50, color="grey", alpha=0.08, zorder=0)
    ax.text(
        48,
        ax.get_ylim()[1] * 0.02,
        "48h",
        fontsize=5,
        ha="center",
        color="grey",
        alpha=0.6,
    )

    # Annotations
    n_trt = len(trt)
    n_ctl = len(ctl)
    rate_trt = ci_trt[-1]
    rate_ctl = ci_ctl[-1]
    ax.text(
        0.97,
        0.97,
        f"IV Mg: {rate_trt:.1f}%\nControl: {rate_ctl:.1f}%\nn={n_trt}",
        transform=ax.transAxes,
        fontsize=5.5,
        va="top",
        ha="right",
        bbox=dict(
            boxstyle="round,pad=0.3",
            facecolor="white",
            edgecolor="grey",
            alpha=0.8,
            linewidth=0.3,
        ),
    )

    # Direction arrow
    if rate_trt < rate_ctl - 1:
        ax.annotate(
            "",
            xy=(140, rate_trt),
            xytext=(140, rate_ctl),
            arrowprops=dict(
                arrowstyle="->", color=WONG["green"], lw=1.0, shrinkA=2, shrinkB=2
            ),
        )
        ax.text(
            145,
            (rate_trt + rate_ctl) / 2,
            "Protection",
            fontsize=5,
            color=WONG["green"],
            va="center",
            rotation=90,
        )
    elif rate_trt > rate_ctl + 1:
        ax.annotate(
            "",
            xy=(140, rate_trt),
            xytext=(140, rate_ctl),
            arrowprops=dict(
                arrowstyle="->", color=WONG["vermil"], lw=1.0, shrinkA=2, shrinkB=2
            ),
        )
        ax.text(
            145,
            (rate_trt + rate_ctl) / 2,
            "Harm",
            fontsize=5,
            color=WONG["vermil"],
            va="center",
            rotation=90,
        )

    ax.set_xlabel("Hours from T₀")
    ax.set_xticks([0, 24, 48, 72, 96, 120, 144, 168])
    ax.set_xticklabels(["0", "24", "48", "72", "96", "120", "144", "7d"])
    ax.set_xlim(-2, 172)

    ax.set_title(stratum, fontweight="bold", pad=6)
    ax.text(
        -0.14,
        1.05,
        panel_label,
        transform=ax.transAxes,
        fontsize=9,
        fontweight="bold",
        va="top",
    )


def main():
    print("── 04c_fig_km_egfr.py: KM-style AKI by eGFR ──\n")

    # ── Load both databases ───────────────────────────────────────
    dfs = {}
    for tag in ["mimic", "eicu"]:
        try:
            dfs[tag] = load_db(tag)
        except Exception as e:
            print(f"  {tag}: {e}")

    # ══════════════════════════════════════════════════════════════
    # FIGURE 1: MIMIC-only (clean, for slides)
    # ══════════════════════════════════════════════════════════════
    if "mimic" in dfs:
        data = dfs["mimic"]
        strata = ["eGFR ≥ 90", "eGFR 60–89", "eGFR < 45"]

        fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.8), sharey=True)

        for i, (stratum, label) in enumerate(zip(strata, ["a", "b", "c"])):
            plot_km_panel(
                axes[i],
                data,
                stratum,
                WONG["blue"],
                "MIMIC-IV",
                label,
                show_ci=True,
                n_boot=200,
            )

        axes[0].set_ylabel("Cumulative AKI incidence (%)")

        # Shared legend
        legend_elements = [
            Line2D([0], [0], color=WONG["blue"], ls="-", lw=1.2, label="IV Mg"),
            Line2D(
                [0], [0], color=WONG["blue"], ls="--", lw=0.9, label="Matched control"
            ),
        ]
        fig.legend(
            handles=legend_elements,
            loc="lower center",
            ncol=2,
            bbox_to_anchor=(0.5, -0.06),
            fontsize=6.5,
        )

        fig.suptitle(
            "Cumulative AKI Incidence by Baseline Renal Function (MIMIC-IV)",
            fontsize=8,
            fontweight="bold",
            y=1.04,
        )

        plt.tight_layout()
        for ext in ["pdf", "png"]:
            out = os.path.join(FIG_DIR, f"fig_km_egfr_mimic.{ext}")
            fig.savefig(out, format=ext, dpi=300 if ext == "png" else None)
            print(f"  Saved: {out}")
        plt.close(fig)

    # ══════════════════════════════════════════════════════════════
    # FIGURE 2: Both databases (2 rows × 3 columns)
    # ══════════════════════════════════════════════════════════════
    if len(dfs) == 2:
        strata = ["eGFR ≥ 90", "eGFR 60–89", "eGFR < 45"]
        fig, axes = plt.subplots(2, 3, figsize=(7.2, 5.0), sharey="row")

        for i, stratum in enumerate(strata):
            plot_km_panel(
                axes[0, i],
                dfs["mimic"],
                stratum,
                WONG["blue"],
                "MIMIC-IV",
                chr(ord("a") + i),
                show_ci=True,
                n_boot=200,
            )
            plot_km_panel(
                axes[1, i],
                dfs["eicu"],
                stratum,
                WONG["vermil"],
                "eICU-CRD",
                chr(ord("d") + i),
                show_ci=True,
                n_boot=200,
            )

        axes[0, 0].set_ylabel("Cumulative AKI (%)\nMIMIC-IV")
        axes[1, 0].set_ylabel("Cumulative AKI (%)\neICU-CRD")

        legend_elements = [
            Line2D([0], [0], color="grey", ls="-", lw=1.2, label="IV Mg (treated)"),
            Line2D([0], [0], color="grey", ls="--", lw=0.9, label="Matched control"),
        ]
        fig.legend(
            handles=legend_elements,
            loc="lower center",
            ncol=2,
            bbox_to_anchor=(0.5, -0.04),
            fontsize=6.5,
        )

        fig.suptitle(
            "Cumulative AKI Incidence by Baseline Renal Function",
            fontsize=9,
            fontweight="bold",
            y=1.02,
        )

        plt.tight_layout()
        for ext in ["pdf", "png"]:
            out = os.path.join(FIG_DIR, f"fig_km_egfr.{ext}")
            fig.savefig(out, format=ext, dpi=300 if ext == "png" else None)
            print(f"  Saved: {out}")
        plt.close(fig)

    print("\nDone.")


if __name__ == "__main__":
    main()
