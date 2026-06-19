#!/usr/bin/env python3
"""
fig_did_timecourse.py — Publication figures for Mg→AKI DiD study

Figure 1: Matching-based DiD time course (PRIMARY)
  Panel a: eICU, Panel b: MIMIC
  Three tolerance bands (±2h, ±4h, ±6h)
  X = hours post-IV-Mg, Y = DiD estimate ± 95% CI

Figure S1: IPTW variant comparison (SUPPLEMENTARY)
  Panel a: eICU, Panel b: MIMIC
  Five methods on common ICU-time anchor
  X = hours from ICU, Y = DiD estimate

Usage: python fig_did_timecourse.py
  Reads from ~/mg_aki/results/
  Outputs to ~/mg_aki/results/
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ── Nature Portfolio style ──────────────────────────────────────────────
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
mpl.rcParams["font.size"] = 7
mpl.rcParams["axes.labelsize"] = 7
mpl.rcParams["axes.titlesize"] = 7
mpl.rcParams["xtick.labelsize"] = 6
mpl.rcParams["ytick.labelsize"] = 6
mpl.rcParams["legend.fontsize"] = 6
mpl.rcParams["legend.title_fontsize"] = 7
mpl.rcParams["axes.linewidth"] = 0.5
mpl.rcParams["xtick.major.width"] = 0.5
mpl.rcParams["ytick.major.width"] = 0.5
mpl.rcParams["xtick.major.size"] = 3
mpl.rcParams["ytick.major.size"] = 3
mpl.rcParams["xtick.direction"] = "out"
mpl.rcParams["ytick.direction"] = "out"
mpl.rcParams["lines.linewidth"] = 1.0
mpl.rcParams["lines.markersize"] = 4
mpl.rcParams["legend.frameon"] = False
mpl.rcParams["axes.grid"] = False
mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["figure.facecolor"] = "white"
mpl.rcParams["savefig.facecolor"] = "white"
mpl.rcParams["savefig.dpi"] = 300
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["savefig.pad_inches"] = 0.02
mpl.rcParams["figure.constrained_layout.use"] = True

# Wong/Okabe-Ito colorblind-safe palette
WONG = {
    "black": "#000000",
    "orange": "#E69F00",
    "skyblue": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermil": "#D55E00",
    "purple": "#CC79A7",
}

RESULTS = os.path.expanduser("~/mg_aki/results")

# Double-column width for 2-panel figures
FIGW = 7.205  # inches (183 mm)
FIGH = 3.0  # inches


# ============================================================================
# FIGURE 1: Matching-based DiD time course (PRIMARY)
# ============================================================================
def fig1_matching():
    print("── Figure 1: Matching-based DiD (primary) ──")

    fig, axes = plt.subplots(1, 2, figsize=(FIGW, FIGH), sharey=True)

    tol_styles = {
        2: {
            "color": WONG["skyblue"],
            "ls": ":",
            "lw": 0.8,
            "alpha": 0.7,
            "marker": "v",
            "ms": 3,
            "label": "±2 h",
        },
        4: {
            "color": WONG["orange"],
            "ls": "--",
            "lw": 0.8,
            "alpha": 0.8,
            "marker": "s",
            "ms": 3,
            "label": "±4 h",
        },
        6: {
            "color": WONG["blue"],
            "ls": "-",
            "lw": 1.2,
            "alpha": 1.0,
            "marker": "o",
            "ms": 4,
            "label": "±6 h (primary)",
        },
    }

    for idx, (db, label) in enumerate([("eicu", "eICU-CRD"), ("mimic", "MIMIC-IV")]):
        ax = axes[idx]
        path = os.path.join(RESULTS, f"did_timing_sweep_{db}.csv")
        if not os.path.exists(path):
            print(f"  {path} not found, skipping")
            continue

        df = pd.read_csv(path)
        df_match = df[df["method"] == "matching"].copy()

        for tol_h, style in tol_styles.items():
            sub = df_match[df_match["tol_h"] == tol_h].sort_values("target_h")
            if len(sub) == 0:
                continue

            x = sub["target_h"].values
            y = sub["did_adj"].values
            lo = sub["ci_lo"].values
            hi = sub["ci_hi"].values
            sig = sub["p_adj"].values < 0.05

            ax.plot(
                x,
                y,
                color=style["color"],
                ls=style["ls"],
                lw=style["lw"],
                alpha=style["alpha"],
                marker=style["marker"],
                ms=style["ms"],
                label=style["label"],
                zorder=3,
            )

            # CI shading for primary tolerance only
            if tol_h == 6:
                ax.fill_between(x, lo, hi, color=style["color"], alpha=0.12, zorder=1)

            # Mark significant points
            if np.any(sig):
                ax.scatter(
                    x[sig],
                    y[sig],
                    color=style["color"],
                    marker="*",
                    s=50,
                    zorder=5,
                    edgecolors="none",
                )

        ax.axhline(0, color="grey", lw=0.5, ls="-", zorder=0)
        ax.set_xlabel("Hours after IV MgSO₄")
        ax.set_xticks(range(6, 37, 6))
        ax.set_title(label, fontweight="bold", fontsize=8)

        # Panel label
        ax.text(
            -0.12,
            1.05,
            chr(ord("a") + idx),
            transform=ax.transAxes,
            fontsize=8,
            fontweight="bold",
            va="top",
        )

    axes[0].set_ylabel("DiD estimate, ΔΔCr (mg/dL)")
    axes[0].legend(loc="lower left", title="Temporal tolerance")

    # Annotation
    fig.text(
        0.5,
        -0.02,
        "Negative values indicate IV Mg protective (less creatinine rise)",
        ha="center",
        fontsize=6,
        fontstyle="italic",
        color="grey",
    )

    out = os.path.join(RESULTS, "fig1_did_matching_timecourse.pdf")
    fig.savefig(out, format="pdf")
    out_png = out.replace(".pdf", ".png")
    fig.savefig(out_png, format="png", dpi=300)
    print(f"  Saved: {out}")
    print(f"  Saved: {out_png}")
    plt.close(fig)


# ============================================================================
# FIGURE S1: IPTW variant comparison (SUPPLEMENTARY)
# ============================================================================
def figs1_iptw():
    print("\n── Figure S1: IPTW variants (supplementary) ──")

    fig, axes = plt.subplots(1, 2, figsize=(FIGW, FIGH), sharey=True)

    method_styles = {
        "sIPTW": {
            "color": WONG["skyblue"],
            "ls": ":",
            "lw": 0.7,
            "marker": "v",
            "ms": 2,
            "label": "sIPTW",
        },
        "sIPTW_t99": {
            "color": WONG["orange"],
            "ls": "--",
            "lw": 0.8,
            "marker": "s",
            "ms": 2,
            "label": "sIPTW (trim 99th)",
        },
        "sIPTW_t95": {
            "color": WONG["green"],
            "ls": "--",
            "lw": 0.8,
            "marker": "^",
            "ms": 2,
            "label": "sIPTW (trim 95th)",
        },
        "sIPTW_DR": {
            "color": WONG["blue"],
            "ls": "-",
            "lw": 1.2,
            "marker": "o",
            "ms": 4,
            "label": "sIPTW-DR (primary)",
        },
        "AIPW": {
            "color": WONG["vermil"],
            "ls": "-",
            "lw": 0.8,
            "marker": "D",
            "ms": 3,
            "label": "AIPW",
        },
    }

    for idx, (db, label) in enumerate([("eicu", "eICU-CRD"), ("mimic", "MIMIC-IV")]):
        ax = axes[idx]
        path = os.path.join(RESULTS, f"did_iptw_{db}.csv")
        if not os.path.exists(path):
            print(f"  {path} not found, skipping")
            continue

        df = pd.read_csv(path)

        for meth, style in method_styles.items():
            sub = df[df["method"] == meth].sort_values("target_h")
            if len(sub) == 0:
                continue

            x = sub["target_h"].values
            y = sub["did"].values
            sig = sub["p"].values < 0.05

            ax.plot(
                x,
                y,
                color=style["color"],
                ls=style["ls"],
                lw=style["lw"],
                marker=style["marker"],
                ms=style["ms"],
                label=style["label"],
                zorder=3,
            )

            # CI shading for primary method only
            if meth == "sIPTW_DR" and "ci_lo" in sub.columns:
                lo = sub["ci_lo"].values
                hi = sub["ci_hi"].values
                ax.fill_between(x, lo, hi, color=style["color"], alpha=0.12, zorder=1)

            if np.any(sig):
                ax.scatter(
                    x[sig],
                    y[sig],
                    color=style["color"],
                    marker="*",
                    s=50,
                    zorder=5,
                    edgecolors="none",
                )

        ax.axhline(0, color="grey", lw=0.5, ls="-", zorder=0)
        ax.set_xlabel("Hours after ICU admission")
        ax.set_xticks(range(6, 37, 6))
        ax.set_title(label, fontweight="bold", fontsize=8)
        ax.text(
            -0.12,
            1.05,
            chr(ord("a") + idx),
            transform=ax.transAxes,
            fontsize=8,
            fontweight="bold",
            va="top",
        )

    axes[0].set_ylabel("DiD estimate, ΔΔCr (mg/dL)")
    axes[0].legend(loc="lower left", title="IPTW variant")

    fig.text(
        0.5,
        -0.02,
        "All methods estimate ATE using ICU-time anchor; "
        "negative values = IV Mg protective",
        ha="center",
        fontsize=6,
        fontstyle="italic",
        color="grey",
    )

    out = os.path.join(RESULTS, "figs1_iptw_timecourse.pdf")
    fig.savefig(out, format="pdf")
    out_png = out.replace(".pdf", ".png")
    fig.savefig(out_png, format="png", dpi=300)
    print(f"  Saved: {out}")
    print(f"  Saved: {out_png}")
    plt.close(fig)


# ============================================================================
if __name__ == "__main__":
    fig1_matching()
    figs1_iptw()
    print("\nDone. Both figures saved to", RESULTS)
