#!/usr/bin/env python3
"""
04_fig_timecourse.py — Time course of IV Mg renoprotective effect

Reads did_riskset_{db}.csv (output of 02_psm.R with dual control pools)
and produces a Nature-style figure showing:
  - DiD (ΔCr treatment effect) from T₀+6h to T₀+36h
  - Yet-untreated (primary) and never-treated (sensitivity) side by side
  - 95% CI bands
  - Two panels: eICU and MIMIC

Usage:
  python 04_fig_timecourse.py            # both databases
  python 04_fig_timecourse.py eicu       # eICU only
  python 04_fig_timecourse.py mimic      # MIMIC only
"""

import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ═══════════════════════════════════════════════════════════════════
# Nature Portfolio rcParams
# ═══════════════════════════════════════════════════════════════════
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

# Wong/Okabe-Ito
BLUE = "#0072B2"
VERMIL = "#D55E00"
SKY = "#56B4E9"
ORANGE = "#E69F00"

RESULTS = os.path.expanduser("~/mg_aki/results")
DB_LABELS = {"eicu": "eICU-CRD", "mimic": "MIMIC-IV"}
PRIMARY_H = 24


def load_results(tag):
    path = os.path.join(RESULTS, f"did_riskset_{tag}.csv")
    if not os.path.exists(path):
        print(f"  WARN: {path} not found, skipping {tag}")
        return None
    df = pd.read_csv(path)
    # Ensure pool column exists (back-compat with old single-pool runs)
    if "pool" not in df.columns:
        df["pool"] = "yet_untreated"
    return df


def plot_timecourse(dbs, out_path):
    """Two-panel (or single-panel) time course figure."""
    n_panels = len(dbs)
    # Nature double-column = 183mm = 7.205in; single = 89mm = 3.504in
    if n_panels == 2:
        fig, axes = plt.subplots(1, 2, figsize=(7.205, 2.8), sharey=True)
    else:
        fig, ax = plt.subplots(1, 1, figsize=(3.504, 2.8))
        axes = [ax]

    for i, (tag, ax) in enumerate(zip(dbs, axes)):
        df = load_results(tag)
        if df is None:
            ax.text(
                0.5,
                0.5,
                f"No data for {tag}",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=7,
                color="gray",
            )
            continue

        # Reference lines
        ax.axhline(0, color="#888888", linewidth=0.4, linestyle="-", zorder=1)
        ax.axvline(PRIMARY_H, color="#cccccc", linewidth=0.4, linestyle=":", zorder=1)

        pool_styles = {
            "yet_untreated": dict(
                color=BLUE,
                marker="o",
                ls="-",
                lw=1.2,
                ms=5,
                zorder=4,
                label="Yet-untreated (primary)",
                fill_alpha=0.12,
            ),
            "never_treated": dict(
                color=VERMIL,
                marker="s",
                ls="--",
                lw=1.0,
                ms=4.5,
                zorder=3,
                label="Never-treated (sensitivity)",
                fill_alpha=0.08,
            ),
        }

        for pool_name, sty in pool_styles.items():
            sub = df[df["pool"] == pool_name].sort_values("target_h")
            if len(sub) == 0:
                continue
            sub = sub.dropna(subset=["did"])
            if len(sub) == 0:
                continue

            h = sub["target_h"].values
            did = sub["did"].values
            ci_lo = sub["ci_lo"].values
            ci_hi = sub["ci_hi"].values
            pvals = sub["p"].values

            # CI band
            ax.fill_between(
                h, ci_lo, ci_hi, alpha=sty["fill_alpha"], color=sty["color"], zorder=2
            )
            # Line + markers
            ax.plot(
                h,
                did,
                color=sty["color"],
                marker=sty["marker"],
                linestyle=sty["ls"],
                linewidth=sty["lw"],
                markersize=sty["ms"],
                markeredgecolor=sty["color"],
                markerfacecolor="white",
                markeredgewidth=0.8,
                zorder=sty["zorder"],
                label=sty["label"],
            )

            # Significance markers
            for j in range(len(h)):
                if not np.isnan(pvals[j]) and pvals[j] < 0.05:
                    ax.plot(
                        h[j],
                        did[j],
                        marker=sty["marker"],
                        markersize=sty["ms"],
                        markeredgecolor=sty["color"],
                        markerfacecolor=sty["color"],
                        zorder=sty["zorder"] + 1,
                    )

        # Formatting
        ax.set_xlabel("Hours from T₀ (first IV Mg)")
        if i == 0:
            ax.set_ylabel("DiD: ΔCr treatment effect (mg/dL)")
        ax.set_xticks([6, 12, 18, 24, 30, 36])
        ax.set_xlim(3, 39)

        # Panel label + title
        panel_letter = chr(ord("a") + i)
        ax.text(
            -0.12,
            1.08,
            panel_letter,
            transform=ax.transAxes,
            fontsize=8,
            fontweight="bold",
            va="top",
            ha="right",
        )
        ax.set_title(DB_LABELS.get(tag, tag), fontsize=7, pad=4)

        # Annotate primary endpoint
        ax.annotate(
            "24h primary",
            xy=(24, 0),
            xytext=(24, ax.get_ylim()[0] * 0.15),
            fontsize=5,
            color="#999999",
            ha="center",
            va="top",
            arrowprops=dict(arrowstyle="-", color="#cccccc", lw=0.3),
        )

        ax.legend(
            loc="lower left",
            fontsize=5.5,
            handlelength=2.5,
            borderpad=0.2,
            labelspacing=0.3,
        )

    # Shared legend note
    fig.text(
        0.5,
        -0.02,
        "Filled markers = P < 0.05. "
        "Bands = 95% CI (HC1 robust SE). "
        "DiD < 0 indicates renoprotection.",
        ha="center",
        fontsize=5,
        color="#666666",
    )

    fig.savefig(out_path.replace(".pdf", ".pdf"), format="pdf", dpi=300)
    fig.savefig(out_path.replace(".pdf", ".png"), format="png", dpi=300)
    plt.close(fig)
    print(f"  ✓ Saved: {out_path} (.pdf + .png)")


if __name__ == "__main__":
    print("=" * 70)
    print("04_fig_timecourse.py — DiD time course (dual control pool)")
    print("=" * 70)

    args = [a.lower() for a in sys.argv[1:]]

    # Detect which databases have results
    available = []
    for tag in ["eicu", "mimic"]:
        if not args or tag in args:
            p = os.path.join(RESULTS, f"did_riskset_{tag}.csv")
            if os.path.exists(p):
                available.append(tag)
            else:
                print(f"  {tag}: no results file, skipping")

    if not available:
        print("  ERROR: no result files found in", RESULTS)
        sys.exit(1)

    out = os.path.join(RESULTS, "fig_timecourse_dual.pdf")
    plot_timecourse(available, out)

    print("\n" + "=" * 70)
    print("Done.")
    print("=" * 70)
