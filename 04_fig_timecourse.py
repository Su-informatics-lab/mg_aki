#!/usr/bin/env python3
"""
04_fig_timecourse.py — 2x2 time course comparison figure

Rows: MIMIC, eICU  |  Cols: primary (no K/Mg), sens_a (all labs)
Each panel: yet-untreated (solid) vs never-treated (dashed), PSM_DR

Usage: python 04_fig_timecourse.py
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

for k, v in {
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 8,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 5.5,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "legend.frameon": False,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.03,
}.items():
    mpl.rcParams[k] = v

BLUE = "#0072B2"
VERMIL = "#D55E00"
RESULTS = os.path.expanduser("~/mg_aki/results")
METHOD = "psm_dr"
DB_LABELS = {"mimic": "MIMIC-IV", "eicu": "eICU-CRD"}
SPEC_LABELS = {
    "primary": "Primary (19 var, no K⁺/Mg)",
    "sens_a": "Sensitivity A (21 var, + K⁺/Mg)",
    "sens_b": "Sensitivity B (19 var, FIRST labs)",
}
POOL_STYLES = {
    "yet_untreated": dict(
        color=BLUE,
        marker="o",
        ls="-",
        lw=1.2,
        ms=4.5,
        label="Yet-untreated",
        fa=0.12,
        z=4,
    ),
    "never_treated": dict(
        color=VERMIL,
        marker="s",
        ls="--",
        lw=1.0,
        ms=4,
        label="Never-treated",
        fa=0.08,
        z=3,
    ),
}


def load(tag):
    p = os.path.join(RESULTS, f"did_riskset_{tag}.csv")
    if not os.path.exists(p):
        return None
    df = pd.read_csv(p)
    return df if "spec" in df.columns else None


def plot_panel(ax, df, spec, letter):
    ax.axhline(0, color="#888888", lw=0.4, zorder=1)
    ax.axvline(24, color="#dddddd", lw=6, alpha=0.5, zorder=0)
    for pn, sty in POOL_STYLES.items():
        sub = (
            df[(df.spec == spec) & (df.pool == pn) & (df.method == METHOD)]
            .sort_values("target_h")
            .dropna(subset=["did"])
        )
        if len(sub) == 0:
            continue
        h, d, lo, hi, pv = [
            sub[c].values for c in ["target_h", "did", "ci_lo", "ci_hi", "p"]
        ]
        ax.fill_between(h, lo, hi, alpha=sty["fa"], color=sty["color"], zorder=2)
        ax.plot(
            h,
            d,
            color=sty["color"],
            marker=sty["marker"],
            ls=sty["ls"],
            lw=sty["lw"],
            ms=sty["ms"],
            markeredgecolor=sty["color"],
            markerfacecolor="white",
            markeredgewidth=0.8,
            zorder=sty["z"],
            label=sty["label"],
        )
        for j in range(len(h)):
            if not np.isnan(pv[j]) and pv[j] < 0.05:
                ax.plot(
                    h[j],
                    d[j],
                    marker=sty["marker"],
                    ms=sty["ms"],
                    markeredgecolor=sty["color"],
                    markerfacecolor=sty["color"],
                    zorder=sty["z"] + 1,
                )
    ax.set_xticks([6, 12, 18, 24, 30, 36, 42, 48])
    ax.set_xlim(3, 51)
    ax.text(
        -0.14,
        1.06,
        letter,
        transform=ax.transAxes,
        fontsize=9,
        fontweight="bold",
        va="top",
    )


def main():
    print("=" * 70)
    print("04_fig_timecourse.py — 2×2 comparison figure")
    print("=" * 70)
    dbs = []
    for tag in ["mimic", "eicu"]:
        df = load(tag)
        if df is not None:
            dbs.append((tag, df))
            print(f"  {tag}: loaded")

    specs = ["primary", "sens_a"]
    nr, nc = len(dbs), len(specs)
    fig, axes = plt.subplots(
        nr,
        nc,
        figsize=(7.2, 2.4 * nr),
        sharey="row",
        sharex=True,
        constrained_layout=True,
    )
    if nr == 1:
        axes = axes.reshape(1, -1)

    letters = "abcdefgh"
    idx = 0
    for row, (tag, df) in enumerate(dbs):
        for col, spec in enumerate(specs):
            ax = axes[row, col]
            plot_panel(ax, df, spec, letters[idx])
            idx += 1
            ax.set_title(
                f"{DB_LABELS.get(tag,tag)}  —  {SPEC_LABELS.get(spec,spec)}",
                fontsize=6.5,
                pad=3,
            )
            if col == 0:
                ax.set_ylabel("DiD: ΔCr (mg/dL)")
            if row == nr - 1:
                ax.set_xlabel("Hours from T₀")
            if row == 0 and col == 0:
                ax.legend(loc="lower left", handlelength=2, borderpad=0.2)

    fig.text(
        0.5,
        -0.01,
        "Filled = P<0.05  |  Bands = 95% CI (HC1)  |  Gray = 24h primary  |  DiD<0 = renoprotective",
        ha="center",
        fontsize=5,
        color="#666666",
    )
    out = os.path.join(RESULTS, "fig_timecourse.pdf")
    fig.savefig(out, format="pdf")
    fig.savefig(out.replace(".pdf", ".png"), format="png")
    plt.close()
    print(f"\n  ✓ {out}")
    print("=" * 70)


if __name__ == "__main__":
    main()
