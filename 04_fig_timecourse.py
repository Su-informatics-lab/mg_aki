#!/usr/bin/env python3
"""
04_fig_timecourse.py — Time course of IV Mg effect on ΔCr

Reads did_riskset_{db}.csv (02_psm.R output with spec/pool/method columns)
4-panel figure: 2 databases × 2 covariate specs, PSM_DR method, both pools.

Usage:
  python 04_fig_timecourse.py
"""

import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ── Nature Portfolio rcParams ─────────────────────────────────────
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
PRIMARY_H = 24
METHOD = "psm_dr"

DB_LABELS = {"mimic": "MIMIC-IV", "eicu": "eICU-CRD"}
SPEC_LABELS = {
    "primary": "Primary (21 var, all labs)",
    "sensitivity": "Sensitivity (19 var, no K⁺/Mg)",
}
POOL_STYLES = {
    "yet_untreated": dict(
        color=BLUE,
        marker="o",
        ls="-",
        lw=1.2,
        ms=4.5,
        label="Yet-untreated",
        fill_alpha=0.12,
        zorder=4,
    ),
    "never_treated": dict(
        color=VERMIL,
        marker="s",
        ls="--",
        lw=1.0,
        ms=4,
        label="Never-treated",
        fill_alpha=0.08,
        zorder=3,
    ),
}


def load_results(tag):
    path = os.path.join(RESULTS, f"did_riskset_{tag}.csv")
    if not os.path.exists(path):
        return None
    return pd.read_csv(path)


def plot_panel(ax, df, db_tag, spec_name, panel_letter):
    """One panel: one database × one spec, two pool lines."""
    ax.axhline(0, color="#888888", lw=0.4, zorder=1)
    ax.axvline(PRIMARY_H, color="#dddddd", lw=6, alpha=0.5, zorder=0)

    for pool_name, sty in POOL_STYLES.items():
        sub = df[
            (df["spec"] == spec_name)
            & (df["pool"] == pool_name)
            & (df["method"] == METHOD)
        ].sort_values("target_h")
        sub = sub.dropna(subset=["did"])
        if len(sub) == 0:
            continue

        h = sub["target_h"].values
        did = sub["did"].values
        ci_lo = sub["ci_lo"].values
        ci_hi = sub["ci_hi"].values
        pvals = sub["p"].values

        ax.fill_between(
            h, ci_lo, ci_hi, alpha=sty["fill_alpha"], color=sty["color"], zorder=2
        )
        ax.plot(
            h,
            did,
            color=sty["color"],
            marker=sty["marker"],
            ls=sty["ls"],
            lw=sty["lw"],
            ms=sty["ms"],
            markeredgecolor=sty["color"],
            markerfacecolor="white",
            markeredgewidth=0.8,
            zorder=sty["zorder"],
            label=sty["label"],
        )

        # Fill significant markers
        for j in range(len(h)):
            if not np.isnan(pvals[j]) and pvals[j] < 0.05:
                ax.plot(
                    h[j],
                    did[j],
                    marker=sty["marker"],
                    ms=sty["ms"],
                    markeredgecolor=sty["color"],
                    markerfacecolor=sty["color"],
                    zorder=sty["zorder"] + 1,
                )

    ax.set_xticks([6, 12, 18, 24, 30, 36, 42, 48])
    ax.set_xlim(3, 51)
    ax.text(
        -0.14,
        1.06,
        panel_letter,
        transform=ax.transAxes,
        fontsize=9,
        fontweight="bold",
        va="top",
    )
    ax.set_title(
        f"{DB_LABELS.get(db_tag, db_tag)}  —  {SPEC_LABELS.get(spec_name, spec_name)}",
        fontsize=6.5,
        pad=3,
    )


def main():
    print("=" * 70)
    print("04_fig_timecourse.py — Canonical time course figure")
    print("=" * 70)

    # First clean up stale files from old runs
    for f in os.listdir(RESULTS):
        if (
            f.startswith("did_balance_")
            or f.startswith("did_pairs_yt_")
            or f.startswith("did_pairs_nt_")
        ):
            path = os.path.join(RESULTS, f)
            os.remove(path)
            print(f"  Cleaned stale: {f}")

    dbs = []
    for tag in ["mimic", "eicu"]:
        df = load_results(tag)
        if df is not None and "spec" in df.columns:
            dbs.append((tag, df))
        else:
            print(f"  {tag}: no valid results, skipping")

    if not dbs:
        print("  ERROR: no results found")
        sys.exit(1)

    specs = ["primary", "sensitivity"]
    n_rows = len(dbs)
    n_cols = len(specs)

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(7.2, 2.4 * n_rows),
        sharey="row",
        sharex=True,
        constrained_layout=True,
    )
    if n_rows == 1:
        axes = axes.reshape(1, -1)

    letters = "abcdefgh"
    idx = 0
    for row, (tag, df) in enumerate(dbs):
        for col, spec in enumerate(specs):
            ax = axes[row, col]
            plot_panel(ax, df, tag, spec, letters[idx])
            idx += 1

            if col == 0:
                ax.set_ylabel("DiD: ΔCr (mg/dL)")
            if row == n_rows - 1:
                ax.set_xlabel("Hours from T₀")
            if row == 0 and col == 0:
                ax.legend(loc="lower left", handlelength=2, borderpad=0.2)

    fig.text(
        0.5,
        -0.01,
        "Filled = P<0.05  |  Bands = 95% CI (HC1)  |  "
        "Gray band = 24h primary endpoint  |  DiD<0 = renoprotective",
        ha="center",
        fontsize=5,
        color="#666666",
    )

    out = os.path.join(RESULTS, "fig_timecourse.pdf")
    fig.savefig(out, format="pdf", dpi=300)
    fig.savefig(out.replace(".pdf", ".png"), format="png", dpi=300)
    plt.close()
    print(f"\n  ✓ {out} (.pdf + .png)")
    print("=" * 70)


if __name__ == "__main__":
    main()
