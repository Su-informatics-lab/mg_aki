#!/usr/bin/env python3
"""
04b_fig_egfr.py — eGFR-stratified dose-response figure

Reads from:  results/egfr_aki_stages_{db}.csv  (output of 03b_egfr_aki_stages.R)
Outputs:     results/fig_egfr_doseresponse.{pdf,png}

Panel a: AKI Stage 1+ by eGFR stratum (both databases)
Panel b: Hospital mortality by eGFR stratum (both databases)
Panel c: AKI Stage 2+ by eGFR stratum (severe AKI)
Panel d: AKI Stage 3+ by eGFR stratum (very severe AKI)

Usage: python 04b_fig_egfr.py
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

# ── Nature Portfolio / JAMA style (consistent with 04_figures.py) ─────
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
mpl.rcParams["font.size"] = 7
mpl.rcParams["axes.labelsize"] = 7
mpl.rcParams["axes.titlesize"] = 8
mpl.rcParams["xtick.labelsize"] = 6
mpl.rcParams["ytick.labelsize"] = 6.5
mpl.rcParams["legend.fontsize"] = 6
mpl.rcParams["axes.linewidth"] = 0.5
mpl.rcParams["xtick.major.width"] = 0.5
mpl.rcParams["ytick.major.width"] = 0.5
mpl.rcParams["xtick.major.size"] = 3
mpl.rcParams["ytick.major.size"] = 3
mpl.rcParams["xtick.direction"] = "out"
mpl.rcParams["ytick.direction"] = "out"
mpl.rcParams["legend.frameon"] = False
mpl.rcParams["axes.grid"] = False
mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["figure.facecolor"] = "white"
mpl.rcParams["savefig.facecolor"] = "white"
mpl.rcParams["savefig.dpi"] = 300
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["savefig.pad_inches"] = 0.02

# Wong/Okabe-Ito colorblind-safe palette (same as 03_did_figures.py)
WONG = {
    "blue": "#0072B2",
    "vermil": "#D55E00",
    "orange": "#E69F00",
    "skyblue": "#56B4E9",
    "green": "#009E73",
    "purple": "#CC79A7",
    "black": "#000000",
}

# Database visual encoding (consistent with existing figures)
DB_STYLE = {
    "MIMIC": {
        "color": WONG["blue"],
        "marker": "o",
        "label": "MIMIC-IV",
        "offset": 0.12,
    },
    "EICU": {
        "color": WONG["vermil"],
        "marker": "s",
        "label": "eICU-CRD",
        "offset": -0.12,
    },
}

RESULTS = os.path.expanduser("~/mg_aki/results")

# eGFR stratum display order (top to bottom: healthy → impaired)
STRATA_ORDER = ["eGFR>=90", "eGFR_60-89", "eGFR_45-59", "eGFR_30-44", "eGFR<30"]
STRATA_LABELS = {
    "eGFR>=90": "eGFR ≥ 90\n(CKD Stage 1)",
    "eGFR_60-89": "eGFR 60–89\n(CKD Stage 2)",
    "eGFR_45-59": "eGFR 45–59\n(CKD Stage 3a)",
    "eGFR_30-44": "eGFR 30–44\n(CKD Stage 3b)",
    "eGFR<30": "eGFR < 30\n(CKD Stage 4–5)",
}


def load_data():
    """Load both databases, return combined DataFrame."""
    frames = []
    for db in ["mimic", "eicu"]:
        path = os.path.join(RESULTS, f"egfr_aki_stages_{db}.csv")
        if os.path.exists(path):
            df = pd.read_csv(path)
            frames.append(df)
            print(f"  Loaded {path}: {len(df)} rows")
        else:
            print(f"  WARNING: {path} not found")
    if not frames:
        raise FileNotFoundError("No data files found")
    return pd.concat(frames, ignore_index=True)


def plot_forest_panel(ax, data, outcome, title, panel_label, xlog=True):
    """Plot one forest panel: eGFR strata on y-axis, OR on x-axis."""

    strata = [s for s in STRATA_ORDER if s in data.stratum.values]
    y_positions = {s: i for i, s in enumerate(reversed(strata))}

    # Reference line at OR=1
    ax.axvline(1, color="grey", lw=0.5, ls="-", zorder=0)

    # Shading for harm zone (OR > 1)
    if xlog:
        ax.axvspan(1, 20, color="#FFEEEE", alpha=0.3, zorder=0)
        ax.axvspan(0.02, 1, color="#EEEEFF", alpha=0.3, zorder=0)

    for db_key, style in DB_STYLE.items():
        sub = data[(data.db == db_key) & (data.outcome == outcome)]
        if len(sub) == 0:
            continue

        for _, row in sub.iterrows():
            if row.stratum not in y_positions:
                continue
            if pd.isna(row["or"]):
                continue

            y = y_positions[row.stratum] + style["offset"]
            or_val = row["or"]
            lo = row["or_lo"]
            hi = row["or_hi"]
            sig = row["p"] < 0.05 if pd.notna(row["p"]) else False

            # CI whiskers
            ax.plot([lo, hi], [y, y], color=style["color"], lw=0.8, zorder=2)

            # Point: filled if significant, open if not
            if sig:
                ax.plot(
                    or_val,
                    y,
                    marker=style["marker"],
                    color=style["color"],
                    ms=5,
                    zorder=3,
                    markeredgewidth=0.5,
                    markeredgecolor="white",
                )
            else:
                ax.plot(
                    or_val,
                    y,
                    marker=style["marker"],
                    color="white",
                    ms=5,
                    zorder=3,
                    markeredgewidth=1.0,
                    markeredgecolor=style["color"],
                )

            # n annotation (right side)
            n_txt = f"n={row.n:.0f}" if pd.notna(row.n) else ""
            if db_key == "MIMIC":
                ax.text(
                    hi * 1.15 if xlog else hi + 0.1,
                    y,
                    n_txt,
                    fontsize=5,
                    color=style["color"],
                    va="center",
                    ha="left",
                )

    # Y-axis labels
    ax.set_yticks(list(y_positions.values()))
    ax.set_yticklabels([STRATA_LABELS.get(s, s) for s in reversed(strata)])

    # X-axis
    if xlog:
        ax.set_xscale("log")
        ax.set_xlim(0.15, 15)
        ax.set_xticks([0.25, 0.5, 1, 2, 4, 8])
        ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.set_xlabel("Odds Ratio (95% CI)")

    # Title and panel label
    ax.set_title(title, fontweight="bold", pad=8)
    ax.text(
        -0.22,
        1.05,
        panel_label,
        transform=ax.transAxes,
        fontsize=9,
        fontweight="bold",
        va="top",
    )

    # Direction annotations
    ax.text(
        0.18,
        -0.08,
        "← Favors IV Mg",
        transform=ax.transAxes,
        fontsize=5,
        color=WONG["blue"],
        fontstyle="italic",
        ha="left",
    )
    ax.text(
        0.82,
        -0.08,
        "Favors control →",
        transform=ax.transAxes,
        fontsize=5,
        color=WONG["vermil"],
        fontstyle="italic",
        ha="right",
    )


def main():
    print("── 04b_fig_egfr.py: eGFR dose-response figure ──")
    data = load_data()

    # ── 4-panel figure ────────────────────────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.5))

    plot_forest_panel(
        axes[0, 0], data, "AKI_Stage1+", "AKI Stage ≥1 (KDIGO, 7-day)", "a"
    )
    plot_forest_panel(axes[0, 1], data, "hosp_mortality", "Hospital Mortality", "b")
    plot_forest_panel(
        axes[1, 0], data, "AKI_Stage2+", "AKI Stage ≥2 (Cr ≥2× baseline)", "c"
    )
    plot_forest_panel(
        axes[1, 1], data, "AKI_Stage3+", "AKI Stage ≥3 (Cr ≥3× baseline)", "d"
    )

    # Legend
    legend_elements = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=WONG["blue"],
            markeredgecolor=WONG["blue"],
            ms=6,
            label="MIMIC-IV",
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            color="w",
            markerfacecolor=WONG["vermil"],
            markeredgecolor=WONG["vermil"],
            ms=6,
            label="eICU-CRD",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=WONG["blue"],
            markeredgecolor="white",
            ms=6,
            label="P < .05 (filled)",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="white",
            markeredgecolor=WONG["blue"],
            ms=6,
            label="P ≥ .05 (open)",
        ),
    ]
    fig.legend(
        handles=legend_elements,
        loc="lower center",
        ncol=4,
        bbox_to_anchor=(0.5, -0.02),
        fontsize=6.5,
    )

    fig.suptitle(
        "eGFR-Stratified Treatment Effect of IV Magnesium\n"
        "on AKI and Mortality After Cardiac Surgery",
        fontsize=9,
        fontweight="bold",
        y=1.02,
    )

    plt.tight_layout()

    # Save
    for ext in ["pdf", "png"]:
        out = os.path.join(RESULTS, f"fig_egfr_doseresponse.{ext}")
        fig.savefig(out, format=ext, dpi=300 if ext == "png" else None)
        print(f"  Saved: {out}")

    plt.close(fig)

    # ── Compact 2-panel version (for slides) ──────────────────────
    fig2, axes2 = plt.subplots(1, 2, figsize=(7.2, 3.2))

    plot_forest_panel(axes2[0], data, "AKI_Stage1+", "AKI Stage ≥1 (7-day)", "a")
    plot_forest_panel(axes2[1], data, "hosp_mortality", "Hospital Mortality", "b")

    fig2.legend(
        handles=legend_elements[:2],
        loc="lower center",
        ncol=2,
        bbox_to_anchor=(0.5, -0.06),
        fontsize=6.5,
    )

    plt.tight_layout()

    for ext in ["pdf", "png"]:
        out = os.path.join(RESULTS, f"fig_egfr_doseresponse_2panel.{ext}")
        fig2.savefig(out, format=ext, dpi=300 if ext == "png" else None)
        print(f"  Saved: {out}")

    plt.close(fig2)
    print("\nDone.")


if __name__ == "__main__":
    main()
