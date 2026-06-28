#!/usr/bin/env python3
"""
04d_fig_secondary_egfr.py — eGFR-stratified multi-endpoint forest + heatmap

Shows the cross-organ eGFR gradient: eGFR≥90 = safe (AKI/mort protective,
rest null), eGFR<45 = multi-organ harm zone.

  Panel A: Forest plot — 6 endpoints × 5 eGFR strata
  Panel B: Heatmap    — OR colored blue/white/red, significance overlay

Run:  python 04d_fig_secondary_egfr.py
"""

import os
import sys

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

for k, v in {
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": "white",
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
}.items():
    mpl.rcParams[k] = v

BLUE = "#0072B2"
VERMIL = "#D55E00"
GRAY = "#999999"
RESULTS = os.path.expanduser("~/mg_aki/results")

ENDPOINTS = [
    ("aki_7d", "7-day AKI"),
    ("hosp_mortality", "Mortality"),
    ("poaf", "POAF (LLM)"),
    ("encephalopathy_delirium", "Enceph/delirium"),
    ("transfusion", "Transfusion"),
    ("reintubation", "Reintubation"),
]

EGFR_STRATA = [
    ("eGFR >= 90", "≥90 (G1)"),
    ("eGFR 60-89", "60–89 (G2)"),
    ("eGFR 45-59", "45–59 (G3a)"),
    ("eGFR 30-44", "30–44 (G3b)"),
    ("eGFR < 30", "<30 (G4–5)"),
]

# Organ-system colors
EP_COLORS = {
    "aki_7d": "#2166AC",  # renal — blue
    "hosp_mortality": "#333333",  # overall — dark
    "poaf": "#D6604D",  # cardiac — red
    "encephalopathy_delirium": "#8856A7",  # neuro — purple
    "transfusion": "#E08214",  # hemato — orange
    "reintubation": "#35978F",  # pulmonary — teal
}


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(RESULTS, f"{name}.{ext}"), format=ext)
    plt.close(fig)
    print(f"  ✓ {name}")


def load_hte():
    p = os.path.join(RESULTS, "did_hte_mimic.csv")
    if not os.path.exists(p):
        print(f"  ERROR: {p} not found. Run 03_hte.R mimic first.")
        sys.exit(1)
    return pd.read_csv(p)


def get_row(df, oc, sg):
    r = df[(df.outcome == oc) & (df.subgroup == sg)]
    if len(r) == 0:
        return None
    return r.iloc[0]


# ====================================================================
# PANEL A: FOREST PLOT
# ====================================================================
def make_forest(ax, df):
    """Multi-endpoint eGFR-stratified forest plot."""
    ax.axvline(1, color="#ddd", lw=0.6, zorder=0)

    y = 0
    yticks, ylabels = [], []
    n_ep = len(ENDPOINTS)

    for ei, (oc, oc_lbl) in enumerate(reversed(ENDPOINTS)):
        color = EP_COLORS[oc]

        for si, (sg, sg_lbl) in enumerate(reversed(EGFR_STRATA)):
            r = get_row(df, oc, sg)
            if r is None or pd.isna(r["or"]) or r["or"] > 20:
                yticks.append(y)
                ylabels.append(sg_lbl if ei == n_ep - 1 else "")
                y += 1
                continue

            est = r["or"]
            lo, hi = r["or_lo"], r["or_hi"]
            sig = not pd.isna(r["p"]) and r["p"] < 0.05

            ax.errorbar(
                est,
                y,
                xerr=[[max(est - lo, 0.001)], [max(hi - est, 0.001)]],
                fmt="o",
                color=color,
                ms=5 if sig else 3.5,
                markerfacecolor=color if sig else "white",
                markeredgecolor=color,
                markeredgewidth=0.7,
                capsize=1.5,
                capthick=0.4,
                lw=0.5,
                zorder=3,
            )
            # n label
            ax.text(
                0.055,
                y,
                f"n={int(r['n'])}",
                fontsize=4.5,
                color=GRAY,
                va="center",
                transform=mpl.transforms.blended_transform_factory(
                    ax.transAxes, ax.transData
                ),
            )
            yticks.append(y)
            ylabels.append(sg_lbl if ei == n_ep - 1 else "")
            y += 1

        # Endpoint label (right side)
        mid_y = y - len(EGFR_STRATA) / 2 - 0.5
        ax.text(
            15,
            mid_y,
            oc_lbl,
            fontsize=6.5,
            fontweight="bold",
            color=color,
            va="center",
            ha="left",
        )

        # Overall
        r_ov = get_row(df, oc, "Overall")
        if r_ov is not None and not pd.isna(r_ov["or"]):
            sig_ov = not pd.isna(r_ov["p"]) and r_ov["p"] < 0.05
            ax.errorbar(
                r_ov["or"],
                y,
                xerr=[
                    [max(r_ov["or"] - r_ov["or_lo"], 0.001)],
                    [max(r_ov["or_hi"] - r_ov["or"], 0.001)],
                ],
                fmt="D",
                color=color,
                ms=5,
                markerfacecolor=color if sig_ov else "white",
                markeredgecolor=color,
                markeredgewidth=0.8,
                capsize=2,
                capthick=0.5,
                lw=0.6,
                zorder=3,
            )
            yticks.append(y)
            ylabels.append("Overall" if ei == n_ep - 1 else "")
            y += 1

        # Gap between endpoints
        y += 0.8

    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=5.5)
    ax.set_xscale("log")
    ax.set_xlim(0.15, 12)
    ax.set_xticks([0.25, 0.5, 1, 2, 4, 8])
    ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.set_xlabel("Odds Ratio (95% CI)")

    ax.text(
        0.01, -0.04, "← Favors IV Mg", transform=ax.transAxes, fontsize=5, color=GRAY
    )
    ax.text(
        0.99,
        -0.04,
        "Favors control →",
        transform=ax.transAxes,
        fontsize=5,
        color=GRAY,
        ha="right",
    )


# ====================================================================
# PANEL B: HEATMAP
# ====================================================================
def make_heatmap(ax, df):
    """eGFR × Endpoint OR heatmap."""
    n_egfr = len(EGFR_STRATA)
    n_ep = len(ENDPOINTS)
    or_mat = np.full((n_egfr, n_ep), np.nan)
    p_mat = np.full((n_egfr, n_ep), np.nan)
    n_mat = np.full((n_egfr, n_ep), 0)

    for ei, (oc, _) in enumerate(ENDPOINTS):
        for si, (sg, _) in enumerate(EGFR_STRATA):
            r = get_row(df, oc, sg)
            if r is not None and not pd.isna(r["or"]):
                or_mat[si, ei] = r["or"]
                p_mat[si, ei] = r["p"]
                n_mat[si, ei] = int(r["n"])

    cmap = plt.get_cmap("RdBu_r")
    vmin, vmax = np.log(0.2), np.log(5.0)
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

    log_or = np.log(np.where(np.isnan(or_mat), 1, or_mat))
    ax.imshow(log_or, cmap=cmap, norm=norm, aspect="auto")

    for si in range(n_egfr):
        for ei in range(n_ep):
            if np.isnan(or_mat[si, ei]):
                ax.text(ei, si, "—", ha="center", va="center", fontsize=6, color=GRAY)
                continue
            orv = or_mat[si, ei]
            pv = p_mat[si, ei]
            nv = n_mat[si, ei]
            sig = "*" if (not np.isnan(pv) and pv < 0.05) else ""
            bg = norm(np.log(orv))
            tc = "white" if (bg < 0.25 or bg > 0.75) else "black"
            ax.text(
                ei,
                si,
                f"{orv:.2f}{sig}",
                ha="center",
                va="center",
                fontsize=7,
                fontweight="bold",
                color=tc,
            )
            ax.text(
                ei,
                si + 0.35,
                f"n={nv}",
                ha="center",
                va="center",
                fontsize=4.5,
                color=tc,
                alpha=0.7,
            )

    ax.set_xticks(range(n_ep))
    ax.set_xticklabels(
        [lbl for _, lbl in ENDPOINTS], fontsize=6, rotation=30, ha="right"
    )
    ax.set_yticks(range(n_egfr))
    ax.set_yticklabels([lbl for _, lbl in EGFR_STRATA], fontsize=6.5)
    ax.set_ylabel("eGFR (mL/min/1.73m²)")

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.3)
        spine.set_color("#999")

    return norm, cmap


# ====================================================================
# MAIN FIGURE
# ====================================================================
def main():
    print("\n── eFig: Secondary endpoints × eGFR strata ──")
    df = load_hte()

    fig = plt.figure(figsize=(7.205, 5.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.3, 1], wspace=0.35)

    # Panel A: Forest
    ax_forest = fig.add_subplot(gs[0])
    make_forest(ax_forest, df)
    ax_forest.text(
        -0.15,
        1.02,
        "a",
        transform=ax_forest.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
    )

    # Panel B: Heatmap
    ax_hm = fig.add_subplot(gs[1])
    norm, cmap = make_heatmap(ax_hm, df)
    ax_hm.text(
        -0.22,
        1.02,
        "b",
        transform=ax_hm.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
    )
    ax_hm.set_title("MIMIC-IV: OR by eGFR × Endpoint", fontsize=7, pad=8)

    # Colorbar
    cbar_ax = fig.add_axes([0.93, 0.15, 0.015, 0.7])
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_ticks([np.log(0.25), np.log(0.5), 0, np.log(2), np.log(4)])
    cbar.set_ticklabels(["0.25", "0.5", "1.0", "2.0", "4.0"])
    cbar.ax.tick_params(labelsize=5.5, width=0.3, length=2)
    cbar.ax.set_ylabel("Odds Ratio", fontsize=6, rotation=270, labelpad=8)
    cbar.outline.set_linewidth(0.3)

    # Footnote
    fig.text(
        0.5,
        -0.02,
        "● P < 0.05  ○ P ≥ 0.05  ◆ Overall  |  "
        "Blue = protective  Red = harmful  |  MIMIC-IV, n = 6,354 pairs",
        ha="center",
        fontsize=5,
        color="#666",
        style="italic",
    )

    save(fig, "efig_secondary_egfr")


if __name__ == "__main__":
    main()
