#!/usr/bin/env python3
"""
04d_fig_secondary_egfr.py — eGFR-stratified multi-endpoint forest + heatmap

Shows the cross-organ eGFR gradient: eGFR≥90 = safe (AKI/mort protective,
rest null), eGFR<45 = multi-organ harm zone.

  Panel A: Forest plot — 6 endpoints × 5 eGFR strata, grouped by endpoint
  Panel B: Heatmap    — OR colored blue/white/red, significance overlay

Run:  python 04d_fig_secondary_egfr.py
"""

import os
import sys

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
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
    ("encephalopathy_delirium", "Enceph / delirium"),
    ("transfusion", "Transfusion"),
    ("reintubation", "Reintubation"),
]

EGFR_STRATA = [
    ("eGFR >= 90", "≥90"),
    ("eGFR 60-89", "60–89"),
    ("eGFR 45-59", "45–59"),
    ("eGFR 30-44", "30–44"),
    ("eGFR < 30", "<30"),
]

EP_COLORS = {
    "aki_7d": "#2166AC",
    "hosp_mortality": "#333333",
    "poaf": "#D6604D",
    "encephalopathy_delirium": "#8856A7",
    "transfusion": "#E08214",
    "reintubation": "#35978F",
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
    return r.iloc[0] if len(r) > 0 else None


# ====================================================================
# PANEL A: FOREST PLOT
# ====================================================================
def make_forest(ax, df):
    ax.axvline(1, color="#ddd", lw=0.6, zorder=0)
    ax.axvspan(1, 12, alpha=0.02, color="#d62728", zorder=0)

    y = 0
    yticks, ylabels = [], []
    group_spans = []

    for ei, (oc, oc_lbl) in enumerate(reversed(ENDPOINTS)):
        color = EP_COLORS[oc]
        y_start = y

        for si, (sg, sg_lbl) in enumerate(reversed(EGFR_STRATA)):
            r = get_row(df, oc, sg)
            if r is not None and not pd.isna(r["or"]) and r["or"] <= 20:
                est, lo, hi = r["or"], r["or_lo"], r["or_hi"]
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
            yticks.append(y)
            ylabels.append(sg_lbl)
            y += 1

        # Overall diamond
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
        ylabels.append("Overall")
        y_end = y

        group_spans.append((y_start, y_end, oc_lbl, color))
        y += 1.5

    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=5.5)
    ax.set_xscale("log")
    ax.set_xlim(0.15, 12)
    ax.set_xticks([0.25, 0.5, 1, 2, 4, 8])
    ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
    ax.set_xlabel("Odds Ratio (95% CI)")

    # Bracket + endpoint labels on right
    for y_start, y_end, lbl, color in group_spans:
        mid = (y_start + y_end) / 2
        x_b = 0.97
        ax.plot(
            [x_b, x_b],
            [y_start - 0.3, y_end + 0.3],
            transform=mpl.transforms.blended_transform_factory(
                ax.transAxes, ax.transData
            ),
            color=color,
            lw=1.5,
            clip_on=False,
            solid_capstyle="round",
        )
        ax.text(
            1.01,
            mid,
            lbl,
            fontsize=6.5,
            fontweight="bold",
            color=color,
            va="center",
            ha="left",
            transform=mpl.transforms.blended_transform_factory(
                ax.transAxes, ax.transData
            ),
            clip_on=False,
        )

    ax.text(
        0.02, -0.04, "← Favors IV Mg", transform=ax.transAxes, fontsize=5, color=GRAY
    )
    ax.text(
        0.90,
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

    # Column labels on top, colored by endpoint
    ep_short = ["AKI", "Mort.", "POAF", "Enceph", "Transfus.", "Reintub."]
    ax.set_xticks(range(n_ep))
    ax.set_xticklabels(ep_short, fontsize=6)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    for ei, (oc, _) in enumerate(ENDPOINTS):
        ax.get_xticklabels()[ei].set_color(EP_COLORS[oc])
        ax.get_xticklabels()[ei].set_fontweight("bold")

    # Row labels — n embedded to avoid colorbar collision
    egfr_labels = [
        f"≥90 (G1)  n={n_mat[0,0]:,}",
        f"60–89 (G2)  n={n_mat[1,0]:,}",
        f"45–59 (G3a)  n={n_mat[2,0]:,}",
        f"30–44 (G3b)  n={n_mat[3,0]:,}",
        f"<30 (G4–5)  n={n_mat[4,0]:,}",
    ]
    ax.set_yticks(range(n_egfr))
    ax.set_yticklabels(egfr_labels, fontsize=5.5)
    ax.set_ylabel("eGFR (mL/min/1.73m²)", fontsize=7)

    # White grid lines
    for i in range(1, n_egfr):
        ax.axhline(i - 0.5, color="white", lw=0.8)
    for i in range(1, n_ep):
        ax.axvline(i - 0.5, color="white", lw=0.8)

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.5)
        spine.set_color("#666")

    return norm, cmap


# ====================================================================
# MAIN
# ====================================================================
def main():
    print("\n── eFig: Secondary endpoints × eGFR strata ──")
    df = load_hte()

    fig = plt.figure(figsize=(7.205, 6.5), constrained_layout=False)
    gs = fig.add_gridspec(
        1,
        2,
        width_ratios=[1.1, 1],
        wspace=0.55,
        left=0.08,
        right=0.88,
        top=0.93,
        bottom=0.07,
    )

    ax_forest = fig.add_subplot(gs[0])
    make_forest(ax_forest, df)
    ax_forest.text(
        -0.12,
        1.03,
        "a",
        transform=ax_forest.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
    )
    ax_forest.set_title("eGFR-stratified ORs by endpoint", fontsize=7.5, pad=10)

    ax_hm = fig.add_subplot(gs[1])
    norm, cmap = make_heatmap(ax_hm, df)
    ax_hm.text(
        -0.28,
        1.03,
        "b",
        transform=ax_hm.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
    )

    cbar_ax = fig.add_axes([0.90, 0.15, 0.015, 0.65])
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_ticks([np.log(0.25), np.log(0.5), 0, np.log(2), np.log(4)])
    cbar.set_ticklabels(["0.25", "0.5", "1.0", "2.0", "4.0"])
    cbar.ax.tick_params(labelsize=5.5, width=0.3, length=2)
    cbar.ax.set_ylabel("Odds Ratio", fontsize=6, rotation=270, labelpad=8)
    cbar.outline.set_linewidth(0.3)

    fig.text(
        0.5,
        0.01,
        "● P < 0.05   ○ P ≥ 0.05   ◆ Overall   |   "
        "Blue = protective   Red = harmful   |   "
        "MIMIC-IV, 6,354 matched pairs   |   * P < 0.05 in heatmap",
        ha="center",
        fontsize=5,
        color="#666",
        style="italic",
    )

    save(fig, "efig_secondary_egfr")


if __name__ == "__main__":
    main()
