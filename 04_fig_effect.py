#!/usr/bin/env python3
"""
04_fig_effect.py — Clean effect-size figure (primary analysis)

Shows DiD ΔCr treatment effect over time for the primary spec (19var,
LAST labs, no K/Mg), PSM_DR, both databases side by side.
Single curve per panel with CI band and annotated values.

Usage: python 04_fig_effect.py
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
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "legend.frameon": False,
    "legend.fontsize": 6.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.04,
}.items():
    mpl.rcParams[k] = v

BLUE = "#0072B2"
BLUE_LT = "#56B4E9"
VERMIL = "#D55E00"
GRAY = "#999999"
RESULTS = os.path.expanduser("~/mg_aki/results")
SPEC = "primary"
METHOD = "psm_dr"
DB_LABELS = {"mimic": "MIMIC-IV (BIDMC)", "eicu": "eICU-CRD (208 hospitals)"}


def load(tag):
    p = os.path.join(RESULTS, f"did_riskset_{tag}.csv")
    if not os.path.exists(p):
        return None
    df = pd.read_csv(p)
    if "spec" not in df.columns:
        return None
    return df


def plot_effect(dbs):
    n = len(dbs)
    fig, axes = plt.subplots(
        1, n, figsize=(3.5 * n, 3.0), sharey=True, constrained_layout=True
    )
    if n == 1:
        axes = [axes]

    for i, (tag, df) in enumerate(dbs):
        ax = axes[i]

        # Reference
        ax.axhline(0, color="#aaaaaa", lw=0.5, ls="-", zorder=1)
        ax.axvspan(22, 26, color="#f0f0f0", zorder=0)  # 24h band

        # Yet-untreated (primary)
        yt = (
            df[(df.spec == SPEC) & (df.pool == "yet_untreated") & (df.method == METHOD)]
            .sort_values("target_h")
            .dropna(subset=["did"])
        )

        if len(yt) > 0:
            h = yt.target_h.values
            d = yt.did.values
            lo = yt.ci_lo.values
            hi = yt.ci_hi.values
            pv = yt.p.values

            ax.fill_between(h, lo, hi, alpha=0.15, color=BLUE, zorder=2)
            ax.plot(
                h,
                d,
                color=BLUE,
                lw=1.5,
                marker="o",
                ms=5,
                markerfacecolor="white",
                markeredgewidth=1,
                markeredgecolor=BLUE,
                zorder=4,
                label="Yet-untreated (primary)",
            )

            for j in range(len(h)):
                if not np.isnan(pv[j]) and pv[j] < 0.05:
                    ax.plot(h[j], d[j], "o", ms=5, color=BLUE, zorder=5)

            # Annotate 24h
            r24 = yt[yt.target_h == 24]
            if len(r24) == 1:
                did24 = r24.did.values[0]
                p24 = r24.p.values[0]
                sig = "*" if p24 < 0.05 else ""
                ax.annotate(
                    f"24h: {did24:+.003f}{sig}\n(P={p24:.3f})",
                    xy=(24, did24),
                    xytext=(34, did24 + 0.005),
                    fontsize=6,
                    color=BLUE,
                    ha="left",
                    arrowprops=dict(arrowstyle="->", color=BLUE, lw=0.6),
                )

            # Annotate peak (most negative)
            peak_idx = np.argmin(d)
            if h[peak_idx] != 24:
                ax.annotate(
                    f"{int(h[peak_idx])}h: {d[peak_idx]:+.003f}*",
                    xy=(h[peak_idx], d[peak_idx]),
                    xytext=(h[peak_idx] - 8, d[peak_idx] - 0.005),
                    fontsize=5.5,
                    color=BLUE,
                    ha="center",
                    arrowprops=dict(arrowstyle="->", color=BLUE, lw=0.5),
                )

        # Never-treated (reference, lighter)
        nt = (
            df[(df.spec == SPEC) & (df.pool == "never_treated") & (df.method == METHOD)]
            .sort_values("target_h")
            .dropna(subset=["did"])
        )
        if len(nt) > 0:
            ax.plot(
                nt.target_h.values,
                nt.did.values,
                color=VERMIL,
                lw=0.8,
                ls="--",
                marker="s",
                ms=3.5,
                markerfacecolor="white",
                markeredgecolor=VERMIL,
                markeredgewidth=0.7,
                alpha=0.7,
                zorder=3,
                label="Never-treated (ref)",
            )

        ax.set_xticks([6, 12, 18, 24, 30, 36, 42, 48])
        ax.set_xlim(3, 51)
        ax.set_xlabel("Hours from T₀ (first IV Mg)")
        if i == 0:
            ax.set_ylabel("Treatment effect: DiD ΔCr (mg/dL)")
            ax.legend(loc="lower left", handlelength=2.5)

        ax.text(
            -0.12,
            1.06,
            chr(ord("a") + i),
            transform=ax.transAxes,
            fontsize=10,
            fontweight="bold",
            va="top",
        )
        ax.set_title(DB_LABELS.get(tag, tag), fontsize=8, pad=6)

    # Footer
    fig.text(
        0.5,
        -0.02,
        "Primary: 19-var PS (no K⁺/Mg), LAST labs, PSM+DR (HC1 SE)  |  "
        "Filled ● = P<0.05  |  DiD<0 = renoprotective",
        ha="center",
        fontsize=5.5,
        color="#666666",
    )

    out = os.path.join(RESULTS, "fig_effect.pdf")
    fig.savefig(out, format="pdf")
    fig.savefig(out.replace(".pdf", ".png"), format="png")
    plt.close()
    print(f"  ✓ {out} (.pdf + .png)")


if __name__ == "__main__":
    print("=" * 70)
    print("04_fig_effect.py — Primary effect-size figure")
    print("=" * 70)

    dbs = []
    for tag in ["mimic", "eicu"]:
        df = load(tag)
        if df is not None:
            dbs.append((tag, df))
            print(f"  {tag}: loaded")
        else:
            print(f"  {tag}: no data")

    if dbs:
        plot_effect(dbs)
    print("=" * 70)
