#!/usr/bin/env python3
"""
fig_did.py — Publication figures for Mg→AKI DiD study (AIPW framework)

Fig 1:  AIPW + sIPTW_DR time course (2 panels: eICU, MIMIC)
Fig 2:  Sensitivity forest plot (5 models × 2 databases)
Fig 3:  AKI subgroup forest plot (MIMIC KDIGO≥1, sIPTW_DR OR)
Fig S1: Specification curve (256 specs, AIPW, 24h)

Reads: ~/mg_aki/results/did_timecourse_*_primary.csv
       ~/mg_aki/results/did_primary_*_primary.csv (for sensitivity)
       ~/mg_aki/results/did_subgroups_full_mimic.csv
       ~/mg_aki/results/did_sweep.csv
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

WONG = {
    "blue": "#0072B2",
    "vermil": "#D55E00",
    "green": "#009E73",
    "orange": "#E69F00",
    "skyblue": "#56B4E9",
    "purple": "#CC79A7",
    "black": "#000000",
}

RESULTS = os.path.expanduser("~/mg_aki/results")
FIGW_2 = 7.205  # double-column (183mm)
FIGW_1 = 3.504  # single-column (89mm)


def save(fig, name):
    for ext in ["pdf", "png"]:
        out = os.path.join(RESULTS, f"{name}.{ext}")
        fig.savefig(out, format=ext, dpi=300)
    print(f"  Saved: {name}.pdf/.png")


# ============================================================================
# FIG 1: AIPW + sIPTW_DR time course
# ============================================================================
def fig1_timecourse():
    print("── Fig 1: AIPW time course ──")
    fig, axes = plt.subplots(1, 2, figsize=(FIGW_2, 3.0), sharey=True)

    for idx, (db, title) in enumerate([("eicu", "eICU-CRD"), ("mimic", "MIMIC-IV")]):
        ax = axes[idx]
        path = os.path.join(RESULTS, f"did_timecourse_{db}_primary.csv")
        if not os.path.exists(path):
            print(f"  {path} not found")
            continue
        df = pd.read_csv(path).dropna(subset=["aipw_did"])

        x = df["target_h"].values
        y = df["aipw_did"].values
        lo = df["aipw_lo"].values
        hi = df["aipw_hi"].values
        sig = df["aipw_p"].values < 0.05

        # AIPW (primary)
        ax.fill_between(x, lo, hi, color=WONG["blue"], alpha=0.12, zorder=1)
        ax.plot(
            x, y, color=WONG["blue"], lw=1.2, marker="o", ms=4, label="AIPW", zorder=3
        )
        if np.any(sig):
            ax.plot(
                x[sig],
                y[sig],
                "o",
                color=WONG["blue"],
                ms=5,
                markeredgecolor="white",
                markeredgewidth=0.5,
                zorder=4,
            )
            ax.plot(
                x[~sig],
                y[~sig],
                "o",
                color=WONG["blue"],
                ms=4,
                fillstyle="none",
                zorder=4,
            )

        # sIPTW_DR (secondary)
        ys = df["siptw_did"].values
        sig_s = df["siptw_p"].values < 0.05
        ax.plot(
            x,
            ys,
            color=WONG["vermil"],
            lw=0.8,
            ls="--",
            marker="s",
            ms=3,
            label="sIPTW-DR",
            alpha=0.8,
            zorder=2,
        )

        ax.axhline(0, color="grey", lw=0.5, ls="-", zorder=0)
        ax.set_xlabel("Hours from ICU admission")
        ax.set_xticks(range(6, 37, 6))
        ax.set_title(title, fontweight="bold", fontsize=8)
        ax.text(
            -0.12,
            1.05,
            chr(97 + idx),
            transform=ax.transAxes,
            fontsize=8,
            fontweight="bold",
            va="top",
        )

    axes[0].set_ylabel("DiD estimate, ΔΔCr (mg/dL)")
    axes[0].legend(loc="lower left")
    fig.text(
        0.5,
        -0.02,
        "Negative values indicate IV Mg protective (less Cr rise)",
        ha="center",
        fontsize=6,
        fontstyle="italic",
        color="grey",
    )
    save(fig, "fig1_timecourse")
    plt.close(fig)


# ============================================================================
# FIG 2: Sensitivity analysis forest plot
# ============================================================================
def fig2_sensitivity():
    print("── Fig 2: Sensitivity forest ──")

    models = [
        ("primary", "Primary (21 covars, labs only)"),
        ("sens_a", "Sens A (18, K+ Ca only)"),
        ("sens_b", "Sens B (20, +lactate, no Mg)"),
        ("sens_c", "Sens C (26, +drugs+steroids)"),
        ("sens_d", "Sens D (16, base only)"),
    ]

    rows = []
    for mtag, mlabel in models:
        for db, dblabel in [("eicu", "eICU"), ("mimic", "MIMIC")]:
            path = os.path.join(RESULTS, f"did_timecourse_{db}_{mtag}.csv")
            if not os.path.exists(path):
                continue
            df = pd.read_csv(path)
            r24 = df[df.target_h == 24]
            if len(r24) == 0:
                continue
            r = r24.iloc[0]
            rows.append(
                {
                    "model": mlabel,
                    "db": dblabel,
                    "did": r["aipw_did"],
                    "lo": r["aipw_lo"],
                    "hi": r["aipw_hi"],
                    "p": r["aipw_p"],
                    "n": r.get("n_trt", 0),
                }
            )
    if not rows:
        print("  No data")
        return

    rdf = pd.DataFrame(rows)

    fig, axes = plt.subplots(1, 2, figsize=(FIGW_2, 2.8), sharey=True)

    for idx, (db, title) in enumerate([("eICU", "eICU-CRD"), ("MIMIC", "MIMIC-IV")]):
        ax = axes[idx]
        sub = rdf[rdf.db == db].reset_index(drop=True)
        y_pos = np.arange(len(sub))

        for i, r in sub.iterrows():
            color = WONG["blue"] if "Primary" in r["model"] else WONG["black"]
            lw = 1.5 if "Primary" in r["model"] else 0.8
            ms = 6 if "Primary" in r["model"] else 4
            ax.errorbar(
                r["did"],
                i,
                xerr=[[r["did"] - r["lo"]], [r["hi"] - r["did"]]],
                fmt="o",
                color=color,
                ms=ms,
                lw=lw,
                capsize=2,
                capthick=lw,
            )
            sig = "**" if r["p"] < 0.01 else "*" if r["p"] < 0.05 else ""
            ax.text(
                r["hi"] + 0.003,
                i,
                f'P={r["p"]:.3f}{sig}',
                va="center",
                fontsize=5.5,
                color="grey",
            )

        ax.axvline(0, color="grey", lw=0.5, ls="--")
        ax.set_yticks(y_pos)
        ax.set_yticklabels([r["model"] for _, r in sub.iterrows()], fontsize=6)
        ax.set_xlabel("AIPW DiD at 24 h (mg/dL)")
        ax.set_title(title, fontweight="bold", fontsize=8)
        ax.text(
            -0.30,
            1.05,
            chr(97 + idx),
            transform=ax.transAxes,
            fontsize=8,
            fontweight="bold",
            va="top",
        )
        ax.invert_yaxis()

    save(fig, "fig2_sensitivity")
    plt.close(fig)


# ============================================================================
# FIG 3: AKI subgroup forest (MIMIC KDIGO≥1)
# ============================================================================
def fig3_subgroups():
    print("── Fig 3: AKI subgroup forest (MIMIC) ──")
    path = os.path.join(RESULTS, "did_subgroups_full_mimic.csv")
    if not os.path.exists(path):
        print(f"  {path} not found")
        return

    df = pd.read_csv(path)
    k1 = df[df.outcome == "AKI KDIGO>=1"].copy()
    k1 = k1[k1["or"].notna() & (k1.subgroup != "Overall")].reset_index(drop=True)

    # Build labels
    k1["label"] = k1.apply(
        lambda r: f"{r['subgroup']}: {r['level']}"
        + (" (ref)" if r.get("ref") == "ref" else ""),
        axis=1,
    )

    fig, ax = plt.subplots(figsize=(FIGW_1 * 1.5, len(k1) * 0.22 + 0.8))
    y = np.arange(len(k1))

    for i, r in k1.iterrows():
        is_ref = r.get("ref") == "ref"
        color = "grey" if is_ref else (WONG["blue"] if r["or"] < 1 else WONG["vermil"])
        ms = 3 if is_ref else 5
        ax.errorbar(
            r["or"],
            i,
            xerr=[[r["or"] - r["or_lo"]], [r["or_hi"] - r["or"]]],
            fmt="o" if not is_ref else "D",
            color=color,
            ms=ms,
            lw=0.8,
            capsize=1.5,
            capthick=0.5,
        )
        # P value
        if not is_ref and r["p"] < 0.05:
            ax.text(
                max(r["or_hi"], 1.05) + 0.05,
                i,
                "*",
                fontsize=7,
                va="center",
                color=WONG["blue"],
                fontweight="bold",
            )

    ax.axvline(1, color="grey", lw=0.5, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(k1["label"].values, fontsize=5.5)
    ax.set_xlabel("OR (95% CI) for AKI KDIGO ≥ 1")
    ax.set_title("MIMIC-IV: Subgroup analysis", fontweight="bold", fontsize=8)
    ax.set_xlim(0, max(k1["or_hi"].max() * 1.1, 2.5))
    ax.invert_yaxis()

    save(fig, "fig3_aki_subgroups_mimic")
    plt.close(fig)


# ============================================================================
# FIG S1: Specification curve (256 specs, AIPW, 24h)
# ============================================================================
def figs1_speccurve():
    print("── Fig S1: Specification curve ──")
    path = os.path.join(RESULTS, "did_sweep.csv")
    if not os.path.exists(path):
        print(f"  {path} not found")
        return

    df = pd.read_csv(path)
    s24 = df[df.target_h == 24].copy()
    s24 = s24.dropna(subset=["e_aipw", "m_aipw"])

    # Sort by MIMIC AIPW estimate
    s24 = s24.sort_values("m_aipw").reset_index(drop=True)
    x = np.arange(len(s24))

    fig, axes = plt.subplots(
        3,
        1,
        figsize=(FIGW_2, 5.5),
        gridspec_kw={"height_ratios": [2, 2, 1.2]},
        sharex=True,
    )

    # Panel a: eICU AIPW estimates
    ax = axes[0]
    both_neg = (s24.e_aipw < 0) & (s24.m_aipw < 0)
    colors_e = [WONG["blue"] if bn else WONG["vermil"] for bn in both_neg]
    ax.scatter(x, s24.e_aipw, c=colors_e, s=4, alpha=0.6, edgecolors="none")
    ax.axhline(0, color="grey", lw=0.5, ls="--")
    ax.set_ylabel("eICU AIPW DiD")
    ax.set_title(
        "Specification curve: 256 covariate specifications, 24 h, AIPW",
        fontsize=7,
        fontweight="bold",
    )
    ax.text(
        -0.06,
        1.05,
        "a",
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        va="top",
    )

    # Panel b: MIMIC AIPW estimates (sorted)
    ax = axes[1]
    colors_m = [WONG["blue"] if v < 0 else WONG["vermil"] for v in s24.m_aipw]
    ax.scatter(x, s24.m_aipw, c=colors_m, s=4, alpha=0.6, edgecolors="none")
    ax.axhline(0, color="grey", lw=0.5, ls="--")
    ax.set_ylabel("MIMIC AIPW DiD")
    ax.text(
        -0.06,
        1.05,
        "b",
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        va="top",
    )

    # Panel c: Toggle indicators
    ax = axes[2]
    toggle_cols = ["D1", "D2", "D3", "D4", "L1", "L2", "L3", "L4"]
    toggle_labels = [
        "Chronic\ndrugs",
        "Steroids",
        "ICU\ndrugs",
        "Beta\nblockers",
        "K+",
        "Ca",
        "Lactate",
        "Mg",
    ]
    for ti, tc in enumerate(toggle_cols):
        on = s24[tc].values == 1
        ax.scatter(
            x[on], np.full(on.sum(), ti), s=1.5, c=WONG["black"], alpha=0.4, marker="|"
        )

    ax.set_yticks(range(len(toggle_cols)))
    ax.set_yticklabels(toggle_labels, fontsize=5)
    ax.set_xlabel("Specification (sorted by MIMIC estimate)")
    ax.set_ylim(-0.5, len(toggle_cols) - 0.5)
    ax.invert_yaxis()
    ax.text(
        -0.06,
        1.05,
        "c",
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        va="top",
    )

    # Summary text
    n_conc = both_neg.sum()
    fig.text(
        0.02,
        -0.03,
        f"Blue = both databases negative ({n_conc}/256, {100*n_conc/256:.0f}%); "
        f"Red = discordant",
        fontsize=6,
        color="grey",
    )

    save(fig, "figs1_speccurve")
    plt.close(fig)


# ============================================================================
if __name__ == "__main__":
    os.makedirs(RESULTS, exist_ok=True)
    fig1_timecourse()
    fig2_sensitivity()
    fig3_subgroups()
    figs1_speccurve()
    print("\nDone.")
