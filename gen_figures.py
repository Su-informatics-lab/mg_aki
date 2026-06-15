#!/usr/bin/env python3
"""
gen_figures.py — All publication figures (Nature/JNO style)

  Figure 2:  Forest plot (AC primary + sensitivity + controls)
  Figure 3:  Surgery-type interaction (cardioplegia hypothesis)
  eFigure 1: Subgroup forest (who benefits: surgery, age, Mg, eGFR, PPI)
  eFigure 2: PS overlap density (both databases)

Reads: results/02_results.csv, results/etables_4_5_subgroups.csv
Writes: figs/fig2_forest.pdf, figs/fig3_interaction.pdf,
        figs/efig1_subgroups.pdf, figs/efig2_ps_overlap.pdf

Run: python gen_figures.py
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import NullLocator

# =====================================================================
# Nature Portfolio rcParams
# =====================================================================
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

# Wong/Okabe-Ito palette
C_BLUE = "#0072B2"  # eICU
C_VERMILLION = "#D55E00"  # MIMIC
C_BLACK = "#000000"  # Pooled
C_SKYBLUE = "#56B4E9"  # secondary
C_GREEN = "#009E73"  # tertiary
C_GRAY = "#999999"

RESULTS = "results"
FIGS = "figs"
os.makedirs(FIGS, exist_ok=True)

# Double-column width for main figures
W_DOUBLE = 7.205  # inches (183mm)
W_SINGLE = 3.504  # inches (89mm)
W_ONEHALF = 4.724  # inches (120mm)


def save(fig, name):
    for ext in ["pdf", "png"]:
        fig.savefig(os.path.join(FIGS, f"{name}.{ext}"), format=ext)
    print(f"  Saved: figs/{name}.pdf + .png")
    plt.close(fig)


# =====================================================================
# FIGURE 2: Forest Plot (primary + sensitivity + controls)
# =====================================================================
def fig2_forest():
    print("Figure 2: Forest plot")
    res = pd.read_csv(os.path.join(RESULTS, "02_results.csv"))
    # Use m=5 (first/primary)
    m_val = res["m"].dropna().min() if "m" in res.columns else None
    if m_val is not None:
        res = res[(res["m"].isna()) | (res["m"] == m_val)]

    # Define sections (top to bottom in the figure)
    sections = [
        (
            "Primary: Active Comparator",
            [
                ("AKI KDIGO ≥1 (Mg+K⁺ vs K⁺-only)", "ac_aki1"),
            ],
        ),
        (
            "Sensitivity: All-Patient AKI",
            [
                ("IPTW", "iptw_aki1"),
                ("Overlap weighting", "ow_aki1"),
                ("PS matching", "psm_aki1"),
            ],
        ),
        (
            "Exploratory",
            [
                ("Hospital mortality", "ow_mort"),
                ("Encephalopathy", "ow_enceph"),
            ],
        ),
        (
            "Control",
            [
                ("Fracture (negative control)", "ow_frac"),
            ],
        ),
    ]

    # Build rows (from bottom to top for matplotlib y-axis)
    rows = []
    y = 0
    for sec_label, outcomes in reversed(sections):
        for label, key in reversed(outcomes):
            e = res[(res["db"] == "eICU") & (res["analysis"] == key)]
            m = res[(res["db"] == "MIMIC") & (res["analysis"] == key)]
            p = res[(res["db"] == "Pooled") & (res["analysis"] == key)]
            if len(p) == 0:
                continue

            # Pooled (diamond)
            rows.append(
                dict(
                    y=y,
                    label="  Pooled",
                    or_=p.or_.values[0] if "or_" in p.columns else p["or"].values[0],
                    lo=p["lo"].values[0],
                    hi=p["hi"].values[0],
                    source="Pooled",
                    section=sec_label,
                    is_header=False,
                )
            )
            y += 1

            # MIMIC
            if len(m) > 0:
                or_val = m["or_"].values[0] if "or_" in m.columns else m["or"].values[0]
                rows.append(
                    dict(
                        y=y,
                        label="  MIMIC-IV",
                        or_=or_val,
                        lo=m["lo"].values[0],
                        hi=m["hi"].values[0],
                        source="MIMIC-IV",
                        section=sec_label,
                        is_header=False,
                    )
                )
                y += 1

            # eICU
            if len(e) > 0:
                or_val = e["or_"].values[0] if "or_" in e.columns else e["or"].values[0]
                rows.append(
                    dict(
                        y=y,
                        label="  eICU-CRD",
                        or_=or_val,
                        lo=e["lo"].values[0],
                        hi=e["hi"].values[0],
                        source="eICU",
                        section=sec_label,
                        is_header=False,
                    )
                )
                y += 1

            # Outcome header (no point)
            rows.append(
                dict(
                    y=y,
                    label=label,
                    or_=np.nan,
                    lo=np.nan,
                    hi=np.nan,
                    source="header",
                    section=sec_label,
                    is_header=True,
                )
            )
            y += 1

        # Section header
        rows.append(
            dict(
                y=y,
                label=sec_label,
                or_=np.nan,
                lo=np.nan,
                hi=np.nan,
                source="section",
                section=sec_label,
                is_header=True,
            )
        )
        y += 1.5  # extra space between sections

    fd = pd.DataFrame(rows)

    fig, ax = plt.subplots(figsize=(W_DOUBLE, 5.5))

    # Null reference line
    ax.axvline(1, color=C_GRAY, linestyle="--", linewidth=0.5, zorder=0)

    # Plot points and CIs
    for _, r in fd.iterrows():
        if np.isnan(r["or_"]):
            continue
        color = {"eICU": C_BLUE, "MIMIC-IV": C_VERMILLION, "Pooled": C_BLACK}.get(
            r["source"], C_GRAY
        )
        marker = "D" if r["source"] == "Pooled" else "s"
        ms = 7 if r["source"] == "Pooled" else 5

        # CI line
        ax.plot(
            [r["lo"], r["hi"]],
            [r["y"], r["y"]],
            color=color,
            linewidth=1.0,
            solid_capstyle="round",
        )
        # Point
        ax.plot(
            r["or_"],
            r["y"],
            marker=marker,
            color=color,
            markersize=ms,
            markeredgewidth=0,
            zorder=5,
        )

        # OR text on right
        ax.text(
            2.8,
            r["y"],
            f"{r['or_']:.2f} ({r['lo']:.2f}–{r['hi']:.2f})",
            va="center",
            ha="left",
            fontsize=5.5,
            color="#333333",
        )

    # Labels on left
    for _, r in fd.iterrows():
        weight = "bold" if r["source"] in ("section", "header", "Pooled") else "normal"
        size = 7 if r["source"] == "section" else 6.5
        ax.text(
            0.30,
            r["y"],
            r["label"],
            va="center",
            ha="left",
            fontsize=size,
            fontweight=weight,
        )

    ax.set_xscale("log")
    ax.set_xlim(0.25, 4.0)
    ax.set_xticks([0.5, 0.75, 1, 1.5, 2])
    ax.set_xticklabels(["0.50", "0.75", "1.00", "1.50", "2.00"])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("Odds ratio (95% CI)")
    ax.set_ylim(-1.5, fd["y"].max() + 1)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)

    # Direction annotations
    ax.text(
        0.55,
        -1.2,
        "← Favors supplementation",
        fontsize=5.5,
        ha="center",
        color=C_GRAY,
        style="italic",
    )
    ax.text(
        1.6,
        -1.2,
        "Favors no supplementation →",
        fontsize=5.5,
        ha="center",
        color=C_GRAY,
        style="italic",
    )

    # Legend
    from matplotlib.lines import Line2D

    legend_elements = [
        Line2D(
            [0],
            [0],
            marker="s",
            color="w",
            markerfacecolor=C_BLUE,
            markersize=6,
            label="eICU-CRD",
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            color="w",
            markerfacecolor=C_VERMILLION,
            markersize=6,
            label="MIMIC-IV",
        ),
        Line2D(
            [0],
            [0],
            marker="D",
            color="w",
            markerfacecolor=C_BLACK,
            markersize=6,
            label="Pooled",
        ),
    ]
    ax.legend(
        handles=legend_elements,
        loc="lower right",
        fontsize=6,
        handletextpad=0.3,
        borderpad=0.2,
    )

    save(fig, "fig2_forest")


# =====================================================================
# FIGURE 3: Surgery-Type Interaction (cardioplegia hypothesis)
# =====================================================================
def efig1_interaction():
    print("eFigure 1: Surgery-type interaction")

    data = pd.DataFrame(
        {
            "database": ["eICU-CRD"] * 2 + ["MIMIC-IV"] * 2,
            "surgery": ["Simple\n(CABG/other)", "Complex\n(valve/combined)"] * 2,
            "or": [1.36, 1.74, 1.01, 1.09],
            "lo": [1.19, 1.45, 0.85, 0.85],
            "hi": [1.55, 2.11, 1.19, 1.41],
        }
    )

    fig, ax = plt.subplots(figsize=(W_SINGLE, W_SINGLE * 0.9))

    ax.axhline(1, color=C_GRAY, linestyle="--", linewidth=0.5, zorder=0)

    offsets = {"eICU-CRD": -0.12, "MIMIC-IV": 0.12}
    colors = {"eICU-CRD": C_BLUE, "MIMIC-IV": C_VERMILLION}
    markers = {"eICU-CRD": "s", "MIMIC-IV": "^"}

    for db in ["eICU-CRD", "MIMIC-IV"]:
        sub = data[data.database == db]
        x = np.arange(len(sub)) + offsets[db]
        ax.errorbar(
            x,
            sub["or"],
            yerr=[sub["or"] - sub["lo"], sub["hi"] - sub["or"]],
            fmt=markers[db],
            color=colors[db],
            markersize=7,
            capsize=3,
            capthick=1.0,
            linewidth=1.0,
            label=db,
            markeredgewidth=0,
        )

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Simple\n(CABG / other)", "Complex\n(valve / combined)"])
    ax.set_ylabel("OR per 1 mg/dL serum Mg increase\n(prognostic association with AKI)")
    ax.set_ylim(0.7, 2.3)
    ax.set_yticks(np.arange(0.8, 2.4, 0.2))
    ax.legend(loc="upper left", fontsize=6)

    # Annotation: cardioplegia arrow
    # eICU Complex is at x = 1 + offset(-0.12) = 0.88, OR = 1.74
    ax.annotate(
        "More cardioplegia →\nhigher Mg + more AKI",
        xy=(0.88, 1.74),
        xytext=(0.40, 2.10),
        fontsize=5.5,
        color=C_GRAY,
        style="italic",
        arrowprops=dict(arrowstyle="->", color=C_GRAY, lw=0.5),
    )

    save(fig, "efig1_interaction")


# =====================================================================
# eFIGURE 1: Subgroup Forest (downstream analysis — who benefits)
# =====================================================================
def fig3_subgroups():
    print("Figure 3: Subgroup forest")

    sub = pd.read_csv(os.path.join(RESULTS, "etables_4_5_subgroups.csv"))
    # Focus on eICU AC results (where signal is clearest)
    ac = sub[(sub["db"] == "eICU") & (sub["analysis"] == "ac_ow")].copy()

    # Define display order and nice labels (top to bottom)
    display = [
        ("5d", "eGFR >=60", "eGFR ≥60"),
        ("5d", "eGFR <60", "eGFR <60"),
        ("5d", "No PPI", "No PPI"),
        ("5d", "PPI user", "PPI user"),
        ("5c", "Mg >=2.0", "Baseline Mg ≥2.0"),
        ("5c", "Mg <2.0 (hypo)", "Baseline Mg <2.0"),
        ("5b", "Age >=60", "Age ≥60"),
        ("5b", "Age <60", "Age <60"),
        ("5", "MDS 0-1", "MDS 0–1"),
        ("5", "MDS >=2", "MDS ≥2"),
        ("4", "Valve", "Valve surgery"),
        ("4", "CABG", "CABG"),
        ("4", "Other cardiac", "Other cardiac"),
    ]

    fig, ax = plt.subplots(figsize=(W_ONEHALF, 4.5))
    ax.axvline(1, color=C_GRAY, linestyle="--", linewidth=0.5, zorder=0)

    y_positions = []
    y_labels = []
    y = 0
    prev_table = None

    for etable, subgroup, nice_label in reversed(display):
        row = ac[(ac["etable"] == etable) & (ac["subgroup"] == subgroup)]

        # Add space between subgroup pairs
        if prev_table is not None and etable != prev_table:
            y += 0.6

        if len(row) == 0:
            y_positions.append(y)
            y_labels.append(nice_label)
            y += 1
            prev_table = etable
            continue

        r = row.iloc[0]
        # Color by direction
        sig = r["p"] < 0.05
        color = C_BLUE if sig else C_GRAY

        ax.plot([r["lo"], r["hi"]], [y, y], color=color, linewidth=1.0)
        ax.plot(r["or"], y, "s", color=color, markersize=6, markeredgewidth=0, zorder=5)

        # OR text with sample size
        star = "*" if sig else ""
        ax.text(
            2.5,
            y,
            f"{r['or']:.2f} ({r['lo']:.2f}–{r['hi']:.2f}){star}  n={int(r['n'])}",
            va="center",
            ha="left",
            fontsize=5.5,
            color="#333333",
        )

        y_positions.append(y)
        y_labels.append(nice_label)
        y += 1
        prev_table = etable

    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels, fontsize=6)
    ax.set_xscale("log")
    ax.set_xlim(0.2, 3.5)
    ax.set_xticks([0.5, 0.75, 1, 1.5, 2])
    ax.set_xticklabels(["0.50", "0.75", "1.00", "1.50", "2.00"])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("Odds ratio (95% CI)")
    ax.set_ylim(-1, y + 0.5)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)

    ax.set_title(
        "eICU Active-Comparator Subgroup Analysis",
        fontsize=7,
        fontweight="bold",
        loc="left",
    )

    ax.text(0.55, -0.7, "← Favors Mg+K⁺", fontsize=5, color=C_GRAY, style="italic")
    ax.text(1.5, -0.7, "Favors K⁺-only →", fontsize=5, color=C_GRAY, style="italic")

    save(fig, "fig3_subgroups")


# =====================================================================
# eFIGURE 2: PS Overlap Density (both databases)
# =====================================================================
def efig2_ps_overlap():
    print("eFigure 2: PS overlap (skipped — needs weighted cohort CSVs with PS column)")
    print("  To generate: add PS column export to 02_analysis.R")


# =====================================================================
# MAIN
# =====================================================================
if __name__ == "__main__":
    print("=" * 50)
    print("Generating publication figures")
    print("=" * 50)
    fig2_forest()
    efig1_interaction()
    fig3_subgroups()
    efig2_ps_overlap()
    print("\nDone. Check figs/ directory.")
