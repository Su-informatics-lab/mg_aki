#!/usr/bin/env python3
"""
gen_figures.py — Publication figures (Nature Portfolio style)

  Figure 1:  Forest plot (AC primary + sensitivity + controls)
             → figs/fig2_forest.pdf
  Figure 2:  Mg-stratified treatment effect (NEW — core finding)
             → figs/fig_mg_stratified.pdf
  Figure 3:  Subgroup forest (eICU AC: surgery, age, Mg, eGFR, PPI)
             → figs/fig3_subgroups.pdf
  eFigure 2: Surgery-type interaction (cardioplegia hypothesis)
             → figs/efig1_interaction.pdf

Reads:  results/02_results.csv
        results/08_mg_stratified.csv
        results/etables_4_5_subgroups.csv

Run: python gen_figures.py
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import NullLocator

# =====================================================================
# Nature Portfolio rcParams (mandatory)
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
C_BLACK = "#000000"  # Pooled / reference
C_SKYBLUE = "#56B4E9"
C_GREEN = "#009E73"
C_GRAY = "#999999"

RESULTS = os.path.expanduser("~/mg_aki/results")
FIGS = os.path.expanduser("~/mg_aki/figs")
os.makedirs(FIGS, exist_ok=True)

# Nature widths (inches)
W_DOUBLE = 7.205  # 183 mm
W_SINGLE = 3.504  # 89 mm
W_ONEHALF = 4.724  # 120 mm


def save(fig, name):
    """Save figure as PDF + PNG."""
    for ext in ["pdf", "png"]:
        fig.savefig(os.path.join(FIGS, f"{name}.{ext}"), format=ext)
    print(f"  Saved: figs/{name}.pdf + .png")
    plt.close(fig)


def add_panel_label(ax, label, x=-0.12, y=1.06):
    """Nature-style panel label: lowercase bold, 8pt."""
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        va="top",
        ha="right",
    )


# =====================================================================
# FIGURE 1: Forest Plot (AC primary + sensitivity + controls)
# Output: figs/fig2_forest.pdf  (matches \includegraphics{fig2_forest})
# =====================================================================
def fig1_forest():
    print("Figure 1: Forest plot (AC + sensitivity + controls)")
    res = pd.read_csv(os.path.join(RESULTS, "02_results.csv"))

    or_col = "or" if "or" in res.columns else "or_"

    # Sections (top to bottom in the figure)
    sections = [
        (
            "Primary: Active Comparator",
            [
                ("AKI KDIGO \u22651 (Mg+K\u207a vs K\u207a-only)", "ac_aki1"),
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

    # Build rows bottom-to-top for matplotlib y-axis
    rows = []
    y = 0
    for sec_label, outcomes in reversed(sections):
        for label, key in reversed(outcomes):
            e = res[(res["db"] == "eICU") & (res["analysis"] == key)]
            m = res[(res["db"] == "MIMIC") & (res["analysis"] == key)]
            p = res[(res["db"] == "Pooled") & (res["analysis"] == key)]
            if len(p) == 0:
                continue

            # Pooled
            rows.append(
                dict(
                    y=y,
                    label="  Pooled",
                    or_=p[or_col].values[0],
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
                rows.append(
                    dict(
                        y=y,
                        label="  MIMIC-IV",
                        or_=m[or_col].values[0],
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
                rows.append(
                    dict(
                        y=y,
                        label="  eICU-CRD",
                        or_=e[or_col].values[0],
                        lo=e["lo"].values[0],
                        hi=e["hi"].values[0],
                        source="eICU",
                        section=sec_label,
                        is_header=False,
                    )
                )
                y += 1
            # Outcome header
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
        y += 1.5

    fd = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(W_DOUBLE, 4.8))

    ax.axvline(1, color=C_GRAY, linestyle="--", linewidth=0.5, zorder=0)

    for _, r in fd.iterrows():
        if np.isnan(r["or_"]):
            continue
        color = {"eICU": C_BLUE, "MIMIC-IV": C_VERMILLION, "Pooled": C_BLACK}.get(
            r["source"], C_GRAY
        )
        marker = "D" if r["source"] == "Pooled" else "s"
        ms = 7 if r["source"] == "Pooled" else 5

        ax.plot(
            [r["lo"], r["hi"]],
            [r["y"], r["y"]],
            color=color,
            linewidth=1.0,
            solid_capstyle="round",
        )
        ax.plot(
            r["or_"],
            r["y"],
            marker=marker,
            color=color,
            markersize=ms,
            markeredgewidth=0,
            zorder=5,
        )
        ax.text(
            2.8,
            r["y"],
            f"{r['or_']:.2f} ({r['lo']:.2f}\u2013{r['hi']:.2f})",
            va="center",
            ha="left",
            fontsize=5.5,
            color="#333333",
        )

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

    ax.text(
        0.55,
        -1.2,
        "\u2190 Favors supplementation",
        fontsize=5.5,
        ha="center",
        color=C_GRAY,
        style="italic",
    )
    ax.text(
        1.6,
        -1.2,
        "Favors no supplementation \u2192",
        fontsize=5.5,
        ha="center",
        color=C_GRAY,
        style="italic",
    )

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
# FIGURE 2: Mg-Stratified Treatment Effect (THE CORE FINDING)
# Output: figs/fig_mg_stratified.pdf
#
# Design:
#   Panel a: 4 serum Mg strata × 2 databases (forest)
#   Panel b: Narrow sub-bands within >2.3 (eICU only) + Khalili ref
# =====================================================================
def fig_mg_stratified():
    print("Figure 2: Mg-stratified forest (core finding)")

    csv_path = os.path.join(RESULTS, "08_mg_stratified.csv")

    # ── Try reading from pipeline CSV ────────────────────────────
    if os.path.exists(csv_path):
        raw = pd.read_csv(csv_path)
        or_col = "or" if "or" in raw.columns else "or_"

        def _get(db, stratum, analysis="all_ow"):
            r = raw[
                (raw["db"] == db)
                & (raw["stratum"] == stratum)
                & (raw["analysis"] == analysis)
            ]
            if len(r) == 0:
                return None
            r = r.iloc[0]
            return dict(or_=r[or_col], lo=r["lo"], hi=r["hi"], p=r["p"])

        strata_eicu = [_get("eICU", s) for s in ["<1.8", "1.8-2.0", "2.0-2.3", ">2.3"]]
        strata_mimic = [
            _get("MIMIC", s) for s in ["<1.8", "1.8-2.0", "2.0-2.3", ">2.3"]
        ]

        # Narrow sub-bands from hospital_re CSV or stratified CSV
        sub_bands = []
        re_path = os.path.join(RESULTS, "08b_hospital_re.csv")
        if os.path.exists(re_path):
            re = pd.read_csv(re_path)
            re_or_col = "or" if "or" in re.columns else "or_"
            for sb in [">2.3:2.3-2.6", ">2.3:2.6-3.0", ">2.3:>3.0"]:
                r = re[(re["stratum"] == sb) & (re["model"] == "fixed_subband")]
                if len(r) > 0:
                    r = r.iloc[0]
                    sub_bands.append(
                        dict(or_=r[re_or_col], lo=r["lo"], hi=r["hi"], p=r["p"])
                    )
                else:
                    sub_bands.append(None)
        else:
            sub_bands = [None, None, None]

        # Interaction P values
        int_eicu = _get("eICU", "interaction", "trt_x_mg")
        int_mimic = _get("MIMIC", "interaction", "trt_x_mg")
        int_p_eicu = int_eicu["p"] if int_eicu else 0.005
        int_p_mimic = int_mimic["p"] if int_mimic else 0.12
    else:
        # ── Fallback: hardcoded from STATE_DUMP_R4 ───────────────
        print("  (using hardcoded values — 08_mg_stratified.csv not found)")
        strata_eicu = [
            dict(or_=0.98, lo=0.66, hi=1.46, p=0.92),
            dict(or_=0.96, lo=0.73, hi=1.26, p=0.76),
            dict(or_=0.91, lo=0.62, hi=1.33, p=0.61),
            dict(or_=0.53, lo=0.35, hi=0.80, p=0.003),
        ]
        strata_mimic = [
            dict(or_=1.07, lo=0.80, hi=1.44, p=0.64),
            dict(or_=0.80, lo=0.54, hi=1.19, p=0.27),
            dict(or_=0.77, lo=0.33, hi=1.80, p=0.55),
            dict(or_=0.92, lo=0.58, hi=1.47, p=0.73),
        ]
        sub_bands = [
            dict(or_=0.83, lo=0.49, hi=1.39, p=0.47),
            dict(or_=0.35, lo=0.16, hi=0.76, p=0.008),
            dict(or_=0.34, lo=0.12, hi=0.97, p=0.044),
        ]
        int_p_eicu = 0.005
        int_p_mimic = 0.12

    # Khalili RCT reference (target ~3.0 mg/dL, adj OR 0.26)
    KHALILI_OR = 0.26

    strata_labels = ["<1.8", "1.8\u20132.0", "2.0\u20132.3", ">2.3"]
    sub_labels = ["2.3\u20132.6", "2.6\u20133.0", ">3.0"]

    # ── Build figure: two panels ─────────────────────────────────
    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(W_DOUBLE, 3.8), gridspec_kw={"width_ratios": [1.3, 1]}
    )

    # ── Panel a: 4 strata × 2 databases ─────────────────────────
    ax = ax_a
    ax.axvline(1, color=C_GRAY, linestyle="--", linewidth=0.5, zorder=0)

    y_positions = []
    y_labels_list = []
    y = 0
    offset_db = 0.18  # vertical offset between databases within a stratum

    for i, (s_label, s_e, s_m) in enumerate(
        zip(reversed(strata_labels), reversed(strata_eicu), reversed(strata_mimic))
    ):

        if i > 0:
            y += 0.4  # gap between strata

        # MIMIC (bottom within stratum)
        if s_m is not None:
            ax.plot(
                [s_m["lo"], s_m["hi"]],
                [y, y],
                color=C_VERMILLION,
                linewidth=1.0,
                solid_capstyle="round",
            )
            ax.plot(
                s_m["or_"],
                y,
                "s",
                color=C_VERMILLION,
                markersize=5,
                markeredgewidth=0,
                zorder=5,
            )
            star = "\u2217" if s_m["p"] < 0.05 else ""
            ax.text(
                2.2,
                y,
                f"{s_m['or_']:.2f} ({s_m['lo']:.2f}\u2013{s_m['hi']:.2f}){star}",
                va="center",
                ha="left",
                fontsize=5,
                color=C_VERMILLION,
            )
        y_mimic = y
        y += 2 * offset_db

        # eICU (top within stratum)
        if s_e is not None:
            ax.plot(
                [s_e["lo"], s_e["hi"]],
                [y, y],
                color=C_BLUE,
                linewidth=1.0,
                solid_capstyle="round",
            )
            ax.plot(
                s_e["or_"],
                y,
                "s",
                color=C_BLUE,
                markersize=5,
                markeredgewidth=0,
                zorder=5,
            )
            star = "\u2217" if s_e["p"] < 0.05 else ""
            ax.text(
                2.2,
                y,
                f"{s_e['or_']:.2f} ({s_e['lo']:.2f}\u2013{s_e['hi']:.2f}){star}",
                va="center",
                ha="left",
                fontsize=5,
                color=C_BLUE,
            )
        y_eicu = y

        # Label at midpoint
        y_mid = (y_mimic + y_eicu) / 2
        y_positions.append(y_mid)
        y_labels_list.append(f"Mg {s_label}")

        y += 1

    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels_list, fontsize=6)
    ax.set_xscale("log")
    ax.set_xlim(0.25, 3.5)
    ax.set_xticks([0.25, 0.5, 1, 2])
    ax.set_xticklabels(["0.25", "0.50", "1.00", "2.00"])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("Odds ratio (95% CI)")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)

    # Interaction P annotation
    int_text = f"Interaction P = {int_p_eicu:.3f} (eICU), " f"{int_p_mimic:.2f} (MIMIC)"
    ax.text(
        0.03,
        0.06,
        int_text,
        transform=ax.transAxes,
        fontsize=5.5,
        color="#555555",
        style="italic",
    )

    # Direction labels
    ax.text(
        0.55,
        -1.0,
        "\u2190 Favors Mg",
        fontsize=5,
        color=C_GRAY,
        style="italic",
        ha="center",
    )
    ax.text(
        1.5,
        -1.0,
        "Favors no Mg \u2192",
        fontsize=5,
        color=C_GRAY,
        style="italic",
        ha="center",
    )

    # Legend
    legend_a = [
        Line2D(
            [0],
            [0],
            marker="s",
            color="w",
            markerfacecolor=C_BLUE,
            markersize=5,
            label="eICU-CRD",
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            color="w",
            markerfacecolor=C_VERMILLION,
            markersize=5,
            label="MIMIC-IV",
        ),
    ]
    ax.legend(
        handles=legend_a,
        loc="upper left",
        fontsize=5.5,
        handletextpad=0.2,
        borderpad=0.2,
    )

    add_panel_label(ax, "a", x=-0.08, y=1.04)

    # ── Panel b: Narrow sub-bands within >2.3 (eICU only) ───────
    ax = ax_b
    ax.axvline(1, color=C_GRAY, linestyle="--", linewidth=0.5, zorder=0)

    # Khalili reference line
    ax.axvline(
        KHALILI_OR, color=C_GREEN, linestyle=":", linewidth=0.8, zorder=0, alpha=0.7
    )
    ax.text(
        KHALILI_OR * 0.85,
        -0.4,
        f"Khalili RCT\nOR {KHALILI_OR}",
        fontsize=5,
        color=C_GREEN,
        ha="right",
        va="top",
        style="italic",
    )

    y_sub = []
    y_sub_labels = []
    y = 0
    for i, (sb_label, sb) in enumerate(zip(reversed(sub_labels), reversed(sub_bands))):
        if sb is None:
            y_sub.append(y)
            y_sub_labels.append(sb_label)
            y += 1.2
            continue

        color = C_BLUE if sb["p"] >= 0.05 else C_BLUE
        lw = 1.5 if sb["p"] < 0.05 else 1.0

        ax.plot(
            [sb["lo"], sb["hi"]],
            [y, y],
            color=C_BLUE,
            linewidth=lw,
            solid_capstyle="round",
        )
        ax.plot(
            sb["or_"], y, "s", color=C_BLUE, markersize=6, markeredgewidth=0, zorder=5
        )

        star = "\u2217" if sb["p"] < 0.05 else ""
        ax.text(
            2.5,
            y,
            f"{sb['or_']:.2f} ({sb['lo']:.2f}\u2013{sb['hi']:.2f}){star}",
            va="center",
            ha="left",
            fontsize=5.5,
            color="#333333",
        )

        y_sub.append(y)
        y_sub_labels.append(f"Mg {sb_label}")
        y += 1.2

    ax.set_yticks(y_sub)
    ax.set_yticklabels(y_sub_labels, fontsize=6)
    ax.set_xscale("log")
    ax.set_xlim(0.08, 4.0)
    ax.set_xticks([0.1, 0.25, 0.5, 1, 2])
    ax.set_xticklabels(["0.10", "0.25", "0.50", "1.00", "2.00"])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("Odds ratio (95% CI)")
    ax.set_ylim(-0.8, y - 0.2)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)

    # Sub-band title
    ax.text(
        0.5,
        1.02,
        "Sub-bands within >2.3 mg/dL (eICU)",
        transform=ax.transAxes,
        fontsize=6,
        ha="center",
        fontweight="bold",
    )

    add_panel_label(ax, "b", x=-0.06, y=1.04)

    fig.align_labels()
    save(fig, "fig_mg_stratified")


# =====================================================================
# FIGURE 3: Subgroup Forest (eICU AC — who benefits)
# Output: figs/fig3_subgroups.pdf
# =====================================================================
def fig3_subgroups():
    print("Figure 3: Subgroup forest (eICU AC)")

    csv_path = os.path.join(RESULTS, "etables_4_5_subgroups.csv")
    if not os.path.exists(csv_path):
        print(f"  SKIP: {csv_path} not found")
        return

    sub = pd.read_csv(csv_path)
    or_col = "or" if "or" in sub.columns else "or_"
    ac = sub[(sub["db"] == "eICU") & (sub["analysis"] == "ac_ow")].copy()

    # Display order (top to bottom on figure, reversed for y-axis)
    display = [
        ("5d", "eGFR >=60", "eGFR \u226560"),
        ("5d", "eGFR <60", "eGFR <60"),
        ("5d", "No PPI", "No PPI"),
        ("5d", "PPI user", "PPI user"),
        ("5c", "Mg >=2.0", "Baseline Mg \u22652.0"),
        ("5c", "Mg <2.0 (hypo)", "Baseline Mg <2.0"),
        ("5b", "Age >=60", "Age \u226560"),
        ("5b", "Age <60", "Age <60"),
        ("5", "MDS 0-1", "MDS 0\u20131"),
        ("5", "MDS >=2", "MDS \u22652"),
        ("4", "Valve", "Valve surgery"),
        ("4", "CABG", "CABG"),
        ("4", "Other cardiac", "Other cardiac"),
    ]

    fig, ax = plt.subplots(figsize=(W_ONEHALF, 4.5))
    ax.axvline(1, color=C_GRAY, linestyle="--", linewidth=0.5, zorder=0)

    y_positions = []
    y_labels_list = []
    y = 0
    prev_table = None

    for etable, subgroup, nice_label in reversed(display):
        row = ac[(ac["etable"] == etable) & (ac["subgroup"] == subgroup)]

        if prev_table is not None and etable != prev_table:
            y += 0.6

        if len(row) == 0:
            y_positions.append(y)
            y_labels_list.append(nice_label)
            y += 1
            prev_table = etable
            continue

        r = row.iloc[0]
        sig = r["p"] < 0.05
        color = C_BLUE if sig else C_GRAY

        ax.plot([r["lo"], r["hi"]], [y, y], color=color, linewidth=1.0)
        ax.plot(
            r[or_col], y, "s", color=color, markersize=6, markeredgewidth=0, zorder=5
        )

        star = "\u2217" if sig else ""
        ax.text(
            2.5,
            y,
            f"{r[or_col]:.2f} ({r['lo']:.2f}\u2013{r['hi']:.2f}){star}  "
            f"n={int(r['n'])}",
            va="center",
            ha="left",
            fontsize=5.5,
            color="#333333",
        )

        y_positions.append(y)
        y_labels_list.append(nice_label)
        y += 1
        prev_table = etable

    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels_list, fontsize=6)
    ax.set_xscale("log")
    ax.set_xlim(0.2, 3.5)
    ax.set_xticks([0.5, 0.75, 1, 1.5, 2])
    ax.set_xticklabels(["0.50", "0.75", "1.00", "1.50", "2.00"])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("Odds ratio (95% CI)")
    ax.set_ylim(-1, y + 0.5)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)

    ax.text(
        0.55, -0.7, "\u2190 Favors Mg+K\u207a", fontsize=5, color=C_GRAY, style="italic"
    )
    ax.text(
        1.5,
        -0.7,
        "Favors K\u207a-only \u2192",
        fontsize=5,
        color=C_GRAY,
        style="italic",
    )

    save(fig, "fig3_subgroups")


# =====================================================================
# eFIGURE 2: Surgery-Type Interaction (cardioplegia hypothesis)
# Output: figs/efig1_interaction.pdf  (matches \includegraphics)
# Hardcoded prognostic ORs from 07b_prognostic.R
# =====================================================================
def efig2_interaction():
    print("eFigure 2: Surgery-type interaction (cardioplegia)")

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
        s = data[data.database == db]
        x = np.arange(len(s)) + offsets[db]
        ax.errorbar(
            x,
            s["or"],
            yerr=[s["or"] - s["lo"], s["hi"] - s["or"]],
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
    ax.set_ylabel(
        "OR per 1 mg/dL serum Mg increase\n" "(prognostic association with AKI)"
    )
    ax.set_ylim(0.7, 2.3)
    ax.set_yticks(np.arange(0.8, 2.4, 0.2))
    ax.legend(loc="upper left", fontsize=6)

    ax.annotate(
        "More cardioplegia \u2192\nhigher Mg + more AKI",
        xy=(0.88, 1.74),
        xytext=(0.40, 2.10),
        fontsize=5.5,
        color=C_GRAY,
        style="italic",
        arrowprops=dict(arrowstyle="->", color=C_GRAY, lw=0.5),
    )

    save(fig, "efig1_interaction")


# =====================================================================
# MAIN
# =====================================================================
if __name__ == "__main__":
    print("=" * 55)
    print("Generating publication figures")
    print("=" * 55)

    fig1_forest()
    fig_mg_stratified()
    fig3_subgroups()
    efig2_interaction()

    print("\nDone. Check figs/ directory.")
