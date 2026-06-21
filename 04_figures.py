#!/usr/bin/env python3
"""
04_figures.py — All publication figures (except CONSORT)

Usage:
  python 04_figures.py              # all figures
  python 04_figures.py primary      # Fig 1: primary outcome forest
  python 04_figures.py hte          # Fig 2: HTE forest
  python 04_figures.py benefit_harm # Fig 3: benefit-harm spectrum
  python 04_figures.py timecourse   # eFig: ΔCr time course
"""

import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ── Nature-style rcParams ─────────────────────────────────────────
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
    "legend.fontsize": 6,
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
    "savefig.pad_inches": 0.04,
}.items():
    mpl.rcParams[k] = v

# Wong/Okabe-Ito
BLUE = "#0072B2"
VERMIL = "#D55E00"
SKY = "#56B4E9"
ORANGE = "#E69F00"
GREEN = "#009E73"
GRAY = "#999999"

RESULTS = os.path.expanduser("~/mg_aki/results")
DB_ORDER = ["mimic", "eicu"]
DB_LABEL = {"mimic": "MIMIC-IV", "eicu": "eICU-CRD"}
DB_COLOR = {"mimic": BLUE, "eicu": VERMIL}
DB_MARKER = {"mimic": "o", "eicu": "s"}

SPEC = "primary"
POOL = "yet_untreated"
METHOD = "psm_dr"


def save(fig, name):
    for ext in ["pdf", "png"]:
        fig.savefig(os.path.join(RESULTS, f"{name}.{ext}"), format=ext, dpi=300)
    plt.close(fig)
    print(f"  ✓ {name}.pdf/.png")


def load_hte(tag):
    p = os.path.join(RESULTS, f"did_hte_{tag}.csv")
    return pd.read_csv(p) if os.path.exists(p) else None


def load_riskset(tag):
    p = os.path.join(RESULTS, f"did_riskset_{tag}.csv")
    return pd.read_csv(p) if os.path.exists(p) else None


def or_ci(row):
    """Compute OR CI from est (=OR) and se (=SE of log-OR)."""
    if pd.isna(row["est"]) or pd.isna(row["se"]) or row["est"] <= 0:
        return row["est"], np.nan, np.nan
    lo = np.exp(np.log(row["est"]) - 1.96 * row["se"])
    hi = np.exp(np.log(row["est"]) + 1.96 * row["se"])
    return row["est"], lo, hi


# ═══════════════════════════════════════════════════════════════════
# FIGURE 1: PRIMARY OUTCOME FOREST
# ═══════════════════════════════════════════════════════════════════
def fig_primary():
    """Forest plot: overall outcomes × 2 DBs + spec sensitivity for 48h AKI."""
    print("\n── Figure 1: Primary outcome forest ──")

    outcomes = [
        ("aki_48h", "48-hour AKI (KDIGO ≥ 1)"),
        ("aki_7d", "7-day AKI (ratio ≥ 1.5)"),
        ("hosp_mortality", "Hospital mortality"),
        ("vent_arrhythmia", "Ventricular arrhythmia"),
    ]

    fig, (ax_main, ax_spec) = plt.subplots(
        1,
        2,
        figsize=(7.2, 3.2),
        gridspec_kw={"width_ratios": [3, 2]},
        constrained_layout=True,
    )

    # ── Panel A: Overall outcomes ──
    ax = ax_main
    ax.axvline(1, color="#cccccc", lw=0.5, zorder=0)
    y_pos = []
    y_labels = []
    y = 0
    for oc_col, oc_label in reversed(outcomes):
        for tag in reversed(DB_ORDER):
            hte = load_hte(tag)
            if hte is None:
                continue
            row = hte[(hte.subgroup == "Overall") & (hte.outcome == oc_col)]
            if len(row) == 0:
                continue
            row = row.iloc[0]
            est, lo, hi = or_ci(row)
            if pd.isna(est):
                continue

            sig = not pd.isna(row["p"]) and row["p"] < 0.05
            ax.errorbar(
                est,
                y,
                xerr=[[est - lo], [hi - est]],
                fmt=DB_MARKER[tag],
                color=DB_COLOR[tag],
                ms=6 if sig else 4.5,
                markerfacecolor=DB_COLOR[tag] if sig else "white",
                markeredgecolor=DB_COLOR[tag],
                markeredgewidth=0.8,
                capsize=2,
                capthick=0.5,
                lw=0.7,
                zorder=3,
            )

            # Annotate OR
            txt = f"{est:.2f}"
            if not pd.isna(row["p"]):
                txt += "*" if row["p"] < 0.05 else ""
            ax.text(
                max(hi, est) + 0.04,
                y,
                txt,
                va="center",
                fontsize=5.5,
                color=DB_COLOR[tag],
            )
            y_pos.append(y)
            y_labels.append(f"{DB_LABEL[tag]}" if tag == DB_ORDER[0] else "")
            y += 1
        y += 0.5  # gap between outcomes

    ax.set_yticks([i * 2.5 + 0.5 for i in range(len(outcomes))])
    ax.set_yticklabels([oc[1] for oc in outcomes], fontsize=6.5)
    ax.set_xlabel("Odds Ratio (95% CI)")
    ax.set_xlim(0.3, 1.8)
    ax.text(
        -0.15,
        1.05,
        "a",
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="top",
    )
    ax.set_title("Overall outcomes", fontsize=7, pad=4)

    # Legend
    for tag in DB_ORDER:
        ax.plot(
            [],
            [],
            DB_MARKER[tag],
            color=DB_COLOR[tag],
            markerfacecolor=DB_COLOR[tag],
            ms=5,
            label=DB_LABEL[tag],
        )
    ax.legend(loc="upper right", fontsize=6, handletextpad=0.3)

    # ── Panel B: 48h AKI across specs ──
    ax = ax_spec
    ax.axvline(1, color="#cccccc", lw=0.5, zorder=0)

    spec_labels = {
        "primary": "Primary\n(19 var, no K⁺/Mg)",
        "sens_a": "Sensitivity A\n(21 var, + K⁺/Mg)",
        "sens_b": "Sensitivity B\n(19 var, FIRST labs)",
    }
    specs = ["primary", "sens_a", "sens_b"]
    y = 0
    for sn in reversed(specs):
        for tag in reversed(DB_ORDER):
            rs = load_riskset(tag)
            if rs is None:
                continue
            # Get 48h AKI from HTE for primary, but for sens_a/sens_b we need
            # to use ΔCr at 24h from riskset as a proxy... Actually we only
            # have HTE for primary spec. For spec comparison, use riskset DiD.
            row = rs[
                (rs.spec == sn)
                & (rs.pool == POOL)
                & (rs.method == METHOD)
                & (rs.target_h == 24)
            ]
            if len(row) == 0:
                continue
            row = row.iloc[0]
            if pd.isna(row["did"]):
                continue

            sig = not pd.isna(row["p"]) and row["p"] < 0.05
            ax.errorbar(
                row["did"],
                y,
                xerr=[[row["did"] - row["ci_lo"]], [row["ci_hi"] - row["did"]]],
                fmt=DB_MARKER[tag],
                color=DB_COLOR[tag],
                ms=6 if sig else 4.5,
                markerfacecolor=DB_COLOR[tag] if sig else "white",
                markeredgecolor=DB_COLOR[tag],
                markeredgewidth=0.8,
                capsize=2,
                capthick=0.5,
                lw=0.7,
                zorder=3,
            )
            y += 1
        y += 0.5

    ax.axvline(0, color="#cccccc", lw=0.5, zorder=0)
    ax.set_yticks([i * 2.5 + 0.5 for i in range(len(specs))])
    ax.set_yticklabels([spec_labels[s] for s in specs], fontsize=6)
    ax.set_xlabel("DiD: ΔCr at 24h (mg/dL)")
    ax.text(
        -0.18,
        1.05,
        "b",
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="top",
    )
    ax.set_title("Sensitivity: covariate specification", fontsize=7, pad=4)

    save(fig, "fig1_primary")


# ═══════════════════════════════════════════════════════════════════
# FIGURE 2: HTE FOREST
# ═══════════════════════════════════════════════════════════════════
def fig_hte():
    """HTE forest: subgroups × 48h AKI, both databases."""
    print("\n── Figure 2: HTE forest ──")

    subgroups = [
        "Overall",
        "Age < 65",
        "Age >= 65",
        "eGFR < 60",
        "eGFR >= 60",
        "Mg < 1.8",
        "Mg >= 1.8",
        "CABG",
        "Non-CABG",
        "Diabetes",
        "No diabetes",
        "CKD",
        "No CKD",
        "Heart failure",
        "No HF",
        "BMI >= 30",
        "BMI < 30",
        "DM + CKD",
        "HF + CABG",
        "Mg<1.8 + CKD",
    ]
    # Group separators (after these indices)
    group_gaps = {1, 3, 5, 7, 9, 11, 13, 15, 16}

    fig, ax = plt.subplots(figsize=(5.5, 7), constrained_layout=True)
    ax.axvline(1, color="#cccccc", lw=0.5, zorder=0)
    ax.axvspan(0.95, 1.05, color="#f5f5f5", zorder=0)

    y = 0
    yticks = []
    ylabels = []
    for i, sg in enumerate(reversed(subgroups)):
        for j, tag in enumerate(reversed(DB_ORDER)):
            hte = load_hte(tag)
            if hte is None:
                continue
            row = hte[(hte.subgroup == sg) & (hte.outcome == "aki_48h")]
            if len(row) == 0:
                continue
            row = row.iloc[0]
            est, lo, hi = or_ci(row)
            if pd.isna(est) or est > 10:
                continue

            sig = not pd.isna(row["p"]) and row["p"] < 0.05
            offset = 0.15 if tag == DB_ORDER[0] else -0.15
            ax.errorbar(
                est,
                y + offset,
                xerr=[[max(est - lo, 0)], [max(hi - est, 0)]],
                fmt=DB_MARKER[tag],
                color=DB_COLOR[tag],
                ms=5 if sig else 3.5,
                markerfacecolor=DB_COLOR[tag] if sig else "white",
                markeredgecolor=DB_COLOR[tag],
                markeredgewidth=0.7,
                capsize=1.5,
                capthick=0.4,
                lw=0.6,
                zorder=3,
            )

        yticks.append(y)
        ylabels.append(sg)
        idx_from_end = len(subgroups) - 1 - i
        if idx_from_end in group_gaps:
            y += 1.5
        else:
            y += 1

    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=6)
    ax.set_xlabel("Odds Ratio for 48h AKI (95% CI)")
    ax.set_xlim(0.05, 3.0)
    ax.set_xscale("log")
    ax.set_xticks([0.1, 0.25, 0.5, 1, 2])
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())

    # Legend
    for tag in DB_ORDER:
        ax.plot(
            [],
            [],
            DB_MARKER[tag],
            color=DB_COLOR[tag],
            markerfacecolor=DB_COLOR[tag],
            ms=5,
            label=DB_LABEL[tag],
        )
    ax.legend(loc="upper right", fontsize=6)

    # Annotations
    ax.text(
        0.08,
        -0.03,
        "← Favors IV Mg",
        transform=ax.transAxes,
        fontsize=5,
        color=GRAY,
        ha="left",
    )
    ax.text(
        0.92,
        -0.03,
        "Favors control →",
        transform=ax.transAxes,
        fontsize=5,
        color=GRAY,
        ha="right",
    )

    save(fig, "fig2_hte")


# ═══════════════════════════════════════════════════════════════════
# FIGURE 3: BENEFIT-HARM SPECTRUM
# ═══════════════════════════════════════════════════════════════════
def fig_benefit_harm():
    """Mg-stratified RD + crossed phenotype panel."""
    print("\n── Figure 3: Benefit-harm spectrum ──")

    fig, (ax_mg, ax_cross) = plt.subplots(
        1,
        2,
        figsize=(7.2, 3.5),
        constrained_layout=True,
    )

    # ── Panel A: Mg-stratified outcomes ──
    ax = ax_mg
    outcomes = [
        ("aki_48h", "48h AKI"),
        ("aki_7d", "7d AKI"),
        ("hosp_mortality", "Mortality"),
    ]
    mg_groups = ["Mg < 1.8", "Mg >= 1.8"]
    x = np.arange(len(outcomes))
    width = 0.18

    for gi, mg_sg in enumerate(mg_groups):
        for ti, tag in enumerate(DB_ORDER):
            hte = load_hte(tag)
            if hte is None:
                continue
            rds = []
            for oc_col, _ in outcomes:
                row = hte[(hte.subgroup == mg_sg) & (hte.outcome == oc_col)]
                if len(row) > 0 and not pd.isna(row.iloc[0]["rd"]):
                    rds.append(row.iloc[0]["rd"] * 100)
                else:
                    rds.append(0)

            offset = (gi * len(DB_ORDER) + ti - 1.5) * width
            color = DB_COLOR[tag]
            alpha = 1.0 if gi == 1 else 0.5  # Mg≥1.8 solid, <1.8 faded
            hatch = "" if gi == 1 else "///"
            bars = ax.bar(
                x + offset,
                rds,
                width * 0.9,
                color=color,
                alpha=alpha,
                hatch=hatch,
                edgecolor=color,
                linewidth=0.5,
            )

            for bar, rd in zip(bars, rds):
                if abs(rd) > 0.5:
                    ax.text(
                        bar.get_x() + bar.get_width() / 2,
                        rd,
                        f"{rd:+.1f}",
                        ha="center",
                        va="bottom" if rd > 0 else "top",
                        fontsize=4.5,
                        color=color,
                    )

    ax.axhline(0, color="#aaa", lw=0.4)
    ax.set_xticks(x)
    ax.set_xticklabels([oc[1] for oc in outcomes], fontsize=6.5)
    ax.set_ylabel("Risk Difference (%)")
    ax.text(
        -0.15,
        1.05,
        "a",
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="top",
    )
    ax.set_title("Mg-stratified treatment effect", fontsize=7, pad=4)

    # Manual legend
    from matplotlib.patches import Patch

    legend_elements = [
        Patch(facecolor=BLUE, alpha=1, label="MIMIC Mg≥1.8"),
        Patch(
            facecolor=BLUE, alpha=0.5, hatch="///", edgecolor=BLUE, label="MIMIC Mg<1.8"
        ),
        Patch(facecolor=VERMIL, alpha=1, label="eICU Mg≥1.8"),
        Patch(
            facecolor=VERMIL,
            alpha=0.5,
            hatch="///",
            edgecolor=VERMIL,
            label="eICU Mg<1.8",
        ),
    ]
    ax.legend(handles=legend_elements, fontsize=5, loc="lower left", ncol=2)

    # ── Panel B: Crossed phenotypes ──
    ax = ax_cross
    ax.axvline(1, color="#cccccc", lw=0.5, zorder=0)

    crossed = [
        ("HF + CABG", "HF + CABG"),
        ("DM + CKD", "DM + CKD"),
        ("Mg<1.8 + CKD", "Mg<1.8 + CKD (HARM?)"),
    ]

    y = 0
    yticks, ylabels = [], []
    for sg_key, sg_label in reversed(crossed):
        for tag in reversed(DB_ORDER):
            hte = load_hte(tag)
            if hte is None:
                continue
            row = hte[(hte.subgroup == sg_key) & (hte.outcome == "aki_48h")]
            if len(row) == 0:
                continue
            row = row.iloc[0]
            est, lo, hi = or_ci(row)
            if pd.isna(est) or est > 20:
                continue

            sig = not pd.isna(row["p"]) and row["p"] < 0.05
            offset = 0.12 if tag == DB_ORDER[0] else -0.12
            ax.errorbar(
                est,
                y + offset,
                xerr=[[max(est - lo, 0)], [max(hi - est, 0)]],
                fmt=DB_MARKER[tag],
                color=DB_COLOR[tag],
                ms=6 if sig else 4,
                markerfacecolor=DB_COLOR[tag] if sig else "white",
                markeredgecolor=DB_COLOR[tag],
                markeredgewidth=0.8,
                capsize=2,
                capthick=0.5,
                lw=0.7,
                zorder=3,
            )

            # n annotation
            n_txt = f"n={row['n_trt']}v{row['n_ctl']}"
            ax.text(
                0.02,
                y + offset,
                n_txt,
                fontsize=4.5,
                color=GRAY,
                transform=mpl.transforms.blended_transform_factory(
                    ax.transAxes, ax.transData
                ),
            )
        yticks.append(y)
        ylabels.append(sg_label)
        y += 1.5

    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=6.5)
    ax.set_xlabel("Odds Ratio for 48h AKI (95% CI)")
    ax.set_xscale("log")
    ax.set_xlim(0.05, 5)
    ax.set_xticks([0.1, 0.25, 0.5, 1, 2, 4])
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.text(
        -0.2,
        1.05,
        "b",
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="top",
    )
    ax.set_title("Crossed phenotypes", fontsize=7, pad=4)

    save(fig, "fig3_benefit_harm")


# ═══════════════════════════════════════════════════════════════════
# eFIGURE: TIME COURSE
# ═══════════════════════════════════════════════════════════════════
def fig_timecourse():
    """ΔCr time course (6-48h), primary spec, both databases."""
    print("\n── eFigure: ΔCr time course ──")

    dbs = []
    for tag in DB_ORDER:
        df = load_riskset(tag)
        if df is not None and "spec" in df.columns:
            dbs.append((tag, df))

    n = len(dbs)
    fig, axes = plt.subplots(
        1, n, figsize=(3.5 * n, 3.0), sharey=True, constrained_layout=True
    )
    if n == 1:
        axes = [axes]

    for i, (tag, df) in enumerate(dbs):
        ax = axes[i]
        ax.axhline(0, color="#aaa", lw=0.5, zorder=1)
        ax.axvspan(22, 26, color="#f0f0f0", zorder=0)

        sub = (
            df[(df.spec == SPEC) & (df.pool == POOL) & (df.method == METHOD)]
            .sort_values("target_h")
            .dropna(subset=["did"])
        )
        if len(sub) == 0:
            continue

        h = sub.target_h.values
        d = sub.did.values
        lo = sub.ci_lo.values
        hi = sub.ci_hi.values
        pv = sub.p.values

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
        )

        for j in range(len(h)):
            if not np.isnan(pv[j]) and pv[j] < 0.05:
                ax.plot(h[j], d[j], "o", ms=5, color=BLUE, zorder=5)

        ax.set_xticks([6, 12, 18, 24, 30, 36, 42, 48])
        ax.set_xlim(3, 51)
        ax.set_xlabel("Hours from T₀ (first IV Mg)")
        if i == 0:
            ax.set_ylabel("DiD: ΔCr (mg/dL)")

        ax.text(
            -0.12,
            1.06,
            chr(ord("a") + i),
            transform=ax.transAxes,
            fontsize=10,
            fontweight="bold",
            va="top",
        )
        ax.set_title(DB_LABEL.get(tag, tag), fontsize=8, pad=6)

    fig.text(
        0.5,
        -0.02,
        "Primary spec (19 var, no K⁺/Mg), PSM+DR  |  "
        "Filled ● = P<0.05  |  Gray = 24h  |  DiD<0 = renoprotective",
        ha="center",
        fontsize=5.5,
        color="#666",
    )

    save(fig, "efig_timecourse")


# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════
FIGURES = {
    "primary": fig_primary,
    "hte": fig_hte,
    "benefit_harm": fig_benefit_harm,
    "timecourse": fig_timecourse,
}

if __name__ == "__main__":
    print("=" * 70)
    print("04_figures.py — Publication figures")
    print("=" * 70)

    args = [a.lower() for a in sys.argv[1:]]
    to_draw = args if args else list(FIGURES.keys())

    for name in to_draw:
        if name in FIGURES:
            FIGURES[name]()
        else:
            print(f"  Unknown figure: {name}. Options: {list(FIGURES.keys())}")

    print("\n" + "=" * 70)
    print("Done.")
    print("=" * 70)
