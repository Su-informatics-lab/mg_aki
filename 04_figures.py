#!/usr/bin/env python3
"""
04_figures.py — Publication figures (aligned to manuscript blueprint)

  fig1_primary     Fig 1: ΔCr DiD 36h + AKI RD (the "does it work" figure)
  fig2_hte         Fig 2: HTE forest (the "who benefits" figure)
  fig3_benefit     Fig 3: Benefit-harm spectrum (the "who to treat" figure)
  efig_timecourse  eFig 3: ΔCr time course (PK plausibility)
  efig_sensitivity eFig 6: Primary vs Sens A (positivity demonstration)

Usage:
  python 04_figures.py                # all
  python 04_figures.py fig1_primary   # just one
"""

import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

# ── Nature rcParams ───────────────────────────────────────────────
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
    "savefig.pad_inches": 0.05,
}.items():
    mpl.rcParams[k] = v

BLUE = "#0072B2"
VERMIL = "#D55E00"
GRAY = "#999999"
RESULTS = os.path.expanduser("~/mg_aki/results")
DBS = ["mimic", "eicu"]
LBL = {"mimic": "MIMIC-IV", "eicu": "eICU-CRD"}
CLR = {"mimic": BLUE, "eicu": VERMIL}
MKR = {"mimic": "o", "eicu": "s"}


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(RESULTS, f"{name}.{ext}"), format=ext)
    plt.close(fig)
    print(f"  ✓ {name}")


def load_hte(tag):
    p = os.path.join(RESULTS, f"did_hte_{tag}.csv")
    return pd.read_csv(p) if os.path.exists(p) else None


def load_rs(tag):
    p = os.path.join(RESULTS, f"did_riskset_{tag}.csv")
    return pd.read_csv(p) if os.path.exists(p) else None


def hte_row(hte, sg, oc):
    """Get one row from HTE data."""
    r = hte[(hte.subgroup == sg) & (hte.outcome == oc)]
    return r.iloc[0] if len(r) > 0 else None


def or_ci(row):
    if row is None or pd.isna(row["est"]) or row["est"] <= 0 or pd.isna(row["se"]):
        return np.nan, np.nan, np.nan
    lo = np.exp(np.log(row["est"]) - 1.96 * row["se"])
    hi = np.exp(np.log(row["est"]) + 1.96 * row["se"])
    return row["est"], lo, hi


# ═══════════════════════════════════════════════════════════════════
# FIGURE 1: PRIMARY RESULT
#   Panel A: ΔCr DiD at 36h (bars)
#   Panel B: 48h AKI + 7d AKI risk difference (bars)
# ═══════════════════════════════════════════════════════════════════
def fig1_primary():
    print("\n── Fig 1: Primary result ──")
    fig, (ax_cr, ax_aki) = plt.subplots(
        1, 2, figsize=(6.5, 3.0), constrained_layout=True
    )
    width = 0.35

    # ── Panel A: ΔCr DiD at 36h ──
    ax = ax_cr
    ax.axhline(0, color="#bbb", lw=0.5)
    for i, tag in enumerate(DBS):
        rs = load_rs(tag)
        if rs is None:
            continue
        r = rs[
            (rs.spec == "primary")
            & (rs.pool == "yet_untreated")
            & (rs.method == "psm_dr")
            & (rs.target_h == 36)
        ]
        if len(r) == 0:
            continue
        r = r.iloc[0]
        x = i * (width + 0.1)
        sig = not pd.isna(r["p"]) and r["p"] < 0.05
        bar = ax.bar(
            x,
            r["did"],
            width,
            color=CLR[tag],
            alpha=0.85,
            edgecolor=CLR[tag],
            linewidth=0.5,
        )
        ax.errorbar(
            x,
            r["did"],
            yerr=[[r["did"] - r["ci_lo"]], [r["ci_hi"] - r["did"]]],
            fmt="none",
            color="black",
            capsize=3,
            capthick=0.6,
            lw=0.6,
        )
        stars = (
            "***"
            if r["p"] < 0.001
            else "**" if r["p"] < 0.01 else "*" if r["p"] < 0.05 else ""
        )
        ax.text(
            x,
            min(r["ci_lo"], r["did"]) - 0.003,
            f'{r["did"]:+.003f}{stars}',
            ha="center",
            va="top",
            fontsize=7,
            fontweight="bold",
            color=CLR[tag],
        )

    ax.set_xticks([i * (width + 0.1) for i in range(len(DBS))])
    ax.set_xticklabels([LBL[t] for t in DBS], fontsize=7)
    ax.set_ylabel("DiD: ΔCr at 36h (mg/dL)")
    ax.set_title("Creatinine change", fontsize=8, pad=4)
    ax.text(
        -0.18,
        1.06,
        "a",
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
    )

    # ── Panel B: AKI Risk Difference ──
    ax = ax_aki
    ax.axhline(0, color="#bbb", lw=0.5)
    outcomes = [
        ("aki_48h", "48h AKI\n(≥0.3 mg/dL)"),
        ("aki_7d", "7d AKI\n(ratio ≥1.5)"),
    ]

    for oi, (oc, oc_label) in enumerate(outcomes):
        for di, tag in enumerate(DBS):
            hte = load_hte(tag)
            if hte is None:
                continue
            r = hte_row(hte, "Overall", oc)
            if r is None or pd.isna(r["rd"]):
                continue
            rd_pct = r["rd"] * 100
            x = oi * 1.2 + di * (width + 0.05)
            sig = not pd.isna(r["p"]) and r["p"] < 0.05
            ax.bar(
                x,
                rd_pct,
                width,
                color=CLR[tag],
                alpha=0.85,
                edgecolor=CLR[tag],
                linewidth=0.5,
            )
            nnt = int(round(1 / abs(r["rd"]))) if abs(r["rd"]) > 0.005 else None
            label = f"{rd_pct:+.1f}%"
            if nnt:
                label += f"\nNNT={nnt}"
            ax.text(
                x,
                rd_pct - 0.3,
                label,
                ha="center",
                va="top",
                fontsize=5.5,
                color=CLR[tag],
                fontweight="bold",
            )

    ax.set_xticks([oi * 1.2 + 0.2 for oi in range(len(outcomes))])
    ax.set_xticklabels([oc[1] for oc in outcomes], fontsize=6.5)
    ax.set_ylabel("Risk Difference (%)")
    ax.set_title("Clinical AKI", fontsize=8, pad=4)
    ax.text(
        -0.18,
        1.06,
        "b",
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
    )

    # Shared legend
    from matplotlib.patches import Patch

    fig.legend(
        handles=[Patch(color=CLR[t], label=LBL[t]) for t in DBS],
        loc="lower center",
        ncol=2,
        fontsize=6.5,
        bbox_to_anchor=(0.5, -0.04),
    )

    save(fig, "fig1_primary")


# ═══════════════════════════════════════════════════════════════════
# FIGURE 2: HTE FOREST
# ═══════════════════════════════════════════════════════════════════
def fig2_hte():
    print("\n── Fig 2: HTE forest ──")

    # Ordered top to bottom as displayed
    groups = [
        ("Overall", None),
        ("Age < 65", "age_ge65"),
        ("Age >= 65", None),
        ("eGFR < 60", "egfr_lt60"),
        ("eGFR >= 60", None),
        ("Mg < 1.8", "mg_lt18"),
        ("Mg >= 1.8", None),
        ("CABG", "surg_cabg"),
        ("Non-CABG", None),
        ("Diabetes", "diabetes"),
        ("No diabetes", None),
        ("CKD", "ckd"),
        ("No CKD", None),
        ("Heart failure", "heart_failure"),
        ("No HF", None),
        ("BMI >= 30", "bmi_ge30"),
        ("BMI < 30", None),
        None,  # separator
        ("DM + CKD", None),
        ("HF + CABG", None),
        ("Mg<1.8 + CKD", None),
    ]

    # Load interaction P from the HTE data (MIMIC)
    # We'll just display them from known values
    htes = {tag: load_hte(tag) for tag in DBS}

    fig, ax = plt.subplots(figsize=(5.5, 7.5), constrained_layout=True)
    ax.axvline(1, color="#ddd", lw=0.6, zorder=0)

    y = 0
    yticks, ylabels = [], []
    prev_was_sep = False

    for item in reversed(groups):
        if item is None:
            y += 0.8
            prev_was_sep = True
            continue

        sg_name, int_var = item
        for di, tag in enumerate(reversed(DBS)):
            hte = htes.get(tag)
            if hte is None:
                continue
            r = hte_row(hte, sg_name, "aki_48h")
            if r is None:
                continue
            est, lo, hi = or_ci(r)
            if pd.isna(est) or est > 15:
                continue

            offset = 0.15 if di == 0 else -0.15
            sig = not pd.isna(r["p"]) and r["p"] < 0.05
            t = tag if di == 0 else DBS[0]  # reversed order
            actual_tag = DBS[1] if di == 0 else DBS[0]
            ax.errorbar(
                est,
                y + offset,
                xerr=[[max(est - lo, 0.001)], [max(hi - est, 0.001)]],
                fmt=MKR[actual_tag],
                color=CLR[actual_tag],
                ms=5 if sig else 3.5,
                markerfacecolor=CLR[actual_tag] if sig else "white",
                markeredgecolor=CLR[actual_tag],
                markeredgewidth=0.7,
                capsize=1.5,
                capthick=0.4,
                lw=0.5,
                zorder=3,
            )

        yticks.append(y)
        label = sg_name
        if sg_name == "Overall":
            label = "Overall"
        ylabels.append(label)

        # Gap between pairs
        idx = len(groups) - 1 - groups.index(item) if item in groups else 0
        y += 1.0

    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=6)
    ax.set_xlabel("Odds Ratio for 48h AKI (95% CI)")
    ax.set_xscale("log")
    ax.set_xlim(0.08, 4.0)
    ax.set_xticks([0.1, 0.25, 0.5, 1, 2])
    ax.xaxis.set_major_formatter(mticker.ScalarFormatter())

    for tag in DBS:
        ax.plot(
            [],
            [],
            MKR[tag],
            color=CLR[tag],
            markerfacecolor=CLR[tag],
            ms=5,
            label=LBL[tag],
        )
    ax.legend(loc="upper right", fontsize=6.5)

    ax.text(
        0.02, -0.025, "← Favors IV Mg", transform=ax.transAxes, fontsize=5.5, color=GRAY
    )
    ax.text(
        0.98,
        -0.025,
        "Favors control →",
        transform=ax.transAxes,
        fontsize=5.5,
        color=GRAY,
        ha="right",
    )

    save(fig, "fig2_hte")


# ═══════════════════════════════════════════════════════════════════
# FIGURE 3: BENEFIT-HARM SPECTRUM
# ═══════════════════════════════════════════════════════════════════
def fig3_benefit():
    print("\n── Fig 3: Benefit-harm spectrum ──")
    fig, (ax_mg, ax_cross) = plt.subplots(
        1,
        2,
        figsize=(7.0, 3.2),
        gridspec_kw={"width_ratios": [3, 2]},
        constrained_layout=True,
    )
    width = 0.28

    # ── Panel A: Mg-stratified RD% for 48h AKI ──
    ax = ax_mg
    ax.axhline(0, color="#bbb", lw=0.5)
    mg_groups = ["Mg < 1.8", "Mg >= 1.8"]
    oc_list = [
        ("aki_48h", "48h AKI"),
        ("aki_7d", "7d AKI"),
        ("hosp_mortality", "Mortality"),
    ]

    for gi, mg_sg in enumerate(mg_groups):
        for oi, (oc, oc_lbl) in enumerate(oc_list):
            for di, tag in enumerate(DBS):
                hte = load_hte(tag)
                if hte is None:
                    continue
                r = hte_row(hte, mg_sg, oc)
                if r is None or pd.isna(r["rd"]):
                    continue
                rd_pct = r["rd"] * 100
                x = oi * 1.6 + (gi * len(DBS) + di) * (width + 0.02)
                alpha = 0.9 if gi == 1 else 0.4
                hatch = "" if gi == 1 else "///"
                ax.bar(
                    x,
                    rd_pct,
                    width,
                    color=CLR[tag],
                    alpha=alpha,
                    hatch=hatch,
                    edgecolor=CLR[tag],
                    linewidth=0.4,
                )
                if abs(rd_pct) > 0.8:
                    va = "bottom" if rd_pct > 0 else "top"
                    ax.text(
                        x,
                        rd_pct,
                        f"{rd_pct:+.1f}",
                        ha="center",
                        va=va,
                        fontsize=4.5,
                        color=CLR[tag],
                    )

    ax.set_xticks([oi * 1.6 + 0.6 for oi in range(len(oc_list))])
    ax.set_xticklabels([oc[1] for oc in oc_list], fontsize=6.5)
    ax.set_ylabel("Risk Difference (%)")
    ax.text(
        -0.14,
        1.06,
        "a",
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
    )
    ax.set_title("Mg-stratified treatment effect", fontsize=7, pad=4)

    from matplotlib.patches import Patch

    ax.legend(
        handles=[
            Patch(color=BLUE, alpha=0.9, label="MIMIC Mg≥1.8"),
            Patch(
                color=BLUE, alpha=0.4, hatch="///", edgecolor=BLUE, label="MIMIC Mg<1.8"
            ),
            Patch(color=VERMIL, alpha=0.9, label="eICU Mg≥1.8"),
            Patch(
                color=VERMIL,
                alpha=0.4,
                hatch="///",
                edgecolor=VERMIL,
                label="eICU Mg<1.8",
            ),
        ],
        fontsize=5,
        loc="lower left",
        ncol=2,
    )

    # ── Panel B: Crossed phenotype forest ──
    ax = ax_cross
    ax.axvline(1, color="#ddd", lw=0.6, zorder=0)

    crossed = [
        ("HF + CABG", "HF + CABG"),
        ("DM + CKD", "DM + CKD"),
        ("Mg<1.8 + CKD", "Mg<1.8 + CKD"),
    ]

    y = 0
    yticks, ylabels = [], []
    for sg_key, sg_label in reversed(crossed):
        for di, tag in enumerate(reversed(DBS)):
            hte = load_hte(tag)
            if hte is None:
                continue
            r = hte_row(hte, sg_key, "aki_48h")
            if r is None:
                continue
            est, lo, hi = or_ci(r)
            if pd.isna(est) or est > 20:
                continue

            actual_tag = DBS[1 - di]
            offset = 0.12 if di == 0 else -0.12
            sig = not pd.isna(r["p"]) and r["p"] < 0.05
            ax.errorbar(
                est,
                y + offset,
                xerr=[[max(est - lo, 0.001)], [max(hi - est, 0.001)]],
                fmt=MKR[actual_tag],
                color=CLR[actual_tag],
                ms=6 if sig else 4,
                markerfacecolor=CLR[actual_tag] if sig else "white",
                markeredgecolor=CLR[actual_tag],
                markeredgewidth=0.8,
                capsize=2,
                capthick=0.5,
                lw=0.6,
                zorder=3,
            )

            # n annotation
            ax.text(
                0.03,
                y + offset,
                f'n={int(r["n_trt"])}v{int(r["n_ctl"])}',
                fontsize=4.5,
                color=GRAY,
                va="center",
                transform=mpl.transforms.blended_transform_factory(
                    ax.transAxes, ax.transData
                ),
            )

        yticks.append(y)
        ylabels.append(sg_label)
        y += 1.3

    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=6.5)
    ax.set_xlabel("OR for 48h AKI (95% CI)")
    ax.set_xscale("log")
    ax.set_xlim(0.05, 6)
    ax.set_xticks([0.1, 0.25, 0.5, 1, 2, 4])
    ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
    ax.text(
        -0.22,
        1.06,
        "b",
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
    )
    ax.set_title("Crossed phenotypes", fontsize=7, pad=4)

    save(fig, "fig3_benefit_harm")


# ═══════════════════════════════════════════════════════════════════
# eFIGURE: TIME COURSE
# ═══════════════════════════════════════════════════════════════════
def efig_timecourse():
    print("\n── eFig: ΔCr time course ──")
    avail = []
    for tag in DBS:
        df = load_rs(tag)
        if df is not None and "spec" in df.columns:
            avail.append((tag, df))

    n = len(avail)
    fig, axes = plt.subplots(
        1, n, figsize=(3.5 * n, 2.8), sharey=True, constrained_layout=True
    )
    if n == 1:
        axes = [axes]

    for i, (tag, df) in enumerate(avail):
        ax = axes[i]
        ax.axhline(0, color="#aaa", lw=0.5)
        ax.axvspan(34, 38, color="#f0f0f0", zorder=0)

        sub = (
            df[
                (df.spec == "primary")
                & (df.pool == "yet_untreated")
                & (df.method == "psm_dr")
            ]
            .sort_values("target_h")
            .dropna(subset=["did"])
        )
        if len(sub) == 0:
            continue

        h, d = sub.target_h.values, sub.did.values
        lo, hi, pv = sub.ci_lo.values, sub.ci_hi.values, sub.p.values

        ax.fill_between(h, lo, hi, alpha=0.15, color=BLUE, zorder=2)
        ax.plot(
            h,
            d,
            color=BLUE,
            lw=1.4,
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
        ax.set_xlabel("Hours from T₀")
        if i == 0:
            ax.set_ylabel("DiD: ΔCr (mg/dL)")
        ax.text(
            -0.14,
            1.06,
            chr(ord("a") + i),
            transform=ax.transAxes,
            fontsize=10,
            fontweight="bold",
            va="top",
        )
        ax.set_title(LBL.get(tag, tag), fontsize=8, pad=5)

    fig.text(
        0.5,
        -0.03,
        "Primary (19 var, no K⁺/Mg), PSM+DR  |  Filled ● = P<0.05  |  "
        "Gray band = 36h primary  |  DiD<0 = renoprotective",
        ha="center",
        fontsize=5.5,
        color="#666",
    )
    save(fig, "efig_timecourse")


# ═══════════════════════════════════════════════════════════════════
# eFIGURE: SENSITIVITY (primary vs sens_a)
# ═══════════════════════════════════════════════════════════════════
def efig_sensitivity():
    print("\n── eFig: Primary vs Sensitivity A ──")
    avail = []
    for tag in DBS:
        df = load_rs(tag)
        if df is not None and "spec" in df.columns:
            avail.append((tag, df))

    specs = [
        ("primary", "Primary (no K⁺/Mg)", BLUE),
        ("sens_a", "Sensitivity A (+K⁺/Mg)", VERMIL),
    ]
    n = len(avail)
    fig, axes = plt.subplots(
        1, n, figsize=(3.5 * n, 2.8), sharey=True, constrained_layout=True
    )
    if n == 1:
        axes = [axes]

    for i, (tag, df) in enumerate(avail):
        ax = axes[i]
        ax.axhline(0, color="#aaa", lw=0.5)
        for sn, slbl, scol in specs:
            sub = (
                df[
                    (df.spec == sn)
                    & (df.pool == "yet_untreated")
                    & (df.method == "psm_dr")
                ]
                .sort_values("target_h")
                .dropna(subset=["did"])
            )
            if len(sub) == 0:
                continue
            h, d = sub.target_h.values, sub.did.values
            lo, hi, pv = sub.ci_lo.values, sub.ci_hi.values, sub.p.values
            ls = "-" if sn == "primary" else "--"
            ax.fill_between(h, lo, hi, alpha=0.10, color=scol)
            ax.plot(
                h,
                d,
                color=scol,
                lw=1.2,
                ls=ls,
                marker="o",
                ms=4,
                markerfacecolor="white",
                markeredgecolor=scol,
                markeredgewidth=0.8,
                label=slbl,
            )
            for j in range(len(h)):
                if not np.isnan(pv[j]) and pv[j] < 0.05:
                    ax.plot(h[j], d[j], "o", ms=4, color=scol, zorder=5)

        ax.set_xticks([6, 12, 18, 24, 30, 36, 42, 48])
        ax.set_xlim(3, 51)
        ax.set_xlabel("Hours from T₀")
        if i == 0:
            ax.set_ylabel("DiD: ΔCr (mg/dL)")
            ax.legend(fontsize=5.5, loc="lower left")
        ax.text(
            -0.14,
            1.06,
            chr(ord("a") + i),
            transform=ax.transAxes,
            fontsize=10,
            fontweight="bold",
            va="top",
        )
        ax.set_title(LBL.get(tag, tag), fontsize=8, pad=5)

    save(fig, "efig_sensitivity")


# ═══════════════════════════════════════════════════════════════════
FIGURES = {
    "fig1_primary": fig1_primary,
    "fig2_hte": fig2_hte,
    "fig3_benefit": fig3_benefit,
    "efig_timecourse": efig_timecourse,
    "efig_sensitivity": efig_sensitivity,
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
            print(f"  Unknown: {name}. Options: {list(FIGURES.keys())}")
    print("=" * 70)
