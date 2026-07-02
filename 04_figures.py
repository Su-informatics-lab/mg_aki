#!/usr/bin/env python3
"""
04_figures.py — Publication figures for Critical Care submission

Main text:
  fig2_primary_forest   — Overall binary outcomes (KDIGO cascade + mortality)
  fig3_egfr_forest      — eGFR-stratified AKI stages + mortality (THE finding)
  fig4_timecourse       — ΔCr time course (PK plausibility)
  fig5_egfr_mg_heatmap  — eGFR × Mg cross-stratification, 48-h AKI (confounding defense)

Supplement:
  efig1_love            — Love plots (primary + earliest-labs sensitivity)
  efig2_sensitivity     — Primary vs earliest-labs sensitivity DiD
  efig4_hte_forest      — Pre-specified subgroup forest
  efig5_crossed_forest  — Crossed phenotype forest
  egfr_mg_heatmap_7d    — eGFR × Mg cross-stratification, 7-d AKI

Usage:
  python 04_figures.py                    # all
  python 04_figures.py fig3_egfr_forest   # single figure
"""

import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

# ═══════════════════════════════════════════════════════════════════
# NATURE / CC STYLE (nature-visualizer skill)
# ═══════════════════════════════════════════════════════════════════
mpl.rcParams.update(
    {
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "legend.title_fontsize": 7,
        "axes.linewidth": 0.5,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "lines.linewidth": 1.0,
        "lines.markersize": 4,
        "legend.frameon": False,
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "figure.constrained_layout.use": True,
    }
)

# Wong/Okabe-Ito palette (Nature Methods 2011)
BLUE = "#0072B2"  # MIMIC
VERMIL = "#D55E00"  # eICU
GREEN = "#009E73"
ORANGE = "#E69F00"
SKY = "#56B4E9"
PURPLE = "#CC79A7"
BLACK = "#000000"
GRAY = "#999999"

# Sizing
W_SINGLE = 3.504  # 89mm
W_HALF = 4.724  # 120mm
W_DOUBLE = 7.205  # 183mm

RESULTS = os.path.expanduser("~/mg_aki/results")
FIG_DIR = os.path.join(RESULTS, "figures")
os.makedirs(FIG_DIR, exist_ok=True)
DBS = ["mimic", "eicu"]
LBL = {"mimic": "MIMIC-IV", "eicu": "eICU-CRD"}
CLR = {"mimic": BLUE, "eicu": VERMIL}
MKR = {"mimic": "o", "eicu": "s"}


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG_DIR, f"{name}.{ext}"), format=ext)
    plt.close(fig)
    print(f"  ✓ {name}")


def panel_label(ax, label, x=-0.12, y=1.06):
    """Nature-style lowercase bold panel label."""
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


# ═══════════════════════════════════════════════════════════════════
# FIG 2: PRIMARY OUTCOMES FOREST
# Two panels (MIMIC | eICU), rows = binary outcomes
# ═══════════════════════════════════════════════════════════════════
def fig2_primary_forest():
    """Overall binary outcomes forest — reads did_binary_{db}.csv.
    AKI outcomes only (mortality removed per editorial decision)."""
    fig, axes = plt.subplots(1, 2, figsize=(W_DOUBLE, W_DOUBLE * 0.30), sharey=True)
    # Row order: AKI outcomes only (mortality kept in Table 2 only)
    row_order = [
        ("aki1_48h", "48-h AKI (KDIGO ≥1)"),
        ("aki1_7d", "7-d AKI (KDIGO ≥1)"),
        ("aki2_7d", "7-d AKI (KDIGO ≥2)"),
        ("aki3_7d", "7-d AKI (KDIGO ≥3)"),
    ]
    for i, (db, ax) in enumerate(zip(DBS, axes)):
        panel_label(ax, chr(ord("a") + i))
        df = pd.read_csv(os.path.join(RESULTS, f"did_binary_{db}.csv"))
        df = df.set_index("outcome")
        y_pos = list(range(len(row_order)))
        for j, (oname, label) in enumerate(row_order):
            if oname not in df.index:
                continue
            r = df.loc[oname]
            if pd.isna(r["or"]):
                ax.text(1.0, j, "n/a", ha="center", va="center", fontsize=5, color=GRAY)
                continue
            or_val, lo, hi = r["or"], r["ci_lo"], r["ci_hi"]
            color = GREEN if or_val < 1 else VERMIL
            ax.plot(or_val, j, MKR[db], color=CLR[db], markersize=5, zorder=3)
            ax.plot([lo, hi], [j, j], color=CLR[db], linewidth=1.0, zorder=2)
            # Annotate OR (CI)
            txt = f"{or_val:.2f} ({lo:.2f}–{hi:.2f})"
            sig = "  *" if not pd.isna(r["p"]) and r["p"] < 0.05 else ""
            ax.annotate(txt + sig, (hi + 0.02, j), fontsize=5, va="center", ha="left")

        ax.axvline(1.0, color=BLACK, linewidth=0.5, linestyle="--", zorder=1)
        ax.set_yticks(y_pos)
        ax.set_yticklabels([lbl for _, lbl in row_order], fontsize=6)
        ax.set_xlabel("Odds ratio (95% CI)")
        ax.set_title(LBL[db], fontweight="bold", fontsize=7)
        ax.set_xlim(0.4, 2.0)
        ax.invert_yaxis()

    save(fig, "primary_forest")


# ═══════════════════════════════════════════════════════════════════
# FIG 3: eGFR-STRATIFIED FOREST (THE finding)
# 3 panels: AKI Stage 1+ | Stage 2+/3+ | Mortality
# Each panel: 5 eGFR strata + Overall, both DBs overlaid
# ═══════════════════════════════════════════════════════════════════
def fig3_egfr_forest():
    """eGFR-stratified AKI stages + mortality (2x3 layout).
    Row 1: 48h primary (AKI >=1, >=2) + hospital mortality
    Row 2: 7d secondary (AKI >=1, >=2, >=3)
    """
    panels_top = [
        ("AKI_48h_Stage1+", "48-h AKI (KDIGO \u22651)"),
        ("AKI_48h_Stage2+", "48-h AKI (KDIGO \u22652)"),
        ("hosp_mortality", "Hospital mortality"),
    ]
    panels_bot = [
        ("AKI_Stage1+", "7-d AKI (KDIGO \u22651)"),
        ("AKI_Stage2+", "7-d AKI (KDIGO \u22652)"),
        ("AKI_Stage3+", "7-d AKI (KDIGO \u22653)"),
    ]
    all_panels = panels_top + panels_bot
    # Grid positions: 2x3, all 6 slots filled (mortality at top-right)
    grid_pos = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]
    strata_order = [
        "Overall",
        "eGFR>=90",
        "eGFR_60-89",
        "eGFR_45-59",
        "eGFR_30-44",
        "eGFR<30",
    ]
    strata_labels = [
        "Overall",
        "eGFR ≥90",
        "eGFR 60–89",
        "eGFR 45–59",
        "eGFR 30–44",
        "eGFR <30",
    ]
    n_strata = len(strata_order)

    fig, axes = plt.subplots(2, 3, figsize=(W_DOUBLE, W_DOUBLE * 0.75), sharey=True)

    for pi, (outcome, title) in enumerate(all_panels):
        row, col = grid_pos[pi]
        ax = axes[row, col]
        panel_label(ax, chr(ord("a") + pi))
        ax.set_title(title, fontsize=7)

        for di, db in enumerate(DBS):
            df = pd.read_csv(os.path.join(RESULTS, f"egfr_aki_stages_{db}.csv"))
            sub = df[df.outcome == outcome].copy()
            sub["stratum"] = pd.Categorical(
                sub.stratum, categories=strata_order, ordered=True
            )
            sub = sub.sort_values("stratum")

            offset = -0.12 if di == 0 else 0.12  # dodge
            for j, strat in enumerate(strata_order):
                row = sub[sub.stratum == strat]
                if len(row) == 0 or pd.isna(row.iloc[0]["or"]):
                    continue
                r = row.iloc[0]
                y = j + offset
                ax.plot(r["or"], y, MKR[db], color=CLR[db], markersize=4, zorder=3)
                ax.plot(
                    [r.or_lo, r.or_hi], [y, y], color=CLR[db], linewidth=0.8, zorder=2
                )

        # Reference line + styling
        ax.axvline(1.0, color=BLACK, linewidth=0.5, linestyle="--", zorder=1)
        ax.set_yticks(range(n_strata))
        if col == 0:
            ax.set_yticklabels(strata_labels, fontsize=6)
        ax.set_xlabel("OR (95% CI)", fontsize=6)
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: f"{x:.1f}" if x >= 1 else f"{x:.2f}")
        )
        ax.set_xlim(0.1, 20) if outcome == "hosp_mortality" else ax.set_xlim(0.15, 12)
        ax.invert_yaxis()

        # Shading for overall row
        ax.axhspan(-0.5, 0.5, color="#f0f0f0", zorder=0)

    # Shared legend
    from matplotlib.lines import Line2D

    handles = [
        Line2D(
            [0],
            [0],
            marker=MKR[db],
            color=CLR[db],
            linestyle="",
            markersize=4,
            label=LBL[db],
        )
        for db in DBS
    ]
    axes[1, 2].legend(handles=handles, loc="lower right", fontsize=6)

    # Interaction P annotation
    for db in DBS:
        df = pd.read_csv(os.path.join(RESULTS, f"did_hte_interact_{db}.csv"))
        egfr_p = df[
            (df.variable.str.contains("eGFR|egfr", case=False))
            & (df.outcome.str.contains("aki_7d", case=False))
        ]
        if len(egfr_p) > 0:
            p = egfr_p.iloc[0].p_interaction
            print(f"    {db} eGFR interaction P (AKI 7d): {p:.6f}")

    save(fig, "egfr_forest")


# ═══════════════════════════════════════════════════════════════════
# FIG 4: ΔCr TIME COURSE (PK plausibility)
# Both DBs, CI bands, peak annotation
# ═══════════════════════════════════════════════════════════════════
def fig4_timecourse():
    """ΔCr DiD time course — reads did_riskset_{db}.csv."""
    fig, ax = plt.subplots(figsize=(W_HALF, W_HALF * 0.65))

    for db in DBS:
        df = pd.read_csv(os.path.join(RESULTS, f"did_riskset_{db}.csv"))
        sub = df[
            (df.spec == "primary")
            & (df.pool == "yet_untreated")
            & (df.method == "psm_dr")
        ].sort_values("target_h")
        if len(sub) == 0:
            continue
        x = np.concatenate([[0], sub.target_h.values])
        y = np.concatenate([[0], sub.did.values])
        lo = np.concatenate([[0], sub.ci_lo.values])
        hi = np.concatenate([[0], sub.ci_hi.values])

        ax.plot(
            x,
            y,
            f"-{MKR[db]}",
            color=CLR[db],
            label=LBL[db],
            markersize=4,
            linewidth=1.0,
        )
        ax.fill_between(x, lo, hi, color=CLR[db], alpha=0.12)

        # Mark significant (skip origin point)
        sig = np.concatenate([[False], (sub.p < 0.05).values])
        ax.scatter(
            x[sig],
            y[sig],
            marker=MKR[db],
            color=CLR[db],
            s=16,
            zorder=5,
            edgecolors="white",
            linewidths=0.3,
        )

    ax.axhline(0, color=BLACK, linewidth=0.5, linestyle="--")
    ax.set_xlabel("Hours from T₀")
    ax.set_ylabel("ΔCr DiD (mg/dL)")
    ax.set_xticks([0, 6, 12, 18, 24, 30, 36, 42, 48])
    ax.set_xlim(0, 50)
    ax.legend(loc="lower left")

    save(fig, "timecourse")


# ═══════════════════════════════════════════════════════════════════
# FIG 5: eGFR × Mg CROSS-STRATIFICATION HEATMAP
# Shows eGFR reversal persists within each Mg stratum
# ═══════════════════════════════════════════════════════════════════
def fig5_egfr_mg_heatmap(aki_outcome="aki_48h", name="egfr_mg_heatmap"):
    """eGFR × Mg OR heatmap — reads mg_strat_{db}.csv.

    aki_outcome: AKI window (default 48h primary; "aki_7d" for suppl).
    name: output filename stem.
    """
    fig, axes = plt.subplots(1, 2, figsize=(W_DOUBLE, W_DOUBLE * 0.32))

    for col_i, db in enumerate(DBS):
        df = pd.read_csv(os.path.join(RESULTS, f"mg_strat_{db}.csv"))

        ax = axes[col_i]
        panel_label(ax, chr(ord("a") + col_i))

        sub = df[(df.outcome == aki_outcome) & (df.egfr_strat != "All")].copy()
        if len(sub) == 0:
            ax.text(
                0.5,
                0.5,
                "no cross-strat data",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=6,
            )
            continue

        mg_bins = [
            b
            for b in [
                "Mg<1.6",
                "Mg_1.6-2.0",
                "Mg>=2.0",
                "Mg<1.8",
                "Mg>=1.8",
                "Mg_1.6-1.8",
                "Mg_1.8-2.0",
                "Mg_2.0-2.3",
                "Mg>=2.3",
            ]
            if b in sub.mg_strat.values
        ]
        egfr_bins = [
            b
            for b in ["eGFR>=90", "eGFR_60-89", "eGFR_45-59", "eGFR<45"]
            if b in sub.egfr_strat.values
        ]

        if not mg_bins or not egfr_bins:
            ax.text(
                0.5,
                0.5,
                "insufficient strata",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=6,
            )
            continue

        mat = np.full((len(mg_bins), len(egfr_bins)), np.nan)
        pmat = np.full_like(mat, np.nan)
        nmat = np.full((len(mg_bins), len(egfr_bins)), 0)
        rmat_t = np.full_like(mat, np.nan)
        rmat_c = np.full_like(mat, np.nan)
        for mi, mg in enumerate(mg_bins):
            for ei, eg in enumerate(egfr_bins):
                row = sub[(sub.mg_strat == mg) & (sub.egfr_strat == eg)]
                if len(row) > 0:
                    r = row.iloc[0]
                    mat[mi, ei] = r["or"]
                    pmat[mi, ei] = r["p"]
                    nmat[mi, ei] = int(r["n"])
                    rmat_t[mi, ei] = r["rate_trt"]
                    rmat_c[mi, ei] = r["rate_ctl"]

        import matplotlib.colors as mcolors
        from matplotlib.patches import Rectangle

        # Log-symmetric norm: OR 0.5 and OR 2.0 get equal visual weight
        mat_clipped = np.clip(mat, 1 / 3.0, 3.0)
        norm = mcolors.LogNorm(vmin=1 / 3.0, vmax=3.0)
        im = ax.imshow(mat_clipped, cmap="RdBu_r", norm=norm, aspect="auto")
        for mi in range(len(mg_bins)):
            for ei in range(len(egfr_bins)):
                v = mat[mi, ei]
                p = pmat[mi, ei]
                n = nmat[mi, ei]
                rt = rmat_t[mi, ei]
                rc = rmat_c[mi, ei]
                if pd.isna(v):
                    continue
                events_trt = n * rt if not pd.isna(rt) else 0
                events_ctl = n * rc if not pd.isna(rc) else 0
                if min(events_trt, events_ctl) < 20:
                    rect = Rectangle(
                        (ei - 0.5, mi - 0.5),
                        1,
                        1,
                        fill=False,
                        hatch="///",
                        linewidth=0,
                        edgecolor="gray",
                        alpha=0.5,
                    )
                    ax.add_patch(rect)
                sig = "*" if not pd.isna(p) and p < 0.05 else ""
                txt_color = "white" if (v < 0.4 or v > 2.5) else "black"
                ax.text(
                    ei,
                    mi,
                    f"{v:.2f}{sig}\nn={n}",
                    ha="center",
                    va="center",
                    fontsize=4.5,
                    color=txt_color,
                    fontweight="bold",
                )

        ax.set_xticks(range(len(egfr_bins)))
        ax.set_xticklabels(
            [b.replace("eGFR", "").replace("_", " ") for b in egfr_bins], fontsize=5
        )
        ax.set_yticks(range(len(mg_bins)))
        ax.set_yticklabels([b.replace("Mg_", "Mg ") for b in mg_bins], fontsize=5)
        ax.set_xlabel("eGFR stratum", fontsize=6)
        if col_i == 0:
            ax.set_ylabel("Baseline Mg", fontsize=6)
        metric = {"aki_48h": "48-h AKI", "aki_7d": "7-d AKI"}.get(
            aki_outcome, aki_outcome
        )
        ax.set_title(f"{LBL[db]} \u2014 {metric}", fontsize=6, fontweight="bold")

    cbar = fig.colorbar(im, ax=axes, shrink=0.6, pad=0.02)
    cbar.set_ticks([1 / 3.0, 0.5, 1.0, 2.0, 3.0])
    cbar.set_ticklabels(["0.33", "0.5", "1.0", "2.0", "3.0"])
    cbar.set_label("OR", fontsize=6)
    cbar.ax.tick_params(labelsize=5)

    save(fig, name)


# ═══════════════════════════════════════════════════════════════════
# eFIG 1: LOVE PLOTS (primary + earliest-labs sensitivity)
# ═══════════════════════════════════════════════════════════════════
def efig1_love():
    """Love plots: primary and earliest-labs sensitivity, separate panels."""
    # Per-spec variable lists matching 02_psm.R SPECS
    PS_BASE = [
        "age",
        "is_female",
        "bmi",
        "surg_cabg",
        "surg_valve",
        "surg_combined",
        "heart_failure",
        "hypertension",
        "diabetes",
        "ckd",
        "copd",
        "pvd",
        "stroke",
        "liver_disease",
        "egfr",
    ]
    SPEC_VARS = {
        "primary_yet_untreated": PS_BASE
        + ["last_calcium", "last_lactate", "last_lactate_missing", "last_heartrate"],
        "sens_b_yet_untreated": PS_BASE
        + [
            "first_calcium",
            "first_lactate",
            "first_lactate_missing",
            "first_heartrate",
        ],
    }
    specs = [
        ("primary_yet_untreated", "Primary (19 covariates)"),
        ("sens_b_yet_untreated", "Sensitivity (earliest labs)"),
    ]
    db = "mimic"
    all_pts = pd.read_csv(os.path.join(RESULTS, f"did_all_{db}.csv"))
    all_pts["_pid_str"] = all_pts.pid.astype(str)

    # Compute time-varying labs from did_labs_all
    labs_path = os.path.join(RESULTS, f"did_labs_all_{db}.csv")
    pid_col_l = "stay_id"
    if os.path.exists(labs_path):
        labs_raw = pd.read_csv(labs_path)
        if "patientunitstayid" in labs_raw.columns:
            pid_col_l = "patientunitstayid"
        labs_raw["pid"] = labs_raw[pid_col_l].astype(str)
        lab_names = ["magnesium", "potassium", "calcium", "lactate", "heartrate"]
        # For each patient, compute last and first values before T0
        mg_map = all_pts.set_index("_pid_str")["mg_offset_h"].to_dict()
        for prefix, ascending in [("last", False), ("first", True)]:
            for ln in lab_names:
                sub = labs_raw[labs_raw.lab_name == ln].copy()
                if len(sub) == 0:
                    continue
                sub["mg_h"] = sub.pid.map(mg_map)
                # Keep only measurements before T0 (or all for controls)
                sub = sub[
                    (sub.offset_h >= 0) & (sub.mg_h.isna() | (sub.offset_h < sub.mg_h))
                ]
                if ascending:
                    sub = sub.sort_values("offset_h").groupby("pid").first()
                else:
                    sub = (
                        sub.sort_values("offset_h", ascending=False)
                        .groupby("pid")
                        .first()
                    )
                col = f"{prefix}_{ln}"
                val_map = sub["value"].to_dict()
                all_pts[col] = all_pts._pid_str.map(val_map)
        # Lactate missing indicator
        for prefix in ["last", "first"]:
            col = f"{prefix}_lactate"
            miss_col = f"{prefix}_lactate_missing"
            all_pts[miss_col] = all_pts[col].isna().astype(float)
        print("    Time-varying labs computed from did_labs_all")

    fig, axes = plt.subplots(1, 2, figsize=(W_HALF, W_HALF * 0.85), sharey=False)

    def compute_smd(x1, x0):
        m1, m0 = x1.mean(), x0.mean()
        sp = np.sqrt((x1.var() + x0.var()) / 2)
        return abs(m1 - m0) / sp if sp > 1e-10 else 0

    raw_trt = all_pts[all_pts.treated == 1]
    raw_ctl = all_pts[all_pts.treated == 0]

    for si, (spec_tag, spec_label) in enumerate(specs):
        ax = axes[si]
        panel_label(ax, chr(ord("a") + si))
        pairs_path = os.path.join(RESULTS, f"did_pairs_{spec_tag}_{db}.csv")
        if not os.path.exists(pairs_path):
            ax.set_title(f"{spec_label}\n(pairs not found)", fontsize=6)
            continue

        pairs = pd.read_csv(pairs_path)
        trt_pids = set(pairs.trt_pid.astype(str))
        ctl_pids = set(pairs.ctl_pid.astype(str))
        trt = all_pts[all_pts._pid_str.isin(trt_pids)]
        ctl = all_pts[all_pts._pid_str.isin(ctl_pids)]

        ps_vars = SPEC_VARS.get(spec_tag, PS_BASE)
        smds_raw, smds_matched, labels = [], [], []
        for v in ps_vars:
            if v not in all_pts.columns:
                continue
            sr = compute_smd(
                raw_trt[v].dropna().astype(float), raw_ctl[v].dropna().astype(float)
            )
            sm = compute_smd(
                trt[v].dropna().astype(float), ctl[v].dropna().astype(float)
            )
            smds_raw.append(sr)
            smds_matched.append(sm)
            labels.append(v)

        idx = np.argsort(smds_matched)[::-1]
        y_pos = range(len(labels))

        ax.scatter(
            [smds_raw[i] for i in idx],
            y_pos,
            marker="x",
            color=GRAY,
            s=12,
            label="Before",
            zorder=2,
        )
        ax.scatter(
            [smds_matched[i] for i in idx],
            y_pos,
            marker="o",
            color=CLR[db],
            s=12,
            label="After",
            zorder=3,
        )
        ax.axvline(0.1, color="red", linewidth=0.5, linestyle=":", zorder=1)
        ax.set_yticks(y_pos)
        ax.set_yticklabels([labels[i] for i in idx], fontsize=5)
        ax.set_xlabel("|SMD|", fontsize=6)
        ax.set_title(spec_label, fontsize=6)
        ax.set_xlim(-0.02, max(max(smds_raw) if smds_raw else 0.5, 0.5) + 0.05)

        n_viol = sum(1 for s in smds_matched if s > 0.1)
        ax.text(
            0.95,
            0.05,
            f"{n_viol}/{len(labels)} > 0.1",
            transform=ax.transAxes,
            fontsize=5,
            ha="right",
            va="bottom",
            color="red" if n_viol > 0 else GREEN,
        )

        if si == 0:
            ax.legend(fontsize=5, loc="upper right")

    save(fig, "love_mimic")


# ═══════════════════════════════════════════════════════════════════
# eFIG 2: SENSITIVITY DiD COMPARISON
# Primary vs earliest-labs, both DBs
# ═══════════════════════════════════════════════════════════════════
def efig2_sensitivity():
    """Primary vs earliest-labs sensitivity DiD across horizons — reads did_riskset_{db}.csv."""
    fig, axes = plt.subplots(1, 2, figsize=(W_DOUBLE, W_DOUBLE * 0.35), sharey=True)
    spec_styles = {
        "primary": {
            "color": BLUE,
            "marker": "o",
            "label": "Primary (labs closest to T\u2080)",
        },
        "sens_b": {
            "color": GREEN,
            "marker": "^",
            "label": "Sensitivity (earliest labs)",
        },
    }
    for i, db in enumerate(DBS):
        ax = axes[i]
        panel_label(ax, chr(ord("a") + i))
        df = pd.read_csv(os.path.join(RESULTS, f"did_riskset_{db}.csv"))
        df = df[(df.pool == "yet_untreated") & (df.method == "psm_dr")]

        for spec, sty in spec_styles.items():
            sub = df[df.spec == spec].sort_values("target_h")
            if len(sub) == 0:
                continue
            x = sub.target_h.values
            y = sub.did.values
            ax.plot(
                x,
                y,
                f'-{sty["marker"]}',
                color=sty["color"],
                label=sty["label"],
                markersize=3,
                linewidth=0.8,
            )

        ax.axhline(0, color=BLACK, linewidth=0.5, linestyle="--")
        ax.set_xlabel("Hours from T₀")
        ax.set_title(LBL[db], fontsize=7, fontweight="bold")
        ax.set_xticks([6, 12, 18, 24, 30, 36, 42, 48])
        if i == 0:
            ax.set_ylabel("ΔCr DiD (mg/dL)")
            ax.legend(fontsize=5, loc="lower left")

    save(fig, "sensitivity_did")


# ═══════════════════════════════════════════════════════════════════
# eFIG 3: SENSITIVITY — BINARY 48h AKI BY eGFR
# Primary vs earliest-labs PS, 48h AKI (KDIGO ≥1), both DBs
# ═══════════════════════════════════════════════════════════════════
def efig3_sensitivity_binary():
    """Sensitivity analysis: 48h AKI ORs by eGFR under primary vs
    earliest-labs PS specifications. Reads egfr_aki_stages_{db}.csv
    (primary) and egfr_aki_stages_sens_b_{db}.csv (sensitivity)."""
    fig, axes = plt.subplots(1, 2, figsize=(W_DOUBLE, W_DOUBLE * 0.40), sharey=True)

    outcome = "AKI_48h_Stage1+"
    strata_order = [
        "Overall",
        "eGFR>=90",
        "eGFR_60-89",
        "eGFR_45-59",
        "eGFR_30-44",
        "eGFR<30",
    ]
    strata_labels = [
        "Overall",
        "eGFR \u226590",
        "eGFR 60\u201389",
        "eGFR 45\u201359",
        "eGFR 30\u201344",
        "eGFR <30",
    ]
    n_strata = len(strata_order)

    spec_styles = {
        "primary": {
            "marker": "o",
            "label": "Primary (labs closest to T\u2080)",
            "size": 5,
        },
        "sens_b": {"marker": "^", "label": "Sensitivity (earliest labs)", "size": 5},
    }
    spec_files = {
        "primary": "egfr_aki_stages_{db}.csv",
        "sens_b": "egfr_aki_stages_sens_b_{db}.csv",
    }

    for i, db in enumerate(DBS):
        ax = axes[i]
        panel_label(ax, chr(ord("a") + i))
        ax.set_title(LBL[db], fontsize=7, fontweight="bold")

        for si, (spec, sty) in enumerate(spec_styles.items()):
            fname = spec_files[spec].replace("{db}", db)
            fpath = os.path.join(RESULTS, fname)
            if not os.path.exists(fpath):
                print(f"    {db} {spec}: {fname} not found, skipping")
                continue
            df = pd.read_csv(fpath)
            sub = df[df.outcome == outcome].copy()
            sub["stratum"] = pd.Categorical(
                sub.stratum, categories=strata_order, ordered=True
            )
            sub = sub.sort_values("stratum")

            offset = -0.12 if si == 0 else 0.12
            for j, strat in enumerate(strata_order):
                row = sub[sub.stratum == strat]
                if len(row) == 0 or pd.isna(row.iloc[0]["or"]):
                    continue
                r = row.iloc[0]
                y = j + offset
                ax.plot(
                    r["or"],
                    y,
                    sty["marker"],
                    color=CLR[db],
                    markersize=sty["size"],
                    zorder=3,
                    fillstyle="full" if si == 0 else "none",
                )
                ax.plot(
                    [r.or_lo, r.or_hi],
                    [y, y],
                    color=CLR[db],
                    linewidth=0.8,
                    zorder=2,
                    alpha=1.0 if si == 0 else 0.6,
                )

        ax.axvline(1.0, color=BLACK, linewidth=0.5, linestyle="--", zorder=1)
        ax.set_yticks(range(n_strata))
        if i == 0:
            ax.set_yticklabels(strata_labels, fontsize=6)
        ax.set_xlabel("OR (95% CI)", fontsize=6)
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: f"{x:.1f}" if x >= 1 else f"{x:.2f}")
        )
        ax.set_xlim(0.2, 6)
        ax.invert_yaxis()
        ax.axhspan(-0.5, 0.5, color="#f0f0f0", zorder=0)

    # Shared legend
    from matplotlib.lines import Line2D

    handles = []
    for spec, sty in spec_styles.items():
        handles.append(
            Line2D(
                [0],
                [0],
                marker=sty["marker"],
                color="grey",
                linestyle="",
                markersize=5,
                label=sty["label"],
                fillstyle="full" if spec == "primary" else "none",
            )
        )
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=2,
        bbox_to_anchor=(0.5, -0.06),
        fontsize=6,
    )

    save(fig, "sensitivity_binary")


# ═══════════════════════════════════════════════════════════════════
# eFIG 4: HTE SUBGROUP FOREST
# Pre-specified subgroups × aki_7d, both DBs
# ═══════════════════════════════════════════════════════════════════
def efig4_hte_forest():
    """Subgroup forest for 7d AKI — reads did_hte_{db}.csv."""
    subgroups = [
        ("Overall", "Overall"),
        ("Age < 65", "Age <65"),
        ("Age >= 65", "Age ≥65"),
        ("eGFR >= 90", "eGFR ≥90 (G1)"),
        ("eGFR 60-89", "eGFR 60–89 (G2)"),
        ("eGFR 45-59", "eGFR 45–59 (G3a)"),
        ("eGFR 30-44", "eGFR 30–44 (G3b)"),
        ("eGFR < 30", "eGFR <30 (G4–5)"),
        ("CABG", "CABG"),
        ("Non-CABG", "Non-CABG"),
        ("Diabetes", "Diabetes"),
        ("No diabetes", "No diabetes"),
        ("Heart failure", "Heart failure"),
        ("No HF", "No HF"),
        ("BMI >= 30", "BMI ≥30"),
        ("BMI < 30", "BMI <30"),
        ("Female", "Female"),
        ("Male", "Male"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(W_DOUBLE, W_DOUBLE * 0.65), sharey=True)

    for i, db in enumerate(DBS):
        ax = axes[i]
        panel_label(ax, chr(ord("a") + i))
        df = pd.read_csv(os.path.join(RESULTS, f"did_hte_{db}.csv"))
        df = df[df.outcome == "aki_7d"].copy()

        y_pos = []
        labels_used = []
        for j, (sub_name, label) in enumerate(subgroups):
            row = df[df.subgroup == sub_name]
            if len(row) == 0:
                continue
            r = row.iloc[0]
            y_pos.append(j)
            labels_used.append(label)
            if pd.isna(r["or"]):
                continue
            ax.plot(r["or"], j, MKR[db], color=CLR[db], markersize=4, zorder=3)
            ax.plot([r.or_lo, r.or_hi], [j, j], color=CLR[db], linewidth=0.8, zorder=2)

        ax.axvline(1.0, color=BLACK, linewidth=0.5, linestyle="--", zorder=1)
        ax.set_yticks(y_pos)
        if i == 0:
            ax.set_yticklabels(labels_used, fontsize=5)
        ax.set_xlabel("OR (95% CI)")
        ax.set_title(LBL[db], fontsize=7, fontweight="bold")
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: f"{x:.1f}" if x >= 1 else f"{x:.2f}")
        )
        ax.set_xlim(0.3, 3.5)
        ax.invert_yaxis()
        # Shade overall row
        ax.axhspan(-0.5, 0.5, color="#f0f0f0", zorder=0)

    save(fig, "hte_forest")


# ═══════════════════════════════════════════════════════════════════
# eFIG 5: CROSSED PHENOTYPE FOREST
# ═══════════════════════════════════════════════════════════════════
def efig5_crossed_forest():
    """Crossed phenotype forest — reads did_hte_crossed_{db}.csv."""
    fig, axes = plt.subplots(1, 2, figsize=(W_DOUBLE, W_DOUBLE * 0.55), sharey=False)

    for i, db in enumerate(DBS):
        ax = axes[i]
        panel_label(ax, chr(ord("a") + i))
        path = os.path.join(RESULTS, f"did_hte_crossed_{db}.csv")
        if not os.path.exists(path):
            ax.set_title(f"{LBL[db]} (not found)", fontsize=6)
            continue
        df = pd.read_csv(path)
        if "outcome" in df.columns:
            df = df[df.outcome == "aki_7d"]
        df = df.sort_values("or", ascending=True).reset_index(drop=True)

        for j in range(len(df)):
            r = df.iloc[j]
            if pd.isna(r["or"]):
                continue
            color = GREEN if r["or"] < 1 else VERMIL
            ax.plot(r["or"], j, MKR[db], color=color, markersize=4, zorder=3)
            ax.plot([r.or_lo, r.or_hi], [j, j], color=color, linewidth=0.8, zorder=2)

        ax.axvline(1.0, color=BLACK, linewidth=0.5, linestyle="--", zorder=1)
        pheno_col = "phenotype" if "phenotype" in df.columns else df.columns[0]
        ax.set_yticks(range(len(df)))
        ax.set_yticklabels(df[pheno_col].values, fontsize=4.5)
        ax.set_xlabel("OR (95% CI)")
        ax.set_title(LBL[db], fontsize=7, fontweight="bold")
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: f"{x:.1f}" if x >= 1 else f"{x:.2f}")
        )

    save(fig, "crossed_forest")


# ═══════════════════════════════════════════════════════════════════
# FIG 5 DUAL: eGFR × Mg HEATMAP with AKI + MORTALITY rows
# ═══════════════════════════════════════════════════════════════════
def fig5_egfr_mg_heatmap_dual(aki_outcome="aki_48h", name="egfr_mg_heatmap"):
    """eGFR x Mg OR heatmap — 2x2: top = AKI, bottom = mortality."""
    import matplotlib.colors as mcolors
    from matplotlib.patches import Rectangle

    outcomes = [
        (
            aki_outcome,
            {"aki_48h": "48-h AKI", "aki_7d": "7-d AKI"}.get(aki_outcome, aki_outcome),
        ),
        ("mortality", "Hospital mortality"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(W_DOUBLE, W_DOUBLE * 0.58))

    im = None
    for row_i, (oc, oc_label) in enumerate(outcomes):
        for col_i, db in enumerate(DBS):
            df = pd.read_csv(os.path.join(RESULTS, f"mg_strat_{db}.csv"))
            ax = axes[row_i, col_i]
            panel_label(ax, chr(ord("a") + row_i * 2 + col_i))
            sub = df[(df.outcome == oc) & (df.egfr_strat != "All")].copy()
            if len(sub) == 0:
                ax.text(
                    0.5,
                    0.5,
                    "no data",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                    fontsize=6,
                )
                continue
            mg_bins = [
                b
                for b in ["Mg<1.6", "Mg_1.6-2.0", "Mg>=2.0", "Mg<1.8", "Mg>=1.8"]
                if b in sub.mg_strat.values
            ]
            egfr_bins = [
                b
                for b in ["eGFR>=90", "eGFR_60-89", "eGFR_45-59", "eGFR<45"]
                if b in sub.egfr_strat.values
            ]
            if not mg_bins or not egfr_bins:
                ax.text(
                    0.5,
                    0.5,
                    "insufficient strata",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                    fontsize=6,
                )
                continue
            mat = np.full((len(mg_bins), len(egfr_bins)), np.nan)
            pmat = np.full_like(mat, np.nan)
            nmat = np.full((len(mg_bins), len(egfr_bins)), 0)
            rmat_t = np.full_like(mat, np.nan)
            rmat_c = np.full_like(mat, np.nan)
            for mi, mg in enumerate(mg_bins):
                for ei, eg in enumerate(egfr_bins):
                    row = sub[(sub.mg_strat == mg) & (sub.egfr_strat == eg)]
                    if len(row) > 0:
                        r = row.iloc[0]
                        mat[mi, ei] = r["or"]
                        pmat[mi, ei] = r["p"]
                        nmat[mi, ei] = int(r["n"])
                        rmat_t[mi, ei] = r["rate_trt"]
                        rmat_c[mi, ei] = r["rate_ctl"]
            mat_clipped = np.clip(mat, 1 / 3.0, 3.0)
            norm = mcolors.LogNorm(vmin=1 / 3.0, vmax=3.0)
            im = ax.imshow(mat_clipped, cmap="RdBu_r", norm=norm, aspect="auto")
            for mi in range(len(mg_bins)):
                for ei in range(len(egfr_bins)):
                    v = mat[mi, ei]
                    p = pmat[mi, ei]
                    n = nmat[mi, ei]
                    rt = rmat_t[mi, ei]
                    rc = rmat_c[mi, ei]
                    if pd.isna(v):
                        continue
                    events_trt = n * rt if not pd.isna(rt) else 0
                    events_ctl = n * rc if not pd.isna(rc) else 0
                    if min(events_trt, events_ctl) < 20:
                        rect = Rectangle(
                            (ei - 0.5, mi - 0.5),
                            1,
                            1,
                            fill=False,
                            hatch="///",
                            linewidth=0,
                            edgecolor="gray",
                            alpha=0.5,
                        )
                        ax.add_patch(rect)
                    sig = "*" if not pd.isna(p) and p < 0.05 else ""
                    txt_color = "white" if (v < 0.4 or v > 2.5) else "black"
                    ax.text(
                        ei,
                        mi,
                        f"{v:.2f}{sig}\nn={n}",
                        ha="center",
                        va="center",
                        fontsize=4.5,
                        color=txt_color,
                        fontweight="bold",
                    )
            ax.set_xticks(range(len(egfr_bins)))
            ax.set_xticklabels(
                [b.replace("eGFR", "").replace("_", " ") for b in egfr_bins], fontsize=5
            )
            ax.set_yticks(range(len(mg_bins)))
            ax.set_yticklabels([b.replace("Mg_", "Mg ") for b in mg_bins], fontsize=5)
            if row_i == 1:
                ax.set_xlabel("eGFR stratum", fontsize=6)
            if col_i == 0:
                ax.set_ylabel("Baseline Mg", fontsize=6)
            ax.set_title(f"{LBL[db]} \u2014 {oc_label}", fontsize=6, fontweight="bold")

    if im is not None:
        cbar = fig.colorbar(im, ax=axes, shrink=0.6, pad=0.02)
        cbar.set_ticks([1 / 3.0, 0.5, 1.0, 2.0, 3.0])
        cbar.set_ticklabels(["0.33", "0.5", "1.0", "2.0", "3.0"])
        cbar.set_label("OR", fontsize=6)
        cbar.ax.tick_params(labelsize=5)
    save(fig, name)


# ═══════════════════════════════════════════════════════════════════
# INTERNAL: ALL ENDPOINTS × eGFR (for group circulation)
# ═══════════════════════════════════════════════════════════════════
def fig_internal_all_endpoints():
    """Comprehensive eGFR-stratified forest for ALL computed endpoints.
    For internal circulation only (not for manuscript).
    """
    outcome_labels = {
        "AKI_48h_Stage1+": "48h AKI \u22651",
        "AKI_48h_Stage2+": "48h AKI \u22652",
        "AKI_48h_Stage3+": "48h AKI \u22653",
        "AKI_Stage1+": "7d AKI \u22651",
        "AKI_Stage2+": "7d AKI \u22652",
        "AKI_Stage3+": "7d AKI \u22653",
        "hosp_mortality": "Hospital mortality",
        "poaf": "POAF",
        "vent_arrhythmia": "Vent. arrhythmia",
    }
    strata_order = [
        "Overall",
        "eGFR>=90",
        "eGFR_60-89",
        "eGFR_45-59",
        "eGFR_30-44",
        "eGFR<30",
    ]
    strata_labels = ["Overall", "\u226590", "60-89", "45-59", "30-44", "<30"]
    n_strata = len(strata_order)

    fig, axes = plt.subplots(
        len(outcome_labels),
        2,
        figsize=(W_DOUBLE, 1.4 * len(outcome_labels)),
        sharey=True,
        sharex="col",
    )

    for col_i, db in enumerate(DBS):
        df = pd.read_csv(os.path.join(RESULTS, f"egfr_aki_stages_{db}.csv"))
        for row_i, (outcome, label) in enumerate(outcome_labels.items()):
            ax = axes[row_i, col_i]
            sub = df[df.outcome == outcome].copy()
            sub["stratum"] = pd.Categorical(
                sub.stratum, categories=strata_order, ordered=True
            )
            sub = sub.sort_values("stratum")
            for j, strat in enumerate(strata_order):
                r = sub[sub.stratum == strat]
                if len(r) == 0 or pd.isna(r.iloc[0]["or"]):
                    continue
                r = r.iloc[0]
                color = CLR[db]
                ax.plot(r["or"], j, MKR[db], color=color, markersize=3.5, zorder=3)
                ax.plot(
                    [r.or_lo, r.or_hi], [j, j], color=color, linewidth=0.6, zorder=2
                )
            ax.axvline(1.0, color=BLACK, linewidth=0.4, linestyle="--", zorder=1)
            ax.axhspan(-0.5, 0.5, color="#f0f0f0", zorder=0)
            ax.set_xscale("log")
            ax.set_xlim(0.08, 25)
            ax.xaxis.set_major_formatter(
                mticker.FuncFormatter(lambda x, _: f"{x:.1f}" if x >= 1 else f"{x:.2f}")
            )
            ax.set_yticks(range(n_strata))
            if col_i == 0:
                ax.set_yticklabels(strata_labels, fontsize=5)
                ax.set_ylabel(label, fontsize=5.5, fontweight="bold")
            if row_i == 0:
                ax.set_title(LBL[db], fontsize=7, fontweight="bold")
            if row_i == len(outcome_labels) - 1:
                ax.set_xlabel("OR (95% CI)", fontsize=6)
            ax.invert_yaxis()

    fig.suptitle(
        "All Endpoints: eGFR-Stratified ORs (Internal)",
        fontsize=9,
        fontweight="bold",
        y=1.01,
    )
    plt.tight_layout()
    save(fig, "internal_all_endpoints")


# ═══════════════════════════════════════════════════════════════════
# MAIN — dispatch
# ═══════════════════════════════════════════════════════════════════
ALL_FIGS = {
    "primary_forest": fig2_primary_forest,
    "egfr_forest": fig3_egfr_forest,
    "timecourse": fig4_timecourse,
    "egfr_mg_heatmap": lambda: fig5_egfr_mg_heatmap_dual("aki_48h", "egfr_mg_heatmap"),
    "egfr_mg_heatmap_7d": lambda: fig5_egfr_mg_heatmap_dual(
        "aki_7d", "egfr_mg_heatmap_7d"
    ),
    "egfr_mg_heatmap_aki_only": fig5_egfr_mg_heatmap,  # legacy AKI-only version
    "love": efig1_love,
    "sensitivity_did": efig2_sensitivity,
    "sensitivity_binary": efig3_sensitivity_binary,
    "hte_forest": efig4_hte_forest,
    "crossed_forest": efig5_crossed_forest,
    "internal_all_endpoints": fig_internal_all_endpoints,
}

if __name__ == "__main__":
    print("=" * 60)
    print("04_figures.py — Critical Care submission figures")
    print("=" * 60)

    args = sys.argv[1:]
    targets = args if args else list(ALL_FIGS.keys())

    for name in targets:
        if name in ALL_FIGS:
            print(f"\n  Generating {name}...")
            try:
                ALL_FIGS[name]()
            except Exception as e:
                print(f"  ✗ {name} FAILED: {e}")
                import traceback

                traceback.print_exc()
        else:
            print(f"  Unknown figure: {name}")
            print(f'  Available: {", ".join(ALL_FIGS.keys())}')

    print(f'\n{"=" * 60}')
    print("DONE")
    print(f'{"=" * 60}')
