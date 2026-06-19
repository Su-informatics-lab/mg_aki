#!/usr/bin/env python3
"""
fig_did.py — All publication figures for Mg→AKI DiD (AIPW framework)

Fig 1:  AIPW + sIPTW_DR time course (from t=0, 2 panels)
Fig 2:  Sensitivity forest (5 models × 2 dbs)
Fig 3:  AKI subgroup forest (MIMIC KDIGO≥1)
Fig 4:  Love plot — raw vs weighted SMD (2 panels)
Fig 5:  Spaghetti Cr trajectories (2 panels per db)
Fig S1: Specification curve (256 specs)

Output: ~/mg_aki/figs/
Reads:  ~/mg_aki/results/
"""

import os
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

warnings.filterwarnings("ignore")

# ── Nature style ────────────────────────────────────────────────────────
for k, v in {
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
}.items():
    mpl.rcParams[k] = v

W = {
    "blue": "#0072B2",
    "vermil": "#D55E00",
    "green": "#009E73",
    "orange": "#E69F00",
    "skyblue": "#56B4E9",
    "purple": "#CC79A7",
    "black": "#000000",
    "grey": "#999999",
}

RESULTS = os.path.expanduser("~/mg_aki/results")
FIGS = os.path.expanduser("~/mg_aki/figs")
FIGW2 = 7.205  # double-column 183mm
FIGW1 = 3.504  # single-column 89mm

PS_COVARS = [
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
    "first_heartrate",
    "first_potassium",
    "first_calcium",
    "first_lactate",
    "lactate_missing",
    "first_mg_value",
]

NICE = {
    "age": "Age",
    "is_female": "Female sex",
    "bmi": "BMI",
    "surg_cabg": "CABG",
    "surg_valve": "Valve",
    "surg_combined": "Combined",
    "heart_failure": "Heart failure",
    "hypertension": "Hypertension",
    "diabetes": "Diabetes",
    "ckd": "CKD",
    "copd": "COPD",
    "pvd": "PVD",
    "stroke": "Stroke",
    "liver_disease": "Liver disease",
    "egfr": "eGFR",
    "first_heartrate": "Heart rate",
    "first_potassium": "Potassium",
    "first_calcium": "Calcium",
    "first_lactate": "Lactate",
    "lactate_missing": "Lactate missing",
    "first_mg_value": "Serum Mg",
}

os.makedirs(FIGS, exist_ok=True)


def save(fig, name):
    for ext in ["pdf", "png"]:
        fig.savefig(os.path.join(FIGS, f"{name}.{ext}"), format=ext, dpi=300)
    print(f"  Saved: {name}.pdf/.png")


def load_combined(tag):
    """Load treated + control, stack, median-impute."""
    trt = pd.read_csv(os.path.join(RESULTS, f"did_treated_{tag}.csv"))
    ctl = pd.read_csv(os.path.join(RESULTS, f"did_control_{tag}.csv"))
    id_col = "patientunitstayid" if "patientunitstayid" in trt.columns else "stay_id"
    trt["pid"] = trt[id_col]
    ctl["pid"] = ctl[id_col]
    shared = list(set(trt.columns) & set(ctl.columns))
    df = pd.concat([trt[shared], ctl[shared]], ignore_index=True)
    avail = [c for c in PS_COVARS if c in df.columns]
    for c in avail:
        if df[c].isna().any():
            df[c] = df[c].fillna(df[c].median())
    return df, avail


def compute_smd(df, var, trt_col="treated"):
    x1 = df.loc[df[trt_col] == 1, var].dropna()
    x0 = df.loc[df[trt_col] == 0, var].dropna()
    sp = np.sqrt((x1.var() + x0.var()) / 2)
    return abs(x1.mean() - x0.mean()) / sp if sp > 1e-10 else 0


def compute_wsmd(df, var, w_col="w", trt_col="treated"):
    t1 = df[trt_col] == 1
    t0 = df[trt_col] == 0
    w1 = df.loc[t1, w_col]
    w0 = df.loc[t0, w_col]
    x1 = df.loc[t1, var]
    x0 = df.loc[t0, var]
    m1 = np.average(x1.fillna(0), weights=w1)
    m0 = np.average(x0.fillna(0), weights=w0)
    sp = np.sqrt((x1.var() + x0.var()) / 2)
    return abs(m1 - m0) / sp if sp > 1e-10 else 0


def compute_iptw(df, avail):
    """Fit PS, compute stabilized IPTW, trim 1/99."""
    X = df[avail].values
    y = df["treated"].values
    lr = LogisticRegression(max_iter=1000, C=1e6, solver="lbfgs")
    lr.fit(X, y)
    ps = lr.predict_proba(X)[:, 1]
    ps = np.clip(ps, 0.01, 0.99)
    prev = y.mean()
    w = np.where(y == 1, prev / ps, (1 - prev) / (1 - ps))
    q01, q99 = np.percentile(w, [1, 99])
    w = np.clip(w, q01, q99)
    df = df.copy()
    df["w"] = w
    df["ps"] = ps
    return df


# ============================================================================
# FIG 1: Time course (AIPW + sIPTW_DR)
# ============================================================================
def fig1():
    print("── Fig 1: Time course ──")
    fig, axes = plt.subplots(1, 2, figsize=(FIGW2, 3.0), sharey=True)

    for idx, (db, title) in enumerate([("eicu", "eICU-CRD"), ("mimic", "MIMIC-IV")]):
        ax = axes[idx]
        path = os.path.join(RESULTS, f"did_timecourse_{db}_primary.csv")
        if not os.path.exists(path):
            continue
        df = pd.read_csv(path).dropna(subset=["aipw_did"])

        # Add t=0 anchor
        anchor = pd.DataFrame(
            {
                "target_h": [0],
                "aipw_did": [0],
                "aipw_lo": [0],
                "aipw_hi": [0],
                "aipw_p": [1],
                "siptw_did": [0],
                "siptw_p": [1],
            }
        )
        df = pd.concat([anchor, df], ignore_index=True).sort_values("target_h")

        x = df["target_h"].values
        y = df["aipw_did"].values
        lo = df["aipw_lo"].values
        hi = df["aipw_hi"].values
        sig = df["aipw_p"].values < 0.05

        ax.fill_between(x, lo, hi, color=W["blue"], alpha=0.12, zorder=1)
        ax.plot(x, y, color=W["blue"], lw=1.2, zorder=3, label="AIPW")
        ax.plot(
            x[sig],
            y[sig],
            "o",
            color=W["blue"],
            ms=5,
            markeredgecolor="white",
            markeredgewidth=0.5,
            zorder=4,
        )
        ax.plot(
            x[~sig], y[~sig], "o", color=W["blue"], ms=3, fillstyle="none", zorder=4
        )

        ys = df["siptw_did"].values
        ax.plot(
            x,
            ys,
            color=W["vermil"],
            lw=0.8,
            ls="--",
            marker="s",
            ms=2.5,
            label="sIPTW-DR",
            alpha=0.8,
            zorder=2,
        )

        ax.axhline(0, color="grey", lw=0.5, zorder=0)
        ax.axvline(0, color="grey", lw=0.3, ls=":", zorder=0)
        ax.set_xlabel("Hours from ICU admission")
        ax.set_xticks(range(0, 37, 6))
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

    axes[0].set_ylabel("AIPW DiD estimate, ΔΔCr (mg/dL)")
    axes[0].legend(loc="lower left")
    fig.text(
        0.5,
        -0.02,
        "Filled circles = P < 0.05; open = not significant. "
        "Negative = IV Mg protective.",
        ha="center",
        fontsize=5.5,
        fontstyle="italic",
        color="grey",
    )
    save(fig, "fig1_timecourse")
    plt.close(fig)


# ============================================================================
# FIG 2: Sensitivity forest
# ============================================================================
def fig2():
    print("── Fig 2: Sensitivity forest ──")
    models = [
        ("primary", "Primary (21, labs only)"),
        ("sens_a", "Sens A (18, K+ Ca)"),
        ("sens_b", "Sens B (20, +lactate)"),
        ("sens_c", "Sens C (26, +drugs)"),
        ("sens_d", "Sens D (16, base only)"),
    ]
    rows = []
    for mtag, mlabel in models:
        for db, dbl in [("eicu", "eICU"), ("mimic", "MIMIC")]:
            p = os.path.join(RESULTS, f"did_timecourse_{db}_{mtag}.csv")
            if not os.path.exists(p):
                continue
            d = pd.read_csv(p)
            r24 = d[d.target_h == 24]
            if len(r24) == 0:
                continue
            r = r24.iloc[0]
            rows.append(
                {
                    "model": mlabel,
                    "db": dbl,
                    "did": r["aipw_did"],
                    "lo": r["aipw_lo"],
                    "hi": r["aipw_hi"],
                    "p": r["aipw_p"],
                }
            )
    if not rows:
        return
    rdf = pd.DataFrame(rows)

    fig, axes = plt.subplots(1, 2, figsize=(FIGW2, 2.8), sharey=True)
    for idx, (db, title) in enumerate([("eICU", "eICU-CRD"), ("MIMIC", "MIMIC-IV")]):
        ax = axes[idx]
        sub = rdf[rdf.db == db].reset_index(drop=True)
        for i, r in sub.iterrows():
            c = W["blue"] if "Primary" in r["model"] else W["black"]
            lw = 1.5 if "Primary" in r["model"] else 0.8
            ms = 6 if "Primary" in r["model"] else 4
            ax.errorbar(
                r["did"],
                i,
                xerr=[[r["did"] - r["lo"]], [r["hi"] - r["did"]]],
                fmt="o",
                color=c,
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
        ax.set_yticks(range(len(sub)))
        ax.set_yticklabels([r["model"] for _, r in sub.iterrows()], fontsize=6)
        ax.set_xlabel("AIPW DiD at 24 h (mg/dL)")
        ax.set_title(title, fontweight="bold", fontsize=8)
        ax.text(
            -0.35,
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
# FIG 3: AKI subgroup forest (MIMIC)
# ============================================================================
def fig3():
    print("── Fig 3: AKI subgroup forest ──")
    p = os.path.join(RESULTS, "did_subgroups_full_mimic.csv")
    if not os.path.exists(p):
        return
    df = pd.read_csv(p)
    k1 = df[(df.outcome == "AKI KDIGO>=1") & df["or"].notna()].copy()
    k1 = k1[k1["or"] < 50].reset_index(drop=True)

    # Build display label
    k1["label"] = k1.apply(
        lambda r: f"{r['subgroup']}: {r['level']}"
        + (" (ref)" if str(r.get("ref", "")) == "ref" else ""),
        axis=1,
    )

    fig, ax = plt.subplots(figsize=(FIGW1 * 1.5, len(k1) * 0.22 + 0.8))
    for i, r in k1.iterrows():
        is_ref = str(r.get("ref", "")) == "ref"
        c = "grey" if is_ref else (W["blue"] if r["or"] < 1 else W["vermil"])
        ms = 3 if is_ref else 5
        marker = "D" if is_ref else "o"
        ax.errorbar(
            r["or"],
            i,
            xerr=[[r["or"] - r["or_lo"]], [r["or_hi"] - r["or"]]],
            fmt=marker,
            color=c,
            ms=ms,
            lw=0.8,
            capsize=1.5,
            capthick=0.5,
        )
        if not is_ref and r["p"] < 0.05:
            ax.text(
                max(r["or_hi"], 1.05) + 0.05,
                i,
                "*",
                fontsize=7,
                va="center",
                color=W["blue"],
                fontweight="bold",
            )
    ax.axvline(1, color="grey", lw=0.5, ls="--")
    ax.set_yticks(range(len(k1)))
    ax.set_yticklabels(k1["label"].values, fontsize=5.5)
    ax.set_xlabel("OR (95% CI) for AKI KDIGO ≥ 1")
    ax.set_title(
        "MIMIC-IV: Subgroup analysis (sIPTW-DR)", fontweight="bold", fontsize=8
    )
    ax.set_xlim(0, min(k1["or_hi"].max() * 1.15, 4))
    ax.invert_yaxis()
    save(fig, "fig3_aki_subgroups")
    plt.close(fig)


# ============================================================================
# FIG 4: Love plot (raw vs weighted SMD)
# ============================================================================
def fig4():
    print("── Fig 4: Love plot ──")
    fig, axes = plt.subplots(1, 2, figsize=(FIGW2, 4.0), sharey=True)

    for idx, (tag, title) in enumerate([("eicu", "eICU-CRD"), ("mimic", "MIMIC-IV")]):
        ax = axes[idx]
        try:
            df, avail = load_combined(tag)
            df = compute_iptw(df, avail)
        except Exception as e:
            print(f"  {tag}: {e}")
            continue

        names_plot = [c for c in avail if c in NICE]
        raw_smds = [compute_smd(df, c) for c in names_plot]
        wt_smds = [compute_wsmd(df, c) for c in names_plot]
        labels = [NICE.get(c, c) for c in names_plot]

        # Sort by raw SMD descending
        order = np.argsort(raw_smds)[::-1]
        y = np.arange(len(names_plot))

        for i, oi in enumerate(order):
            ax.plot([raw_smds[oi], wt_smds[oi]], [i, i], color="grey", lw=0.3, zorder=1)
        ax.scatter(
            [raw_smds[o] for o in order],
            y,
            c=W["vermil"],
            s=15,
            marker="o",
            label="Raw",
            zorder=2,
            alpha=0.8,
        )
        ax.scatter(
            [wt_smds[o] for o in order],
            y,
            c=W["blue"],
            s=15,
            marker="s",
            label="sIPTW weighted",
            zorder=3,
        )

        ax.axvline(0.1, color="grey", lw=0.5, ls="--")
        ax.axvline(0.05, color="grey", lw=0.3, ls=":")
        ax.set_yticks(y)
        ax.set_yticklabels([labels[o] for o in order], fontsize=5.5)
        ax.set_xlabel("Absolute SMD")
        ax.set_title(title, fontweight="bold", fontsize=8)
        ax.legend(loc="lower right", fontsize=5.5)
        ax.text(
            -0.30,
            1.03,
            chr(97 + idx),
            transform=ax.transAxes,
            fontsize=8,
            fontweight="bold",
            va="top",
        )
        ax.set_xlim(-0.01, max(max(raw_smds) * 1.1, 0.2))

    save(fig, "fig4_loveplot")
    plt.close(fig)


# ============================================================================
# FIG 5: Spaghetti Cr trajectories (ICU-time anchor)
# ============================================================================
def fig5():
    print("── Fig 5: Spaghetti Cr trajectories ──")

    for tag, title in [("eicu", "eICU-CRD"), ("mimic", "MIMIC-IV")]:
        cr_path = os.path.join(RESULTS, f"did_cr_all_{tag}.csv")
        if not os.path.exists(cr_path):
            continue

        try:
            df, avail = load_combined(tag)
        except:
            continue

        cr = pd.read_csv(cr_path)
        id_col = "patientunitstayid" if "patientunitstayid" in cr.columns else "stay_id"
        cr["pid"] = cr[id_col]
        if "labresultoffset" not in cr.columns:
            cr["labresultoffset"] = cr["offset_min"]
        cr["offset_h"] = cr["labresultoffset"] / 60

        # Cr_pre: first Cr within 0-6h
        pre = cr[(cr.offset_h >= 0) & (cr.offset_h <= 6)].sort_values(
            ["pid", "offset_h"]
        )
        pre = pre.drop_duplicates("pid", keep="first")[["pid", "labresult"]].rename(
            columns={"labresult": "cr_pre"}
        )

        # Merge
        cr2 = cr.merge(pre, on="pid").merge(df[["pid", "treated"]], on="pid")
        cr2["dcr"] = cr2["labresult"] - cr2["cr_pre"]
        cr2 = cr2[(cr2.offset_h >= 0) & (cr2.offset_h <= 48)]
        cr2["dcr"] = cr2["dcr"].clip(-2, 3)  # winsorize for display
        cr2["group"] = cr2["treated"].map({1: "IV Mg", 0: "Control"})

        fig, axes = plt.subplots(1, 2, figsize=(FIGW2, 3.5), sharey=True)

        for gi, (grp, color, marker) in enumerate(
            [("IV Mg", W["blue"], "o"), ("Control", W["vermil"], "s")]
        ):
            ax = axes[gi]
            gd = cr2[cr2.group == grp]

            # Subsample for readability
            pids = gd.pid.unique()
            np.random.seed(42)
            show_pids = np.random.choice(pids, min(400, len(pids)), replace=False)
            gd_show = gd[gd.pid.isin(show_pids)]

            for pid, pdf in gd_show.groupby("pid"):
                pdf = pdf.sort_values("offset_h")
                ax.plot(pdf.offset_h, pdf.dcr, color=color, alpha=0.04, lw=0.2)

            # Mean trend (3h bins)
            gd["tbin"] = (gd.offset_h / 3).round() * 3
            mean_t = gd.groupby("tbin").dcr.agg(["mean", "sem"]).reset_index()
            ax.fill_between(
                mean_t.tbin,
                mean_t["mean"] - 1.96 * mean_t["sem"],
                mean_t["mean"] + 1.96 * mean_t["sem"],
                color=color,
                alpha=0.25,
            )
            ax.plot(mean_t.tbin, mean_t["mean"], color=color, lw=1.5)
            ax.plot(mean_t.tbin, mean_t["mean"], marker, color=color, ms=2.5)

            ax.axhline(0, color="grey", lw=0.5, ls="--")
            ax.set_xlabel("Hours from ICU admission")
            ax.set_xticks(range(0, 49, 6))
            ax.set_title(f"{grp} (n={len(pids):,})", fontweight="bold", fontsize=8)
            ax.text(
                -0.12,
                1.05,
                chr(97 + gi),
                transform=ax.transAxes,
                fontsize=8,
                fontweight="bold",
                va="top",
            )

        axes[0].set_ylabel("ΔCr from baseline (mg/dL)")
        fig.suptitle(
            f"{title}: Individual Cr trajectories", fontsize=8, fontweight="bold"
        )
        save(fig, f"fig5_spaghetti_{tag}")
        plt.close(fig)


# ============================================================================
# FIG S1: Specification curve
# ============================================================================
def figs1():
    print("── Fig S1: Spec curve ──")
    path = os.path.join(RESULTS, "did_sweep.csv")
    if not os.path.exists(path):
        return

    df = pd.read_csv(path)
    s24 = df[df.target_h == 24].dropna(subset=["e_aipw", "m_aipw"]).copy()
    s24 = s24.sort_values("m_aipw").reset_index(drop=True)
    x = np.arange(len(s24))
    both_neg = (s24.e_aipw < 0) & (s24.m_aipw < 0)

    fig, axes = plt.subplots(
        3,
        1,
        figsize=(FIGW2, 5.5),
        gridspec_kw={"height_ratios": [2, 2, 1.2]},
        sharex=True,
    )

    # Panel a: eICU
    ax = axes[0]
    c_e = [W["blue"] if bn else W["vermil"] for bn in both_neg]
    ax.scatter(x, s24.e_aipw, c=c_e, s=4, alpha=0.6, edgecolors="none")
    ax.axhline(0, color="grey", lw=0.5, ls="--")
    ax.set_ylabel("eICU AIPW DiD")
    ax.set_title(
        "Specification curve: 256 covariate sets, 24 h from ICU, AIPW",
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

    # Panel b: MIMIC (sorted)
    ax = axes[1]
    c_m = [W["blue"] if v < 0 else W["vermil"] for v in s24.m_aipw]
    ax.scatter(x, s24.m_aipw, c=c_m, s=4, alpha=0.6, edgecolors="none")
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
    tcols = ["D1", "D2", "D3", "D4", "L1", "L2", "L3", "L4"]
    tlabels = [
        "Chronic\ndrugs",
        "Steroids",
        "ICU\ndrugs",
        "Beta\nblockers",
        "K+",
        "Ca",
        "Lactate",
        "Mg",
    ]
    for ti, tc in enumerate(tcols):
        on = s24[tc].values == 1
        ax.scatter(
            x[on], np.full(on.sum(), ti), s=1.5, c=W["black"], alpha=0.4, marker="|"
        )
    ax.set_yticks(range(len(tcols)))
    ax.set_yticklabels(tlabels, fontsize=5)
    ax.set_xlabel("Specification (sorted by MIMIC estimate)")
    ax.set_ylim(-0.5, len(tcols) - 0.5)
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

    n_conc = both_neg.sum()
    fig.text(
        0.02,
        -0.03,
        f"Blue = both negative ({n_conc}/256, {100*n_conc/256:.0f}%); "
        f"Red = discordant. AIPW both-sig-negative: "
        f"{((s24.e_aipw<0)&(s24.e_aipw_p<0.05)&(s24.m_aipw<0)&(s24.m_aipw_p<0.05)).sum()}/256",
        fontsize=5.5,
        color="grey",
    )
    save(fig, "figs1_speccurve")
    plt.close(fig)


# ============================================================================
if __name__ == "__main__":
    fig1()
    fig2()
    fig3()
    fig4()
    fig5()
    figs1()
    print(f"\nAll figures saved to {FIGS}")
