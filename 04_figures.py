#!/usr/bin/env python3
"""
04_figures.py — Publication figures (v4: fixed love plot)

  fig1_primary     Fig 1: ΔCr bar + KM AKI incidence + AKI rates
  fig2_hte         Fig 2: HTE forest (7d AKI — primary outcome)
  fig3_benefit     Fig 3: Benefit-harm spectrum (7d AKI in Panel B)
  efig_timecourse  eFig: ΔCr time course
  efig_sensitivity eFig: Primary vs Sens A
  efig_hte_48h     eFig: 48h AKI HTE forest (supplement)
  efig_love        eFig: Love plot — ALL 19 PS covariates + supp Mg/K

Changes from v3:
  - efig_love: shows all 19 PS covariates (was 15, missing 4 labs)
  - efig_love: sorted by matched SMD descending (violation at top)
  - efig_love: adds Mg/K as supplementary vars below separator

Usage:
  python 04_figures.py                # all
  python 04_figures.py efig_love      # just the love plot
"""

import os
import sys

import matplotlib as mpl
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
    print(f"  \u2713 {name}")


def load_hte(tag):
    p = os.path.join(RESULTS, f"did_hte_{tag}.csv")
    if not os.path.exists(p):
        return None
    df = pd.read_csv(p)
    # Compat: map new 03_hte.R columns to old names expected by figures
    if "or" in df.columns and "est" not in df.columns:
        df["est"] = df["or"]
    if "rate_trt" in df.columns and "rd" not in df.columns:
        df["rd"] = df["rate_trt"] - df["rate_ctl"]
    if "or_lo" in df.columns and "se" not in df.columns:
        df["se"] = (np.log(df["or_hi"].clip(1e-6)) - np.log(df["or_lo"].clip(1e-6))) / (
            2 * 1.96
        )
    if "n" in df.columns and "n_trt" not in df.columns:
        df["n_trt"] = df.get("events_trt", df["n"] // 2)
        df["n_ctl"] = df.get("events_ctl", df["n"] // 2)
    return df


def load_rs(tag):
    p = os.path.join(RESULTS, f"did_riskset_{tag}.csv")
    return pd.read_csv(p) if os.path.exists(p) else None


def hte_row(hte, sg, oc):
    r = hte[(hte.subgroup == sg) & (hte.outcome == oc)]
    return r.iloc[0] if len(r) > 0 else None


def or_ci(row):
    if row is None or pd.isna(row["est"]) or row["est"] <= 0 or pd.isna(row["se"]):
        return np.nan, np.nan, np.nan
    lo = np.exp(np.log(row["est"]) - 1.96 * row["se"])
    hi = np.exp(np.log(row["est"]) + 1.96 * row["se"])
    return row["est"], lo, hi


# ====================================================================
# KM COMPUTATION
# ====================================================================


def ensure_hte_data(tag):
    """Generate did_hte_data_{db}.csv from pairs + did_all if missing."""
    cache = os.path.join(RESULTS, f"did_hte_data_{tag}.csv")
    if os.path.exists(cache):
        return
    print(f"    Generating did_hte_data_{tag}.csv from matched pairs...")
    pairs_path = os.path.join(RESULTS, f"did_pairs_primary_yet_untreated_{tag}.csv")
    all_path = os.path.join(RESULTS, f"did_all_{tag}.csv")
    if not os.path.exists(pairs_path) or not os.path.exists(all_path):
        print(f"    Missing input files for {tag}")
        return
    pairs = pd.read_csv(pairs_path)
    all_pts = pd.read_csv(all_path)
    pid_map = all_pts.set_index("pid")
    rows = []
    for _, pr in pairs.iterrows():
        for pid, trt_val in [(pr.trt_pid, 1), (pr.ctl_pid, 0)]:
            if pid in pid_map.index:
                row = pid_map.loc[pid].to_dict()
                row["pid"] = pid
                row["treated"] = trt_val
                row["t_mg"] = pr.t_mg
                rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(cache, index=False)
    print(
        f"    → {len(df)} rows ({df.treated.sum():.0f} trt, {(1-df.treated).sum():.0f} ctl)"
    )


def compute_km_data(tag):
    cache = os.path.join(RESULTS, f"km_aki_{tag}.csv")
    if os.path.exists(cache):
        df = pd.read_csv(cache)
        if len(df) > 0 and "event" in df.columns:
            return df
        os.remove(cache)

    ensure_hte_data(tag)
    print(f"    Computing KM data for {tag} (one-time)...")
    hte_data = pd.read_csv(os.path.join(RESULTS, f"did_hte_data_{tag}.csv"))
    cr_all = pd.read_csv(os.path.join(RESULTS, f"did_cr_all_{tag}.csv"))

    cr_id = "patientunitstayid" if "patientunitstayid" in cr_all.columns else "stay_id"
    cr_all["pid"] = cr_all[cr_id].astype(int).astype(str)
    if "offset_h" not in cr_all.columns:
        cr_all["offset_h"] = cr_all["labresultoffset"] / 60

    cr_by_pid = {pid: g.sort_values("offset_h") for pid, g in cr_all.groupby("pid")}
    print(f"    Cr data: {len(cr_by_pid)} unique patients")

    rows = []
    n_no_pre = n_no_cr = n_ok = 0
    for _, r in hte_data.iterrows():
        pid = str(int(r["pid"]))
        tmg = r["t_mg"]
        cr_pt = cr_by_pid.get(pid)
        if cr_pt is None:
            n_no_cr += 1
            continue
        pre = cr_pt[(cr_pt.offset_h >= 0) & (cr_pt.offset_h < tmg)]
        if len(pre) == 0:
            n_no_pre += 1
            continue
        cr_pre = pre.iloc[-1].labresult
        if pd.isna(cr_pre) or cr_pre <= 0:
            n_no_pre += 1
            continue
        post = cr_pt[cr_pt.offset_h > tmg].copy()
        post["hours_from_t0"] = post.offset_h - tmg
        post = post[post.hours_from_t0 <= 168]
        time_aki = np.nan
        for _, cr in post.iterrows():
            h = cr.hours_from_t0
            delta = cr.labresult - cr_pre
            ratio = cr.labresult / cr_pre if cr_pre > 0 else 0
            if (h <= 48 and delta >= 0.3) or ratio >= 1.5:
                time_aki = h
                break
        if not np.isnan(time_aki):
            rows.append(
                {"treated": int(r.treated), "time": time_aki, "event": 1, "db": tag}
            )
        else:
            censor_t = post.hours_from_t0.max() if len(post) > 0 else 0
            if censor_t > 0:
                rows.append(
                    {"treated": int(r.treated), "time": censor_t, "event": 0, "db": tag}
                )
        n_ok += 1

    print(
        f"    Processed: {n_ok} ok, {n_no_pre} no Cr_pre, {n_no_cr} pid not in Cr data"
    )
    if len(rows) == 0:
        return None
    df = pd.DataFrame(rows)
    df.to_csv(cache, index=False)
    print(f"    \u2192 {len(df)} patients, {int(df.event.sum())} AKI events")
    return df


def kaplan_meier(times, events):
    n = len(times)
    order = np.argsort(times)
    t_sorted = times[order]
    e_sorted = events[order]
    unique_t = np.unique(t_sorted[e_sorted == 1])
    km_t, km_surv, km_lo, km_hi = [0.0], [1.0], [1.0], [1.0]
    var_sum = 0.0
    surv = 1.0
    at_risk = n
    for ut in unique_t:
        n_censor = np.sum(
            (t_sorted < ut)
            & (e_sorted == 0)
            & (t_sorted > (km_t[-1] if len(km_t) > 1 else 0))
        )
        at_risk -= n_censor
        d = np.sum((t_sorted == ut) & (e_sorted == 1))
        if at_risk <= 0:
            break
        surv *= 1 - d / at_risk
        if at_risk > d:
            var_sum += d / (at_risk * (at_risk - d))
        se = surv * np.sqrt(var_sum) if var_sum > 0 else 0
        km_t.append(ut)
        km_surv.append(surv)
        km_lo.append(max(0, surv - 1.96 * se))
        km_hi.append(min(1, surv + 1.96 * se))
        at_risk -= d
    km_t = np.array(km_t)
    return km_t, 1 - np.array(km_surv), 1 - np.array(km_hi), 1 - np.array(km_lo)


# ====================================================================
# FIGURE 1: PRIMARY RESULT (3 panels)
# ====================================================================
def fig1_primary():
    print("\n\u2500\u2500 Fig 1: Primary result \u2500\u2500")
    fig = plt.figure(figsize=(7.2, 3.2), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.2, 2.2, 1.5])
    ax_cr = fig.add_subplot(gs[0])
    ax_km = fig.add_subplot(gs[1])
    ax_bar = fig.add_subplot(gs[2])
    width = 0.32

    # Panel A: dCr DiD at 36h
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
        ax.bar(
            i,
            r["did"],
            width * 1.5,
            color=CLR[tag],
            alpha=0.85,
            edgecolor=CLR[tag],
            linewidth=0.5,
        )
        ax.errorbar(
            i,
            r["did"],
            yerr=[[r["did"] - r["ci_lo"]], [r["ci_hi"] - r["did"]]],
            fmt="none",
            color="black",
            capsize=3,
            capthick=0.6,
            lw=0.6,
        )
        sig = (
            "***"
            if r["p"] < 0.001
            else "**" if r["p"] < 0.01 else "*" if r["p"] < 0.05 else ""
        )
        ax.text(
            i,
            min(r["ci_lo"], r["did"]) - 0.003,
            f'{r["did"]:+.003f}{sig}',
            ha="center",
            va="top",
            fontsize=6.5,
            fontweight="bold",
            color=CLR[tag],
        )
    ax.set_xticks(range(len(DBS)))
    ax.set_xticklabels([LBL[t] for t in DBS], fontsize=6)
    ax.set_ylabel("DiD: \u0394Cr at 36h\n(mg/dL)")
    ax.text(
        -0.25,
        1.06,
        "a",
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
    )

    # Panel B: KM cumulative AKI incidence
    ax = ax_km
    for tag in DBS:
        km_data = compute_km_data(tag)
        if km_data is None or len(km_data) == 0:
            continue
        for trt_val, trt_lbl, ls, alpha in [
            (1, "IV Mg", "-", 0.15),
            (0, "Control", "--", 0.08),
        ]:
            sub = km_data[km_data.treated == trt_val]
            if len(sub) < 20:
                continue
            t, ci, ci_lo, ci_hi = kaplan_meier(sub.time.values, sub.event.values)
            ax.step(
                t,
                ci * 100,
                where="post",
                color=CLR[tag],
                ls=ls,
                lw=1.2,
                label=f"{LBL[tag]} {trt_lbl}",
            )
            ax.fill_between(
                t, ci_lo * 100, ci_hi * 100, step="post", alpha=alpha, color=CLR[tag]
            )
    ax.axvline(48, color="#ddd", lw=4, alpha=0.5, zorder=0)
    ax.text(
        48,
        ax.get_ylim()[1] * 0.02 if ax.get_ylim()[1] > 0 else 1,
        "48h",
        fontsize=5,
        color=GRAY,
        ha="center",
        va="bottom",
    )
    ax.set_xlabel("Hours from T\u2080")
    ax.set_ylabel("Cumulative AKI\nincidence (%)")
    ax.set_xlim(0, 168)
    ax.set_xticks([0, 24, 48, 72, 96, 120, 144, 168])
    ax.set_xticklabels(["0", "24", "48", "72", "96", "120", "144", "7d"])
    ax.legend(fontsize=5, loc="upper left", ncol=2, handlelength=1.5, columnspacing=0.8)
    ax.text(
        -0.12,
        1.06,
        "b",
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
    )

    # Panel C: AKI rate bars
    ax = ax_bar
    ax.axhline(0, color="#bbb", lw=0.5)
    outcomes = [("aki_48h", "48h AKI"), ("aki_7d", "7d AKI")]
    for oi, (oc, oc_lbl) in enumerate(outcomes):
        for di, tag in enumerate(DBS):
            hte = load_hte(tag)
            if hte is None:
                continue
            r = hte_row(hte, "Overall", oc)
            if r is None or pd.isna(r["rd"]):
                continue
            rd_pct = r["rd"] * 100
            x = oi * 1.2 + di * (width + 0.05)
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
                rd_pct - 0.2,
                label,
                ha="center",
                va="top",
                fontsize=5,
                color=CLR[tag],
                fontweight="bold",
            )
    ax.set_xticks([oi * 1.2 + 0.2 for oi in range(len(outcomes))])
    ax.set_xticklabels([oc[1] for oc in outcomes], fontsize=6.5)
    ax.set_ylabel("Risk Difference (%)")
    ax.text(
        -0.22,
        1.06,
        "c",
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
    )

    save(fig, "fig1_primary")


# ====================================================================
# FIGURE 2: HTE FOREST — 7-day AKI
# ====================================================================
def _hte_forest(outcome, xlabel, figname):
    subgroups = [
        "Overall",
        None,
        "Age < 65",
        "Age >= 65",
        None,
        "eGFR < 60",
        "eGFR >= 60",
        None,
        "Mg < 1.8",
        "Mg >= 1.8",
        None,
        "CABG",
        "Non-CABG",
        None,
        "Diabetes",
        "No diabetes",
        None,
        "CKD",
        "No CKD",
        None,
        "Heart failure",
        "No HF",
        None,
        "BMI >= 30",
        "BMI < 30",
        None,
        "DM + CKD",
        "HF + CABG",
        "Mg<1.8 + CKD",
    ]
    htes = {tag: load_hte(tag) for tag in DBS}
    fig, ax = plt.subplots(figsize=(5.5, 7.5), constrained_layout=True)
    ax.axvline(1, color="#ddd", lw=0.6, zorder=0)
    y = 0
    yticks, ylabels = [], []
    for item in reversed(subgroups):
        if item is None:
            y += 0.6
            continue
        for di, tag in enumerate(DBS):
            hte = htes.get(tag)
            if hte is None:
                continue
            r = hte_row(hte, item, outcome)
            if r is None:
                continue
            est, lo, hi = or_ci(r)
            if pd.isna(est) or est > 15:
                continue
            sig = not pd.isna(r["p"]) and r["p"] < 0.05
            offset = 0.13 * (1 - 2 * di)
            ax.errorbar(
                est,
                y + offset,
                xerr=[[max(est - lo, 0.001)], [max(hi - est, 0.001)]],
                fmt=MKR[tag],
                color=CLR[tag],
                ms=5 if sig else 3.5,
                markerfacecolor=CLR[tag] if sig else "white",
                markeredgecolor=CLR[tag],
                markeredgewidth=0.7,
                capsize=1.5,
                capthick=0.4,
                lw=0.5,
                zorder=3,
            )
        yticks.append(y)
        ylabels.append(item)
        y += 1.0
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=6)
    ax.set_xlabel(xlabel)
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
        0.02,
        -0.02,
        "\u2190 Favors IV Mg",
        transform=ax.transAxes,
        fontsize=5.5,
        color=GRAY,
    )
    ax.text(
        0.98,
        -0.02,
        "Favors control \u2192",
        transform=ax.transAxes,
        fontsize=5.5,
        color=GRAY,
        ha="right",
    )
    save(fig, figname)


def fig2_hte():
    print("\n\u2500\u2500 Fig 2: HTE forest \u2500\u2500")
    _hte_forest("aki_7d", "Odds Ratio for 7-day AKI (95% CI)", "fig2_hte")


def efig_hte_48h():
    print("\n\u2500\u2500 eFig: HTE forest (48h AKI, supplement) \u2500\u2500")
    _hte_forest("aki_48h", "Odds Ratio for 48-hour AKI (95% CI)", "efig_hte_48h")


# ====================================================================
# FIGURE 3: BENEFIT-HARM SPECTRUM
# ====================================================================
def fig3_benefit():
    print("\n\u2500\u2500 Fig 3: Benefit-harm spectrum \u2500\u2500")
    from matplotlib.patches import Patch

    fig, (ax_mg, ax_cross) = plt.subplots(
        1,
        2,
        figsize=(7.0, 3.2),
        gridspec_kw={"width_ratios": [3, 2]},
        constrained_layout=True,
    )
    width = 0.28

    # Panel A: Mg-stratified RD%
    ax = ax_mg
    ax.axhline(0, color="#bbb", lw=0.5)
    mg_groups = ["Mg < 1.8", "Mg >= 1.8"]
    oc_list = [
        ("aki_48h", "48h AKI"),
        ("aki_7d", "7d AKI"),
        ("hosp_mortality", "Mortality"),
    ]
    for gi, mg in enumerate(mg_groups):
        for oi, (oc, _) in enumerate(oc_list):
            for di, tag in enumerate(DBS):
                hte = load_hte(tag)
                if hte is None:
                    continue
                r = hte_row(hte, mg, oc)
                if r is None or pd.isna(r["rd"]):
                    continue
                rd = r["rd"] * 100
                x = oi * 1.6 + (gi * 2 + di) * (width + 0.02)
                alpha = 0.9 if gi == 1 else 0.4
                hatch = "" if gi == 1 else "///"
                ax.bar(
                    x,
                    rd,
                    width,
                    color=CLR[tag],
                    alpha=alpha,
                    hatch=hatch,
                    edgecolor=CLR[tag],
                    linewidth=0.4,
                )
                if abs(rd) > 0.8:
                    ax.text(
                        x,
                        rd,
                        f"{rd:+.1f}",
                        ha="center",
                        va="bottom" if rd > 0 else "top",
                        fontsize=4.5,
                        color=CLR[tag],
                    )
    ax.set_xticks([oi * 1.6 + 0.6 for oi in range(len(oc_list))])
    ax.set_xticklabels([o[1] for o in oc_list], fontsize=6.5)
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
    ax.legend(
        handles=[
            Patch(color=BLUE, alpha=0.9, label="MIMIC Mg\u22651.8"),
            Patch(
                color=BLUE, alpha=0.4, hatch="///", edgecolor=BLUE, label="MIMIC Mg<1.8"
            ),
            Patch(color=VERMIL, alpha=0.9, label="eICU Mg\u22651.8"),
            Patch(
                color=VERMIL,
                alpha=0.4,
                hatch="///",
                edgecolor=VERMIL,
                label="eICU Mg<1.8",
            ),
        ],
        fontsize=4.5,
        loc="lower left",
        ncol=2,
    )

    # Panel B: Crossed phenotype forest — 7-day AKI
    ax = ax_cross
    ax.axvline(1, color="#ddd", lw=0.6, zorder=0)
    crossed = [
        ("HF + CABG", "HF + CABG"),
        ("DM + CKD", "DM + CKD"),
        ("Mg<1.8 + CKD", "Mg<1.8 + CKD"),
    ]
    y = 0
    yticks, ylabels = [], []
    for sg_key, sg_lbl in reversed(crossed):
        for di, tag in enumerate(DBS):
            hte = load_hte(tag)
            if hte is None:
                continue
            r = hte_row(hte, sg_key, "aki_7d")
            if r is None:
                continue
            est, lo, hi = or_ci(r)
            if pd.isna(est) or est > 20:
                continue
            sig = not pd.isna(r["p"]) and r["p"] < 0.05
            offset = 0.12 * (1 - 2 * di)
            ax.errorbar(
                est,
                y + offset,
                xerr=[[max(est - lo, 0.001)], [max(hi - est, 0.001)]],
                fmt=MKR[tag],
                color=CLR[tag],
                ms=6 if sig else 4,
                markerfacecolor=CLR[tag] if sig else "white",
                markeredgecolor=CLR[tag],
                markeredgewidth=0.8,
                capsize=2,
                capthick=0.5,
                lw=0.6,
                zorder=3,
            )
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
        ylabels.append(sg_lbl)
        y += 1.3
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=6.5)
    ax.set_xlabel("OR for 7-day AKI")
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


# ====================================================================
# eFIGURE: TIME COURSE
# ====================================================================
def efig_timecourse():
    print("\n\u2500\u2500 eFig: \u0394Cr time course \u2500\u2500")
    avail = [(t, load_rs(t)) for t in DBS if load_rs(t) is not None]
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
        ax.fill_between(h, lo, hi, alpha=0.15, color=BLUE)
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
        ax.set_xlabel("Hours from T\u2080")
        if i == 0:
            ax.set_ylabel("DiD: \u0394Cr (mg/dL)")
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
        "Primary (19 var), PSM+DR  |  \u25cf = P<0.05  |  Gray = 36h primary  |  DiD<0 = renoprotective",
        ha="center",
        fontsize=5.5,
        color="#666",
    )
    save(fig, "efig_timecourse")


# ====================================================================
# eFIGURE: SENSITIVITY (primary vs sens_a)
# ====================================================================
def efig_sensitivity():
    print("\n\u2500\u2500 eFig: Primary vs Sensitivity A \u2500\u2500")
    avail = [(t, load_rs(t)) for t in DBS if load_rs(t) is not None]
    n = len(avail)
    fig, axes = plt.subplots(
        1, n, figsize=(3.5 * n, 2.8), sharey=True, constrained_layout=True
    )
    if n == 1:
        axes = [axes]
    specs = [
        ("primary", "Primary (no K\u207a/Mg)", BLUE),
        ("sens_a", "+K\u207a/Mg (positivity issue)", VERMIL),
    ]
    for i, (tag, df) in enumerate(avail):
        ax = axes[i]
        ax.axhline(0, color="#aaa", lw=0.5)
        for sn, slbl, sc in specs:
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
            ax.fill_between(h, lo, hi, alpha=0.10, color=sc)
            ax.plot(
                h,
                d,
                color=sc,
                lw=1.2,
                ls=ls,
                marker="o",
                ms=4,
                markerfacecolor="white",
                markeredgecolor=sc,
                markeredgewidth=0.8,
                label=slbl,
            )
            for j in range(len(h)):
                if not np.isnan(pv[j]) and pv[j] < 0.05:
                    ax.plot(h[j], d[j], "o", ms=4, color=sc, zorder=5)
        ax.set_xticks([6, 12, 18, 24, 30, 36, 42, 48])
        ax.set_xlim(3, 51)
        ax.set_xlabel("Hours from T\u2080")
        if i == 0:
            ax.set_ylabel("DiD: \u0394Cr (mg/dL)")
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


# ====================================================================
# eFIGURE: LOVE PLOT — ALL 19 PS covariates + supplementary Mg/K
#   FIXED in v4: was 15 vars (missing 4 labs), now 19 + 2 supp
#   Sorted by matched SMD descending (violation visible at top)
# ====================================================================
# All 19 primary PS covariates
PS_VARS = [
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
    "last_calcium",
    "last_lactate",
    "last_lactate_missing",
    "last_heartrate",
]
# Supplementary (not in PS model, shown below separator)
SUPP_VARS = ["last_magnesium", "last_potassium"]

NICE = {
    "age": "Age",
    "is_female": "Female sex",
    "bmi": "BMI",
    "surg_cabg": "CABG",
    "surg_valve": "Valve surgery",
    "surg_combined": "Combined surgery",
    "heart_failure": "Heart failure",
    "hypertension": "Hypertension",
    "diabetes": "Diabetes",
    "ckd": "CKD",
    "copd": "COPD",
    "pvd": "PVD",
    "stroke": "Stroke/TIA",
    "liver_disease": "Liver disease",
    "egfr": "eGFR",
    "last_calcium": "Calcium*",
    "last_lactate": "Lactate*",
    "last_lactate_missing": "Lactate missing*",
    "last_heartrate": "Heart rate*",
    "last_magnesium": "Magnesium\u2020",
    "last_potassium": "Potassium\u2020",
}


def _smd(x1, x0):
    x1, x0 = x1.dropna(), x0.dropna()
    if len(x1) < 5 or len(x0) < 5:
        return np.nan
    m1, m0 = x1.mean(), x0.mean()
    sp = np.sqrt((x1.var() + x0.var()) / 2)
    return abs(m1 - m0) / sp if sp > 1e-10 else 0.0


def _compute_smds(df, var_list):
    t1, t0 = df[df.treated == 1], df[df.treated == 0]
    return {
        v: _smd(t1[v].astype(float), t0[v].astype(float))
        for v in var_list
        if v in df.columns
    }


def efig_love():
    print("\n\u2500\u2500 eFig: Love plot (covariate balance) \u2500\u2500")
    from matplotlib.lines import Line2D

    n = len(DBS)
    fig, axes = plt.subplots(
        1, n, figsize=(3.8 * n, 5.8), sharey=False, constrained_layout=True
    )
    if n == 1:
        axes = [axes]

    for i, tag in enumerate(DBS):
        ax = axes[i]
        all_path = os.path.join(RESULTS, f"did_all_{tag}.csv")
        hte_path = os.path.join(RESULTS, f"did_hte_data_{tag}.csv")
        ensure_hte_data(tag)
        if not os.path.exists(all_path) or not os.path.exists(hte_path):
            print(f"    SKIP {tag}: missing data files")
            continue

        raw_df = pd.read_csv(all_path)
        matched_df = pd.read_csv(hte_path)

        # Matched SMDs: pair-level (controls duplicated = correct for replacement)
        matched_smds = _compute_smds(matched_df, PS_VARS + SUPP_VARS)

        # Raw SMDs: baseline vars from did_all, lab vars from unique matched patients
        baseline = [v for v in PS_VARS if not v.startswith("last_")]
        raw_smds = _compute_smds(raw_df, baseline)
        unique_matched = matched_df.drop_duplicates(subset="pid")
        lab_raw = _compute_smds(
            unique_matched, [v for v in PS_VARS if v.startswith("last_")] + SUPP_VARS
        )
        raw_smds.update(lab_raw)

        # ── Build rows for PS covariates, sorted by matched SMD desc ──
        ps_rows = []
        for v in PS_VARS:
            m = matched_smds.get(v, np.nan)
            r = raw_smds.get(v, np.nan)
            if np.isnan(m):
                continue
            ps_rows.append((v, NICE.get(v, v), r, m))
        ps_rows.sort(key=lambda x: -x[3])  # sort by matched SMD, largest first

        # Supplementary vars
        supp_rows = []
        for v in SUPP_VARS:
            m = matched_smds.get(v, np.nan)
            r = raw_smds.get(v, np.nan)
            if np.isnan(m):
                continue
            supp_rows.append((v, NICE.get(v, v), r, m))

        all_rows = ps_rows + supp_rows
        n_ps = len(ps_rows)
        n_all = len(all_rows)

        # Y positions: PS at top, gap, supplementary at bottom
        ys = []
        for j in range(n_all):
            if j < n_ps:
                ys.append(n_all - j)
            else:
                ys.append(n_all - j - 0.7)  # gap before supplementary

        # Threshold line
        ax.axvline(0.10, color="#999", lw=0.5, ls="--", zorder=0)

        # Separator between PS and supplementary
        if n_ps < n_all:
            sep_y = ys[n_ps - 1] - 0.5
            ax.axhline(sep_y, color="#ccc", lw=0.4, ls="-")

        # Connecting lines
        for j, (v, lbl, r_val, m_val) in enumerate(all_rows):
            if not np.isnan(r_val):
                ax.plot([r_val, m_val], [ys[j], ys[j]], color="#ddd", lw=0.5, zorder=1)

        # Before-matching crosses
        for j, (v, lbl, r_val, m_val) in enumerate(all_rows):
            if np.isnan(r_val):
                continue
            c = VERMIL if j < n_ps else GRAY
            ax.scatter(
                r_val, ys[j], marker="x", s=22, color=c, linewidths=0.7, zorder=3
            )

        # After-matching dots
        for j, (v, lbl, r_val, m_val) in enumerate(all_rows):
            c = BLUE if j < n_ps else GRAY
            ax.scatter(m_val, ys[j], marker="o", s=20, color=c, linewidths=0, zorder=4)

        ax.set_yticks(ys)
        ax.set_yticklabels([r[1] for r in all_rows], fontsize=5.5)
        ax.set_xlabel("Absolute SMD")
        ax.set_xlim(-0.01, None)
        ax.text(
            -0.02,
            1.03,
            chr(ord("a") + i),
            transform=ax.transAxes,
            fontsize=10,
            fontweight="bold",
            va="top",
        )
        ax.set_title(LBL[tag], fontsize=8, pad=5)

        # Print summary
        viol = [(lbl, m) for _, lbl, _, m in ps_rows if m > 0.10]
        print(
            f"    {LBL[tag]}: {len(ps_rows)} PS covariates, "
            f"{len(viol)} violations"
            + (f" ({', '.join(f'{l}={m:.3f}' for l, m in viol)})" if viol else "")
        )

    # Shared legend
    fig.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker="x",
                color=VERMIL,
                ls="None",
                ms=5,
                markeredgewidth=0.8,
                label="Before matching",
            ),
            Line2D(
                [0],
                [0],
                marker="o",
                color=BLUE,
                ls="None",
                ms=4,
                label="After matching",
            ),
        ],
        loc="lower center",
        ncol=2,
        frameon=False,
        fontsize=6,
        bbox_to_anchor=(0.5, -0.01),
    )

    # Footnote
    fig.text(
        0.5,
        -0.035,
        "* Lab value closest to T\u2080 (PS covariate)   "
        "\u2020 Not in PS model (shown for reference)",
        ha="center",
        fontsize=5,
        style="italic",
        color="#666",
    )

    save(fig, "efig_love_plot")


# ====================================================================
FIGURES = {
    "fig1_primary": fig1_primary,
    "fig2_hte": fig2_hte,
    "fig3_benefit": fig3_benefit,
    "efig_timecourse": efig_timecourse,
    "efig_sensitivity": efig_sensitivity,
    "efig_hte_48h": efig_hte_48h,
    "efig_love": efig_love,
}

if __name__ == "__main__":
    print("=" * 70)
    print("04_figures.py \u2014 Publication figures")
    print("=" * 70)
    args = [a.lower() for a in sys.argv[1:]]
    to_draw = args if args else list(FIGURES.keys())
    for name in to_draw:
        if name in FIGURES:
            FIGURES[name]()
        else:
            print(f"  Unknown: {name}")
    print("=" * 70)
