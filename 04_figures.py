#!/usr/bin/env python3
"""
04_figures.py — Publication figures (v3: KM for AKI)

  fig1_primary     Fig 1: ΔCr bar + KM AKI incidence + AKI rates
  fig2_hte         Fig 2: HTE forest
  fig3_benefit     Fig 3: Benefit-harm spectrum
  efig_timecourse  eFig: ΔCr time course
  efig_sensitivity eFig: Primary vs Sens A

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
    r = hte[(hte.subgroup == sg) & (hte.outcome == oc)]
    return r.iloc[0] if len(r) > 0 else None


def or_ci(row):
    if row is None or pd.isna(row["est"]) or row["est"] <= 0 or pd.isna(row["se"]):
        return np.nan, np.nan, np.nan
    lo = np.exp(np.log(row["est"]) - 1.96 * row["se"])
    hi = np.exp(np.log(row["est"]) + 1.96 * row["se"])
    return row["est"], lo, hi


# ═══════════════════════════════════════════════════════════════════
# KM COMPUTATION: time-to-first-AKI from raw Cr
# ═══════════════════════════════════════════════════════════════════
def compute_km_data(tag):
    """Compute time-to-first-AKI for each matched pair member."""
    cache = os.path.join(RESULTS, f"km_aki_{tag}.csv")
    if os.path.exists(cache):
        df = pd.read_csv(cache)
        if len(df) > 0 and "event" in df.columns:
            return df
        os.remove(cache)  # stale/empty cache

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
    n_no_pre = 0
    n_no_cr = 0
    n_ok = 0
    for _, r in hte_data.iterrows():
        pid = str(int(r["pid"]))
        tmg = r["t_mg"]

        cr_pt = cr_by_pid.get(pid)
        if cr_pt is None:
            n_no_cr += 1
            continue

        # Cr_pre = last Cr before T₀
        pre = cr_pt[(cr_pt.offset_h >= 0) & (cr_pt.offset_h < tmg)]
        if len(pre) == 0:
            n_no_pre += 1
            continue
        cr_pre = pre.iloc[-1].labresult  # already sorted by offset_h
        if pd.isna(cr_pre) or cr_pre <= 0:
            n_no_pre += 1
            continue

        # Post-T₀ Cr measurements
        post = cr_pt[cr_pt.offset_h > tmg].copy()
        post["hours_from_t0"] = post.offset_h - tmg
        post = post[post.hours_from_t0 <= 168]  # 7d max

        # Find first AKI (KDIGO temporal criteria)
        #   ΔCr ≥ 0.3: within 48h only
        #   ratio ≥ 1.5: within 7 days
        time_aki = np.nan
        for _, cr in post.iterrows():
            h = cr.hours_from_t0
            delta = cr.labresult - cr_pre
            ratio = cr.labresult / cr_pre if cr_pre > 0 else 0
            aki_absolute = (h <= 48) and (delta >= 0.3)
            aki_ratio = ratio >= 1.5
            if aki_absolute or aki_ratio:
                time_aki = h
                break

        if not np.isnan(time_aki):
            rows.append(
                {"treated": int(r.treated), "time": time_aki, "event": 1, "db": tag}
            )
        else:
            # Censored at last available Cr or 168h
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
        print("    WARN: no KM rows — returning None")
        return None
    df = pd.DataFrame(rows)
    df.to_csv(cache, index=False)
    print(f"    → {len(df)} patients, {int(df.event.sum())} AKI events")
    return df


def kaplan_meier(times, events):
    """Compute KM cumulative incidence (1 - survival) with Greenwood CI."""
    n = len(times)
    order = np.argsort(times)
    t_sorted = times[order]
    e_sorted = events[order]

    unique_t = np.unique(t_sorted[e_sorted == 1])
    km_t = [0.0]
    km_surv = [1.0]
    km_lo = [1.0]
    km_hi = [1.0]
    var_sum = 0.0
    surv = 1.0
    at_risk = n

    for ut in unique_t:
        # Events and censorings before this time
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
    cum_inc = 1 - np.array(km_surv)
    ci_lo = 1 - np.array(km_hi)  # inverted for cumulative incidence
    ci_hi = 1 - np.array(km_lo)
    return km_t, cum_inc, ci_lo, ci_hi


# ═══════════════════════════════════════════════════════════════════
# FIGURE 1: PRIMARY RESULT (3 panels)
#   A: ΔCr DiD bar at 36h
#   B: KM cumulative AKI incidence (0-7d)
#   C: AKI rates (48h + 7d) by DB
# ═══════════════════════════════════════════════════════════════════
def fig1_primary():
    print("\n── Fig 1: Primary result ──")
    fig = plt.figure(figsize=(7.2, 3.2), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.2, 2.2, 1.5])
    ax_cr = fig.add_subplot(gs[0])
    ax_km = fig.add_subplot(gs[1])
    ax_bar = fig.add_subplot(gs[2])

    width = 0.32

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
    ax.set_ylabel("DiD: ΔCr at 36h\n(mg/dL)")
    ax.text(
        -0.25,
        1.06,
        "a",
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
    )

    # ── Panel B: KM cumulative AKI incidence ──
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
            color = CLR[tag]
            ax.step(
                t,
                ci * 100,
                where="post",
                color=color,
                ls=ls,
                lw=1.2,
                label=f"{LBL[tag]} {trt_lbl}",
            )
            ax.fill_between(
                t, ci_lo * 100, ci_hi * 100, step="post", alpha=alpha, color=color
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
    ax.set_xlabel("Hours from T₀")
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

    # ── Panel C: AKI rate bars ──
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


# ═══════════════════════════════════════════════════════════════════
# FIGURE 2: HTE FOREST
# ═══════════════════════════════════════════════════════════════════
def fig2_hte():
    print("\n── Fig 2: HTE forest ──")
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
    yticks = []
    ylabels = []
    for item in reversed(subgroups):
        if item is None:
            y += 0.6
            continue
        for di, tag in enumerate(DBS):
            hte = htes.get(tag)
            if hte is None:
                continue
            r = hte_row(hte, item, "aki_48h")
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
        0.02, -0.02, "← Favors IV Mg", transform=ax.transAxes, fontsize=5.5, color=GRAY
    )
    ax.text(
        0.98,
        -0.02,
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
    from matplotlib.patches import Patch

    # ── Panel A: Mg-stratified RD% ──
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
        fontsize=4.5,
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
    yticks = []
    ylabels = []
    for sg_key, sg_lbl in reversed(crossed):
        for di, tag in enumerate(DBS):
            hte = load_hte(tag)
            if hte is None:
                continue
            r = hte_row(hte, sg_key, "aki_48h")
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
    ax.set_xlabel("OR for 48h AKI")
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
        "Primary (19 var), PSM+DR  |  ● = P<0.05  |  "
        "Gray = 36h primary  |  DiD<0 = renoprotective",
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
    avail = [(t, load_rs(t)) for t in DBS if load_rs(t) is not None]
    n = len(avail)
    fig, axes = plt.subplots(
        1, n, figsize=(3.5 * n, 2.8), sharey=True, constrained_layout=True
    )
    if n == 1:
        axes = [axes]
    specs = [
        ("primary", "Primary (no K⁺/Mg)", BLUE),
        ("sens_a", "+K⁺/Mg (positivity issue)", VERMIL),
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
            print(f"  Unknown: {name}")
    print("=" * 70)
