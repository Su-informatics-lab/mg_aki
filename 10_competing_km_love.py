#!/usr/bin/env python3
"""
10_competing_km_love.py — Competing risk probe + KM + Love plots

Run on Tempest:
  python 10_competing_km_love.py          # everything
  python 10_competing_km_love.py death    # death probe only
  python 10_competing_km_love.py km       # KM only
  python 10_competing_km_love.py love     # Love only
"""

import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

import matplotlib as mpl
import matplotlib.pyplot as plt

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

C_BLUE = "#0072B2"
C_VERMILLION = "#D55E00"
C_BLACK = "#000000"
C_GRAY = "#999999"
C_GREEN = "#009E73"
C_SKYBLUE = "#56B4E9"

W_DOUBLE = 7.205
W_SINGLE = 3.504

RESULTS = os.path.expanduser("~/mg_aki/results")
FIGS = os.path.expanduser("~/mg_aki/figs")
os.makedirs(FIGS, exist_ok=True)

LANDMARK_MIN = 360
FOLLOWUP_MAX_MIN = 7 * 24 * 60

RENAME_MAP = {
    "mg_supplementation": "mg_supp",
    "hosp_mortality": "hospital_mortality",
    "age_num": "age",
    "hx_chf": "heart_failure",
    "hx_hypertension": "hypertension",
    "hx_diabetes": "diabetes",
    "hx_ckd": "ckd",
    "hx_copd": "copd",
    "hx_pvd": "pvd",
    "hx_stroke": "stroke",
    "hx_liver": "liver_disease",
    "baseline_cr": "baseline_creatinine",
    "baseline_egfr": "egfr",
    "nephrotox_loop_diuretic": "loop_diuretics",
    "nephrotox_nsaid": "nsaids",
    "nephrotox_acei_arb": "acei_arb",
    "nephrotox_ppi": "ppi",
    "has_betablocker": "beta_blockers",
    "has_steroid": "steroids",
    "preop_antiarrhythmic": "antiarrhythmics",
    "first_k_value": "first_potassium",
    "first_ca_value": "first_calcium",
    "first_hr": "first_heartrate",
    "has_vasopressor": "vasopressor_6h",
    "nc_fracture": "fracture",
    "neuro_encephalopathy": "encephalopathy",
}

# ◆ 30 PS covariates (egfr replaces baseline_creatinine)
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
    "loop_diuretics",
    "nsaids",
    "acei_arb",
    "ppi",
    "beta_blockers",
    "steroids",
    "antiarrhythmics",
    "first_potassium",
    "first_calcium",
    "first_heartrate",
    "vasopressor_6h",
    "transfusion_6h",
    "first_mg_value",
    "first_lactate",
    "lactate_missing",
]

COVAR_LABELS = {
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
    "loop_diuretics": "Loop diuretics",
    "nsaids": "NSAIDs",
    "acei_arb": "ACEi/ARB",
    "ppi": "PPI",
    "beta_blockers": "\u03b2-blockers",
    "steroids": "Steroids",
    "antiarrhythmics": "Antiarrhythmics",
    "first_potassium": "First K\u207a",
    "first_calcium": "First Ca\u00b2\u207a",
    "first_heartrate": "First HR",
    "vasopressor_6h": "Vasopressor",
    "transfusion_6h": "Transfusion (0\u20136 h)",
    "first_mg_value": "First serum Mg",
    "first_lactate": "First lactate",
    "lactate_missing": "Lactate missing",
}


def standardize(d):
    for old, new in RENAME_MAP.items():
        if old in d.columns and new not in d.columns:
            d = d.rename(columns={old: new})
    if "age" in d.columns and d["age"].dtype == object:
        d["age"] = pd.to_numeric(
            d["age"].astype(str).str.replace(">", ""), errors="coerce"
        ).fillna(90)
    if "surgery_type" in d.columns:
        d["surg_cabg"] = (d.surgery_type == "cabg").astype(int)
        d["surg_valve"] = (d.surgery_type == "valve").astype(int)
        d["surg_combined"] = (d.surgery_type == "combined").astype(int)
    # ◆ Lactate: create missing indicator BEFORE median imputation
    if "first_lactate" in d.columns:
        d["lactate_missing"] = d["first_lactate"].isna().astype(int)
    # Median-impute for PS
    for v in PS_COVARS:
        if v in d.columns and d[v].isna().any():
            d[v] = d[v].fillna(d[v].median())
    return d


def load_cohorts():
    print("Loading cohorts...")
    d_e = standardize(pd.read_csv(os.path.join(RESULTS, "01_analysis_a_cohort.csv")))
    d_m = standardize(pd.read_csv(os.path.join(RESULTS, "04_mimic_cohort.csv")))
    print(f"  eICU: N={len(d_e)}, trt={d_e.mg_supp.sum()}")
    print(f"  MIMIC: N={len(d_m)}, trt={d_m.mg_supp.sum()}")
    return d_e, d_m


# =====================================================================
# SECTION 1: DEATH RATE PROBE
# =====================================================================
def probe_death(d, db_name):
    print(f"\n{'='*65}")
    print(f"DEATH RATE PROBE: {db_name}")
    print(f"{'='*65}")
    has_death_offset = "death_offset_min" in d.columns
    has_aki_time = "aki_time_offset" in d.columns
    for mort_col in ["icu_mortality", "hospital_mortality"]:
        if mort_col in d.columns:
            n = d[mort_col].sum()
            print(f"  {mort_col}: {n}/{len(d)} ({100*n/len(d):.1f}%)")
    n7 = 0
    if has_death_offset:
        d7 = d[d.death_offset_min.notna() & (d.death_offset_min <= FOLLOWUP_MAX_MIN)]
        n7 = len(d7)
        print(f"\n  7-day mortality: {n7}/{len(d)} ({100*n7/len(d):.1f}%)")
        for g, label in [(1, "Mg supp"), (0, "No supp")]:
            sub = d[d.mg_supp == g]
            sub7 = sub[
                sub.death_offset_min.notna()
                & (sub.death_offset_min <= FOLLOWUP_MAX_MIN)
            ]
            print(
                f"    {label}: {len(sub7)}/{len(sub)} ({100*len(sub7)/max(len(sub),1):.1f}%)"
            )
        if "ac_group" in d.columns:
            for g, label in [("mg_k", "Mg+K"), ("k_only", "K-only")]:
                sub = d[d.ac_group == g]
                if len(sub) == 0:
                    continue
                sub7 = sub[
                    sub.death_offset_min.notna()
                    & (sub.death_offset_min <= FOLLOWUP_MAX_MIN)
                ]
                print(
                    f"    AC {label}: {len(sub7)}/{len(sub)} ({100*len(sub7)/max(len(sub),1):.1f}%)"
                )
        if has_aki_time:
            died_no_aki = d[
                (d.death_offset_min.notna())
                & (d.death_offset_min <= FOLLOWUP_MAX_MIN)
                & (d.aki_kdigo1 == 0)
            ]
            aki_pts = d[d.aki_kdigo1 == 1]
            print(f"\n  Competing events:")
            print(f"    Died within 7d, no AKI: {len(died_no_aki)}")
            print(f"    AKI events: {aki_pts.aki_kdigo1.sum()}")
        mort_rate = n7 / len(d)
        if mort_rate >= 0.10:
            verdict = "⚠ COMPETING RISK CONCERN"
        elif mort_rate >= 0.05:
            verdict = "⚡ MODERATE"
        else:
            verdict = "✓ LOW"
        print(f"\n  VERDICT: 7-day mortality = {100*mort_rate:.1f}% → {verdict}")
    return {
        "db": db_name,
        "N": len(d),
        "mort_7d": n7,
        "mort_7d_pct": round(100 * n7 / len(d), 1) if n7 else None,
    }


# =====================================================================
# SECTION 2: KM CURVES
# =====================================================================
def kaplan_meier(time, event):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=int)
    order = np.argsort(time)
    t, e = time[order], event[order]
    event_times = np.unique(t[e == 1])
    if len(event_times) == 0:
        return np.array([0]), np.array([1.0]), np.array([1.0]), np.array([1.0])
    n_risk, n_event = [], []
    for ut in event_times:
        n_risk.append(np.sum(t >= ut))
        n_event.append(np.sum((t == ut) & (e == 1)))
    n_risk = np.array(n_risk)
    n_event = np.array(n_event)
    surv_steps = (n_risk - n_event) / n_risk
    survival = np.cumprod(surv_steps)
    var_terms = np.cumsum(n_event / (n_risk * (n_risk - n_event + 1e-10)))
    se = survival * np.sqrt(var_terms)
    ci_lo = np.maximum(0, survival - 1.96 * se)
    ci_hi = np.minimum(1, survival + 1.96 * se)
    event_times = np.concatenate([[0], event_times])
    survival = np.concatenate([[1.0], survival])
    ci_lo = np.concatenate([[1.0], ci_lo])
    ci_hi = np.concatenate([[1.0], ci_hi])
    return event_times, survival, ci_lo, ci_hi


def log_rank_test(time1, event1, time2, event2):
    t_all = np.concatenate([time1, time2])
    e_all = np.concatenate([event1, event2])
    event_times = np.sort(np.unique(t_all[e_all == 1]))
    if len(event_times) == 0:
        return 0, 1.0
    O1 = E1 = V = 0
    for ut in event_times:
        ar1 = np.sum(time1 >= ut)
        ar2 = np.sum(time2 >= ut)
        n_tot = ar1 + ar2
        if n_tot == 0:
            continue
        ev1 = np.sum((time1 == ut) & (event1 == 1))
        ev2 = np.sum((time2 == ut) & (event2 == 1))
        d = ev1 + ev2
        O1 += ev1
        E1 += d * ar1 / n_tot
        if n_tot > 1:
            V += (d * ar1 * ar2 * (n_tot - d)) / (n_tot**2 * (n_tot - 1))
    if V < 1e-10:
        return 0, 1.0
    chi2 = (O1 - E1) ** 2 / V
    return chi2, 1 - stats.chi2.cdf(chi2, df=1)


def construct_tte_eicu(d, event_col="aki_kdigo1", time_col="aki_time_offset"):
    mask = np.ones(len(d), dtype=bool)
    if "unitdischargeoffset" in d.columns:
        mask &= d.unitdischargeoffset > LANDMARK_MIN
    if "death_offset_min" in d.columns:
        mask &= ~(d.death_offset_min.notna() & (d.death_offset_min <= LANDMARK_MIN))
    if event_col == "aki_kdigo1" and time_col in d.columns:
        mask &= ~(
            (d[event_col] == 1) & d[time_col].notna() & (d[time_col] <= LANDMARK_MIN)
        )
    d_lm = d[mask].copy()
    if event_col == "aki_kdigo1" and time_col in d.columns:
        evt_min = np.where(
            (d_lm[event_col] == 1) & d_lm[time_col].notna(),
            d_lm[time_col] - LANDMARK_MIN,
            np.nan,
        )
    elif event_col == "death_7d":
        evt_min = np.where(
            d_lm.death_offset_min.notna()
            & (d_lm.death_offset_min <= FOLLOWUP_MAX_MIN)
            & (d_lm.death_offset_min > LANDMARK_MIN),
            d_lm.death_offset_min - LANDMARK_MIN,
            np.nan,
        )
    else:
        evt_min = np.full(len(d_lm), np.nan)
    censor_min = (
        d_lm.unitdischargeoffset.values.clip(max=FOLLOWUP_MAX_MIN) - LANDMARK_MIN
        if "unitdischargeoffset" in d_lm.columns
        else np.full(len(d_lm), FOLLOWUP_MAX_MIN - LANDMARK_MIN)
    )
    has_event = ~np.isnan(evt_min)
    time_min = np.where(has_event, evt_min, censor_min)
    return np.maximum(time_min, 0) / (60 * 24), has_event.astype(int), d_lm


def construct_tte_mimic(d, event_col="aki_kdigo1"):
    d = d.copy()
    if "intime" in d.columns and "outtime" in d.columns:
        d["intime"] = pd.to_datetime(d["intime"])
        d["outtime"] = pd.to_datetime(d["outtime"])
        d["los_min"] = (d.outtime - d.intime).dt.total_seconds() / 60
    elif "los_min" not in d.columns:
        d["los_min"] = FOLLOWUP_MAX_MIN
    mask = d.los_min > LANDMARK_MIN
    if "death_offset_min" in d.columns:
        mask &= ~(d.death_offset_min.notna() & (d.death_offset_min <= LANDMARK_MIN))
    d_lm = d[mask].copy()
    if event_col == "death_7d":
        evt_min = np.where(
            d_lm.death_offset_min.notna()
            & (d_lm.death_offset_min <= FOLLOWUP_MAX_MIN)
            & (d_lm.death_offset_min > LANDMARK_MIN),
            d_lm.death_offset_min - LANDMARK_MIN,
            np.nan,
        )
        censor_min = d_lm.los_min.values.clip(max=FOLLOWUP_MAX_MIN) - LANDMARK_MIN
        has_event = ~np.isnan(evt_min)
        time_min = np.where(has_event, evt_min, censor_min)
        return np.maximum(time_min, 0) / (60 * 24), has_event.astype(int), d_lm
    else:
        censor_min = d_lm.los_min.values.clip(max=FOLLOWUP_MAX_MIN) - LANDMARK_MIN
        has_aki = (d_lm.aki_kdigo1 == 1).values
        time_min = np.where(has_aki, censor_min * 0.5, censor_min)
        return np.maximum(time_min, 0) / (60 * 24), has_aki.astype(int), d_lm


def plot_km(
    ax,
    time_days,
    event,
    trt,
    label_trt="Mg supp",
    label_ctrl="No supp",
    color_trt=C_VERMILLION,
    color_ctrl=C_BLUE,
    ylabel="Event-free probability",
):
    mask_t = trt == 1
    mask_c = trt == 0
    t_t, s_t, lo_t, hi_t = kaplan_meier(time_days[mask_t], event[mask_t])
    t_c, s_c, lo_c, hi_c = kaplan_meier(time_days[mask_c], event[mask_c])
    ax.step(
        t_t,
        s_t,
        where="post",
        color=color_trt,
        linewidth=1.2,
        label=f"{label_trt} (n={mask_t.sum()}, events={event[mask_t].sum()})",
    )
    ax.fill_between(t_t, lo_t, hi_t, step="post", alpha=0.15, color=color_trt)
    ax.step(
        t_c,
        s_c,
        where="post",
        color=color_ctrl,
        linewidth=1.2,
        label=f"{label_ctrl} (n={mask_c.sum()}, events={event[mask_c].sum()})",
    )
    ax.fill_between(t_c, lo_c, hi_c, step="post", alpha=0.15, color=color_ctrl)
    chi2, p = log_rank_test(
        time_days[mask_t], event[mask_t], time_days[mask_c], event[mask_c]
    )
    p_str = f"P < .001" if p < 0.001 else f"P = {p:.3f}"
    ax.text(
        0.98,
        0.03,
        f"Log-rank {p_str}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=6,
        color="#555555",
    )
    ax.set_xlabel("Days from 6 h landmark")
    ax.set_ylabel(ylabel)
    ax.set_xlim(0, 7)
    ax.set_ylim(0, 1.02)
    ax.legend(loc="lower left", fontsize=6, handlelength=1.5)
    return chi2, p


def run_km(d_e, d_m):
    print(f"\n{'='*65}\nKM CURVES (6 h landmark)\n{'='*65}")
    fig, axes = plt.subplots(2, 3, figsize=(W_DOUBLE, 5.5))
    fig.subplots_adjust(hspace=0.35, wspace=0.3)
    time_d, evt, d_lm = construct_tte_eicu(d_e, "aki_kdigo1", "aki_time_offset")
    print(f"  eICU landmark: {len(d_lm)} patients")
    plot_km(axes[0, 0], time_d, evt, d_lm.mg_supp.values, ylabel="AKI-free probability")
    axes[0, 0].set_title("eICU: AKI (all-patient)", fontsize=7, fontweight="bold")
    time_d2, evt2, d_lm2 = construct_tte_eicu(d_e, "death_7d", "death_offset_min")
    plot_km(
        axes[0, 1], time_d2, evt2, d_lm2.mg_supp.values, ylabel="Survival probability"
    )
    axes[0, 1].set_title("eICU: Death (all-patient)", fontsize=7, fontweight="bold")
    d_ac_e = d_e[d_e.ac_group.isin(["mg_k", "k_only"])].copy()
    d_ac_e["ac_trt"] = (d_ac_e.ac_group == "mg_k").astype(int)
    time_d3, evt3, d_lm3 = construct_tte_eicu(d_ac_e, "aki_kdigo1", "aki_time_offset")
    plot_km(
        axes[0, 2],
        time_d3,
        evt3,
        d_lm3.ac_trt.values,
        label_trt="Mg+K\u207a",
        label_ctrl="K\u207a-only",
        ylabel="AKI-free probability",
    )
    axes[0, 2].set_title("eICU: AKI (active comparator)", fontsize=7, fontweight="bold")
    time_d4, evt4, d_lm4 = construct_tte_mimic(d_m, "aki_kdigo1")
    print(f"  MIMIC landmark: {len(d_lm4)} patients")
    plot_km(
        axes[1, 0], time_d4, evt4, d_lm4.mg_supp.values, ylabel="AKI-free probability"
    )
    axes[1, 0].set_title("MIMIC: AKI (all-patient)*", fontsize=7, fontweight="bold")
    time_d5, evt5, d_lm5 = construct_tte_mimic(d_m, "death_7d")
    plot_km(
        axes[1, 1], time_d5, evt5, d_lm5.mg_supp.values, ylabel="Survival probability"
    )
    axes[1, 1].set_title("MIMIC: Death (all-patient)", fontsize=7, fontweight="bold")
    d_ac_m = d_m[d_m.ac_group.isin(["mg_k", "k_only"])].copy()
    d_ac_m["ac_trt"] = (d_ac_m.ac_group == "mg_k").astype(int)
    time_d6, evt6, d_lm6 = construct_tte_mimic(d_ac_m, "aki_kdigo1")
    plot_km(
        axes[1, 2],
        time_d6,
        evt6,
        d_lm6.ac_trt.values,
        label_trt="Mg+K\u207a",
        label_ctrl="K\u207a-only",
        ylabel="AKI-free probability",
    )
    axes[1, 2].set_title(
        "MIMIC: AKI (active comparator)*", fontsize=7, fontweight="bold"
    )
    for i, ax in enumerate(axes.flat):
        ax.text(
            -0.12,
            1.06,
            chr(ord("a") + i),
            transform=ax.transAxes,
            fontsize=8,
            fontweight="bold",
            va="top",
            ha="right",
        )
    fig.text(
        0.5,
        0.01,
        "*MIMIC AKI timing approximate",
        ha="center",
        fontsize=5,
        color=C_GRAY,
    )
    path = os.path.join(FIGS, "fig_km_panels.pdf")
    fig.savefig(path, format="pdf")
    fig.savefig(path.replace(".pdf", ".png"), format="png")
    plt.close(fig)
    print(f"  Saved: {path}")


# =====================================================================
# SECTION 3: LOVE PLOTS
# =====================================================================
def compute_smd(x, trt, w=None):
    x = np.asarray(x, dtype=float)
    trt = np.asarray(trt, dtype=int)
    mask = ~np.isnan(x)
    x, trt = x[mask], trt[mask]
    if w is not None:
        w = np.asarray(w, dtype=float)[mask]
    if w is None:
        m1 = np.mean(x[trt == 1])
        m0 = np.mean(x[trt == 0])
    else:
        w1 = w[trt == 1]
        w0 = w[trt == 0]
        m1 = np.average(x[trt == 1], weights=w1) if w1.sum() > 0 else 0
        m0 = np.average(x[trt == 0], weights=w0) if w0.sum() > 0 else 0
    s1 = np.std(x[trt == 1], ddof=1) if (trt == 1).sum() > 1 else 0
    s0 = np.std(x[trt == 0], ddof=1) if (trt == 0).sum() > 1 else 0
    sp = np.sqrt((s1**2 + s0**2) / 2)
    return abs(m1 - m0) / sp if sp > 1e-10 else 0.0


def fit_ps_weights(d, trt_col, covars):
    from sklearn.linear_model import LogisticRegression

    available = [c for c in covars if c in d.columns]
    X = d[available].values.astype(float)
    y = d[trt_col].values.astype(float)
    mask = ~np.isnan(X).any(axis=1) & ~np.isnan(y)
    X_cc, y_cc = X[mask], y[mask]
    lr = LogisticRegression(max_iter=2000, C=1e6, solver="lbfgs")
    lr.fit(X_cc, y_cc)
    ps = np.clip(lr.predict_proba(X_cc)[:, 1], 0.01, 0.99)
    ow = np.where(y_cc == 1, 1 - ps, ps)
    prev = y_cc.mean()
    iptw = np.where(y_cc == 1, prev / ps, (1 - prev) / (1 - ps))
    q01, q99 = np.percentile(iptw, [1, 99])
    iptw = np.clip(iptw, q01, q99)
    full_ow = np.full(len(d), np.nan)
    full_iptw = np.full(len(d), np.nan)
    full_ow[mask] = ow
    full_iptw[mask] = iptw
    return full_ow, full_iptw, mask


def plot_love(ax, d, trt_col, covars, title=""):
    available = [c for c in covars if c in d.columns]
    ow, iptw, mask = fit_ps_weights(d, trt_col, available)
    d_cc = d[mask].copy()
    trt = d_cc[trt_col].values
    ow_cc = ow[mask]
    iptw_cc = iptw[mask]
    smd_raw, smd_ow, smd_iptw, labels = [], [], [], []
    for v in available:
        if not np.issubdtype(d_cc[v].dtype, np.number):
            continue
        smd_raw.append(compute_smd(d_cc[v].values, trt))
        smd_ow.append(compute_smd(d_cc[v].values, trt, ow_cc))
        smd_iptw.append(compute_smd(d_cc[v].values, trt, iptw_cc))
        labels.append(COVAR_LABELS.get(v, v))
    smd_raw = np.array(smd_raw)
    smd_ow = np.array(smd_ow)
    smd_iptw = np.array(smd_iptw)
    labels = np.array(labels)
    order = np.argsort(smd_raw)
    smd_raw, smd_ow, smd_iptw, labels = (
        smd_raw[order],
        smd_ow[order],
        smd_iptw[order],
        labels[order],
    )
    y = np.arange(len(labels))
    ax.axvline(0.10, color=C_GRAY, linestyle="--", linewidth=0.5, alpha=0.7)
    ax.axvline(0.05, color=C_GRAY, linestyle=":", linewidth=0.3, alpha=0.5)
    ax.scatter(
        smd_raw,
        y,
        marker="o",
        facecolors="none",
        edgecolors=C_VERMILLION,
        s=20,
        linewidths=0.7,
        label="Unweighted",
        zorder=3,
    )
    ax.scatter(
        smd_iptw,
        y,
        marker="^",
        facecolors=C_GREEN,
        edgecolors=C_GREEN,
        s=18,
        linewidths=0.5,
        label="IPTW",
        zorder=4,
    )
    ax.scatter(
        smd_ow,
        y,
        marker="s",
        facecolors=C_BLUE,
        edgecolors=C_BLUE,
        s=18,
        linewidths=0.5,
        label="OW",
        zorder=5,
    )
    for i in range(len(y)):
        ax.plot(
            [smd_raw[i], smd_iptw[i]],
            [y[i], y[i]],
            color=C_GRAY,
            linewidth=0.2,
            zorder=1,
        )
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=5)
    ax.set_xlabel("|Standardized mean difference|")
    ax.set_xlim(-0.01, max(0.25, max(smd_raw) * 1.15))
    ax.legend(loc="lower right", fontsize=5, markerscale=0.8, handletextpad=0.3)
    if title:
        ax.set_title(title, fontsize=7, fontweight="bold")
    print(f"    {title}")
    print(
        f"      {len(labels)} covariates: raw max={max(smd_raw):.3f}, OW max={max(smd_ow):.4f}"
    )


def run_love(d_e, d_m):
    print(f"\n{'='*65}\nLOVE PLOTS (raw vs IPTW vs OW)\n{'='*65}")

    fig, axes = plt.subplots(2, 2, figsize=(W_DOUBLE, 8.0))
    fig.subplots_adjust(hspace=0.35, wspace=0.50)
    plot_love(axes[0, 0], d_e, "mg_supp", PS_COVARS, "eICU: All-patient")
    d_ac_e = d_e[d_e.ac_group.isin(["mg_k", "k_only"])].copy()
    d_ac_e["ac_trt"] = (d_ac_e.ac_group == "mg_k").astype(int)
    plot_love(axes[0, 1], d_ac_e, "ac_trt", PS_COVARS, "eICU: Active comparator")
    plot_love(axes[1, 0], d_m, "mg_supp", PS_COVARS, "MIMIC: All-patient")
    d_ac_m = d_m[d_m.ac_group.isin(["mg_k", "k_only"])].copy()
    d_ac_m["ac_trt"] = (d_ac_m.ac_group == "mg_k").astype(int)
    plot_love(axes[1, 1], d_ac_m, "ac_trt", PS_COVARS, "MIMIC: Active comparator")
    for i, ax in enumerate(axes.flat):
        ax.text(
            -0.25,
            1.04,
            chr(ord("a") + i),
            transform=ax.transAxes,
            fontsize=8,
            fontweight="bold",
            va="top",
            ha="right",
        )
    path = os.path.join(FIGS, "fig_love_plots.pdf")
    fig.savefig(path, format="pdf")
    fig.savefig(path.replace(".pdf", ".png"), format="png")
    plt.close(fig)
    print(f"  Saved: {path}")


# =====================================================================
# MAIN
# =====================================================================
if __name__ == "__main__":
    args = [a.lower() for a in sys.argv[1:]]
    run_all = len(args) == 0
    d_e, d_m = load_cohorts()
    if run_all or "death" in args:
        rows = [probe_death(d_e, "eICU"), probe_death(d_m, "MIMIC")]
        pd.DataFrame(rows).to_csv(
            os.path.join(RESULTS, "10_death_rates.csv"), index=False
        )
        print(f"\n  Saved: results/10_death_rates.csv")
    if run_all or "km" in args:
        run_km(d_e, d_m)
    if run_all or "love" in args:
        run_love(d_e, d_m)
    print(f"\n{'='*65}\n10_competing_km_love.py COMPLETE\n{'='*65}")
