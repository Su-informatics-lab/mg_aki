#!/usr/bin/env python3
"""
fig_interactive.py — Interactive Plotly exploration of Mg→AKI patient data

NOTE: AIPW is NOT matching. There are no 1:1 pairs.
      AIPW weights every patient to achieve covariate balance.
      This script shows individual patient ΔCr with all attributes on hover.

Outputs: ~/mg_aki/figs/interactive_{db}.html (one per database)

Run: python fig_interactive.py
"""

import os
import warnings

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from sklearn.linear_model import LogisticRegression

warnings.filterwarnings("ignore")

RESULTS = os.path.expanduser("~/mg_aki/results")
FIGS = os.path.expanduser("~/mg_aki/figs")
os.makedirs(FIGS, exist_ok=True)

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
    "is_female": "Female",
    "bmi": "BMI",
    "surg_cabg": "CABG",
    "surg_valve": "Valve",
    "surg_combined": "Combined",
    "heart_failure": "HF",
    "hypertension": "HTN",
    "diabetes": "DM",
    "ckd": "CKD",
    "copd": "COPD",
    "pvd": "PVD",
    "stroke": "Stroke",
    "liver_disease": "Liver",
    "egfr": "eGFR",
    "first_heartrate": "HR",
    "first_potassium": "K+",
    "first_calcium": "Ca",
    "first_lactate": "Lactate",
    "lactate_missing": "Lactate miss",
    "first_mg_value": "Serum Mg",
}


def build_dataset(tag):
    """Load data, build 24h ICU-anchor ΔCr, compute IPTW weights."""
    trt = pd.read_csv(os.path.join(RESULTS, f"did_treated_{tag}.csv"))
    ctl = pd.read_csv(os.path.join(RESULTS, f"did_control_{tag}.csv"))
    cr_all = pd.read_csv(os.path.join(RESULTS, f"did_cr_all_{tag}.csv"))

    id_col = "patientunitstayid" if "patientunitstayid" in trt.columns else "stay_id"
    trt["pid"] = trt[id_col]
    ctl["pid"] = ctl[id_col]
    cr_id = "patientunitstayid" if "patientunitstayid" in cr_all.columns else "stay_id"
    cr_all["pid"] = cr_all[cr_id]
    if "labresultoffset" not in cr_all.columns:
        cr_all["labresultoffset"] = cr_all["offset_min"]
    cr_all["offset_h"] = cr_all["labresultoffset"] / 60

    # Stack
    shared = list(set(trt.columns) & set(ctl.columns))
    full = pd.concat([trt[shared], ctl[shared]], ignore_index=True)
    avail = [c for c in PS_COVARS if c in full.columns]
    for c in avail:
        if full[c].isna().any():
            full[c] = full[c].fillna(full[c].median())

    # Build ΔCr at 24h ICU-anchor
    pre = cr_all[(cr_all.offset_h >= 0) & (cr_all.offset_h <= 6)].sort_values(
        ["pid", "offset_h"]
    )
    pre = pre.drop_duplicates("pid", keep="first")[["pid", "labresult"]].rename(
        columns={"labresult": "cr_pre"}
    )
    post = cr_all[(cr_all.offset_h >= 18) & (cr_all.offset_h <= 30)].copy()
    post["dist"] = abs(post.offset_h - 24)
    post = post.sort_values(["pid", "dist"]).drop_duplicates("pid", keep="first")
    post = post[["pid", "labresult", "offset_h"]].rename(
        columns={"labresult": "cr_post", "offset_h": "cr_post_h"}
    )
    m = pre.merge(post, on="pid")
    m = m[m.cr_post_h > 0]  # ensure post > pre
    m["delta_cr"] = m["cr_post"] - m["cr_pre"]
    m["aki"] = (m.delta_cr >= 0.3).astype(int)

    df = full.merge(
        m[["pid", "cr_pre", "cr_post", "delta_cr", "aki", "cr_post_h"]], on="pid"
    )

    # IPTW weights
    X = df[avail].values
    y = df["treated"].values
    lr = LogisticRegression(max_iter=1000, C=1e6, solver="lbfgs")
    lr.fit(X, y)
    ps = np.clip(lr.predict_proba(X)[:, 1], 0.01, 0.99)
    prev = y.mean()
    df["w"] = np.where(y == 1, prev / ps, (1 - prev) / (1 - ps))
    q01, q99 = np.percentile(df["w"], [1, 99])
    df["w"] = np.clip(df["w"], q01, q99)
    df["ps"] = ps

    # Derived categories
    df["group"] = df["treated"].map({1: "IV Mg", 0: "Control"})
    df["age_cat"] = pd.cut(
        df["age"], [0, 55, 65, 75, 100], labels=["<55", "55-64", "65-74", "≥75"]
    )
    df["egfr_cat"] = pd.cut(
        df["egfr"], [-1, 45, 60, 90, 999], labels=["<45", "45-59", "60-89", "≥90"]
    )
    df["mg_cat"] = pd.cut(
        df["first_mg_value"],
        [-1, 1.8, 2.0, 2.3, 99],
        labels=["<1.8", "1.8-2.0", "2.0-2.3", ">2.3"],
    )
    df["bmi_cat"] = pd.cut(df["bmi"], [0, 25, 30, 99], labels=["<25", "25-30", "≥30"])
    df["cr_pre_cat"] = pd.cut(
        df["cr_pre"],
        [0, 0.8, 1.0, 1.2, 1.5, 99],
        labels=["≤0.8", "0.8-1.0", "1.0-1.2", "1.2-1.5", ">1.5"],
    )
    surg_map = {(1, 0, 0): "CABG", (0, 1, 0): "Valve", (0, 0, 1): "Combined"}
    df["surgery"] = df.apply(
        lambda r: surg_map.get(
            (r.get("surg_cabg", 0), r.get("surg_valve", 0), r.get("surg_combined", 0)),
            "Other",
        ),
        axis=1,
    )

    # Hover text
    df["hover"] = df.apply(
        lambda r: f"Age: {r['age']:.0f} | {'F' if r['is_female'] else 'M'} | BMI: {r['bmi']:.1f}<br>"
        f"Surgery: {r['surgery']} | eGFR: {r['egfr']:.0f}<br>"
        f"Cr pre: {r['cr_pre']:.2f} → post: {r['cr_post']:.2f} (Δ{r['delta_cr']:+.2f})<br>"
        f"Mg: {r['first_mg_value']:.1f} | K: {r.get('first_potassium',0):.1f} | Ca: {r.get('first_calcium',0):.1f}<br>"
        f"HF:{int(r.get('heart_failure',0))} DM:{int(r.get('diabetes',0))} CKD:{int(r.get('ckd',0))}<br>"
        f"IPTW weight: {r['w']:.2f} | PS: {r['ps']:.3f}<br>"
        f"AKI: {'YES' if r['aki'] else 'no'}",
        axis=1,
    )

    return df, avail


def make_interactive(tag, title):
    print(f"\n── Building interactive: {title} ──")
    df, avail = build_dataset(tag)
    nt = df.treated.sum()
    nc = len(df) - nt
    print(f"  {nt} treated + {nc} control = {len(df)} total")

    colors = {"IV Mg": "#0072B2", "Control": "#D55E00"}

    # ═══════════════════════════════════════════════════════════════════
    # Panel 1: ΔCr strip plot by treatment group, jittered
    # ═══════════════════════════════════════════════════════════════════
    fig1 = go.Figure()
    for grp, color in colors.items():
        sub = df[df.group == grp].sample(
            min(2000, len(df[df.group == grp])), random_state=42
        )
        fig1.add_trace(
            go.Scatter(
                x=sub["group"] + np.random.normal(0, 0.08, len(sub)).astype(str),
                y=sub["delta_cr"],
                mode="markers",
                marker=dict(size=3, color=color, opacity=0.3),
                text=sub["hover"],
                hoverinfo="text",
                name=grp,
                showlegend=True,
            )
        )
    # Add box overlay
    for grp, color in colors.items():
        sub = df[df.group == grp]
        fig1.add_trace(
            go.Box(
                y=sub["delta_cr"],
                name=grp,
                marker_color=color,
                boxmean=True,
                opacity=0.5,
                showlegend=False,
            )
        )
    fig1.add_hline(y=0, line_dash="dash", line_color="grey", line_width=0.5)
    fig1.add_hline(
        y=0.3,
        line_dash="dot",
        line_color="red",
        line_width=0.5,
        annotation_text="KDIGO ≥1 threshold",
    )
    fig1.update_layout(
        title=f"{title}: ΔCr distribution by treatment group (hover for details)",
        yaxis_title="ΔCr from baseline (mg/dL)",
        height=500,
        template="plotly_white",
    )

    # ═══════════════════════════════════════════════════════════════════
    # Panel 2: ΔCr scatter by attribute (dropdown to select x-axis)
    # ═══════════════════════════════════════════════════════════════════
    continuous_vars = [
        "age",
        "egfr",
        "bmi",
        "cr_pre",
        "first_mg_value",
        "first_potassium",
        "first_calcium",
        "first_heartrate",
        "ps",
        "w",
    ]
    nice_cont = {
        "age": "Age",
        "egfr": "eGFR",
        "bmi": "BMI",
        "cr_pre": "Baseline Cr",
        "first_mg_value": "Serum Mg",
        "first_potassium": "K+",
        "first_calcium": "Ca",
        "first_heartrate": "HR",
        "ps": "Propensity Score",
        "w": "IPTW Weight",
    }

    # Subsample for performance
    sub = df.sample(min(3000, len(df)), random_state=42)

    fig2 = go.Figure()
    for i, var in enumerate(continuous_vars):
        if var not in sub.columns:
            continue
        for grp, color in colors.items():
            g = sub[sub.group == grp]
            fig2.add_trace(
                go.Scatter(
                    x=g[var],
                    y=g["delta_cr"],
                    mode="markers",
                    marker=dict(size=3, color=color, opacity=0.3),
                    text=g["hover"],
                    hoverinfo="text",
                    name=grp,
                    visible=(i == 0),
                    legendgroup=grp,
                    showlegend=(i == 0),
                )
            )

    # Dropdown buttons
    buttons = []
    for i, var in enumerate(continuous_vars):
        if var not in sub.columns:
            continue
        vis = [False] * len(fig2.data)
        for j in range(2):  # 2 traces per variable (treated, control)
            idx = i * 2 + j
            if idx < len(vis):
                vis[idx] = True
        buttons.append(
            dict(
                label=nice_cont.get(var, var),
                method="update",
                args=[{"visible": vis}, {"xaxis.title": nice_cont.get(var, var)}],
            )
        )

    fig2.update_layout(
        updatemenus=[
            dict(buttons=buttons, direction="down", x=0.02, y=1.15, showactive=True)
        ],
        title=f"{title}: ΔCr vs patient attributes (dropdown to change x-axis)",
        xaxis_title=nice_cont.get(continuous_vars[0], continuous_vars[0]),
        yaxis_title="ΔCr (mg/dL)",
        height=550,
        template="plotly_white",
    )
    fig2.add_hline(y=0, line_dash="dash", line_color="grey", line_width=0.5)
    fig2.add_hline(y=0.3, line_dash="dot", line_color="red", line_width=0.5)

    # ═══════════════════════════════════════════════════════════════════
    # Panel 3: Subgroup mean ΔCr heatmap (treated vs control)
    # ═══════════════════════════════════════════════════════════════════
    cat_vars = [
        ("age_cat", "Age"),
        ("egfr_cat", "eGFR"),
        ("mg_cat", "Serum Mg"),
        ("surgery", "Surgery"),
        ("cr_pre_cat", "Baseline Cr"),
        ("bmi_cat", "BMI"),
    ]
    binary_vars = [
        ("diabetes", "Diabetes"),
        ("ckd", "CKD"),
        ("heart_failure", "HF"),
        ("is_female", "Female"),
    ]

    rows_heat = []
    for var, label in cat_vars + binary_vars:
        if var not in df.columns:
            continue
        for val in sorted(df[var].dropna().unique()):
            sub_t = df[(df[var] == val) & (df.treated == 1)]
            sub_c = df[(df[var] == val) & (df.treated == 0)]
            if len(sub_t) < 10 or len(sub_c) < 10:
                continue
            rows_heat.append(
                {
                    "Subgroup": f"{label}: {val}",
                    "Treated ΔCr": sub_t.delta_cr.mean(),
                    "Control ΔCr": sub_c.delta_cr.mean(),
                    "Diff (T-C)": sub_t.delta_cr.mean() - sub_c.delta_cr.mean(),
                    "Treated AKI%": 100 * sub_t.aki.mean(),
                    "Control AKI%": 100 * sub_c.aki.mean(),
                    "AKI Diff%": 100 * (sub_t.aki.mean() - sub_c.aki.mean()),
                    "n_treated": len(sub_t),
                    "n_control": len(sub_c),
                }
            )

    hdf = pd.DataFrame(rows_heat)
    if len(hdf) > 0:
        fig3 = make_subplots(
            rows=1,
            cols=2,
            subplot_titles=["Mean ΔCr by subgroup", "AKI rate (%) by subgroup"],
            horizontal_spacing=0.15,
        )

        # ΔCr bars
        fig3.add_trace(
            go.Bar(
                y=hdf["Subgroup"],
                x=hdf["Treated ΔCr"],
                name="IV Mg",
                orientation="h",
                marker_color=colors["IV Mg"],
                opacity=0.7,
                text=[f"{v:.3f}" for v in hdf["Treated ΔCr"]],
                textposition="auto",
            ),
            row=1,
            col=1,
        )
        fig3.add_trace(
            go.Bar(
                y=hdf["Subgroup"],
                x=hdf["Control ΔCr"],
                name="Control",
                orientation="h",
                marker_color=colors["Control"],
                opacity=0.7,
                text=[f"{v:.3f}" for v in hdf["Control ΔCr"]],
                textposition="auto",
            ),
            row=1,
            col=1,
        )

        # AKI rate bars
        fig3.add_trace(
            go.Bar(
                y=hdf["Subgroup"],
                x=hdf["Treated AKI%"],
                name="IV Mg",
                orientation="h",
                marker_color=colors["IV Mg"],
                opacity=0.7,
                showlegend=False,
                text=[f"{v:.1f}%" for v in hdf["Treated AKI%"]],
                textposition="auto",
            ),
            row=1,
            col=2,
        )
        fig3.add_trace(
            go.Bar(
                y=hdf["Subgroup"],
                x=hdf["Control AKI%"],
                name="Control",
                orientation="h",
                marker_color=colors["Control"],
                opacity=0.7,
                showlegend=False,
                text=[f"{v:.1f}%" for v in hdf["Control AKI%"]],
                textposition="auto",
            ),
            row=1,
            col=2,
        )

        fig3.update_layout(
            title=f"{title}: Subgroup comparison (unweighted descriptive)",
            barmode="group",
            height=max(400, len(hdf) * 25 + 100),
            template="plotly_white",
        )

    # ═══════════════════════════════════════════════════════════════════
    # Combine into single HTML
    # ═══════════════════════════════════════════════════════════════════
    out = os.path.join(FIGS, f"interactive_{tag}.html")
    with open(out, "w") as f:
        f.write(
            f"<html><head><title>{title} — Interactive Explorer</title></head><body>\n"
        )
        f.write(f"<h2>{title}: AIPW Patient-Level Explorer</h2>\n")
        f.write(
            f'<p style="color:grey">n = {nt:,} treated + {nc:,} control. '
            f"AIPW = weighting (no 1:1 pairs). Hover for patient details.</p>\n"
        )
        f.write("<h3>1. ΔCr Distribution</h3>\n")
        f.write(fig1.to_html(full_html=False, include_plotlyjs="cdn"))
        f.write("\n<h3>2. ΔCr vs Attributes (use dropdown)</h3>\n")
        f.write(fig2.to_html(full_html=False, include_plotlyjs=False))
        if len(hdf) > 0:
            f.write("\n<h3>3. Subgroup Comparison</h3>\n")
            f.write(fig3.to_html(full_html=False, include_plotlyjs=False))
        f.write("\n</body></html>")

    print(f"  Saved: {out}")


if __name__ == "__main__":
    for tag, title in [("eicu", "eICU-CRD"), ("mimic", "MIMIC-IV")]:
        try:
            make_interactive(tag, title)
        except Exception as e:
            print(f"  {tag} failed: {e}")
    print(f"\nDone. Open HTML files in browser: {FIGS}/interactive_*.html")
