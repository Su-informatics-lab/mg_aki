#!/usr/bin/env python3
"""
fig_raw_interaction.py — Raw data + subgroup interaction figures

Fig A: Raw ΔCr violin+strip at 12/18/24/30/36h (treated vs control)
Fig B: Subgroup interaction heatmap (AIPW RD across all attributes, both dbs)

Output: ~/mg_aki/figs/
"""

import os
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

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
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.direction": "out",
    "ytick.direction": "out",
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
    "black": "#000000",
    "grey": "#999999",
}

RESULTS = os.path.expanduser("~/mg_aki/results")
FIGS = os.path.expanduser("~/mg_aki/figs")
os.makedirs(FIGS, exist_ok=True)
FIGW2 = 7.205


def save(fig, name):
    for ext in ["pdf", "png"]:
        fig.savefig(os.path.join(FIGS, f"{name}.{ext}"), format=ext, dpi=300)
    print(f"  Saved: {name}.pdf/.png")


# ============================================================================
# FIG A: Raw ΔCr violin + strip
# ============================================================================
def fig_raw_strip():
    print("── Fig A: Raw ΔCr strip plot ──")
    TIME_PTS = [12, 18, 24, 30, 36]
    WINDOW = 6

    for tag, title in [("eicu", "eICU-CRD"), ("mimic", "MIMIC-IV")]:
        cr_path = os.path.join(RESULTS, f"did_cr_all_{tag}.csv")
        trt_path = os.path.join(RESULTS, f"did_treated_{tag}.csv")
        ctl_path = os.path.join(RESULTS, f"did_control_{tag}.csv")
        if not all(os.path.exists(p) for p in [cr_path, trt_path, ctl_path]):
            continue

        trt = pd.read_csv(trt_path)
        ctl = pd.read_csv(ctl_path)
        cr = pd.read_csv(cr_path)
        id_col = (
            "patientunitstayid" if "patientunitstayid" in trt.columns else "stay_id"
        )
        trt["pid"] = trt[id_col]
        ctl["pid"] = ctl[id_col]
        cr_id = "patientunitstayid" if "patientunitstayid" in cr.columns else "stay_id"
        cr["pid"] = cr[cr_id]
        if "labresultoffset" not in cr.columns:
            cr["labresultoffset"] = cr["offset_min"]
        cr["offset_h"] = cr["labresultoffset"] / 60

        # Treatment map
        trt_set = set(trt.pid)
        ctl_set = set(ctl.pid)

        # Cr_pre: first Cr within 0-6h
        pre = cr[(cr.offset_h >= 0) & (cr.offset_h <= 6)].sort_values(
            ["pid", "offset_h"]
        )
        pre = pre.drop_duplicates("pid", keep="first")[["pid", "labresult"]].rename(
            columns={"labresult": "cr_pre"}
        )

        fig, axes = plt.subplots(1, len(TIME_PTS), figsize=(FIGW2, 3.5), sharey=True)

        for ti, th in enumerate(TIME_PTS):
            ax = axes[ti]

            # Cr_post closest to th within ±WINDOW
            post = cr[
                (cr.offset_h >= th - WINDOW) & (cr.offset_h <= th + WINDOW)
            ].copy()
            post["dist"] = abs(post.offset_h - th)
            post = post.sort_values(["pid", "dist"]).drop_duplicates(
                "pid", keep="first"
            )

            m = post.merge(pre, on="pid")
            m["dcr"] = m["labresult"] - m["cr_pre"]
            m["group"] = "Other"
            m.loc[m.pid.isin(trt_set), "group"] = "IV Mg"
            m.loc[m.pid.isin(ctl_set), "group"] = "Control"
            m = m[m.group.isin(["IV Mg", "Control"])]
            m["dcr"] = m["dcr"].clip(-2, 3)

            # Violin
            for gi, (grp, color, xoff) in enumerate(
                [("IV Mg", W["blue"], -0.2), ("Control", W["vermil"], 0.2)]
            ):
                gd = m[m.group == grp]["dcr"].dropna().values
                if len(gd) < 10:
                    continue

                # Violin body
                parts = ax.violinplot(
                    [gd],
                    positions=[xoff],
                    widths=0.35,
                    showextrema=False,
                    showmedians=False,
                )
                for pc in parts["bodies"]:
                    pc.set_facecolor(color)
                    pc.set_alpha(0.25)
                    pc.set_edgecolor("none")

                # Jittered strip (subsample for readability)
                np.random.seed(42 + gi)
                show = np.random.choice(gd, min(200, len(gd)), replace=False)
                jitter = np.random.normal(0, 0.05, len(show))
                ax.scatter(
                    jitter + xoff,
                    show,
                    c=color,
                    s=0.5,
                    alpha=0.15,
                    edgecolors="none",
                    rasterized=True,
                )

                # Median + IQR
                q25, med, q75 = np.percentile(gd, [25, 50, 75])
                ax.plot([xoff - 0.08, xoff + 0.08], [med, med], color=color, lw=1.5)
                ax.plot([xoff, xoff], [q25, q75], color=color, lw=1)

                # N label
                ax.text(
                    xoff, -1.9, f"n={len(gd):,}", ha="center", fontsize=4.5, color=color
                )

            ax.axhline(0, color="grey", lw=0.3, ls="--")
            ax.set_xticks([-0.2, 0.2])
            ax.set_xticklabels(["Mg", "Ctl"], fontsize=5)
            ax.set_title(f"{th}h", fontsize=7, fontweight="bold")
            if ti == 0:
                ax.set_ylabel("ΔCr from baseline (mg/dL)")

            # Mean difference annotation
            mg_med = np.median(m.loc[m.group == "IV Mg", "dcr"].dropna())
            ct_med = np.median(m.loc[m.group == "Control", "dcr"].dropna())
            diff = mg_med - ct_med
            ax.text(
                0, 2.7, f"Δmed={diff:+.3f}", ha="center", fontsize=4.5, color=W["black"]
            )

        fig.suptitle(
            f"{title}: Raw ΔCr distributions by treatment group",
            fontsize=8,
            fontweight="bold",
        )

        # Legend
        from matplotlib.patches import Patch

        fig.legend(
            [Patch(fc=W["blue"], alpha=0.4), Patch(fc=W["vermil"], alpha=0.4)],
            ["IV Mg", "Control"],
            loc="lower center",
            ncol=2,
            fontsize=6,
            bbox_to_anchor=(0.5, -0.02),
        )

        save(fig, f"figA_raw_strip_{tag}")
        plt.close(fig)


# ============================================================================
# FIG B: Subgroup interaction heatmap
# ============================================================================
def fig_interaction_heatmap():
    print("── Fig B: Subgroup interaction heatmap ──")

    # Load subgroup results from both databases
    dfs = {}
    for tag in ["eicu", "mimic"]:
        p = os.path.join(RESULTS, f"did_subgroups_full_{tag}.csv")
        if os.path.exists(p):
            d = pd.read_csv(p)
            d["db"] = tag
            dfs[tag] = d

    if not dfs:
        print("  No subgroup data")
        return

    # Focus on AKI KDIGO>=1 and key secondary outcomes
    outcomes_show = [
        "AKI KDIGO>=1",
        "Hospital mortality",
        "POAF",
        "Encephalopathy",
        "Vent arrhythmia",
    ]

    for db_tag, db_title in [("mimic", "MIMIC-IV"), ("eicu", "eICU-CRD")]:
        if db_tag not in dfs:
            continue
        df = dfs[db_tag]

        # Filter to AKI KDIGO>=1 (primary binary), exclude Overall
        k1 = df[(df.outcome == "AKI KDIGO>=1") & (df.subgroup != "Overall")].copy()
        if "rd" not in k1.columns:
            print(f"  {db_tag}: no rd column")
            continue
        k1 = k1[k1.rd.notna()].copy()

        # Build label
        k1["label"] = k1.apply(
            lambda r: f"{r['subgroup']}: {r['level']}"
            + (" ●" if str(r.get("ref", "")) == "ref" else ""),
            axis=1,
        )

        # Sort by RD (most protective first)
        k1 = k1.sort_values("rd").reset_index(drop=True)

        fig, ax = plt.subplots(figsize=(FIGW2 * 0.7, len(k1) * 0.28 + 1.5))

        y = np.arange(len(k1))
        colors = []
        for _, r in k1.iterrows():
            is_ref = str(r.get("ref", "")) == "ref"
            if is_ref:
                colors.append(W["grey"])
            elif r["rd"] < 0 and r["p"] < 0.05:
                colors.append(W["blue"])
            elif r["rd"] > 0 and r["p"] < 0.05:
                colors.append(W["vermil"])
            elif r["rd"] < 0:
                colors.append(W["skyblue"])
            else:
                colors.append(W["orange"])

        # Horizontal bars
        bars = ax.barh(
            y, k1["rd"].values, color=colors, height=0.6, alpha=0.8, edgecolor="none"
        )

        # Error bars
        for i, (_, r) in enumerate(k1.iterrows()):
            if pd.notna(r.get("rd_lo")) and pd.notna(r.get("rd_hi")):
                ax.plot(
                    [r["rd_lo"], r["rd_hi"]], [i, i], color="black", lw=0.4, zorder=3
                )

        # Annotations: rate_trt vs rate_ctl, n, P, NNT
        for i, (_, r) in enumerate(k1.iterrows()):
            is_ref = str(r.get("ref", "")) == "ref"
            # N
            n_str = f"n={int(r['n_trt'])}+{int(r['n_ctl'])}"
            ax.text(
                max(k1["rd_hi"].max(), 5) + 0.5,
                i,
                n_str,
                fontsize=4.5,
                va="center",
                color="grey",
            )
            # Rates
            if pd.notna(r.get("rate_trt")):
                rate_str = f"{r['rate_trt']:.0f}% vs {r['rate_ctl']:.0f}%"
                ax.text(
                    max(k1["rd_hi"].max(), 5) + 6,
                    i,
                    rate_str,
                    fontsize=4.5,
                    va="center",
                    color="grey",
                )
            # Significance star
            if not is_ref and pd.notna(r["p"]) and r["p"] < 0.05:
                ax.text(
                    r["rd"] - 0.3 if r["rd"] < 0 else r["rd"] + 0.3,
                    i,
                    "*",
                    fontsize=7,
                    va="center",
                    ha="center",
                    color="white" if abs(r["rd"]) > 2 else W["blue"],
                    fontweight="bold",
                )

        ax.axvline(0, color="black", lw=0.5)
        ax.set_yticks(y)
        ax.set_yticklabels(k1["label"].values, fontsize=5.5)
        ax.set_xlabel("AIPW risk difference, % (negative = protective)")
        ax.set_title(
            f"{db_title}: AKI KDIGO≥1 subgroup stratification\n"
            f"Sorted by effect size (most protective → least)",
            fontsize=7,
            fontweight="bold",
        )
        ax.invert_yaxis()

        # Legend
        from matplotlib.patches import Patch

        legend_elements = [
            Patch(fc=W["blue"], alpha=0.8, label="Protective (P<0.05)"),
            Patch(fc=W["skyblue"], alpha=0.8, label="Protective (NS)"),
            Patch(fc=W["orange"], alpha=0.8, label="Harmful (NS)"),
            Patch(fc=W["vermil"], alpha=0.8, label="Harmful (P<0.05)"),
            Patch(fc=W["grey"], alpha=0.8, label="Reference group"),
        ]
        ax.legend(
            handles=legend_elements,
            loc="lower right",
            fontsize=5,
            title="Direction + significance",
            title_fontsize=5.5,
        )

        save(fig, f"figB_interaction_{db_tag}")
        plt.close(fig)

    # ── Cross-database comparison heatmap ──────────────────────────────
    if len(dfs) == 2:
        print("  Building cross-database heatmap...")
        e = dfs["eicu"]
        m = dfs["mimic"]
        k1_e = e[(e.outcome == "AKI KDIGO>=1")].copy()
        k1_m = m[(m.outcome == "AKI KDIGO>=1")].copy()

        # Build unified key
        k1_e["key"] = k1_e["subgroup"] + "|" + k1_e["level"].astype(str)
        k1_m["key"] = k1_m["subgroup"] + "|" + k1_m["level"].astype(str)
        shared_keys = sorted(set(k1_e.key) & set(k1_m.key))

        if len(shared_keys) < 5:
            print("  Too few shared subgroups for heatmap")
            return

        rd_e = k1_e.set_index("key")["rd"].to_dict()
        rd_m = k1_m.set_index("key")["rd"].to_dict()
        p_e = k1_e.set_index("key")["p"].to_dict()
        p_m = k1_m.set_index("key")["p"].to_dict()

        # Sort by mean RD across databases
        mean_rd = {k: (rd_e.get(k, 0) + rd_m.get(k, 0)) / 2 for k in shared_keys}
        shared_keys = sorted(shared_keys, key=lambda k: mean_rd[k])

        labels = [k.replace("|", ": ") for k in shared_keys]
        mat = np.array(
            [[rd_e.get(k, np.nan), rd_m.get(k, np.nan)] for k in shared_keys]
        )

        fig, ax = plt.subplots(figsize=(4, len(shared_keys) * 0.25 + 1))
        im = ax.imshow(mat, cmap="RdBu_r", aspect="auto", vmin=-15, vmax=15)

        # Significance annotations
        for i, k in enumerate(shared_keys):
            for j, (rd_dict, p_dict) in enumerate([(rd_e, p_e), (rd_m, p_m)]):
                rd_val = rd_dict.get(k, np.nan)
                p_val = p_dict.get(k, np.nan)
                if pd.notna(rd_val):
                    txt = f"{rd_val:+.1f}"
                    if pd.notna(p_val) and p_val < 0.05:
                        txt += "*"
                    color = "white" if abs(rd_val) > 7 else "black"
                    ax.text(
                        j, i, txt, ha="center", va="center", fontsize=4.5, color=color
                    )

        ax.set_xticks([0, 1])
        ax.set_xticklabels(["eICU", "MIMIC"], fontsize=6)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=5)
        ax.set_title(
            "AKI KDIGO≥1: AIPW risk difference (%) by subgroup",
            fontsize=7,
            fontweight="bold",
        )

        cbar = plt.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
        cbar.set_label("RD (%)", fontsize=6)
        cbar.ax.tick_params(labelsize=5)

        save(fig, "figB_heatmap_crossdb")
        plt.close(fig)


# ============================================================================
if __name__ == "__main__":
    fig_raw_strip()
    fig_interaction_heatmap()
    print(f"\nDone. Figures in {FIGS}")
