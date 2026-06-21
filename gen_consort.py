#!/usr/bin/env python3
"""
eFigure 1: CONSORT Flow Diagram — Risk-set PSM cohort construction
Both databases (eICU-CRD + MIMIC-IV) in parallel columns.
Numbers sourced from authoritative ETL output (01_etl.py).

Output: efig1_consort.pdf / .png
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# ── Nature-style rcParams ──
for k, v in {
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.linewidth": 0.5,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
}.items():
    mpl.rcParams[k] = v

# ── Authoritative numbers from ETL ──
E = {  # eICU-CRD
    "total_icu": "200,859",
    "cardiac": "26,725",
    "eskd_excl": "903",
    "post_eskd": "25,822",
    "ivmg": "6,761",
    "ivmg_pct": "26.2%",
    "no_ivmg": "19,218",
    "cr_pre_icu": "4,025",
    "no_cr_pre": "2,736",
    "hosp_fallback": "923",
    "cr_pre_total": "4,948",
    "excl_cr_high": "103",
    "treated_final": "4,845",
    "ctrl_has_2cr": "14,582",
    "ctrl_final": "14,244",
    "eligible": "19,089",
    "matched": "4,834",
    "match_pct": "99.8%",
    "max_smd": "0.12",
}
M = {  # MIMIC-IV
    "total_icu": "94,458",
    "cardiac": "13,404",
    "eskd_excl": "429",
    "post_eskd": "12,975",
    "ivmg": "6,539",
    "ivmg_pct": "50.4%",
    "no_ivmg": "6,436",
    "cr_pre_icu": "6,381",
    "no_cr_pre": "158",
    "hosp_fallback": "5",
    "cr_pre_total": "6,386",
    "excl_cr_high": "20",
    "treated_final": "6,366",
    "ctrl_has_2cr": "6,272",
    "ctrl_final": "6,211",
    "eligible": "12,577",
    "matched": "6,354",
    "match_pct": "99.8%",
    "max_smd": "0.13",
}

# ── Layout constants ──
FIG_W, FIG_H = 7.5, 10.5  # double-column width, tall
BOX_W = 2.4  # box width
BOX_H = 0.48  # standard box height
BOX_H_TALL = 0.68  # for boxes with more text
EXCL_W = 1.3  # exclusion side-box width
EXCL_H = 0.36

# Column centers — push apart to leave room for exclusion boxes
CX_E = 2.0  # eICU column center
CX_M = 5.5  # MIMIC column center
CX_MID = (CX_E + CX_M) / 2  # midpoint for shared labels

# Colors
CLR_BOX = "#F7F7F7"
CLR_BORDER = "#333333"
CLR_EXCL = "#FFF3E0"
CLR_EXCL_BD = "#E65100"
CLR_FINAL = "#E8F5E9"
CLR_FINAL_BD = "#2E7D32"
CLR_MATCH = "#E3F2FD"
CLR_MATCH_BD = "#1565C0"
CLR_ARROW = "#555555"
CLR_HEADER = "#455A64"

LW = 0.6  # line width for borders


def add_box(
    ax,
    cx,
    cy,
    w,
    h,
    text,
    fc=CLR_BOX,
    ec=CLR_BORDER,
    fontsize=6.5,
    fontweight="normal",
    ha="center",
    va="center",
):
    """Draw a rounded box with centered text."""
    box = FancyBboxPatch(
        (cx - w / 2, cy - h / 2),
        w,
        h,
        boxstyle="round,pad=0.03",
        facecolor=fc,
        edgecolor=ec,
        linewidth=LW,
        zorder=2,
    )
    ax.add_patch(box)
    ax.text(
        cx,
        cy,
        text,
        fontsize=fontsize,
        fontweight=fontweight,
        ha=ha,
        va=va,
        zorder=3,
        linespacing=1.3,
    )
    return cy


def arrow(ax, x1, y1, x2, y2, color=CLR_ARROW):
    """Draw a straight arrow."""
    ax.annotate(
        "",
        xy=(x2, y2),
        xytext=(x1, y1),
        arrowprops=dict(arrowstyle="-|>", color=color, lw=0.8, mutation_scale=8),
        zorder=1,
    )


def side_arrow(ax, x_from, y_from, x_to, y_to, color=CLR_EXCL_BD):
    """Draw a horizontal arrow to an exclusion box."""
    ax.annotate(
        "",
        xy=(x_to, y_to),
        xytext=(x_from, y_from),
        arrowprops=dict(arrowstyle="-|>", color=color, lw=0.6, mutation_scale=6),
        zorder=1,
    )


def build_column(ax, cx, db, nums, x_excl_side):
    """Build one database column of the CONSORT flow."""
    excl_sign = 1 if x_excl_side > cx else -1
    x_excl = cx + excl_sign * (BOX_W / 2 + EXCL_W / 2 + 0.12)

    y = 9.4

    # Row 1: Total ICU
    add_box(
        ax,
        cx,
        y,
        BOX_W,
        BOX_H,
        f"{db}\nTotal ICU admissions\nn = {nums['total_icu']}",
        fontsize=6.5,
        fontweight="bold",
    )

    # Arrow + exclusion: non-cardiac / non-adult / repeat stays
    y_next = y - 0.85
    arrow(ax, cx, y - BOX_H / 2, cx, y_next + BOX_H / 2)
    n_excl_cardiac = int(nums["total_icu"].replace(",", "")) - int(
        nums["cardiac"].replace(",", "")
    )
    add_box(
        ax,
        x_excl,
        y - 0.42,
        EXCL_W,
        EXCL_H,
        f"Excluded: {n_excl_cardiac:,}\n(non-cardiac, <18,\nrepeat stays)",
        fc=CLR_EXCL,
        ec=CLR_EXCL_BD,
        fontsize=5.5,
    )
    side_arrow(
        ax,
        cx + excl_sign * BOX_W / 2,
        y - 0.30,
        x_excl - excl_sign * EXCL_W / 2,
        y - 0.42,
    )

    # Row 2: Cardiac surgery
    y = y_next
    add_box(
        ax,
        cx,
        y,
        BOX_W,
        BOX_H,
        f"Cardiac surgery, adult,\nfirst ICU stay\nn = {nums['cardiac']}",
    )

    # Arrow + exclusion: ESKD
    y_next = y - 0.80
    arrow(ax, cx, y - BOX_H / 2, cx, y_next + BOX_H / 2)
    add_box(
        ax,
        x_excl,
        y - 0.40,
        EXCL_W,
        EXCL_H,
        f"Excluded: {nums['eskd_excl']}\n(ESKD)",
        fc=CLR_EXCL,
        ec=CLR_EXCL_BD,
        fontsize=5.5,
    )
    side_arrow(
        ax,
        cx + excl_sign * BOX_W / 2,
        y - 0.28,
        x_excl - excl_sign * EXCL_W / 2,
        y - 0.40,
    )

    # Row 3: Post-ESKD
    y = y_next
    add_box(ax, cx, y, BOX_W, BOX_H, f"Post-ESKD exclusion\nn = {nums['post_eskd']}")

    # Split into treated / control
    y_split = y - 0.72
    x_trt = cx - 0.75
    x_ctl = cx + 0.75
    bw_half = 1.2  # half-box width for split

    # Split arrows
    arrow(ax, cx, y - BOX_H / 2, x_trt, y_split + BOX_H / 2)
    arrow(ax, cx, y - BOX_H / 2, x_ctl, y_split + BOX_H / 2)

    # Row 4a: IV Mg group
    add_box(
        ax,
        x_trt,
        y_split,
        bw_half,
        BOX_H,
        f"IV Mg received\nn = {nums['ivmg']}\n({nums['ivmg_pct']})",
        fontsize=5.5,
    )

    # Row 4b: No IV Mg group
    add_box(
        ax,
        x_ctl,
        y_split,
        bw_half,
        BOX_H,
        f"No IV Mg\nn = {nums['no_ivmg']}",
        fontsize=5.5,
    )

    # --- Treated arm filters ---
    y_trt1 = y_split - 0.80
    arrow(ax, x_trt, y_split - BOX_H / 2, x_trt, y_trt1 + BOX_H_TALL / 2)

    add_box(
        ax,
        x_trt,
        y_trt1,
        bw_half,
        BOX_H_TALL,
        (
            f"≥1 Cr before T₀\n"
            f"ICU: {nums['cr_pre_icu']}\n"
            f"+ hosp fallback: {nums['hosp_fallback']}\n"
            f"Excl Cr ≥ 4.0: {nums['excl_cr_high']}"
        ),
        fontsize=5,
        fc="#FAFAFA",
    )

    # Treated final
    y_trt2 = y_trt1 - 0.72
    arrow(ax, x_trt, y_trt1 - BOX_H_TALL / 2, x_trt, y_trt2 + BOX_H / 2)
    add_box(
        ax,
        x_trt,
        y_trt2,
        bw_half,
        BOX_H,
        f"Treated final\nn = {nums['treated_final']}",
        fc=CLR_FINAL,
        ec=CLR_FINAL_BD,
        fontweight="bold",
        fontsize=6,
    )

    # --- Control arm filters ---
    y_ctl1 = y_split - 0.80
    arrow(ax, x_ctl, y_split - BOX_H / 2, x_ctl, y_ctl1 + BOX_H / 2)

    n_no_2cr = int(nums["no_ivmg"].replace(",", "")) - int(
        nums["ctrl_has_2cr"].replace(",", "")
    )
    n_ctrl_excl = int(nums["ctrl_has_2cr"].replace(",", "")) - int(
        nums["ctrl_final"].replace(",", "")
    )
    add_box(
        ax,
        x_ctl,
        y_ctl1,
        bw_half,
        BOX_H_TALL,
        (
            f"≥2 postop Cr: {nums['ctrl_has_2cr']}\n"
            f"< 2 Cr: {n_no_2cr:,}\n"
            f"Excl Cr ≥ 4.0: {n_ctrl_excl:,}"
        ),
        fontsize=5,
        fc="#FAFAFA",
    )

    # Control final
    y_ctl2 = y_ctl1 - 0.72
    arrow(ax, x_ctl, y_ctl1 - BOX_H_TALL / 2, x_ctl, y_ctl2 + BOX_H / 2)
    add_box(
        ax,
        x_ctl,
        y_ctl2,
        bw_half,
        BOX_H,
        f"Control final\nn = {nums['ctrl_final']}",
        fc=CLR_FINAL,
        ec=CLR_FINAL_BD,
        fontweight="bold",
        fontsize=6,
    )

    # --- Merge into eligible cohort ---
    y_elig = y_trt2 - 0.72
    arrow(ax, x_trt, y_trt2 - BOX_H / 2, cx, y_elig + BOX_H / 2)
    arrow(ax, x_ctl, y_ctl2 - BOX_H / 2, cx, y_elig + BOX_H / 2)
    add_box(
        ax,
        cx,
        y_elig,
        BOX_W,
        BOX_H,
        f"Eligible cohort\nn = {nums['eligible']}",
        fontsize=6.5,
        fontweight="bold",
    )

    # --- Risk-set PSM ---
    y_match = y_elig - 0.72
    arrow(ax, cx, y_elig - BOX_H / 2, cx, y_match + BOX_H / 2)
    add_box(
        ax,
        cx,
        y_match,
        BOX_W,
        BOX_H,
        f"Risk-set PSM (1:1)\n{nums['matched']} matched pairs\n({nums['match_pct']} matched, max SMD {nums['max_smd']})",
        fc=CLR_MATCH,
        ec=CLR_MATCH_BD,
        fontweight="bold",
        fontsize=6,
    )

    return y_match


# ── Build figure ──
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
ax.set_xlim(0, FIG_W)
ax.set_ylim(0, FIG_H)
ax.set_aspect("equal")
ax.axis("off")

# Title
ax.text(
    CX_MID,
    FIG_H - 0.15,
    "eFigure 1. Cohort Selection and Risk-Set Propensity Score Matching",
    fontsize=8.5,
    fontweight="bold",
    ha="center",
    va="top",
    color=CLR_HEADER,
)

# Column headers
ax.text(
    CX_E,
    FIG_H - 0.55,
    "eICU-CRD\n(208 US hospitals, 2014–2015)",
    fontsize=6.5,
    fontweight="bold",
    ha="center",
    va="top",
    color=CLR_HEADER,
)
ax.text(
    CX_M,
    FIG_H - 0.55,
    "MIMIC-IV\n(BIDMC, 2008–2019)",
    fontsize=6.5,
    fontweight="bold",
    ha="center",
    va="top",
    color=CLR_HEADER,
)

# Vertical separator
ax.axvline(CX_MID, ymin=0.02, ymax=0.92, color="#E0E0E0", lw=0.5, ls="--", zorder=0)

# Build both columns
# eICU: exclusion boxes go LEFT (negative side)
build_column(ax, CX_E, "eICU-CRD", E, x_excl_side=CX_E - 2)
# MIMIC: exclusion boxes go RIGHT (positive side)
build_column(ax, CX_M, "MIMIC-IV", M, x_excl_side=CX_M + 2)

# Footer
ax.text(
    CX_MID,
    0.3,
    (
        "T₀ = time of first IV magnesium administration (patient-specific).  "
        "Risk-set matching: 1:1 with replacement, caliper 0.2 SD, "
        "19 PS covariates, MICE m = 20."
    ),
    fontsize=5.5,
    ha="center",
    va="bottom",
    color="#666666",
    style="italic",
)

# Save
for ext in ("pdf", "png"):
    fig.savefig(
        f"/home/g91p721/mg_aki/results/efig1_consort.{ext}",
        format=ext,
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.08,
    )
plt.close(fig)
print("✓ Saved: efig1_consort.pdf / .png")
