#!/usr/bin/env python3
"""
Generate Figure 2: Forest plot for mg_aki manuscript v2.
- AKI outcomes (4 rows)
- Active comparator (1 row)
- Hospital mortality exploratory (1 row)
- Control outcomes (2 rows)
- Significance annotations (* for P<.05)
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Data ──────────────────────────────────────────────────────────────
# Each entry: (label, eicu_or, eicu_lo, eicu_hi, eicu_p,
#              mimic_or, mimic_lo, mimic_hi, mimic_p,
#              pool_or, pool_lo, pool_hi, pool_p)

sections = [
    (
        "AKI Outcomes",
        [
            (
                "KDIGO stage ≥1 (primary)",
                0.77,
                0.61,
                0.98,
                0.032,
                0.91,
                0.74,
                1.13,
                0.40,
                0.85,
                0.72,
                0.99,
                0.040,
            ),
            (
                "Creatinine ratio ≥1.5×",
                0.80,
                0.61,
                1.05,
                0.11,
                0.90,
                0.71,
                1.13,
                0.36,
                0.85,
                0.72,
                1.02,
                0.078,
            ),
            (
                "KDIGO stage ≥2",
                0.95,
                0.73,
                1.25,
                0.73,
                0.87,
                0.63,
                1.20,
                0.39,
                0.92,
                0.75,
                1.13,
                0.42,
            ),
            (
                "AKI within 48 h",
                0.68,
                0.46,
                0.99,
                0.045,
                0.95,
                0.70,
                1.27,
                0.71,
                0.83,
                0.66,
                1.05,
                0.12,
            ),
        ],
    ),
    (
        "Active Comparator (Mg+K⁺ vs K⁺-only)",
        [
            (
                "KDIGO stage ≥1",
                0.71,
                0.48,
                1.05,
                0.083,
                0.94,
                0.68,
                1.30,
                0.71,
                0.84,
                0.66,
                1.08,
                0.168,
            ),
        ],
    ),
    (
        "Mortality (exploratory)",
        [
            (
                "Hospital mortality",
                0.94,
                0.64,
                1.37,
                0.73,
                0.65,
                0.45,
                0.95,
                0.025,
                0.78,
                0.60,
                1.02,
                0.069,
            ),
        ],
    ),
    (
        "Control Outcomes",
        [
            (
                "Fracture (negative control)",
                0.107e1,
                0.55,
                2.09,
                0.83,
                1.01,
                0.56,
                1.80,
                0.99,
                1.03,
                0.67,
                1.60,
                0.88,
            ),
            (
                "Encephalopathy (exploratory)",
                0.72,
                0.41,
                1.26,
                0.25,
                0.47,
                0.26,
                0.86,
                0.014,
                0.59,
                0.39,
                0.89,
                0.012,
            ),
        ],
    ),
]

# Fix fracture eICU value
sections[3][1][0] = (
    "Fracture (negative control)",
    1.07,
    0.55,
    2.09,
    0.83,
    1.01,
    0.56,
    1.80,
    0.99,
    1.03,
    0.67,
    1.60,
    0.88,
)

# ── Layout ────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(12, 8.5))

y_pos = 0
y_ticks = []
y_labels = []
spacing = 1.0
section_gap = 0.6
header_gap = 0.3

eicu_color = "#2166AC"
mimic_color = "#B2182B"
pool_color = "#333333"

# Process bottom-up so y=0 is bottom
all_rows = []
for sec_name, rows in reversed(sections):
    for row in reversed(rows):
        all_rows.append(("row", row))
    all_rows.append(("header", sec_name))
    all_rows.append(("gap", None))

# Remove trailing gap
if all_rows and all_rows[-1][0] == "gap":
    all_rows = all_rows[:-1]

y = 0
plot_data = []

for item_type, item in all_rows:
    if item_type == "gap":
        y += section_gap
    elif item_type == "header":
        plot_data.append(("header", item, y))
        y += header_gap
    elif item_type == "row":
        plot_data.append(("row", item, y))
        y += spacing

# ── Draw ──────────────────────────────────────────────────────────────
ax.axvline(x=1.0, color="#999999", linestyle="--", linewidth=0.8, zorder=0)

for ptype, pdata, py in plot_data:
    if ptype == "header":
        ax.text(
            0.25,
            py,
            pdata,
            fontsize=10.5,
            fontweight="bold",
            fontfamily="sans-serif",
            va="center",
            transform=ax.get_yaxis_transform(),
        )
    elif ptype == "row":
        label, e_or, e_lo, e_hi, e_p, m_or, m_lo, m_hi, m_p, p_or, p_lo, p_hi, p_p = (
            pdata
        )

        # eICU square
        ax.plot(e_or, py + 0.22, "s", color=eicu_color, markersize=7, zorder=3)
        ax.plot(
            [e_lo, e_hi],
            [py + 0.22, py + 0.22],
            "-",
            color=eicu_color,
            linewidth=1.2,
            zorder=2,
        )

        # MIMIC square
        ax.plot(m_or, py, "s", color=mimic_color, markersize=7, zorder=3)
        ax.plot([m_lo, m_hi], [py, py], "-", color=mimic_color, linewidth=1.2, zorder=2)

        # Pooled diamond
        diamond_h = 0.12
        diamond_w = (p_hi - p_lo) * 0.05  # visual width
        ax.plot(p_or, py - 0.22, "D", color=pool_color, markersize=8, zorder=3)
        ax.plot(
            [p_lo, p_hi],
            [py - 0.22, py - 0.22],
            "-",
            color=pool_color,
            linewidth=1.5,
            zorder=2,
        )

        # Significance star
        sig_marker = ""
        if p_p < 0.05:
            sig_marker = "*"

        # Labels on left
        indent = "    "
        ax.text(
            0.02,
            py,
            f"{indent}{label}",
            fontsize=9.5,
            va="center",
            fontfamily="sans-serif",
            transform=ax.get_yaxis_transform(),
        )

        # Values on right — format with sig annotation
        def fmt(or_val, lo, hi, p, color):
            s = f"{or_val:.2f} ({lo:.2f}–{hi:.2f})"
            if p < 0.05:
                s += "*"
            return s

        e_str = fmt(e_or, e_lo, e_hi, e_p, eicu_color)
        m_str = fmt(m_or, m_lo, m_hi, m_p, mimic_color)
        p_str = fmt(p_or, p_lo, p_hi, p_p, pool_color)

        # Right-side annotation column — use axes x-coord so labels don't clip
        import matplotlib.transforms as mtrans

        trans_r = mtrans.blended_transform_factory(ax.transAxes, ax.transData)
        rx = 0.83
        ax.text(
            rx,
            py + 0.22,
            e_str,
            fontsize=8,
            va="center",
            color=eicu_color,
            fontfamily="monospace",
            fontweight="bold" if e_p < 0.05 else "normal",
            transform=trans_r,
            clip_on=False,
        )
        ax.text(
            rx,
            py,
            m_str,
            fontsize=8,
            va="center",
            color=mimic_color,
            fontfamily="monospace",
            fontweight="bold" if m_p < 0.05 else "normal",
            transform=trans_r,
            clip_on=False,
        )
        ax.text(
            rx,
            py - 0.22,
            p_str,
            fontsize=8,
            va="center",
            color=pool_color,
            fontfamily="monospace",
            fontweight="bold" if p_p < 0.05 else "normal",
            transform=trans_r,
            clip_on=False,
        )

# ── Axes ──────────────────────────────────────────────────────────────
ax.set_xlim(0.35, 2.15)
ax.set_xscale("log")
from matplotlib.ticker import FixedFormatter, FixedLocator

ax.xaxis.set_major_locator(FixedLocator([0.50, 0.75, 1.00, 1.50, 2.00]))
ax.xaxis.set_major_formatter(FixedFormatter(["0.50", "0.75", "1.00", "1.50", "2.00"]))
ax.xaxis.set_minor_locator(FixedLocator([]))  # kill minor ticks
ax.tick_params(axis="x", labelsize=10)
ax.set_xlabel("Odds Ratio (95% CI)", fontsize=11, fontfamily="sans-serif")

max_y = max(py for _, _, py in plot_data) + 1.0
ax.set_ylim(-0.8, max_y)
ax.set_yticks([])

# Favors labels
ax.text(
    0.68,
    -0.55,
    "← Favors supplementation",
    fontsize=9,
    ha="center",
    color="#555555",
    fontfamily="sans-serif",
)
ax.text(
    1.45,
    -0.55,
    "Favors no supplementation →",
    fontsize=9,
    ha="center",
    color="#555555",
    fontfamily="sans-serif",
)

# Legend
from matplotlib.lines import Line2D

legend_elements = [
    Line2D(
        [0],
        [0],
        marker="s",
        color="w",
        markerfacecolor=eicu_color,
        markersize=8,
        label="eICU-CRD",
    ),
    Line2D(
        [0],
        [0],
        marker="s",
        color="w",
        markerfacecolor=mimic_color,
        markersize=8,
        label="MIMIC-IV",
    ),
    Line2D(
        [0],
        [0],
        marker="D",
        color="w",
        markerfacecolor=pool_color,
        markersize=8,
        label="Pooled",
    ),
]
ax.legend(
    handles=legend_elements,
    loc="upper right",
    fontsize=9,
    framealpha=0.9,
    edgecolor="#cccccc",
)

# Title
ax.set_title(
    "Figure 2. Pooled Odds Ratios for AKI, Mortality, and Control Outcomes",
    fontsize=12,
    fontweight="bold",
    fontfamily="sans-serif",
    pad=12,
)

# Footnote
fig.text(
    0.12,
    0.01,
    "* P < .05. Bold values indicate statistical significance.",
    fontsize=8,
    color="#666666",
    fontfamily="sans-serif",
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)

plt.tight_layout(rect=[0, 0.03, 0.82, 0.97])
plt.savefig("/home/claude/mg_aki/fig2_forest.pdf", dpi=300, bbox_inches="tight")
plt.savefig("/home/claude/mg_aki/fig2_forest.png", dpi=300, bbox_inches="tight")
print("Saved fig2_forest.pdf and .png")
