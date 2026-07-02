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
