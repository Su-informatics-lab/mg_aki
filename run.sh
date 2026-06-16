#!/usr/bin/env bash
# ============================================================================
# run.sh — Full pipeline for mg_aki
#
# Usage:
#   bash run.sh              # all steps
#   bash run.sh 1 2 7        # subset
#
# Steps:
#   1   → 01_etl.py             Cohort construction (eICU + MIMIC)
#   2   → 02_analysis.R         Primary analysis (MICE + OW/IPTW/PSM/AC)
#   3   → 03_augment.py         Tier 2 variables (vent, IABP, etc.)
#   4   → 04_subgroups.R        Subgroup + safety analyses
#   5   → 05_tables.R           Publication Tables 1–2
#   6   → 06_figures.R          R-based figures (CONSORT CSV)
#   7   → 07_sensitivity.R      E-values + prognostic + MICE + AC table1
#   8   → 08_stratified.R       Mg-stratified + hospital RE
#   fig → gen_figures.py        Python figures
# ============================================================================
set -euo pipefail
cd ~/mg_aki

if [ -f .venv/bin/activate ]; then source .venv/bin/activate; fi
if command -v module &>/dev/null; then
  module purge 2>/dev/null && module load R/4.4.2-gfbf-2024a 2>/dev/null || true
fi

SECONDS=0
mkdir -p results figs
step() { echo -e "\n──── Step $1: $2 ────"; }

if [ $# -eq 0 ]; then
  STEPS=(1 2 3 4 5 6 7 8 fig)
else
  STEPS=("$@")
fi

for s in "${STEPS[@]}"; do
  case $s in
    1)   step 1   "ETL (eICU + MIMIC)"; python 01_etl.py ;;
    2)   step 2   "Primary Analysis";   Rscript 02_analysis.R ;;
    3)   step 3   "Augment Cohorts";    python 03_augment.py ;;
    4)   step 4   "Subgroup/Safety";    Rscript 04_subgroups.R ;;
    5)   step 5   "Tables";            Rscript 05_tables.R ;;
    6)   step 6   "Figures (R)";       Rscript 06_figures.R ;;
    7)   step 7   "Sensitivity";       Rscript 07_sensitivity.R ;;
    8)   step 8   "Stratified + RE";   Rscript 08_stratified.R ;;
    fig) step fig "Figures (Python)";  python gen_figures.py ;;
    *)   echo "Unknown step: $s"; echo "Valid: 1 2 3 4 5 6 7 8 fig"; exit 1 ;;
  esac
done

echo -e "\n✓ Done ($(( SECONDS/60 ))m$(( SECONDS%60 ))s)"
ls -lh results/*.csv figs/*.pdf 2>/dev/null || true
