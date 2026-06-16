#!/usr/bin/env bash
# ============================================================================
# run.sh — Full pipeline for mg_aki project (v4: threshold narrative)
#
# Usage:
#   bash run.sh              # all steps
#   bash run.sh 1 1b 2       # subset
#   bash run.sh 8 8b fig     # just new analyses + figures
#
# Steps:
#   1   → 01_etl.py             eICU cohort construction
#   1b  → 01b_mimic_etl.py      MIMIC-IV cohort construction
#   8e  → 08e_fix_mimic_poaf.py Fix MIMIC POAF phenotype (run before 2)
#   2   → 02_analysis.R         Primary analysis (MICE + OW/IPTW/PSM/AC)
#   3   → 03_augment.py         Tier 2 variables (vent, IABP, etc.)
#   4   → 04_subgroups.R        Subgroup + safety analyses
#   5   → 05_tables.R           Publication Tables 1–2
#   6   → 06_figures.R          R-based figures (CONSORT CSV)
#   7   → 07_evalues.py         E-value sensitivity analysis
#   7b  → 07b_prognostic.R      Prognostic Mg→AKI by severity
#   7c  → 07c_mice_stability.R  MICE m=10,20 stability
#   7d  → 07d_ac_table1.R       AC baseline characteristics
#   8   → 08_mg_stratified.R    Mg-stratified treatment effects
#   8b  → 08b_hospital_re.R     Hospital random effects
#   fig → gen_figures.py         Python figures (Nature style)
# ============================================================================
set -euo pipefail
cd ~/mg_aki

# Activate Python venv
if [ -f .venv/bin/activate ]; then
  source .venv/bin/activate
fi

# Load R module (Tempest HPC)
if command -v module &>/dev/null; then
  module purge 2>/dev/null && module load R/4.4.2-gfbf-2024a 2>/dev/null || true
fi

SECONDS=0
LOGDIR="results"
mkdir -p "$LOGDIR"

step() { echo -e "\n──── Step $1: $2 ────"; }

# Default: all steps in order
if [ $# -eq 0 ]; then
  STEPS=(1 1b 8e 2 3 4 5 6 7 7b 7c 7d 8 8b fig)
else
  STEPS=("$@")
fi

for s in "${STEPS[@]}"; do
  case $s in
    1)   step 1   "eICU ETL";             python 01_etl.py ;;
    1b)  step 1b  "MIMIC ETL";            python 01b_mimic_etl.py ;;
    8e)  step 8e  "Fix MIMIC POAF";       python 08e_fix_mimic_poaf.py ;;
    2)   step 2   "Primary Analysis";     Rscript 02_analysis.R ;;
    3)   step 3   "Augment Cohorts";      python 03_augment.py ;;
    4)   step 4   "Subgroup/Safety";      Rscript 04_subgroups.R ;;
    5)   step 5   "Tables";              Rscript 05_tables.R ;;
    6)   step 6   "Figures (R)";         Rscript 06_figures.R ;;
    7)   step 7   "E-values";            python 07_evalues.py ;;
    7b)  step 7b  "Prognostic";          Rscript 07b_prognostic.R ;;
    7c)  step 7c  "MICE stability";      Rscript 07c_mice_stability.R ;;
    7d)  step 7d  "AC Table 1";          Rscript 07d_ac_table1.R ;;
    8)   step 8   "Mg-stratified";       Rscript 08_mg_stratified.R ;;
    8b)  step 8b  "Hospital RE";         Rscript 08b_hospital_re.R ;;
    fig) step fig "Figures (Python)";    python gen_figures.py ;;
    *)   echo "Unknown step: $s"; echo "Valid: 1 1b 8e 2 3 4 5 6 7 7b 7c 7d 8 8b fig"; exit 1 ;;
  esac
done

echo -e "\n✓ Done ($(( SECONDS/60 ))m$(( SECONDS%60 ))s)"
ls -lh results/*.csv figs/*.pdf 2>/dev/null || true
