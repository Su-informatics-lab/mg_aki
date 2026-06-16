#!/usr/bin/env bash
# ============================================================================
# run.sh — Full pipeline for mg_aki (v5: AC primary, threshold exploratory)
#
# Usage:
#   bash run.sh              # all steps (1-9 + fig)
#   bash run.sh 1 2 7        # subset
#   bash run.sh qc           # QC probes only
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
#   9   → 09_robustness.R       Lactate sensitivity + QBA
#   fig → gen_figures.py         Python figures (Nature style)
#   qc  → probe scripts          BMI audit, APACHE, achieved-Mg, weighted SMDs
# ============================================================================
set -euo pipefail
cd ~/mg_aki

if [ -f .venv/bin/activate ]; then source .venv/bin/activate; fi
if command -v module &>/dev/null; then
  module purge 2>/dev/null && module load R/4.4.2-gfbf-2024a 2>/dev/null || true
fi

SECONDS=0
LOGFILE="results/run_$(date +%Y%m%d_%H%M%S).log"
mkdir -p results figs

step() { echo -e "\n──── Step $1: $2 ────"; }

if [ $# -eq 0 ]; then
  STEPS=(1 2 3 4 5 6 7 8 9 fig)
else
  STEPS=("$@")
fi

{
  echo "Pipeline started: $(date)"
  echo "Steps: ${STEPS[*]}"
  echo ""

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
      9)   step 9   "Robustness";        Rscript 09_robustness.R ;;
      fig) step fig "Figures (Python)";  python gen_figures.py ;;
      qc)  step qc  "QC Probes"
           python probe_bmi_apache.py
           python probe_bmi_rootcause.py
           Rscript probe_v5_experiments.R
           ;;
      *)   echo "Unknown step: $s"
           echo "Valid: 1 2 3 4 5 6 7 8 9 fig qc"
           exit 1 ;;
    esac
  done

  echo -e "\n✓ Done ($(( SECONDS/60 ))m$(( SECONDS%60 ))s)"
  ls -lh results/*.csv figs/*.pdf 2>/dev/null || true

} 2>&1 | tee "$LOGFILE"

echo "Log: $LOGFILE"
