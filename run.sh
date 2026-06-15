#!/usr/bin/env bash
set -euo pipefail
cd ~/mg_aki
source .venv/bin/activate
module purge && module load R/4.4.2-gfbf-2024a
SECONDS=0

step() { echo -e "\n──── Step $1: $2 ────"; }

STEPS=("${@:-1 2 3 4 5 6}")

for s in ${STEPS[@]}; do
  case $s in
    1)  step 1 "eICU ETL";          python 01_etl.py ;;
    1b) step 1b "MIMIC ETL";        python 01b_mimic_etl.py ;;
    2)  step 2 "Primary Analysis";   Rscript 02_analysis.R ;;
    3)  step 3 "Augment Cohorts";    python 03_augment.py ;;
    4)  step 4 "Subgroup/Safety";    Rscript 04_subgroups.R ;;
    5)  step 5 "Tables";            Rscript 05_tables.R ;;
    6)  step 6 "Figures";           Rscript 06_figures.R ;;
    7)  step 7 "E-values";          python 07_evalues.py ;;
    *)  echo "Unknown step: $s (valid: 1 1b 2 3 4 5 6 7)"; exit 1 ;;
  esac
done

echo -e "\n✓ Done ($(( SECONDS/60 ))m$(( SECONDS%60 ))s)"
ls -lh results/*.csv results/*.pdf 2>/dev/null
