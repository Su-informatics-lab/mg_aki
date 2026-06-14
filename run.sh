#!/usr/bin/env bash
set -euo pipefail
cd ~/mg_aki
source .venv/bin/activate
module purge && module load R/4.4.2-gfbf-2024a
SECONDS=0

step() { echo -e "\n──── Step $1/6: $2 ────"; }

STEPS=("${@:-1 2 3 4 5 6}")

for s in ${STEPS[@]}; do
  case $s in
    1) step 1 "eICU ETL";      python 01_etl.py ;;
    2) step 2 "eICU PS/IPTW";  Rscript 02_psm.R ;;
    3) step 3 "eICU Models";   Rscript 03_models.R ;;
    4) step 4 "MIMIC ETL";     python 04_mimic_validation.py ;;
    5) step 5 "MIMIC TTE";     Rscript 05_mimic_tte.R ;;
    6) step 6 "Meta-analysis"; Rscript 06_meta.R ;;
    *) echo "Unknown step: $s"; exit 1 ;;
  esac
done

echo -e "\n✓ Done ($(( SECONDS/60 ))m$(( SECONDS%60 ))s)"
ls -lh results/*.csv results/*.pdf 2>/dev/null
