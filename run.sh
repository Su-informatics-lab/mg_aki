#!/bin/bash
# run.sh — Full pipeline for Mg → AKI (Critical Care submission)
#
# Usage:
#   bash run.sh              # full run, both databases
#   bash run.sh 2 3          # steps 2 and 3 only
#   bash run.sh clean        # just clean results/
set -euo pipefail
cd "$(dirname "$0")"

RESULTS="results"
SEP="======================================================================"
log() { echo -e "\n$SEP\n  STEP $1: $2\n$SEP"; }
elapsed() { echo "  ⏱  $(( SECONDS - $1 ))s"; }

STEPS=("$@")
if [ ${#STEPS[@]} -eq 0 ]; then STEPS=(clean 1 2 3 4 5 6 7 8); fi

for step in "${STEPS[@]}"; do
case $step in

clean)
  log "0" "Clean results/"
  mkdir -p $RESULTS
  rm -f $RESULTS/did_*.csv $RESULTS/table1_*.csv $RESULTS/etable*.csv
  rm -f $RESULTS/egfr_*.csv $RESULTS/mg_strat_*.csv $RESULTS/mg_did_*.csv
  rm -f $RESULTS/fig*.pdf $RESULTS/fig*.png $RESULTS/efig*.pdf $RESULTS/efig*.png
  echo "  ✓ Cleaned"
  ;;

1)
  log "1" "ETL — cohort + RRT + death_offset + Table1 labs"
  t0=$SECONDS; python 01_etl.py; elapsed $t0
  ;;

2)
  log "2" "PSM — risk-set matching + DiD + KDIGO staging + binary outcomes"
  t0=$SECONDS
  Rscript 02_psm.R mimic
  Rscript 02_psm.R eicu
  elapsed $t0
  ;;

3)
  log "3" "HTE — pre-specified subgroups"
  t0=$SECONDS
  Rscript 03_hte.R mimic
  Rscript 03_hte.R eicu
  elapsed $t0
  ;;

4)
  log "4" "eGFR × AKI staging (KDIGO ≥1/≥2/≥3)"
  t0=$SECONDS
  Rscript 03b_egfr_aki_stages.R mimic
  Rscript 03b_egfr_aki_stages.R eicu
  elapsed $t0
  ;;

5)
  log "5" "Mg-stratified analysis + eGFR×Mg cross-strat"
  t0=$SECONDS
  Rscript 03c_mg_strat.R mimic
  Rscript 03c_mg_strat.R eicu
  elapsed $t0
  ;;

6)
  log "6" "Figures"
  t0=$SECONDS
  python 04_figures.py
  python 04b_fig_egfr.py
  python 04c_fig_km_egfr.py
  python 04d_fig_secondary_egfr.py
  elapsed $t0
  ;;

7)
  log "7" "Tables + CONSORT"
  t0=$SECONDS
  python gen_table1.py
  python gen_etables.py
  python gen_consort.py
  elapsed $t0
  ;;

8)
  log "8" "Supplement (optional, errors non-fatal)"
  for s in 03d_mg_did.R 03e_pair_dcr.R gen_mortality_table.R; do
    for d in mimic eicu; do
      Rscript "$s" "$d" 2>/dev/null || echo "  $s $d skipped"
    done
  done
  ;;

*) echo "  Unknown step: $step (valid: clean, 1-8)" ;;
esac
done

echo -e "\n$SEP"
echo "  DONE"
echo "  CSVs: $(ls $RESULTS/*.csv 2>/dev/null | wc -l)"
echo "  PDFs: $(ls $RESULTS/*.pdf 2>/dev/null | wc -l)"
echo "  Key: did_binary_*.csv (KDIGO ORs), did_riskset_*.csv (DiD)"
echo "  NEXT: git tag v6.0-cc-submission"
echo "$SEP"
