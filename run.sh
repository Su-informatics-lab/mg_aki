#!/bin/bash
# ============================================================================
# run.sh — Full pipeline rerun for Mg → AKI study
#
# Usage:
#   bash run.sh          # full clean rerun (all steps)
#   bash run.sh 2 3 4    # only steps 2, 3, 4
#   bash run.sh clean     # just clean results/
#
# Prerequisites:
#   - Python 3.12+ with pandas, numpy, matplotlib, duckdb
#   - R 4.4+ with mice, sandwich, lmtest
#   - Raw data at ~/mg_aki/eicu-crd-2.0/ and ~/mg_aki/mimic-iv-3.1/
# ============================================================================
set -euo pipefail

cd "$(dirname "$0")"
RESULTS="results"
SEP="======================================================================"

log() { echo -e "\n$SEP\n  STEP $1: $2\n$SEP"; }
elapsed() { echo "  ⏱  $(( SECONDS - $1 ))s"; }

# ── Parse args ──
STEPS=("$@")
RUN_ALL=false
if [ ${#STEPS[@]} -eq 0 ]; then
    RUN_ALL=true
    STEPS=(clean 1 2 3 4 5 6 7)
fi

for step in "${STEPS[@]}"; do
case $step in

# ── CLEAN ──────────────────────────────────────────────────────────
clean)
    log "0" "Clean results/"
    mkdir -p $RESULTS
    rm -f $RESULTS/did_*.csv
    rm -f $RESULTS/table1_*.csv
    rm -f $RESULTS/etable*.csv
    rm -f $RESULTS/km_aki_*.csv
    rm -f $RESULTS/fig*.pdf $RESULTS/fig*.png
    rm -f $RESULTS/efig*.pdf $RESULTS/efig*.png
    echo "  ✓ Cleaned"
    ;;

# ── STEP 1: ETL ───────────────────────────────────────────────────
1)
    log "1" "ETL (01_etl.py) — cohort construction"
    t0=$SECONDS
    python 01_etl.py
    elapsed $t0
    echo "  Outputs: did_all_{db}.csv, did_cr_all_{db}.csv, did_labs_all_{db}.csv"
    ;;

# ── STEP 2: PSM ───────────────────────────────────────────────────
2)
    log "2" "PSM (02_psm.R) — risk-set matching [set.seed(2026)]"
    t0=$SECONDS
    Rscript 02_psm.R mimic
    Rscript 02_psm.R eicu
    elapsed $t0
    echo "  Outputs: did_riskset_{db}.csv, did_pairs_*_{db}.csv"
    ;;

# ── STEP 3: HTE ───────────────────────────────────────────────────
3)
    log "3" "HTE (03_hte.R) — subgroup analysis"
    t0=$SECONDS
    Rscript 03_hte.R mimic
    Rscript 03_hte.R eicu
    elapsed $t0
    echo "  Outputs: did_hte_{db}.csv, did_hte_data_{db}.csv"
    ;;

# ── STEP 4: FIGURES ────────────────────────────────────────────────
4)
    log "4" "Figures (04_figures.py) — 7 publication figures"
    t0=$SECONDS
    python 04_figures.py
    elapsed $t0
    echo "  Outputs: fig1_primary, fig2_hte, fig3_benefit_harm,"
    echo "           efig_timecourse, efig_sensitivity, efig_hte_48h, efig_love_plot"
    ;;

# ── STEP 5: TABLE 1 ───────────────────────────────────────────────
5)
    log "5" "Table 1 (gen_table1.py)"
    t0=$SECONDS
    python gen_table1.py
    elapsed $t0
    echo "  Outputs: table1_{db}.csv, table1_combined.csv"
    ;;

# ── STEP 6: eTABLES ───────────────────────────────────────────────
6)
    log "6" "eTables (gen_etables.py)"
    t0=$SECONDS
    python gen_etables.py
    elapsed $t0
    echo "  Outputs: etable2_balance.csv, etable3_hte_matrix.csv"
    ;;

# ── STEP 7: CONSORT ───────────────────────────────────────────────
7)
    log "7" "CONSORT (gen_consort.py)"
    t0=$SECONDS
    python gen_consort.py
    elapsed $t0
    echo "  Outputs: efig1_consort.pdf/.png"
    ;;

*)
    echo "  Unknown step: $step (valid: clean, 1-7)"
    ;;

esac
done

echo -e "\n$SEP"
echo "  DONE — $(ls $RESULTS/*.pdf 2>/dev/null | wc -l) PDFs in $RESULTS/"
echo "$SEP"
