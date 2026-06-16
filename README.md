# Postoperative Magnesium Supplementation and AKI After Cardiac Surgery

Target trial emulation evaluating the association between postoperative
IV magnesium supplementation and AKI after cardiac surgery.
Active-comparator design (Mg+K⁺ vs K⁺-only) as primary analysis,
with exploratory serum magnesium–stratified analysis.
eICU-CRD (primary) and MIMIC-IV (replication).

## Key finding

In the active-comparator analysis, magnesium-plus-potassium repletion
was associated with lower AKI than potassium-only repletion (eICU
OR 0.75, P=.02; MIMIC-IV OR 0.93, P=.64). A complexity-specific
negative control (perioperative transfusion) validated the
active-comparator design: confounded in the all-patient comparison
(OR 0.76, P=.009), null within the AC framework (OR 0.95, P=.66).

Exploratory analysis revealed effect modification by baseline serum
magnesium (interaction P=.005), with the association concentrated
above 2.3 mg/dL (OR 0.53, P=.003). However, baseline serum magnesium
is not achieved serum magnesium — only 25% of supplemented patients
achieved follow-up Mg >2.3 mg/dL. The threshold observation is
hypothesis-generating and warrants a randomized trial targeting
achieved serum magnesium levels.

## Data

Both require credentialed access via [PhysioNet](https://physionet.org):

- [eICU-CRD v2.0](https://physionet.org/content/eicu-crd/2.0/)
- [MIMIC-IV v3.1](https://physionet.org/content/mimiciv/3.1/)

## Dependencies

**Python 3.12:** `pandas`, `numpy`, `matplotlib`

**R 4.4:** `mice`, `sandwich`, `lmtest`, `MatchIt`, `lme4`,
`tableone`, `ggplot2`, `gridExtra`

```bash
pip install pandas numpy matplotlib
Rscript -e 'install.packages(c("mice","sandwich","lmtest","MatchIt",
  "lme4","tableone","ggplot2","gridExtra"))'
```

## Setup

Edit data paths in `00_config.py` (eICU) and the `MIMIC_ROOT`
constant at the top of `01_etl.py` (MIMIC-IV).

## Pipeline

```bash
bash run.sh          # all steps
bash run.sh 1 2 7    # subset
```

| Step | Script | Output |
|------|--------|--------|
| 1 | `01_etl.py` | `results/01_analysis_a_cohort.csv`, `results/04_mimic_cohort.csv` |
| 2 | `02_analysis.R` | `results/02_results.csv` |
| 3 | `03_augment.py` | `*_enriched.csv` |
| 4 | `04_subgroups.R` | `etables_4_5_subgroups.csv`, `etable6_safety.csv` |
| 5 | `05_tables.R` | `05_table1.csv`, `05_table2.csv` |
| 6 | `06_figures.R` | `06_fig*.pdf`, `06_fig1_consort.csv` |
| 7 | `07_sensitivity.R` | `07_evalues.csv`, `07b_prognostic.csv`, `07c_mice_stability.csv`, `07d_ac_table1.csv` |
| 8 | `08_stratified.R` | `08_mg_stratified.csv`, `08b_hospital_re.csv` |
| 9 | `09_robustness.R` | `09_lactate_sensitivity.csv`, `09_qba.csv` |
| fig | `gen_figures.py` | `figs/fig2_forest.pdf`, `figs/fig_mg_stratified.pdf`, `figs/fig3_subgroups.pdf` |

All phenotype definitions, ICD codes, and lab ranges are in
`00_config.py`.

## Standalone probes (not in pipeline)

| Script | Purpose |
|--------|---------|
| `probe_complexity_nc.R` | Transfusion negative control (results in gen_figures.py) |
| `probe_bmi_apache.py` | BMI QC audit + APACHE IV extraction + achieved-Mg distribution |
| `probe_bmi_rootcause.py` | BMI outlier root cause classification (height/weight unit errors) |
| `probe_v5_experiments.R` | BMI sensitivity, APACHE IV PS probe, weighted SMD verification |

## QC notes

- **BMI**: eICU raw data contains 22 outliers (0.27%) from height
  recorded in inches and data entry errors across 16 hospitals.
  Capped at [10, 80] in ETL (matches MIMIC). Probe confirmed
  |ΔOR| = 0.0009 — zero practical effect on results.
- **APACHE IV**: 87% availability in eICU. Adding to PS model
  changed no estimate by more than |ΔOR| = 0.009.
- **Weighted SMDs**: All OW-weighted SMDs = 0.0000 in both databases
  for the AC analysis (verified by probe_v5_experiments.R).

## License

MIT. Datasets governed by PhysioNet data use agreements.
