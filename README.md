# Postoperative Magnesium Supplementation and AKI After Cardiac Surgery

Target trial emulation evaluating the association between postoperative
IV magnesium supplementation and AKI after cardiac surgery, with
serum magnesium–stratified analysis revealing a pharmacological
threshold. eICU-CRD (208 US hospitals, primary) and MIMIC-IV (BIDMC,
validation).

## Key finding

Magnesium supplementation was associated with reduced AKI only when
serum Mg exceeded ~2.3 mg/dL (OR 0.53, P=.003; interaction P=.005),
with the strongest signal at 2.6–3.0 mg/dL (AKI 12% vs 26%, OR 0.35).
Standard supplementation doses may be insufficient to reach the
renal-protective threshold.

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
bash run.sh
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
| fig | `gen_figures.py` | `figs/fig2_forest.pdf`, `figs/fig_mg_stratified.pdf`, `figs/fig3_subgroups.pdf` |

All phenotype definitions, ICD codes, and lab ranges are in
`00_config.py`.

## Key results

Across 11,855 cardiac surgery patients (eICU: 8,109; MIMIC-IV: 3,746):

| Analysis | eICU OR (95% CI) | P |
|---|---|---|
| Active comparator (Mg+K⁺ vs K⁺-only) | 0.75 (0.58–0.96) | .02 |
| All-patient IPTW | 0.76 (0.61–0.96) | .02 |
| Mg-stratified >2.3 mg/dL (OW) | 0.53 (0.35–0.80) | .003 |
| Sub-band 2.6–3.0 mg/dL | 0.35 (0.16–0.76) | .008 |
| Fracture negative control | 0.98 (0.55–1.76) | .95 |


## License

MIT. Datasets governed by PhysioNet data use agreements.
