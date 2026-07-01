# mg_aki

Retrospective cohort study of postoperative IV magnesium and AKI after cardiac surgery, with pre-specified eGFR×treatment interaction.
Risk-set propensity score matching + difference-in-differences across MIMIC-IV and eICU-CRD.

## Data

Both require credentialed access via [PhysioNet](https://physionet.org):

- [MIMIC-IV v3.1](https://physionet.org/content/mimiciv/3.1/)
- [eICU-CRD v2.0](https://physionet.org/content/eicu-crd/2.0/)

## Dependencies

**R 4.5.1:** `mice`, `sandwich`, `lmtest`, `MatchIt`, `lme4`

**Python 3.10.8:** `pandas`, `numpy`, `matplotlib`

## Pipeline

Scripts are numbered in execution order. Each reads from `~/mg_aki/results/` and writes back to it.

| Script | What it does |
|--------|--------------|
| `00_config.py` | Shared constants: item IDs, ICD codes, lab ranges, file paths |
| `01_etl.py` | Cohort extraction from raw MIMIC-IV / eICU-CRD tables |
| `02_psm.R` | Risk-set propensity score matching (19 covariates, 1:1, MICE m=20). Primary + sensitivity (earliest-labs) specs |
| `03_hte.R` | Pre-specified subgroup ORs and interaction tests |
| `03b_egfr_aki_stages.R` | eGFR-stratified binary AKI stages (48h + 7d). Accepts optional `[spec]` argument for sensitivity |
| `03c_mg_strat.R` | eGFR × baseline Mg cross-stratification |
| `03d_mg_did.R` | Continuous ΔCr difference-in-differences at 6h intervals |
| `03e_pair_dcr.R` | Pair-level ΔCr export for cumulative incidence figures |
| `03f_hosp_re.R` | Hospital random-effects sensitivity (eICU-CRD only; requires `lme4`) |
| `04_figures.py` | All publication figures (forests, heatmaps, love plots, sensitivity) |
| `04c_fig_cuminc_egfr.py` | Cumulative AKI incidence by eGFR stratum (combined 48h+7d) |
| `05_qc.R` | Post-hoc balance and outcome QC checks |

### Generators (tables and supplementary materials)

| Script | Output |
|--------|--------|
| `gen_table1.py` | Table 1 (baseline characteristics) |
| `gen_etables.py` | eTables for supplement |
| `gen_mortality_table.R` | Mortality summary (Table 2 supplement) |
| `gen_consort.py` | CONSORT flow numbers |

### Probes (standalone, not in pipeline)

| Script | Purpose |
|--------|---------|
| `06_replicate_xiong_koh.py` | Reproduce Xiong 2023 Mg-AKI association in both databases |
| `probe_etl_qc.py` | ETL sanity checks (counts, distributions) |
| `probe_km_vs_binary.py` | KM vs naive cumulative incidence comparison |
| `probe_competing_risks.py` | Competing risk (death before AKI) assessment |

## Reproduction

```bash
# on Tempest — ETL (Python, run first)
module load Python/3.10.8-GCCcore-12.2.0
python 01_etl.py

# R steps
module purge
module load R/4.5.1-gfbf-2025a
Rscript 02_psm.R mimic
Rscript 02_psm.R eicu
Rscript 03_hte.R mimic
Rscript 03_hte.R eicu
Rscript 03b_egfr_aki_stages.R mimic
Rscript 03b_egfr_aki_stages.R eicu
Rscript 03c_mg_strat.R mimic
Rscript 03c_mg_strat.R eicu
Rscript 03d_mg_did.R mimic
Rscript 03d_mg_did.R eicu
Rscript 03e_pair_dcr.R mimic
Rscript 03e_pair_dcr.R eicu
Rscript 03f_hosp_re.R

# figures (Python)
module purge
module load Python/3.10.8-GCCcore-12.2.0
python 04_figures.py
python 04c_fig_cuminc_egfr.py
```


## License

MIT. Datasets governed by PhysioNet data use agreements.

## Citation

```bibtex
[PLACEHOLDER]
```

## Contact

- [Jing Su](mailto:su1@iu.edu) for general questions.
- [Haining Wang](mailto:hw56@iu.edu) for reproduction.

Su Lab in Biomedical Informatics, Biostatistics & Health Data Science · Indiana University School of Medicine
