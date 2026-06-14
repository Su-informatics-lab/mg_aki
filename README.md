# Postoperative Magnesium Supplementation and AKI After Cardiac Surgery

Target trial emulation of postoperative IV magnesium supplementation
on AKI and mortality after cardiac surgery. eICU-CRD (208 US hospitals,
primary) and MIMIC-IV (BIDMC, validation), pooled via fixed-effects
meta-analysis.

## Data

Both require credentialed access via [PhysioNet](https://physionet.org):

- [eICU-CRD v2.0](https://physionet.org/content/eicu-crd/2.0/)
- [MIMIC-IV v3.1](https://physionet.org/content/mimiciv/3.1/)

## Dependencies

**Python 3.12.3:** `pandas`, `numpy`

**R 4.4.2:** `tidyverse 2.0.0`, `survey 4.5`, `WeightIt 1.7.0`,
`cobalt 4.6.3`, `MatchIt 4.7.2`, `survival 3.8-6`, `broom 1.0.13`,
`sandwich 3.1-1`, `lmtest 0.9-40`, `tableone 0.13.2`, `EValue 4.1.4`

## Setup

1. Edit data paths in `00_config.py` and at the top of
   `04_mimic_validation.py`.

2. Install dependencies:

   ```bash
   python -m venv .venv && source .venv/bin/activate
   pip install pandas numpy
   Rscript -e 'install.packages(c("tidyverse","survey","broom","sandwich",
     "lmtest","WeightIt","cobalt","MatchIt","survival","tableone","EValue"))'
   ```

## Run

```bash
bash run.sh        # full pipeline (steps 1-6)
bash run.sh 4 5 6  # subset
```

| Step | Script | Output |
|------|--------|--------|
| 1 | `01_etl.py` | `results/01_analysis_a_cohort.csv` |
| 2 | `02_psm.R` | `results/02e_all_iptw.csv` |
| 3 | `03_models.R` | `results/03_results_summary.csv` |
| 4 | `04_mimic_validation.py` | `results/04_mimic_cohort.csv` |
| 5 | `05_mimic_tte.R` | `results/05_mimic_results_summary.csv` |
| 6 | `06_meta.R` | `results/06_meta_results.csv` |

All phenotype definitions, ICD codes, and lab ranges are in `00_config.py`.

## Key results

Across 11,670 cardiac surgery patients (eICU: 8,109; MIMIC-IV: 3,746):

| Outcome | Pooled OR (95% CI) | *P* | I² |
|---|---|---|---|
| AKI (KDIGO ≥ 1) | 0.85 (0.72–0.99) | .037 | 7% |
| Hospital mortality | 0.66 (0.51–0.85) | .002 | 0% |
| Fracture (negative control) | 1.04 (0.67–1.61) | .872 | 0% |

## License

MIT. Datasets governed by PhysioNet data use agreements.
