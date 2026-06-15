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
A full pipeline log is included at `results/run_20260614_123038.log`.

bash run.sh        # full pipeline (steps 1-6)
bash run.sh 4 5 6  # subset
```

| Step | Script | Output |
|------|--------|--------|
| 1 | `01_etl.py` | `results/01_analysis_a_cohort.csv` (incl. K⁺ supp, death timing) |
| 2 | `02_psm.R` | `results/02e_all_iptw.csv`, `results/02g_ac_iptw.csv` |
| 3 | `03_models.R` | `results/03_results_summary.csv` (IPTW + AC + landmark + overlap) |
| 4 | `04_mimic_validation.py` | `results/04_mimic_cohort.csv` (incl. K⁺ supp, death timing) |
| 5 | `05_mimic_tte.R` | `results/05_mimic_results_summary.csv` (IPTW + AC + landmark) |
| 6 | `06_meta.R` | `results/06_meta_results.csv` (all-patient + active comparator) |

All phenotype definitions, ICD codes, and lab ranges are in `00_config.py`.

## Key results

Across 11,670 cardiac surgery patients (eICU: 8,109; MIMIC-IV: 3,746):

| Analysis | Outcome | eICU OR | Pooled OR (95% CI) | *P* |
|---|---|---|---|---|
| All-patient IPTW | AKI KDIGO ≥1 | 0.78 | 0.84 (0.73–0.96) | .012 |
| Active comparator (Mg+K⁺ vs K⁺) | AKI KDIGO ≥1 | 0.71 | 0.79 (0.65–0.96) | .018 |
| All-patient IPTW | Mortality (exploratory) | 0.89 | — | NS |
| All-patient IPTW | Fracture (negative control) | 0.97 | 1.03 (0.67–1.60) | .88 |

The previously reported association between higher serum magnesium and
worse AKI (Xiong et al. 2023) reflects confounding by cardioplegia volume,
not magnesium harm.

## License

MIT. Datasets governed by PhysioNet data use agreements.
