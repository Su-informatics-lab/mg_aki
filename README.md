# mg_aki

Target trial emulation: postoperative IV magnesium supplementation
and AKI/mortality after cardiac surgery. eICU-CRD (primary, 208
hospitals) + MIMIC-IV (validation). Targeting JAMA Network Open.

## Data

Requires PhysioNet credentialed access (username: hw56).

```
~/mg_aki/eicu-crd-2.0/     # eICU-CRD v2.0 CSVs
~/mg_aki/mimic-iv-3.1/     # MIMIC-IV v3.1 CSVs
```

## Environment

IU Tempest HPC, user g91p721.

```bash
module purge && module load R/4.4.2-gfbf-2024a
python -m venv .venv
source .venv/bin/activate
pip install pandas numpy
```

### R packages (~/R/libs)

```r
install.packages(c(
  "tidyverse", "survey", "broom", "sandwich", "lmtest",
  "WeightIt", "cobalt", "MatchIt", "survival", "tableone",
  "EValue", "metafor"
), lib = "~/R/libs")
```

Tested versions: R 4.4.2, Python 3.12.3, pandas 2.x, numpy 1.x.

## Pipeline

```
00_config.py              constants, ICD codes, lab ranges
01_etl.py                 eICU ETL → results/01_analysis_a_cohort.csv
02_psm.R                  eICU PS estimation, IPTW, matching
03_models.R               eICU outcomes → results/03_results_summary.csv
04_mimic_validation.py    MIMIC-IV ETL → results/04_mimic_cohort.csv
05_mimic_tte.R            MIMIC outcomes → results/05_mimic_results_summary.csv
06_meta.R                 pooled meta-analysis → results/06_meta_results.csv
```

## Run

```bash
bash run.sh        # all steps
bash run.sh 4 5 6  # MIMIC + meta only
```

## Results

| Outcome | Pooled OR (95% CI) | p | I² |
|---|---|---|---|
| KDIGO ≥1 AKI | 0.85 (0.72–0.99) | .037 | 7% |
| Hospital mortality | 0.66 (0.51–0.85) | .002 | 0% |
| Fracture (neg ctrl) | 1.04 (0.67–1.61) | .872 | 0% |
