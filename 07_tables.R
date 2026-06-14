#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────
# 07_tables.R — Publication tables
#   Table 1: Baseline characteristics by treatment group (both DBs)
#   Table 2: Primary/secondary outcomes (eICU, MIMIC, pooled)
# ─────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyverse); library(tableone)
})

RESULTS <- path.expand("~/mg_aki/results")

# =====================================================================
# TABLE 1: Baseline Characteristics
# =====================================================================
cat("TABLE 1: Baseline Characteristics\n")

make_table1 <- function(path, db_label) {
  d <- read_csv(path, show_col_types = FALSE)
  d$trt <- as.integer(d$mg_supplementation)
  d$trt_label <- ifelse(d$trt == 1, "Mg Supplementation", "No Supplementation")

  vars <- c(
    "age_num", "is_female", "bmi", "surgery_type",
    "hx_chf", "hx_hypertension", "hx_diabetes", "hx_ckd",
    "hx_copd", "hx_pvd", "hx_stroke", "hx_liver",
    "baseline_cr", "baseline_egfr",
    "nephrotox_loop_diuretic", "nephrotox_nsaid",
    "nephrotox_acei_arb", "nephrotox_ppi",
    "has_betablocker", "has_steroid", "preop_antiarrhythmic",
    "has_vasopressor",
    "first_mg_value", "first_k_value", "first_ca_value", "first_hr",
    "aki_kdigo1", "hosp_mortality"
  )
  # Use age instead of age_num for MIMIC
  if (!"age_num" %in% names(d) && "age" %in% names(d)) {
    d$age_num <- d$age
  }

  vars <- intersect(vars, names(d))

  cat_vars <- intersect(c("surgery_type", "is_female",
    "hx_chf", "hx_hypertension", "hx_diabetes", "hx_ckd",
    "hx_copd", "hx_pvd", "hx_stroke", "hx_liver",
    "nephrotox_loop_diuretic", "nephrotox_nsaid",
    "nephrotox_acei_arb", "nephrotox_ppi",
    "has_betablocker", "has_steroid", "preop_antiarrhythmic",
    "has_vasopressor", "aki_kdigo1", "hosp_mortality"), vars)

  t1 <- CreateTableOne(vars = vars, strata = "trt_label",
                       factorVars = cat_vars, data = d, test = FALSE)
  out <- print(t1, smd = TRUE, printToggle = FALSE)

  # Convert to data frame and add database label
  df <- as.data.frame(out)
  df$Variable <- rownames(df)
  df$Database <- db_label
  df
}

t1_eicu <- make_table1(file.path(RESULTS, "01_analysis_a_cohort.csv"), "eICU")
t1_mimic <- make_table1(file.path(RESULTS, "04_mimic_cohort.csv"), "MIMIC-IV")

# Combine and save
t1_combined <- bind_rows(t1_eicu, t1_mimic) %>%
  select(Database, Variable, everything())

out_path <- file.path(RESULTS, "07_table1.csv")
write_csv(t1_combined, out_path)
cat(sprintf("  Saved: %s\n", out_path))

# Print for review
cat("\n  eICU:\n")
print(t1_eicu[, c("Variable", names(t1_eicu)[1:3])])
cat("\n  MIMIC-IV:\n")
print(t1_mimic[, c("Variable", names(t1_mimic)[1:3])])


# =====================================================================
# TABLE 2: Primary and Secondary Outcomes
# =====================================================================
cat("\n\nTABLE 2: Primary and Secondary Outcomes\n")

eicu  <- read_csv(file.path(RESULTS, "03_results_summary.csv"), show_col_types = FALSE)
mimic <- read_csv(file.path(RESULTS, "05_mimic_results_summary.csv"), show_col_types = FALSE)
meta  <- read_csv(file.path(RESULTS, "06_meta_results.csv"), show_col_types = FALSE)

fmt_or <- function(est, lo, hi, p) {
  sprintf("%.2f (%.2f-%.2f), P=%s",
          est, lo, hi,
          ifelse(p < 0.001, "<.001", sprintf("%.3f", p)))
}

# Define outcome rows for Table 2
outcome_map <- list(
  list(label = "AKI, KDIGO stage ≥1 (primary)",
       eicu = "ALL_aki_KDIGO >=1", mimic = "ALL_aki_KDIGO >=1",
       meta_row = "KDIGO Stage >=1"),
  list(label = "AKI, creatinine ratio ≥1.5×",
       eicu = "ALL_aki_Ratio >=1.5x", mimic = "ALL_aki_Ratio >=1.5x",
       meta_row = "Ratio >=1.5x"),
  list(label = "AKI, stage ≥2",
       eicu = "ALL_aki_Stage >=2", mimic = "ALL_aki_Stage >=2",
       meta_row = "Stage >=2"),
  list(label = "AKI, stage ≥3",
       eicu = "ALL_aki_Stage >=3", mimic = "ALL_aki_Stage >=3",
       meta_row = "Stage >=3"),
  list(label = "AKI within 48 h",
       eicu = "ALL_tw_AKI 1.5x <=48h", mimic = "ALL_tw_AKI 1.5x <=48h",
       meta_row = "AKI 1.5x <=48h"),
  list(label = "Hospital mortality",
       eicu = "ALL_sec_Hospital mortality", mimic = "ALL_sec_Hospital mortality",
       meta_row = "Mortality (all+all)"),
  list(label = "Hospital mortality (hypo+all)",
       eicu = "HYPO_sec_Hospital mortality", mimic = "ALL_sec_Hospital mortality",
       meta_row = "Mortality (hypo+all)"),
  list(label = "Fracture (negative control)",
       eicu = "ALL_nc_Fracture", mimic = "ALL_nc_Fracture",
       meta_row = "Fracture (neg ctrl)"),
  list(label = "Encephalopathy (exploratory)",
       eicu = "ALL_neuro_Encephalopathy", mimic = "ALL_neuro_Encephalopathy",
       meta_row = "Encephalopathy")
)

t2_rows <- list()
for (om in outcome_map) {
  e <- eicu %>% filter(model == om$eicu)
  m <- mimic %>% filter(model == om$mimic)
  p <- meta %>% filter(outcome == om$meta_row)

  t2_rows[[length(t2_rows) + 1]] <- tibble(
    Outcome = om$label,
    `eICU OR (95% CI)` = if (nrow(e) > 0) fmt_or(e$estimate[1], e$conf.low[1], e$conf.high[1], e$p.value[1]) else "—",
    `eICU Events` = if (nrow(e) > 0) e$n_events[1] else NA,
    `MIMIC OR (95% CI)` = if (nrow(m) > 0) fmt_or(m$estimate[1], m$conf.low[1], m$conf.high[1], m$p.value[1]) else "—",
    `MIMIC Events` = if (nrow(m) > 0) m$n_events[1] else NA,
    `Pooled OR (95% CI)` = if (nrow(p) > 0) fmt_or(p$pooled_or[1], p$pooled_lo[1], p$pooled_hi[1], p$pooled_p[1]) else "—",
    `I²` = if (nrow(p) > 0) sprintf("%.0f%%", p$I2[1]) else "—"
  )
}

t2 <- bind_rows(t2_rows)
out_path2 <- file.path(RESULTS, "07_table2.csv")
write_csv(t2, out_path2)
cat(sprintf("  Saved: %s\n", out_path2))

# Print
print(t2, width = Inf)

cat("\n07_tables.R COMPLETE\n")
