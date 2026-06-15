#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────
# 05_tables.R — Publication tables
#   Table 1: Baseline characteristics by treatment group (both DBs)
#   Table 2: Primary (AC) + sensitivity + exploratory + control outcomes
#
# Reads: results/01_analysis_a_cohort.csv, results/04_mimic_cohort.csv,
#        results/02_results.csv
# ─────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({ library(tableone) })
RESULTS <- path.expand("~/mg_aki/results")

# =====================================================================
# TABLE 1
# =====================================================================
cat("TABLE 1: Baseline Characteristics\n")

make_table1 <- function(path, db_label) {
  d <- read.csv(file.path(RESULTS, path), stringsAsFactors = FALSE)
  # Standardize treatment column
  if ("mg_supplementation" %in% names(d)) d$trt <- d$mg_supplementation
  d$trt_label <- ifelse(d$trt == 1, "Mg Supplementation", "No Supplementation")
  # Standardize age
  if (!"age_num" %in% names(d) && "age" %in% names(d)) d$age_num <- d$age

  vars <- intersect(c(
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
  ), names(d))

  cat_vars <- intersect(c("surgery_type", "is_female",
    "hx_chf", "hx_hypertension", "hx_diabetes", "hx_ckd",
    "hx_copd", "hx_pvd", "hx_stroke", "hx_liver",
    "nephrotox_loop_diuretic", "nephrotox_nsaid",
    "nephrotox_acei_arb", "nephrotox_ppi",
    "has_betablocker", "has_steroid", "preop_antiarrhythmic",
    "has_vasopressor", "aki_kdigo1", "hosp_mortality"), vars)

  t1 <- CreateTableOne(vars=vars, strata="trt_label", factorVars=cat_vars,
                       data=d, test=FALSE)
  out <- print(t1, smd=TRUE, printToggle=FALSE)
  df <- as.data.frame(out); df$Variable <- rownames(df); df$Database <- db_label
  df
}

t1_eicu <- make_table1("01_analysis_a_cohort.csv", "eICU")
t1_mimic <- make_table1("04_mimic_cohort.csv", "MIMIC-IV")
t1_combined <- rbind(t1_eicu, t1_mimic)
write.csv(t1_combined, file.path(RESULTS, "05_table1.csv"), row.names=FALSE)
cat("  Saved: 05_table1.csv\n")
cat("\n  eICU:\n"); print(t1_eicu[, c("Variable", names(t1_eicu)[1:3])])
cat("\n  MIMIC:\n"); print(t1_mimic[, c("Variable", names(t1_mimic)[1:3])])

# =====================================================================
# TABLE 2
# =====================================================================
cat("\n\nTABLE 2: Primary and Sensitivity Outcomes\n")

res <- read.csv(file.path(RESULTS, "02_results.csv"), stringsAsFactors=FALSE)
# Use m=5 (first M_VALUE saved as primary)
res <- res[is.na(res$m) | res$m == min(res$m, na.rm=TRUE), ]

fmt_or <- function(or, lo, hi, p) {
  sprintf("%.2f (%.2f-%.2f), P=%s", or, lo, hi,
          ifelse(p < 0.001, "<.001", sprintf("%.3f", p)))
}

get_est <- function(db_name, analysis_name) {
  r <- res[res$db == db_name & res$analysis == analysis_name, ]
  if (nrow(r) == 0) return("—")
  fmt_or(r$or[1], r$lo[1], r$hi[1], r$p[1])
}

get_i2 <- function(analysis_name) {
  r <- res[res$db == "Pooled" & res$analysis == analysis_name, ]
  if (nrow(r) == 0 || is.na(r$I2[1])) return("—")
  sprintf("%.0f%%", r$I2[1])
}

# Define Table 2 rows
t2_spec <- list(
  list(label="Primary: AC (Mg+K vs K-only, OW)", a="ac_aki1"),
  list(label="Sensitivity: IPTW", a="iptw_aki1"),
  list(label="Sensitivity: Overlap weighting", a="ow_aki1"),
  list(label="Sensitivity: PS matching", a="psm_aki1"),
  list(label="Hospital mortality (OW)", a="ow_mort"),
  list(label="Encephalopathy (OW)", a="ow_enceph"),
  list(label="Fracture neg. control (OW)", a="ow_frac")
)

t2_rows <- lapply(t2_spec, function(s) {
  data.frame(
    Outcome = s$label,
    eICU = get_est("eICU", s$a),
    MIMIC = get_est("MIMIC", s$a),
    Pooled = get_est("Pooled", s$a),
    I2 = get_i2(s$a),
    stringsAsFactors = FALSE
  )
})

t2 <- do.call(rbind, t2_rows)
write.csv(t2, file.path(RESULTS, "05_table2.csv"), row.names=FALSE)
cat("  Saved: 05_table2.csv\n\n")
print(t2, right=FALSE)

cat("\n05_tables.R COMPLETE\n")
