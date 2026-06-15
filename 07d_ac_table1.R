#!/usr/bin/env Rscript
# 07d_ac_table1.R — eTable 3: AC subgroup baseline characteristics
# Saves CreateTableOne for the K+-repleted population to CSV.
# Run: Rscript 07d_ac_table1.R

suppressPackageStartupMessages({ library(tableone) })
RESULTS <- path.expand("~/mg_aki/results")

stdz <- function(d) {
  rmap <- c(mg_supplementation="mg_supp", age_num="age",
    baseline_cr="baseline_creatinine", baseline_egfr="egfr",
    hx_chf="heart_failure", hx_hypertension="hypertension",
    hx_diabetes="diabetes", hx_ckd="ckd", hx_copd="copd",
    hx_pvd="pvd", hx_stroke="stroke", hx_liver="liver_disease",
    nephrotox_loop_diuretic="loop_diuretics", nephrotox_nsaid="nsaids",
    nephrotox_acei_arb="acei_arb", nephrotox_ppi="ppi",
    has_betablocker="beta_blockers", has_steroid="steroids",
    preop_antiarrhythmic="antiarrhythmics",
    first_k_value="first_potassium", first_ca_value="first_calcium",
    first_hr="first_heartrate", has_vasopressor="vasopressor_6h",
    hosp_mortality="hospital_mortality")
  for (old in names(rmap)) {
    new <- rmap[[old]]
    if (old %in% names(d) && !new %in% names(d)) names(d)[names(d)==old] <- new
  }
  if (is.character(d$age)) { d$age <- suppressWarnings(as.numeric(d$age)); d$age[is.na(d$age)] <- 90 }
  d
}

make_ac_table <- function(path, db_label) {
  d <- stdz(read.csv(file.path(RESULTS, path), stringsAsFactors=FALSE))
  if (!"ac_group" %in% names(d)) { cat(sprintf("  %s: no ac_group\n", db_label)); return(NULL) }
  d_ac <- d[d$ac_group %in% c("mg_k","k_only"),]
  d_ac$trt_label <- ifelse(d_ac$ac_group=="mg_k", "Mg+K", "K-only")

  vars <- intersect(c("age","is_female","bmi","surgery_type",
    "heart_failure","hypertension","diabetes","ckd","copd","pvd","stroke","liver_disease",
    "baseline_creatinine","egfr",
    "loop_diuretics","nsaids","acei_arb","ppi",
    "beta_blockers","steroids","antiarrhythmics","vasopressor_6h",
    "first_mg_value","first_potassium","first_calcium","first_heartrate",
    "aki_kdigo1","hospital_mortality"), names(d_ac))

  cat_vars <- intersect(c("surgery_type","is_female",
    "heart_failure","hypertension","diabetes","ckd","copd","pvd","stroke","liver_disease",
    "loop_diuretics","nsaids","acei_arb","ppi",
    "beta_blockers","steroids","antiarrhythmics","vasopressor_6h",
    "aki_kdigo1","hospital_mortality"), vars)

  cat(sprintf("\n  %s AC population: N=%d (Mg+K=%d, K-only=%d)\n",
      db_label, nrow(d_ac), sum(d_ac$trt_label=="Mg+K"), sum(d_ac$trt_label=="K-only")))

  t1 <- CreateTableOne(vars=vars, strata="trt_label", factorVars=cat_vars,
                       data=d_ac, test=FALSE)
  out <- print(t1, smd=TRUE, printToggle=FALSE)
  cat("\n"); print(print(t1, smd=TRUE, printToggle=TRUE))

  df <- as.data.frame(out)
  df$Variable <- rownames(df)
  df$Database <- db_label
  df
}

cat(strrep("=",60), "\n07d: AC Subgroup Baseline Characteristics (eTable 3)\n",
    strrep("=",60), "\n")

t1_e <- make_ac_table("01_analysis_a_cohort.csv", "eICU")
t1_m <- make_ac_table("04_mimic_cohort.csv", "MIMIC-IV")

combined <- rbind(t1_e, t1_m)
outpath <- file.path(RESULTS, "07d_ac_table1.csv")
write.csv(combined, outpath, row.names=FALSE)
cat(sprintf("\nSaved: %s\n", outpath))
