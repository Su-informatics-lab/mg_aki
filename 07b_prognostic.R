#!/usr/bin/env Rscript
# 07b_prognostic.R — Prognostic Mg–AKI association by severity
# Fills eTable 1: OR per 1 mg/dL first postop Mg, adjusted for 27 PS covariates
# Run: Rscript 07b_prognostic.R

suppressPackageStartupMessages({ library(sandwich); library(lmtest) })
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
    first_hr="first_heartrate", has_vasopressor="vasopressor_6h")
  for (old in names(rmap)) {
    new <- rmap[[old]]
    if (old %in% names(d) && !new %in% names(d)) names(d)[names(d)==old] <- new
  }
  if (is.character(d$age)) { d$age <- suppressWarnings(as.numeric(d$age)); d$age[is.na(d$age)] <- 90 }
  if ("surgery_type" %in% names(d)) {
    d$surg_cabg <- as.integer(d$surgery_type=="cabg")
    d$surg_valve <- as.integer(d$surgery_type=="valve")
    d$surg_combined <- as.integer(d$surgery_type=="combined")
  }
  # Median-impute PS vars
  for (v in c("bmi","first_heartrate","first_calcium","first_potassium"))
    if (v %in% names(d) && any(is.na(d[[v]])))
      d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm=TRUE)
  # Lactate
  if ("first_lactate" %in% names(d)) {
    d$lactate_missing <- as.integer(is.na(d$first_lactate))
    d$first_lactate[is.na(d$first_lactate)] <- median(d$first_lactate, na.rm=TRUE)
  }
  d
}

# Covariates (same 27 as PS model, minus first_mg_value since it's the exposure here)
adj_vars <- c("age","is_female","bmi",
  "surg_cabg","surg_valve","surg_combined",
  "heart_failure","hypertension","diabetes","ckd","copd","pvd","stroke","liver_disease",
  "baseline_creatinine","egfr",
  "loop_diuretics","nsaids","acei_arb","ppi",
  "beta_blockers","steroids","antiarrhythmics",
  "first_potassium","first_calcium","first_heartrate","vasopressor_6h")

aki_defs <- c(
  "KDIGO stage >=1"      = "aki_kdigo1",
  "Cr ratio >=1.5x"      = "aki_primary",
  "Delta Cr >=0.3/48h"   = "aki_delta03",
  "KDIGO stage >=2"      = "aki_stage2",
  "KDIGO stage >=3"      = "aki_stage3"
)

run_prognostic <- function(d, db_name) {
  cat(sprintf("\n%s (N=%d)\n", db_name, nrow(d)))
  avail <- intersect(adj_vars, names(d))
  if ("first_lactate" %in% names(d)) avail <- c(avail, "first_lactate", "lactate_missing")

  rows <- list()
  for (lbl in names(aki_defs)) {
    outcome <- aki_defs[[lbl]]
    if (!outcome %in% names(d)) { cat(sprintf("  %-25s MISSING\n", lbl)); next }
    fml <- as.formula(paste(outcome, "~ first_mg_value +", paste(avail, collapse="+")))
    d_cc <- d[complete.cases(d[, c("first_mg_value", avail, outcome)]),]
    tryCatch({
      fit <- glm(fml, data=d_cc, family=binomial())
      ct <- coeftest(fit, vcov.=vcovHC(fit, type="HC1"))
      mg_row <- which(rownames(ct)=="first_mg_value")
      or <- exp(ct[mg_row,1]); lo <- exp(ct[mg_row,1]-1.96*ct[mg_row,2])
      hi <- exp(ct[mg_row,1]+1.96*ct[mg_row,2])
      p <- 2*pnorm(-abs(ct[mg_row,1]/ct[mg_row,2]))
      sig <- ifelse(p < 0.05, " *", "")
      cat(sprintf("  %-25s OR %.3f (%.3f-%.3f) P=%.4f%s\n", lbl, or, lo, hi, p, sig))
      rows[[length(rows)+1]] <- data.frame(db=db_name, outcome=lbl, outcome_col=outcome,
        surgery="All", or=round(or,3), lo=round(lo,3), hi=round(hi,3), p=round(p,4))
    }, error=function(e) cat(sprintf("  %-25s FAILED: %s\n", lbl, e$message)))
  }

  # Surgery-type interaction (KDIGO >=1 only)
  if ("aki_kdigo1" %in% names(d) && "surgery_type" %in% names(d)) {
    cat("  Surgery-type stratification (KDIGO >=1):\n")
    d$complex <- as.integer(d$surgery_type %in% c("valve","combined"))
    for (stype in c("Simple","Complex")) {
      sub <- if (stype=="Simple") d[d$complex==0,] else d[d$complex==1,]
      if (nrow(sub) < 50) next
      fml <- as.formula(paste("aki_kdigo1 ~ first_mg_value +", paste(avail, collapse="+")))
      sub_cc <- sub[complete.cases(sub[, c("first_mg_value", avail, "aki_kdigo1")]),]
      tryCatch({
        fit <- glm(fml, data=sub_cc, family=binomial())
        ct <- coeftest(fit, vcov.=vcovHC(fit, type="HC1"))
        mg_row <- which(rownames(ct)=="first_mg_value")
        or <- exp(ct[mg_row,1]); lo <- exp(ct[mg_row,1]-1.96*ct[mg_row,2])
        hi <- exp(ct[mg_row,1]+1.96*ct[mg_row,2])
        p <- 2*pnorm(-abs(ct[mg_row,1]/ct[mg_row,2]))
        cat(sprintf("    %-22s OR %.3f (%.3f-%.3f) P=%.4f\n", stype, or, lo, hi, p))
        rows[[length(rows)+1]] <- data.frame(db=db_name, outcome="KDIGO stage >=1",
          outcome_col="aki_kdigo1", surgery=stype,
          or=round(or,3), lo=round(lo,3), hi=round(hi,3), p=round(p,4))
      }, error=function(e) cat(sprintf("    %-22s FAILED: %s\n", stype, e$message)))
    }
  }

  do.call(rbind, rows)
}

cat(strrep("=", 60), "\n07b: Prognostic Mg-AKI Association\n", strrep("=", 60), "\n")
dat_e <- stdz(read.csv(file.path(RESULTS, "01_analysis_a_cohort.csv"), stringsAsFactors=FALSE))
dat_m <- stdz(read.csv(file.path(RESULTS, "04_mimic_cohort.csv"), stringsAsFactors=FALSE))

res_e <- run_prognostic(dat_e, "eICU")
res_m <- run_prognostic(dat_m, "MIMIC")
out <- rbind(res_e, res_m)
write.csv(out, file.path(RESULTS, "07b_prognostic.csv"), row.names=FALSE)
cat(sprintf("\nSaved: results/07b_prognostic.csv (%d rows)\n", nrow(out)))
