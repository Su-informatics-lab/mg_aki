#!/usr/bin/env Rscript
# ============================================================================
# 02c_strict_sensitivity.R — Lenient vs Strict baseline comparison
#
# Lenient (primary):  all patients (current cohort)
# Strict (sensitivity): exclude pre_lab_mg_supp == 1
#
# Runs key analyses on both, produces side-by-side comparison.
# Output: results/02c_strict_comparison.csv
# ============================================================================

suppressPackageStartupMessages({
  library(sandwich); library(lmtest)
})
RESULTS <- path.expand("~/mg_aki/results")

stdz <- function(d) {
  rmap <- c(mg_supplementation="mg_supp", hosp_mortality="hospital_mortality",
    age_num="age", baseline_cr="baseline_creatinine", baseline_egfr="egfr",
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
    if (old %in% names(d) && !new %in% names(d))
      names(d)[names(d) == old] <- new
  }
  if (is.character(d$age)) { d$age <- suppressWarnings(as.numeric(d$age)); d$age[is.na(d$age)] <- 90 }
  if ("surgery_type" %in% names(d)) {
    d$surg_cabg <- as.integer(d$surgery_type=="cabg")
    d$surg_valve <- as.integer(d$surgery_type=="valve")
    d$surg_combined <- as.integer(d$surgery_type=="combined")
  }
  if ("first_lactate" %in% names(d)) {
    d$lactate_missing <- as.integer(is.na(d$first_lactate))
    d$first_lactate[is.na(d$first_lactate)] <- median(d$first_lactate, na.rm=TRUE)
  }
  for (v in c("bmi","first_heartrate","first_calcium","first_potassium"))
    if (v %in% names(d) && any(is.na(d[[v]])))
      d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm=TRUE)
  d
}

wglm <- function(formula, data, w, cluster=NULL) {
  data$.w <- w
  fit <- glm(formula, data=data, weights=.w, family=quasibinomial())
  vc <- if (!is.null(cluster) && length(unique(cluster))>1)
    vcovCL(fit, cluster=cluster) else vcovHC(fit, type="HC1")
  ct <- coeftest(fit, vcov.=vc)
  tr <- 2
  list(or=exp(ct[tr,1]), lo=exp(ct[tr,1]-1.96*ct[tr,2]),
       hi=exp(ct[tr,1]+1.96*ct[tr,2]), p=2*pnorm(-abs(ct[tr,1]/ct[tr,2])))
}

run_key <- function(dat, db, version, has_cluster) {
  cat(sprintf("\n  %s %s [N=%d, trt=%d]\n", db, version, nrow(dat), sum(dat$mg_supp)))

  ps_covars <- intersect(c("age","is_female","bmi",
    "surg_cabg","surg_valve","surg_combined",
    "heart_failure","hypertension","diabetes","ckd",
    "copd","pvd","stroke","liver_disease",
    "baseline_creatinine","egfr",
    "loop_diuretics","nsaids","acei_arb","ppi",
    "beta_blockers","steroids","antiarrhythmics",
    "first_potassium","first_calcium","first_heartrate",
    "vasopressor_6h","first_mg_value",
    "first_lactate","lactate_missing"), names(dat))

  d <- dat[complete.cases(dat[,ps_covars]),]
  ps_fit <- glm(as.formula(paste("mg_supp ~", paste(ps_covars, collapse="+"))),
                data=d, family=binomial())
  d$ps <- pmax(pmin(fitted(ps_fit), 0.99), 0.01)
  d$ow <- ifelse(d$mg_supp==1, 1-d$ps, d$ps)
  prev <- mean(d$mg_supp)
  d$iptw <- pmax(pmin(ifelse(d$mg_supp==1, prev/d$ps, (1-prev)/(1-d$ps)),
                       quantile(ifelse(d$mg_supp==1, prev/d$ps, (1-prev)/(1-d$ps)), 0.99)),
                 quantile(ifelse(d$mg_supp==1, prev/d$ps, (1-prev)/(1-d$ps)), 0.01))

  cl <- if (has_cluster && "hospitalid" %in% names(d)) d$hospitalid else NULL

  rows <- list()
  # OW
  r <- wglm(aki_kdigo1 ~ mg_supp, d, d$ow, cl)
  cat(sprintf("    OW:   OR %.3f (%.3f-%.3f) P=%.4f\n", r$or, r$lo, r$hi, r$p))
  rows[[1]] <- data.frame(db=db, version=version, analysis="ow", n=nrow(d),
                           or=round(r$or,3), lo=round(r$lo,3), hi=round(r$hi,3), p=round(r$p,4))
  # IPTW
  r <- wglm(aki_kdigo1 ~ mg_supp, d, d$iptw, cl)
  cat(sprintf("    IPTW: OR %.3f (%.3f-%.3f) P=%.4f\n", r$or, r$lo, r$hi, r$p))
  rows[[2]] <- data.frame(db=db, version=version, analysis="iptw", n=nrow(d),
                           or=round(r$or,3), lo=round(r$lo,3), hi=round(r$hi,3), p=round(r$p,4))

  # Mg >2.3 stratum
  d_hi <- d[!is.na(d$first_mg_value) & d$first_mg_value > 2.3,]
  if (sum(d_hi$mg_supp) >= 10) {
    ps_no_mg <- ps_covars[ps_covars != "first_mg_value"]
    ps2 <- glm(as.formula(paste("mg_supp ~", paste(ps_no_mg, collapse="+"))),
               data=d_hi, family=binomial())
    d_hi$ow2 <- ifelse(d_hi$mg_supp==1, 1-pmax(pmin(fitted(ps2),0.99),0.01),
                        pmax(pmin(fitted(ps2),0.99),0.01))
    r <- wglm(aki_kdigo1 ~ mg_supp, d_hi, d_hi$ow2, if(has_cluster) d_hi$hospitalid else NULL)
    cat(sprintf("    >2.3: OR %.3f (%.3f-%.3f) P=%.4f  [N=%d trt=%d]\n",
                r$or, r$lo, r$hi, r$p, nrow(d_hi), sum(d_hi$mg_supp)))
    rows[[3]] <- data.frame(db=db, version=version, analysis="mg_gt23_ow", n=nrow(d_hi),
                             or=round(r$or,3), lo=round(r$lo,3), hi=round(r$hi,3), p=round(r$p,4))
  }

  # Interaction
  tryCatch({
    int_fit <- glm(as.formula(paste("aki_kdigo1 ~ mg_supp * first_mg_value +",
                                     paste(ps_covars[ps_covars!="first_mg_value"], collapse="+"))),
                   data=d, family=binomial())
    ct <- coeftest(int_fit, vcov.=vcovHC(int_fit, type="HC1"))
    ir <- grep("mg_supp:first_mg_value", rownames(ct))
    if (length(ir)>0) {
      cat(sprintf("    Interaction: OR %.3f P=%.4f\n", exp(ct[ir,1]), ct[ir,4]))
      rows[[length(rows)+1]] <- data.frame(db=db, version=version, analysis="interaction",
        n=nrow(d), or=round(exp(ct[ir,1]),3), lo=NA, hi=NA, p=round(ct[ir,4],4))
    }
  }, error=function(e){})

  # Composite
  if ("aki_or_death" %in% names(d)) {
    r <- wglm(aki_or_death ~ mg_supp, d, d$ow, cl)
    cat(sprintf("    Composite: OR %.3f (%.3f-%.3f) P=%.4f\n", r$or, r$lo, r$hi, r$p))
    rows[[length(rows)+1]] <- data.frame(db=db, version=version, analysis="composite_ow",
      n=nrow(d), or=round(r$or,3), lo=round(r$lo,3), hi=round(r$hi,3), p=round(r$p,4))
  }

  do.call(rbind, rows)
}

# ── Load + run ────────────────────────────────────────────────────────
cat("=" , rep("=",59), "\n")
cat("LENIENT vs STRICT COMPARISON\n")
cat(rep("=",60), "\n")

all_rows <- list()
for (db_info in list(
  list(file="01_analysis_a_cohort.csv", db="eICU", cl=TRUE),
  list(file="04_mimic_cohort.csv", db="MIMIC", cl=FALSE)
)) {
  dat <- stdz(read.csv(file.path(RESULTS, db_info$file), stringsAsFactors=FALSE))
  if (!"pre_lab_mg_supp" %in% names(dat)) {
    cat(sprintf("  %s: no pre_lab_mg_supp column — run 01_etl.py first\n", db_info$db))
    next
  }

  n_pre <- sum(dat$pre_lab_mg_supp)
  cat(sprintf("\n%s: N=%d, pre-lab Mg supp=%d (%.1f%%)\n",
              db_info$db, nrow(dat), n_pre, 100*n_pre/nrow(dat)))

  # Lenient (primary)
  all_rows[[length(all_rows)+1]] <- run_key(dat, db_info$db, "lenient", db_info$cl)

  # Strict (sensitivity)
  dat_strict <- dat[dat$pre_lab_mg_supp == 0,]
  all_rows[[length(all_rows)+1]] <- run_key(dat_strict, db_info$db, "strict", db_info$cl)
}

res <- do.call(rbind, all_rows)
write.csv(res, file.path(RESULTS, "02c_strict_comparison.csv"), row.names=FALSE)

cat(sprintf("\n%s\nSIDE-BY-SIDE COMPARISON\n%s\n", strrep("=",60), strrep("=",60)))
for (a in unique(res$analysis)) {
  cat(sprintf("  %s:\n", a))
  for (db in unique(res$db)) {
    len <- res[res$db==db & res$analysis==a & res$version=="lenient",]
    str <- res[res$db==db & res$analysis==a & res$version=="strict",]
    if (nrow(len)==0) next
    cat(sprintf("    %s lenient: OR %.3f P=%.4f (N=%d)\n", db, len$or, len$p, len$n))
    if (nrow(str)>0) {
      delta <- abs(len$or - str$or)
      cat(sprintf("    %s strict:  OR %.3f P=%.4f (N=%d)  |ΔOR|=%.3f\n",
                  db, str$or, str$p, str$n, delta))
    }
  }
}

cat(sprintf("\n✓ Saved: results/02c_strict_comparison.csv\n"))
cat("  If |ΔOR| < 0.05 across all analyses → pre-lab contamination is negligible\n")
