#!/usr/bin/env Rscript
# probe_clean_comparator.R — Treated vs Never-Received-Mg
# Excludes "untreated" patients who actually received Mg after 6h
# Run after: python probe_clean_comparator.py (adds comparator_status column)

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
    preop_antiarrhythmic="antiarrhythmics", first_k_value="first_potassium",
    first_ca_value="first_calcium", first_hr="first_heartrate",
    has_vasopressor="vasopressor_6h")
  for (old in names(rmap)) {
    new <- rmap[[old]]
    if (old %in% names(d) && !new %in% names(d))
      names(d)[names(d) == old] <- new
  }
  if (is.character(d$age)) { d$age <- suppressWarnings(as.numeric(d$age)); d$age[is.na(d$age)] <- 90 }
  if ("surgery_type" %in% names(d)) {
    d$surg_cabg <- as.integer(d$surgery_type == "cabg")
    d$surg_valve <- as.integer(d$surgery_type == "valve")
    d$surg_combined <- as.integer(d$surgery_type == "combined")
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

ps_covars <- c("age","is_female","bmi","surg_cabg","surg_valve","surg_combined",
  "heart_failure","hypertension","diabetes","ckd","copd","pvd","stroke",
  "liver_disease","baseline_creatinine","egfr","loop_diuretics","nsaids",
  "acei_arb","ppi","beta_blockers","steroids","antiarrhythmics",
  "first_potassium","first_calcium","first_heartrate","vasopressor_6h",
  "first_mg_value")

wglm <- function(fml, data, w, cl=NULL) {
  data$.w <- w
  fit <- glm(fml, data=data, weights=.w, family=quasibinomial())
  vc <- if (!is.null(cl) && length(unique(cl))>1) vcovCL(fit, cluster=cl) else vcovHC(fit, type="HC1")
  ct <- coeftest(fit, vcov.=vc)
  tr <- 2
  list(or=exp(ct[tr,1]), lo=exp(ct[tr,1]-1.96*ct[tr,2]),
       hi=exp(ct[tr,1]+1.96*ct[tr,2]), p=2*pnorm(-abs(ct[tr,1]/ct[tr,2])))
}

run_analysis <- function(d, label, covs, has_cluster=FALSE) {
  cv <- intersect(covs, names(d))
  if ("first_lactate" %in% names(d)) cv <- c(cv, "first_lactate", "lactate_missing")
  d <- d[complete.cases(d[, cv]), ]
  ps_fit <- glm(as.formula(paste("mg_supp ~", paste(cv, collapse="+"))),
                data=d, family=binomial())
  d$ps <- pmax(pmin(fitted(ps_fit), 0.99), 0.01)
  cl <- if (has_cluster && "hospitalid" %in% names(d)) d$hospitalid else NULL

  # OW
  d$ow <- ifelse(d$mg_supp==1, 1-d$ps, d$ps)
  r_ow <- wglm(aki_kdigo1 ~ mg_supp, d, d$ow, cl)

  # IPTW
  prev <- mean(d$mg_supp)
  d$iptw <- ifelse(d$mg_supp==1, prev/d$ps, (1-prev)/(1-d$ps))
  d$iptw <- pmax(pmin(d$iptw, quantile(d$iptw,.99)), quantile(d$iptw,.01))
  r_iptw <- wglm(aki_kdigo1 ~ mg_supp, d, d$iptw, cl)

  cat(sprintf("  %-25s OW:   OR=%.3f (%.3f-%.3f) P=%.4f  n=%d trt=%d\n",
      label, r_ow$or, r_ow$lo, r_ow$hi, r_ow$p, nrow(d), sum(d$mg_supp)))
  cat(sprintf("  %-25s IPTW: OR=%.3f (%.3f-%.3f) P=%.4f\n",
      label, r_iptw$or, r_iptw$lo, r_iptw$hi, r_iptw$p))
  list(ow=r_ow, iptw=r_iptw, n=nrow(d), n_trt=sum(d$mg_supp))
}

# ── Load data ────────────────────────────────────────────────────
eicu <- stdz(read.csv(file.path(RESULTS, "01_analysis_a_cohort.csv"), stringsAsFactors=FALSE))
mimic <- stdz(read.csv(file.path(RESULTS, "04_mimic_cohort.csv"), stringsAsFactors=FALSE))

for (db_name in c("MIMIC-IV", "eICU")) {
  d <- if (db_name == "MIMIC-IV") mimic else eicu
  hc <- db_name == "eICU"

  cat(sprintf("\n%s\n%s: CLEAN COMPARATOR ANALYSIS\n%s\n", strrep("=",60), db_name, strrep("=",60)))

  if (!"comparator_status" %in% names(d)) {
    cat("  ERROR: comparator_status column not found. Run probe_clean_comparator.py first.\n")
    next
  }

  n_trt <- sum(d$comparator_status == "treated")
  n_contam <- sum(d$comparator_status == "contaminated")
  n_clean <- sum(d$comparator_status == "clean_untreated")
  cat(sprintf("  Treated: %d | Contaminated controls: %d | Clean controls: %d\n",
      n_trt, n_contam, n_clean))
  cat(sprintf("  Contamination rate: %.0f%% of untreated received Mg after 6h\n\n",
      100*n_contam/(n_contam+n_clean)))

  # AKI rates
  for (s in c("treated", "contaminated", "clean_untreated")) {
    sub <- d[d$comparator_status == s, ]
    cat(sprintf("  %-20s AKI=%.1f%% (n=%d)\n", s, 100*mean(sub$aki_kdigo1), nrow(sub)))
  }
  cat("\n")

  # ── Standard: treated vs ALL untreated ──
  cat("  ── Standard (treated vs all untreated) ──\n")
  d_std <- d[d$comparator_status %in% c("treated", "contaminated", "clean_untreated"), ]
  run_analysis(d_std, "Standard", ps_covars, hc)

  # ── Clean: treated vs NEVER-treated only ──
  cat("\n  ── Clean (treated vs never-received-Mg) ──\n")
  d_clean <- d[d$comparator_status %in% c("treated", "clean_untreated"), ]
  run_analysis(d_clean, "Clean comparator", ps_covars, hc)

  # ── AC within clean ──
  if ("ac_group" %in% names(d_clean)) {
    cat("\n  ── AC within clean comparator ──\n")
    d_ac <- d_clean[d_clean$ac_group %in% c("mg_k", "k_only"), ]
    # Need to check if k_only patients are also contaminated
    # Some K-only patients may have received Mg after 6h
    d_ac$mg_supp <- as.integer(d_ac$ac_group == "mg_k")
    if (sum(d_ac$mg_supp) >= 10 && sum(d_ac$mg_supp==0) >= 10) {
      run_analysis(d_ac, "AC clean", ps_covars, FALSE)
    } else {
      cat("  AC clean: insufficient sample after exclusion\n")
    }
  }

  # ── Mg > 2.3 stratum within clean ──
  cat("\n  ── >2.3 stratum within clean comparator ──\n")
  d_hi <- d_clean[!is.na(d_clean$first_mg_value) & d_clean$first_mg_value > 2.3, ]
  if (sum(d_hi$mg_supp) >= 15) {
    covs_nomg <- ps_covars[ps_covars != "first_mg_value"]
    run_analysis(d_hi, ">2.3 clean", covs_nomg, hc)
  } else {
    cat("  >2.3 clean: insufficient treated\n")
  }
}

cat(sprintf("\n%s\nVERDICT\n%s\n", strrep("=",60), strrep("=",60)))
cat("If MIMIC clean-comparator OR < standard OR:
  → contamination was diluting the signal
  → 'MIMIC null' is partially artifactual
  → major finding for the paper

If MIMIC clean-comparator OR ≈ standard OR:
  → contamination is not the explanation
  → MIMIC null remains real
\n")
