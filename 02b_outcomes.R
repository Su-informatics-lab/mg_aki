#!/usr/bin/env Rscript
# ============================================================================
# 02b_outcomes.R — Additional outcomes: POAF, composite AKI-or-death
#
# Runs the same MICE + OW/IPTW framework as 02_analysis.R but for:
#   1. POAF (new-onset postoperative atrial fibrillation)
#   2. AKI-or-death composite (competing risk sensitivity)
#   3. AC versions of both
#
# Run AFTER 02_analysis.R:  Rscript 02b_outcomes.R
# Output: results/02b_additional_outcomes.csv
# ============================================================================

cat("======================================================================\n")
cat("02b: ADDITIONAL OUTCOMES (POAF + Composite AKI-or-Death)\n")
cat("======================================================================\n")

needed <- c("mice", "sandwich", "lmtest")
miss <- needed[!sapply(needed, requireNamespace, quietly = TRUE)]
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org", quiet = TRUE)
suppressPackageStartupMessages({
  library(mice); library(sandwich); library(lmtest)
})

SEED <- 42
M_IMP <- 5
RESULTS <- path.expand("~/mg_aki/results")

# ── Helpers ───────────────────────────────────────────────────────────
pool_simple <- function(ests, ses) {
  m <- length(ests); qbar <- mean(ests); ubar <- mean(ses^2)
  b <- var(ests); tv <- ubar + (1 + 1/m) * b; se <- sqrt(tv)
  list(or = exp(qbar), lo = exp(qbar - 1.96*se),
       hi = exp(qbar + 1.96*se), p = 2*pnorm(-abs(qbar/se)))
}

wglm <- function(formula, data, w, cluster = NULL) {
  data$.w <- w
  fit <- glm(formula, data = data, weights = .w, family = quasibinomial())
  vc <- if (!is.null(cluster) && length(unique(cluster)) > 1)
    vcovCL(fit, cluster = cluster) else vcovHC(fit, type = "HC1")
  ct <- coeftest(fit, vcov. = vc)
  tr <- which(rownames(ct) == "mg_supp")
  if (length(tr) == 0) tr <- 2
  list(logOR = ct[tr, 1], se = ct[tr, 2])
}

# ── Standardize (same as 02_analysis.R) ───────────────────────────────
stdz <- function(d) {
  rmap <- c(mg_supplementation="mg_supp", hosp_mortality="hospital_mortality",
    age_num="age", hx_chf="heart_failure", hx_hypertension="hypertension",
    hx_diabetes="diabetes", hx_ckd="ckd", hx_copd="copd",
    hx_pvd="pvd", hx_stroke="stroke", hx_liver="liver_disease",
    baseline_cr="baseline_creatinine", baseline_egfr="egfr",
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
  if (is.character(d$age)) {
    d$age <- suppressWarnings(as.numeric(d$age))
    d$age[is.na(d$age)] <- 90
  }
  if ("surgery_type" %in% names(d)) {
    d$surg_cabg     <- as.integer(d$surgery_type == "cabg")
    d$surg_valve    <- as.integer(d$surgery_type == "valve")
    d$surg_combined <- as.integer(d$surgery_type == "combined")
  }
  if ("first_lactate" %in% names(d)) {
    d$lactate_missing <- as.integer(is.na(d$first_lactate))
    d$first_lactate[is.na(d$first_lactate)] <- median(d$first_lactate, na.rm=TRUE)
  }
  d
}

# ── Load ──────────────────────────────────────────────────────────────
cat("Loading cohorts...\n")
dat_e <- stdz(read.csv(file.path(RESULTS, "01_analysis_a_cohort.csv"), stringsAsFactors=FALSE))
dat_m <- stdz(read.csv(file.path(RESULTS, "04_mimic_cohort.csv"), stringsAsFactors=FALSE))
cat(sprintf("  eICU: N=%d  MIMIC: N=%d\n", nrow(dat_e), nrow(dat_m)))

# ── Analysis function ─────────────────────────────────────────────────
run_additional <- function(dat, db_name, has_cluster = FALSE) {
  cat(sprintf("\n%s  [N=%d]\n", db_name, nrow(dat)))

  # PS covariates
  ps_covars <- intersect(c("age","is_female","bmi",
    "surg_cabg","surg_valve","surg_combined",
    "heart_failure","hypertension","diabetes","ckd",
    "copd","pvd","stroke","liver_disease",
    "baseline_creatinine","egfr",
    "loop_diuretics","nsaids","acei_arb","ppi",
    "beta_blockers","steroids","antiarrhythmics",
    "first_potassium","first_calcium","first_heartrate",
    "vasopressor_6h","transfusion_6h","first_mg_value",
    "first_lactate","lactate_missing"), names(dat))

  ps_formula <- as.formula(paste("mg_supp ~", paste(ps_covars, collapse="+")))

  # MICE
  mice_vars <- intersect(c("bmi","first_heartrate","first_calcium","first_potassium"),
                          names(dat))
  any_missing <- any(sapply(mice_vars, function(v) sum(is.na(dat[[v]])) > 0))
  if (any_missing) {
    imp_preds <- unique(c(mice_vars, "age","is_female","baseline_creatinine",
                          "mg_supp","aki_kdigo1"))
    imp_preds <- imp_preds[imp_preds %in% names(dat)]
    imp <- mice(dat[, imp_preds], m=M_IMP, method="pmm",
                seed=SEED, printFlag=FALSE, maxit=10)
    n_imp <- M_IMP
  } else {
    n_imp <- 1
  }

  # Outcomes to run
  outcomes <- list(
    list(name="ow_poaf",      outcome="poaf",        trt="mg_supp",  filter_na=TRUE),
    list(name="ow_composite",  outcome="aki_or_death", trt="mg_supp",  filter_na=FALSE),
    list(name="iptw_composite",outcome="aki_or_death", trt="mg_supp",  filter_na=FALSE),
    list(name="ac_poaf",       outcome="poaf",        trt="ac_trt",   filter_na=TRUE),
    list(name="ac_composite",  outcome="aki_or_death", trt="ac_trt",   filter_na=FALSE),
    # Transfusion negative control (complexity-specific)
    list(name="ow_transfusion",  outcome="nc_transfusion", trt="mg_supp", filter_na=FALSE),
    list(name="ac_transfusion",  outcome="nc_transfusion", trt="ac_trt",  filter_na=FALSE)
  )

  stores <- lapply(outcomes, function(x) list(est=numeric(0), se=numeric(0)))
  names(stores) <- sapply(outcomes, `[[`, "name")

  for (i in seq_len(n_imp)) {
    d <- dat
    if (any_missing && n_imp > 1) {
      imputed <- complete(imp, i)
      for (v in mice_vars)
        if (v %in% names(imputed)) d[[v]] <- imputed[[v]]
    }

    # Median-impute remaining
    for (v in ps_covars)
      if (v %in% names(d) && any(is.na(d[[v]])))
        d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm=TRUE)

    d_ps <- d[complete.cases(d[, ps_covars]),]

    # All-patient PS + weights
    ps_fit <- glm(ps_formula, data=d_ps, family=binomial())
    d_ps$ps <- pmax(pmin(fitted(ps_fit), 0.99), 0.01)
    d_ps$ow <- ifelse(d_ps$mg_supp==1, 1-d_ps$ps, d_ps$ps)
    prev <- mean(d_ps$mg_supp)
    d_ps$iptw <- ifelse(d_ps$mg_supp==1, prev/d_ps$ps, (1-prev)/(1-d_ps$ps))
    q01 <- quantile(d_ps$iptw, 0.01); q99 <- quantile(d_ps$iptw, 0.99)
    d_ps$iptw <- pmax(pmin(d_ps$iptw, q99), q01)

    cl <- if (has_cluster && "hospitalid" %in% names(d_ps)) d_ps$hospitalid else NULL

    # AC setup
    if ("ac_group" %in% names(d_ps)) {
      d_ac <- d_ps[d_ps$ac_group %in% c("mg_k","k_only"),]
      d_ac$ac_trt <- as.integer(d_ac$ac_group == "mg_k")
      if (nrow(d_ac) > 50 && sum(d_ac$ac_trt) > 10) {
        ac_fml <- as.formula(paste("ac_trt ~", paste(ps_covars, collapse="+")))
        tryCatch({
          ac_fit <- glm(ac_fml, data=d_ac, family=binomial())
          d_ac$ac_ps <- pmax(pmin(fitted(ac_fit), 0.99), 0.01)
          d_ac$ac_ow <- ifelse(d_ac$ac_trt==1, 1-d_ac$ac_ps, d_ac$ac_ps)
        }, error=function(e) { d_ac <<- NULL })
      } else { d_ac <- NULL }
    } else { d_ac <- NULL }

    # Run each outcome
    for (oc in outcomes) {
      oc_name <- oc$name
      oc_var  <- oc$outcome
      oc_trt  <- oc$trt

      # Select population
      if (oc_trt == "ac_trt") {
        if (is.null(d_ac)) next
        dd <- d_ac
        wt <- dd$ac_ow
        cl_use <- NULL
      } else if (grepl("iptw", oc_name)) {
        dd <- d_ps; wt <- dd$iptw; cl_use <- cl
      } else {
        dd <- d_ps; wt <- dd$ow; cl_use <- cl
      }

      if (!oc_var %in% names(dd)) next

      # Filter NA (for POAF: exclude preexisting AF)
      if (oc$filter_na) dd <- dd[!is.na(dd[[oc_var]]),]
      if (sum(dd[[oc_var]], na.rm=TRUE) < 5) next

      tryCatch({
        dd$.w <- if (oc_trt == "ac_trt") dd$ac_ow else wt
        fml <- as.formula(paste(oc_var, "~", oc_trt))
        fit <- glm(fml, data=dd, weights=.w, family=quasibinomial())
        vc <- if (!is.null(cl_use) && length(unique(cl_use[seq_len(nrow(dd))])) > 1)
          vcovCL(fit, cluster=cl_use[seq_len(nrow(dd))]) else vcovHC(fit, type="HC1")
        ct <- coeftest(fit, vcov.=vc)
        stores[[oc_name]]$est <- c(stores[[oc_name]]$est, ct[2,1])
        stores[[oc_name]]$se  <- c(stores[[oc_name]]$se,  ct[2,2])
      }, error=function(e) {})
    }
  }

  # Pool results
  out_rows <- list()
  for (nm in names(stores)) {
    s <- stores[[nm]]
    if (length(s$est) == 0) next
    if (length(s$est) == 1) {
      r <- list(or=exp(s$est), lo=exp(s$est-1.96*s$se),
                hi=exp(s$est+1.96*s$se), p=2*pnorm(-abs(s$est/s$se)))
    } else {
      r <- pool_simple(s$est, s$se)
    }
    sig <- ifelse(r$p < 0.05, " *", "")
    cat(sprintf("  %-18s OR=%.3f (%.3f–%.3f) P=%.4f%s\n",
                nm, r$or, r$lo, r$hi, r$p, sig))
    out_rows[[length(out_rows)+1]] <- data.frame(
      db=db_name, analysis=nm,
      or=round(r$or,3), lo=round(r$lo,3), hi=round(r$hi,3),
      p=round(r$p,4), stringsAsFactors=FALSE)
  }
  do.call(rbind, out_rows)
}

# ── Run ───────────────────────────────────────────────────────────────
res_e <- run_additional(dat_e, "eICU", has_cluster=TRUE)
res_m <- run_additional(dat_m, "MIMIC", has_cluster=FALSE)

all_res <- rbind(res_e, res_m)
write.csv(all_res, file.path(RESULTS, "02b_additional_outcomes.csv"), row.names=FALSE)
cat(sprintf("\n✓ Saved: results/02b_additional_outcomes.csv (%d rows)\n", nrow(all_res)))

cat("\nKey results for manuscript:\n")
cat("  POAF: exploratory outcome (AF prevention is known Mg mechanism)\n")
cat("  AKI-or-death composite: competing risk sensitivity\n")
cat("    If composite ≈ AKI-only → competing risk is minimal\n")
cat("    If composite diverges → need cause-specific hazard\n")
cat("  Transfusion negative control (complexity-specific):\n")
cat("    ow_transfusion: expect confounded (OR<1 = Mg pts had simpler surgery)\n")
cat("    ac_transfusion: expect NULL (AC controls complexity)\n")
cat("    If ow confounded + ac null → design works ✓\n")
