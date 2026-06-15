#!/usr/bin/env Rscript
# 07c_mice_stability.R — MICE stability across m=5, 10, 20
# Fills eTable 2. Reuses the PS + OW machinery from 02_analysis.R.
# Run: Rscript 07c_mice_stability.R   (takes ~5 min on Tempest)

suppressPackageStartupMessages({
  library(mice); library(sandwich); library(lmtest)
})
RESULTS <- path.expand("~/mg_aki/results")
SEED <- 42

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
  if ("first_lactate" %in% names(d)) {
    d$lactate_missing <- as.integer(is.na(d$first_lactate))
    d$first_lactate[is.na(d$first_lactate)] <- median(d$first_lactate, na.rm=TRUE)
  }
  d
}

pool_simple <- function(ests, ses) {
  m <- length(ests); qbar <- mean(ests); ubar <- mean(ses^2)
  b <- var(ests); tv <- ubar + (1+1/m)*b; se <- sqrt(tv)
  list(or=exp(qbar), lo=exp(qbar-1.96*se), hi=exp(qbar+1.96*se),
       p=2*pnorm(-abs(qbar/se)))
}

wglm <- function(formula, data, w, cluster=NULL) {
  data$.w <- w
  fit <- glm(formula, data=data, weights=.w, family=quasibinomial())
  vc <- if (!is.null(cluster) && length(unique(cluster))>1)
    vcovCL(fit, cluster=cluster) else vcovHC(fit, type="HC1")
  ct <- coeftest(fit, vcov.=vc)
  tr <- which(rownames(ct)=="mg_supp"); if (length(tr)==0) tr <- 2
  list(logOR=ct[tr,1], se=ct[tr,2])
}

run_m <- function(dat, m_val, db_name, has_cluster) {
  cat(sprintf("  %s m=%d...\n", db_name, m_val))
  ps_covars <- intersect(c("age","is_female","bmi",
    "surg_cabg","surg_valve","surg_combined",
    "heart_failure","hypertension","diabetes","ckd","copd","pvd","stroke","liver_disease",
    "baseline_creatinine","egfr","loop_diuretics","nsaids","acei_arb","ppi",
    "beta_blockers","steroids","antiarrhythmics",
    "first_potassium","first_calcium","first_heartrate",
    "vasopressor_6h","first_mg_value","first_lactate","lactate_missing"), names(dat))

  mice_vars <- intersect(c("bmi","first_heartrate","first_calcium","first_potassium"), names(dat))
  any_missing <- any(sapply(mice_vars, function(v) sum(is.na(dat[[v]]))>0))

  imp_preds <- unique(c(mice_vars,"age","is_female","baseline_creatinine","mg_supp","aki_kdigo1"))
  imp_preds <- imp_preds[imp_preds %in% names(dat)]

  if (any_missing) {
    imp <- mice(dat[, imp_preds], m=m_val, method="pmm", seed=SEED, printFlag=FALSE, maxit=10)
  }

  ps_formula <- as.formula(paste("mg_supp ~", paste(ps_covars, collapse="+")))

  analyses <- list(
    ow_aki1=list(outcome="aki_kdigo1", weight="ow"),
    iptw_aki1=list(outcome="aki_kdigo1", weight="iptw"),
    ac_aki1=list(outcome="aki_kdigo1", weight="ac_ow")
  )

  results <- list()
  for (aname in names(analyses)) {
    ests <- ses <- numeric(0)
    for (i in seq_len(m_val)) {
      d <- dat
      if (any_missing) {
        imputed <- complete(imp, i)
        for (v in mice_vars) if (v %in% names(imputed)) d[[v]] <- imputed[[v]]
      }
      d_ps <- d[complete.cases(d[, ps_covars]),]
      ps_fit <- glm(ps_formula, data=d_ps, family=binomial())
      d_ps$ps <- pmax(pmin(fitted(ps_fit), 0.99), 0.01)

      if (aname == "ac_aki1" && "ac_group" %in% names(d_ps)) {
        d_ac <- d_ps[d_ps$ac_group %in% c("mg_k","k_only"),]
        d_ac$ac_trt <- as.integer(d_ac$ac_group=="mg_k")
        ac_fml <- as.formula(paste("ac_trt ~", paste(ps_covars, collapse="+")))
        tryCatch({
          ac_fit <- glm(ac_fml, data=d_ac, family=binomial())
          d_ac$ac_ps <- pmax(pmin(fitted(ac_fit), 0.99), 0.01)
          d_ac$ac_ow <- ifelse(d_ac$ac_trt==1, 1-d_ac$ac_ps, d_ac$ac_ps)
          d_ac$.w <- d_ac$ac_ow
          fit <- glm(aki_kdigo1 ~ ac_trt, data=d_ac, weights=.w, family=quasibinomial())
          ct <- coeftest(fit, vcov.=vcovHC(fit, type="HC1"))
          ests <- c(ests, ct[2,1]); ses <- c(ses, ct[2,2])
        }, error=function(e) {})
      } else {
        d_ps$ow <- ifelse(d_ps$mg_supp==1, 1-d_ps$ps, d_ps$ps)
        prev <- mean(d_ps$mg_supp)
        d_ps$iptw <- ifelse(d_ps$mg_supp==1, prev/d_ps$ps, (1-prev)/(1-d_ps$ps))
        q01 <- quantile(d_ps$iptw, 0.01); q99 <- quantile(d_ps$iptw, 0.99)
        d_ps$iptw <- pmax(pmin(d_ps$iptw, q99), q01)

        w <- if (aname=="iptw_aki1") d_ps$iptw else d_ps$ow
        cl <- if (has_cluster && "hospitalid" %in% names(d_ps)) d_ps$hospitalid else NULL
        res <- wglm(aki_kdigo1 ~ mg_supp, d_ps, w, cl)
        ests <- c(ests, res$logOR); ses <- c(ses, res$se)
      }
    }
    if (length(ests)>0) {
      r <- pool_simple(ests, ses)
      results[[aname]] <- data.frame(db=db_name, analysis=aname, m=m_val,
        or=round(r$or,3), lo=round(r$lo,3), hi=round(r$hi,3), p=round(r$p,4))
    }
  }
  do.call(rbind, results)
}

cat(strrep("=",60), "\n07c: MICE Stability (m=10, 20)\n", strrep("=",60), "\n")
dat_e <- stdz(read.csv(file.path(RESULTS,"01_analysis_a_cohort.csv"), stringsAsFactors=FALSE))
dat_m <- stdz(read.csv(file.path(RESULTS,"04_mimic_cohort.csv"), stringsAsFactors=FALSE))

all_res <- list()
for (m_val in c(10, 20)) {
  all_res[[length(all_res)+1]] <- run_m(dat_e, m_val, "eICU", TRUE)
  all_res[[length(all_res)+1]] <- run_m(dat_m, m_val, "MIMIC", FALSE)
}
out <- do.call(rbind, all_res)

# Append m=5 from existing results
res5 <- read.csv(file.path(RESULTS,"02_results.csv"), stringsAsFactors=FALSE)
if ("m" %in% names(res5)) {
  m5 <- res5[!is.na(res5$m) & res5$m==min(res5$m,na.rm=TRUE) &
             res5$analysis %in% c("ow_aki1","iptw_aki1","ac_aki1") &
             res5$db != "Pooled",]
  if (nrow(m5)>0) {
    m5$m <- 5
    m5 <- m5[, intersect(names(out), names(m5))]
    out <- rbind(m5, out)
  }
}

outpath <- file.path(RESULTS, "07c_mice_stability.csv")
write.csv(out, outpath, row.names=FALSE)
cat(sprintf("\nSaved: %s (%d rows)\n", outpath, nrow(out)))
print(out)
