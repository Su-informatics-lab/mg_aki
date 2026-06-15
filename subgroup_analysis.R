#!/usr/bin/env Rscript
# ============================================================================
# subgroup_analysis.R — Downstream subgroup + safety analyses
#
# Approach: PS + OW estimated ONCE on full sample, then subgroup analyses
# use existing weights (standard clinical approach, per Yan).
#
# Outputs:
#   eTable 3:  AC subgroup baseline characteristics
#   eTable 4:  Treatment effect by surgery type
#   eTable 5:  Treatment effect by MDS score
#   eTable 5b: Treatment effect by age (≥60 vs <60)
#   eTable 5c: Treatment effect by baseline Mg
#   eTable 6:  Safety outcomes
#
# Run: Rscript subgroup_analysis.R
# ============================================================================

suppressPackageStartupMessages({
  library(sandwich); library(lmtest); library(tableone)
})

RESULTS <- path.expand("~/mg_aki/results")

# ── Load ──────────────────────────────────────────────────────────────────
load_cohort <- function(enriched, original) {
  f <- if (file.exists(file.path(RESULTS, enriched))) enriched else original
  cat(sprintf("  Loading %s\n", f))
  read.csv(file.path(RESULTS, f), stringsAsFactors = FALSE)
}

cat("Loading cohorts...\n")
dat_e <- load_cohort("01_analysis_a_cohort_enriched.csv", "01_analysis_a_cohort.csv")
dat_m <- load_cohort("04_mimic_cohort_enriched.csv", "04_mimic_cohort.csv")

# ── Standardize ───────────────────────────────────────────────────────────
stdz <- function(d) {
  rmap <- c(mg_supplementation="mg_supp", hosp_mortality="hospital_mortality",
            age_num="age", baseline_egfr="egfr", baseline_cr="baseline_creatinine",
            nephrotox_loop_diuretic="loop_diuretics", nephrotox_ppi="ppi",
            has_betablocker="beta_blockers", has_steroid="steroids",
            first_hr="first_heartrate", has_vasopressor="vasopressor_6h",
            nephrotox_nsaid="nsaids", nephrotox_acei_arb="acei_arb",
            preop_antiarrhythmic="antiarrhythmics",
            first_k_value="first_potassium", first_ca_value="first_calcium",
            hx_chf="heart_failure", hx_hypertension="hypertension",
            hx_diabetes="diabetes", hx_ckd="ckd", hx_copd="copd",
            hx_pvd="pvd", hx_stroke="stroke", hx_liver="liver_disease",
            nc_fracture="fracture", neuro_encephalopathy="encephalopathy")
  for (old in names(rmap)) {
    new <- rmap[[old]]
    if (old %in% names(d) && !new %in% names(d)) names(d)[names(d)==old] <- new
  }
  if (is.character(d$age)) { d$age <- suppressWarnings(as.numeric(d$age)); d$age[is.na(d$age)] <- 90 }
  d
}
dat_e <- stdz(dat_e); dat_m <- stdz(dat_m)

# ── Helpers ───────────────────────────────────────────────────────────────
fmt_or <- function(r) sprintf("OR %.2f (%.2f-%.2f) P=%.3f", r$or, r$lo, r$hi, r$p)
fmt_d  <- function(r) sprintf("diff %.2f (%.2f to %.2f) P=%.3f", r$est, r$lo, r$hi, r$p)

smd_w <- function(x, trt, w) {
  m1 <- weighted.mean(x[trt==1], w[trt==1], na.rm=TRUE)
  m0 <- weighted.mean(x[trt==0], w[trt==0], na.rm=TRUE)
  sp <- sqrt((var(x[trt==1],na.rm=TRUE) + var(x[trt==0],na.rm=TRUE))/2)
  if (sp < 1e-10) return(0)
  abs(m1 - m0) / sp
}

wglm_sub <- function(outcome, trt_var, d, w, cluster=NULL) {
  if (!outcome %in% names(d)) return(NULL)
  d <- d[!is.na(d[[outcome]]),]
  if (sum(d[[outcome]],na.rm=TRUE) < 5) return(NULL)
  d$.w <- w[!is.na(d[[outcome]])]  # align weights
  # Safer: re-index
  d$.w <- d[[paste0(trt_var, "_ow")]]
  if (is.null(d$.w) || all(is.na(d$.w))) return(NULL)
  tryCatch({
    fit <- glm(as.formula(paste(outcome, "~", trt_var)), data=d, weights=.w,
               family=quasibinomial())
    vc <- if (!is.null(cluster) && length(unique(cluster[!is.na(d[[outcome]])])) > 1)
      vcovCL(fit, cluster=cluster[!is.na(d[[outcome]])]) else vcovHC(fit, type="HC1")
    ct <- coeftest(fit, vcov.=vc)
    tr <- 2
    list(or=exp(ct[tr,1]), lo=exp(ct[tr,1]-1.96*ct[tr,2]),
         hi=exp(ct[tr,1]+1.96*ct[tr,2]), p=2*pnorm(-abs(ct[tr,1]/ct[tr,2])))
  }, error=function(e) NULL)
}

wlm_sub <- function(outcome, trt_var, d, cluster=NULL) {
  if (!outcome %in% names(d)) return(NULL)
  d <- d[!is.na(d[[outcome]]),]
  d$.w <- d[[paste0(trt_var, "_ow")]]
  if (is.null(d$.w)) return(NULL)
  tryCatch({
    fit <- lm(as.formula(paste(outcome, "~", trt_var)), data=d, weights=.w)
    vc <- if (!is.null(cluster) && length(unique(cluster[!is.na(d[[outcome]])])) > 1)
      vcovCL(fit, cluster=cluster[!is.na(d[[outcome]])]) else vcovHC(fit, type="HC1")
    ct <- coeftest(fit, vcov.=vc)
    tr <- 2
    list(est=ct[tr,1], lo=ct[tr,1]-1.96*ct[tr,2], hi=ct[tr,1]+1.96*ct[tr,2],
         p=2*pnorm(-abs(ct[tr,1]/ct[tr,2])))
  }, error=function(e) NULL)
}

# ============================================================================
# STEP 1: Estimate PS + OW ONCE on full sample (per database)
# ============================================================================
cat("\n", strrep("=",60), "\n")
cat("STEP 1: Full-sample PS + OW estimation\n")
cat(strrep("=",60), "\n")

fit_full_ps <- function(d, db_name) {
  ps_vars <- intersect(c("age","is_female","bmi",
    "heart_failure","hypertension","diabetes","ckd","copd","pvd","stroke","liver_disease",
    "baseline_creatinine","egfr","loop_diuretics","nsaids","acei_arb","ppi",
    "beta_blockers","steroids","antiarrhythmics",
    "first_potassium","first_calcium","first_heartrate",
    "vasopressor_6h","first_mg_value"), names(d))
  if ("surgery_type" %in% names(d)) {
    d$s_cabg <- as.integer(d$surgery_type=="cabg")
    d$s_valve <- as.integer(d$surgery_type=="valve")
    d$s_combined <- as.integer(d$surgery_type=="combined")
    ps_vars <- c(ps_vars, "s_cabg","s_valve","s_combined")
  }
  # Median-impute for PS estimation
  for (v in ps_vars) if (any(is.na(d[[v]])))
    d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm=TRUE)

  # All-patient PS
  fml <- as.formula(paste("mg_supp ~", paste(ps_vars, collapse="+")))
  d <- d[complete.cases(d[,ps_vars]),]
  ps_fit <- glm(fml, data=d, family=binomial())
  d$ps_all <- pmax(pmin(fitted(ps_fit), 0.99), 0.01)
  d$mg_supp_ow <- ifelse(d$mg_supp==1, 1-d$ps_all, d$ps_all)

  smds <- sapply(ps_vars, function(v) if(is.numeric(d[[v]])) smd_w(d[[v]], d$mg_supp, d$mg_supp_ow) else NA)
  cat(sprintf("  %s all-patient: N=%d, max SMD=%.4f\n", db_name, nrow(d), max(smds,na.rm=TRUE)))

  # AC PS (within K+-repleted)
  if ("ac_group" %in% names(d)) {
    d$ac_trt <- NA_integer_
    d$ac_trt[d$ac_group == "mg_k"] <- 1L
    d$ac_trt[d$ac_group == "k_only"] <- 0L
    d_ac <- d[!is.na(d$ac_trt),]
    ac_fml <- as.formula(paste("ac_trt ~", paste(ps_vars, collapse="+")))
    tryCatch({
      ac_fit <- glm(ac_fml, data=d_ac, family=binomial())
      d_ac$ps_ac <- pmax(pmin(fitted(ac_fit), 0.99), 0.01)
      d_ac$ac_trt_ow <- ifelse(d_ac$ac_trt==1, 1-d_ac$ps_ac, d_ac$ps_ac)
      # Write AC weights back
      d$ac_trt_ow <- NA_real_
      d$ac_trt_ow[match(rownames(d_ac), rownames(d))] <- d_ac$ac_trt_ow
      d$ps_ac <- NA_real_
      d$ps_ac[match(rownames(d_ac), rownames(d))] <- d_ac$ps_ac
      ac_smds <- sapply(ps_vars, function(v) if(is.numeric(d_ac[[v]])) smd_w(d_ac[[v]], d_ac$ac_trt, d_ac$ac_trt_ow) else NA)
      cat(sprintf("  %s AC: N=%d (trt=%d), max SMD=%.4f, Mg SMD=%.4f\n",
                  db_name, nrow(d_ac), sum(d_ac$ac_trt), max(ac_smds,na.rm=TRUE),
                  ac_smds["first_mg_value"]))
    }, error=function(e) cat(sprintf("  AC PS failed: %s\n", e$message)))
  }

  # MDS score
  d$mds <- 0L
  if ("loop_diuretics" %in% names(d)) d$mds <- d$mds + d$loop_diuretics
  if ("ppi" %in% names(d)) d$mds <- d$mds + d$ppi
  if ("alcohol_history" %in% names(d)) d$mds <- d$mds + d$alcohol_history
  if ("egfr" %in% names(d)) d$mds <- d$mds + ifelse(!is.na(d$egfr) & d$egfr<60, 2L,
                                                      ifelse(!is.na(d$egfr) & d$egfr<90, 1L, 0L))
  d$mds_high <- as.integer(d$mds >= 2)
  cat(sprintf("  %s MDS >=2: %d/%d (%.1f%%)\n", db_name, sum(d$mds_high), nrow(d), 100*mean(d$mds_high)))

  # ICU LOS
  if (!"icu_los_h" %in% names(d)) {
    if ("unitdischargeoffset" %in% names(d)) d$icu_los_h <- d$unitdischargeoffset / 60
    else if ("outtime" %in% names(d) && "intime" %in% names(d)) {
      d$icu_los_h <- as.numeric(difftime(as.POSIXct(d$outtime), as.POSIXct(d$intime), units="hours"))
    }
  }

  d
}

dat_e <- fit_full_ps(dat_e, "eICU")
dat_m <- fit_full_ps(dat_m, "MIMIC")

# ============================================================================
# HELPER: Run downstream subgroup (uses existing weights, checks balance)
# ============================================================================
run_subgroup <- function(d, subset_idx, trt_var, outcome, db_name, cluster=NULL) {
  d_sub <- d[subset_idx,]
  n_total <- nrow(d_sub)
  n_trt <- sum(d_sub[[trt_var]], na.rm=TRUE)
  if (n_trt < 10 || (n_total - n_trt) < 10) return(NULL)

  wt_col <- paste0(trt_var, "_ow")
  if (!wt_col %in% names(d_sub) || all(is.na(d_sub[[wt_col]]))) return(NULL)

  # Check balance within subgroup
  d_bal <- d_sub[!is.na(d_sub[[wt_col]]),]
  check_vars <- intersect(c("age","is_female","bmi","baseline_creatinine","egfr",
                             "first_mg_value","first_potassium"), names(d_bal))
  smds <- sapply(check_vars, function(v)
    if (is.numeric(d_bal[[v]])) smd_w(d_bal[[v]], d_bal[[trt_var]], d_bal[[wt_col]]) else NA)
  max_smd <- max(smds, na.rm=TRUE)

  # Outcome model
  d_sub$.w <- d_sub[[wt_col]]
  d_sub <- d_sub[!is.na(d_sub$.w) & !is.na(d_sub[[outcome]]),]
  if (sum(d_sub[[outcome]]) < 5) return(NULL)

  cl <- if (!is.null(cluster)) cluster[subset_idx][!is.na(d[subset_idx, wt_col]) & !is.na(d[subset_idx, outcome])] else NULL

  tryCatch({
    fit <- glm(as.formula(paste(outcome, "~", trt_var)), data=d_sub, weights=.w,
               family=quasibinomial())
    vc <- if (!is.null(cl) && length(unique(cl))>1) vcovCL(fit, cluster=cl) else vcovHC(fit, type="HC1")
    ct <- coeftest(fit, vcov.=vc)
    list(or=exp(ct[2,1]), lo=exp(ct[2,1]-1.96*ct[2,2]), hi=exp(ct[2,1]+1.96*ct[2,2]),
         p=2*pnorm(-abs(ct[2,1]/ct[2,2])), n=n_total, n_trt=n_trt, max_smd=max_smd)
  }, error=function(e) NULL)
}

# ============================================================================
# eTABLE 3: AC Subgroup Baseline
# ============================================================================
cat("\n", strrep("=",60), "\n")
cat("eTABLE 3: AC Subgroup Baseline\n")
cat(strrep("=",60), "\n")

for (db in list(list(d=dat_e, nm="eICU"), list(d=dat_m, nm="MIMIC"))) {
  d <- db$d
  if (!"ac_trt" %in% names(d)) next
  d_ac <- d[!is.na(d$ac_trt),]
  d_ac$trt_label <- ifelse(d_ac$ac_trt==1, "Mg+K", "K-only")
  vars <- intersect(c("age","is_female","bmi","surgery_type",
    "heart_failure","hypertension","diabetes","ckd","baseline_creatinine","egfr",
    "loop_diuretics","ppi","first_mg_value","mds","aki_kdigo1","hospital_mortality"), names(d_ac))
  cat(sprintf("\n  %s AC (N=%d):\n", db$nm, nrow(d_ac)))
  tryCatch({
    t1 <- CreateTableOne(vars=vars, strata="trt_label", data=d_ac, test=FALSE)
    print(t1, smd=TRUE, printToggle=TRUE)
  }, error=function(e) cat(sprintf("  Failed: %s\n", e$message)))
}

# ============================================================================
# eTABLES 4-5c: All subgroup analyses (downstream, existing weights)
# ============================================================================

subgroups <- list(
  # eTable 4: Surgery type
  list(tbl="4", var="surgery_type", cuts=list(
    list(label="CABG", expr=quote(surgery_type=="cabg")),
    list(label="Valve", expr=quote(surgery_type=="valve")),
    list(label="Combined", expr=quote(surgery_type=="combined")),
    list(label="Other cardiac", expr=quote(surgery_type=="other_cardiac"))
  )),
  # eTable 5: MDS
  list(tbl="5", var="mds_high", cuts=list(
    list(label="MDS 0-1", expr=quote(mds_high==0)),
    list(label="MDS >=2", expr=quote(mds_high==1))
  )),
  # eTable 5b: Age
  list(tbl="5b", var="age", cuts=list(
    list(label="Age <60", expr=quote(age < 60)),
    list(label="Age >=60", expr=quote(age >= 60))
  )),
  # eTable 5c: Baseline Mg
  list(tbl="5c", var="first_mg_value", cuts=list(
    list(label="Mg <2.0 (hypo)", expr=quote(first_mg_value < 2.0)),
    list(label="Mg >=2.0", expr=quote(first_mg_value >= 2.0))
  )),
  # eTable 5d: Individual MDS components (Yan request)
  list(tbl="5d", var="depletion_risk", cuts=list(
    list(label="PPI user", expr=quote(ppi == 1)),
    list(label="No PPI", expr=quote(ppi == 0)),
    list(label="Loop diuretic", expr=quote(loop_diuretics == 1)),
    list(label="No loop diuretic", expr=quote(loop_diuretics == 0)),
    list(label="eGFR <60", expr=quote(egfr < 60)),
    list(label="eGFR >=60", expr=quote(egfr >= 60))
  ))
)

all_sub_rows <- list()

for (sg in subgroups) {
  cat(sprintf("\n%s\neTABLE %s: Subgroup by %s\n%s\n", strrep("=",60), sg$tbl, sg$var, strrep("=",60)))
  for (cut in sg$cuts) {
    for (db in list(list(d=dat_e, nm="eICU", cl=TRUE), list(d=dat_m, nm="MIMIC", cl=FALSE))) {
      idx <- tryCatch(eval(cut$expr, db$d), error=function(e) rep(FALSE, nrow(db$d)))
      idx[is.na(idx)] <- FALSE
      cl <- if (db$cl && "hospitalid" %in% names(db$d)) db$d$hospitalid else NULL

      # All-patient OW
      r <- run_subgroup(db$d, idx, "mg_supp", "aki_kdigo1", db$nm, cl)
      if (!is.null(r)) {
        cat(sprintf("  %s %-20s OW:  N=%d trt=%d %s [bal %.3f]\n",
                    db$nm, cut$label, r$n, r$n_trt, fmt_or(r), r$max_smd))
        all_sub_rows[[length(all_sub_rows)+1]] <- data.frame(
          etable=sg$tbl, db=db$nm, subgroup=cut$label, analysis="all_ow",
          n=r$n, n_trt=r$n_trt, or=round(r$or,3), lo=round(r$lo,3),
          hi=round(r$hi,3), p=round(r$p,4), max_smd=round(r$max_smd,4))
      }

      # AC OW (downstream)
      r_ac <- run_subgroup(db$d, idx & !is.na(db$d$ac_trt), "ac_trt", "aki_kdigo1", db$nm, NULL)
      if (!is.null(r_ac)) {
        cat(sprintf("  %s %-20s AC:  N=%d trt=%d %s [bal %.3f]\n",
                    db$nm, cut$label, r_ac$n, r_ac$n_trt, fmt_or(r_ac), r_ac$max_smd))
        all_sub_rows[[length(all_sub_rows)+1]] <- data.frame(
          etable=sg$tbl, db=db$nm, subgroup=cut$label, analysis="ac_ow",
          n=r_ac$n, n_trt=r_ac$n_trt, or=round(r_ac$or,3), lo=round(r_ac$lo,3),
          hi=round(r_ac$hi,3), p=round(r_ac$p,4), max_smd=round(r_ac$max_smd,4))
      }
    }
  }
}

if (length(all_sub_rows) > 0) {
  sub_df <- do.call(rbind, all_sub_rows)
  write.csv(sub_df, file.path(RESULTS, "etables_4_5_subgroups.csv"), row.names=FALSE)
  cat(sprintf("\n✓ Saved etables_4_5_subgroups.csv (%d rows)\n", nrow(sub_df)))
}

# ============================================================================
# eTABLE 6: Safety Outcomes
# ============================================================================
cat(sprintf("\n%s\neTABLE 6: Safety Outcomes\n%s\n", strrep("=",60), strrep("=",60)))

t6_rows <- list()
for (db in list(list(d=dat_e, nm="eICU", cl=TRUE), list(d=dat_m, nm="MIMIC", cl=FALSE))) {
  d <- db$d
  cl <- if (db$cl && "hospitalid" %in% names(d)) d$hospitalid else NULL
  cat(sprintf("\n  %s (N=%d):\n", db$nm, nrow(d)))

  # Continuous outcomes (OW-weighted linear model)
  for (oc in c("icu_los_h", "vent_duration_h", "post_hr_6h")) {
    if (!oc %in% names(d) || sum(!is.na(d[[oc]])) < 50) next
    d_oc <- d[!is.na(d[[oc]]) & !is.na(d$mg_supp_ow),]
    d_oc$.w <- d_oc$mg_supp_ow
    tryCatch({
      fit <- lm(as.formula(paste(oc, "~ mg_supp")), data=d_oc, weights=.w)
      vc <- if (!is.null(cl)) vcovCL(fit, cluster=cl[!is.na(d[[oc]]) & !is.na(d$mg_supp_ow)]) else vcovHC(fit, type="HC1")
      ct <- coeftest(fit, vcov.=vc)
      r <- list(est=ct[2,1], lo=ct[2,1]-1.96*ct[2,2], hi=ct[2,1]+1.96*ct[2,2],
                p=2*pnorm(-abs(ct[2,1]/ct[2,2])))
      cat(sprintf("    %-20s %s\n", oc, fmt_d(r)))
      # Descriptive
      trt_v <- d_oc[[oc]][d_oc$mg_supp==1]; ctrl_v <- d_oc[[oc]][d_oc$mg_supp==0]
      cat(sprintf("      Trt: median %.1f, Ctrl: median %.1f\n", median(trt_v), median(ctrl_v)))
      t6_rows[[length(t6_rows)+1]] <- data.frame(
        db=db$nm, outcome=oc, est=round(r$est,2), lo=round(r$lo,2),
        hi=round(r$hi,2), p=round(r$p,4), type="continuous")
    }, error=function(e) cat(sprintf("    %-20s failed: %s\n", oc, e$message)))
  }

  # Max post-treatment Mg (descriptive safety)
  if ("max_posttreat_mg" %in% names(d)) {
    d_mg <- d[!is.na(d$max_posttreat_mg),]
    trt_mg <- d_mg$max_posttreat_mg[d_mg$mg_supp==1]
    ctrl_mg <- d_mg$max_posttreat_mg[d_mg$mg_supp==0]
    cat(sprintf("    Max post-Mg:       Trt median %.2f (max %.2f), Ctrl median %.2f (max %.2f)\n",
                median(trt_mg), max(trt_mg), median(ctrl_mg), max(ctrl_mg)))
    cat(sprintf("      >4.8: Trt %d/%d (%.1f%%), Ctrl %d/%d (%.1f%%)\n",
                sum(trt_mg>4.8), length(trt_mg), 100*mean(trt_mg>4.8),
                sum(ctrl_mg>4.8), length(ctrl_mg), 100*mean(ctrl_mg>4.8)))
    t6_rows[[length(t6_rows)+1]] <- data.frame(
      db=db$nm, outcome="max_mg_gt4.8_trt_pct", est=round(100*mean(trt_mg>4.8),1),
      lo=length(trt_mg), hi=sum(trt_mg>4.8), p=NA, type="descriptive")
  }

  # Binary outcomes (OW-weighted)
  for (oc in c("has_iabp", "has_ecmo")) {
    if (!oc %in% names(d) || sum(d[[oc]],na.rm=TRUE) < 3) next
    r <- run_subgroup(d, rep(TRUE, nrow(d)), "mg_supp", oc, db$nm, cl)
    if (!is.null(r)) {
      cat(sprintf("    %-20s %s [events=%d]\n", oc, fmt_or(r), sum(d[[oc]],na.rm=TRUE)))
      t6_rows[[length(t6_rows)+1]] <- data.frame(
        db=db$nm, outcome=oc, est=round(r$or,3), lo=round(r$lo,3),
        hi=round(r$hi,3), p=round(r$p,4), type="binary")
    }
  }
}

if (length(t6_rows) > 0) {
  t6 <- do.call(rbind, t6_rows)
  write.csv(t6, file.path(RESULTS, "etable6_safety.csv"), row.names=FALSE)
  cat("\n✓ Saved etable6_safety.csv\n")
}

# ============================================================================
# SUMMARY
# ============================================================================
cat(sprintf("\n%s\nSUMMARY\n%s\n", strrep("=",60), strrep("=",60)))
cat("  All subgroup analyses use full-sample OW weights (downstream).\n")
cat("  Balance checked within each subgroup (max_smd column).\n\n")
cat("  Files:\n")
cat("    etables_4_5_subgroups.csv — surgery, MDS, age, baseline Mg, PPI/diuretic/eGFR\n")
cat("    etable6_safety.csv — ICU LOS, vent, HR, max Mg, IABP/ECMO\n\n")
cat("  Key predictions:\n")
cat("    Mg <2.0 effect > Mg >=2.0 → 'repletion benefits the depleted'\n")
cat("    Age >=60 effect > Age <60 → age-related Mg depletion\n")
cat("    MDS >=2 effect > MDS 0-1 → composite risk guides therapy\n")
cat("    Safety null across the board → 'no harm at these doses'\n")
cat("\nsubgroup_analysis.R COMPLETE\n")
