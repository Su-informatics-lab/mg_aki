#!/usr/bin/env Rscript
# ============================================================================
# subgroup_analysis.R — Subgroup + Safety eTables
#
#   eTable 3: AC subgroup baseline characteristics
#   eTable 4: Treatment effect by surgery type
#   eTable 5: Treatment effect by MDS score (modified, 0-4)
#   eTable 6: Safety outcomes (ICU LOS, vent duration, max Mg, HR, IABP/ECMO)
#
# Run: Rscript subgroup_analysis.R
# Input: results/*_enriched.csv (from subgroup_extract.py)
#        Falls back to results/01_analysis_a_cohort.csv if enriched not found
# ============================================================================

suppressPackageStartupMessages({
  library(sandwich); library(lmtest); library(tableone)
})

RESULTS <- path.expand("~/mg_aki/results")

# ── Load (prefer enriched, fall back to original) ─────────────────────────
load_cohort <- function(enriched, original) {
  f <- if (file.exists(file.path(RESULTS, enriched))) enriched else original
  cat(sprintf("  Loading %s\n", f))
  read.csv(file.path(RESULTS, f), stringsAsFactors = FALSE)
}

cat("Loading cohorts...\n")
dat_e <- load_cohort("01_analysis_a_cohort_enriched.csv", "01_analysis_a_cohort.csv")
dat_m <- load_cohort("04_mimic_cohort_enriched.csv", "04_mimic_cohort.csv")

# ── Standardize (same mapping as test_r2.R) ───────────────────────────────
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

# ── Compute MDS score ─────────────────────────────────────────────────────
compute_mds <- function(d) {
  d$mds <- 0L
  if ("loop_diuretics" %in% names(d)) d$mds <- d$mds + d$loop_diuretics
  if ("ppi" %in% names(d)) d$mds <- d$mds + d$ppi
  if ("alcohol_history" %in% names(d)) d$mds <- d$mds + d$alcohol_history
  if ("egfr" %in% names(d)) {
    d$mds <- d$mds + ifelse(!is.na(d$egfr) & d$egfr < 60, 2L,
                            ifelse(!is.na(d$egfr) & d$egfr < 90, 1L, 0L))
  }
  d$mds_high <- as.integer(d$mds >= 2)
  cat(sprintf("    MDS: mean=%.1f, >=2: %d/%d (%.1f%%)\n",
              mean(d$mds), sum(d$mds_high), nrow(d), 100*mean(d$mds_high)))
  d
}

cat("\neICU MDS:\n"); dat_e <- compute_mds(dat_e)
cat("MIMIC MDS:\n"); dat_m <- compute_mds(dat_m)

# ── Weighted GLM helper ───────────────────────────────────────────────────
wglm <- function(fml, d, w, cluster=NULL) {
  d$.w <- w
  fit <- glm(fml, data=d, weights=.w, family=quasibinomial())
  vc <- if (!is.null(cluster) && length(unique(cluster))>1) vcovCL(fit, cluster=cluster)
        else vcovHC(fit, type="HC1")
  ct <- coeftest(fit, vcov.=vc)
  tr <- 2  # treatment coefficient row
  list(logOR=ct[tr,1], se=ct[tr,2], or=exp(ct[tr,1]),
       lo=exp(ct[tr,1]-1.96*ct[tr,2]), hi=exp(ct[tr,1]+1.96*ct[tr,2]),
       p=2*pnorm(-abs(ct[tr,1]/ct[tr,2])))
}

wlm <- function(fml, d, w, cluster=NULL) {
  d$.w <- w
  fit <- lm(fml, data=d, weights=.w)
  vc <- if (!is.null(cluster) && length(unique(cluster))>1) vcovCL(fit, cluster=cluster)
        else vcovHC(fit, type="HC1")
  ct <- coeftest(fit, vcov.=vc)
  tr <- 2
  list(est=ct[tr,1], se=ct[tr,2],
       lo=ct[tr,1]-1.96*ct[tr,2], hi=ct[tr,1]+1.96*ct[tr,2],
       p=2*pnorm(-abs(ct[tr,1]/ct[tr,2])))
}

fmt_or <- function(r) sprintf("OR %.2f (%.2f-%.2f) P=%.3f", r$or, r$lo, r$hi, r$p)
fmt_d  <- function(r) sprintf("diff %.2f (%.2f to %.2f) P=%.3f", r$est, r$lo, r$hi, r$p)

# ── Quick PS + OW for a given subset ──────────────────────────────────────
run_ow <- function(d, trt_var, outcome, db_name, cluster=NULL) {
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
  # Median-impute remaining NAs for quick subgroup analysis
  for (v in ps_vars) if (v %in% names(d) && any(is.na(d[[v]])))
    d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm=TRUE)

  fml <- as.formula(paste(trt_var, "~", paste(ps_vars, collapse="+")))
  d <- d[complete.cases(d[,c(ps_vars, trt_var, outcome)]),]
  if (sum(d[[trt_var]]) < 10 || sum(d[[trt_var]]==0) < 10) return(NULL)
  tryCatch({
    ps_fit <- glm(fml, data=d, family=binomial())
    d$ps <- pmax(pmin(fitted(ps_fit), 0.99), 0.01)
    d$ow <- ifelse(d[[trt_var]]==1, 1-d$ps, d$ps)
    wglm(as.formula(paste(outcome, "~", trt_var)), d, d$ow, cluster)
  }, error=function(e) NULL)
}

# ============================================================================
# eTABLE 3: AC Subgroup Baseline Characteristics
# ============================================================================
cat("\n", strrep("=",60), "\n")
cat("eTABLE 3: AC Subgroup Baseline\n")
cat(strrep("=",60), "\n")

for (db in list(list(d=dat_e, nm="eICU"), list(d=dat_m, nm="MIMIC"))) {
  d <- db$d
  if (!"ac_group" %in% names(d)) next
  d_ac <- d[d$ac_group %in% c("mg_k","k_only"),]
  d_ac$trt_label <- ifelse(d_ac$mg_supp==1, "Mg+K", "K-only")
  vars <- intersect(c("age","is_female","bmi","surgery_type",
    "heart_failure","hypertension","diabetes","ckd","baseline_creatinine","egfr",
    "loop_diuretics","ppi","first_mg_value","aki_kdigo1","hospital_mortality"), names(d_ac))
  cat(sprintf("\n  %s AC subgroup (N=%d):\n", db$nm, nrow(d_ac)))
  tryCatch({
    t1 <- CreateTableOne(vars=vars, strata="trt_label", data=d_ac, test=FALSE)
    print(t1, smd=TRUE, printToggle=TRUE)
  }, error=function(e) cat(sprintf("  Table failed: %s\n", e$message)))
}

# ============================================================================
# eTABLE 4: Treatment Effect by Surgery Type
# ============================================================================
cat("\n", strrep("=",60), "\n")
cat("eTABLE 4: Effect by Surgery Type\n")
cat(strrep("=",60), "\n")

t4_rows <- list()
for (stype in c("cabg","valve","combined","other_cardiac")) {
  for (db in list(list(d=dat_e, nm="eICU", cl=TRUE), list(d=dat_m, nm="MIMIC", cl=FALSE))) {
    d_sub <- db$d[db$d$surgery_type == stype,]
    n_trt <- sum(d_sub$mg_supp, na.rm=TRUE)
    n_total <- nrow(d_sub)
    if (n_trt < 10 || n_total < 50) {
      cat(sprintf("  %s %s: N=%d (trt=%d) — skipped\n", db$nm, stype, n_total, n_trt))
      next
    }
    cl <- if (db$cl && "hospitalid" %in% names(d_sub)) d_sub$hospitalid else NULL
    r <- run_ow(d_sub, "mg_supp", "aki_kdigo1", db$nm, cl)
    if (!is.null(r)) {
      cat(sprintf("  %s %s (N=%d, trt=%d): %s\n", db$nm, stype, n_total, n_trt, fmt_or(r)))
      t4_rows[[length(t4_rows)+1]] <- data.frame(
        db=db$nm, surgery=stype, n=n_total, n_trt=n_trt,
        or=round(r$or,3), lo=round(r$lo,3), hi=round(r$hi,3), p=round(r$p,4))
    }
  }
}
if (length(t4_rows) > 0) {
  t4 <- do.call(rbind, t4_rows)
  write.csv(t4, file.path(RESULTS, "etable4_surgery_subgroups.csv"), row.names=FALSE)
  cat("  Saved etable4_surgery_subgroups.csv\n")
}

# ============================================================================
# eTABLE 5: Treatment Effect by MDS Score
# ============================================================================
cat("\n", strrep("=",60), "\n")
cat("eTABLE 5: Effect by MDS Score\n")
cat(strrep("=",60), "\n")

t5_rows <- list()
for (mds_level in list(list(label="MDS 0-1", val=0), list(label="MDS >=2", val=1))) {
  for (db in list(list(d=dat_e, nm="eICU", cl=TRUE), list(d=dat_m, nm="MIMIC", cl=FALSE))) {
    d_sub <- db$d[db$d$mds_high == mds_level$val,]
    n_trt <- sum(d_sub$mg_supp, na.rm=TRUE)

    # All-patient OW
    cl <- if (db$cl && "hospitalid" %in% names(d_sub)) d_sub$hospitalid else NULL
    r <- run_ow(d_sub, "mg_supp", "aki_kdigo1", db$nm, cl)
    if (!is.null(r)) {
      cat(sprintf("  %s %s (N=%d, trt=%d): %s\n", db$nm, mds_level$label, nrow(d_sub), n_trt, fmt_or(r)))
      t5_rows[[length(t5_rows)+1]] <- data.frame(
        db=db$nm, mds=mds_level$label, analysis="all_ow", n=nrow(d_sub), n_trt=n_trt,
        or=round(r$or,3), lo=round(r$lo,3), hi=round(r$hi,3), p=round(r$p,4))
    }

    # AC within subgroup
    if ("ac_group" %in% names(d_sub)) {
      d_ac <- d_sub[d_sub$ac_group %in% c("mg_k","k_only"),]
      d_ac$ac_trt <- as.integer(d_ac$ac_group == "mg_k")
      r_ac <- run_ow(d_ac, "ac_trt", "aki_kdigo1", db$nm, NULL)
      if (!is.null(r_ac)) {
        cat(sprintf("    AC: %s\n", fmt_or(r_ac)))
        t5_rows[[length(t5_rows)+1]] <- data.frame(
          db=db$nm, mds=mds_level$label, analysis="ac_ow", n=nrow(d_ac),
          n_trt=sum(d_ac$ac_trt), or=round(r_ac$or,3), lo=round(r_ac$lo,3),
          hi=round(r_ac$hi,3), p=round(r_ac$p,4))
      }
    }
  }
}
if (length(t5_rows) > 0) {
  t5 <- do.call(rbind, t5_rows)
  write.csv(t5, file.path(RESULTS, "etable5_mds_subgroups.csv"), row.names=FALSE)
  cat("  Saved etable5_mds_subgroups.csv\n")
}

# ============================================================================
# eTABLE 5b: Treatment Effect by Age (≥60 vs <60)
# ============================================================================
cat("\n", strrep("=",60), "\n")
cat("eTABLE 5b: Effect by Age\n")
cat(strrep("=",60), "\n")

t5b_rows <- list()
for (age_cut in list(list(label="Age <60", lo=-Inf, hi=60),
                     list(label="Age >=60", lo=60, hi=Inf))) {
  for (db in list(list(d=dat_e, nm="eICU", cl=TRUE), list(d=dat_m, nm="MIMIC", cl=FALSE))) {
    d_sub <- db$d[db$d$age >= age_cut$lo & db$d$age < age_cut$hi,]
    n_trt <- sum(d_sub$mg_supp, na.rm=TRUE)
    cl <- if (db$cl && "hospitalid" %in% names(d_sub)) d_sub$hospitalid else NULL
    r <- run_ow(d_sub, "mg_supp", "aki_kdigo1", db$nm, cl)
    if (!is.null(r)) {
      cat(sprintf("  %s %s (N=%d, trt=%d): %s\n", db$nm, age_cut$label, nrow(d_sub), n_trt, fmt_or(r)))
      t5b_rows[[length(t5b_rows)+1]] <- data.frame(
        db=db$nm, subgroup=age_cut$label, analysis="all_ow", n=nrow(d_sub),
        n_trt=n_trt, or=round(r$or,3), lo=round(r$lo,3), hi=round(r$hi,3), p=round(r$p,4))
    }
    # AC within age subgroup
    if ("ac_group" %in% names(d_sub)) {
      d_ac <- d_sub[d_sub$ac_group %in% c("mg_k","k_only"),]
      d_ac$ac_trt <- as.integer(d_ac$ac_group == "mg_k")
      r_ac <- run_ow(d_ac, "ac_trt", "aki_kdigo1", db$nm, NULL)
      if (!is.null(r_ac)) {
        cat(sprintf("    AC: %s\n", fmt_or(r_ac)))
        t5b_rows[[length(t5b_rows)+1]] <- data.frame(
          db=db$nm, subgroup=age_cut$label, analysis="ac_ow", n=nrow(d_ac),
          n_trt=sum(d_ac$ac_trt), or=round(r_ac$or,3), lo=round(r_ac$lo,3),
          hi=round(r_ac$hi,3), p=round(r_ac$p,4))
      }
    }
  }
}
if (length(t5b_rows) > 0) {
  t5b <- do.call(rbind, t5b_rows)
  write.csv(t5b, file.path(RESULTS, "etable5b_age_subgroups.csv"), row.names=FALSE)
  cat("  Saved etable5b_age_subgroups.csv\n")
}

# ============================================================================
# eTABLE 5c: Treatment Effect by Baseline Serum Mg
# ============================================================================
cat("\n", strrep("=",60), "\n")
cat("eTABLE 5c: Effect by Baseline Serum Mg\n")
cat(strrep("=",60), "\n")

t5c_rows <- list()
for (db in list(list(d=dat_e, nm="eICU", cl=TRUE), list(d=dat_m, nm="MIMIC", cl=FALSE))) {
  d <- db$d
  if (!"first_mg_value" %in% names(d)) next
  mg_med <- median(d$first_mg_value, na.rm=TRUE)
  cat(sprintf("  %s: Mg median = %.2f mg/dL\n", db$nm, mg_med))

  for (mg_cut in list(
    list(label=sprintf("Mg <%.1f (below median)", mg_med), lo=-Inf, hi=mg_med),
    list(label=sprintf("Mg >=%.1f (above median)", mg_med), lo=mg_med, hi=Inf),
    list(label="Mg <2.0 (hypoMg)", lo=-Inf, hi=2.0),
    list(label="Mg >=2.0 (normal+)", lo=2.0, hi=Inf)
  )) {
    d_sub <- d[d$first_mg_value >= mg_cut$lo & d$first_mg_value < mg_cut$hi,]
    n_trt <- sum(d_sub$mg_supp, na.rm=TRUE)
    cl <- if (db$cl && "hospitalid" %in% names(d_sub)) d_sub$hospitalid else NULL
    r <- run_ow(d_sub, "mg_supp", "aki_kdigo1", db$nm, cl)
    if (!is.null(r)) {
      cat(sprintf("  %s %s (N=%d, trt=%d): %s\n", db$nm, mg_cut$label, nrow(d_sub), n_trt, fmt_or(r)))
      t5c_rows[[length(t5c_rows)+1]] <- data.frame(
        db=db$nm, subgroup=mg_cut$label, n=nrow(d_sub),
        n_trt=n_trt, or=round(r$or,3), lo=round(r$lo,3), hi=round(r$hi,3), p=round(r$p,4))
    }
  }
}
if (length(t5c_rows) > 0) {
  t5c <- do.call(rbind, t5c_rows)
  write.csv(t5c, file.path(RESULTS, "etable5c_mg_subgroups.csv"), row.names=FALSE)
  cat("  Saved etable5c_mg_subgroups.csv\n")
}

# ============================================================================
# eTABLE 6: Safety Outcomes
# ============================================================================
cat("\n", strrep("=",60), "\n")
cat("eTABLE 6: Safety Outcomes\n")
cat(strrep("=",60), "\n")

t6_rows <- list()

for (db in list(list(d=dat_e, nm="eICU", cl=TRUE), list(d=dat_m, nm="MIMIC", cl=FALSE))) {
  d <- db$d
  cl <- if (db$cl && "hospitalid" %in% names(d)) d$hospitalid else NULL

  # Quick PS + OW for all patients (median impute for speed)
  ps_vars <- intersect(c("age","is_female","bmi","heart_failure","hypertension","diabetes","ckd",
    "copd","pvd","stroke","liver_disease","baseline_creatinine","egfr",
    "loop_diuretics","nsaids","acei_arb","ppi","beta_blockers","steroids","antiarrhythmics",
    "first_potassium","first_calcium","first_heartrate","vasopressor_6h","first_mg_value"), names(d))
  if ("surgery_type" %in% names(d)) {
    d$s_cabg <- as.integer(d$surgery_type=="cabg")
    d$s_valve <- as.integer(d$surgery_type=="valve")
    d$s_combined <- as.integer(d$surgery_type=="combined")
    ps_vars <- c(ps_vars, "s_cabg","s_valve","s_combined")
  }
  for (v in ps_vars) if (v %in% names(d) && any(is.na(d[[v]])))
    d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm=TRUE)
  d <- d[complete.cases(d[,ps_vars]),]

  tryCatch({
    fml <- as.formula(paste("mg_supp ~", paste(ps_vars, collapse="+")))
    ps_fit <- glm(fml, data=d, family=binomial())
    d$ps <- pmax(pmin(fitted(ps_fit), 0.99), 0.01)
    d$ow <- ifelse(d$mg_supp==1, 1-d$ps, d$ps)
  }, error=function(e) { cat(sprintf("  %s PS failed: %s\n", db$nm, e$message)); return() })

  cat(sprintf("\n  %s (N=%d):\n", db$nm, nrow(d)))

  # ── ICU LOS (continuous, OW-weighted) ───────────────────────────
  if ("icu_los_h" %in% names(d)) {
    d_los <- d[!is.na(d$icu_los_h),]
    r <- wlm(icu_los_h ~ mg_supp, d_los, d_los$ow, cl)
    cat(sprintf("    ICU LOS: %s\n", fmt_d(r)))
    # Descriptive
    trt <- d_los$icu_los_h[d_los$mg_supp==1]; ctrl <- d_los$icu_los_h[d_los$mg_supp==0]
    cat(sprintf("      Trt: median %.1fh (IQR %.1f-%.1f), Ctrl: %.1fh (%.1f-%.1f)\n",
                median(trt), quantile(trt,.25), quantile(trt,.75),
                median(ctrl), quantile(ctrl,.25), quantile(ctrl,.75)))
    t6_rows[[length(t6_rows)+1]] <- data.frame(
      db=db$nm, outcome="ICU LOS (h)", est=round(r$est,2),
      lo=round(r$lo,2), hi=round(r$hi,2), p=round(r$p,4), type="continuous")
  }

  # ── Vent duration ───────────────────────────────────────────────
  if ("vent_duration_h" %in% names(d)) {
    d_v <- d[!is.na(d$vent_duration_h),]
    if (nrow(d_v) > 50) {
      r <- wlm(vent_duration_h ~ mg_supp, d_v, d_v$ow, cl)
      cat(sprintf("    Vent duration: %s\n", fmt_d(r)))
      t6_rows[[length(t6_rows)+1]] <- data.frame(
        db=db$nm, outcome="Vent duration (h)", est=round(r$est,2),
        lo=round(r$lo,2), hi=round(r$hi,2), p=round(r$p,4), type="continuous")
    }
  }

  # ── Post-treatment HR ───────────────────────────────────────────
  if ("post_hr_6h" %in% names(d)) {
    d_hr <- d[!is.na(d$post_hr_6h),]
    if (nrow(d_hr) > 50) {
      r <- wlm(post_hr_6h ~ mg_supp, d_hr, d_hr$ow, cl)
      cat(sprintf("    Post-6h HR: %s\n", fmt_d(r)))
      t6_rows[[length(t6_rows)+1]] <- data.frame(
        db=db$nm, outcome="Post-6h HR (bpm)", est=round(r$est,2),
        lo=round(r$lo,2), hi=round(r$hi,2), p=round(r$p,4), type="continuous")
    }
  }

  # ── Max post-treatment Mg (descriptive + safety flag) ───────────
  if ("max_posttreat_mg" %in% names(d)) {
    d_mg <- d[!is.na(d$max_posttreat_mg),]
    trt_mg <- d_mg$max_posttreat_mg[d_mg$mg_supp==1]
    ctrl_mg <- d_mg$max_posttreat_mg[d_mg$mg_supp==0]
    cat(sprintf("    Max post-Mg: Trt median %.2f (max %.2f), Ctrl median %.2f (max %.2f)\n",
                median(trt_mg), max(trt_mg), median(ctrl_mg), max(ctrl_mg)))
    cat(sprintf("      >4.8 mg/dL: Trt %d/%d, Ctrl %d/%d\n",
                sum(trt_mg>4.8), length(trt_mg), sum(ctrl_mg>4.8), length(ctrl_mg)))
    cat(sprintf("      >6.0 mg/dL: Trt %d/%d, Ctrl %d/%d\n",
                sum(trt_mg>6.0), length(trt_mg), sum(ctrl_mg>6.0), length(ctrl_mg)))
    t6_rows[[length(t6_rows)+1]] <- data.frame(
      db=db$nm, outcome="Max Mg >4.8 mg/dL (trt)", est=sum(trt_mg>4.8),
      lo=length(trt_mg), hi=round(100*mean(trt_mg>4.8),1), p=NA, type="count")
  }

  # ── IABP / ECMO ────────────────────────────────────────────────
  for (oc in c("has_iabp", "has_ecmo")) {
    if (oc %in% names(d) && sum(d[[oc]], na.rm=TRUE) >= 5) {
      r <- wglm(as.formula(paste(oc, "~ mg_supp")), d, d$ow, cl)
      cat(sprintf("    %s: %s  [events=%d]\n", oc, fmt_or(r), sum(d[[oc]])))
      t6_rows[[length(t6_rows)+1]] <- data.frame(
        db=db$nm, outcome=oc, est=round(r$or,3),
        lo=round(r$lo,3), hi=round(r$hi,3), p=round(r$p,4), type="binary")
    }
  }
}

if (length(t6_rows) > 0) {
  t6 <- do.call(rbind, t6_rows)
  write.csv(t6, file.path(RESULTS, "etable6_safety.csv"), row.names=FALSE)
  cat("\n  Saved etable6_safety.csv\n")
}

# ── Summary ───────────────────────────────────────────────────────────────
cat(sprintf("\n%s\nSUMMARY\n%s\n", strrep("=",60), strrep("=",60)))
cat("  eTable 3: AC baseline characteristics — printed above\n")
cat("  eTable 4: Surgery-type subgroups — etable4_surgery_subgroups.csv\n")
cat("  eTable 5: MDS subgroups — etable5_mds_subgroups.csv\n")
cat("  eTable 5b: Age subgroups — etable5b_age_subgroups.csv\n")
cat("  eTable 5c: Baseline Mg subgroups — etable5c_mg_subgroups.csv\n")
cat("  eTable 6: Safety outcomes — etable6_safety.csv\n")
cat("\n  Key clinical predictions:\n")
cat("  - If Mg <2.0 shows stronger effect → 'repletion benefits the depleted'\n")
cat("  - If age >=60 shows stronger effect → 'age-related Mg depletion' story\n")
cat("  - If MDS >=2 shows stronger effect → composite risk score guides therapy\n")
cat("  - If none differ → effect generalizes (also fine)\n")
cat("\n  Manuscript sentence if all safety outcomes null:\n")
cat("  'Subgroup analyses by surgery type and magnesium depletion risk\n")
cat("   did not identify significant effect modification (eTables 4-5).\n")
cat("   Magnesium supplementation was not associated with prolonged ICU\n")
cat("   stay, prolonged mechanical ventilation, or bradycardia, and\n")
cat("   post-treatment serum magnesium remained below symptomatic\n")
cat("   thresholds in all patients (eTable 6).'\n")
cat("\nsubgroup_analysis.R COMPLETE\n")
