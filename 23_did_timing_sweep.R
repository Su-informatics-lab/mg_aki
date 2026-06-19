#!/usr/bin/env Rscript
# ============================================================================
# 23_did_timing_sweep.R — Cr timing sweep + secondary outcomes
#
# Part 1: Sweep closest_Xh for X = 12, 18, 24, 36, 48, 72
#   Shows how the DiD evolves over time post-IV-Mg
#   Re-uses r=1, caliper=0.2, replace=T PS match
#
# Part 2: Secondary binary outcomes (not DiD — just matched comparison)
#   Mortality, encephalopathy, POAF, arrhythmia
#   Merged from old cohort CSVs (01_analysis_a_cohort.csv, 04_mimic_cohort.csv)
#
# Run:  Rscript 23_did_timing_sweep.R eicu
#       Rscript 23_did_timing_sweep.R mimic
# ============================================================================

suppressPackageStartupMessages({
  library(MatchIt)
  library(sandwich)
  library(lmtest)
})

RESULTS <- path.expand("~/mg_aki/results")
CALIPER <- 0.2

# Timing targets (hours after IV Mg)
TIMING_TARGETS <- c(12, 18, 24, 36, 48, 72)

PS_COVARS <- c(
  "age", "is_female", "bmi",
  "surg_cabg", "surg_valve", "surg_combined",
  "heart_failure", "hypertension", "diabetes", "ckd",
  "copd", "pvd", "stroke", "liver_disease",
  "egfr",
  "loop_diuretics", "nsaids", "acei_arb", "ppi",
  "beta_blockers", "steroids", "antiarrhythmics",
  "first_potassium", "first_calcium", "first_heartrate",
  "first_mg_value", "first_lactate", "lactate_missing"
)

# ── Helpers ──────────────────────────────────────────────────────────────
median_impute <- function(d, vars) {
  for (v in vars)
    if (v %in% names(d) && any(is.na(d[[v]])))
      d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm = TRUE)
  d
}

compute_smds <- function(d, vars) {
  sapply(vars, function(v) {
    if (!v %in% names(d)) return(NA)
    x1 <- d[[v]][d$treated==1]; x0 <- d[[v]][d$treated==0]
    sp <- sqrt((var(x1,na.rm=T)+var(x0,na.rm=T))/2)
    if (is.na(sp)||sp<1e-10) NA else abs(mean(x1,na.rm=T)-mean(x0,na.rm=T))/sp
  })
}

did_robust <- function(df, ps_vars) {
  if (sum(df$treated==1)<10 || sum(df$treated==0)<10) return(NULL)
  smds <- compute_smds(df, ps_vars)
  adjust <- names(smds[!is.na(smds) & smds > 0.05])
  adjust <- intersect(adjust, names(df))
  adjust <- adjust[sapply(adjust, function(v) var(df[[v]],na.rm=T)>1e-10)]

  fit0 <- lm(delta_cr ~ treated, data=df)
  cl0 <- tryCatch(if(length(unique(df$match_pair_id))>1)
    vcovCL(fit0,cluster=df$match_pair_id) else vcovHC(fit0,type="HC1"),
    error=function(e) vcovHC(fit0,type="HC1"))
  ct0 <- coeftest(fit0, vcov.=cl0)

  if (length(adjust)>0) {
    fml <- as.formula(paste("delta_cr ~ treated +", paste(adjust,collapse="+")))
    fit1 <- tryCatch(lm(fml, data=df), error=function(e) NULL)
    if (!is.null(fit1)) {
      cl1 <- tryCatch(if(length(unique(df$match_pair_id))>1)
        vcovCL(fit1,cluster=df$match_pair_id) else vcovHC(fit1,type="HC1"),
        error=function(e) vcovHC(fit1,type="HC1"))
      ct1 <- coeftest(fit1, vcov.=cl1)
    } else ct1 <- ct0
  } else ct1 <- ct0

  list(n_trt=sum(df$treated==1), n_ctl=sum(df$treated==0),
       did_unadj=ct0["treated","Estimate"], p_unadj=ct0["treated","Pr(>|t|)"],
       did_adj=ct1["treated","Estimate"], se_adj=ct1["treated","Std. Error"],
       p_adj=ct1["treated","Pr(>|t|)"],
       ci_lo=ct1["treated","Estimate"]-1.96*ct1["treated","Std. Error"],
       ci_hi=ct1["treated","Estimate"]+1.96*ct1["treated","Std. Error"])
}

# Binary outcome: logistic regression in matched sample, doubly robust
binary_robust <- function(df, outcome_col, ps_vars) {
  if (!outcome_col %in% names(df)) return(NULL)
  df$y <- df[[outcome_col]]
  df <- df[!is.na(df$y), ]
  n1 <- sum(df$treated==1); n0 <- sum(df$treated==0)
  if (n1<20 || n0<20) return(NULL)
  rate1 <- mean(df$y[df$treated==1]); rate0 <- mean(df$y[df$treated==0])

  # Unadjusted
  fit0 <- tryCatch(glm(y ~ treated, data=df, family=quasibinomial()),
                   error=function(e) NULL)
  if (is.null(fit0)) return(NULL)
  cl0 <- tryCatch(if(length(unique(df$match_pair_id))>1)
    vcovCL(fit0,cluster=df$match_pair_id) else vcovHC(fit0,type="HC1"),
    error=function(e) vcovHC(fit0,type="HC1"))
  ct0 <- coeftest(fit0, vcov.=cl0)

  # Adjusted
  smds <- compute_smds(df, ps_vars)
  adjust <- names(smds[!is.na(smds) & smds > 0.05])
  adjust <- intersect(adjust, names(df))
  adjust <- adjust[sapply(adjust, function(v) var(df[[v]],na.rm=T)>1e-10)]

  if (length(adjust)>0) {
    fml <- as.formula(paste("y ~ treated +", paste(adjust,collapse="+")))
    fit1 <- tryCatch(glm(fml, data=df, family=quasibinomial()),
                     error=function(e) NULL)
    if (!is.null(fit1)) {
      cl1 <- tryCatch(if(length(unique(df$match_pair_id))>1)
        vcovCL(fit1,cluster=df$match_pair_id) else vcovHC(fit1,type="HC1"),
        error=function(e) vcovHC(fit1,type="HC1"))
      ct1 <- coeftest(fit1, vcov.=cl1)
    } else ct1 <- ct0
  } else ct1 <- ct0

  or <- exp(ct1["treated","Estimate"])
  or_lo <- exp(ct1["treated","Estimate"]-1.96*ct1["treated","Std. Error"])
  or_hi <- exp(ct1["treated","Estimate"]+1.96*ct1["treated","Std. Error"])

  list(outcome=outcome_col, n_trt=n1, n_ctl=n0,
       rate_trt=round(rate1,4), rate_ctl=round(rate0,4),
       or=round(or,3), or_lo=round(or_lo,3), or_hi=round(or_hi,3),
       p=round(ct1["treated","Pr(>|t|)"],4))
}

# ============================================================================
run_timing <- function(db) {
  tag <- tolower(db)
  SEP <- paste(rep("=",70), collapse="")
  cat(sprintf("\n%s\n%s: Cr Timing Sweep + Secondary Outcomes\n%s\n", SEP, db, SEP))

  # ── Load DiD data ────────────────────────────────────────────────────────
  trt <- read.csv(file.path(RESULTS, sprintf("20_did_treated_%s.csv",tag)), stringsAsFactors=F)
  ctl <- read.csv(file.path(RESULTS, sprintf("20_did_control_%s.csv",tag)), stringsAsFactors=F)
  cr_all <- read.csv(file.path(RESULTS, sprintf("20_did_cr_all_%s.csv",tag)), stringsAsFactors=F)

  id_col <- if ("patientunitstayid" %in% names(trt)) "patientunitstayid" else "stay_id"
  trt$pid <- trt[[id_col]]; ctl$pid <- ctl[[id_col]]
  cr_id <- if ("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
  cr_all$pid <- cr_all[[cr_id]]
  if (!"labresultoffset" %in% names(cr_all) && "offset_min" %in% names(cr_all))
    cr_all$labresultoffset <- cr_all$offset_min

  ps_vars <- intersect(PS_COVARS, intersect(names(trt), names(ctl)))
  stack_cols <- intersect(unique(c("pid","treated",ps_vars)),
                          intersect(names(trt), names(ctl)))
  combined <- rbind(trt[,stack_cols], ctl[,stack_cols])
  rownames(combined) <- seq_len(nrow(combined))
  combined <- median_impute(combined, ps_vars)
  cat(sprintf("  %d treated + %d control\n", sum(combined$treated==1), sum(combined$treated==0)))

  # Pre-index Cr
  cr_list <- split(cr_all[,c("pid","labresult","labresultoffset")], cr_all$pid)

  # Cr_pre for treated
  cr_pre_map <- setNames(trt$cr_pre, trt$pid)
  cr_pre_off_map <- setNames(trt$cr_pre_offset_min, trt$pid)
  mg_times <- setNames(trt$mg_offset_min, trt$pid)

  # All post-IV-Mg Cr for treated
  cr_trt <- cr_all[cr_all$pid %in% trt$pid, ]
  cr_trt$mg_off <- mg_times[as.character(cr_trt$pid)]
  cr_trt$post_h <- (cr_trt$labresultoffset - cr_trt$mg_off) / 60

  # ── PS matching (r=1, cached) ────────────────────────────────────────────
  cat("\n  PS matching (r=1, caliper=0.2, replace=T)...\n")
  ps_formula <- as.formula(paste("treated ~", paste(ps_vars, collapse="+")))
  m <- suppressWarnings(matchit(ps_formula, data=combined, method="nearest",
                                 distance="glm", ratio=1, caliper=CALIPER, replace=TRUE))
  md <- match.data(m)
  cat(sprintf("  Matched: %d trt, %d ctl\n", sum(md$treated==1), sum(md$treated==0)))

  # Extract pairs
  mm <- m$match.matrix
  trt_idx <- as.integer(rownames(mm))
  pairs <- list()
  for (i in seq_len(nrow(mm)))
    for (j in seq_len(ncol(mm))) {
      ci <- mm[i,j]
      if (!is.na(ci))
        pairs[[length(pairs)+1]] <- c(combined$pid[trt_idx[i]], combined$pid[as.integer(ci)])
    }
  pairs_df <- as.data.frame(do.call(rbind, pairs))
  names(pairs_df) <- c("trt_pid","ctl_pid")

  # Covariate lookup for merge
  covar_want <- c("pid", ps_vars, "surgery_type", "hosp_mortality")
  covar_trt <- trt[, intersect(covar_want, names(trt)), drop=F]
  covar_ctl <- ctl[, intersect(covar_want, names(ctl)), drop=F]
  shared <- intersect(names(covar_trt), names(covar_ctl))
  covar_all <- rbind(covar_trt[,shared], covar_ctl[,shared])
  covar_all <- covar_all[!duplicated(covar_all$pid),]

  # ════════════════════════════════════════════════════════════════════════
  # PART 1: Timing sweep
  # ════════════════════════════════════════════════════════════════════════
  cat(sprintf("\n%s\nPART 1: Cr Timing Sweep (closest_Xh, X = %s)\n%s\n",
              paste(rep("─",60),collapse=""),
              paste(TIMING_TARGETS,collapse=","),
              paste(rep("─",60),collapse="")))

  tol_min <- 6 * 60  # ±6h tolerance throughout
  timing_results <- list()

  for (target_h in TIMING_TARGETS) {
    cat(sprintf("\n── closest_%dh ──\n", target_h))

    # Build Cr_post for this target
    # Window: 0 to 2×target (generous) or at least target±12h
    win_lo <- max(0, target_h - 12)
    win_hi <- target_h + 12
    cr_win <- cr_trt[cr_trt$post_h >= win_lo & cr_trt$post_h <= win_hi, ]
    cr_win$dist <- abs(cr_win$post_h - target_h)
    strat <- cr_win[order(cr_win$pid, cr_win$dist), ]
    strat <- strat[!duplicated(strat$pid), ]

    strat$cr_pre <- cr_pre_map[as.character(strat$pid)]
    strat$cr_pre_offset_min <- cr_pre_off_map[as.character(strat$pid)]
    strat$cr_post <- strat$labresult
    strat$cr_post_offset_min <- strat$labresultoffset
    strat$delta_cr <- strat$cr_post - strat$cr_pre

    n_with_post <- nrow(strat)
    cat(sprintf("  Treated with Cr_post: %d\n", n_with_post))
    if (n_with_post > 0) {
      cat(sprintf("  Actual post_h: median=%.1f, IQR=[%.1f–%.1f]\n",
                  median(strat$post_h), quantile(strat$post_h,0.25),
                  quantile(strat$post_h,0.75)))
    }

    # Temporal matching
    trt_valid <- strat[, c("pid","cr_pre","cr_pre_offset_min","cr_post",
                            "cr_post_offset_min","delta_cr")]
    trt_valid$treated <- 1
    pairs_valid <- pairs_df[pairs_df$trt_pid %in% trt_valid$pid, ]

    n_valid <- 0; ctl_rows <- list()
    for (r in seq_len(nrow(pairs_valid))) {
      tpid <- pairs_valid$trt_pid[r]; cpid <- pairs_valid$ctl_pid[r]
      tidx <- which(trt_valid$pid==tpid)[1]
      if (is.na(tidx)) next
      t_pre <- trt_valid$cr_pre_offset_min[tidx]
      t_post <- trt_valid$cr_post_offset_min[tidx]
      ctl_cr <- cr_list[[as.character(cpid)]]
      if (is.null(ctl_cr)||nrow(ctl_cr)<2) next

      pd <- abs(ctl_cr$labresultoffset - t_pre)
      pi <- which.min(pd)
      if (pd[pi] > tol_min) next

      pc <- ctl_cr[ctl_cr$labresultoffset > ctl_cr$labresultoffset[pi], ]
      if (nrow(pc)==0) next
      pod <- abs(pc$labresultoffset - t_post)
      poi <- which.min(pod)
      if (pod[poi] > tol_min) next

      n_valid <- n_valid+1
      ctl_rows[[n_valid]] <- data.frame(
        pid=cpid, cr_pre=ctl_cr$labresult[pi],
        cr_pre_offset_min=ctl_cr$labresultoffset[pi],
        cr_post=pc$labresult[poi],
        cr_post_offset_min=pc$labresultoffset[poi],
        delta_cr=pc$labresult[poi]-ctl_cr$labresult[pi],
        treated=0, match_pair_id=tidx, stringsAsFactors=F)
    }

    if (n_valid==0) {
      cat("  No temporal matches\n")
      timing_results[[as.character(target_h)]] <- data.frame(
        target_h=target_h, n_trt=0, n_ctl=0,
        did_adj=NA, p_adj=NA, ci_lo=NA, ci_hi=NA, stringsAsFactors=F)
      next
    }

    ctl_matched <- do.call(rbind, ctl_rows)
    valid_pairs <- sort(unique(ctl_matched$match_pair_id))
    trt_out <- trt_valid[valid_pairs, ]
    trt_out$match_pair_id <- seq_along(valid_pairs)
    ctl_matched$match_pair_id <- match(ctl_matched$match_pair_id, valid_pairs)
    matched_df <- rbind(trt_out, ctl_matched)

    # Merge covariates
    matched_df <- merge(matched_df, covar_all, by="pid", all.x=T, suffixes=c("",".c"))
    matched_df <- median_impute(matched_df, ps_vars)

    cat(sprintf("  Matched: %d trt, %d ctl\n", nrow(trt_out), n_valid))

    res <- did_robust(matched_df, ps_vars)
    if (!is.null(res)) {
      sig <- if (res$p_adj < 0.05) " *" else ""
      cat(sprintf("  DiD = %+.4f (P=%.4f, CI [%+.4f,%+.4f])%s\n",
                  res$did_adj, res$p_adj, res$ci_lo, res$ci_hi, sig))
      timing_results[[as.character(target_h)]] <- data.frame(
        target_h=target_h, n_trt=res$n_trt, n_ctl=res$n_ctl,
        did_unadj=round(res$did_unadj,4), did_adj=round(res$did_adj,4),
        se_adj=round(res$se_adj,4), p_adj=round(res$p_adj,4),
        ci_lo=round(res$ci_lo,4), ci_hi=round(res$ci_hi,4),
        stringsAsFactors=F)
    }
  }

  # Timing summary table
  timing_df <- do.call(rbind, timing_results)
  cat(sprintf("\n%s\nTIMING SWEEP SUMMARY\n%s\n",
              paste(rep("─",60),collapse=""), paste(rep("─",60),collapse="")))
  cat("\n  target_h  n_trt  n_ctl  DiD_unadj  DiD_adj    P       95% CI\n")
  cat("  ────────  ─────  ─────  ─────────  ────────  ──────  ──────────────\n")
  for (i in seq_len(nrow(timing_df))) {
    r <- timing_df[i,]
    sig <- if (!is.na(r$p_adj) && r$p_adj < 0.05) " *" else ""
    cat(sprintf("    %3dh    %5d  %5d   %+.4f   %+.4f  %.4f  [%+.4f,%+.4f]%s\n",
                r$target_h, r$n_trt, r$n_ctl,
                ifelse(is.na(r$did_unadj),0,r$did_unadj),
                ifelse(is.na(r$did_adj),0,r$did_adj),
                ifelse(is.na(r$p_adj),1,r$p_adj),
                ifelse(is.na(r$ci_lo),0,r$ci_lo),
                ifelse(is.na(r$ci_hi),0,r$ci_hi), sig))
  }

  write.csv(timing_df, file.path(RESULTS, sprintf("23_timing_sweep_%s.csv",tag)), row.names=F)

  # ════════════════════════════════════════════════════════════════════════
  # PART 2: Secondary binary outcomes
  # ════════════════════════════════════════════════════════════════════════
  cat(sprintf("\n%s\nPART 2: Secondary Outcomes (matched comparison)\n%s\n",
              paste(rep("─",60),collapse=""), paste(rep("─",60),collapse="")))

  # Load the primary matched dataset (closest_24h)
  primary_path <- file.path(RESULTS, sprintf("21_matched_%s_sens_at24.csv",tag))
  if (!file.exists(primary_path)) {
    primary_path <- file.path(RESULTS, sprintf("21_matched_%s_primary.csv",tag))
  }
  if (!file.exists(primary_path)) {
    cat("  No matched dataset found, skipping\n")
    return(timing_df)
  }
  matched <- read.csv(primary_path, stringsAsFactors=F)
  matched <- median_impute(matched, ps_vars)
  cat(sprintf("  Using matched set: %s (%d rows)\n", basename(primary_path), nrow(matched)))

  # Merge secondary outcomes from old cohort CSVs
  cat("  Merging secondary outcomes from old cohort...\n")
  old_eicu <- file.path(RESULTS, "01_analysis_a_cohort.csv")
  old_mimic <- file.path(RESULTS, "04_mimic_cohort.csv")

  if (tag == "eicu" && file.exists(old_eicu)) {
    old <- read.csv(old_eicu, stringsAsFactors=F)
    old$pid <- old$patientunitstayid
    # Available outcomes
    outcome_cols <- intersect(c("hosp_mortality", "icu_mortality",
                                 "neuro_encephalopathy", "nc_poaf",
                                 "vent_arrhythmia", "aki_kdigo1",
                                 "nc_transfusion"),
                              names(old))
    cat(sprintf("  Old cohort outcomes available: %s\n", paste(outcome_cols, collapse=", ")))
    old_sub <- old[, c("pid", outcome_cols)]
    matched <- merge(matched, old_sub, by="pid", all.x=T, suffixes=c("",".old"))
  } else if (tag == "mimic" && file.exists(old_mimic)) {
    old <- read.csv(old_mimic, stringsAsFactors=F)
    old$pid <- old$stay_id
    outcome_cols <- intersect(c("hosp_mortality", "icu_mortality",
                                 "neuro_encephalopathy", "nc_poaf",
                                 "vent_arrhythmia", "aki_kdigo1",
                                 "nc_transfusion"),
                              names(old))
    # Rename if needed
    rename_map <- c(mg_supplementation="mg_supp")
    for (nm in names(rename_map))
      if (nm %in% names(old) && !rename_map[nm] %in% names(old))
        names(old)[names(old)==nm] <- rename_map[nm]
    outcome_cols <- intersect(c("hosp_mortality",
                                 "neuro_encephalopathy", "nc_poaf",
                                 "vent_arrhythmia", "aki_kdigo1"),
                              names(old))
    cat(sprintf("  Old cohort outcomes available: %s\n", paste(outcome_cols, collapse=", ")))
    old_sub <- old[, c("pid", outcome_cols)]
    matched <- merge(matched, old_sub, by="pid", all.x=T, suffixes=c("",".old"))
  } else {
    cat("  Old cohort CSV not found\n")
    outcome_cols <- character(0)
  }

  # Also use hosp_mortality already in matched set
  if (!"hosp_mortality" %in% outcome_cols && "hosp_mortality" %in% names(matched))
    outcome_cols <- c("hosp_mortality", outcome_cols)

  # Run binary outcome analyses
  outcome_labels <- c(
    hosp_mortality = "Hospital mortality",
    icu_mortality = "ICU mortality",
    neuro_encephalopathy = "Encephalopathy",
    nc_poaf = "POAF",
    vent_arrhythmia = "Ventricular arrhythmia",
    aki_kdigo1 = "AKI (KDIGO ≥1, reference)"
  )

  cat("\n  outcome                n_trt  n_ctl  rate_trt  rate_ctl  OR (95% CI)        P\n")
  cat("  ───────────────────── ────── ────── ──────── ──────── ─────────────────── ──────\n")

  sec_results <- list()
  for (oc in outcome_cols) {
    res <- binary_robust(matched, oc, ps_vars)
    if (is.null(res)) next
    label <- ifelse(oc %in% names(outcome_labels), outcome_labels[oc], oc)
    sig <- if (res$p < 0.05) " *" else ""
    cat(sprintf("  %-23s %5d  %5d    %.1f%%    %.1f%%  %.2f (%.2f–%.2f)  %.4f%s\n",
                label, res$n_trt, res$n_ctl,
                100*res$rate_trt, 100*res$rate_ctl,
                res$or, res$or_lo, res$or_hi, res$p, sig))
    sec_results[[oc]] <- as.data.frame(res, stringsAsFactors=F)
  }

  if (length(sec_results) > 0) {
    sec_df <- do.call(rbind, sec_results)
    write.csv(sec_df, file.path(RESULTS, sprintf("23_secondary_%s.csv",tag)), row.names=F)
    cat(sprintf("\n  Saved: 23_secondary_%s.csv\n", tag))
  }

  return(timing_df)
}

# ============================================================================
cat("======================================================================\n")
cat("23_did_timing_sweep.R — Cr timing + secondary outcomes\n")
cat(sprintf("  Timing targets: closest_%sh\n", paste(TIMING_TARGETS, collapse="h, closest_")))
cat("  Secondary: mortality, encephalopathy, POAF, arrhythmia\n")
cat("======================================================================\n")

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) { cat("Usage: Rscript 23_did_timing_sweep.R eicu|mimic\n"); quit(status=1) }
for (a in args) run_timing(toupper(a))

cat("\n======================================================================\n")
cat("Key: where does the DiD peak? If at 24-36h → consistent with Cr kinetics\n")
cat("     If monotonically increasing → treatment effect accumulates over time\n")
cat("======================================================================\n")
