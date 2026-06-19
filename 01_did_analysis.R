#!/usr/bin/env Rscript
# ============================================================================
# 01_did_analysis.R — DiD analysis with model selection
#
# Usage:  Rscript 01_did_analysis.R <db> [model]
#   db:    eicu | mimic
#   model: v3 (default) | yan | original
#
# Output files include model tag: did_timecourse_eicu_v3.csv etc.
# ============================================================================

suppressPackageStartupMessages({
  library(MatchIt)
  library(sandwich)
  library(lmtest)
})

source(file.path(path.expand("~/mg_aki"), "did_covars.R"))

RESULTS   <- path.expand("~/mg_aki/results")
CALIPER   <- 0.2
TARGETS   <- seq(6, 36, by = 3)
TOLERANCES <- c(2, 4, 6)
PRIMARY_TARGET <- 24
PRIMARY_TOL    <- 6

MG_BINS <- list("<1.8"=c(0,1.8), "1.8-2.0"=c(1.8,2.0),
                "2.0-2.3"=c(2.0,2.3), ">2.3"=c(2.3,99))
SURG_TYPES <- c("cabg","valve","combined","other_cardiac")
SECONDARY_OUTCOMES <- c("hosp_mortality","poaf","encephalopathy","vent_arrhythmia")
OUTCOME_LABELS <- c(hosp_mortality="Hospital mortality", poaf="POAF (new-onset)",
                     encephalopathy="Encephalopathy", vent_arrhythmia="Ventricular arrhythmia")

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

did_robust <- function(df, ps_vars, pair_col="match_pair_id") {
  nt <- sum(df$treated==1); nc <- sum(df$treated==0)
  if (nt<10||nc<10) return(NULL)
  smds <- compute_smds(df, ps_vars)
  adj <- names(smds[!is.na(smds) & smds > 0.05])
  adj <- intersect(adj, names(df))
  adj <- adj[vapply(adj, function(v) var(df[[v]],na.rm=T)>1e-10, logical(1))]

  get_vcov <- function(fit) {
    if (pair_col %in% names(df) && length(unique(df[[pair_col]]))>1)
      tryCatch(vcovCL(fit,cluster=df[[pair_col]]),
               error=function(e) vcovHC(fit,type="HC1"))
    else vcovHC(fit, type="HC1")
  }

  fit0 <- lm(delta_cr ~ treated, data=df)
  ct0 <- coeftest(fit0, vcov.=get_vcov(fit0))
  if (length(adj)>0) {
    fml <- as.formula(paste("delta_cr ~ treated +", paste(adj,collapse="+")))
    fit1 <- tryCatch(lm(fml,data=df), error=function(e) NULL)
    if (!is.null(fit1)) ct1 <- coeftest(fit1, vcov.=get_vcov(fit1))
    else ct1 <- ct0
  } else ct1 <- ct0

  list(n_trt=nt, n_ctl=nc,
       did_unadj=ct0["treated","Estimate"], p_unadj=ct0["treated","Pr(>|t|)"],
       did_adj=ct1["treated","Estimate"], se_adj=ct1["treated","Std. Error"],
       p_adj=ct1["treated","Pr(>|t|)"],
       ci_lo=ct1["treated","Estimate"]-1.96*ct1["treated","Std. Error"],
       ci_hi=ct1["treated","Estimate"]+1.96*ct1["treated","Std. Error"],
       n_adjust=length(adj))
}

binary_outcome <- function(df, outcome_col, ps_vars, pair_col="match_pair_id") {
  if (!outcome_col %in% names(df)) return(NULL)
  df$y <- as.numeric(df[[outcome_col]])
  df <- df[!is.na(df$y), ]
  nt <- sum(df$treated==1); nc <- sum(df$treated==0)
  if (nt<20||nc<20) return(NULL)
  rate1 <- mean(df$y[df$treated==1]); rate0 <- mean(df$y[df$treated==0])
  if (rate1==0 && rate0==0) return(NULL)

  fit0 <- tryCatch(glm(y ~ treated, data=df, family=quasibinomial()), error=function(e) NULL)
  if (is.null(fit0) || !"treated" %in% names(coef(fit0))) return(NULL)

  smds <- compute_smds(df, ps_vars)
  adj <- names(smds[!is.na(smds) & smds > 0.05])
  adj <- intersect(adj, names(df))
  adj <- adj[vapply(adj, function(v) var(df[[v]],na.rm=T)>1e-10, logical(1))]

  if (length(adj)>0) {
    fml <- as.formula(paste("y ~ treated +", paste(adj,collapse="+")))
    fit <- tryCatch(glm(fml, data=df, family=quasibinomial()), error=function(e) NULL)
    if (is.null(fit) || !"treated" %in% names(coef(fit))) fit <- fit0
  } else fit <- fit0

  vc <- tryCatch({
    if (pair_col %in% names(df) && length(unique(df[[pair_col]]))>1)
      vcovCL(fit, cluster=df[[pair_col]])
    else vcovHC(fit, type="HC1")
  }, error=function(e) tryCatch(vcovHC(fit, type="HC1"), error=function(e2) vcov(fit)))

  ct <- tryCatch(coeftest(fit, vcov.=vc), error=function(e) NULL)
  if (is.null(ct)) return(NULL)
  trt_row <- which(rownames(ct) == "treated")
  if (length(trt_row)==0) return(NULL)
  p_col <- grep("^Pr", colnames(ct))
  if (length(p_col)==0) return(NULL)

  or <- exp(ct[trt_row, "Estimate"])
  list(outcome=outcome_col, n_trt=nt, n_ctl=nc,
       rate_trt=round(rate1,4), rate_ctl=round(rate0,4), or=round(or,3),
       or_lo=round(exp(ct[trt_row,"Estimate"]-1.96*ct[trt_row,"Std. Error"]),3),
       or_hi=round(exp(ct[trt_row,"Estimate"]+1.96*ct[trt_row,"Std. Error"]),3),
       p=round(ct[trt_row, p_col],4))
}

temporal_align <- function(trt_valid, pairs_df, cr_list, tol_min) {
  n_valid <- 0; ctl_rows <- list()
  for (r in seq_len(nrow(pairs_df))) {
    tpid <- pairs_df$trt_pid[r]; cpid <- pairs_df$ctl_pid[r]
    tidx <- which(trt_valid$pid==tpid)[1]
    if (is.na(tidx)) next
    t_pre <- trt_valid$cr_pre_offset_min[tidx]
    t_post <- trt_valid$cr_post_offset_min[tidx]
    cc <- cr_list[[as.character(cpid)]]
    if (is.null(cc)||nrow(cc)<2) next
    pd <- abs(cc$labresultoffset-t_pre); pi <- which.min(pd)
    if (pd[pi]>tol_min) next
    pc <- cc[cc$labresultoffset > cc$labresultoffset[pi],]
    if (nrow(pc)==0) next
    pod <- abs(pc$labresultoffset-t_post); poi <- which.min(pod)
    if (pod[poi]>tol_min) next
    n_valid <- n_valid+1
    ctl_rows[[n_valid]] <- data.frame(
      pid=cpid, cr_pre=cc$labresult[pi],
      cr_pre_offset_min=cc$labresultoffset[pi],
      cr_post=pc$labresult[poi],
      cr_post_offset_min=pc$labresultoffset[poi],
      delta_cr=pc$labresult[poi]-cc$labresult[pi],
      treated=0, match_pair_id=tidx, stringsAsFactors=F)
  }
  if (n_valid==0) return(NULL)
  do.call(rbind, ctl_rows)
}

build_matched <- function(trt_valid, ctl_matched, covar_lu, ps_vars) {
  vp <- sort(unique(ctl_matched$match_pair_id))
  t_out <- trt_valid[vp,]; t_out$match_pair_id <- seq_along(vp)
  ctl_matched$match_pair_id <- match(ctl_matched$match_pair_id, vp)
  mdf <- rbind(t_out, ctl_matched)
  mdf <- merge(mdf, covar_lu, by="pid", all.x=T, suffixes=c("",".c"))
  median_impute(mdf, ps_vars)
}


# ============================================================================
# MAIN
# ============================================================================
run_analysis <- function(db, model_name, ps_vars) {
  tag <- tolower(db)
  mtag <- tolower(model_name)
  ftag <- sprintf("%s_%s", tag, mtag)  # for output filenames
  SEP <- paste(rep("=",70),collapse="")
  cat(sprintf("\n%s\n%s [%s]: %d PS covariates, +/-6h tolerance\n%s\n",
              SEP, db, mtag, length(ps_vars), SEP))

  trt <- read.csv(file.path(RESULTS,sprintf("did_treated_%s.csv",tag)),stringsAsFactors=F)
  ctl <- read.csv(file.path(RESULTS,sprintf("did_control_%s.csv",tag)),stringsAsFactors=F)
  cr_all <- read.csv(file.path(RESULTS,sprintf("did_cr_all_%s.csv",tag)),stringsAsFactors=F)

  id_col <- if("patientunitstayid" %in% names(trt)) "patientunitstayid" else "stay_id"
  trt$pid <- trt[[id_col]]; ctl$pid <- ctl[[id_col]]
  cr_id <- if("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
  cr_all$pid <- cr_all[[cr_id]]
  if(!"labresultoffset" %in% names(cr_all)) cr_all$labresultoffset <- cr_all$offset_min

  ps_avail <- intersect(ps_vars, intersect(names(trt),names(ctl)))
  cat(sprintf("  PS covariates: %d requested, %d available\n", length(ps_vars), length(ps_avail)))
  missing_v <- setdiff(ps_vars, names(trt))
  if (length(missing_v)>0) cat(sprintf("  WARNING missing: %s\n", paste(missing_v,collapse=", ")))

  stack_cols <- intersect(unique(c("pid","treated",ps_avail)), intersect(names(trt),names(ctl)))
  combined <- rbind(trt[,stack_cols], ctl[,stack_cols])
  rownames(combined) <- seq_len(nrow(combined))
  combined <- median_impute(combined, ps_avail)
  cat(sprintf("  Loaded: %d treated + %d control\n",
              sum(combined$treated==1), sum(combined$treated==0)))

  # ── A. PS matching ──────────────────────────────────────────────────────
  cat("  PS matching (r=1, caliper=0.2, replace=T)...\n")
  ps_fml <- as.formula(paste("treated ~", paste(ps_avail,collapse="+")))
  m <- suppressWarnings(matchit(ps_fml, data=combined, method="nearest",
                                 distance="glm", ratio=1, caliper=CALIPER, replace=TRUE))
  md <- match.data(m)
  smds_raw <- compute_smds(combined, ps_avail)
  smds_matched <- compute_smds(md, ps_avail)
  cat(sprintf("  Raw max SMD: %.3f | Matched max SMD: %.3f | n>0.1: %d/%d\n",
              max(smds_raw,na.rm=T), max(smds_matched,na.rm=T),
              sum(smds_matched>0.1,na.rm=T), length(ps_avail)))
  worst <- sort(smds_matched[smds_matched > 0.1], decreasing=TRUE)
  if (length(worst)>0) {
    cat("  Covariates with matched SMD > 0.1:\n")
    for (nm in names(worst)) cat(sprintf("    %-25s %.3f\n", nm, worst[nm]))
  }

  mm <- m$match.matrix; trt_idx <- as.integer(rownames(mm))
  pairs <- list()
  for(i in seq_len(nrow(mm))) for(j in seq_len(ncol(mm))) {
    ci <- mm[i,j]; if(!is.na(ci))
      pairs[[length(pairs)+1]] <- c(combined$pid[trt_idx[i]], combined$pid[as.integer(ci)])
  }
  pairs_df <- as.data.frame(do.call(rbind,pairs)); names(pairs_df) <- c("trt_pid","ctl_pid")
  cat(sprintf("  Pairs: %d\n", nrow(pairs_df)))

  cr_list <- split(cr_all[,c("pid","labresult","labresultoffset")], cr_all$pid)
  mg_times <- setNames(trt$mg_offset_min, trt$pid)
  cr_pre_map <- setNames(trt$cr_pre, trt$pid)
  cr_pre_off_map <- setNames(trt$cr_pre_offset_min, trt$pid)

  cr_trt <- cr_all[cr_all$pid %in% trt$pid,]
  cr_trt$mg_off <- mg_times[as.character(cr_trt$pid)]
  cr_trt$post_h <- (cr_trt$labresultoffset - cr_trt$mg_off) / 60

  covar_want <- unique(c("pid", ps_avail, "surgery_type", "hosp_mortality",
                          "poaf", "encephalopathy", "vent_arrhythmia", "prior_af"))
  cov_t <- trt[,intersect(covar_want,names(trt)),drop=F]
  cov_c <- ctl[,intersect(covar_want,names(ctl)),drop=F]
  sh <- intersect(names(cov_t),names(cov_c))
  covar_lu <- rbind(cov_t[,sh], cov_c[,sh])
  covar_lu <- covar_lu[!duplicated(covar_lu$pid),]

  # ── B. Time course ──────────────────────────────────────────────────────
  cat(sprintf("\n  TIME COURSE (6-36h, +/-2/4/6h)\n"))
  cat("  target  tol   n_trt  n_ctl    DiD      P       95% CI\n")
  cat("  ------  ---  -----  -----  --------  ------  --------------\n")

  tc_results <- list(); tidx <- 0

  for (target_h in TARGETS) {
    win_lo <- max(0, target_h-12); win_hi <- target_h+12
    cw <- cr_trt[cr_trt$post_h>=win_lo & cr_trt$post_h<=win_hi,]
    cw$dist <- abs(cw$post_h - target_h)
    strat <- cw[order(cw$pid, cw$dist),]; strat <- strat[!duplicated(strat$pid),]
    strat$cr_pre <- cr_pre_map[as.character(strat$pid)]
    strat$cr_pre_offset_min <- cr_pre_off_map[as.character(strat$pid)]
    strat$cr_post <- strat$labresult
    strat$cr_post_offset_min <- strat$labresultoffset
    strat$delta_cr <- strat$cr_post - strat$cr_pre
    tv <- strat[,c("pid","cr_pre","cr_pre_offset_min","cr_post","cr_post_offset_min","delta_cr")]
    tv$treated <- 1
    pv <- pairs_df[pairs_df$trt_pid %in% tv$pid,]

    for (tol_h in TOLERANCES) {
      ctl_m <- temporal_align(tv, pv, cr_list, tol_h*60)
      if (is.null(ctl_m) || nrow(ctl_m)<10) {
        tidx <- tidx+1
        tc_results[[tidx]] <- data.frame(target_h=target_h, tol_h=tol_h,
          n_trt=0, n_ctl=0, did_adj=NA, se_adj=NA, p_adj=NA,
          ci_lo=NA, ci_hi=NA, stringsAsFactors=F)
        next
      }
      mdf <- build_matched(tv, ctl_m, covar_lu, ps_avail)
      res <- did_robust(mdf, ps_avail)
      tidx <- tidx+1
      if (!is.null(res)) {
        tc_results[[tidx]] <- data.frame(target_h=target_h, tol_h=tol_h,
          n_trt=res$n_trt, n_ctl=res$n_ctl,
          did_adj=round(res$did_adj,4), se_adj=round(res$se_adj,4),
          p_adj=round(res$p_adj,4),
          ci_lo=round(res$ci_lo,4), ci_hi=round(res$ci_hi,4),
          stringsAsFactors=F)
        primary <- (res$n_trt>0 && target_h==PRIMARY_TARGET && tol_h==PRIMARY_TOL)
        sig <- if(res$p_adj<0.05) " *" else ""
        tag_str <- if(primary) " << PRIMARY" else ""
        cat(sprintf("  %4dh   +/-%dh  %5d  %5d  %+.4f  %.4f  [%+.4f,%+.4f]%s%s\n",
                    target_h,tol_h,res$n_trt,res$n_ctl,res$did_adj,res$p_adj,
                    res$ci_lo,res$ci_hi,sig,tag_str))
      }
    }
  }

  tc <- do.call(rbind, tc_results)
  write.csv(tc, file.path(RESULTS,sprintf("did_timecourse_%s.csv",ftag)), row.names=F)

  # ── C. Primary (24h +/-6h) ─────────────────────────────────────────────
  cat(sprintf("\n  PRIMARY DiD (24h +/-6h)\n"))
  cw24 <- cr_trt[cr_trt$post_h>=12 & cr_trt$post_h<=36,]
  cw24$dist <- abs(cw24$post_h - 24)
  s24 <- cw24[order(cw24$pid,cw24$dist),]; s24 <- s24[!duplicated(s24$pid),]
  s24$cr_pre <- cr_pre_map[as.character(s24$pid)]
  s24$cr_pre_offset_min <- cr_pre_off_map[as.character(s24$pid)]
  s24$cr_post <- s24$labresult; s24$cr_post_offset_min <- s24$labresultoffset
  s24$delta_cr <- s24$cr_post - s24$cr_pre
  tv24 <- s24[,c("pid","cr_pre","cr_pre_offset_min","cr_post","cr_post_offset_min","delta_cr")]
  tv24$treated <- 1
  pv24 <- pairs_df[pairs_df$trt_pid %in% tv24$pid,]
  ctl24 <- temporal_align(tv24, pv24, cr_list, PRIMARY_TOL*60)
  primary_matched <- build_matched(tv24, ctl24, covar_lu, ps_avail)

  res24 <- did_robust(primary_matched, ps_avail)
  cat(sprintf("  Unadjusted: DiD=%+.4f, P=%.4f\n", res24$did_unadj, res24$p_unadj))
  cat(sprintf("  Doubly robust: DiD=%+.4f (SE=%.4f), P=%.4f\n",
              res24$did_adj, res24$se_adj, res24$p_adj))
  cat(sprintf("  95%% CI: [%+.4f, %+.4f]\n", res24$ci_lo, res24$ci_hi))
  cat(sprintf("  n = %d/%d, %d covariates adjusted\n",
              res24$n_trt, res24$n_ctl, res24$n_adjust))

  write.csv(primary_matched,
            file.path(RESULTS,sprintf("did_matched_%s_24h.csv",ftag)), row.names=F)

  # ── D. Subgroups ─────────────────────────────────────────────────────────
  cat(sprintf("\n  SUBGROUPS\n"))
  sg_results <- list(); sidx <- 0
  pm <- primary_matched
  pair_ids <- pm$match_pair_id[pm$treated==1]

  cat("  Mg strata:\n")
  for (mg_nm in names(MG_BINS)) {
    lo <- MG_BINS[[mg_nm]][1]; hi <- MG_BINS[[mg_nm]][2]
    trt_mg <- pm$first_mg_value[pm$treated==1]
    in_st <- !is.na(trt_mg) & trt_mg>=lo & trt_mg<hi
    sub <- pm[pm$match_pair_id %in% pair_ids[in_st],]
    nt <- sum(sub$treated==1); nc <- sum(sub$treated==0)
    res <- if(nt>=10&&nc>=10) did_robust(sub, ps_avail) else NULL
    sig <- if(!is.null(res) && res$p_adj<0.05) " *" else ""
    if(!is.null(res))
      cat(sprintf("    %s: n=%d/%d, DiD=%+.4f, P=%.4f%s\n",
                  mg_nm,nt,nc,res$did_adj,res$p_adj,sig))
    else cat(sprintf("    %s: n=%d/%d -- too few\n", mg_nm,nt,nc))
    sidx <- sidx+1
    sg_results[[sidx]] <- data.frame(
      subgroup="Mg_stratum", stratum=mg_nm, n_trt=nt, n_ctl=nc,
      did_adj=if(!is.null(res)) round(res$did_adj,4) else NA,
      p_adj=if(!is.null(res)) round(res$p_adj,4) else NA, stringsAsFactors=F)
  }

  avail_adj <- intersect(ps_avail[ps_avail!="first_mg_value"], names(pm))
  avail_adj <- avail_adj[vapply(avail_adj, function(v) var(pm[[v]],na.rm=T)>1e-10, logical(1))]
  int_fml <- as.formula(paste("delta_cr ~ treated * first_mg_value +",
                               paste(avail_adj,collapse="+")))
  int_fit <- tryCatch(lm(int_fml, data=pm), error=function(e) NULL)
  if(!is.null(int_fit)) {
    int_cl <- tryCatch(vcovCL(int_fit,cluster=pm$match_pair_id),
                       error=function(e) vcovHC(int_fit,type="HC1"))
    int_ct <- coeftest(int_fit, vcov.=int_cl)
    ir <- grep("treated:first_mg_value", rownames(int_ct))
    if(length(ir)>0) cat(sprintf("  Interaction (trt x Mg): b=%+.4f, P=%.4f\n",
                                  int_ct[ir,"Estimate"], int_ct[ir,"Pr(>|t|)"]))
  }

  cat("\n  Surgery type:\n")
  if ("surgery_type" %in% names(pm)) {
    for (st in SURG_TYPES) {
      trt_st <- pm$surgery_type[pm$treated==1]
      in_st <- !is.na(trt_st) & trt_st==st
      sub <- pm[pm$match_pair_id %in% pair_ids[in_st],]
      nt <- sum(sub$treated==1); nc <- sum(sub$treated==0)
      res <- if(nt>=10&&nc>=10) did_robust(sub, ps_avail) else NULL
      sig <- if(!is.null(res) && res$p_adj<0.05) " *" else ""
      if(!is.null(res))
        cat(sprintf("    %s: n=%d/%d, DiD=%+.4f, P=%.4f%s\n",
                    st,nt,nc,res$did_adj,res$p_adj,sig))
      else cat(sprintf("    %s: n=%d/%d -- too few\n", st,nt,nc))
      sidx <- sidx+1
      sg_results[[sidx]] <- data.frame(
        subgroup="Surgery", stratum=st, n_trt=nt, n_ctl=nc,
        did_adj=if(!is.null(res)) round(res$did_adj,4) else NA,
        p_adj=if(!is.null(res)) round(res$p_adj,4) else NA, stringsAsFactors=F)
    }
  }

  sg_df <- do.call(rbind, sg_results)
  write.csv(sg_df, file.path(RESULTS,sprintf("did_subgroups_%s.csv",ftag)), row.names=F)

  # ── E. Secondary outcomes ────────────────────────────────────────────────
  cat(sprintf("\n  SECONDARY OUTCOMES\n"))
  cat("  outcome                n_trt  n_ctl  rate_trt  rate_ctl  OR (95% CI)       P\n")
  cat("  --------------------- ------ ------ -------- -------- ------------------ ------\n")
  sec_results <- list()
  for (oc in SECONDARY_OUTCOMES) {
    res <- binary_outcome(pm, oc, ps_avail)
    if (is.null(res)) {
      cat(sprintf("  %-23s  -- unavailable\n",
                  ifelse(oc %in% names(OUTCOME_LABELS), OUTCOME_LABELS[oc], oc))); next
    }
    label <- ifelse(oc %in% names(OUTCOME_LABELS), OUTCOME_LABELS[oc], oc)
    sig <- if(res$p<0.05) " *" else ""
    cat(sprintf("  %-23s %5d  %5d    %.1f%%    %.1f%%  %.2f (%.2f-%.2f)  %.4f%s\n",
                label,res$n_trt,res$n_ctl,100*res$rate_trt,100*res$rate_ctl,
                res$or,res$or_lo,res$or_hi,res$p,sig))
    sec_results[[oc]] <- as.data.frame(res, stringsAsFactors=F)
  }
  if (length(sec_results)>0) {
    sec_df <- do.call(rbind, sec_results)
    write.csv(sec_df, file.path(RESULTS,sprintf("did_secondary_%s.csv",ftag)), row.names=F)
  }

  cat(sprintf("\n%s\n%s [%s]: COMPLETE\n", SEP, db, mtag))
  cat(sprintf("  did_timecourse_%s.csv | did_matched_%s_24h.csv\n", ftag, ftag))
  cat(sprintf("  did_subgroups_%s.csv  | did_secondary_%s.csv\n", ftag, ftag))
}

# ============================================================================
# CLI
# ============================================================================
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript 01_did_analysis.R <db> [model]\n")
  cat("  db:    eicu | mimic\n")
  cat("  model: v3 (default) | yan | original\n")
  cat("\nExamples:\n")
  cat("  Rscript 01_did_analysis.R eicu v3\n")
  cat("  Rscript 01_did_analysis.R mimic yan\n")
  quit(status=1)
}

db_arg <- toupper(args[1])
model_arg <- if (length(args) >= 2) tolower(args[2]) else "v3"

cat("======================================================================\n")
cat(sprintf("01_did_analysis.R — model=%s, +/-6h tolerance\n", model_arg))
ps_vars <- select_model(model_arg)
cat("======================================================================\n")

run_analysis(db_arg, model_arg, ps_vars)
