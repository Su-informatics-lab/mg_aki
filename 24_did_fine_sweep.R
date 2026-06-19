#!/usr/bin/env Rscript
# ============================================================================
# 24_did_fine_sweep.R — Fine timing grid + OW parallel analysis
#
# Part A: Matching-based DiD — sweep target 6-36h (every 3h) × tol ±2/4/6h
# Part B: OW-DiD — no matching, common ICU-time anchor, OW weights
#   For each target hour X: ΔCr = Cr(closest to Xh from ICU) - Cr_pre
#   OW-weighted regression: ΔCr ~ treated + [covariates]
#
# Run:  Rscript 24_did_fine_sweep.R eicu
#       Rscript 24_did_fine_sweep.R mimic
# ============================================================================

suppressPackageStartupMessages({
  library(MatchIt)
  library(sandwich)
  library(lmtest)
})

RESULTS <- path.expand("~/mg_aki/results")
CALIPER <- 0.2

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

TARGETS <- seq(6, 36, by = 3)         # 6,9,12,...,36
TOLERANCES <- c(2, 4, 6)              # hours

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

robust_did <- function(df, ps_vars, pair_col="match_pair_id") {
  if (sum(df$treated==1)<10||sum(df$treated==0)<10) return(NULL)
  smds <- compute_smds(df, ps_vars)
  adj <- names(smds[!is.na(smds) & smds > 0.05])
  adj <- intersect(adj, names(df))
  adj <- adj[sapply(adj, function(v) var(df[[v]],na.rm=T)>1e-10)]

  fit0 <- lm(delta_cr ~ treated, data=df)
  if (pair_col %in% names(df) && length(unique(df[[pair_col]]))>1) {
    cl <- tryCatch(vcovCL(fit0,cluster=df[[pair_col]]),
                   error=function(e) vcovHC(fit0,type="HC1"))
  } else cl <- vcovHC(fit0, type="HC1")
  ct0 <- coeftest(fit0, vcov.=cl)

  if (length(adj)>0) {
    fml <- as.formula(paste("delta_cr ~ treated +", paste(adj,collapse="+")))
    fit1 <- tryCatch(lm(fml,data=df), error=function(e) NULL)
    if (!is.null(fit1)) {
      cl1 <- if (pair_col %in% names(df) && length(unique(df[[pair_col]]))>1)
        tryCatch(vcovCL(fit1,cluster=df[[pair_col]]),
                 error=function(e) vcovHC(fit1,type="HC1"))
      else vcovHC(fit1, type="HC1")
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

# OW-weighted regression (no pairs)
ow_did <- function(df, ps_vars) {
  if (sum(df$treated==1)<10||sum(df$treated==0)<10) return(NULL)
  # PS model
  avail <- intersect(ps_vars, names(df))
  avail <- avail[sapply(avail, function(v) var(df[[v]],na.rm=T)>1e-10)]
  fml_ps <- as.formula(paste("treated ~", paste(avail, collapse="+")))
  ps_fit <- tryCatch(glm(fml_ps, data=df, family=binomial()),
                     error=function(e) NULL)
  if (is.null(ps_fit)) return(NULL)
  df$ps <- predict(ps_fit, type="response")

  # OW weights: treated get (1-PS), controls get PS
  df$ow <- ifelse(df$treated==1, 1-df$ps, df$ps)
  df$ow <- df$ow / mean(df$ow)  # normalize

  # Weighted balance check
  w_smds <- sapply(avail, function(v) {
    x1 <- df[[v]][df$treated==1]; w1 <- df$ow[df$treated==1]
    x0 <- df[[v]][df$treated==0]; w0 <- df$ow[df$treated==0]
    m1 <- weighted.mean(x1, w1, na.rm=T)
    m0 <- weighted.mean(x0, w0, na.rm=T)
    sp <- sqrt((var(x1,na.rm=T)+var(x0,na.rm=T))/2)
    if (is.na(sp)||sp<1e-10) NA else abs(m1-m0)/sp
  })
  max_wsmd <- max(w_smds, na.rm=T)

  # OW-weighted regression
  fit <- lm(delta_cr ~ treated, data=df, weights=ow)
  cl <- vcovHC(fit, type="HC1")
  ct <- coeftest(fit, vcov.=cl)

  # Doubly robust: add covariates with weighted SMD > 0.05
  adj <- names(w_smds[!is.na(w_smds) & w_smds > 0.05])
  if (length(adj)>0) {
    fml <- as.formula(paste("delta_cr ~ treated +", paste(adj,collapse="+")))
    fit1 <- tryCatch(lm(fml, data=df, weights=ow), error=function(e) NULL)
    if (!is.null(fit1)) {
      cl1 <- vcovHC(fit1, type="HC1")
      ct1 <- coeftest(fit1, vcov.=cl1)
    } else ct1 <- ct
  } else ct1 <- ct

  list(n_trt=sum(df$treated==1), n_ctl=sum(df$treated==0),
       max_wsmd=round(max_wsmd,4),
       did_ow=ct["treated","Estimate"], p_ow=ct["treated","Pr(>|t|)"],
       did_owdr=ct1["treated","Estimate"], se_owdr=ct1["treated","Std. Error"],
       p_owdr=ct1["treated","Pr(>|t|)"],
       ci_lo=ct1["treated","Estimate"]-1.96*ct1["treated","Std. Error"],
       ci_hi=ct1["treated","Estimate"]+1.96*ct1["treated","Std. Error"])
}

# ============================================================================
run_fine <- function(db) {
  tag <- tolower(db)
  SEP <- paste(rep("=",70),collapse="")
  cat(sprintf("\n%s\n%s: Fine Timing Sweep + OW\n%s\n", SEP, db, SEP))

  trt <- read.csv(file.path(RESULTS,sprintf("20_did_treated_%s.csv",tag)),stringsAsFactors=F)
  ctl <- read.csv(file.path(RESULTS,sprintf("20_did_control_%s.csv",tag)),stringsAsFactors=F)
  cr_all <- read.csv(file.path(RESULTS,sprintf("20_did_cr_all_%s.csv",tag)),stringsAsFactors=F)

  id_col <- if("patientunitstayid" %in% names(trt)) "patientunitstayid" else "stay_id"
  trt$pid <- trt[[id_col]]; ctl$pid <- ctl[[id_col]]
  cr_id <- if("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
  cr_all$pid <- cr_all[[cr_id]]
  if(!"labresultoffset" %in% names(cr_all) && "offset_min" %in% names(cr_all))
    cr_all$labresultoffset <- cr_all$offset_min

  ps_vars <- intersect(PS_COVARS, intersect(names(trt),names(ctl)))
  stack_cols <- intersect(unique(c("pid","treated",ps_vars)),
                          intersect(names(trt),names(ctl)))

  combined <- rbind(trt[,stack_cols], ctl[,stack_cols])
  rownames(combined) <- seq_len(nrow(combined))
  combined <- median_impute(combined, ps_vars)
  cat(sprintf("  %d treated + %d control\n",
              sum(combined$treated==1), sum(combined$treated==0)))

  cr_list <- split(cr_all[,c("pid","labresult","labresultoffset")], cr_all$pid)
  mg_times <- setNames(trt$mg_offset_min, trt$pid)
  cr_pre_map <- setNames(trt$cr_pre, trt$pid)
  cr_pre_off_map <- setNames(trt$cr_pre_offset_min, trt$pid)

  # All post-Mg Cr for treated
  cr_trt <- cr_all[cr_all$pid %in% trt$pid, ]
  cr_trt$mg_off <- mg_times[as.character(cr_trt$pid)]
  cr_trt$post_h <- (cr_trt$labresultoffset - cr_trt$mg_off) / 60

  # ── PS match (r=1, cached for Part A) ────────────────────────────────
  cat("  PS matching (r=1)...\n")
  ps_formula <- as.formula(paste("treated ~", paste(ps_vars,collapse="+")))
  m <- suppressWarnings(matchit(ps_formula,data=combined,method="nearest",
                                 distance="glm",ratio=1,caliper=CALIPER,replace=TRUE))
  mm <- m$match.matrix
  trt_idx <- as.integer(rownames(mm))
  pairs <- list()
  for(i in seq_len(nrow(mm)))
    for(j in seq_len(ncol(mm))) {
      ci <- mm[i,j]
      if(!is.na(ci)) pairs[[length(pairs)+1]] <- c(combined$pid[trt_idx[i]],
                                                     combined$pid[as.integer(ci)])
    }
  pairs_df <- as.data.frame(do.call(rbind,pairs))
  names(pairs_df) <- c("trt_pid","ctl_pid")
  cat(sprintf("  Matched pairs: %d\n", nrow(pairs_df)))

  covar_want <- c("pid", ps_vars)
  covar_trt <- trt[,intersect(covar_want,names(trt)),drop=F]
  covar_ctl <- ctl[,intersect(covar_want,names(ctl)),drop=F]
  sh <- intersect(names(covar_trt),names(covar_ctl))
  covar_lu <- rbind(covar_trt[,sh],covar_ctl[,sh])
  covar_lu <- covar_lu[!duplicated(covar_lu$pid),]

  # ════════════════════════════════════════════════════════════════════════
  # PART A: Matching-based DiD sweep
  # ════════════════════════════════════════════════════════════════════════
  cat(sprintf("\n%s\nPART A: Matching-based DiD (r=1, targets %d-%dh, tol ±%s h)\n%s\n",
              paste(rep("─",60),collapse=""), min(TARGETS), max(TARGETS),
              paste(TOLERANCES,collapse="/"), paste(rep("─",60),collapse="")))

  match_results <- list(); midx <- 0

  for (target_h in TARGETS) {
    # Build Cr_post for this target
    win_lo <- max(0, target_h - 12); win_hi <- target_h + 12
    cw <- cr_trt[cr_trt$post_h >= win_lo & cr_trt$post_h <= win_hi, ]
    cw$dist <- abs(cw$post_h - target_h)
    strat <- cw[order(cw$pid, cw$dist),]; strat <- strat[!duplicated(strat$pid),]
    strat$cr_pre <- cr_pre_map[as.character(strat$pid)]
    strat$cr_pre_offset_min <- cr_pre_off_map[as.character(strat$pid)]
    strat$cr_post <- strat$labresult
    strat$cr_post_offset_min <- strat$labresultoffset
    strat$delta_cr <- strat$cr_post - strat$cr_pre

    trt_valid <- strat[,c("pid","cr_pre","cr_pre_offset_min","cr_post",
                           "cr_post_offset_min","delta_cr")]
    trt_valid$treated <- 1
    pv <- pairs_df[pairs_df$trt_pid %in% trt_valid$pid, ]

    for (tol_h in TOLERANCES) {
      tol_min <- tol_h * 60
      n_valid <- 0; ctl_rows <- list()

      for (r in seq_len(nrow(pv))) {
        tpid <- pv$trt_pid[r]; cpid <- pv$ctl_pid[r]
        tidx <- which(trt_valid$pid==tpid)[1]
        if(is.na(tidx)) next
        t_pre <- trt_valid$cr_pre_offset_min[tidx]
        t_post <- trt_valid$cr_post_offset_min[tidx]
        cc <- cr_list[[as.character(cpid)]]
        if(is.null(cc)||nrow(cc)<2) next
        pd <- abs(cc$labresultoffset-t_pre); pi <- which.min(pd)
        if(pd[pi]>tol_min) next
        pc <- cc[cc$labresultoffset > cc$labresultoffset[pi],]
        if(nrow(pc)==0) next
        pod <- abs(pc$labresultoffset-t_post); poi <- which.min(pod)
        if(pod[poi]>tol_min) next
        n_valid <- n_valid+1
        ctl_rows[[n_valid]] <- data.frame(
          pid=cpid,cr_pre=cc$labresult[pi],
          cr_pre_offset_min=cc$labresultoffset[pi],
          cr_post=pc$labresult[poi],
          cr_post_offset_min=pc$labresultoffset[poi],
          delta_cr=pc$labresult[poi]-cc$labresult[pi],
          treated=0,match_pair_id=tidx,stringsAsFactors=F)
      }

      if(n_valid<10) {
        midx <- midx+1
        match_results[[midx]] <- data.frame(method="matching",target_h=target_h,
          tol_h=tol_h,n_trt=0,n_ctl=0,did_adj=NA,p_adj=NA,ci_lo=NA,ci_hi=NA,
          stringsAsFactors=F)
        next
      }

      ctl_m <- do.call(rbind,ctl_rows)
      vp <- sort(unique(ctl_m$match_pair_id))
      t_out <- trt_valid[vp,]; t_out$match_pair_id <- seq_along(vp)
      ctl_m$match_pair_id <- match(ctl_m$match_pair_id,vp)
      mdf <- rbind(t_out, ctl_m)
      mdf <- merge(mdf, covar_lu, by="pid", all.x=T, suffixes=c("",".c"))
      mdf <- median_impute(mdf, ps_vars)

      res <- robust_did(mdf, ps_vars)
      midx <- midx+1
      if(!is.null(res)) {
        match_results[[midx]] <- data.frame(method="matching",target_h=target_h,
          tol_h=tol_h,n_trt=res$n_trt,n_ctl=res$n_ctl,
          did_adj=round(res$did_adj,4),p_adj=round(res$p_adj,4),
          ci_lo=round(res$ci_lo,4),ci_hi=round(res$ci_hi,4),
          stringsAsFactors=F)
      }
    }
  }

  mr <- do.call(rbind, match_results)

  # ════════════════════════════════════════════════════════════════════════
  # PART B: OW-DiD (no matching, common ICU-time anchor)
  # ════════════════════════════════════════════════════════════════════════
  cat(sprintf("\n%s\nPART B: Overlap Weighting DiD (common ICU-time anchor)\n%s\n",
              paste(rep("─",60),collapse=""), paste(rep("─",60),collapse="")))

  # For OW: everyone uses Cr closest to Xh from ICU admission
  # Cr_pre = first postop Cr (close to ICU admission) for everyone
  # Treated: cr_pre already defined (before IV Mg)
  # Controls: first_postop_cr already in their CSV

  ow_results <- list(); oidx <- 0

  for (target_h in TARGETS) {
    target_min <- target_h * 60

    # Treated: Cr closest to target_h from ICU admission (must be after Cr_pre)
    t_cr <- cr_all[cr_all$pid %in% trt$pid, ]
    t_cr$dist <- abs(t_cr$labresultoffset - target_min)
    # Must be after Cr_pre offset
    t_cr_valid <- merge(t_cr, data.frame(pid=trt$pid, cr_pre_off=trt$cr_pre_offset_min),
                        by="pid")
    t_cr_valid <- t_cr_valid[t_cr_valid$labresultoffset > t_cr_valid$cr_pre_off, ]
    t_post <- t_cr_valid[order(t_cr_valid$pid, t_cr_valid$dist), ]
    t_post <- t_post[!duplicated(t_post$pid), ]

    t_df <- merge(trt[,c("pid","cr_pre","treated")],
                  t_post[,c("pid","labresult")], by="pid")
    names(t_df)[names(t_df)=="labresult"] <- "cr_post"
    t_df$delta_cr <- t_df$cr_post - t_df$cr_pre

    # Controls: Cr closest to target_h from ICU
    # Cr_pre = first_postop_cr
    c_cr <- cr_all[cr_all$pid %in% ctl$pid, ]
    c_cr$dist <- abs(c_cr$labresultoffset - target_min)

    if ("first_cr_offset_min" %in% names(ctl)) {
      c_cr_valid <- merge(c_cr, data.frame(pid=ctl$pid,
                          cr_pre_off=ctl$first_cr_offset_min), by="pid")
      c_cr_valid <- c_cr_valid[c_cr_valid$labresultoffset > c_cr_valid$cr_pre_off, ]
    } else {
      c_cr_valid <- c_cr[c_cr$labresultoffset > 0, ]
    }
    c_post <- c_cr_valid[order(c_cr_valid$pid, c_cr_valid$dist), ]
    c_post <- c_post[!duplicated(c_post$pid), ]

    cr_pre_col <- if("first_postop_cr" %in% names(ctl)) "first_postop_cr" else "cr_pre"
    c_df <- merge(ctl[,c("pid",cr_pre_col,"treated")],
                  c_post[,c("pid","labresult")], by="pid")
    names(c_df)[names(c_df)==cr_pre_col] <- "cr_pre"
    names(c_df)[names(c_df)=="labresult"] <- "cr_post"
    c_df$delta_cr <- c_df$cr_post - c_df$cr_pre

    ow_df <- rbind(t_df[,c("pid","cr_pre","cr_post","delta_cr","treated")],
                   c_df[,c("pid","cr_pre","cr_post","delta_cr","treated")])
    ow_df <- merge(ow_df, covar_lu, by="pid", all.x=T, suffixes=c("",".c"))
    ow_df <- median_impute(ow_df, ps_vars)

    n_t <- sum(ow_df$treated==1); n_c <- sum(ow_df$treated==0)
    if (n_t < 20 || n_c < 20) next

    res <- ow_did(ow_df, ps_vars)
    oidx <- oidx+1
    if (!is.null(res)) {
      sig <- if(res$p_owdr < 0.05) " *" else ""
      cat(sprintf("  %3dh: n=%d/%d, OW max_wSMD=%.4f, DiD=%+.4f P=%.4f%s\n",
                  target_h, n_t, n_c, res$max_wsmd, res$did_owdr, res$p_owdr, sig))
      ow_results[[oidx]] <- data.frame(method="OW",target_h=target_h,tol_h=NA,
        n_trt=res$n_trt,n_ctl=res$n_ctl,max_wsmd=res$max_wsmd,
        did_adj=round(res$did_owdr,4),p_adj=round(res$p_owdr,4),
        ci_lo=round(res$ci_lo,4),ci_hi=round(res$ci_hi,4),
        stringsAsFactors=F)
    }
  }

  owr <- if(length(ow_results)>0) do.call(rbind, ow_results) else NULL

  # ════════════════════════════════════════════════════════════════════════
  # COMBINED SUMMARY
  # ════════════════════════════════════════════════════════════════════════
  cat(sprintf("\n%s\n%s: MATCHING-BASED DiD RESULTS\n%s\n", SEP, db, SEP))
  cat("\n  target  tol   n_trt  n_ctl    DiD      P       95% CI\n")
  cat("  ──────  ───  ─────  ─────  ────────  ──────  ──────────────\n")
  for(i in seq_len(nrow(mr))) {
    r <- mr[i,]
    if(is.na(r$did_adj)) {
      cat(sprintf("  %4dh   ±%dh      —      —       —       —\n",
                  r$target_h, r$tol_h))
    } else {
      sig <- if(r$p_adj<0.05) " *" else ""
      cat(sprintf("  %4dh   ±%dh  %5d  %5d  %+.4f  %.4f  [%+.4f,%+.4f]%s\n",
                  r$target_h,r$tol_h,r$n_trt,r$n_ctl,r$did_adj,r$p_adj,
                  r$ci_lo,r$ci_hi,sig))
    }
  }

  if(!is.null(owr)) {
    cat(sprintf("\n%s\n%s: OVERLAP WEIGHTING DiD RESULTS\n%s\n", SEP, db, SEP))
    cat("\n  target  n_trt   n_ctl   wSMD    DiD      P       95% CI\n")
    cat("  ──────  ──────  ──────  ──────  ────────  ──────  ──────────────\n")
    for(i in seq_len(nrow(owr))) {
      r <- owr[i,]
      sig <- if(r$p_adj<0.05) " *" else ""
      cat(sprintf("  %4dh   %5d   %5d  %.4f  %+.4f  %.4f  [%+.4f,%+.4f]%s\n",
                  r$target_h,r$n_trt,r$n_ctl,r$max_wsmd,r$did_adj,r$p_adj,
                  r$ci_lo,r$ci_hi,sig))
    }
  }

  # Save all
  all_res <- mr
  if(!is.null(owr)) {
    owr_aligned <- owr
    if(!"tol_h" %in% names(owr_aligned)) owr_aligned$tol_h <- NA
    if(!"max_wsmd" %in% names(mr)) mr$max_wsmd <- NA
    common <- intersect(names(mr), names(owr_aligned))
    all_res <- rbind(mr[,common], owr_aligned[,common])
  }
  write.csv(all_res, file.path(RESULTS,sprintf("24_fine_sweep_%s.csv",tag)), row.names=F)
  cat(sprintf("\n  Saved: 24_fine_sweep_%s.csv\n", tag))
}

# ============================================================================
cat("======================================================================\n")
cat("24_did_fine_sweep.R — Fine timing + OW parallel\n")
cat(sprintf("  Targets: %s h (post-IV-Mg for matching, post-ICU for OW)\n",
            paste(TARGETS,collapse=",")))
cat(sprintf("  Tolerances: ±%s h\n", paste(TOLERANCES,collapse=",")))
cat("======================================================================\n")

args <- commandArgs(trailingOnly=TRUE)
if(length(args)==0){cat("Usage: Rscript 24_did_fine_sweep.R eicu|mimic\n");quit(status=1)}
for(a in args) run_fine(toupper(a))

cat("\n======================================================================\n")
cat("Compare matching vs OW: if they agree, the result is robust.\n")
cat("OW gives exact balance (wSMD→0) without temporal matching.\n")
cat("======================================================================\n")
