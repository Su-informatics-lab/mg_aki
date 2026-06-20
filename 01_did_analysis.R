#!/usr/bin/env Rscript
# ============================================================================
# 01_did_analysis.R — PSM + MICE (JNO-grade)
#
# Primary:       MICE m=20, 21 PS covars, PSM 1:1, Rubin's rules
# Sensitivity 1: Reduced model (<30% missing covars only)
# Sensitivity 2: Complete case
#
# Usage: Rscript 01_did_analysis.R eicu [primary|reduced|complete]
# ============================================================================

suppressPackageStartupMessages({
  library(MatchIt); library(sandwich); library(lmtest); library(mice)
})
source(file.path(path.expand("~/mg_aki"), "did_covars.R"))

RESULTS   <- path.expand("~/mg_aki/results")
CALIPER   <- 0.2
M_IMP     <- 20
PRIMARY_T <- 24
WINDOW    <- 6

SEC_OC <- c("hosp_mortality","poaf","encephalopathy","vent_arrhythmia")
SEC_LB <- c(hosp_mortality="Hospital mortality",poaf="POAF (new-onset)",
            encephalopathy="Encephalopathy",vent_arrhythmia="Ventricular arrhythmia")

# ── Helpers ──────────────────────────────────────────────────────────────
compute_smds <- function(d, vars) {
  sapply(vars, function(v) {
    if (!v %in% names(d)) return(NA)
    x1 <- d[[v]][d$treated==1]; x0 <- d[[v]][d$treated==0]
    sp <- sqrt((var(x1,na.rm=T)+var(x0,na.rm=T))/2)
    if (is.na(sp)||sp<1e-10) NA else abs(mean(x1,na.rm=T)-mean(x0,na.rm=T))/sp
  })
}

did_dr <- function(df, ps_vars) {
  nt <- sum(df$treated==1); nc <- sum(df$treated==0)
  if (nt<10||nc<10) return(NULL)
  smds <- compute_smds(df, ps_vars)
  adj <- names(smds[!is.na(smds) & smds>0.05])
  adj <- intersect(adj, names(df))
  adj <- adj[vapply(adj, function(v) var(df[[v]],na.rm=T)>1e-10, logical(1))]
  if (length(adj)>0)
    fml <- as.formula(paste("delta_cr ~ treated +", paste(adj,collapse="+")))
  else fml <- delta_cr ~ treated
  fit <- lm(fml, data=df)
  ct <- coeftest(fit, vcov.=tryCatch(vcovHC(fit,type="HC1"),error=function(e) vcov(fit)))
  list(did=ct["treated","Estimate"], se=ct["treated","Std. Error"],
       p=ct["treated","Pr(>|t|)"], n_trt=nt, n_ctl=nc)
}

build_dcr <- function(data_with_pid, cr_all, target_h) {
  pre <- cr_all[cr_all$offset_h>=0 & cr_all$offset_h<=6,]
  pre <- pre[order(pre$pid,pre$offset_h),]; pre <- pre[!duplicated(pre$pid),]
  post <- cr_all[cr_all$offset_h>=(target_h-WINDOW) & cr_all$offset_h<=(target_h+WINDOW),]
  post$dist <- abs(post$offset_h - target_h)
  post <- post[order(post$pid,post$dist),]; post <- post[!duplicated(post$pid),]
  m <- merge(pre[,c("pid","labresult","offset_h")],
             post[,c("pid","labresult","offset_h")],
             by="pid", suffixes=c("_pre","_post"))
  m <- m[m$offset_h_post > m$offset_h_pre,]
  m$delta_cr <- m$labresult_post - m$labresult_pre
  merge(data_with_pid, m[,c("pid","delta_cr")], by="pid")
}

rubin_pool <- function(estimates, ses) {
  m <- length(estimates)
  theta <- mean(estimates)
  W <- mean(ses^2)
  B <- var(estimates)
  T_var <- W + (1 + 1/m) * B
  se <- sqrt(T_var)
  r <- (1 + 1/m) * B / W
  df_old <- (m - 1) * (1 + 1/r)^2
  p <- 2 * pt(-abs(theta/se), df=max(df_old, 2))
  list(theta=theta, se=se, p=p,
       ci_lo=theta-1.96*se, ci_hi=theta+1.96*se,
       fmi=round((r+2/(df_old+3))/(r+1), 3))
}

run_single_psm <- function(data, cr_all, ps_vars, target_h=PRIMARY_T) {
  ps_avail <- intersect(ps_vars, names(data))
  ps_fml <- as.formula(paste("treated ~", paste(ps_avail,collapse="+")))
  m <- tryCatch(
    suppressWarnings(matchit(ps_fml, data=data, method="nearest",
                              distance="glm", ratio=1, caliper=CALIPER, replace=TRUE)),
    error=function(e) { cat("    matchit error:", conditionMessage(e), "\n"); NULL })
  if (is.null(m)) return(NULL)
  md <- match.data(m)
  smds_match <- compute_smds(md, ps_avail)
  dcr <- build_dcr(md, cr_all, target_h)
  if (nrow(dcr) < 40) return(NULL)
  res <- did_dr(dcr, ps_avail)
  if (is.null(res)) return(NULL)
  list(did=res$did, se=res$se, p=res$p,
       n_trt=res$n_trt, n_ctl=res$n_ctl,
       max_smd=max(smds_match, na.rm=TRUE),
       n_above_01=sum(smds_match>0.1, na.rm=TRUE))
}

# ============================================================================
run_analysis <- function(db, mode="primary") {
  tag <- tolower(db)
  SEP <- paste(rep("=",70), collapse="")
  cat(sprintf("\n%s\n%s [%s]: PSM + MICE (m=%d)\n%s\n", SEP, db, mode, M_IMP, SEP))

  # ── Load ────────────────────────────────────────────────────────────
  trt <- read.csv(file.path(RESULTS,sprintf("did_treated_%s.csv",tag)),stringsAsFactors=F)
  ctl <- read.csv(file.path(RESULTS,sprintf("did_control_%s.csv",tag)),stringsAsFactors=F)
  cr_all <- read.csv(file.path(RESULTS,sprintf("did_cr_all_%s.csv",tag)),stringsAsFactors=F)

  id_col <- if("patientunitstayid" %in% names(trt)) "patientunitstayid" else "stay_id"
  trt$pid <- trt[[id_col]]; ctl$pid <- ctl[[id_col]]
  cr_id <- if("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
  cr_all$pid <- cr_all[[cr_id]]
  if(!"labresultoffset" %in% names(cr_all)) cr_all$labresultoffset <- cr_all$offset_min
  cr_all$offset_h <- cr_all$labresultoffset / 60

  all_want <- unique(c("pid","treated",PS_PRIMARY,"surgery_type","first_mg_value",SEC_OC))
  sh <- intersect(all_want, intersect(names(trt),names(ctl)))
  combined <- rbind(trt[,sh], ctl[,sh])

  cat(sprintf("  Loaded: %d treated + %d control\n",
              sum(combined$treated==1), sum(combined$treated==0)))

  # ── Select PS vars ─────────────────────────────────────────────────
  if (mode == "reduced") {
    miss_rates <- sapply(PS_PRIMARY, function(v)
      if(v %in% names(combined)) mean(is.na(combined[[v]])) else 1)
    ps_vars <- names(miss_rates[miss_rates < 0.30])
    cat(sprintf("  Reduced: %d covars (dropped: %s)\n",
                length(ps_vars), paste(setdiff(PS_PRIMARY,ps_vars),collapse=", ")))
  } else {
    ps_vars <- PS_PRIMARY
  }

  # ── Complete case ──────────────────────────────────────────────────
  if (mode == "complete") {
    cat("\n  ── COMPLETE CASE ANALYSIS ──\n")
    ps_avail <- intersect(ps_vars, names(combined))
    cc <- combined[complete.cases(combined[, ps_avail]), ]
    cat(sprintf("  Complete: %d/%d (%.1f%%), treated=%d, control=%d\n",
                nrow(cc), nrow(combined), 100*nrow(cc)/nrow(combined),
                sum(cc$treated==1), sum(cc$treated==0)))
    res <- run_single_psm(cc, cr_all, ps_vars)
    if (!is.null(res))
      cat(sprintf("  PSM-DR: DiD=%+.4f (SE=%.4f), P=%.4f, 95%%CI [%+.4f,%+.4f]\n",
                  res$did, res$se, res$p, res$did-1.96*res$se, res$did+1.96*res$se))
    return(invisible(NULL))
  }

  # ── MICE ───────────────────────────────────────────────────────────
  cat(sprintf("\n  ── MICE (m=%d) ──\n", M_IMP))

  # Build imputation data: PS covars + treated (NO delta_cr_aux)
  imp_vars <- unique(c(intersect(ps_vars, names(combined)), "treated"))
  imp_data <- combined[, imp_vars]

  # Set methods
  ini <- mice(imp_data, maxit=0, print=FALSE)
  meth <- ini$method
  meth["treated"] <- ""  # don't impute treatment
  for (v in names(meth)) {
    if (meth[v] == "") next
    if (all(imp_data[[v]] %in% c(0,1,NA))) meth[v] <- "logreg"
    else meth[v] <- "pmm"
  }

  cat(sprintf("  Methods: %s\n",
              paste(sprintf("%s=%s", names(meth[meth!=""]), meth[meth!=""]), collapse=", ")))
  cat("  Running MICE...")
  imp <- mice(imp_data, m=M_IMP, method=meth, maxit=10, seed=42, printFlag=FALSE)
  cat(" done.\n")
  cat(sprintf("  Logged events: %d\n", nrow(imp$loggedEvents)))

  # ── PSM on each imputed dataset ────────────────────────────────────
  cat(sprintf("\n  ── PSM 1:1 on %d imputed datasets ──\n", M_IMP))

  dids <- ses <- smds_v <- n_trts <- n_ctls <- numeric(M_IMP)
  valid <- 0

  for (i in seq_len(M_IMP)) {
    imp_i <- complete(imp, i)

    # Add back pid and non-PS columns
    imp_i$pid <- combined$pid
    for (col in setdiff(names(combined), names(imp_i))) {
      imp_i[[col]] <- combined[[col]]
    }

    res <- tryCatch(run_single_psm(imp_i, cr_all, ps_vars), error=function(e) {
      cat(sprintf("    m=%d error: %s\n", i, conditionMessage(e))); NULL })

    if (!is.null(res)) {
      valid <- valid + 1
      dids[valid] <- res$did; ses[valid] <- res$se
      smds_v[valid] <- res$max_smd
      n_trts[valid] <- res$n_trt; n_ctls[valid] <- res$n_ctl
      if (i <= 3 || i == M_IMP)
        cat(sprintf("    m=%2d: DiD=%+.4f, SE=%.4f, maxSMD=%.3f, n=%d/%d\n",
                    i, res$did, res$se, res$max_smd, res$n_trt, res$n_ctl))
      else if (i == 4) cat("    ...\n")
    } else {
      cat(sprintf("    m=%2d: FAILED\n", i))
    }
  }

  if (valid < 3) { cat("  Too few valid imputations\n"); return(NULL) }
  dids <- dids[1:valid]; ses <- ses[1:valid]

  # ── Pool ───────────────────────────────────────────────────────────
  cat(sprintf("\n  ── POOLED (Rubin's rules, %d/%d valid) ──\n", valid, M_IMP))
  pool <- rubin_pool(dids, ses)

  cat(sprintf("  DiD = %+.4f (SE=%.4f), P=%.4f\n", pool$theta, pool$se, pool$p))
  cat(sprintf("  95%% CI: [%+.4f, %+.4f]\n", pool$ci_lo, pool$ci_hi))
  cat(sprintf("  FMI: %.3f\n", pool$fmi))
  cat(sprintf("  Mean maxSMD: %.3f, Mean n_trt: %d\n",
              mean(smds_v[1:valid]), round(mean(n_trts[1:valid]))))

  # ── Time course (pooled) ───────────────────────────────────────────
  cat(sprintf("\n  ── TIME COURSE ──\n"))
  cat("  target   DiD        SE      P        95%% CI               FMI\n")
  cat("  ──────  ────────  ──────  ──────  ────────────────────  ─────\n")

  tc_rows <- list()
  for (th in seq(6, 36, by=6)) {
    tc_d <- tc_s <- numeric(valid); tc_v <- 0
    for (i in seq_len(valid)) {
      imp_i <- complete(imp, i)
      imp_i$pid <- combined$pid
      for (col in setdiff(names(combined), names(imp_i)))
        imp_i[[col]] <- combined[[col]]
      r <- tryCatch(run_single_psm(imp_i, cr_all, ps_vars, target_h=th),
                    error=function(e) NULL)
      if (!is.null(r)) { tc_v <- tc_v+1; tc_d[tc_v] <- r$did; tc_s[tc_v] <- r$se }
    }
    if (tc_v >= 3) {
      tp <- rubin_pool(tc_d[1:tc_v], tc_s[1:tc_v])
      sig <- if(tp$p<0.05) " *" else ""
      tag_str <- if(th==PRIMARY_T) " << PRIMARY" else ""
      cat(sprintf("  %4dh   %+.4f  %.4f  %.4f  [%+.4f,%+.4f]  %.3f%s%s\n",
                  th, tp$theta, tp$se, tp$p, tp$ci_lo, tp$ci_hi, tp$fmi, sig, tag_str))
      tc_rows[[length(tc_rows)+1]] <- data.frame(
        target_h=th, did=round(tp$theta,5), se=round(tp$se,5),
        p=round(tp$p,4), ci_lo=round(tp$ci_lo,5), ci_hi=round(tp$ci_hi,5),
        fmi=tp$fmi, stringsAsFactors=F)
    }
  }
  if (length(tc_rows)>0)
    write.csv(do.call(rbind, tc_rows),
              file.path(RESULTS, sprintf("did_timecourse_%s_mice.csv",tag)), row.names=F)

  # ── Subgroups (on primary 24h, first imputed dataset) ──────────────
  cat(sprintf("\n  ── SUBGROUPS (24h, imputation m=1) ──\n"))
  imp1 <- complete(imp, 1)
  imp1$pid <- combined$pid
  for (col in setdiff(names(combined), names(imp1))) imp1[[col]] <- combined[[col]]

  ps_avail <- intersect(ps_vars, names(imp1))
  ps_fml <- as.formula(paste("treated ~", paste(ps_avail,collapse="+")))
  m1 <- suppressWarnings(matchit(ps_fml, data=imp1, method="nearest",
                                  distance="glm", ratio=1, caliper=CALIPER, replace=TRUE))
  md1 <- match.data(m1)
  d24 <- build_dcr(md1, cr_all, PRIMARY_T)

  if ("first_mg_value" %in% names(d24)) {
    cat("  Mg strata:\n")
    for (mg_nm in c("<1.8","1.8-2.0","2.0-2.3",">2.3")) {
      bins <- list("<1.8"=c(0,1.8),"1.8-2.0"=c(1.8,2.0),"2.0-2.3"=c(2.0,2.3),">2.3"=c(2.3,99))
      lo <- bins[[mg_nm]][1]; hi <- bins[[mg_nm]][2]
      sub <- d24[!is.na(d24$first_mg_value) & d24$first_mg_value>=lo & d24$first_mg_value<hi,]
      r <- if(sum(sub$treated==1)>=10 && sum(sub$treated==0)>=10) did_dr(sub,ps_avail) else NULL
      if(!is.null(r)) cat(sprintf("    %s: n=%d/%d, DiD=%+.4f, P=%.4f%s\n",
                                   mg_nm,r$n_trt,r$n_ctl,r$did,r$p,if(r$p<0.05)" *"else""))
    }
  }
  if ("surgery_type" %in% names(d24)) {
    cat("  Surgery:\n")
    for (st in c("cabg","valve","combined","other_cardiac")) {
      sub <- d24[!is.na(d24$surgery_type) & d24$surgery_type==st,]
      r <- if(sum(sub$treated==1)>=10 && sum(sub$treated==0)>=10) did_dr(sub,ps_avail) else NULL
      if(!is.null(r)) cat(sprintf("    %s: n=%d/%d, DiD=%+.4f, P=%.4f%s\n",
                                   st,r$n_trt,r$n_ctl,r$did,r$p,if(r$p<0.05)" *"else""))
    }
  }

  # ── Secondary outcomes (matched set from m=1) ─────────────────────
  cat(sprintf("\n  ── SECONDARY OUTCOMES (matched m=1) ──\n"))
  cat("  outcome                n_trt  n_ctl  rate_trt rate_ctl  OR (95%CI)         P\n")
  cat("  --------------------- ------ ------ -------- -------- -----------------  ------\n")
  for (oc in SEC_OC) {
    if (!oc %in% names(md1)) next
    md1$y <- as.numeric(md1[[oc]]); md_c <- md1[!is.na(md1$y),]
    nt <- sum(md_c$treated==1); nc <- sum(md_c$treated==0)
    if (nt<20||nc<20) next
    r1 <- mean(md_c$y[md_c$treated==1]); r0 <- mean(md_c$y[md_c$treated==0])
    fit <- tryCatch(glm(y~treated, data=md_c, family=quasibinomial()), error=function(e) NULL)
    if (is.null(fit)||!"treated" %in% names(coef(fit))) next
    ct <- tryCatch(coeftest(fit, vcov.=vcovHC(fit,type="HC1")), error=function(e) NULL)
    if (is.null(ct)) next
    tr <- which(rownames(ct)=="treated")
    if (length(tr)==0) next
    or <- exp(ct[tr,"Estimate"])
    or_lo <- exp(ct[tr,"Estimate"]-1.96*ct[tr,"Std. Error"])
    or_hi <- exp(ct[tr,"Estimate"]+1.96*ct[tr,"Std. Error"])
    p <- ct[tr,grep("^Pr",colnames(ct))]
    sig <- if(p<0.05) " *" else ""
    label <- ifelse(oc %in% names(SEC_LB), SEC_LB[oc], oc)
    cat(sprintf("  %-23s %5d  %5d   %.1f%%   %.1f%%  %.2f (%.2f-%.2f)  %.4f%s\n",
                label,nt,nc,100*r1,100*r0,or,or_lo,or_hi,p,sig))
  }

  # ── Save ───────────────────────────────────────────────────────────
  write.csv(data.frame(database=db, mode=mode, did=round(pool$theta,5),
    se=round(pool$se,5), p=round(pool$p,4), ci_lo=round(pool$ci_lo,5),
    ci_hi=round(pool$ci_hi,5), fmi=pool$fmi,
    mean_smd=round(mean(smds_v[1:valid]),3), m_valid=valid,
    stringsAsFactors=F),
    file.path(RESULTS, sprintf("did_mice_%s_%s.csv",tag,mode)), row.names=F)

  # Save matched set from m=1 for downstream
  write.csv(d24, file.path(RESULTS, sprintf("did_matched_%s_mice_24h.csv",tag)), row.names=F)

  # Save balance from m=1
  smds_raw <- compute_smds(combined, ps_avail)
  smds_m <- compute_smds(md1, ps_avail)
  write.csv(data.frame(covariate=names(smds_raw),
    raw_smd=round(as.numeric(smds_raw),4),
    matched_smd=round(as.numeric(smds_m),4), stringsAsFactors=F),
    file.path(RESULTS, sprintf("did_balance_%s_mice.csv",tag)), row.names=F)

  cat(sprintf("\n%s\n%s [%s]: COMPLETE\n%s\n", SEP, db, mode, SEP))
}

# ============================================================================
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  cat("Usage: Rscript 01_did_analysis.R <db> [primary|reduced|complete]\n")
  quit(status=1)
}
db_arg <- toupper(args[1])
mode_arg <- if(length(args)>=2) tolower(args[2]) else "primary"
cat("======================================================================\n")
cat(sprintf("01_did_analysis.R — MICE (m=%d) + PSM, mode=%s\n", M_IMP, mode_arg))
cat("======================================================================\n")
run_analysis(db_arg, mode_arg)
