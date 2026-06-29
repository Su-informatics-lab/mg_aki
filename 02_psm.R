#!/usr/bin/env Rscript
# ============================================================================
# 02_psm.R — Canonical Risk-Set PSM (v5.1 — yet-untreated only)
#
# Specs:
#   PRIMARY:  19 var, LAST labs, no K/Mg (best balance, avoids positivity issue)
#   SENS_A:   21 var, LAST labs, + K/Mg  (over-adjustment check)
#   SENS_B:   19 var, FIRST labs, no K/Mg (lab-timing sensitivity)
#
# Pool:     yet-untreated only (sequential trial estimand)
#           Never-treated pool removed: healthy-control bias, does not
#           answer "should we give Mg now?" estimand.
# Methods:  PSM | PSM+DR (adjust SMD>0.1)
# Horizons: 6–48h
# MICE:     m=20 averaged
# Match:    1:1 with replacement, caliper 0.2 SD, HC1 SE
#
# Usage: Rscript 02_psm.R mimic
#        Rscript 02_psm.R eicu
# ============================================================================

suppressPackageStartupMessages({ library(sandwich); library(lmtest); library(mice) })

RESULTS    <- path.expand("~/mg_aki/results")
PRIMARY_H  <- 24
CR_WINDOW  <- 12
CALIPER_SD <- 0.2
M_IMP      <- 20
TARGETS    <- c(6, 12, 18, 24, 30, 36, 42, 48)

PS_BASE <- c("age","is_female","bmi",
             "surg_cabg","surg_valve","surg_combined",
             "heart_failure","hypertension","diabetes","ckd",
             "copd","pvd","stroke","liver_disease","egfr")

SPECS <- list(
  primary = list(
    vars  = c(PS_BASE, "last_calcium","last_lactate","last_lactate_missing","last_heartrate"),
    label = "19var LAST (no K/Mg)"),
  sens_a  = list(
    vars  = c(PS_BASE, "last_magnesium","last_potassium","last_calcium",
              "last_lactate","last_lactate_missing","last_heartrate"),
    label = "21var LAST (all labs)"),
  sens_b  = list(
    vars  = c(PS_BASE, "first_calcium","first_lactate","first_lactate_missing","first_heartrate"),
    label = "19var FIRST (no K/Mg)")
)

LAB_BASES <- c("magnesium","potassium","calcium","lactate","heartrate")

# ═══════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════
find_cr <- function(cr_pt, target_h, window = CR_WINDOW) {
  if (is.null(cr_pt) || nrow(cr_pt) == 0) return(c(NA, NA))
  cand <- cr_pt[cr_pt$offset_h >= (target_h - window) & cr_pt$offset_h <= (target_h + window), ]
  if (nrow(cand) == 0) return(c(NA, NA))
  best <- cand[which.min(abs(cand$offset_h - target_h)), ]
  c(best$labresult, best$offset_h)
}
find_cr_pre <- function(cr_pt, t_h) {
  if (is.null(cr_pt) || nrow(cr_pt) == 0) return(c(NA, NA))
  cand <- cr_pt[cr_pt$offset_h >= 0 & cr_pt$offset_h < t_h, ]
  if (nrow(cand) == 0) return(c(NA, NA))
  best <- cand[which.max(cand$offset_h), ]
  c(best$labresult, best$offset_h)
}
safe_coeftest <- function(fit) {
  ct <- tryCatch(suppressWarnings(coeftest(fit, vcov. = vcovHC(fit, type = "HC1"))),
                 error = function(e) NULL)
  if (!is.null(ct) && is.matrix(ct) && "treated" %in% rownames(ct) &&
      ncol(ct) >= 4 && !any(is.nan(ct["treated", ]))) return(ct)
  tryCatch(coeftest(fit), error = function(e) NULL)
}

extract_labs <- function(labs, all_pts, lab_bases, prefix, descending) {
  for (ln in lab_bases) {
    sub <- labs[labs$lab_name == ln, ]
    if (nrow(sub) == 0) next
    sub$mg_offset_h <- all_pts$mg_offset_h[match(sub$pid, all_pts$pid)]
    sub <- sub[sub$offset_h >= 0 & (is.na(sub$mg_offset_h) | sub$offset_h < sub$mg_offset_h), ]
    if (nrow(sub) == 0) next
    if (descending) s <- sub[order(-sub$offset_h), ] else s <- sub[order(sub$offset_h), ]
    s <- s[!duplicated(s$pid), ]
    idx <- match(all_pts$pid, s$pid)
    col <- paste0(prefix, "_", ln)
    all_pts[[col]] <- s$value[idx]
  }
  lac_col <- paste0(prefix, "_lactate")
  miss_col <- paste0(prefix, "_lactate_missing")
  all_pts[[miss_col]] <- as.integer(is.na(all_pts[[lac_col]]))
  all_pts
}

# ═══════════════════════════════════════════════════════════════════
# run_spec_pool
# ═══════════════════════════════════════════════════════════════════
run_spec_pool <- function(spec_name, spec_obj, pool_name,
                          all_pts, trt_idx, risk_sets,
                          cr_list, trt_tmg, caliper) {
  ps_vars <- intersect(spec_obj$vars, names(all_pts))
  ps_vars <- ps_vars[vapply(ps_vars, function(v) {
    x <- all_pts[[v]]; !all(is.na(x)) && var(x, na.rm=TRUE) > 1e-10
  }, logical(1))]
  n_trt <- length(trt_idx); trt_pids <- all_pts$pid[trt_idx]

  sep <- paste(rep("-", 60), collapse = "")
  cat(sprintf("\n%s\n  %s [%s] | %s | %d covars\n%s\n",
              sep, toupper(spec_name), spec_obj$label, toupper(pool_name),
              length(ps_vars), sep))

  ps_fml <- as.formula(paste("treated ~", paste(ps_vars, collapse = "+")))
  ps_fit <- suppressWarnings(glm(ps_fml, data = all_pts, family = binomial()))
  all_pts$ps <- predict(ps_fit, type = "response")

  match_trt <- integer(n_trt); match_ctl <- integer(n_trt); matched <- logical(n_trt)
  for (k in seq_len(n_trt)) {
    rs <- risk_sets[[k]]; if (length(rs) == 0) next
    pd <- abs(all_pts$ps[rs] - all_pts$ps[trt_idx[k]])
    wc <- which(pd <= caliper); if (length(wc) == 0) next
    best <- wc[which.min(pd[wc])]
    match_trt[k] <- trt_idx[k]; match_ctl[k] <- rs[best]; matched[k] <- TRUE
  }
  nm <- sum(matched)
  cat(sprintf("  Matched: %d/%d (%.1f%%)\n", nm, n_trt, 100*nm/n_trt))
  if (nm < 50) return(NULL)

  trt_m <- match_trt[matched]; ctl_m <- match_ctl[matched]
  smds <- sapply(ps_vars, function(v) {
    x1 <- all_pts[[v]][trt_m]; x0 <- all_pts[[v]][ctl_m]
    sp <- sqrt((var(x1,na.rm=T)+var(x0,na.rm=T))/2)
    if (is.na(sp)||sp<1e-10) 0 else abs(mean(x1,na.rm=T)-mean(x0,na.rm=T))/sp
  })
  adj_vars <- names(smds[smds > 0.1])
  cat(sprintf("  Balance: max=%.3f, viol=%d/%d\n", max(smds), sum(smds>0.1), length(ps_vars)))
  if (length(adj_vars)>0) cat(sprintf("    DR adjusts: %s\n", paste(adj_vars, collapse=", ")))

  results <- list()
  for (target_h in TARGETS) {
    dcr_t <- dcr_c <- numeric(nm); valid <- logical(nm); idx <- 0
    for (kk in which(matched)) {
      idx <- idx+1
      tp <- as.character(all_pts$pid[match_trt[kk]])
      cp <- as.character(all_pts$pid[match_ctl[kk]])
      tmg <- trt_tmg[kk]
      # Yet-untreated: censor if control received Mg before outcome window
      cmg <- all_pts$mg_offset_h[match_ctl[kk]]
      if (!is.na(cmg) && cmg < tmg+target_h) { valid[idx]<-FALSE; next }
      pt <- find_cr_pre(cr_list[[tp]], tmg); pc <- find_cr_pre(cr_list[[cp]], tmg)
      qt <- find_cr(cr_list[[tp]], tmg+target_h); qc <- find_cr(cr_list[[cp]], tmg+target_h)
      if (any(is.na(c(pt[1],pc[1],qt[1],qc[1])))) { valid[idx]<-FALSE; next }
      dcr_t[idx] <- qt[1]-pt[1]; dcr_c[idx] <- qc[1]-pc[1]; valid[idx] <- TRUE
    }
    nv <- sum(valid)
    if (nv < 30) {
      for (mt in c("psm","psm_dr"))
        results[[length(results)+1]] <- data.frame(spec=spec_name, pool=pool_name,
          target_h=target_h, method=mt, n=nv, did=NA, se=NA, p=NA, ci_lo=NA, ci_hi=NA,
          max_smd=max(smds), n_viol=sum(smds>0.1))
      next
    }
    pdf <- data.frame(delta_cr=c(dcr_t[valid],dcr_c[valid]), treated=rep(c(1,0),each=nv))
    tr <- match_trt[matched][valid]; cr <- match_ctl[matched][valid]
    for (av in adj_vars) if (av %in% names(all_pts))
      pdf[[av]] <- c(all_pts[[av]][tr], all_pts[[av]][cr])

    # PSM plain
    fp <- lm(delta_cr ~ treated, data=pdf); cp <- safe_coeftest(fp)
    results[[length(results)+1]] <- data.frame(spec=spec_name, pool=pool_name,
      target_h=target_h, method="psm", n=nv,
      did=cp["treated",1], se=cp["treated",2], p=cp["treated",4],
      ci_lo=cp["treated",1]-1.96*cp["treated",2], ci_hi=cp["treated",1]+1.96*cp["treated",2],
      max_smd=max(smds), n_viol=sum(smds>0.1))

    # PSM+DR
    ua <- intersect(adj_vars, names(pdf))
    ua <- ua[vapply(ua, function(v) var(pdf[[v]],na.rm=T)>1e-10, logical(1))]
    fd <- if (length(ua)>0) tryCatch(lm(as.formula(paste("delta_cr~treated+",paste(ua,collapse="+"))),
           data=pdf), error=function(e) fp) else fp
    cd <- safe_coeftest(fd)
    results[[length(results)+1]] <- data.frame(spec=spec_name, pool=pool_name,
      target_h=target_h, method="psm_dr", n=nv,
      did=cd["treated",1], se=cd["treated",2], p=cd["treated",4],
      ci_lo=cd["treated",1]-1.96*cd["treated",2], ci_hi=cd["treated",1]+1.96*cd["treated",2],
      max_smd=max(smds), n_viol=sum(smds>0.1))
  }
  res <- do.call(rbind, results)

  cat(sprintf("  ── PSM_DR [%s | %s] ──\n", spec_name, pool_name))
  sub <- res[res$method=="psm_dr",]
  for (i in seq_len(nrow(sub))) { r<-sub[i,]
    if (is.na(r$did)) { cat(sprintf("  %3dh (n=%d)\n",r$target_h,r$n)); next }
    sig<-if(!is.na(r$p)&&r$p<0.05)" *" else "  "
    pri<-if(r$target_h==PRIMARY_H)" <<" else ""
    cat(sprintf("  %3dh  %+.4f  P=%.4f%s  [%+.4f,%+.4f]  n=%d%s\n",
                r$target_h,r$did,r$p,sig,r$ci_lo,r$ci_hi,r$n,pri))
  }
  list(results=res, pairs=data.frame(trt_pid=all_pts$pid[trt_m],ctl_pid=all_pts$pid[ctl_m],
       t_mg=trt_tmg[matched]), n_matched=nm, max_smd=max(smds), n_viol=sum(smds>0.1))
}

# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<1) { cat("Usage: Rscript 02_psm.R <db>\n"); quit(status=1) }
db <- toupper(args[1]); tag <- tolower(db)

SEP <- paste(rep("=",70),collapse="")
cat(sprintf("\n%s\n02_psm.R — %s\n  PRIMARY: 19var LAST (no K/Mg)\n  SENS_A:  21var LAST (all labs)\n  SENS_B:  19var FIRST (no K/Mg)\n  Pool: yet-untreated only | Methods: PSM + DR | Horizons: 6-48h\n%s\n",SEP,db,SEP))

all_pts <- read.csv(file.path(RESULTS,sprintf("did_all_%s.csv",tag)),stringsAsFactors=FALSE)
cr_all <- read.csv(file.path(RESULTS,sprintf("did_cr_all_%s.csv",tag)),stringsAsFactors=FALSE)
cr_id <- if("patientunitstayid"%in%names(cr_all)) "patientunitstayid" else "stay_id"
cr_all$pid <- cr_all[[cr_id]]
if(!"offset_h"%in%names(cr_all)) cr_all$offset_h <- cr_all$labresultoffset/60
cr_all <- cr_all[order(cr_all$pid, cr_all$offset_h),]

cat("  Loading labs...\n")
labs_raw <- read.csv(file.path(RESULTS,sprintf("did_labs_all_%s.csv",tag)),stringsAsFactors=FALSE)
pcl <- if("patientunitstayid"%in%names(labs_raw)) "patientunitstayid" else "stay_id"
if(pcl%in%names(labs_raw)) labs_raw$pid <- labs_raw[[pcl]]
le <- labs_raw[labs_raw$lab_name%in%c("magnesium","potassium","calcium","lactate"),]
lh <- labs_raw[labs_raw$lab_name=="heartrate",]
if(nrow(lh)>500000) {
  lh$hb<-floor(lh$offset_h); lh<-lh[order(lh$pid,lh$offset_h),]
  lh<-lh[!duplicated(paste(lh$pid,lh$hb)),]; lh$hb<-NULL
}
labs<-rbind(le,lh); rm(labs_raw,le,lh); gc()

N<-nrow(all_pts)
cat(sprintf("  Patients: %d (%d trt, %d ctl)\n\n",N,sum(all_pts$treated==1),sum(all_pts$treated==0)))

# ── Compute LAST and FIRST labs ───────────────────────────────────
cat("  LAST labs (closest to T0)...\n")
all_pts <- extract_labs(labs, all_pts, LAB_BASES, "last", descending=TRUE)
for(ln in LAB_BASES) { col<-paste0("last_",ln); nf<-sum(!is.na(all_pts[[col]]))
  cat(sprintf("    %s: %d (%.0f%%)\n",col,nf,100*nf/N)) }

cat("  FIRST labs (earliest postop)...\n")
all_pts <- extract_labs(labs, all_pts, LAB_BASES, "first", descending=FALSE)
for(ln in LAB_BASES) { col<-paste0("first_",ln); nf<-sum(!is.na(all_pts[[col]]))
  cat(sprintf("    %s: %d (%.0f%%)\n",col,nf,100*nf/N)) }

# ── Indices + Cr + AKI ────────────────────────────────────────────
trt_idx<-which(all_pts$treated==1&!is.na(all_pts$mg_offset_h))
n_trt<-length(trt_idx); trt_pids<-all_pts$pid[trt_idx]; trt_tmg<-all_pts$mg_offset_h[trt_idx]

cat("  Cr lists + AKI...\n")
cr_list<-split(cr_all[,c("labresult","offset_h")],cr_all$pid)
ec<-sapply(cr_list,function(x)min(x$offset_h,na.rm=T))
all_pts$earliest_cr_h<-ec[as.character(all_pts$pid)]
all_pts$earliest_cr_h[is.na(all_pts$earliest_cr_h)]<-Inf
faki<-sapply(cr_list,function(cr){
  if(nrow(cr)<2)return(NA_real_); cr<-cr[order(cr$offset_h),]; bl<-cr$labresult[1]
  if(is.na(bl)||bl<=0)return(NA_real_)
  for(k in 2:nrow(cr)){d<-cr$labresult[k]-bl;r<-cr$labresult[k]/bl
    if(!is.na(d)&&(d>=0.3||(!is.na(r)&&r>=1.5)))return(cr$offset_h[k])}
  NA_real_})
all_pts$first_aki_h<-faki[as.character(all_pts$pid)]

# ── Risk sets (yet-untreated only) ────────────────────────────────
cat("  Risk sets (yet-untreated)...\n")
rs_yt<-vector("list",n_trt)
for(k in seq_len(n_trt)){tm<-trt_tmg[k]
  rs_yt[[k]]<-which(all_pts$icu_discharge_h>tm&(is.na(all_pts$first_aki_h)|all_pts$first_aki_h>tm)&
    all_pts$earliest_cr_h<=tm&(is.na(all_pts$mg_offset_h)|all_pts$mg_offset_h>tm+PRIMARY_H)&
    all_pts$pid!=trt_pids[k])
}
cat(sprintf("    YT median risk set: %.0f\n",median(sapply(rs_yt,length))))

# ── MICE on ALL candidate vars ────────────────────────────────────
all_cand <- unique(unlist(lapply(SPECS, `[[`, "vars")))
all_cand <- intersect(all_cand, names(all_pts))
all_cand <- all_cand[vapply(all_cand,function(v){x<-all_pts[[v]];!all(is.na(x))&&var(x,na.rm=T)>1e-10},logical(1))]
to_imp <- all_cand[vapply(all_cand,function(v)any(is.na(all_pts[[v]])),logical(1))]
set.seed(2026)  # reproducibility: MICE + matching
cat(sprintf("\n  MICE m=%d on %d vars: %s\n",M_IMP,length(to_imp),paste(to_imp,collapse=", ")))
if(length(to_imp)>0){
  md<-all_pts[,c("treated",all_cand)]; mt<-rep("",ncol(md)); names(mt)<-names(md)
  for(v in to_imp) mt[v]<-"pmm"
  imp<-mice(md,m=M_IMP,method=mt,printFlag=FALSE,maxit=10)
  cat(sprintf("  MICE done. Logged: %d\n",nrow(imp$loggedEvents)))
  for(v in to_imp){vals<-sapply(1:M_IMP,function(m)complete(imp,m)[[v]])
    all_pts[[v]]<-rowMeans(vals,na.rm=TRUE)}
}

# Caliper from primary PS
pv<-intersect(SPECS$primary$vars,names(all_pts))
pfml<-as.formula(paste("treated~",paste(pv,collapse="+")))
pfit<-suppressWarnings(glm(pfml,data=all_pts,family=binomial()))
caliper<-CALIPER_SD*sd(predict(pfit,type="response"),na.rm=TRUE)
cat(sprintf("  Caliper: %.4f\n",caliper))

# ═══════════════════════════════════════════════════════════════════
# RUN 3 SPECS × 1 POOL (yet-untreated only)
# ═══════════════════════════════════════════════════════════════════
all_res<-list(); all_pairs<-list()
for(sn in names(SPECS)){
  out<-run_spec_pool(sn,SPECS[[sn]],"yet_untreated",all_pts,trt_idx,rs_yt,cr_list,trt_tmg,caliper)
  if(!is.null(out)){all_res[[length(all_res)+1]]<-out$results
    all_pairs[[paste(sn,"yet_untreated",sep="_")]]<-out$pairs}
}

res_all<-do.call(rbind,all_res); res_all$db<-db
write.csv(res_all,file.path(RESULTS,sprintf("did_riskset_%s.csv",tag)),row.names=FALSE)
for(nm in names(all_pairs))
  write.csv(all_pairs[[nm]],file.path(RESULTS,sprintf("did_pairs_%s_%s.csv",nm,tag)),row.names=FALSE)

# ── Summary table ─────────────────────────────────────────────────
cat(sprintf("\n%s\n  SUMMARY: %s (PSM_DR, yet-untreated)\n%s\n",SEP,db,SEP))
cat("  spec               24h DiD       P      n\n")
for(sn in names(SPECS)){
  r<-res_all[res_all$spec==sn&res_all$pool=="yet_untreated"&res_all$method=="psm_dr"&res_all$target_h==PRIMARY_H,]
  if(nrow(r)==1&&!is.na(r$did))
    cat(sprintf("  %-16s  %+.4f  %.3f  %4d\n",sn,r$did,r$p,r$n))
  else cat(sprintf("  %-16s     --      --    --\n",sn))
}

cat(sprintf("\n%s\n02_psm.R -- %s DONE (%d rows)\n%s\n",SEP,db,nrow(res_all),SEP))

# ═══════════════════════════════════════════════════════════════════
# BINARY OUTCOMES: KDIGO staging + time-windowed mortality
# Runs on PRIMARY spec pairs only (the reported results)
# ═══════════════════════════════════════════════════════════════════
cat(sprintf("\n%s\n  BINARY OUTCOMES: %s (primary spec)\n%s\n", SEP, db, SEP))

pri_pairs <- all_pairs[["primary_yet_untreated"]]
if (!is.null(pri_pairs) && nrow(pri_pairs) > 0) {
  npairs <- nrow(pri_pairs)
  cat(sprintf("  Pairs: %d\n", npairs))

  # Pre-index patient info (NOT Cr — baseline must be computed per-pair from T0)
  rrt_map <- if ("rrt_offset_h" %in% names(all_pts))
    setNames(all_pts$rrt_offset_h, as.character(all_pts$pid)) else NULL
  death_map <- if ("death_offset_h" %in% names(all_pts))
    setNames(all_pts$death_offset_h, as.character(all_pts$pid)) else NULL
  mort_map <- setNames(all_pts$hosp_mortality, as.character(all_pts$pid))

  # For each pair, compute outcomes for both arms
  outcome_rows <- vector("list", npairs)
  for (k in seq_len(npairs)) {
    tp <- as.character(pri_pairs$trt_pid[k])
    cp <- as.character(pri_pairs$ctl_pid[k])
    t0 <- pri_pairs$t_mg[k]

    # Correct baseline: last Cr before T0 (same as DiD and 03_hte.R)
    bl_t <- find_cr_pre(cr_list[[tp]], t0)[1]
    bl_c <- find_cr_pre(cr_list[[cp]], t0)[1]

    cr_t <- cr_list[[tp]]; cr_c <- cr_list[[cp]]
    cr_t48 <- if (!is.null(cr_t)) cr_t[cr_t$offset_h > t0 & cr_t$offset_h <= t0 + 48, ] else data.frame()
    cr_c48 <- if (!is.null(cr_c)) cr_c[cr_c$offset_h > t0 & cr_c$offset_h <= t0 + 48, ] else data.frame()
    cr_t7d <- if (!is.null(cr_t)) cr_t[cr_t$offset_h > t0 & cr_t$offset_h <= t0 + 168, ] else data.frame()
    cr_c7d <- if (!is.null(cr_c)) cr_c[cr_c$offset_h > t0 & cr_c$offset_h <= t0 + 168, ] else data.frame()

    max_t48 <- if (nrow(cr_t48) > 0) max(cr_t48$labresult) else NA
    max_c48 <- if (nrow(cr_c48) > 0) max(cr_c48$labresult) else NA
    max_t7d <- if (nrow(cr_t7d) > 0) max(cr_t7d$labresult) else NA
    max_c7d <- if (nrow(cr_c7d) > 0) max(cr_c7d$labresult) else NA

    # RRT within windows
    rrt_t <- if (!is.null(rrt_map)) rrt_map[tp] else NA
    rrt_c <- if (!is.null(rrt_map)) rrt_map[cp] else NA
    rrt_t48 <- !is.na(rrt_t) && rrt_t > t0 && rrt_t <= t0 + 48
    rrt_c48 <- !is.na(rrt_c) && rrt_c > t0 && rrt_c <= t0 + 48
    rrt_t7d <- !is.na(rrt_t) && rrt_t > t0 && rrt_t <= t0 + 168
    rrt_c7d <- !is.na(rrt_c) && rrt_c > t0 && rrt_c <= t0 + 168

    # Death within windows
    d_t <- if (!is.null(death_map)) death_map[tp] else NA
    d_c <- if (!is.null(death_map)) death_map[cp] else NA
    d_t_from_t0 <- if (!is.na(d_t)) d_t - t0 else NA
    d_c_from_t0 <- if (!is.na(d_c)) d_c - t0 else NA

    # Helper: KDIGO staging
    kdigo <- function(max_cr, bl, rrt_in_window) {
      if (is.na(max_cr) || is.na(bl) || bl <= 0) return(c(NA, NA, NA))
      delta <- max_cr - bl; ratio <- max_cr / bl
      stg1 <- as.integer(delta >= 0.3 | ratio >= 1.5 | rrt_in_window)
      stg2 <- as.integer(ratio >= 2.0 | rrt_in_window)
      stg3 <- as.integer(ratio >= 3.0 | max_cr >= 4.0 | rrt_in_window)
      c(stg1, stg2, stg3)
    }

    kt48 <- kdigo(max_t48, bl_t, rrt_t48)
    kc48 <- kdigo(max_c48, bl_c, rrt_c48)
    kt7d <- kdigo(max_t7d, bl_t, rrt_t7d)
    kc7d <- kdigo(max_c7d, bl_c, rrt_c7d)

    outcome_rows[[k]] <- data.frame(
      trt_pid = tp, ctl_pid = cp, t0 = t0,
      # KDIGO 48h
      aki1_48h_trt = kt48[1], aki1_48h_ctl = kc48[1],
      aki2_48h_trt = kt48[2], aki2_48h_ctl = kc48[2],
      aki3_48h_trt = kt48[3], aki3_48h_ctl = kc48[3],
      # KDIGO 7d
      aki1_7d_trt = kt7d[1], aki1_7d_ctl = kc7d[1],
      aki2_7d_trt = kt7d[2], aki2_7d_ctl = kc7d[2],
      aki3_7d_trt = kt7d[3], aki3_7d_ctl = kc7d[3],
      # Mortality
      hosp_mort_trt = as.integer(mort_map[tp] == 1),
      hosp_mort_ctl = as.integer(mort_map[cp] == 1),
      death_7d_trt  = as.integer(!is.na(d_t_from_t0) && d_t_from_t0 >= 0 && d_t_from_t0 <= 168),
      death_7d_ctl  = as.integer(!is.na(d_c_from_t0) && d_c_from_t0 >= 0 && d_c_from_t0 <= 168),
      death_14d_trt = as.integer(!is.na(d_t_from_t0) && d_t_from_t0 >= 0 && d_t_from_t0 <= 336),
      death_14d_ctl = as.integer(!is.na(d_c_from_t0) && d_c_from_t0 >= 0 && d_c_from_t0 <= 336),
      stringsAsFactors = FALSE
    )
  }
  bo <- do.call(rbind, outcome_rows)

  # GLM for each binary outcome
  outcome_names <- c("aki1_48h", "aki2_48h", "aki3_48h",
                     "aki1_7d",  "aki2_7d",  "aki3_7d",
                     "hosp_mort", "death_7d", "death_14d")
  binary_res <- list()
  for (oname in outcome_names) {
    y_t <- bo[[paste0(oname, "_trt")]]
    y_c <- bo[[paste0(oname, "_ctl")]]
    valid <- !is.na(y_t) & !is.na(y_c)
    nv <- sum(valid)
    if (nv < 30 || sum(y_t[valid]) + sum(y_c[valid]) < 5) {
      binary_res[[length(binary_res) + 1]] <- data.frame(
        db = db, outcome = oname, n = nv,
        rate_trt = NA, rate_ctl = NA, or = NA, ci_lo = NA, ci_hi = NA, p = NA)
      next
    }
    y <- c(y_t[valid], y_c[valid])
    x <- rep(c(1L, 0L), each = nv)
    fit <- tryCatch(
      glm(y ~ x, family = quasibinomial()),
      error = function(e) NULL)
    if (is.null(fit)) {
      binary_res[[length(binary_res) + 1]] <- data.frame(
        db = db, outcome = oname, n = nv,
        rate_trt = mean(y_t[valid]), rate_ctl = mean(y_c[valid]),
        or = NA, ci_lo = NA, ci_hi = NA, p = NA)
      next
    }
    ct <- safe_coeftest(fit)
    log_or <- ct["x", 1]; se <- ct["x", 2]; pv <- ct["x", 4]
    binary_res[[length(binary_res) + 1]] <- data.frame(
      db = db, outcome = oname, n = nv,
      rate_trt = mean(y_t[valid]), rate_ctl = mean(y_c[valid]),
      or = exp(log_or),
      ci_lo = exp(log_or - 1.96 * se),
      ci_hi = exp(log_or + 1.96 * se),
      p = pv)
  }
  br <- do.call(rbind, binary_res)

  # Print summary
  cat("\n  outcome         rate_trt  rate_ctl  OR (95% CI)            P\n")
  cat("  ", paste(rep("-", 66), collapse=""), "\n", sep="")
  for (i in seq_len(nrow(br))) {
    r <- br[i, ]
    if (is.na(r$or)) {
      cat(sprintf("  %-16s  %5.1f%%    %5.1f%%    --                    --  (n=%d)\n",
                  r$outcome, 100*r$rate_trt, 100*r$rate_ctl, r$n))
    } else {
      sig <- if (!is.na(r$p) && r$p < 0.05) " *" else "  "
      cat(sprintf("  %-16s  %5.1f%%    %5.1f%%    %.2f (%.2f-%.2f)  %.4f%s  n=%d\n",
                  r$outcome, 100*r$rate_trt, 100*r$rate_ctl,
                  r$or, r$ci_lo, r$ci_hi, r$p, sig, r$n))
    }
  }

  # Save
  write.csv(br, file.path(RESULTS, sprintf("did_binary_%s.csv", tag)), row.names = FALSE)
  write.csv(bo, file.path(RESULTS, sprintf("did_binary_pairs_%s.csv", tag)), row.names = FALSE)
  cat(sprintf("\n  ✓ did_binary_%s.csv: %d outcomes\n", tag, nrow(br)))
  cat(sprintf("  ✓ did_binary_pairs_%s.csv: %d pairs × %d outcomes\n", tag, nrow(bo), length(outcome_names)))
} else {
  cat("  WARN: No primary pairs found, skipping binary outcomes\n")
}
cat(sprintf("\n%s\n02_psm.R -- %s COMPLETE\n%s\n", SEP, db, SEP))
