#!/usr/bin/env Rscript
# ============================================================================
# 03_hte.R — HTE with post-T₀ secondary outcomes (v5.1)
#
# AKI definition: Consolidated KDIGO (2012)
#   ≤48h: ΔCr ≥ 0.3 OR ratio ≥ 1.5 (both criteria active)
#   >48h: ratio ≥ 1.5 only (absolute criterion expired)
#   UO criterion not used (unreliable in cardiac surgery ICU)
#
# Outcomes (all T₀-anchored):
#   ΔCr 48h:        Cr at T₀+48h − Cr_pre
#   48h AKI:        consolidated KDIGO within 48h of T₀
#   7d AKI:         consolidated KDIGO within 7d of T₀
#   Mortality:      per-stay (patient alive at T₀ by construction)
#   POAF/Enceph/VA: eICU: post-T₀ from raw dx table (7d window)
#                   MIMIC: hospitalization-level ICD (no T₀ filter)
#
# Usage: Rscript 03_hte.R mimic
#        Rscript 03_hte.R eicu
# ============================================================================

suppressPackageStartupMessages({ library(sandwich); library(lmtest) })

RESULTS   <- path.expand("~/mg_aki/results")
EICU_ROOT <- path.expand("~/mg_aki/eicu-crd-2.0")
CR_WINDOW <- 12
WINDOW_7D_MIN <- 7 * 24 * 60   # 7 days in minutes
LAB_BASES <- c("magnesium","potassium","calcium","lactate","heartrate")

# Secondary outcome patterns (eICU diagnosis text matching)
AF_PAT  <- c("atrial fibrillation","atrial flutter","a-fib","afib","new onset af")
AF_PRIOR<- c("atrial fibrillation","atrial flutter","a-fib","afib","chronic af")
ENC_PAT <- c("encephalopathy","delirium","altered mental","acute confusional","metabolic encephalopathy")
VA_PAT  <- c("ventricular tachycardia","ventricular fibrillation","v-tach","v-fib","vtach","vfib","cardiac arrest")

# ── Subgroups ──
# Both sides of each comparison are needed: the HTE IS the contrast
# between subgroup and complement (e.g., CKD vs No CKD).
SUBGROUPS <- list(
  list(name="Overall",var=NULL),
  list(name="Age < 65",var="age",op="<",val=65),
  list(name="Age >= 65",var="age",op=">=",val=65),
  list(name="eGFR < 60",var="egfr",op="<",val=60),
  list(name="eGFR >= 60",var="egfr",op=">=",val=60),
  list(name="Mg < 1.8",var="last_magnesium",op="<",val=1.8),
  list(name="Mg >= 1.8",var="last_magnesium",op=">=",val=1.8),
  list(name="CABG",var="surg_cabg",op="==",val=1),
  list(name="Non-CABG",var="surg_cabg",op="==",val=0),
  list(name="Diabetes",var="diabetes",op="==",val=1),
  list(name="No diabetes",var="diabetes",op="==",val=0),
  list(name="CKD",var="ckd",op="==",val=1),
  list(name="No CKD",var="ckd",op="==",val=0),
  list(name="Heart failure",var="heart_failure",op="==",val=1),
  list(name="No HF",var="heart_failure",op="==",val=0),
  list(name="BMI >= 30",var="bmi",op=">=",val=30),
  list(name="BMI < 30",var="bmi",op="<",val=30),
  list(name="DM + CKD",var=c("diabetes","ckd"),op=c("==","=="),val=c(1,1)),
  list(name="HF + CABG",var=c("heart_failure","surg_cabg"),op=c("==","=="),val=c(1,1)),
  list(name="Mg<1.8 + CKD",var=c("last_magnesium","ckd"),op=c("<","=="),val=c(1.8,1))
)

# ═══════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════
find_cr <- function(cp,th,w=CR_WINDOW){if(is.null(cp)||nrow(cp)==0)return(NA)
  cd<-cp[cp$offset_h>=(th-w)&cp$offset_h<=(th+w),];if(nrow(cd)==0)return(NA)
  cd$labresult[which.min(abs(cd$offset_h-th))]}
find_cr_pre <- function(cp,th){if(is.null(cp)||nrow(cp)==0)return(NA)
  cd<-cp[cp$offset_h>=0&cp$offset_h<th,];if(nrow(cd)==0)return(NA)
  cd$labresult[which.max(cd$offset_h)]}
max_cr_win <- function(cp,t0,hrs){if(is.null(cp)||nrow(cp)==0)return(NA)
  cd<-cp[cp$offset_h>=t0&cp$offset_h<=(t0+hrs),];if(nrow(cd)==0)return(NA)
  max(cd$labresult,na.rm=TRUE)}

safe_ct <- function(fit){
  ct<-tryCatch(suppressWarnings(coeftest(fit,vcov.=vcovHC(fit,type="HC1"))),error=function(e)NULL)
  if(!is.null(ct)&&is.matrix(ct)&&"treated"%in%rownames(ct)&&ncol(ct)>=4&&!any(is.nan(ct["treated",]))) ct
  else tryCatch(coeftest(fit),error=function(e)NULL)}

matches_any <- function(s, pats) {
  s <- tolower(as.character(s)); m <- rep(FALSE, length(s))
  for (p in pats) m <- m | grepl(tolower(p), s, fixed=TRUE)
  m
}

test_sg <- function(df, sg, oc, otype) {
  if(!is.null(sg$var)){mask<-rep(TRUE,nrow(df))
    for(i in seq_along(sg$var)){v<-sg$var[i];o<-sg$op[i];val<-sg$val[i]
      if(!v%in%names(df))return(NULL); x<-df[[v]]
      mask<-mask&!is.na(x)&switch(o,"<"=x<val,">="=x>=val,"=="=x==val)}
    sub<-df[mask,]} else sub<-df
  nt<-sum(sub$treated==1);nc<-sum(sub$treated==0)
  if(nt<15||nc<15) return(data.frame(subgroup=sg$name,outcome=oc,type=otype,
    n_trt=nt,n_ctl=nc,est=NA,se=NA,p=NA,rate_trt=NA,rate_ctl=NA,rd=NA,nnt=NA,stringsAsFactors=FALSE))
  if(otype=="continuous"){sub<-sub[!is.na(sub[[oc]]),];if(nrow(sub)<30)return(NULL)
    fit<-lm(as.formula(paste(oc,"~treated")),data=sub);ct<-safe_ct(fit);if(is.null(ct))return(NULL)
    data.frame(subgroup=sg$name,outcome=oc,type="continuous",n_trt=sum(sub$treated==1),
      n_ctl=sum(sub$treated==0),est=ct["treated",1],se=ct["treated",2],p=ct["treated",4],
      rate_trt=NA,rate_ctl=NA,rd=NA,nnt=NA,stringsAsFactors=FALSE)
  } else {sub<-sub[!is.na(sub[[oc]]),]
    r1<-mean(sub[[oc]][sub$treated==1]);r0<-mean(sub[[oc]][sub$treated==0])
    rd<-r1-r0;nnt<-if(abs(rd)>0.001)round(1/abs(rd))else NA
    fit<-tryCatch(glm(as.formula(paste(oc,"~treated")),data=sub,family=quasibinomial()),error=function(e)NULL)
    if(is.null(fit)||!"treated"%in%names(coef(fit)))return(NULL)
    ct<-safe_ct(fit);if(is.null(ct))return(NULL)
    data.frame(subgroup=sg$name,outcome=oc,type="binary",n_trt=sum(sub$treated==1),
      n_ctl=sum(sub$treated==0),est=exp(ct["treated",1]),se=ct["treated",2],p=ct["treated",4],
      rate_trt=r1,rate_ctl=r0,rd=rd,nnt=nnt,stringsAsFactors=FALSE)}
}

test_int <- function(df,iv,oc){if(!iv%in%names(df)||!oc%in%names(df))return(NA)
  d<-df[!is.na(df[[iv]])&!is.na(df[[oc]]),];if(nrow(d)<60)return(NA)
  fit<-tryCatch(lm(as.formula(sprintf("%s~treated*%s",oc,iv)),data=d),error=function(e)NULL)
  if(is.null(fit))return(NA)
  ct<-tryCatch(coeftest(fit,vcov.=vcovHC(fit,type="HC1")),error=function(e)NULL)
  if(is.null(ct))return(NA);ir<-grep(paste0("treated:",iv),rownames(ct))
  if(length(ir)==0)NA else ct[ir[1],4]}

# ═══════════════════════════════════════════════════════════════════
# eICU: post-T₀ secondary outcomes from raw diagnosis table
# ═══════════════════════════════════════════════════════════════════
extract_eicu_post_t0 <- function(pairs, all_pts) {
  cat("  Loading eICU raw diagnosis table for post-T₀ outcomes...\n")
  dx_path <- file.path(EICU_ROOT, "diagnosis.csv.gz")
  if (!file.exists(dx_path)) dx_path <- file.path(EICU_ROOT, "diagnosis.csv")
  if (!file.exists(dx_path)) { cat("    WARN: diagnosis table not found\n"); return(NULL) }

  dx <- read.csv(dx_path, stringsAsFactors=FALSE)
  dx$diagnosisstring <- tolower(dx$diagnosisstring)
  all_pids <- unique(c(pairs$trt_pid, pairs$ctl_pid))
  dx <- dx[dx$patientunitstayid %in% all_pids, ]
  cat(sprintf("    Loaded %d diagnosis rows for %d patients\n", nrow(dx), length(all_pids)))

  if ("prior_af" %in% names(all_pts)) {
    prior_af_pids <- all_pts$pid[all_pts$prior_af == 1]
  } else {
    prior_af_pids <- character(0)
    ph_path <- file.path(EICU_ROOT, "pastHistory.csv.gz")
    if (!file.exists(ph_path)) ph_path <- file.path(EICU_ROOT, "pastHistory.csv")
    if (file.exists(ph_path)) {
      cat("    Loading pastHistory for prior AF...\n")
      ph <- read.csv(ph_path, stringsAsFactors=FALSE)
      ph <- ph[ph$patientunitstayid %in% all_pids, ]
      if ("pasthistorypath" %in% names(ph)) {
        prior_af_pids <- unique(ph$patientunitstayid[matches_any(ph$pasthistorypath, AF_PRIOR)])
        cat(sprintf("    Prior AF: %d patients\n", length(prior_af_pids)))
      }
    }
  }

  disch_map <- setNames(all_pts$icu_discharge_h * 60, all_pts$pid)
  n <- nrow(pairs)
  poaf_trt <- poaf_ctl <- enc_trt <- enc_ctl <- va_trt <- va_ctl <- integer(n)

  for (i in seq_len(n)) {
    tp <- pairs$trt_pid[i]; cp <- pairs$ctl_pid[i]
    t0_min <- pairs$t_mg[i] * 60
    end_trt <- min(t0_min + WINDOW_7D_MIN, disch_map[as.character(tp)], na.rm=TRUE)
    end_ctl <- min(t0_min + WINDOW_7D_MIN, disch_map[as.character(cp)], na.rm=TRUE)

    for (info in list(list(pid=tp,end=end_trt,side="trt"),
                      list(pid=cp,end=end_ctl,side="ctl"))) {
      pdx <- dx[dx$patientunitstayid == info$pid &
                dx$diagnosisoffset > t0_min &
                dx$diagnosisoffset <= info$end, ]
      if (nrow(pdx) == 0) next
      ds <- pdx$diagnosisstring
      if (!info$pid %in% prior_af_pids && any(matches_any(ds, AF_PAT))) {
        if (info$side=="trt") poaf_trt[i]<-1L else poaf_ctl[i]<-1L
      }
      if (any(matches_any(ds, ENC_PAT))) {
        if (info$side=="trt") enc_trt[i]<-1L else enc_ctl[i]<-1L
      }
      if (any(matches_any(ds, VA_PAT))) {
        if (info$side=="trt") va_trt[i]<-1L else va_ctl[i]<-1L
      }
    }
  }

  data.frame(pair_id = seq_len(n),
    poaf_trt=poaf_trt, poaf_ctl=poaf_ctl,
    enc_trt=enc_trt, enc_ctl=enc_ctl,
    va_trt=va_trt, va_ctl=va_ctl)
}

# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════
args<-commandArgs(trailingOnly=TRUE)
if(length(args)<1){cat("Usage: Rscript 03_hte.R <db>\n");quit(status=1)}
db<-toupper(args[1]); tag<-tolower(db); is_eicu <- tag == "eicu"

SEP<-paste(rep("=",70),collapse="")
cat(sprintf("\n%s\n03_hte.R — HTE: %s (post-T0 outcomes)\n%s\n",SEP,db,SEP))

pairs<-read.csv(file.path(RESULTS,sprintf("did_pairs_primary_yet_untreated_%s.csv",tag)),stringsAsFactors=FALSE)
all_pts<-read.csv(file.path(RESULTS,sprintf("did_all_%s.csv",tag)),stringsAsFactors=FALSE)
cr_all<-read.csv(file.path(RESULTS,sprintf("did_cr_all_%s.csv",tag)),stringsAsFactors=FALSE)
cr_id<-if("patientunitstayid"%in%names(cr_all))"patientunitstayid" else "stay_id"
cr_all$pid<-cr_all[[cr_id]]
if(!"offset_h"%in%names(cr_all)) cr_all$offset_h<-cr_all$labresultoffset/60
cr_all<-cr_all[order(cr_all$pid,cr_all$offset_h),]
cr_list<-split(cr_all[,c("labresult","offset_h")],cr_all$pid)
cat(sprintf("  Pairs: %d | Patients: %d\n",nrow(pairs),nrow(all_pts)))

# LAST labs
cat("  LAST labs...\n")
labs_raw<-read.csv(file.path(RESULTS,sprintf("did_labs_all_%s.csv",tag)),stringsAsFactors=FALSE)
pcl<-if("patientunitstayid"%in%names(labs_raw))"patientunitstayid" else "stay_id"
if(pcl%in%names(labs_raw)) labs_raw$pid<-labs_raw[[pcl]]
le<-labs_raw[labs_raw$lab_name%in%c("magnesium","potassium","calcium","lactate"),]
lh<-labs_raw[labs_raw$lab_name=="heartrate",]
if(nrow(lh)>500000){lh$hb<-floor(lh$offset_h);lh<-lh[order(lh$pid,lh$offset_h),]
  lh<-lh[!duplicated(paste(lh$pid,lh$hb)),];lh$hb<-NULL}
labs<-rbind(le,lh); rm(labs_raw,le,lh)
for(ln in LAB_BASES){sub<-labs[labs$lab_name==ln,];if(nrow(sub)==0)next
  sub$mg_offset_h<-all_pts$mg_offset_h[match(sub$pid,all_pts$pid)]
  sub<-sub[sub$offset_h>=0&(is.na(sub$mg_offset_h)|sub$offset_h<sub$mg_offset_h),]
  if(nrow(sub)==0)next;s<-sub[order(-sub$offset_h),];s<-s[!duplicated(s$pid),]
  all_pts[[paste0("last_",ln)]]<-s$value[match(all_pts$pid,s$pid)]}

# eICU: post-T₀ secondary outcomes
eicu_sec <- NULL
if (is_eicu) {
  eicu_sec <- extract_eicu_post_t0(pairs, all_pts)
  if (!is.null(eicu_sec))
    cat(sprintf("    Post-T0 POAF: %d/%d trt, %d/%d ctl\n",
                sum(eicu_sec$poaf_trt), nrow(eicu_sec),
                sum(eicu_sec$poaf_ctl), nrow(eicu_sec)))
}

# Build outcome dataset
cat("  Building outcomes...\n")
rows<-vector("list",nrow(pairs)*2); ri<-0
for(i in seq_len(nrow(pairs))){
  tp<-as.character(pairs$trt_pid[i]);cp<-as.character(pairs$ctl_pid[i]);tmg<-pairs$t_mg[i]
  cr_pre_t<-find_cr_pre(cr_list[[tp]],tmg);cr_pre_c<-find_cr_pre(cr_list[[cp]],tmg)
  dcr48_t<-{v<-find_cr(cr_list[[tp]],tmg+48);if(!is.na(v)&&!is.na(cr_pre_t))v-cr_pre_t else NA}
  dcr48_c<-{v<-find_cr(cr_list[[cp]],tmg+48);if(!is.na(v)&&!is.na(cr_pre_c))v-cr_pre_c else NA}
  m48t<-max_cr_win(cr_list[[tp]],tmg,48);m48c<-max_cr_win(cr_list[[cp]],tmg,48)
  m7t<-max_cr_win(cr_list[[tp]],tmg,168);m7c<-max_cr_win(cr_list[[cp]],tmg,168)
  aki48_t<-as.integer(!is.na(m48t)&&!is.na(cr_pre_t)&&
    ((m48t-cr_pre_t)>=0.3||(cr_pre_t>0&&m48t/cr_pre_t>=1.5)))
  aki48_c<-as.integer(!is.na(m48c)&&!is.na(cr_pre_c)&&
    ((m48c-cr_pre_c)>=0.3||(cr_pre_c>0&&m48c/cr_pre_c>=1.5)))
  aki7_t<-as.integer(!is.na(cr_pre_t)&&cr_pre_t>0&&(
    (!is.na(m48t)&&(m48t-cr_pre_t)>=0.3)||(!is.na(m7t)&&m7t/cr_pre_t>=1.5)))
  aki7_c<-as.integer(!is.na(cr_pre_c)&&cr_pre_c>0&&(
    (!is.na(m48c)&&(m48c-cr_pre_c)>=0.3)||(!is.na(m7c)&&m7c/cr_pre_c>=1.5)))
  ti<-which(all_pts$pid==pairs$trt_pid[i])[1];ci<-which(all_pts$pid==pairs$ctl_pid[i])[1]
  if(is.na(ti)||is.na(ci))next

  if (!is.null(eicu_sec)) {
    poaf_t<-eicu_sec$poaf_trt[i]; poaf_c<-eicu_sec$poaf_ctl[i]
    enc_t<-eicu_sec$enc_trt[i];   enc_c<-eicu_sec$enc_ctl[i]
    va_t<-eicu_sec$va_trt[i];     va_c<-eicu_sec$va_ctl[i]
  } else {
    poaf_t<-all_pts$poaf[ti];     poaf_c<-all_pts$poaf[ci]
    enc_t<-all_pts$encephalopathy[ti]; enc_c<-all_pts$encephalopathy[ci]
    va_t<-all_pts$vent_arrhythmia[ti]; va_c<-all_pts$vent_arrhythmia[ci]
  }

  covs<-c("age","is_female","bmi","egfr","surg_cabg","surg_valve","diabetes","ckd",
          "heart_failure","last_magnesium")
  mk<-function(idx,trt,dcr,a48,a7,crp,pf,enc,va){
    r<-data.frame(pair_id=i,treated=trt,pid=all_pts$pid[idx],t_mg=tmg,
      dcr_48h=dcr,aki_48h=a48,aki_7d=a7,cr_pre=crp,
      hosp_mortality=all_pts$hosp_mortality[idx],
      poaf=pf,encephalopathy=enc,vent_arrhythmia=va,stringsAsFactors=FALSE)
    for(cv in covs) r[[cv]]<-if(cv%in%names(all_pts)) all_pts[[cv]][idx] else NA
    r}
  ri<-ri+1;rows[[ri]]<-mk(ti,1,dcr48_t,aki48_t,aki7_t,cr_pre_t,poaf_t,enc_t,va_t)
  ri<-ri+1;rows[[ri]]<-mk(ci,0,dcr48_c,aki48_c,aki7_c,cr_pre_c,poaf_c,enc_c,va_c)
}
df<-do.call(rbind,rows[1:ri])
cat(sprintf("  Dataset: %d rows (%d pairs)\n",nrow(df),nrow(df)/2))

if (is_eicu && !is.null(eicu_sec))
  cat("  NOTE: POAF/encephalopathy/vent_arrhythmia use post-T0 extraction (7d window)\n")
if (is_eicu && is.null(eicu_sec))
  cat("  WARN: eICU post-T0 extraction FAILED - using ICU-stay flags (BIASED)\n")
if (!is_eicu)
  cat("  NOTE: POAF/encephalopathy/vent_arrhythmia are hospitalization-level ICD (no T0 filter)\n")

# Overall rates
cat("\n  Overall:\n")
for(oc in c("dcr_48h","aki_48h","aki_7d","hosp_mortality","poaf","encephalopathy","vent_arrhythmia")){
  if(!oc%in%names(df))next;r1<-mean(df[[oc]][df$treated==1],na.rm=T);r0<-mean(df[[oc]][df$treated==0],na.rm=T)
  if(oc=="dcr_48h") cat(sprintf("    %s: trt=%+.4f ctl=%+.4f DiD=%+.4f\n",oc,r1,r0,r1-r0))
  else cat(sprintf("    %s: trt=%.1f%% ctl=%.1f%% RD=%+.1f%%\n",oc,100*r1,100*r0,100*(r1-r0)))}

# ═══════════════════════════════════════════════════════════════════
# SWEEP
# ═══════════════════════════════════════════════════════════════════
OC<-list(list(c="dcr_48h",t="continuous",l="dCr 48h"),
  list(c="aki_48h",t="binary",l="48h AKI (KDIGO, both criteria)"),
  list(c="aki_7d",t="binary",l="7d AKI (KDIGO consolidated)"),
  list(c="hosp_mortality",t="binary",l="Mortality"),list(c="poaf",t="binary",l="POAF"),
  list(c="encephalopathy",t="binary",l="Encephalopathy"),
  list(c="vent_arrhythmia",t="binary",l="Vent arrhythmia"))

hte_rows<-list()
for(sg in SUBGROUPS) for(oc in OC){r<-test_sg(df,sg,oc$c,oc$t)
  if(!is.null(r)) hte_rows[[length(hte_rows)+1]]<-r}
hte<-do.call(rbind,hte_rows)

for(oc in OC){cat(sprintf("\n  == %s ==\n",oc$l));sub<-hte[hte$outcome==oc$c,]
  if(oc$t=="continuous"){
    cat(sprintf("  %-20s %5s %5s %9s %8s\n","Subgroup","nT","nC","DiD","P"))
    for(i in 1:nrow(sub)){r<-sub[i,];sig<-if(!is.na(r$p)&&r$p<0.05)" *" else "  "
      cat(sprintf("  %-20s %5d %5d %+.4f  %.4f%s\n",r$subgroup,r$n_trt,r$n_ctl,
                  ifelse(is.na(r$est),0,r$est),ifelse(is.na(r$p),1,r$p),sig))}
  } else {
    cat(sprintf("  %-20s %5s %5s %5s %5s %5s %5s %8s\n","Subgroup","nT","nC","T%","C%","RD%","OR","P"))
    for(i in 1:nrow(sub)){r<-sub[i,];sig<-if(!is.na(r$p)&&r$p<0.05)" *" else "  "
      cat(sprintf("  %-20s %5d %5d %5.1f %5.1f %+4.1f %5.2f  %.4f%s\n",r$subgroup,r$n_trt,r$n_ctl,
                  100*ifelse(is.na(r$rate_trt),0,r$rate_trt),100*ifelse(is.na(r$rate_ctl),0,r$rate_ctl),
                  100*ifelse(is.na(r$rd),0,r$rd),ifelse(is.na(r$est),1,r$est),
                  ifelse(is.na(r$p),1,r$p),sig))}
  }}

# Interactions
cat(sprintf("\n  == INTERACTIONS (dCr 48h) ==\n"))
df$age_ge65<-as.integer(df$age>=65);df$egfr_lt60<-as.integer(df$egfr<60)
df$mg_lt18<-as.integer(!is.na(df$last_magnesium)&df$last_magnesium<1.8)
df$bmi_ge30<-as.integer(!is.na(df$bmi)&df$bmi>=30)
for(iv in c("age_ge65","egfr_lt60","mg_lt18","surg_cabg","diabetes","ckd","heart_failure","bmi_ge30")){
  p<-test_int(df,iv,"dcr_48h");sig<-if(!is.na(p)&&p<0.05)" *" else "  "
  cat(sprintf("    treated x %-16s P=%.4f%s\n",iv,ifelse(is.na(p),NA,p),sig))}

write.csv(hte,file.path(RESULTS,sprintf("did_hte_%s.csv",tag)),row.names=FALSE)
write.csv(df,file.path(RESULTS,sprintf("did_hte_data_%s.csv",tag)),row.names=FALSE)
cat(sprintf("\n%s\n03_hte.R — %s DONE (%d rows)\n%s\n",SEP,db,nrow(hte),SEP))
