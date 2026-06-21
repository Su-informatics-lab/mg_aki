#!/usr/bin/env Rscript
# ============================================================================
# 03_hte.R — Heterogeneous Treatment Effects (v5)
#
# Input:  primary matched pairs (19var LAST, no K/Mg, yet-untreated)
# Outcomes:
#   Continuous: ΔCr at 48h (wrt T₀)
#   AKI binary: 48h AKI (ΔCr≥0.3 wrt T₀), 7d AKI (ratio≥1.5 wrt T₀)
#   Clinical:   mortality, POAF, encephalopathy, ventricular arrhythmia
# Subgroups: eGFR, Mg, surgery, DM, CKD, HF, BMI, age + crossed
#
# Usage: Rscript 03_hte.R mimic
#        Rscript 03_hte.R eicu
# ============================================================================

suppressPackageStartupMessages({ library(sandwich); library(lmtest) })

RESULTS   <- path.expand("~/mg_aki/results")
CR_WINDOW <- 12
LAB_BASES <- c("magnesium","potassium","calcium","lactate","heartrate")

SUBGROUPS <- list(
  list(name="Overall",         var=NULL),
  list(name="Age < 65",        var="age",            op="<",  val=65),
  list(name="Age >= 65",       var="age",            op=">=", val=65),
  list(name="eGFR < 60",       var="egfr",           op="<",  val=60),
  list(name="eGFR >= 60",      var="egfr",           op=">=", val=60),
  list(name="Mg < 1.8",        var="last_magnesium", op="<",  val=1.8),
  list(name="Mg >= 1.8",       var="last_magnesium", op=">=", val=1.8),
  list(name="CABG",            var="surg_cabg",      op="==", val=1),
  list(name="Non-CABG",        var="surg_cabg",      op="==", val=0),
  list(name="Diabetes",        var="diabetes",       op="==", val=1),
  list(name="No diabetes",     var="diabetes",       op="==", val=0),
  list(name="CKD",             var="ckd",            op="==", val=1),
  list(name="No CKD",          var="ckd",            op="==", val=0),
  list(name="Heart failure",   var="heart_failure",  op="==", val=1),
  list(name="No HF",           var="heart_failure",  op="==", val=0),
  list(name="BMI >= 30",       var="bmi",            op=">=", val=30),
  list(name="BMI < 30",        var="bmi",            op="<",  val=30),
  list(name="DM + CKD",        var=c("diabetes","ckd"), op=c("==","=="), val=c(1,1)),
  list(name="HF + CABG",       var=c("heart_failure","surg_cabg"), op=c("==","=="), val=c(1,1)),
  list(name="Mg<1.8 + CKD",    var=c("last_magnesium","ckd"), op=c("<","=="), val=c(1.8,1))
)

# ═══════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════
find_cr <- function(cr_pt, target_h, window=CR_WINDOW) {
  if(is.null(cr_pt)||nrow(cr_pt)==0) return(NA)
  cand<-cr_pt[cr_pt$offset_h>=(target_h-window)&cr_pt$offset_h<=(target_h+window),]
  if(nrow(cand)==0) return(NA)
  cand$labresult[which.min(abs(cand$offset_h-target_h))]
}
find_cr_pre <- function(cr_pt, t_h) {
  if(is.null(cr_pt)||nrow(cr_pt)==0) return(NA)
  cand<-cr_pt[cr_pt$offset_h>=0&cr_pt$offset_h<t_h,]
  if(nrow(cand)==0) return(NA)
  cand$labresult[which.max(cand$offset_h)]
}
max_cr_window <- function(cr_pt, t0, hours) {
  if(is.null(cr_pt)||nrow(cr_pt)==0) return(NA)
  cand<-cr_pt[cr_pt$offset_h>=t0&cr_pt$offset_h<=(t0+hours),]
  if(nrow(cand)==0) return(NA)
  max(cand$labresult, na.rm=TRUE)
}

safe_ct <- function(fit) {
  ct<-tryCatch(suppressWarnings(coeftest(fit,vcov.=vcovHC(fit,type="HC1"))),error=function(e)NULL)
  if(!is.null(ct)&&is.matrix(ct)&&"treated"%in%rownames(ct)&&ncol(ct)>=4&&!any(is.nan(ct["treated",]))) ct
  else tryCatch(coeftest(fit),error=function(e)NULL)
}

test_subgroup <- function(df, sg, outcome_col, outcome_type) {
  if (!is.null(sg$var)) {
    mask <- rep(TRUE, nrow(df))
    for (i in seq_along(sg$var)) {
      v<-sg$var[i]; o<-sg$op[i]; val<-sg$val[i]
      if (!v %in% names(df)) return(NULL)
      x<-df[[v]]
      mask<-mask & !is.na(x) & switch(o,"<"=x<val,">="=x>=val,"=="=x==val,">"=x>val,"<="=x<=val)
    }
    sub<-df[mask,]
  } else sub<-df

  nt<-sum(sub$treated==1); nc<-sum(sub$treated==0)
  if(nt<15||nc<15) return(data.frame(subgroup=sg$name,outcome=outcome_col,
    type=outcome_type,n_trt=nt,n_ctl=nc,est=NA,se=NA,p=NA,
    rate_trt=NA,rate_ctl=NA,rd=NA,nnt=NA,stringsAsFactors=FALSE))

  if (outcome_type=="continuous") {
    sub<-sub[!is.na(sub[[outcome_col]]),]; if(nrow(sub)<30) return(NULL)
    fit<-lm(as.formula(paste(outcome_col,"~treated")),data=sub)
    ct<-safe_ct(fit); if(is.null(ct)) return(NULL)
    return(data.frame(subgroup=sg$name,outcome=outcome_col,type="continuous",
      n_trt=sum(sub$treated==1),n_ctl=sum(sub$treated==0),
      est=ct["treated",1],se=ct["treated",2],p=ct["treated",4],
      rate_trt=NA,rate_ctl=NA,rd=NA,nnt=NA,stringsAsFactors=FALSE))
  } else {
    sub<-sub[!is.na(sub[[outcome_col]]),]
    r1<-mean(sub[[outcome_col]][sub$treated==1]); r0<-mean(sub[[outcome_col]][sub$treated==0])
    rd<-r1-r0; nnt<-if(abs(rd)>0.001) round(1/abs(rd)) else NA
    fit<-tryCatch(glm(as.formula(paste(outcome_col,"~treated")),data=sub,family=quasibinomial()),
                  error=function(e)NULL)
    if(is.null(fit)||!"treated"%in%names(coef(fit))) return(NULL)
    ct<-safe_ct(fit); if(is.null(ct)) return(NULL)
    return(data.frame(subgroup=sg$name,outcome=outcome_col,type="binary",
      n_trt=sum(sub$treated==1),n_ctl=sum(sub$treated==0),
      est=exp(ct["treated",1]),se=ct["treated",2],p=ct["treated",4],
      rate_trt=r1,rate_ctl=r0,rd=rd,nnt=nnt,stringsAsFactors=FALSE))
  }
}

test_interaction <- function(df, sg_var, outcome_col) {
  if(!sg_var%in%names(df)||!outcome_col%in%names(df)) return(NA)
  d<-df[!is.na(df[[sg_var]])&!is.na(df[[outcome_col]]),]; if(nrow(d)<60) return(NA)
  fml<-as.formula(sprintf("%s ~ treated * %s",outcome_col,sg_var))
  fit<-tryCatch(lm(fml,data=d),error=function(e)NULL); if(is.null(fit)) return(NA)
  ct<-tryCatch(coeftest(fit,vcov.=vcovHC(fit,type="HC1")),error=function(e)NULL)
  if(is.null(ct)) return(NA)
  ir<-grep(paste0("treated:",sg_var),rownames(ct)); if(length(ir)==0) return(NA)
  ct[ir[1],4]
}

# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════
args<-commandArgs(trailingOnly=TRUE)
if(length(args)<1){cat("Usage: Rscript 03_hte.R <db>\n");quit(status=1)}
db<-toupper(args[1]); tag<-tolower(db)

SEP<-paste(rep("=",70),collapse="")
cat(sprintf("\n%s\n03_hte.R — HTE: %s\n  Primary matched set (19var LAST, no K/Mg, yet-untreated)\n%s\n",SEP,db,SEP))

pairs<-read.csv(file.path(RESULTS,sprintf("did_pairs_primary_yet_untreated_%s.csv",tag)),stringsAsFactors=FALSE)
all_pts<-read.csv(file.path(RESULTS,sprintf("did_all_%s.csv",tag)),stringsAsFactors=FALSE)
cr_all<-read.csv(file.path(RESULTS,sprintf("did_cr_all_%s.csv",tag)),stringsAsFactors=FALSE)
cr_id<-if("patientunitstayid"%in%names(cr_all))"patientunitstayid" else "stay_id"
cr_all$pid<-cr_all[[cr_id]]
if(!"offset_h"%in%names(cr_all)) cr_all$offset_h<-cr_all$labresultoffset/60
cr_all<-cr_all[order(cr_all$pid,cr_all$offset_h),]
cr_list<-split(cr_all[,c("labresult","offset_h")],cr_all$pid)
cat(sprintf("  Pairs: %d | Patients: %d\n",nrow(pairs),nrow(all_pts)))

# LAST labs for subgroup defs
cat("  Computing LAST labs...\n")
labs_raw<-read.csv(file.path(RESULTS,sprintf("did_labs_all_%s.csv",tag)),stringsAsFactors=FALSE)
pcl<-if("patientunitstayid"%in%names(labs_raw))"patientunitstayid" else "stay_id"
if(pcl%in%names(labs_raw)) labs_raw$pid<-labs_raw[[pcl]]
le<-labs_raw[labs_raw$lab_name%in%c("magnesium","potassium","calcium","lactate"),]
lh<-labs_raw[labs_raw$lab_name=="heartrate",]
if(nrow(lh)>500000){lh$hb<-floor(lh$offset_h);lh<-lh[order(lh$pid,lh$offset_h),]
  lh<-lh[!duplicated(paste(lh$pid,lh$hb)),];lh$hb<-NULL}
labs<-rbind(le,lh); rm(labs_raw,le,lh)
for(ln in LAB_BASES){
  sub<-labs[labs$lab_name==ln,]; if(nrow(sub)==0) next
  sub$mg_offset_h<-all_pts$mg_offset_h[match(sub$pid,all_pts$pid)]
  sub<-sub[sub$offset_h>=0&(is.na(sub$mg_offset_h)|sub$offset_h<sub$mg_offset_h),]
  if(nrow(sub)==0) next
  s<-sub[order(-sub$offset_h),]; s<-s[!duplicated(s$pid),]
  all_pts[[paste0("last_",ln)]]<-s$value[match(all_pts$pid,s$pid)]
}

# Build outcome dataset
cat("  Building outcome dataset...\n")
rows<-vector("list",nrow(pairs)*2)
ri<-0
for(i in seq_len(nrow(pairs))){
  tp<-as.character(pairs$trt_pid[i]); cp<-as.character(pairs$ctl_pid[i]); tmg<-pairs$t_mg[i]
  cr_pre_t<-find_cr_pre(cr_list[[tp]],tmg); cr_pre_c<-find_cr_pre(cr_list[[cp]],tmg)
  cr48_t<-find_cr(cr_list[[tp]],tmg+48); cr48_c<-find_cr(cr_list[[cp]],tmg+48)
  dcr48_t<-if(!is.na(cr48_t)&&!is.na(cr_pre_t)) cr48_t-cr_pre_t else NA
  dcr48_c<-if(!is.na(cr48_c)&&!is.na(cr_pre_c)) cr48_c-cr_pre_c else NA
  max48_t<-max_cr_window(cr_list[[tp]],tmg,48); max48_c<-max_cr_window(cr_list[[cp]],tmg,48)
  aki48_t<-as.integer(!is.na(max48_t)&&!is.na(cr_pre_t)&&(max48_t-cr_pre_t)>=0.3)
  aki48_c<-as.integer(!is.na(max48_c)&&!is.na(cr_pre_c)&&(max48_c-cr_pre_c)>=0.3)
  max7d_t<-max_cr_window(cr_list[[tp]],tmg,168); max7d_c<-max_cr_window(cr_list[[cp]],tmg,168)
  aki7d_t<-as.integer(!is.na(max7d_t)&&!is.na(cr_pre_t)&&cr_pre_t>0&&(max7d_t/cr_pre_t)>=1.5)
  aki7d_c<-as.integer(!is.na(max7d_c)&&!is.na(cr_pre_c)&&cr_pre_c>0&&(max7d_c/cr_pre_c)>=1.5)
  ti<-which(all_pts$pid==pairs$trt_pid[i])[1]; ci<-which(all_pts$pid==pairs$ctl_pid[i])[1]
  if(is.na(ti)||is.na(ci)) next
  base<-data.frame(pair_id=i,t_mg=tmg,stringsAsFactors=FALSE)
  covs<-c("age","is_female","bmi","egfr","surg_cabg","surg_valve","diabetes","ckd",
          "heart_failure","hosp_mortality","poaf","encephalopathy","vent_arrhythmia","last_magnesium")
  for(cv in covs) base[[cv]]<-NA
  mk_row<-function(idx,trt,dcr,a48,a7d,crp){
    r<-base; r$treated<-trt; r$pid<-all_pts$pid[idx]; r$dcr_48h<-dcr
    r$aki_48h<-a48; r$aki_7d<-a7d; r$cr_pre<-crp
    for(cv in covs) if(cv%in%names(all_pts)) r[[cv]]<-all_pts[[cv]][idx]
    r}
  ri<-ri+1; rows[[ri]]<-mk_row(ti,1,dcr48_t,aki48_t,aki7d_t,cr_pre_t)
  ri<-ri+1; rows[[ri]]<-mk_row(ci,0,dcr48_c,aki48_c,aki7d_c,cr_pre_c)
}
df<-do.call(rbind,rows[1:ri])
cat(sprintf("  Outcome dataset: %d rows (%d pairs)\n",nrow(df),nrow(df)/2))

# Overall rates
cat("\n  Overall outcome rates:\n")
for(oc in c("dcr_48h","aki_48h","aki_7d","hosp_mortality","poaf","encephalopathy","vent_arrhythmia")){
  if(!oc%in%names(df)) next
  r1<-mean(df[[oc]][df$treated==1],na.rm=T); r0<-mean(df[[oc]][df$treated==0],na.rm=T)
  if(oc=="dcr_48h") cat(sprintf("    %s: trt=%+.4f, ctl=%+.4f, DiD=%+.4f\n",oc,r1,r0,r1-r0))
  else cat(sprintf("    %s: trt=%.1f%%, ctl=%.1f%%, RD=%+.1f%%\n",oc,100*r1,100*r0,100*(r1-r0)))
}

# ═══════════════════════════════════════════════════════════════════
# SUBGROUP × OUTCOME SWEEP
# ═══════════════════════════════════════════════════════════════════
OUTCOMES<-list(
  list(col="dcr_48h",type="continuous",label="dCr 48h"),
  list(col="aki_48h",type="binary",label="48h AKI (>=0.3)"),
  list(col="aki_7d",type="binary",label="7d AKI (ratio>=1.5)"),
  list(col="hosp_mortality",type="binary",label="Mortality"),
  list(col="poaf",type="binary",label="POAF"),
  list(col="encephalopathy",type="binary",label="Encephalopathy"),
  list(col="vent_arrhythmia",type="binary",label="Vent arrhythmia")
)

all_hte<-list()
for(sg in SUBGROUPS) for(oc in OUTCOMES){
  r<-test_subgroup(df,sg,oc$col,oc$type)
  if(!is.null(r)) all_hte[[length(all_hte)+1]]<-r
}
hte<-do.call(rbind,all_hte)

# Print
for(oc in OUTCOMES){
  cat(sprintf("\n  == %s ==\n",oc$label))
  sub<-hte[hte$outcome==oc$col,]
  if(oc$type=="continuous"){
    cat(sprintf("  %-20s  %5s %5s  %9s  %8s\n","Subgroup","nT","nC","DiD","P"))
    for(i in seq_len(nrow(sub))){r<-sub[i,]
      sig<-if(!is.na(r$p)&&r$p<0.05)" *" else "  "
      cat(sprintf("  %-20s  %5d %5d  %+.4f   %.4f%s\n",r$subgroup,r$n_trt,r$n_ctl,
                  ifelse(is.na(r$est),0,r$est),ifelse(is.na(r$p),1,r$p),sig))}
  } else {
    cat(sprintf("  %-20s  %5s %5s  %5s %5s %6s  %5s  %8s\n",
                "Subgroup","nT","nC","T%","C%","RD%","OR","P"))
    for(i in seq_len(nrow(sub))){r<-sub[i,]
      sig<-if(!is.na(r$p)&&r$p<0.05)" *" else "  "
      cat(sprintf("  %-20s  %5d %5d  %5.1f %5.1f %+5.1f  %5.2f  %.4f%s\n",r$subgroup,r$n_trt,r$n_ctl,
                  100*ifelse(is.na(r$rate_trt),0,r$rate_trt),100*ifelse(is.na(r$rate_ctl),0,r$rate_ctl),
                  100*ifelse(is.na(r$rd),0,r$rd),ifelse(is.na(r$est),1,r$est),
                  ifelse(is.na(r$p),1,r$p),sig))}
  }
}

# Interactions
cat(sprintf("\n  == INTERACTION TESTS (dCr 48h) ==\n"))
df$age_ge65<-as.integer(df$age>=65)
df$egfr_lt60<-as.integer(df$egfr<60)
df$mg_lt18<-as.integer(!is.na(df$last_magnesium)&df$last_magnesium<1.8)
df$bmi_ge30<-as.integer(!is.na(df$bmi)&df$bmi>=30)
for(iv in c("age_ge65","egfr_lt60","mg_lt18","surg_cabg","diabetes","ckd","heart_failure","bmi_ge30")){
  p<-test_interaction(df,iv,"dcr_48h")
  sig<-if(!is.na(p)&&p<0.05)" *" else "  "
  cat(sprintf("    treated x %-16s  P = %.4f%s\n",iv,ifelse(is.na(p),NA,p),sig))
}

# Save
write.csv(hte,file.path(RESULTS,sprintf("did_hte_%s.csv",tag)),row.names=FALSE)
write.csv(df,file.path(RESULTS,sprintf("did_hte_data_%s.csv",tag)),row.names=FALSE)
cat(sprintf("\n%s\n03_hte.R — %s DONE\n  did_hte_%s.csv (%d rows)\n  did_hte_data_%s.csv\n%s\n",
            SEP,db,tag,nrow(hte),tag,SEP))
