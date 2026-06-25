#!/usr/bin/env Rscript
# ============================================================================
# 03b_egfr_aki_stages.R — eGFR-stratified AKI stages + secondary outcomes
#
# Runs on existing primary matched pairs from 02_psm.R.
# No new matching — just outcome computation by eGFR stratum.
#
# Outputs:
#   results/egfr_aki_stages_{db}.csv  — OR, CI, P for all strata × outcomes
#
# Usage: Rscript 03b_egfr_aki_stages.R mimic
#        Rscript 03b_egfr_aki_stages.R eicu
# ============================================================================

suppressPackageStartupMessages({ library(sandwich); library(lmtest) })

RESULTS <- path.expand("~/mg_aki/results")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) { cat("Usage: Rscript 03b_egfr_aki_stages.R <db>\n"); quit(status=1) }
tag <- tolower(args[1]); db <- toupper(tag)

SEP <- paste(rep("=", 60), collapse = "")
cat(sprintf("\n%s\n03b_egfr_aki_stages.R — %s\n%s\n", SEP, db, SEP))

# ── Load data ─────────────────────────────────────────────────────
all_pts <- read.csv(file.path(RESULTS, sprintf("did_all_%s.csv", tag)), stringsAsFactors=FALSE)
pairs   <- read.csv(file.path(RESULTS, sprintf("did_pairs_primary_yet_untreated_%s.csv", tag)),
                     stringsAsFactors=FALSE)
cr_all  <- read.csv(file.path(RESULTS, sprintf("did_cr_all_%s.csv", tag)), stringsAsFactors=FALSE)

cr_id <- if ("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
cr_all$pid <- cr_all[[cr_id]]
if (!"offset_h" %in% names(cr_all)) cr_all$offset_h <- cr_all$labresultoffset / 60
cr_all <- cr_all[order(cr_all$pid, cr_all$offset_h), ]
cr_list <- split(cr_all[, c("labresult","offset_h")], cr_all$pid)

trt_rows <- match(pairs$trt_pid, all_pts$pid)
ctl_rows <- match(pairs$ctl_pid, all_pts$pid)
n_pairs  <- nrow(pairs)

cat(sprintf("  Pairs: %d | Patients: %d\n", n_pairs, nrow(all_pts)))

# ── OR helper ─────────────────────────────────────────────────────
run_or <- function(outcome_trt, outcome_ctl) {
  valid <- !is.na(outcome_trt) & !is.na(outcome_ctl)
  ot <- outcome_trt[valid]; oc <- outcome_ctl[valid]
  n <- sum(valid); nt <- sum(ot); nc <- sum(oc)
  if (n < 30 || (nt + nc) == 0)
    return(data.frame(or=NA, or_lo=NA, or_hi=NA, p=NA,
                      rate_trt=NA, rate_ctl=NA, n=n, events_trt=nt, events_ctl=nc))
  r1 <- mean(ot); r0 <- mean(oc)
  df <- data.frame(outcome=c(ot,oc), treated=rep(c(1,0), each=sum(valid)))
  fit <- tryCatch(glm(outcome ~ treated, data=df, family=quasibinomial()),
                  error=function(e) NULL)
  if (is.null(fit)) return(data.frame(or=NA, or_lo=NA, or_hi=NA, p=NA,
                                       rate_trt=r1, rate_ctl=r0, n=n, events_trt=nt, events_ctl=nc))
  ct <- tryCatch(coeftest(fit, vcov.=vcovHC(fit, type="HC1")),
                 error=function(e) tryCatch(coeftest(fit), error=function(e2) NULL))
  if (is.null(ct)) return(data.frame(or=NA, or_lo=NA, or_hi=NA, p=NA,
                                      rate_trt=r1, rate_ctl=r0, n=n, events_trt=nt, events_ctl=nc))
  or <- exp(ct["treated","Estimate"])
  lo <- exp(ct["treated","Estimate"] - 1.96*ct["treated","Std. Error"])
  hi <- exp(ct["treated","Estimate"] + 1.96*ct["treated","Std. Error"])
  p  <- ct["treated", ncol(ct)]
  data.frame(or=round(or,4), or_lo=round(lo,4), or_hi=round(hi,4), p=round(p,6),
             rate_trt=round(r1,4), rate_ctl=round(r0,4), n=n, events_trt=nt, events_ctl=nc)
}

# ── Compute AKI stages ────────────────────────────────────────────
compute_aki_stages <- function(pid, t_mg) {
  cr <- cr_list[[as.character(pid)]]
  if (is.null(cr) || nrow(cr) < 1) return(c(NA, NA, NA))
  pre <- cr[cr$offset_h >= 0 & cr$offset_h < t_mg, ]
  if (nrow(pre) == 0) return(c(NA, NA, NA))
  bl <- pre$labresult[which.max(pre$offset_h)]
  if (is.na(bl) || bl <= 0) return(c(NA, NA, NA))
  post <- cr[cr$offset_h >= t_mg & cr$offset_h <= (t_mg + 168), ]
  if (nrow(post) == 0) return(c(0, 0, 0))
  stage1 <- 0; stage2 <- 0; stage3 <- 0
  for (i in seq_len(nrow(post))) {
    h <- post$offset_h[i] - t_mg; val <- post$labresult[i]
    delta <- val - bl; ratio <- val / bl
    if (h <= 48 && (delta >= 0.3 || ratio >= 1.5)) stage1 <- 1
    if (h > 48 && ratio >= 1.5) stage1 <- 1
    if (ratio >= 2.0) stage2 <- 1
    if (ratio >= 3.0 || val >= 4.0) stage3 <- 1
  }
  c(stage1, stage2, stage3)
}

cat("  Computing AKI stages...\n")
aki_trt <- aki_ctl <- matrix(NA, n_pairs, 3)
for (i in seq_len(n_pairs)) {
  aki_trt[i,] <- compute_aki_stages(pairs$trt_pid[i], pairs$t_mg[i])
  aki_ctl[i,] <- compute_aki_stages(pairs$ctl_pid[i], pairs$t_mg[i])
}

# ── eGFR strata ───────────────────────────────────────────────────
egfr_trt <- all_pts$egfr[trt_rows]
ckd_stage <- rep(NA_character_, n_pairs)
ckd_stage[!is.na(egfr_trt) & egfr_trt >= 90]                  <- "eGFR>=90"
ckd_stage[!is.na(egfr_trt) & egfr_trt >= 60 & egfr_trt < 90]  <- "eGFR_60-89"
ckd_stage[!is.na(egfr_trt) & egfr_trt >= 45 & egfr_trt < 60]  <- "eGFR_45-59"
ckd_stage[!is.na(egfr_trt) & egfr_trt >= 30 & egfr_trt < 45]  <- "eGFR_30-44"
ckd_stage[!is.na(egfr_trt) & egfr_trt < 30]                    <- "eGFR<30"

strata_order <- c("Overall","eGFR>=90","eGFR_60-89","eGFR_45-59","eGFR_30-44","eGFR<30")

cat(sprintf("  eGFR distribution: %s\n",
            paste(names(table(ckd_stage)), table(ckd_stage), sep="=", collapse=", ")))

# ── Run all combinations ──────────────────────────────────────────
aki_labels <- c("AKI_Stage1+","AKI_Stage2+","AKI_Stage3+")
outcome_cols <- c("hosp_mortality","encephalopathy","vent_arrhythmia","poaf")

results <- list()
row_idx <- 0

# AKI stages × eGFR strata
for (aki_col in 1:3) {
  for (stg in strata_order) {
    if (stg == "Overall") { idx <- seq_len(n_pairs) }
    else { idx <- which(ckd_stage == stg) }
    if (length(idx) < 30) next
    res <- run_or(aki_trt[idx, aki_col], aki_ctl[idx, aki_col])
    res$outcome <- aki_labels[aki_col]; res$stratum <- stg; res$db <- db
    row_idx <- row_idx + 1; results[[row_idx]] <- res
  }
}

# Secondary outcomes × eGFR strata
for (oc in outcome_cols) {
  if (!(oc %in% names(all_pts))) next
  for (stg in strata_order) {
    if (stg == "Overall") { idx <- seq_len(n_pairs) }
    else { idx <- which(ckd_stage == stg) }
    if (length(idx) < 30) next
    ot <- all_pts[[oc]][trt_rows[idx]]
    oc_val <- all_pts[[oc]][ctl_rows[idx]]
    res <- run_or(ot, oc_val)
    res$outcome <- oc; res$stratum <- stg; res$db <- db
    row_idx <- row_idx + 1; results[[row_idx]] <- res
  }
}

res_df <- do.call(rbind, results)
res_df$stratum <- factor(res_df$stratum, levels = strata_order)

# ── Print summary ─────────────────────────────────────────────────
cat(sprintf("\n── Results (%d rows) ──\n", nrow(res_df)))
for (oc in unique(res_df$outcome)) {
  cat(sprintf("\n  [%s]\n", oc))
  sub <- res_df[res_df$outcome == oc, ]
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    if (is.na(r$or)) { cat(sprintf("    %-15s  n=%d (skip)\n", r$stratum, r$n)); next }
    sig <- if(!is.na(r$p) && r$p < 0.05) " *" else "  "
    cat(sprintf("    %-15s  OR=%.3f [%.3f,%.3f]  P=%.4f%s  %.1f%% vs %.1f%%  n=%d\n",
                r$stratum, r$or, r$or_lo, r$or_hi, r$p, sig,
                100*r$rate_trt, 100*r$rate_ctl, r$n))
  }
}

# ── Save ──────────────────────────────────────────────────────────
outpath <- file.path(RESULTS, sprintf("egfr_aki_stages_%s.csv", tag))
write.csv(res_df, outpath, row.names = FALSE)
cat(sprintf("\n  Saved: %s (%d rows)\n", outpath, nrow(res_df)))

cat(sprintf("\n%s\n03b done: %s\n%s\n", SEP, db, SEP))
