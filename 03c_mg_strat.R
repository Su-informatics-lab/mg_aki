#!/usr/bin/env Rscript
# ============================================================================
# 03c_mg_strat.R — Mg-stratified AKI + eGFR × Mg cross-stratification
#
# Answers Dr. Su's question: is the treatment effect driven by
# SUPPLEMENTING Mg, or by HAVING different baseline Mg levels?
#
# If effect persists within Mg strata → supplementation is causal
# If effect disappears → confounded by baseline Mg level
# If eGFR reversal persists within Mg strata → eGFR is real modifier
#
# Pair-preserving: subsets by TREATED patient's baseline Mg/eGFR.
#
# Baseline Mg = first serum Mg measurement before IV Mg time (treated)
#               from did_labs_all_{db}.csv (lab_name == "magnesium").
#
# Outputs:
#   results/mg_strat_{db}.csv          — all results
#
# Usage: Rscript 03c_mg_strat.R mimic
#        Rscript 03c_mg_strat.R eicu
# ============================================================================

suppressPackageStartupMessages({ library(sandwich); library(lmtest) })

RESULTS <- path.expand("~/mg_aki/results")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) { cat("Usage: Rscript 03c_mg_strat.R <db>\n"); quit(status = 1) }
tag <- tolower(args[1]); db <- toupper(tag)

SEP <- paste(rep("=", 60), collapse = "")
cat(sprintf("\n%s\n03c_mg_strat.R — %s\n  Q: Is treatment effect confounded by baseline Mg level?\n%s\n",
            SEP, db, SEP))

# ── Load ──────────────────────────────────────────────────────────
all_pts <- read.csv(file.path(RESULTS, sprintf("did_all_%s.csv", tag)), stringsAsFactors = FALSE)
pairs   <- read.csv(file.path(RESULTS, sprintf("did_pairs_primary_yet_untreated_%s.csv", tag)),
                     stringsAsFactors = FALSE)
cr_all  <- read.csv(file.path(RESULTS, sprintf("did_cr_all_%s.csv", tag)), stringsAsFactors = FALSE)
labs_all <- read.csv(file.path(RESULTS, sprintf("did_labs_all_%s.csv", tag)), stringsAsFactors = FALSE)

# Identify pid column in labs_all (patientunitstayid for eICU, stay_id for MIMIC)
labs_pid_col <- if ("patientunitstayid" %in% names(labs_all)) "patientunitstayid" else "stay_id"
labs_all$pid <- labs_all[[labs_pid_col]]

cr_id <- if ("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
cr_all$pid <- cr_all[[cr_id]]
if (!"offset_h" %in% names(cr_all)) cr_all$offset_h <- cr_all$labresultoffset / 60
cr_all <- cr_all[order(cr_all$pid, cr_all$offset_h), ]
cr_list <- split(cr_all[, c("labresult", "offset_h")], cr_all$pid)

trt_rows <- match(pairs$trt_pid, all_pts$pid)
ctl_rows <- match(pairs$ctl_pid, all_pts$pid)
n_pairs  <- nrow(pairs)

cat(sprintf("  Pairs: %d\n", n_pairs))

# ── Compute baseline Mg from did_labs_all ─────────────────────────
# Baseline Mg = first serum Mg before IV Mg time (t_mg) for treated patients
cat("  Computing baseline Mg from did_labs_all...\n")
mg_labs <- labs_all[labs_all$lab_name == "magnesium" & !is.na(labs_all$value), ]
mg_labs <- mg_labs[order(mg_labs$pid, mg_labs$offset_h), ]
mg_list <- split(mg_labs[, c("value", "offset_h")], mg_labs$pid)

cat(sprintf("  Mg measurements in labs_all: %d across %d patients\n",
            nrow(mg_labs), length(unique(mg_labs$pid))))

mg_trt <- rep(NA_real_, n_pairs)
for (i in seq_len(n_pairs)) {
  pid <- as.character(pairs$trt_pid[i])
  t_mg <- pairs$t_mg[i]  # IV Mg onset in hours
  mg_i <- mg_list[[pid]]
  if (is.null(mg_i) || nrow(mg_i) == 0) next
  # Baseline = last Mg measurement before IV Mg time
  pre <- mg_i[mg_i$offset_h < t_mg, , drop = FALSE]
  if (nrow(pre) == 0) next
  mg_trt[i] <- pre$value[nrow(pre)]  # last pre-treatment value
}

mg_avail <- sum(!is.na(mg_trt))
cat(sprintf("  Mg available: %d/%d (%.1f%%)\n", mg_avail, n_pairs, 100*mg_avail/n_pairs))
if (mg_avail > 0) {
  cat(sprintf("  Mg distribution: mean=%.2f, median=%.2f, IQR=[%.2f,%.2f]\n",
              mean(mg_trt, na.rm=TRUE), median(mg_trt, na.rm=TRUE),
              quantile(mg_trt, 0.25, na.rm=TRUE), quantile(mg_trt, 0.75, na.rm=TRUE)))
} else {
  cat("  WARNING: No baseline Mg found — check did_labs_all for magnesium rows\n")
}

# ── OR helper ─────────────────────────────────────────────────────
run_or <- function(ot, oc) {
  valid <- !is.na(ot) & !is.na(oc)
  ot <- ot[valid]; oc <- oc[valid]
  n <- sum(valid); et <- sum(ot); ec <- sum(oc)
  if (n < 30 || (et + ec) == 0)
    return(data.frame(or=NA, or_lo=NA, or_hi=NA, p=NA,
                      rate_trt=NA, rate_ctl=NA, n=n))
  r1 <- mean(ot); r0 <- mean(oc)
  df <- data.frame(outcome=c(ot,oc), treated=rep(c(1,0),each=sum(valid)))
  fit <- tryCatch(glm(outcome~treated, data=df, family=quasibinomial()),
                  error=function(e) NULL)
  if (is.null(fit)) return(data.frame(or=NA,or_lo=NA,or_hi=NA,p=NA,
                                       rate_trt=r1,rate_ctl=r0,n=n))
  ct <- tryCatch(coeftest(fit, vcov.=vcovHC(fit, type="HC1")),
                 error=function(e) tryCatch(coeftest(fit), error=function(e2) NULL))
  if (is.null(ct)) return(data.frame(or=NA,or_lo=NA,or_hi=NA,p=NA,
                                      rate_trt=r1,rate_ctl=r0,n=n))
  or <- exp(ct["treated","Estimate"])
  lo <- exp(ct["treated","Estimate"] - 1.96*ct["treated","Std. Error"])
  hi <- exp(ct["treated","Estimate"] + 1.96*ct["treated","Std. Error"])
  p  <- ct["treated", ncol(ct)]
  data.frame(or=round(or,4), or_lo=round(lo,4), or_hi=round(hi,4), p=round(p,6),
             rate_trt=round(r1,4), rate_ctl=round(r0,4), n=n)
}

# ── Compute 7d AKI ────────────────────────────────────────────────
compute_aki_7d <- function(pid, t_mg) {
  cr <- cr_list[[as.character(pid)]]
  if (is.null(cr) || nrow(cr) < 1) return(NA)
  pre <- cr[cr$offset_h >= 0 & cr$offset_h < t_mg, ]
  if (nrow(pre) == 0) return(NA)
  bl <- pre$labresult[which.max(pre$offset_h)]
  if (is.na(bl) || bl <= 0) return(NA)
  post <- cr[cr$offset_h >= t_mg & cr$offset_h <= (t_mg + 168), ]
  if (nrow(post) == 0) return(0)
  for (i in seq_len(nrow(post))) {
    h <- post$offset_h[i] - t_mg; val <- post$labresult[i]
    delta <- val - bl; ratio <- val / bl
    if (h <= 48 && (delta >= 0.3 || ratio >= 1.5)) return(1)
    if (h > 48 && ratio >= 1.5) return(1)
  }
  return(0)
}

cat("  Computing 7d AKI...\n")
aki_trt <- aki_ctl <- integer(n_pairs)
for (i in seq_len(n_pairs)) {
  aki_trt[i] <- compute_aki_7d(pairs$trt_pid[i], pairs$t_mg[i])
  aki_ctl[i] <- compute_aki_7d(pairs$ctl_pid[i], pairs$t_mg[i])
}

# Also load mortality
mort_trt <- all_pts$hosp_mortality[trt_rows]
mort_ctl <- all_pts$hosp_mortality[ctl_rows]

# ── Treated patient covariates ────────────────────────────────────
egfr_trt <- all_pts$egfr[trt_rows]

# ── Define strata ─────────────────────────────────────────────────
mg_strat <- rep(NA_character_, n_pairs)
mg_strat[!is.na(mg_trt) & mg_trt < 1.6]                    <- "Mg<1.6"
mg_strat[!is.na(mg_trt) & mg_trt >= 1.6 & mg_trt < 1.8]   <- "Mg_1.6-1.8"
mg_strat[!is.na(mg_trt) & mg_trt >= 1.8 & mg_trt < 2.0]   <- "Mg_1.8-2.0"
mg_strat[!is.na(mg_trt) & mg_trt >= 2.0 & mg_trt < 2.3]   <- "Mg_2.0-2.3"
mg_strat[!is.na(mg_trt) & mg_trt >= 2.3]                   <- "Mg>=2.3"

egfr_strat <- rep(NA_character_, n_pairs)
egfr_strat[!is.na(egfr_trt) & egfr_trt >= 90]                 <- "eGFR>=90"
egfr_strat[!is.na(egfr_trt) & egfr_trt >= 60 & egfr_trt < 90] <- "eGFR_60-89"
egfr_strat[!is.na(egfr_trt) & egfr_trt >= 45 & egfr_trt < 60] <- "eGFR_45-59"
egfr_strat[!is.na(egfr_trt) & egfr_trt < 45]                  <- "eGFR<45"

cat(sprintf("  Mg strata: %s\n", paste(names(table(mg_strat)), table(mg_strat), sep="=", collapse=", ")))
cat(sprintf("  eGFR strata: %s\n", paste(names(table(egfr_strat)), table(egfr_strat), sep="=", collapse=", ")))

# ══════════════════════════════════════════════════════════════════
# SECTION 1: AKI BY MG STRATUM (answers: supplementation or level?)
# ══════════════════════════════════════════════════════════════════
cat(sprintf("\n%s\nSECTION 1: AKI 7d by Mg stratum\n%s\n",
            paste(rep("-",50),collapse=""), paste(rep("-",50),collapse="")))
cat("  If effect persists → supplementation is causal\n")
cat("  If disappears → confounded by baseline Mg\n\n")

mg_order <- c("Overall","Mg<1.6","Mg_1.6-1.8","Mg_1.8-2.0","Mg_2.0-2.3","Mg>=2.3")
results <- list(); ridx <- 0

for (stg in mg_order) {
  if (stg == "Overall") idx <- seq_len(n_pairs)
  else idx <- which(mg_strat == stg)
  if (length(idx) < 30) { cat(sprintf("  %-15s  n=%d (skip)\n", stg, length(idx))); next }

  res_aki <- run_or(aki_trt[idx], aki_ctl[idx])
  res_aki$outcome <- "aki_7d"; res_aki$mg_strat <- stg; res_aki$egfr_strat <- "All"; res_aki$db <- db
  ridx <- ridx+1; results[[ridx]] <- res_aki
  sig <- if(!is.na(res_aki$p) && res_aki$p < 0.05) " *" else "  "
  cat(sprintf("  %-15s  AKI: OR=%.3f [%.3f,%.3f] P=%.4f%s  %.1f%% vs %.1f%%  n=%d\n",
              stg, res_aki$or, res_aki$or_lo, res_aki$or_hi, res_aki$p, sig,
              100*res_aki$rate_trt, 100*res_aki$rate_ctl, res_aki$n))
}

# ══════════════════════════════════════════════════════════════════
# SECTION 2: MORTALITY BY MG STRATUM
# ══════════════════════════════════════════════════════════════════
cat(sprintf("\n%s\nSECTION 2: Mortality by Mg stratum\n%s\n",
            paste(rep("-",50),collapse=""), paste(rep("-",50),collapse="")))

for (stg in mg_order) {
  if (stg == "Overall") idx <- seq_len(n_pairs)
  else idx <- which(mg_strat == stg)
  if (length(idx) < 30) next
  res <- run_or(mort_trt[idx], mort_ctl[idx])
  res$outcome <- "mortality"; res$mg_strat <- stg; res$egfr_strat <- "All"; res$db <- db
  ridx <- ridx+1; results[[ridx]] <- res
  if (is.na(res$or)) next
  sig <- if(!is.na(res$p) && res$p < 0.05) " *" else "  "
  cat(sprintf("  %-15s  Mort: OR=%.3f [%.3f,%.3f] P=%.4f%s  %.1f%% vs %.1f%%  n=%d\n",
              stg, res$or, res$or_lo, res$or_hi, res$p, sig,
              100*res$rate_trt, 100*res$rate_ctl, res$n))
}

# ══════════════════════════════════════════════════════════════════
# SECTION 3: eGFR × Mg CROSS-STRATIFICATION (the key test)
# ══════════════════════════════════════════════════════════════════
cat(sprintf("\n%s\nSECTION 3: eGFR × Mg cross-stratification (AKI 7d)\n",
            paste(rep("=",60),collapse="")))
cat("  Key Q: Does eGFR reversal persist within each Mg stratum?\n")
cat(sprintf("%s\n", paste(rep("=",60),collapse="")))

egfr_order <- c("eGFR>=90","eGFR_60-89","eGFR_45-59","eGFR<45")
mg_coarse <- c("Mg<1.8","Mg>=1.8")

# Coarse Mg split for cross-strat (more power)
mg_coarse_strat <- rep(NA_character_, n_pairs)
mg_coarse_strat[!is.na(mg_trt) & mg_trt < 1.8]  <- "Mg<1.8"
mg_coarse_strat[!is.na(mg_trt) & mg_trt >= 1.8] <- "Mg>=1.8"

cat(sprintf("\n  Mg coarse: %s\n",
            paste(names(table(mg_coarse_strat)), table(mg_coarse_strat), sep="=", collapse=", ")))

cat(sprintf("\n  %15s | %15s | %15s | %15s | %15s\n",
            "", "eGFR>=90", "eGFR 60-89", "eGFR 45-59", "eGFR<45"))
cat(paste(rep("-", 85), collapse=""), "\n")

for (mg_s in mg_coarse) {
  row_str <- sprintf("  %15s |", mg_s)
  for (eg_s in egfr_order) {
    idx <- which(mg_coarse_strat == mg_s & egfr_strat == eg_s)
    if (length(idx) < 30) {
      row_str <- paste0(row_str, sprintf(" %15s |", sprintf("n=%d", length(idx))))
      next
    }
    res <- run_or(aki_trt[idx], aki_ctl[idx])
    res$outcome <- "aki_7d"; res$mg_strat <- mg_s; res$egfr_strat <- eg_s; res$db <- db
    ridx <- ridx+1; results[[ridx]] <- res
    if (is.na(res$or)) {
      row_str <- paste0(row_str, sprintf(" %15s |", "NA"))
    } else {
      sig <- if(!is.na(res$p) && res$p < 0.05) "*" else " "
      row_str <- paste0(row_str, sprintf(" OR=%.2f P=%.3f%s|", res$or, res$p, sig))
    }
  }
  cat(row_str, "\n")
}

# ── Same for mortality ────────────────────────────────────────────
cat(sprintf("\n  %15s | %15s | %15s | %15s | %15s\n",
            "MORTALITY", "eGFR>=90", "eGFR 60-89", "eGFR 45-59", "eGFR<45"))
cat(paste(rep("-", 85), collapse=""), "\n")

for (mg_s in mg_coarse) {
  row_str <- sprintf("  %15s |", mg_s)
  for (eg_s in egfr_order) {
    idx <- which(mg_coarse_strat == mg_s & egfr_strat == eg_s)
    if (length(idx) < 30) {
      row_str <- paste0(row_str, sprintf(" %15s |", sprintf("n=%d", length(idx))))
      next
    }
    res <- run_or(mort_trt[idx], mort_ctl[idx])
    res$outcome <- "mortality"; res$mg_strat <- mg_s; res$egfr_strat <- eg_s; res$db <- db
    ridx <- ridx+1; results[[ridx]] <- res
    if (is.na(res$or)) {
      row_str <- paste0(row_str, sprintf(" %15s |", "NA"))
    } else {
      sig <- if(!is.na(res$p) && res$p < 0.05) "*" else " "
      row_str <- paste0(row_str, sprintf(" OR=%.2f P=%.3f%s|", res$or, res$p, sig))
    }
  }
  cat(row_str, "\n")
}

# ══════════════════════════════════════════════════════════════════
# SECTION 3b: eGFR × Mg 3-BIN cross-stratification
#   Mg<1.6 (severe) | Mg 1.6-2.0 (sweet spot) | Mg≥2.0 (normal+high)
# ══════════════════════════════════════════════════════════════════
cat(sprintf("\n%s\nSECTION 3b: eGFR × Mg 3-bin cross-stratification\n%s\n",
            paste(rep("=",60),collapse=""), paste(rep("=",60),collapse="")))

mg3_names <- c("Mg<1.6", "Mg_1.6-2.0", "Mg>=2.0")
mg3_strat <- rep(NA_character_, n_pairs)
mg3_strat[!is.na(mg_trt) & mg_trt < 1.6]                 <- "Mg<1.6"
mg3_strat[!is.na(mg_trt) & mg_trt >= 1.6 & mg_trt < 2.0] <- "Mg_1.6-2.0"
mg3_strat[!is.na(mg_trt) & mg_trt >= 2.0]                 <- "Mg>=2.0"

cat(sprintf("  Mg 3-bin: %s\n",
            paste(names(table(mg3_strat)), table(mg3_strat), sep="=", collapse=", ")))

for (outcome_pair in list(
  list(label="AKI", ot=aki_trt, oc=aki_ctl, oname="aki_7d"),
  list(label="MORTALITY", ot=mort_trt, oc=mort_ctl, oname="mortality")
)) {
  cat(sprintf("\n  %10s | %15s | %15s | %15s | %15s\n",
              outcome_pair$label, "eGFR>=90", "eGFR 60-89", "eGFR 45-59", "eGFR<45"))
  cat(paste(rep("-", 80), collapse=""), "\n")

  for (mg_s in mg3_names) {
    row_str <- sprintf("  %10s |", mg_s)
    for (eg_s in egfr_order) {
      idx <- which(mg3_strat == mg_s & egfr_strat == eg_s)
      if (length(idx) < 30) {
        row_str <- paste0(row_str, sprintf(" %15s |", sprintf("n=%d", length(idx))))
        next
      }
      res <- run_or(outcome_pair$ot[idx], outcome_pair$oc[idx])
      res$outcome <- outcome_pair$oname; res$mg_strat <- mg_s
      res$egfr_strat <- eg_s; res$db <- db
      ridx <- ridx+1; results[[ridx]] <- res
      if (is.na(res$or)) {
        row_str <- paste0(row_str, sprintf(" %15s |", "NA"))
      } else {
        sig <- if(!is.na(res$p) && res$p < 0.05) "*" else " "
        row_str <- paste0(row_str, sprintf(" OR=%.2f P=%.3f%s|", res$or, res$p, sig))
      }
    }
    cat(row_str, "\n")
  }
}

# ── Interaction test: Mg stratum × treatment ──────────────────────
cat(sprintf("\n%s\nSECTION 4: Interaction P-values\n%s\n",
            paste(rep("-",50),collapse=""), paste(rep("-",50),collapse="")))

# Mg × treatment interaction (continuous)
valid_mg <- !is.na(mg_trt) & !is.na(aki_trt) & !is.na(aki_ctl)
if (sum(valid_mg) >= 100) {
  df_int <- data.frame(
    outcome = c(aki_trt[valid_mg], aki_ctl[valid_mg]),
    treated = rep(c(1,0), each = sum(valid_mg)),
    mg_bl   = rep(mg_trt[valid_mg], 2)
  )
  fit_int <- tryCatch(
    glm(outcome ~ treated * mg_bl, data=df_int, family=quasibinomial()),
    error = function(e) NULL
  )
  if (!is.null(fit_int)) {
    ct_int <- tryCatch(coeftest(fit_int, vcov.=vcovHC(fit_int, type="HC1")),
                       error = function(e) coeftest(fit_int))
    int_row <- grep("treated:mg_bl", rownames(ct_int))
    if (length(int_row) > 0) {
      cat(sprintf("  Mg × treatment interaction (continuous): P=%.4f\n",
                  ct_int[int_row, ncol(ct_int)]))
    }
  }
}

# eGFR × treatment interaction (continuous) — for reference
valid_eg <- !is.na(egfr_trt) & !is.na(aki_trt) & !is.na(aki_ctl)
if (sum(valid_eg) >= 100) {
  df_int2 <- data.frame(
    outcome = c(aki_trt[valid_eg], aki_ctl[valid_eg]),
    treated = rep(c(1,0), each = sum(valid_eg)),
    egfr    = rep(egfr_trt[valid_eg], 2)
  )
  fit_int2 <- tryCatch(
    glm(outcome ~ treated * egfr, data=df_int2, family=quasibinomial()),
    error = function(e) NULL
  )
  if (!is.null(fit_int2)) {
    ct_int2 <- tryCatch(coeftest(fit_int2, vcov.=vcovHC(fit_int2, type="HC1")),
                        error = function(e) coeftest(fit_int2))
    int_row2 <- grep("treated:egfr", rownames(ct_int2))
    if (length(int_row2) > 0) {
      cat(sprintf("  eGFR × treatment interaction (continuous): P=%.4f\n",
                  ct_int2[int_row2, ncol(ct_int2)]))
    }
  }
}

# ── Save ──────────────────────────────────────────────────────────
if (ridx > 0) {
  res_df <- do.call(rbind, results)
  outpath <- file.path(RESULTS, sprintf("mg_strat_%s.csv", tag))
  write.csv(res_df, outpath, row.names = FALSE)
  cat(sprintf("\n  Saved: %s (%d rows)\n", outpath, nrow(res_df)))
} else {
  cat("\n  WARNING: No results to save\n")
}

cat(sprintf("\n%s\n03c_mg_strat.R — %s DONE\n%s\n", SEP, db, SEP))
