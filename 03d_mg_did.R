#!/usr/bin/env Rscript
# ============================================================================
# 03d_mg_did.R — Stratified DiD: eGFR + baseline Mg + Overall
#
# Post-treatment Cr = peak (max) Cr within [T₀, T₀+h] at 9 horizons:
#   6, 12, 18, 24, 30, 36, 42, 48h, and 7d (168h)
# DiD = mean(ΔCr_treated) − mean(ΔCr_control)
# where ΔCr = peak Cr in window − baseline Cr
#
# Baseline Mg = first postoperative serum Mg before T₀
#               (consistent with 03c_mg_strat.R and manuscript)
#
# Outputs:
#   results/did_stratified_{db}.csv — DiD by strat_var × stratum × time
#
# Usage: Rscript 03d_mg_did.R mimic
#        Rscript 03d_mg_did.R eicu
# ============================================================================

suppressPackageStartupMessages({ library(sandwich); library(lmtest) })

RESULTS <- path.expand("~/mg_aki/results")
TARGETS <- c(6, 12, 18, 24, 30, 36, 42, 48, 168)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) { cat("Usage: Rscript 03d_mg_did.R <db>\n"); quit(status = 1) }
tag <- tolower(args[1]); db <- toupper(tag)

SEP <- paste(rep("=", 60), collapse = "")
cat(sprintf("\n%s\n03d_mg_did.R — %s\n  Stratified DiD: eGFR + Mg + Overall\n%s\n", SEP, db, SEP))

# ── Load ──────────────────────────────────────────────────────────
all_pts  <- read.csv(file.path(RESULTS, sprintf("did_all_%s.csv", tag)), stringsAsFactors = FALSE)
pairs    <- read.csv(file.path(RESULTS, sprintf("did_pairs_primary_yet_untreated_%s.csv", tag)),
                      stringsAsFactors = FALSE)
cr_all   <- read.csv(file.path(RESULTS, sprintf("did_cr_all_%s.csv", tag)), stringsAsFactors = FALSE)
labs_all <- read.csv(file.path(RESULTS, sprintf("did_labs_all_%s.csv", tag)), stringsAsFactors = FALSE)

cr_id <- if ("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
cr_all$pid <- cr_all[[cr_id]]
if (!"offset_h" %in% names(cr_all)) cr_all$offset_h <- cr_all$labresultoffset / 60
cr_all <- cr_all[order(cr_all$pid, cr_all$offset_h), ]
cr_list <- split(cr_all[, c("labresult", "offset_h")], cr_all$pid)

labs_pid <- if ("patientunitstayid" %in% names(labs_all)) "patientunitstayid" else "stay_id"
labs_all$pid <- labs_all[[labs_pid]]

trt_rows <- match(pairs$trt_pid, all_pts$pid)
n_pairs  <- nrow(pairs)
cat(sprintf("  Pairs: %d\n", n_pairs))

# ── Cr helpers ────────────────────────────────────────────────────
find_cr_pre <- function(cr_pt, t_h) {
  if (is.null(cr_pt) || nrow(cr_pt) == 0) return(NA)
  cand <- cr_pt[cr_pt$offset_h >= 0 & cr_pt$offset_h < t_h, ]
  if (nrow(cand) == 0) return(NA)
  cand$labresult[which.max(cand$offset_h)]  # last Cr before T0 = baseline
}
find_max_cr <- function(cr_pt, t_start, t_end) {
  if (is.null(cr_pt) || nrow(cr_pt) == 0) return(NA)
  cand <- cr_pt[cr_pt$offset_h >= t_start & cr_pt$offset_h <= t_end, ]
  if (nrow(cand) == 0) return(NA)
  max(cand$labresult, na.rm = TRUE)
}

# ── Baseline Mg = first postop Mg before T₀ (03c-consistent) ─────
cat("  Computing baseline Mg (first postop before T0)...\n")
mg_labs <- labs_all[labs_all$lab_name == "magnesium" & !is.na(labs_all$value), ]
mg_labs <- mg_labs[order(mg_labs$pid, mg_labs$offset_h), ]
mg_list <- split(mg_labs[, c("value", "offset_h")], mg_labs$pid)

mg_trt <- rep(NA_real_, n_pairs)
for (i in seq_len(n_pairs)) {
  pid <- as.character(pairs$trt_pid[i])
  t_mg <- pairs$t_mg[i]
  mg_i <- mg_list[[pid]]
  if (is.null(mg_i) || nrow(mg_i) == 0) next
  pre <- mg_i[mg_i$offset_h < t_mg, , drop = FALSE]
  if (nrow(pre) == 0) next
  mg_trt[i] <- pre$value[1]  # FIRST postop value (consistent with 03c)
}
cat(sprintf("  Mg available: %d/%d (%.1f%%)\n",
            sum(!is.na(mg_trt)), n_pairs, 100*sum(!is.na(mg_trt))/n_pairs))

# ── eGFR strata (treated patient's baseline, preserves pairs) ────
egfr_trt <- all_pts$egfr[trt_rows]
egfr_strat <- rep(NA_character_, n_pairs)
egfr_strat[!is.na(egfr_trt) & egfr_trt >= 90]                 <- "eGFR>=90"
egfr_strat[!is.na(egfr_trt) & egfr_trt >= 60 & egfr_trt < 90] <- "eGFR_60-89"
egfr_strat[!is.na(egfr_trt) & egfr_trt >= 45 & egfr_trt < 60] <- "eGFR_45-59"
egfr_strat[!is.na(egfr_trt) & egfr_trt >= 30 & egfr_trt < 45] <- "eGFR_30-44"
egfr_strat[!is.na(egfr_trt) & egfr_trt < 30]                  <- "eGFR<30"

cat(sprintf("  eGFR strata: %s\n",
            paste(names(table(egfr_strat)), table(egfr_strat), sep="=", collapse=", ")))

# ── Mg strata ─────────────────────────────────────────────────────
mg_strat <- rep(NA_character_, n_pairs)
mg_strat[!is.na(mg_trt) & mg_trt < 1.6]                  <- "Mg<1.6"
mg_strat[!is.na(mg_trt) & mg_trt >= 1.6 & mg_trt < 1.8]  <- "Mg_1.6-1.8"
mg_strat[!is.na(mg_trt) & mg_trt >= 1.8 & mg_trt < 2.0]  <- "Mg_1.8-2.0"
mg_strat[!is.na(mg_trt) & mg_trt >= 2.0 & mg_trt < 2.3]  <- "Mg_2.0-2.3"
mg_strat[!is.na(mg_trt) & mg_trt >= 2.3]                  <- "Mg>=2.3"

cat(sprintf("  Mg strata: %s\n",
            paste(names(table(mg_strat)), table(mg_strat), sep="=", collapse=", ")))

# ── Compute pair-level peak ΔCr at 9 horizons ────────────────────
cat("  Computing pair-level peak ΔCr (9 horizons incl. 7d)...\n")
dcr_trt <- dcr_ctl <- matrix(NA_real_, nrow = n_pairs, ncol = length(TARGETS))
colnames(dcr_trt) <- colnames(dcr_ctl) <- paste0("h", TARGETS)

for (i in seq_len(n_pairs)) {
  t_mg_i <- pairs$t_mg[i]
  tp <- as.character(pairs$trt_pid[i])
  cp <- as.character(pairs$ctl_pid[i])
  pre_t <- find_cr_pre(cr_list[[tp]], t_mg_i)
  pre_c <- find_cr_pre(cr_list[[cp]], t_mg_i)
  if (is.na(pre_t) || is.na(pre_c)) next
  for (j in seq_along(TARGETS)) {
    max_t <- find_max_cr(cr_list[[tp]], t_mg_i, t_mg_i + TARGETS[j])
    max_c <- find_max_cr(cr_list[[cp]], t_mg_i, t_mg_i + TARGETS[j])
    if (!is.na(max_t)) dcr_trt[i, j] <- max_t - pre_t
    if (!is.na(max_c)) dcr_ctl[i, j] <- max_c - pre_c
  }
}

# ── Generic DiD runner ────────────────────────────────────────────
run_did_block <- function(strat_var_name, strat_vec, strat_order) {
  tlabels <- ifelse(TARGETS < 168, paste0(TARGETS, "h"), "7d")
  cat(sprintf("\n── %s ──\n", strat_var_name))
  cat(sprintf("  %-15s  ", "Stratum"))
  for (tl in tlabels) cat(sprintf(" %7s", tl))
  cat("\n", paste(rep("-", 90), collapse = ""), "\n")

  block <- list()
  for (stg in strat_order) {
    if (stg == "Overall") idx <- seq_len(n_pairs)
    else idx <- which(strat_vec == stg)

    row_str <- sprintf("  %-15s  ", stg)
    for (j in seq_along(TARGETS)) {
      v <- idx[!is.na(dcr_trt[idx, j]) & !is.na(dcr_ctl[idx, j])]
      if (length(v) < 30) {
        block[[length(block) + 1]] <- data.frame(
          db = db, strat_var = strat_var_name, stratum = stg,
          target_h = TARGETS[j], did = NA, se = NA, p = NA,
          ci_lo = NA, ci_hi = NA, n = length(v))
        row_str <- paste0(row_str, "      -- ")
        next
      }
      df <- data.frame(
        delta_cr = c(dcr_trt[v, j], dcr_ctl[v, j]),
        treated = rep(c(1, 0), each = length(v))
      )
      fit <- lm(delta_cr ~ treated, data = df)
      ct <- tryCatch(
        coeftest(fit, vcov. = vcovHC(fit, type = "HC1")),
        error = function(e) coeftest(fit)
      )
      est <- ct["treated", 1]; se <- ct["treated", 2]; p <- ct["treated", 4]
      block[[length(block) + 1]] <- data.frame(
        db = db, strat_var = strat_var_name, stratum = stg,
        target_h = TARGETS[j],
        did = round(est, 5), se = round(se, 5), p = round(p, 5),
        ci_lo = round(est - 1.96 * se, 5), ci_hi = round(est + 1.96 * se, 5),
        n = length(v))
      sig <- if (!is.na(p) && p < 0.05) "*" else " "
      row_str <- paste0(row_str, sprintf(" %+.4f%s", est, sig))
    }
    cat(row_str, "\n")
  }
  do.call(rbind, block)
}

# ── Run all three stratifications ─────────────────────────────────
res_overall <- run_did_block("Overall", rep("Overall", n_pairs),
                              "Overall")

res_egfr <- run_did_block("eGFR", egfr_strat,
                           c("Overall", "eGFR>=90", "eGFR_60-89",
                             "eGFR_45-59", "eGFR_30-44", "eGFR<30"))

res_mg <- run_did_block("Mg", mg_strat,
                         c("Overall", "Mg<1.6", "Mg_1.6-1.8",
                           "Mg_1.8-2.0", "Mg_2.0-2.3", "Mg>=2.3"))

# ── Save ──────────────────────────────────────────────────────────
res_df <- rbind(res_overall, res_egfr, res_mg)
outpath <- file.path(RESULTS, sprintf("did_stratified_%s.csv", tag))
write.csv(res_df, outpath, row.names = FALSE)
cat(sprintf("\n  Saved: %s (%d rows)\n", outpath, nrow(res_df)))

# Backward-compat: also save Mg-only slice as mg_did_{db}.csv
mg_slice <- res_df[res_df$strat_var == "Mg", ]
names(mg_slice)[names(mg_slice) == "stratum"] <- "mg_strat"
mg_slice$strat_var <- NULL
write.csv(mg_slice, file.path(RESULTS, sprintf("mg_did_%s.csv", tag)), row.names = FALSE)
cat(sprintf("  Saved: mg_did_%s.csv (backward compat, %d rows)\n", tag, nrow(mg_slice)))

cat(sprintf("\n%s\n03d_mg_did.R — %s DONE\n%s\n", SEP, db, SEP))
