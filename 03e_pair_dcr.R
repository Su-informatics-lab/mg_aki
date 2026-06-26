#!/usr/bin/env Rscript
# ============================================================================
# 03e_pair_dcr.R — Save pair-level ΔCr for distribution plots
#
# Outputs pair_dcr_{db}.csv with per-pair:
#   dcr_36h_trt/ctl  — ΔCr at 36h (same as DiD estimand)
#   max_dcr_7d_trt/ctl — max ΔCr over 7d (same as AKI scan)
#   mg_strat — treated patient's baseline Mg stratum
#
# Usage: Rscript 03e_pair_dcr.R mimic
#        Rscript 03e_pair_dcr.R eicu
# ============================================================================

RESULTS <- path.expand("~/mg_aki/results")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) { cat("Usage: Rscript 03e_pair_dcr.R <db>\n"); quit(status = 1) }
tag <- tolower(args[1]); db <- toupper(tag)

cat(sprintf("\n03e_pair_dcr.R — %s\n", db))

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
n_pairs <- nrow(pairs)

# ── Cr helpers ────────────────────────────────────────────────────
find_cr_pre <- function(cr_pt, t_h) {
  if (is.null(cr_pt) || nrow(cr_pt) == 0) return(NA)
  cand <- cr_pt[cr_pt$offset_h >= 0 & cr_pt$offset_h < t_h, ]
  if (nrow(cand) == 0) return(NA)
  cand$labresult[which.max(cand$offset_h)]
}

# ── Baseline Mg ───────────────────────────────────────────────────
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
  mg_trt[i] <- pre$value[nrow(pre)]
}

mg_strat <- rep(NA_character_, n_pairs)
mg_strat[!is.na(mg_trt) & mg_trt < 1.6]                  <- "Mg<1.6"
mg_strat[!is.na(mg_trt) & mg_trt >= 1.6 & mg_trt < 1.8]  <- "Mg_1.6-1.8"
mg_strat[!is.na(mg_trt) & mg_trt >= 1.8 & mg_trt < 2.0]  <- "Mg_1.8-2.0"
mg_strat[!is.na(mg_trt) & mg_trt >= 2.0 & mg_trt < 2.3]  <- "Mg_2.0-2.3"
mg_strat[!is.na(mg_trt) & mg_trt >= 2.3]                  <- "Mg>=2.3"

# ── Compute pair-level ΔCr ────────────────────────────────────────
cat("  Computing pair-level ΔCr...\n")
dcr_36h_trt <- dcr_36h_ctl <- rep(NA_real_, n_pairs)
max_dcr_7d_trt <- max_dcr_7d_ctl <- rep(NA_real_, n_pairs)

for (i in seq_len(n_pairs)) {
  t_mg_i <- pairs$t_mg[i]
  tp <- as.character(pairs$trt_pid[i])
  cp <- as.character(pairs$ctl_pid[i])

  pre_t <- find_cr_pre(cr_list[[tp]], t_mg_i)
  pre_c <- find_cr_pre(cr_list[[cp]], t_mg_i)

  # Max ΔCr within 0–36h (matches revised DiD estimand)
  if (!is.na(pre_t) && !is.null(cr_t)) {
    post36 <- cr_t[cr_t$offset_h >= t_mg_i & cr_t$offset_h <= (t_mg_i + 36), ]
    if (nrow(post36) > 0) dcr_36h_trt[i] <- max(post36$labresult, na.rm = TRUE) - pre_t
  }
  if (!is.na(pre_c) && !is.null(cr_c)) {
    post36 <- cr_c[cr_c$offset_h >= t_mg_i & cr_c$offset_h <= (t_mg_i + 36), ]
    if (nrow(post36) > 0) dcr_36h_ctl[i] <- max(post36$labresult, na.rm = TRUE) - pre_c
  }

  # Max ΔCr over 7d (AKI-relevant)
  cr_t <- cr_list[[tp]]; cr_c <- cr_list[[cp]]
  if (!is.na(pre_t) && !is.null(cr_t)) {
    post <- cr_t[cr_t$offset_h >= t_mg_i & cr_t$offset_h <= (t_mg_i + 168), ]
    if (nrow(post) > 0) max_dcr_7d_trt[i] <- max(post$labresult - pre_t, na.rm = TRUE)
  }
  if (!is.na(pre_c) && !is.null(cr_c)) {
    post <- cr_c[cr_c$offset_h >= t_mg_i & cr_c$offset_h <= (t_mg_i + 168), ]
    if (nrow(post) > 0) max_dcr_7d_ctl[i] <- max(post$labresult - pre_c, na.rm = TRUE)
  }
}

# ── Save ──────────────────────────────────────────────────────────
out <- data.frame(
  trt_pid     = pairs$trt_pid,
  ctl_pid     = pairs$ctl_pid,
  mg_strat    = mg_strat,
  mg_value    = round(mg_trt, 2),
  dcr_36h_trt = round(dcr_36h_trt, 4),
  dcr_36h_ctl = round(dcr_36h_ctl, 4),
  max_dcr_7d_trt = round(max_dcr_7d_trt, 4),
  max_dcr_7d_ctl = round(max_dcr_7d_ctl, 4),
  stringsAsFactors = FALSE
)

outpath <- file.path(RESULTS, sprintf("pair_dcr_%s.csv", tag))
write.csv(out, outpath, row.names = FALSE)

v36 <- sum(!is.na(dcr_36h_trt) & !is.na(dcr_36h_ctl))
v7d <- sum(!is.na(max_dcr_7d_trt) & !is.na(max_dcr_7d_ctl))
cat(sprintf("  Saved: %s (%d pairs, 36h valid=%d, 7d valid=%d)\n", outpath, n_pairs, v36, v7d))
