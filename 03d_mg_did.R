#!/usr/bin/env Rscript
# ============================================================================
# 03d_mg_did.R — Mg-stratified DiD on continuous ΔCr (max-Cr version)
#
# KEY CHANGE: post-treatment Cr = MAX Cr within [t_mg, t_mg + target_h],
# not the closest Cr to target_h ± 12h. This matches how AKI is defined
# (peak rise triggers KDIGO), making DiD and AKI comparable on the same
# scale. (Dr. Su's suggestion: "如果用0-36的最高CRT呢？")
#
# Uses existing primary pairs — no new matching.
# Baseline Mg from did_labs_all (same source as 03c_mg_strat.R).
#
# Outputs:
#   results/mg_did_{db}.csv — DiD by Mg stratum × time point
#
# Usage: Rscript 03d_mg_did.R mimic
#        Rscript 03d_mg_did.R eicu
# ============================================================================

suppressPackageStartupMessages({ library(sandwich); library(lmtest) })

RESULTS   <- path.expand("~/mg_aki/results")
TARGETS   <- c(6, 12, 18, 24, 30, 36, 42, 48)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) { cat("Usage: Rscript 03d_mg_did.R <db>\n"); quit(status = 1) }
tag <- tolower(args[1]); db <- toupper(tag)

SEP <- paste(rep("=", 60), collapse = "")
cat(sprintf("\n%s\n03d_mg_did.R — %s\n  Mg-stratified DiD on continuous ΔCr\n%s\n", SEP, db, SEP))

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

# ── Cr helpers (exact 02_psm.R logic) ─────────────────────────────
find_cr_pre <- function(cr_pt, t_h) {
  if (is.null(cr_pt) || nrow(cr_pt) == 0) return(NA)
  cand <- cr_pt[cr_pt$offset_h >= 0 & cr_pt$offset_h < t_h, ]
  if (nrow(cand) == 0) return(NA)
  cand$labresult[which.max(cand$offset_h)]
}

find_cr_post_max <- function(cr_pt, t_mg, target_h) {
  # MAX Cr from t_mg to t_mg + target_h (matches AKI scan logic)
  if (is.null(cr_pt) || nrow(cr_pt) == 0) return(NA)
  cand <- cr_pt[cr_pt$offset_h >= t_mg & cr_pt$offset_h <= (t_mg + target_h), ]
  if (nrow(cand) == 0) return(NA)
  max(cand$labresult, na.rm = TRUE)
}

# ── Baseline Mg from did_labs_all (same as 03c) ──────────────────
cat("  Computing baseline Mg...\n")
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

mg_avail <- sum(!is.na(mg_trt))
cat(sprintf("  Mg available: %d/%d (%.1f%%)\n", mg_avail, n_pairs, 100*mg_avail/n_pairs))

# ── Define strata ─────────────────────────────────────────────────
mg_strat <- rep(NA_character_, n_pairs)
mg_strat[!is.na(mg_trt) & mg_trt < 1.6]                  <- "Mg<1.6"
mg_strat[!is.na(mg_trt) & mg_trt >= 1.6 & mg_trt < 1.8]  <- "Mg_1.6-1.8"
mg_strat[!is.na(mg_trt) & mg_trt >= 1.8 & mg_trt < 2.0]  <- "Mg_1.8-2.0"
mg_strat[!is.na(mg_trt) & mg_trt >= 2.0 & mg_trt < 2.3]  <- "Mg_2.0-2.3"
mg_strat[!is.na(mg_trt) & mg_trt >= 2.3]                  <- "Mg>=2.3"

cat(sprintf("  Mg strata: %s\n",
            paste(names(table(mg_strat)), table(mg_strat), sep="=", collapse=", ")))

# ── Compute pair-level ΔCr at each time point ─────────────────────
cat("  Computing pair-level ΔCr at 8 time points...\n")

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
    post_t <- find_cr_post_max(cr_list[[tp]], t_mg_i, TARGETS[j])
    post_c <- find_cr_post_max(cr_list[[cp]], t_mg_i, TARGETS[j])
    if (!is.na(post_t)) dcr_trt[i, j] <- post_t - pre_t
    if (!is.na(post_c)) dcr_ctl[i, j] <- post_c - pre_c
  }
}

# ── DiD within each stratum × time point ──────────────────────────
cat(sprintf("\n  %-15s  ", "Stratum"))
for (h in TARGETS) cat(sprintf(" %6dh", h))
cat("\n", paste(rep("-", 80), collapse = ""), "\n")

mg_order <- c("Overall", "Mg<1.6", "Mg_1.6-1.8", "Mg_1.8-2.0", "Mg_2.0-2.3", "Mg>=2.3")
results <- list()

for (stg in mg_order) {
  if (stg == "Overall") idx <- seq_len(n_pairs)
  else idx <- which(mg_strat == stg)

  row_str <- sprintf("  %-15s  ", stg)
  for (j in seq_along(TARGETS)) {
    v <- idx[!is.na(dcr_trt[idx, j]) & !is.na(dcr_ctl[idx, j])]
    if (length(v) < 30) {
      results[[length(results) + 1]] <- data.frame(
        db = db, mg_strat = stg, target_h = TARGETS[j],
        did = NA, se = NA, p = NA, ci_lo = NA, ci_hi = NA, n = length(v))
      row_str <- paste0(row_str, "     -- ")
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
    est <- ct["treated", 1]
    se  <- ct["treated", 2]
    p   <- ct["treated", 4]
    results[[length(results) + 1]] <- data.frame(
      db = db, mg_strat = stg, target_h = TARGETS[j],
      did = round(est, 5), se = round(se, 5), p = round(p, 5),
      ci_lo = round(est - 1.96 * se, 5), ci_hi = round(est + 1.96 * se, 5),
      n = length(v))

    sig <- if (!is.na(p) && p < 0.05) "*" else " "
    row_str <- paste0(row_str, sprintf(" %+.3f%s", est, sig))
  }
  cat(row_str, "\n")
}

# ── Save ──────────────────────────────────────────────────────────
res_df <- do.call(rbind, results)
outpath <- file.path(RESULTS, sprintf("mg_did_%s.csv", tag))
write.csv(res_df, outpath, row.names = FALSE)
cat(sprintf("\n  Saved: %s (%d rows)\n", outpath, nrow(res_df)))

cat(sprintf("\n%s\n03d_mg_did.R — %s DONE\n%s\n", SEP, db, SEP))
