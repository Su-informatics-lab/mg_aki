#!/usr/bin/env Rscript
# ============================================================================
# 05_qc.R — Quantitative consistency checks for Mg → AKI study
#
# Addresses two reviewer-anticipating concerns:
#
#   A. Simpson's Paradox: Why does overall DiD ΔCr ≈ 0.01 mg/dL coexist
#      with AKI OR = 0.87?  Answer: eGFR-stratified ΔCr shows protection
#      (eGFR≥90, −0.07) and harm (eGFR<45, +0.34) that cancel on average,
#      but AKI's binary threshold lets the larger protective group dominate.
#
#   B. Surveillance Bias: Do treated patients get more Cr draws, inflating
#      AKI detection?  Answer: within pair-matched common observation
#      windows, monitoring intensity is balanced across eGFR strata.
#
# Runs on existing matched pairs from 02_psm.R.  No new matching.
#
# Outputs:
#   results/qc_simpson_{db}.csv     — eGFR-stratified ΔCr vs AKI consistency
#   results/qc_surveillance_{db}.csv — monitoring intensity by eGFR stratum
#   results/qc_kdigo_{db}.csv       — KDIGO arm decomposition
#   results/qc_threshold_{db}.csv   — ΔCr distribution around 0.3
#
# Usage: Rscript 05_qc.R mimic
#        Rscript 05_qc.R eicu
# ============================================================================

suppressPackageStartupMessages({ library(sandwich); library(lmtest) })

RESULTS <- path.expand("~/mg_aki/results")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) { cat("Usage: Rscript 05_qc.R <db>\n"); quit(status = 1) }
tag <- tolower(args[1]); db <- toupper(tag)

SEP  <- paste(rep("=", 60), collapse = "")
SEP2 <- paste(rep("-", 60), collapse = "")

cat(sprintf("\n%s\n05_qc.R — %s — Quantitative Consistency Checks\n%s\n", SEP, db, SEP))

# ═══════════════════════════════════════════════════════════════════
# LOAD DATA
# ═══════════════════════════════════════════════════════════════════
all_pts  <- read.csv(file.path(RESULTS, sprintf("did_all_%s.csv", tag)), stringsAsFactors = FALSE)
pairs    <- read.csv(file.path(RESULTS, sprintf("did_pairs_primary_yet_untreated_%s.csv", tag)),
                      stringsAsFactors = FALSE)
cr_all   <- read.csv(file.path(RESULTS, sprintf("did_cr_all_%s.csv", tag)), stringsAsFactors = FALSE)
labs_all <- read.csv(file.path(RESULTS, sprintf("did_labs_all_%s.csv", tag)), stringsAsFactors = FALSE)

# Standardize pid columns
cr_id <- if ("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
cr_all$pid <- cr_all[[cr_id]]
if (!"offset_h" %in% names(cr_all)) cr_all$offset_h <- cr_all$labresultoffset / 60
cr_all <- cr_all[order(cr_all$pid, cr_all$offset_h), ]

labs_pid <- if ("patientunitstayid" %in% names(labs_all)) "patientunitstayid" else "stay_id"
labs_all$pid <- labs_all[[labs_pid]]

# Pre-split for speed
cr_list   <- split(cr_all[, c("labresult", "offset_h")], cr_all$pid)
labs_split <- split(labs_all[, c("pid", "lab_name", "offset_h")], labs_all$pid)

trt_rows <- match(pairs$trt_pid, all_pts$pid)
ctl_rows <- match(pairs$ctl_pid, all_pts$pid)
n_pairs  <- nrow(pairs)

egfr_trt <- all_pts$egfr[trt_rows]
los_trt  <- all_pts$icu_discharge_h[trt_rows]
los_ctl  <- all_pts$icu_discharge_h[ctl_rows]

cat(sprintf("  Pairs: %d  |  Labs: %d rows (%s)\n",
            n_pairs, nrow(labs_all),
            paste(sort(unique(labs_all$lab_name)), collapse = ", ")))

# ── eGFR bins ─────────────────────────────────────────────────────
egfr_bins <- list(
  "eGFR>=90"   = !is.na(egfr_trt) & egfr_trt >= 90,
  "eGFR 60-89" = !is.na(egfr_trt) & egfr_trt >= 60 & egfr_trt < 90,
  "eGFR 45-59" = !is.na(egfr_trt) & egfr_trt >= 45 & egfr_trt < 60,
  "eGFR <45"   = !is.na(egfr_trt) & egfr_trt < 45
)

# ═══════════════════════════════════════════════════════════════════
# COMPUTE PER-PATIENT 7d Cr TRAJECTORY
# ═══════════════════════════════════════════════════════════════════
cat("  Computing 7d Cr trajectories...\n")

compute_traj <- function(pid, t_mg) {
  out <- list(bl = NA, max_dcr_48h = NA, max_dcr_7d = NA, max_ratio_7d = NA,
              aki_absolute = NA, aki_ratio = NA, aki_any = NA, n_post_cr = 0)
  cr <- cr_list[[as.character(pid)]]
  if (is.null(cr) || nrow(cr) < 1) return(out)
  pre <- cr[cr$offset_h >= 0 & cr$offset_h < t_mg, ]
  if (nrow(pre) == 0) return(out)
  bl <- pre$labresult[which.max(pre$offset_h)]
  if (is.na(bl) || bl <= 0) return(out)
  out$bl <- bl
  post <- cr[cr$offset_h >= t_mg & cr$offset_h <= (t_mg + 168), ]
  out$n_post_cr <- nrow(post)
  if (nrow(post) == 0) {
    out$max_dcr_48h <- 0; out$max_dcr_7d <- 0; out$max_ratio_7d <- 1
    out$aki_absolute <- 0; out$aki_ratio <- 0; out$aki_any <- 0
    return(out)
  }
  dcr <- post$labresult - bl
  ratio <- post$labresult / bl
  hrs <- post$offset_h - t_mg
  in48 <- hrs <= 48
  out$max_dcr_7d   <- max(dcr, na.rm = TRUE)
  out$max_dcr_48h  <- if (any(in48)) max(dcr[in48], na.rm = TRUE) else 0
  out$max_ratio_7d <- max(ratio, na.rm = TRUE)
  out$aki_absolute <- as.integer(any(in48 & dcr >= 0.3, na.rm = TRUE))
  out$aki_ratio    <- as.integer(any(ratio >= 1.5, na.rm = TRUE))
  out$aki_any      <- as.integer(out$aki_absolute == 1 | out$aki_ratio == 1)
  return(out)
}

trt_traj <- ctl_traj <- vector("list", n_pairs)
for (i in seq_len(n_pairs)) {
  trt_traj[[i]] <- compute_traj(pairs$trt_pid[i], pairs$t_mg[i])
  ctl_traj[[i]] <- compute_traj(pairs$ctl_pid[i], pairs$t_mg[i])
}

vec <- function(lst, f) sapply(lst, function(x) x[[f]])
dcr7d_t <- vec(trt_traj, "max_dcr_7d");  dcr7d_c <- vec(ctl_traj, "max_dcr_7d")
dcr48_t <- vec(trt_traj, "max_dcr_48h"); dcr48_c <- vec(ctl_traj, "max_dcr_48h")
aki_abs_t <- vec(trt_traj, "aki_absolute"); aki_abs_c <- vec(ctl_traj, "aki_absolute")
aki_rat_t <- vec(trt_traj, "aki_ratio");    aki_rat_c <- vec(ctl_traj, "aki_ratio")
aki_any_t <- vec(trt_traj, "aki_any");      aki_any_c <- vec(ctl_traj, "aki_any")
npost_cr_t <- vec(trt_traj, "n_post_cr");   npost_cr_c <- vec(ctl_traj, "n_post_cr")

valid <- !is.na(aki_any_t) & !is.na(aki_any_c)

# ═══════════════════════════════════════════════════════════════════
# A. SIMPSON'S PARADOX DEMONSTRATION
# ═══════════════════════════════════════════════════════════════════
cat(sprintf("\n%s\n  A. SIMPSON'S PARADOX: Continuous ΔCr vs Binary AKI\n%s\n", SEP, SEP))

# ── A1: eGFR-stratified ΔCr vs AKI ───────────────────────────────
cat(sprintf("\n  A1. eGFR-stratified consistency (max ΔCr 7d vs AKI 7d)\n%s\n", SEP2))
cat(sprintf("  %-12s  %7s  %7s  %8s  | %6s  %6s  %8s  %5s\n",
            "Stratum", "dCr_trt", "dCr_ctl", "dCr_diff", "AKI_t%", "AKI_c%", "AKI_diff", "n"))
cat(paste(rep("-", 78), collapse = ""), "\n")

simpson_rows <- list()
for (nm in c("Overall", names(egfr_bins))) {
  idx <- if (nm == "Overall") which(valid)
         else which(valid & egfr_bins[[nm]])
  if (length(idx) < 30) next
  mt <- mean(dcr7d_t[idx], na.rm = TRUE)
  mc <- mean(dcr7d_c[idx], na.rm = TRUE)
  at <- 100 * mean(aki_any_t[idx])
  ac <- 100 * mean(aki_any_c[idx])
  cat(sprintf("  %-12s  %+7.4f  %+7.4f  %+8.4f  | %5.1f%%  %5.1f%%  %+7.1f pp  %5d\n",
              nm, mt, mc, mt - mc, at, ac, at - ac, length(idx)))
  simpson_rows[[length(simpson_rows) + 1]] <- data.frame(
    db = db, stratum = nm, n = length(idx),
    dcr7d_trt = round(mt, 4), dcr7d_ctl = round(mc, 4), dcr7d_diff = round(mt - mc, 4),
    aki_rate_trt = round(at, 1), aki_rate_ctl = round(ac, 1), aki_diff_pp = round(at - ac, 1),
    stringsAsFactors = FALSE)
}

cat("\n  Interpretation:\n")
cat("    Overall ΔCr ≈ 0 because eGFR≥90 protection (−) and eGFR<45 harm (+) cancel.\n")
cat("    Binary AKI threshold + larger eGFR≥90 group → overall OR < 1.\n")

# ── A2: Threshold density ─────────────────────────────────────────
cat(sprintf("\n  A2. ΔCr distribution around KDIGO threshold (48h max ΔCr)\n%s\n", SEP2))

bins <- list(c(0.0,0.1), c(0.1,0.2), c(0.2,0.3), c(0.3,0.4), c(0.4,0.5),
             c(0.5,1.0), c(1.0,Inf))
nv <- sum(valid)
thresh_rows <- list()

for (b in bins) {
  nt <- sum(dcr48_t[valid] >= b[1] & dcr48_t[valid] < b[2], na.rm = TRUE)
  nc <- sum(dcr48_c[valid] >= b[1] & dcr48_c[valid] < b[2], na.rm = TRUE)
  mark <- if (b[1] == 0.2 | b[1] == 0.3) " ◀" else ""
  cat(sprintf("  ΔCr [%.1f, %.1f):  trt=%d (%.1f%%)  ctl=%d (%.1f%%)  diff=%+d%s\n",
              b[1], b[2], nt, 100*nt/nv, nc, 100*nc/nv, nt - nc, mark))
  thresh_rows[[length(thresh_rows) + 1]] <- data.frame(
    db = db, bin_lo = b[1], bin_hi = b[2],
    n_trt = nt, pct_trt = round(100*nt/nv, 1),
    n_ctl = nc, pct_ctl = round(100*nc/nv, 1),
    diff = nt - nc, stringsAsFactors = FALSE)
}

# Marginal AKI
near_t <- sum(dcr48_t[valid] >= 0.30 & dcr48_t[valid] < 0.35, na.rm = TRUE)
near_c <- sum(dcr48_c[valid] >= 0.30 & dcr48_c[valid] < 0.35, na.rm = TRUE)
cat(sprintf("\n  Marginal AKI (ΔCr 0.30–0.35): trt=%d, ctl=%d\n", near_t, near_c))
cat(sprintf("  As %% of absolute-arm AKI: trt=%.0f%%, ctl=%.0f%%\n",
            100*near_t/max(sum(aki_abs_t[valid]),1),
            100*near_c/max(sum(aki_abs_c[valid]),1)))

# ═══════════════════════════════════════════════════════════════════
# B. KDIGO ARM DECOMPOSITION
# ═══════════════════════════════════════════════════════════════════
cat(sprintf("\n%s\n  B. KDIGO ARM DECOMPOSITION\n%s\n", SEP, SEP))
cat("  AKI = absolute (ΔCr≥0.3 within 48h) OR ratio (≥1.5× within 7d)\n\n")

both_t <- aki_abs_t[valid] == 1 & aki_rat_t[valid] == 1
both_c <- aki_abs_c[valid] == 1 & aki_rat_c[valid] == 1
abs_only_t <- aki_abs_t[valid] == 1 & aki_rat_t[valid] == 0
abs_only_c <- aki_abs_c[valid] == 1 & aki_rat_c[valid] == 0
rat_only_t <- aki_abs_t[valid] == 0 & aki_rat_t[valid] == 1
rat_only_c <- aki_abs_c[valid] == 0 & aki_rat_c[valid] == 1

kdigo_rows <- list()
for (lab_row in list(
  list("AKI_any",           aki_any_t[valid],  aki_any_c[valid]),
  list("absolute_only_48h", abs_only_t,        abs_only_c),
  list("ratio_only_7d",     rat_only_t,        rat_only_c),
  list("both_arms",         both_t,            both_c)
)) {
  nt <- sum(lab_row[[2]]); nc <- sum(lab_row[[3]])
  cat(sprintf("  %-25s  trt=%d (%.1f%%)  ctl=%d (%.1f%%)  diff=%+d\n",
              lab_row[[1]], nt, 100*nt/nv, nc, 100*nc/nv, nt - nc))
  kdigo_rows[[length(kdigo_rows) + 1]] <- data.frame(
    db = db, arm = lab_row[[1]],
    n_trt = nt, pct_trt = round(100*nt/nv, 1),
    n_ctl = nc, pct_ctl = round(100*nc/nv, 1),
    diff = nt - nc, stringsAsFactors = FALSE)
}

# ═══════════════════════════════════════════════════════════════════
# C. SURVEILLANCE BIAS
# ═══════════════════════════════════════════════════════════════════
cat(sprintf("\n%s\n  C. SURVEILLANCE BIAS\n%s\n", SEP, SEP))

# ── C1: ICU LOS ──────────────────────────────────────────────────
cat(sprintf("  C1. ICU length of stay (hours)\n%s\n", SEP2))

v_los <- !is.na(los_trt) & !is.na(los_ctl)
cat(sprintf("  %-12s  %7s  %7s  %5s  %5s\n", "Stratum", "LOS_trt", "LOS_ctl", "Ratio", "n"))
cat(paste(rep("-", 48), collapse = ""), "\n")
for (nm in c("Overall", names(egfr_bins))) {
  idx <- if (nm == "Overall") which(v_los) else which(v_los & egfr_bins[[nm]])
  if (length(idx) < 30) next
  mt <- mean(los_trt[idx]); mc <- mean(los_ctl[idx])
  cat(sprintf("  %-12s  %7.1f  %7.1f  %5.2f  %5d\n", nm, mt, mc, mt/mc, length(idx)))
}

# ── C2: Pair-matched common-window Cr rate ────────────────────────
cat(sprintf("\n  C2. Pair-matched Cr rate (common observation window)\n%s\n", SEP2))
cat("  For each pair: window = min(LOS_trt, LOS_ctl, t_mg+168) − t_mg\n")
cat("  Pairs with window < 24h excluded (cannot assess monitoring rate)\n\n")

remaining_trt <- los_trt - pairs$t_mg
remaining_ctl <- los_ctl - pairs$t_mg
common_window <- pmin(remaining_trt, remaining_ctl, 168, na.rm = TRUE)
common_valid  <- !is.na(common_window) & common_window >= 24

# Count Cr within common window per pair
cr_cw_trt <- cr_cw_ctl <- numeric(n_pairs)
for (i in seq_len(n_pairs)) {
  if (!common_valid[i]) next
  t_mg_i <- pairs$t_mg[i]; cw <- common_window[i]
  cr_t <- cr_list[[as.character(pairs$trt_pid[i])]]
  cr_c <- cr_list[[as.character(pairs$ctl_pid[i])]]
  if (!is.null(cr_t))
    cr_cw_trt[i] <- sum(cr_t$offset_h >= t_mg_i & cr_t$offset_h <= (t_mg_i + cw))
  if (!is.null(cr_c))
    cr_cw_ctl[i] <- sum(cr_c$offset_h >= t_mg_i & cr_c$offset_h <= (t_mg_i + cw))
}

cr_rate_trt <- cr_cw_trt / (common_window / 24)
cr_rate_ctl <- cr_cw_ctl / (common_window / 24)

cat(sprintf("  %-12s  %9s  %9s  %5s  %5s  %s\n",
            "Stratum", "Cr/24h_t", "Cr/24h_c", "Ratio", "n", "Interpretation"))
cat(paste(rep("-", 72), collapse = ""), "\n")

surv_rows <- list()
for (nm in c("Overall", names(egfr_bins))) {
  idx <- if (nm == "Overall") which(common_valid)
         else which(common_valid & egfr_bins[[nm]])
  if (length(idx) < 30) next
  rt <- mean(cr_rate_trt[idx], na.rm = TRUE)
  rc <- mean(cr_rate_ctl[idx], na.rm = TRUE)
  ratio <- rt / max(rc, 0.001)
  interp <- if (abs(ratio - 1) < 0.10) "balanced"
            else if (ratio < 1) "trt LESS monitored"
            else "trt MORE monitored"
  cat(sprintf("  %-12s  %9.2f  %9.2f  %5.2f  %5d  %s\n",
              nm, rt, rc, ratio, length(idx), interp))
  surv_rows[[length(surv_rows) + 1]] <- data.frame(
    db = db, stratum = nm, measure = "cr_per_24h_common_window",
    trt = round(rt, 2), ctl = round(rc, 2), ratio = round(ratio, 2),
    n = length(idx), stringsAsFactors = FALSE)
}

# ── C3: All-lab rate in common window ─────────────────────────────
cat(sprintf("\n  C3. All-lab monitoring in common window\n%s\n", SEP2))

# Count all labs within common window
all_cw_trt <- all_cw_ctl <- numeric(n_pairs)
for (i in seq_len(n_pairs)) {
  if (!common_valid[i]) next
  t_mg_i <- pairs$t_mg[i]; cw <- common_window[i]
  lb_t <- labs_split[[as.character(pairs$trt_pid[i])]]
  lb_c <- labs_split[[as.character(pairs$ctl_pid[i])]]
  if (!is.null(lb_t))
    all_cw_trt[i] <- sum(lb_t$offset_h >= t_mg_i & lb_t$offset_h <= (t_mg_i + cw))
  if (!is.null(lb_c))
    all_cw_ctl[i] <- sum(lb_c$offset_h >= t_mg_i & lb_c$offset_h <= (t_mg_i + cw))
}

all_rate_trt <- all_cw_trt / (common_window / 24)
all_rate_ctl <- all_cw_ctl / (common_window / 24)

cat(sprintf("  %-12s  %9s  %9s  %5s  %5s\n",
            "Stratum", "Lab/24h_t", "Lab/24h_c", "Ratio", "n"))
cat(paste(rep("-", 48), collapse = ""), "\n")
for (nm in c("Overall", names(egfr_bins))) {
  idx <- if (nm == "Overall") which(common_valid)
         else which(common_valid & egfr_bins[[nm]])
  if (length(idx) < 30) next
  rt <- mean(all_rate_trt[idx], na.rm = TRUE)
  rc <- mean(all_rate_ctl[idx], na.rm = TRUE)
  cat(sprintf("  %-12s  %9.2f  %9.2f  %5.2f  %5d\n",
              nm, rt, rc, rt/max(rc, 0.001), length(idx)))
  surv_rows[[length(surv_rows) + 1]] <- data.frame(
    db = db, stratum = nm, measure = "all_lab_per_24h_common_window",
    trt = round(rt, 2), ctl = round(rc, 2), ratio = round(rt/max(rc, 0.001), 2),
    n = length(idx), stringsAsFactors = FALSE)
}

# ── C4: Cr vs K⁺ correlation (same panel draw?) ──────────────────
cat(sprintf("\n  C4. Cr vs K⁺ count correlation\n%s\n", SEP2))

# Count K within common window
k_cw_trt <- k_cw_ctl <- numeric(n_pairs)
for (i in seq_len(n_pairs)) {
  if (!common_valid[i]) next
  t_mg_i <- pairs$t_mg[i]; cw <- common_window[i]
  lb_t <- labs_split[[as.character(pairs$trt_pid[i])]]
  lb_c <- labs_split[[as.character(pairs$ctl_pid[i])]]
  if (!is.null(lb_t)) {
    sub <- lb_t[lb_t$lab_name == "potassium" & lb_t$offset_h >= t_mg_i & lb_t$offset_h <= (t_mg_i + cw), ]
    k_cw_trt[i] <- nrow(sub)
  }
  if (!is.null(lb_c)) {
    sub <- lb_c[lb_c$lab_name == "potassium" & lb_c$offset_h >= t_mg_i & lb_c$offset_h <= (t_mg_i + cw), ]
    k_cw_ctl[i] <- nrow(sub)
  }
}

all_cr_counts <- c(cr_cw_trt[common_valid], cr_cw_ctl[common_valid])
all_k_counts  <- c(k_cw_trt[common_valid], k_cw_ctl[common_valid])
vv <- all_cr_counts > 0 & all_k_counts > 0
if (sum(vv) > 100) {
  r <- cor(all_cr_counts[vv], all_k_counts[vv])
  cat(sprintf("  Pearson r(Cr, K⁺ draws) = %.3f  (n=%d)\n", r, sum(vv)))
  cat(sprintf("  Median Cr/K ratio       = %.2f\n", median(all_cr_counts[vv] / all_k_counts[vv])))
  if (r > 0.9) cat("  → Same-panel draws confirmed\n")
  else if (r > 0.7) cat("  → Strongly correlated (mostly same panels)\n")
}

# ═══════════════════════════════════════════════════════════════════
# D. SUMMARY
# ═══════════════════════════════════════════════════════════════════
cat(sprintf("\n%s\n  D. SUMMARY: Surveillance bias by eGFR stratum\n%s\n", SEP, SEP))
cat(sprintf("  %-12s  %7s  %10s  %12s  %s\n",
            "Stratum", "Cr_rate", "AKI_effect", "Bias_dir", "Conclusion"))
cat(paste(rep("-", 72), collapse = ""), "\n")

for (nm in names(egfr_bins)) {
  idx <- which(common_valid & egfr_bins[[nm]])
  if (length(idx) < 30) next
  rt <- mean(cr_rate_trt[idx], na.rm = TRUE)
  rc <- mean(cr_rate_ctl[idx], na.rm = TRUE)
  ratio <- rt / max(rc, 0.001)
  aki_dir <- if (nm %in% c("eGFR>=90", "eGFR 60-89")) "protective" else "harmful"

  if (ratio < 0.90) {
    bias <- "under-detect trt AKI"
    concl <- if (aki_dir == "protective") "CONSERVATIVE" else "attenuates harm"
  } else if (ratio > 1.10) {
    bias <- "over-detect trt AKI"
    concl <- if (aki_dir == "harmful") "MAY INFLATE" else "attenuates protection"
  } else {
    bias <- "negligible"; concl <- "unbiased"
  }
  cat(sprintf("  %-12s  %7.2f  %10s  %12s  %s\n", nm, ratio, aki_dir, bias, concl))
}

# ═══════════════════════════════════════════════════════════════════
# SAVE CSVs
# ═══════════════════════════════════════════════════════════════════
save_df <- function(rows, name) {
  if (length(rows) > 0) {
    df <- do.call(rbind, rows)
    path <- file.path(RESULTS, sprintf("%s_%s.csv", name, tag))
    write.csv(df, path, row.names = FALSE)
    cat(sprintf("  ✓ %s (%d rows)\n", basename(path), nrow(df)))
  }
}

cat(sprintf("\n%s\n  Saving...\n", SEP2))
save_df(simpson_rows,  "qc_simpson")
save_df(thresh_rows,   "qc_threshold")
save_df(kdigo_rows,    "qc_kdigo")
save_df(surv_rows,     "qc_surveillance")

cat(sprintf("\n%s\n05_qc.R — %s DONE\n%s\n", SEP, db, SEP))
