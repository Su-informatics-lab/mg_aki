#!/usr/bin/env Rscript
# ============================================================================
# gen_mortality_table.R — Mortality + ΔCr summary by Mg × eGFR strata
#
# For each stratum:
#   Treated:  n, deaths (%), max ΔCr 7d median [IQR]
#   Control:  n, deaths (%), max ΔCr 7d median [IQR]
#   Effect:   OR [95% CI], P
#
# Usage: Rscript gen_mortality_table.R mimic
#        Rscript gen_mortality_table.R eicu
# ============================================================================

RESULTS <- path.expand("~/mg_aki/results")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) { cat("Usage: Rscript gen_mortality_table.R <db>\n"); quit(status = 1) }
tag <- tolower(args[1]); db <- toupper(tag)

cat(sprintf("\n══ Mortality & ΔCr Summary — %s ══\n", db))

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
ctl_rows <- match(pairs$ctl_pid, all_pts$pid)
n_pairs  <- nrow(pairs)

# ── Covariates ────────────────────────────────────────────────────
egfr_trt    <- all_pts$egfr[trt_rows]
mort_trt    <- all_pts$hosp_mortality[trt_rows]
mort_ctl    <- all_pts$hosp_mortality[ctl_rows]

# Baseline Mg
mg_labs <- labs_all[labs_all$lab_name == "magnesium" & !is.na(labs_all$value), ]
mg_labs <- mg_labs[order(mg_labs$pid, mg_labs$offset_h), ]
mg_list <- split(mg_labs[, c("value", "offset_h")], mg_labs$pid)
mg_trt <- rep(NA_real_, n_pairs)
for (i in seq_len(n_pairs)) {
  pid <- as.character(pairs$trt_pid[i]); t_mg <- pairs$t_mg[i]
  mg_i <- mg_list[[pid]]
  if (is.null(mg_i) || nrow(mg_i) == 0) next
  pre <- mg_i[mg_i$offset_h < t_mg, , drop = FALSE]
  if (nrow(pre) > 0) mg_trt[i] <- pre$value[nrow(pre)]
}

# Max ΔCr 7d
find_cr_pre <- function(cr_pt, t_h) {
  if (is.null(cr_pt) || nrow(cr_pt) == 0) return(NA)
  cand <- cr_pt[cr_pt$offset_h >= 0 & cr_pt$offset_h < t_h, ]
  if (nrow(cand) == 0) return(NA)
  cand$labresult[which.max(cand$offset_h)]
}
max_dcr_7d_trt <- max_dcr_7d_ctl <- rep(NA_real_, n_pairs)
for (i in seq_len(n_pairs)) {
  t_mg_i <- pairs$t_mg[i]
  tp <- as.character(pairs$trt_pid[i]); cp <- as.character(pairs$ctl_pid[i])
  cr_t <- cr_list[[tp]]; cr_c <- cr_list[[cp]]
  pre_t <- find_cr_pre(cr_t, t_mg_i); pre_c <- find_cr_pre(cr_c, t_mg_i)
  if (!is.na(pre_t) && !is.null(cr_t)) {
    post <- cr_t[cr_t$offset_h >= t_mg_i & cr_t$offset_h <= (t_mg_i + 168), ]
    if (nrow(post) > 0) max_dcr_7d_trt[i] <- max(post$labresult - pre_t, na.rm = TRUE)
  }
  if (!is.na(pre_c) && !is.null(cr_c)) {
    post <- cr_c[cr_c$offset_h >= t_mg_i & cr_c$offset_h <= (t_mg_i + 168), ]
    if (nrow(post) > 0) max_dcr_7d_ctl[i] <- max(post$labresult - pre_c, na.rm = TRUE)
  }
}

# ── Strata ────────────────────────────────────────────────────────
mg3 <- rep(NA_character_, n_pairs)
mg3[!is.na(mg_trt) & mg_trt < 1.6]                 <- "Mg<1.6"
mg3[!is.na(mg_trt) & mg_trt >= 1.6 & mg_trt < 2.0] <- "Mg 1.6-2.0"
mg3[!is.na(mg_trt) & mg_trt >= 2.0]                 <- "Mg>=2.0"

egfr_s <- rep(NA_character_, n_pairs)
egfr_s[!is.na(egfr_trt) & egfr_trt >= 90]                 <- "eGFR>=90 (G1)"
egfr_s[!is.na(egfr_trt) & egfr_trt >= 60 & egfr_trt < 90] <- "eGFR 60-89 (G2)"
egfr_s[!is.na(egfr_trt) & egfr_trt >= 45 & egfr_trt < 60] <- "eGFR 45-59 (G3a)"
egfr_s[!is.na(egfr_trt) & egfr_trt < 45]                  <- "eGFR<45 (G3b-5)"

# ── Helper: median [IQR] string ───────────────────────────────────
miq <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 5) return("—")
  sprintf("%.3f [%.3f, %.3f]", median(x), quantile(x, 0.25), quantile(x, 0.75))
}

# ── OR helper ─────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(sandwich); library(lmtest) })
or_str <- function(ot, oc) {
  valid <- !is.na(ot) & !is.na(oc); ot <- ot[valid]; oc <- oc[valid]
  n <- sum(valid)
  if (n < 30 || (sum(ot) + sum(oc)) == 0) return("—")
  df <- data.frame(o = c(ot, oc), t = rep(c(1,0), each = n))
  fit <- tryCatch(glm(o ~ t, data = df, family = quasibinomial()), error = function(e) NULL)
  if (is.null(fit)) return("—")
  ct <- tryCatch(coeftest(fit, vcov. = vcovHC(fit, type = "HC1")),
                 error = function(e) tryCatch(coeftest(fit), error = function(e2) NULL))
  if (is.null(ct)) return("—")
  or <- exp(ct["t", 1]); lo <- exp(ct["t",1] - 1.96*ct["t",2])
  hi <- exp(ct["t",1] + 1.96*ct["t",2]); p <- ct["t", ncol(ct)]
  sig <- if (!is.na(p) && p < 0.05) " *" else ""
  sprintf("%.2f [%.2f, %.2f] P=%.3f%s", or, lo, hi, p, sig)
}

# ── Print table ───────────────────────────────────────────────────
print_block <- function(label, idx) {
  mt <- mort_trt[idx]; mc <- mort_ctl[idx]
  dt <- max_dcr_7d_trt[idx]; dc <- max_dcr_7d_ctl[idx]
  n <- length(idx)
  nt_d <- sum(mt == 1, na.rm = TRUE); nc_d <- sum(mc == 1, na.rm = TRUE)
  cat(sprintf("  %-22s  n=%-5d\n", label, n))
  cat(sprintf("    Treated:   deaths %d/%d (%.1f%%)  max ΔCr 7d: %s\n",
              nt_d, n, 100*nt_d/n, miq(dt)))
  cat(sprintf("    Control:   deaths %d/%d (%.1f%%)  max ΔCr 7d: %s\n",
              nc_d, n, 100*nc_d/n, miq(dc)))
  cat(sprintf("    Mortality OR: %s\n", or_str(mt, mc)))
  cat("\n")
}

cat("\n── OVERALL ──\n")
print_block("All pairs", seq_len(n_pairs))

cat("── BY Mg STRATUM (3-bin) ──\n")
for (s in c("Mg<1.6", "Mg 1.6-2.0", "Mg>=2.0")) {
  idx <- which(mg3 == s)
  if (length(idx) >= 30) print_block(s, idx)
}

cat("── BY eGFR STRATUM ──\n")
for (s in c("eGFR>=90 (G1)", "eGFR 60-89 (G2)", "eGFR 45-59 (G3a)", "eGFR<45 (G3b-5)")) {
  idx <- which(egfr_s == s)
  if (length(idx) >= 30) print_block(s, idx)
}

cat("── eGFR × Mg (key cells) ──\n")
for (eg in c("eGFR>=90 (G1)", "eGFR<45 (G3b-5)")) {
  for (mg in c("Mg<1.6", "Mg 1.6-2.0", "Mg>=2.0")) {
    idx <- which(egfr_s == eg & mg3 == mg)
    if (length(idx) >= 20) print_block(sprintf("%s × %s", eg, mg), idx)
  }
}

# ── Save as CSV ───────────────────────────────────────────────────
rows <- list()
add_row <- function(stratum, idx) {
  mt <- mort_trt[idx]; mc <- mort_ctl[idx]
  dt <- max_dcr_7d_trt[idx]; dc <- max_dcr_7d_ctl[idx]
  n <- length(idx)
  rows[[length(rows) + 1]] <<- data.frame(
    db = db, stratum = stratum, n = n,
    mort_trt_n = sum(mt == 1, na.rm = TRUE),
    mort_trt_pct = round(100 * mean(mt, na.rm = TRUE), 1),
    mort_ctl_n = sum(mc == 1, na.rm = TRUE),
    mort_ctl_pct = round(100 * mean(mc, na.rm = TRUE), 1),
    dcr7d_trt_median = round(median(dt, na.rm = TRUE), 3),
    dcr7d_trt_q25 = round(quantile(dt, 0.25, na.rm = TRUE), 3),
    dcr7d_trt_q75 = round(quantile(dt, 0.75, na.rm = TRUE), 3),
    dcr7d_ctl_median = round(median(dc, na.rm = TRUE), 3),
    dcr7d_ctl_q25 = round(quantile(dc, 0.25, na.rm = TRUE), 3),
    dcr7d_ctl_q75 = round(quantile(dc, 0.75, na.rm = TRUE), 3),
    stringsAsFactors = FALSE)
}

add_row("Overall", seq_len(n_pairs))
for (s in c("Mg<1.6", "Mg 1.6-2.0", "Mg>=2.0"))
  if (sum(mg3 == s, na.rm = TRUE) >= 30) add_row(s, which(mg3 == s))
for (s in c("eGFR>=90 (G1)", "eGFR 60-89 (G2)", "eGFR 45-59 (G3a)", "eGFR<45 (G3b-5)"))
  if (sum(egfr_s == s, na.rm = TRUE) >= 30) add_row(s, which(egfr_s == s))

outpath <- file.path(RESULTS, sprintf("mortality_table_%s.csv", tag))
write.csv(do.call(rbind, rows), outpath, row.names = FALSE)
cat(sprintf("  Saved: %s\n", outpath))
