#!/usr/bin/env Rscript
# ============================================================================
# 03f_hosp_re.R — Hospital random-effects sensitivity (eICU-CRD only)
#
# MIMIC-IV is single-center; only eICU-CRD (208 hospitals) needs this.
# Compares primary quasibinomial + sandwich SE with GLMM (1|hospitalid).
#
# Outputs:
#   results/hosp_re_eicu.csv  — side-by-side OR comparison
#
# Usage: Rscript 03f_hosp_re.R
# ============================================================================

suppressPackageStartupMessages({
  library(lme4)
  library(sandwich)
  library(lmtest)
})

set.seed(2026)

RESULTS   <- path.expand("~/mg_aki/results")
EICU_ROOT <- if (dir.exists(path.expand("~/mg_aki/eicu-crd-2.0")))
               path.expand("~/mg_aki/eicu-crd-2.0") else
               path.expand("~/mg_aki/eicu-crd-demo")

SEP <- paste(rep("=", 60), collapse = "")
cat(sprintf("\n%s\n03f_hosp_re.R — Hospital RE sensitivity (eICU-CRD)\n%s\n", SEP, SEP))

# ── Load matched pairs + patient data ─────────────────────────────
all_pts <- read.csv(file.path(RESULTS, "did_all_eicu.csv"), stringsAsFactors = FALSE)
pairs   <- read.csv(file.path(RESULTS, "did_pairs_primary_yet_untreated_eicu.csv"),
                     stringsAsFactors = FALSE)
cr_all  <- read.csv(file.path(RESULTS, "did_cr_all_eicu.csv"), stringsAsFactors = FALSE)

cat(sprintf("  Pairs: %d | Patients: %d\n", nrow(pairs), nrow(all_pts)))

# ── Load hospitalid from eICU patient table ───────────────────────
cat("  Loading hospitalid from patient table...\n")
patient <- read.csv(gzfile(file.path(EICU_ROOT, "patient.csv.gz")),
                    stringsAsFactors = FALSE)
hosp_map <- setNames(patient$hospitalid,
                     as.character(patient$patientunitstayid))
all_pts$hospitalid <- hosp_map[as.character(all_pts$pid)]
n_hosp <- length(unique(all_pts$hospitalid[!is.na(all_pts$hospitalid)]))
cat(sprintf("  Hospitals mapped: %d unique\n", n_hosp))

# ── Creatinine setup (same as 03b) ───────────────────────────────
cr_all$pid <- cr_all$patientunitstayid
if (!"offset_h" %in% names(cr_all)) cr_all$offset_h <- cr_all$labresultoffset / 60
cr_all <- cr_all[order(cr_all$pid, cr_all$offset_h), ]
cr_list <- split(cr_all[, c("labresult", "offset_h")], cr_all$pid)

n_pairs  <- nrow(pairs)
trt_rows <- match(pairs$trt_pid, all_pts$pid)
ctl_rows <- match(pairs$ctl_pid, all_pts$pid)

# ── AKI computation (identical to 03b) ────────────────────────────
compute_aki <- function(pid, t_mg, window_h) {
  cr <- cr_list[[as.character(pid)]]
  if (is.null(cr) || nrow(cr) < 1) return(NA)
  pre <- cr[cr$offset_h >= 0 & cr$offset_h < t_mg, ]
  if (nrow(pre) == 0) return(NA)
  bl <- pre$labresult[which.max(pre$offset_h)]
  if (is.na(bl) || bl <= 0) return(NA)
  post <- cr[cr$offset_h > t_mg & cr$offset_h <= (t_mg + window_h), ]
  if (nrow(post) == 0) return(0L)
  for (i in seq_len(nrow(post))) {
    h <- post$offset_h[i] - t_mg; val <- post$labresult[i]
    delta <- val - bl; ratio <- val / bl
    if (window_h <= 48) {
      if (delta >= 0.3 || ratio >= 1.5) return(1L)
    } else {
      if (h <= 48 && (delta >= 0.3 || ratio >= 1.5)) return(1L)
      if (h > 48 && ratio >= 1.5) return(1L)
    }
  }
  return(0L)
}

cat("  Computing AKI outcomes (48h + 7d)...\n")
aki48_trt <- aki48_ctl <- aki7d_trt <- aki7d_ctl <- rep(NA_integer_, n_pairs)
for (i in seq_len(n_pairs)) {
  aki48_trt[i] <- compute_aki(pairs$trt_pid[i], pairs$t_mg[i], 48)
  aki48_ctl[i] <- compute_aki(pairs$ctl_pid[i], pairs$t_mg[i], 48)
  aki7d_trt[i] <- compute_aki(pairs$trt_pid[i], pairs$t_mg[i], 168)
  aki7d_ctl[i] <- compute_aki(pairs$ctl_pid[i], pairs$t_mg[i], 168)
}
cat(sprintf("  48h AKI valid: %d pairs | 7d AKI valid: %d pairs\n",
            sum(!is.na(aki48_trt) & !is.na(aki48_ctl)),
            sum(!is.na(aki7d_trt) & !is.na(aki7d_ctl))))

# ── eGFR strata (by treated patient, same as 03b) ────────────────
egfr_trt <- all_pts$egfr[trt_rows]
ckd_stage <- rep(NA_character_, n_pairs)
ckd_stage[!is.na(egfr_trt) & egfr_trt >= 90]                  <- "eGFR>=90"
ckd_stage[!is.na(egfr_trt) & egfr_trt >= 60 & egfr_trt < 90]  <- "eGFR_60-89"
ckd_stage[!is.na(egfr_trt) & egfr_trt >= 45 & egfr_trt < 60]  <- "eGFR_45-59"
ckd_stage[!is.na(egfr_trt) & egfr_trt >= 30 & egfr_trt < 45]  <- "eGFR_30-44"
ckd_stage[!is.na(egfr_trt) & egfr_trt < 30]                    <- "eGFR<30"

# ── Hospital IDs for each pair ────────────────────────────────────
hosp_trt <- all_pts$hospitalid[trt_rows]
hosp_ctl <- all_pts$hospitalid[ctl_rows]

# ── Model comparison function ─────────────────────────────────────
compare_models <- function(outcome_trt, outcome_ctl, hosp_t, hosp_c, label, stratum) {
  valid <- !is.na(outcome_trt) & !is.na(outcome_ctl) &
           !is.na(hosp_t) & !is.na(hosp_c)
  ot <- outcome_trt[valid]; oc <- outcome_ctl[valid]
  ht <- hosp_t[valid];     hc <- hosp_c[valid]
  n <- sum(valid); nt <- sum(ot); nc <- sum(oc)

  if (n < 30 || (nt + nc) == 0)
    return(data.frame(outcome = label, stratum = stratum, n = n,
                      or_primary = NA, lo_primary = NA, hi_primary = NA, p_primary = NA,
                      or_glmm = NA, lo_glmm = NA, hi_glmm = NA, p_glmm = NA,
                      n_hosp = NA, var_hosp = NA, note = "skip"))

  # Long-format data: stack treated + control
  df <- data.frame(
    outcome  = c(ot, oc),
    treated  = rep(c(1L, 0L), each = n),
    hosp     = factor(c(ht, hc))
  )
  n_h <- nlevels(df$hosp)

  # ── (A) Primary: quasibinomial + sandwich SE ──
  fit_a <- tryCatch(glm(outcome ~ treated, data = df, family = quasibinomial()),
                    error = function(e) NULL)
  if (!is.null(fit_a)) {
    ct_a <- tryCatch(coeftest(fit_a, vcov. = vcovHC(fit_a, type = "HC1")),
                     error = function(e) NULL)
  } else ct_a <- NULL
  if (!is.null(ct_a)) {
    or_a <- exp(ct_a["treated", "Estimate"])
    lo_a <- exp(ct_a["treated", "Estimate"] - 1.96 * ct_a["treated", "Std. Error"])
    hi_a <- exp(ct_a["treated", "Estimate"] + 1.96 * ct_a["treated", "Std. Error"])
    p_a  <- ct_a["treated", ncol(ct_a)]
  } else { or_a <- lo_a <- hi_a <- p_a <- NA }

  # ── (B) GLMM: binomial + (1|hospital) ──
  note <- ""
  fit_b <- tryCatch(
    glmer(outcome ~ treated + (1 | hosp), data = df, family = binomial(),
          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))),
    warning = function(w) {
      note <<- paste0("warn:", conditionMessage(w))
      suppressWarnings(
        glmer(outcome ~ treated + (1 | hosp), data = df, family = binomial(),
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
      )
    },
    error = function(e) { note <<- paste0("fail:", conditionMessage(e)); NULL }
  )

  if (!is.null(fit_b)) {
    cf <- summary(fit_b)$coefficients
    or_b <- exp(cf["treated", "Estimate"])
    se_b <- cf["treated", "Std. Error"]
    lo_b <- exp(cf["treated", "Estimate"] - 1.96 * se_b)
    hi_b <- exp(cf["treated", "Estimate"] + 1.96 * se_b)
    p_b  <- cf["treated", "Pr(>|z|)"]
    var_h <- as.numeric(VarCorr(fit_b)$hosp)
  } else { or_b <- lo_b <- hi_b <- p_b <- var_h <- NA }

  data.frame(
    outcome    = label,
    stratum    = stratum,
    n          = n,
    or_primary = round(or_a, 4),
    lo_primary = round(lo_a, 4),
    hi_primary = round(hi_a, 4),
    p_primary  = round(p_a, 6),
    or_glmm    = round(or_b, 4),
    lo_glmm    = round(lo_b, 4),
    hi_glmm    = round(hi_b, 4),
    p_glmm     = round(p_b, 6),
    n_hosp     = n_h,
    var_hosp   = round(var_h, 6),
    note       = note,
    stringsAsFactors = FALSE
  )
}

# ── Run all strata ────────────────────────────────────────────────
strata_order <- c("Overall", "eGFR>=90", "eGFR_60-89",
                   "eGFR_45-59", "eGFR_30-44", "eGFR<30")
outcomes <- list(
  list(trt = aki48_trt, ctl = aki48_ctl, label = "AKI_48h_Stage1+"),
  list(trt = aki7d_trt, ctl = aki7d_ctl, label = "AKI_7d_Stage1+")
)

# Also add hospital mortality (already in all_pts)
mort_trt <- all_pts$hosp_mortality[trt_rows]
mort_ctl <- all_pts$hosp_mortality[ctl_rows]
outcomes[[3]] <- list(trt = mort_trt, ctl = mort_ctl, label = "hosp_mortality")

results <- list()
row_idx <- 0

for (oc in outcomes) {
  for (stg in strata_order) {
    if (stg == "Overall") { idx <- seq_len(n_pairs) }
    else { idx <- which(ckd_stage == stg) }
    if (length(idx) < 30) next

    cat(sprintf("  [%s x %s] n=%d ... ", oc$label, stg, length(idx)))
    res <- compare_models(
      oc$trt[idx], oc$ctl[idx],
      hosp_trt[idx], hosp_ctl[idx],
      oc$label, stg
    )
    row_idx <- row_idx + 1
    results[[row_idx]] <- res

    if (!is.na(res$or_primary) && !is.na(res$or_glmm)) {
      cat(sprintf("primary=%.3f  glmm=%.3f  var_h=%.4f  %s\n",
                  res$or_primary, res$or_glmm, res$var_hosp, res$note))
    } else {
      cat(sprintf("skip (%s)\n", res$note))
    }
  }
}

res_df <- do.call(rbind, results)
res_df$stratum <- factor(res_df$stratum, levels = strata_order)

# ── Summary ───────────────────────────────────────────────────────
cat(sprintf("\n%s\nSummary: %d rows\n%s\n", SEP, nrow(res_df), SEP))
for (oc in unique(res_df$outcome)) {
  cat(sprintf("\n  [%s]\n", oc))
  sub <- res_df[res_df$outcome == oc, ]
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    if (is.na(r$or_primary)) { cat(sprintf("    %-15s  n=%d (skip)\n", r$stratum, r$n)); next }
    cat(sprintf("    %-15s  Primary: %.3f [%.3f,%.3f]  GLMM: %.3f [%.3f,%.3f]  var_h=%.4f  n_hosp=%d\n",
                r$stratum,
                r$or_primary, r$lo_primary, r$hi_primary,
                r$or_glmm, r$lo_glmm, r$hi_glmm,
                r$var_hosp, r$n_hosp))
  }
}

# ── Save ──────────────────────────────────────────────────────────
outpath <- file.path(RESULTS, "hosp_re_eicu.csv")
write.csv(res_df, outpath, row.names = FALSE)
cat(sprintf("\n  Saved: %s\n", outpath))

cat(sprintf("\n%s\n03f_hosp_re.R done\n%s\n", SEP, SEP))
