#!/usr/bin/env Rscript
# ============================================================================
# 00b_missing_assessment.R — Missingness diagnostics + MICE validation
#
# For JNO: demonstrate that missingness is handled properly
#   1. Missingness rates by variable × treatment group
#   2. Differential missingness tests (chi-squared)
#   3. Missingness pattern analysis
#   4. Assess MCAR vs MAR (logistic regression of missingness indicators)
#
# Run: Rscript 00b_missing_assessment.R eicu
#      Rscript 00b_missing_assessment.R mimic
# ============================================================================

source(file.path(path.expand("~/mg_aki"), "did_covars.R"))

RESULTS   <- path.expand("~/mg_aki/results")
PS_COVARS <- PS_PRIMARY

run_assessment <- function(db) {
  tag <- tolower(db)
  SEP <- paste(rep("=",70), collapse="")
  cat(sprintf("\n%s\n%s: Missingness Assessment\n%s\n", SEP, db, SEP))

  trt <- read.csv(file.path(RESULTS,sprintf("did_treated_%s.csv",tag)),stringsAsFactors=F)
  ctl <- read.csv(file.path(RESULTS,sprintf("did_control_%s.csv",tag)),stringsAsFactors=F)
  id_col <- if("patientunitstayid" %in% names(trt)) "patientunitstayid" else "stay_id"
  trt$pid <- trt[[id_col]]; ctl$pid <- ctl[[id_col]]
  sh <- intersect(names(trt), names(ctl))
  full <- rbind(trt[,sh], ctl[,sh])

  # Variables to assess
  vars_to_check <- c("bmi", "first_potassium", "first_calcium",
                       "first_lactate", "first_mg_value", "first_heartrate", "egfr")

  # ── 1. Missingness rates ────────────────────────────────────────────
  cat("\n  1. MISSINGNESS RATES\n")
  cat("  variable            overall   treated   control   diff     chi2_P    concern\n")
  cat("  ──────────────────  ────────  ────────  ────────  ─────── ─────────  ───────\n")

  miss_rows <- list()
  for (v in vars_to_check) {
    if (!v %in% names(full)) next
    n_all <- nrow(full); n_trt <- sum(full$treated==1); n_ctl <- sum(full$treated==0)
    miss_all <- sum(is.na(full[[v]])); rate_all <- miss_all/n_all
    miss_trt <- sum(is.na(full[[v]][full$treated==1])); rate_trt <- miss_trt/n_trt
    miss_ctl <- sum(is.na(full[[v]][full$treated==0])); rate_ctl <- miss_ctl/n_ctl
    diff_pct <- rate_trt - rate_ctl

    # Chi-squared test for differential missingness
    tbl <- matrix(c(miss_trt, n_trt-miss_trt, miss_ctl, n_ctl-miss_ctl), nrow=2)
    chi_p <- tryCatch(chisq.test(tbl)$p.value, error=function(e) NA)

    concern <- if(rate_all > 0.3) "HIGH" else if(rate_all > 0.1) "moderate" else "low"
    if (!is.na(chi_p) && chi_p < 0.001 && abs(diff_pct) > 0.05) concern <- paste(concern, "+ DIFFERENTIAL")

    cat(sprintf("  %-20s %6.1f%%   %6.1f%%   %6.1f%%   %+5.1f%%  %.2e   %s\n",
                v, 100*rate_all, 100*rate_trt, 100*rate_ctl, 100*diff_pct,
                ifelse(is.na(chi_p), NA, chi_p), concern))

    miss_rows[[v]] <- data.frame(variable=v, rate_all=round(100*rate_all,1),
      rate_trt=round(100*rate_trt,1), rate_ctl=round(100*rate_ctl,1),
      diff_pct=round(100*diff_pct,1), chi2_p=round(chi_p,6), concern=concern,
      stringsAsFactors=F)
  }

  # ── 2. Missingness predictors (MAR assessment) ─────────────────────
  cat("\n  2. MAR ASSESSMENT: Does missingness depend on observed variables?\n")
  cat("     (Logistic regression: missing_indicator ~ treatment + age + sex + surgery + comorbidities)\n\n")

  predictors <- c("treated","age","is_female","surg_cabg","surg_valve",
                    "heart_failure","diabetes","ckd","egfr")
  predictors <- intersect(predictors, names(full))

  for (v in vars_to_check) {
    if (!v %in% names(full)) next
    if (sum(is.na(full[[v]])) < 10) next

    full$miss_ind <- as.integer(is.na(full[[v]]))
    pred_avail <- setdiff(predictors, v)  # don't predict from itself
    pred_avail <- pred_avail[sapply(pred_avail, function(p) !all(is.na(full[[p]])))]

    fml <- as.formula(paste("miss_ind ~", paste(pred_avail, collapse="+")))
    fit <- tryCatch(glm(fml, data=full, family=binomial()), error=function(e) NULL)
    if (is.null(fit)) next

    # Overall model significance (likelihood ratio test)
    null_fit <- glm(miss_ind ~ 1, data=full, family=binomial())
    lr_test <- anova(null_fit, fit, test="Chisq")
    lr_p <- lr_test$`Pr(>Chi)`[2]

    # Treatment effect on missingness
    ct <- summary(fit)$coefficients
    trt_row <- which(rownames(ct)=="treated")
    trt_or <- if(length(trt_row)>0) exp(ct[trt_row,"Estimate"]) else NA
    trt_p <- if(length(trt_row)>0) ct[trt_row,"Pr(>|z|)"] else NA

    mechanism <- if(!is.na(lr_p) && lr_p < 0.001) {
      if(!is.na(trt_p) && trt_p < 0.05) "MAR (treatment-dependent)" else "MAR (covariate-dependent)"
    } else "Consistent with MCAR"

    cat(sprintf("  %-18s LR_P=%.2e  Trt_OR=%.2f  Trt_P=%.4f  → %s\n",
                v, lr_p, ifelse(is.na(trt_or),1,trt_or),
                ifelse(is.na(trt_p),1,trt_p), mechanism))
  }

  # ── 3. Missingness patterns ────────────────────────────────────────
  cat("\n  3. MISSINGNESS PATTERNS (top 10)\n")
  pattern_vars <- intersect(vars_to_check, names(full))
  miss_mat <- !is.na(full[,pattern_vars])
  pattern_str <- apply(miss_mat, 1, function(x) paste(as.integer(x), collapse=""))
  pattern_tab <- sort(table(pattern_str), decreasing=TRUE)
  cat(sprintf("     %s  | n      | %%\n", paste(sprintf("%-5s", pattern_vars), collapse=" ")))
  cat(sprintf("     %s  | ───── | ─────\n", paste(rep("─────", length(pattern_vars)), collapse=" ")))
  for (i in seq_len(min(10, length(pattern_tab)))) {
    p <- names(pattern_tab)[i]
    bits <- as.integer(strsplit(p, "")[[1]])
    labels <- ifelse(bits==1, "  ✓  ", " miss")
    cat(sprintf("     %s  | %5d | %5.1f%%\n",
                paste(labels, collapse=" "), pattern_tab[i],
                100*pattern_tab[i]/nrow(full)))
  }
  complete <- sum(complete.cases(full[,pattern_vars]))
  cat(sprintf("\n  Complete cases: %d/%d (%.1f%%)\n", complete, nrow(full), 100*complete/nrow(full)))

  # ── 4. Recommendation ──────────────────────────────────────────────
  cat(sprintf("\n  4. RECOMMENDATION\n"))

  has_high_miss <- any(sapply(miss_rows, function(r) r$rate_all > 30))
  has_differential <- any(sapply(miss_rows, function(r) grepl("DIFFERENTIAL", r$concern)))

  if (has_high_miss) {
    cat("  ⚠ HIGH MISSINGNESS (>30%) detected. MICE required.\n")
    cat("    Variables with >30%: ")
    high <- sapply(miss_rows, function(r) if(r$rate_all>30) r$variable else NA)
    cat(paste(na.omit(high), collapse=", "), "\n")
  }
  if (has_differential) {
    cat("  ⚠ DIFFERENTIAL MISSINGNESS detected (rate differs by treatment).\n")
    cat("    → Must include treatment in MICE imputation model.\n")
  }
  cat("  → Primary: MICE (m=10, PMM for continuous, logreg for binary)\n")
  cat("  → Sensitivity: complete case analysis\n")
  cat("  → Include outcome (delta_cr) and treatment in imputation model\n")

  # Save
  miss_df <- do.call(rbind, miss_rows)
  write.csv(miss_df, file.path(RESULTS, sprintf("missingness_%s.csv", tag)), row.names=F)
  cat(sprintf("\n  Saved: missingness_%s.csv\n", tag))
}

# ============================================================================
cat("======================================================================\n")
cat("00b_missing_assessment.R — Missingness diagnostics for JNO\n")
cat("======================================================================\n")
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) { cat("Usage: Rscript 00b_missing_assessment.R eicu|mimic\n"); quit(status=1) }
for (a in args) run_assessment(toupper(a))
