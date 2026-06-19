#!/usr/bin/env Rscript
# ============================================================================
# did_subgroups.R — Mg-stratified & surgery subgroups across all configs
#
# For each matched dataset (excluding peak):
#   - Mg strata: <1.8, 1.8-2.0, 2.0-2.3, >2.3
#   - Surgery: CABG, valve, combined, other_cardiac
#   - Interaction test: treatment × Mg (continuous)
#
# Stratification is by the TREATED patient's value; matched controls
# follow their partner regardless of their own stratum.
#
# Run:  Rscript did_subgroups.R eicu
#       Rscript did_subgroups.R mimic
# ============================================================================

suppressPackageStartupMessages({
  library(sandwich)
  library(lmtest)
})

RESULTS <- path.expand("~/mg_aki/results")

# Configs to analyze (all except peak)
CONFIGS <- c("primary", "sens_t4", "sens_r2", "sens_r4t4", "sens_at24")

CONFIG_LABELS <- c(
  primary  = "r=1 ±6h first_6_24h",
  sens_t4  = "r=1 ±4h first_6_24h",
  sens_r2  = "r=2 ±6h first_6_24h",
  sens_r4t4= "r=4 ±4h first_6_24h",
  sens_at24= "r=1 ±6h closest_24h"
)

PS_COVARS <- c(
  "age", "is_female", "bmi",
  "surg_cabg", "surg_valve", "surg_combined",
  "heart_failure", "hypertension", "diabetes", "ckd",
  "copd", "pvd", "stroke", "liver_disease",
  "egfr",
  "loop_diuretics", "nsaids", "acei_arb", "ppi",
  "beta_blockers", "steroids", "antiarrhythmics",
  "first_potassium", "first_calcium", "first_heartrate",
  "first_mg_value", "first_lactate", "lactate_missing"
)

MG_BINS <- list(
  "<1.8"    = c(0, 1.8),
  "1.8-2.0" = c(1.8, 2.0),
  "2.0-2.3" = c(2.0, 2.3),
  ">2.3"    = c(2.3, 99)
)

SURG_TYPES <- c("cabg", "valve", "combined", "other_cardiac")

# ── Helpers ──────────────────────────────────────────────────────────────
median_impute <- function(d, vars) {
  for (v in vars)
    if (v %in% names(d) && any(is.na(d[[v]])))
      d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm = TRUE)
  d
}

did_est <- function(df, ps_vars) {
  # Doubly robust DiD with cluster SEs
  n_trt <- sum(df$treated == 1)
  n_ctl <- sum(df$treated == 0)
  if (n_trt < 10 || n_ctl < 10) return(NULL)

  # Find covariates with SMD > 0.05
  adjust <- character(0)
  for (v in ps_vars) {
    if (!v %in% names(df)) next
    x1 <- df[[v]][df$treated==1]; x0 <- df[[v]][df$treated==0]
    sp <- sqrt((var(x1,na.rm=T)+var(x0,na.rm=T))/2)
    if (!is.na(sp) && sp > 1e-10) {
      smd <- abs(mean(x1,na.rm=T)-mean(x0,na.rm=T))/sp
      if (!is.na(smd) && smd > 0.05) adjust <- c(adjust, v)
    }
  }

  # Unadjusted
  fit0 <- tryCatch(lm(delta_cr ~ treated, data=df), error=function(e) NULL)
  if (is.null(fit0)) return(NULL)
  cl0 <- tryCatch(
    if (length(unique(df$match_pair_id)) > 1)
      vcovCL(fit0, cluster=df$match_pair_id) else vcovHC(fit0, type="HC1"),
    error = function(e) vcovHC(fit0, type="HC1"))
  ct0 <- coeftest(fit0, vcov.=cl0)

  # Adjusted
  if (length(adjust) > 0) {
    avail <- intersect(adjust, names(df))
    # Remove near-constant vars
    avail <- avail[sapply(avail, function(v) var(df[[v]], na.rm=T) > 1e-10)]
    if (length(avail) > 0) {
      fml <- as.formula(paste("delta_cr ~ treated +", paste(avail, collapse="+")))
      fit1 <- tryCatch(lm(fml, data=df), error=function(e) NULL)
      if (is.null(fit1)) { ct1 <- ct0 } else {
        cl1 <- tryCatch(
          if (length(unique(df$match_pair_id)) > 1)
            vcovCL(fit1, cluster=df$match_pair_id) else vcovHC(fit1, type="HC1"),
          error = function(e) vcovHC(fit1, type="HC1"))
        ct1 <- coeftest(fit1, vcov.=cl1)
      }
    } else { ct1 <- ct0 }
  } else { ct1 <- ct0 }

  list(
    n_trt = n_trt, n_ctl = n_ctl,
    did_unadj = ct0["treated","Estimate"],
    p_unadj = ct0["treated","Pr(>|t|)"],
    did_adj = ct1["treated","Estimate"],
    se_adj = ct1["treated","Std. Error"],
    p_adj = ct1["treated","Pr(>|t|)"],
    ci_lo = ct1["treated","Estimate"] - 1.96*ct1["treated","Std. Error"],
    ci_hi = ct1["treated","Estimate"] + 1.96*ct1["treated","Std. Error"]
  )
}

interaction_test <- function(df, ps_vars) {
  # Treatment × continuous Mg interaction
  if (!"first_mg_value" %in% names(df)) return(list(int_p=NA))
  avail <- intersect(ps_vars[ps_vars != "first_mg_value"], names(df))
  avail <- avail[sapply(avail, function(v) var(df[[v]], na.rm=T) > 1e-10)]
  fml <- as.formula(paste("delta_cr ~ treated * first_mg_value +",
                          paste(avail, collapse="+")))
  fit <- tryCatch(lm(fml, data=df), error=function(e) NULL)
  if (is.null(fit)) return(list(int_p=NA))
  cl <- tryCatch(
    if (length(unique(df$match_pair_id)) > 1)
      vcovCL(fit, cluster=df$match_pair_id) else vcovHC(fit, type="HC1"),
    error = function(e) vcovHC(fit, type="HC1"))
  ct <- coeftest(fit, vcov.=cl)
  int_row <- grep("treated:first_mg_value", rownames(ct))
  if (length(int_row) == 0) return(list(int_p=NA))
  list(
    int_est = ct[int_row, "Estimate"],
    int_se = ct[int_row, "Std. Error"],
    int_p = ct[int_row, "Pr(>|t|)"]
  )
}

# ============================================================================
run_subgroups <- function(db) {
  tag <- tolower(db)
  SEP <- paste(rep("=", 70), collapse="")
  cat(sprintf("\n%s\n%s: Subgroup Analysis (Mg strata + surgery type)\n%s\n", SEP, db, SEP))

  all_rows <- list()
  idx <- 0

  for (cfg in CONFIGS) {
    path <- file.path(RESULTS, sprintf("did_matched_%s_%s.csv", tag, cfg))
    if (!file.exists(path)) {
      cat(sprintf("  %s: file not found, skipping\n", cfg))
      next
    }

    df <- read.csv(path, stringsAsFactors=FALSE)
    df <- median_impute(df, PS_COVARS)
    n_trt <- sum(df$treated==1)
    n_ctl <- sum(df$treated==0)

    cat(sprintf("\n%s\n  Config: %s (%s)\n  N: %d trt, %d ctl\n%s\n",
                paste(rep("─",60), collapse=""), cfg, CONFIG_LABELS[cfg],
                n_trt, n_ctl, paste(rep("─",60), collapse="")))

    # ── Overall (for reference) ────────────────────────────────────────
    res_all <- did_est(df, PS_COVARS)
    if (!is.null(res_all)) {
      cat(sprintf("  Overall: DiD=%+.4f, P=%.4f, 95%%CI [%+.4f,%+.4f]\n",
                  res_all$did_adj, res_all$p_adj, res_all$ci_lo, res_all$ci_hi))
      idx <- idx+1
      all_rows[[idx]] <- data.frame(config=cfg, subgroup="Overall", stratum="all",
                                     n_trt=res_all$n_trt, n_ctl=res_all$n_ctl,
                                     did_adj=round(res_all$did_adj,4),
                                     p_adj=round(res_all$p_adj,4),
                                     ci_lo=round(res_all$ci_lo,4),
                                     ci_hi=round(res_all$ci_hi,4),
                                     stringsAsFactors=F)
    }

    # ── Mg strata ──────────────────────────────────────────────────────
    cat("\n  Mg strata (by treated patient's first serum Mg):\n")

    # Assign strata by treated patient's Mg (controls follow their partner)
    trt_mg <- df$first_mg_value[df$treated == 1]
    pair_ids <- df$match_pair_id[df$treated == 1]

    for (mg_name in names(MG_BINS)) {
      lo <- MG_BINS[[mg_name]][1]; hi <- MG_BINS[[mg_name]][2]
      # Treated patients in this stratum
      in_stratum <- !is.na(trt_mg) & trt_mg >= lo & trt_mg < hi
      pairs_in <- pair_ids[in_stratum]

      sub <- df[df$match_pair_id %in% pairs_in, ]
      nt <- sum(sub$treated==1); nc <- sum(sub$treated==0)

      if (nt < 10 || nc < 10) {
        cat(sprintf("    %s: n_trt=%d — too few\n", mg_name, nt))
        idx <- idx+1
        all_rows[[idx]] <- data.frame(config=cfg, subgroup="Mg_stratum", stratum=mg_name,
                                       n_trt=nt, n_ctl=nc,
                                       did_adj=NA, p_adj=NA, ci_lo=NA, ci_hi=NA,
                                       stringsAsFactors=F)
        next
      }

      res <- did_est(sub, PS_COVARS)
      if (!is.null(res)) {
        sig <- if (res$p_adj < 0.05) " *" else ""
        cat(sprintf("    %s: n=%d/%d, DiD=%+.4f, P=%.4f%s\n",
                    mg_name, nt, nc, res$did_adj, res$p_adj, sig))
        idx <- idx+1
        all_rows[[idx]] <- data.frame(config=cfg, subgroup="Mg_stratum", stratum=mg_name,
                                       n_trt=res$n_trt, n_ctl=res$n_ctl,
                                       did_adj=round(res$did_adj,4),
                                       p_adj=round(res$p_adj,4),
                                       ci_lo=round(res$ci_lo,4),
                                       ci_hi=round(res$ci_hi,4),
                                       stringsAsFactors=F)
      }
    }

    # Interaction test
    int <- interaction_test(df, PS_COVARS)
    if (!is.na(int$int_p)) {
      cat(sprintf("    Interaction (trt × Mg): est=%+.4f, P=%.4f\n",
                  int$int_est, int$int_p))
      idx <- idx+1
      all_rows[[idx]] <- data.frame(config=cfg, subgroup="Mg_interaction",
                                     stratum="trt×Mg",
                                     n_trt=n_trt, n_ctl=n_ctl,
                                     did_adj=round(int$int_est,4),
                                     p_adj=round(int$int_p,4),
                                     ci_lo=NA, ci_hi=NA,
                                     stringsAsFactors=F)
    }

    # ── Surgery type ───────────────────────────────────────────────────
    cat("\n  Surgery type:\n")

    # Need surgery_type column — check if available
    if ("surgery_type" %in% names(df)) {
      for (st in SURG_TYPES) {
        # Use treated patient's surgery type
        trt_surg <- df$surgery_type[df$treated == 1]
        in_st <- !is.na(trt_surg) & trt_surg == st
        pairs_st <- pair_ids[in_st]
        sub <- df[df$match_pair_id %in% pairs_st, ]
        nt <- sum(sub$treated==1); nc <- sum(sub$treated==0)

        if (nt < 10 || nc < 10) {
          cat(sprintf("    %s: n_trt=%d — too few\n", st, nt))
          idx <- idx+1
          all_rows[[idx]] <- data.frame(config=cfg, subgroup="Surgery", stratum=st,
                                         n_trt=nt, n_ctl=nc,
                                         did_adj=NA, p_adj=NA, ci_lo=NA, ci_hi=NA,
                                         stringsAsFactors=F)
          next
        }

        res <- did_est(sub, PS_COVARS)
        if (!is.null(res)) {
          sig <- if (res$p_adj < 0.05) " *" else ""
          cat(sprintf("    %s: n=%d/%d, DiD=%+.4f, P=%.4f%s\n",
                      st, nt, nc, res$did_adj, res$p_adj, sig))
          idx <- idx+1
          all_rows[[idx]] <- data.frame(config=cfg, subgroup="Surgery", stratum=st,
                                         n_trt=res$n_trt, n_ctl=res$n_ctl,
                                         did_adj=round(res$did_adj,4),
                                         p_adj=round(res$p_adj,4),
                                         ci_lo=round(res$ci_lo,4),
                                         ci_hi=round(res$ci_hi,4),
                                         stringsAsFactors=F)
        }
      }
    } else {
      cat("    surgery_type not in matched data — skipped\n")
    }
  }

  # ── Save & Summary ────────────────────────────────────────────────────
  results <- do.call(rbind, all_rows)
  out_path <- file.path(RESULTS, sprintf("did_subgroups_%s.csv", tag))
  write.csv(results, out_path, row.names=FALSE)
  cat(sprintf("\n  Saved: %s (%d rows)\n", basename(out_path), nrow(results)))

  # ── Formatted summary: Mg strata across configs ────────────────────────
  cat(sprintf("\n%s\n%s: Mg-STRATIFIED SUMMARY (doubly robust DiD)\n%s\n",
              SEP, db, SEP))

  mg_res <- results[results$subgroup == "Mg_stratum", ]
  int_res <- results[results$subgroup == "Mg_interaction", ]

  if (nrow(mg_res) > 0) {
    cat("\n  config          stratum   n_trt  n_ctl    DiD_adj    P       95% CI\n")
    cat("  ──────────────  ────────  ─────  ─────  ─────────  ──────  ──────────────\n")
    for (i in seq_len(nrow(mg_res))) {
      r <- mg_res[i, ]
      sig <- if (!is.na(r$p_adj) && r$p_adj < 0.05) " *" else ""
      did_str <- if (is.na(r$did_adj)) "     NA" else sprintf("%+.4f", r$did_adj)
      p_str <- if (is.na(r$p_adj)) "    NA" else sprintf("%.4f", r$p_adj)
      ci_str <- if (is.na(r$ci_lo)) "      NA" else
        sprintf("[%+.4f,%+.4f]", r$ci_lo, r$ci_hi)
      cat(sprintf("  %-14s  %-8s  %5d  %5d  %s  %s  %s%s\n",
                  r$config, r$stratum, r$n_trt, r$n_ctl,
                  did_str, p_str, ci_str, sig))
    }
  }

  if (nrow(int_res) > 0) {
    cat("\n  Interaction tests (trt × continuous Mg):\n")
    for (i in seq_len(nrow(int_res))) {
      r <- int_res[i, ]
      sig <- if (!is.na(r$p_adj) && r$p_adj < 0.05) " *" else ""
      cat(sprintf("    %s: β=%+.4f, P=%.4f%s\n",
                  r$config, r$did_adj, r$p_adj, sig))
    }
  }

  # Surgery summary
  cat(sprintf("\n%s\n%s: SURGERY TYPE SUMMARY\n%s\n", SEP, db, SEP))
  sg_res <- results[results$subgroup == "Surgery", ]
  if (nrow(sg_res) > 0) {
    cat("\n  config          surgery         n_trt  n_ctl    DiD_adj    P\n")
    cat("  ──────────────  ──────────────  ─────  ─────  ─────────  ──────\n")
    for (i in seq_len(nrow(sg_res))) {
      r <- sg_res[i, ]
      sig <- if (!is.na(r$p_adj) && r$p_adj < 0.05) " *" else ""
      did_str <- if (is.na(r$did_adj)) "     NA" else sprintf("%+.4f", r$did_adj)
      p_str <- if (is.na(r$p_adj)) "    NA" else sprintf("%.4f", r$p_adj)
      cat(sprintf("  %-14s  %-14s  %5d  %5d  %s  %s%s\n",
                  r$config, r$stratum, r$n_trt, r$n_ctl,
                  did_str, p_str, sig))
    }
  }

  return(results)
}

# ============================================================================
cat("======================================================================\n")
cat("did_subgroups.R — Mg-stratified & surgery subgroups\n")
cat(sprintf("  Configs: %s\n", paste(CONFIGS, collapse=", ")))
cat(sprintf("  Mg strata: %s\n", paste(names(MG_BINS), collapse=", ")))
cat("======================================================================\n")

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) { cat("Usage: Rscript did_subgroups.R eicu|mimic\n"); quit(status=1) }
for (a in args) run_subgroups(toupper(a))

cat("\n======================================================================\n")
cat("Done. Key question: does the >2.3 Mg stratum show the concentration-\n")
cat("dependent effect seen in the old TTE?\n")
cat("======================================================================\n")
