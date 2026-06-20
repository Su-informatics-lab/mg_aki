#!/usr/bin/env Rscript
# ============================================================================
# 01d_hte_mice.R — Heterogeneous Treatment Effect Analysis
#
# Reads MICE-imputed PSM matched set (m=1, 24h) from 01_did_analysis.R.
# Performs:
#   1. Risk-based quintile (ΔCr) and tertile (AKI) HTE  [Kent et al. BMJ 2018]
#   2. Pre-specified clinical subgroup forest (binary splits, dual outcome)
#   3. Crossed-subgroup benefit + harm detection
#   4. Overall AKI staging summary (KDIGO ≥1 / ≥2 / ≥3)
#   5. Individual trajectory data for time-course spaghetti plots
#
# Usage: Rscript 01d_hte_mice.R eicu
#        Rscript 01d_hte_mice.R mimic
# ============================================================================

suppressPackageStartupMessages({ library(sandwich); library(lmtest) })

RESULTS <- path.expand("~/mg_aki/results")
WINDOW  <- 6

# ── Helpers ────────────────────────────────────────────────────────────────

did_dr <- function(df, covars = NULL, outcome = "delta_cr") {
  if (!outcome %in% names(df)) return(NULL)
  d <- df[!is.na(df[[outcome]]),]
  nt <- sum(d$treated == 1); nc <- sum(d$treated == 0)
  if (nt < 10 || nc < 10) return(NULL)

  adj <- character(0)
  if (length(covars) > 0) {
    avail <- intersect(covars, names(d))
    for (v in avail) {
      x1 <- d[[v]][d$treated == 1]; x0 <- d[[v]][d$treated == 0]
      sp <- sqrt((var(x1, na.rm = T) + var(x0, na.rm = T)) / 2)
      if (!is.na(sp) && sp > 1e-10) {
        smd <- abs(mean(x1, na.rm = T) - mean(x0, na.rm = T)) / sp
        if (!is.na(smd) && smd > 0.05) adj <- c(adj, v)
      }
    }
    adj <- adj[vapply(adj, function(v) var(d[[v]], na.rm = T) > 1e-10, logical(1))]
  }

  fml <- if (length(adj) > 0)
    as.formula(paste(outcome, "~ treated +", paste(adj, collapse = "+")))
  else as.formula(paste(outcome, "~ treated"))

  fit <- tryCatch(lm(fml, data = d), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  ct <- tryCatch(coeftest(fit, vcov. = vcovHC(fit, type = "HC1")),
                 error = function(e) tryCatch(coeftest(fit), error = function(e2) NULL))
  if (is.null(ct) || !"treated" %in% rownames(ct)) return(NULL)

  est <- ct["treated", "Estimate"]; se <- ct["treated", "Std. Error"]
  list(did = est, se = se, p = ct["treated", "Pr(>|t|)"],
       ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se,
       n_trt = nt, n_ctl = nc)
}

aki_rd <- function(df, outcome = "aki1") {
  if (!outcome %in% names(df)) return(NULL)
  d <- df[!is.na(df[[outcome]]),]
  nt <- sum(d$treated == 1); nc <- sum(d$treated == 0)
  if (nt < 10 || nc < 10) return(NULL)

  r1 <- mean(d[[outcome]][d$treated == 1])
  r0 <- mean(d[[outcome]][d$treated == 0])
  rd <- r1 - r0

  fit <- tryCatch(lm(as.formula(paste(outcome, "~ treated")), data = d),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  ct <- tryCatch(coeftest(fit, vcov. = vcovHC(fit, type = "HC1")),
                 error = function(e) tryCatch(coeftest(fit), error = function(e2) NULL))
  if (is.null(ct) || !"treated" %in% rownames(ct)) return(NULL)

  se <- ct["treated", "Std. Error"]; p <- ct["treated", "Pr(>|t|)"]
  nnt <- if (rd < 0 && p < 0.05) round(1 / abs(rd)) else NA
  nnh <- if (rd > 0 && p < 0.05) round(1 / abs(rd)) else NA
  list(rate_trt = r1, rate_ctl = r0, rd = rd, se = se, p = p,
       ci_lo = rd - 1.96 * se, ci_hi = rd + 1.96 * se,
       nnt = nnt, nnh = nnh, n_trt = nt, n_ctl = nc)
}

build_dcr <- function(pids, cr_all, target_h) {
  cr <- cr_all[cr_all$pid %in% pids, ]
  pre  <- cr[cr$offset_h >= 0 & cr$offset_h <= 6, ]
  pre  <- pre[order(pre$pid, pre$offset_h), ]; pre <- pre[!duplicated(pre$pid), ]
  post <- cr[cr$offset_h >= (target_h - WINDOW) & cr$offset_h <= (target_h + WINDOW), ]
  post$dist <- abs(post$offset_h - target_h)
  post <- post[order(post$pid, post$dist), ]; post <- post[!duplicated(post$pid), ]
  m <- merge(pre[, c("pid","labresult","offset_h")],
             post[, c("pid","labresult","offset_h")],
             by = "pid", suffixes = c("_pre","_post"))
  m <- m[m$offset_h_post > m$offset_h_pre, ]
  m$delta_cr <- m$labresult_post - m$labresult_pre
  m$cr_ratio <- m$labresult_post / m$labresult_pre
  m
}

interaction_p <- function(dat, moderator, outcome = "delta_cr") {
  d <- dat[!is.na(dat[[outcome]]) & !is.na(dat[[moderator]]), ]
  if (nrow(d) < 40) return(NA)
  fml <- as.formula(paste(outcome, "~ treated *", moderator))
  fit <- tryCatch(lm(fml, data = d), error = function(e) NULL)
  if (is.null(fit)) return(NA)
  ct <- tryCatch(coeftest(fit, vcov. = vcovHC(fit, type = "HC1")),
                 error = function(e) tryCatch(coeftest(fit), error = function(e2) NULL))
  if (is.null(ct)) return(NA)
  ir <- grep(paste0("treated:", moderator), rownames(ct))
  if (length(ir) == 0) ir <- grep(paste0(moderator, ":treated"), rownames(ct))
  if (length(ir) == 0) return(NA)
  ct[ir[1], "Pr(>|t|)"]
}

fmt_sig  <- function(p) if (!is.na(p) && p < 0.05) " *" else "  "
fmt_nnt  <- function(a) {
  if (is.null(a)) return("\u2014")
  if (!is.na(a$nnt)) return(as.character(a$nnt))
  if (!is.na(a$nnh)) return(paste0("NNH=", a$nnh))
  "\u2014"
}

# ============================================================================
# MAIN
# ============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) { cat("Usage: Rscript 01d_hte_mice.R <db>\n"); quit(status = 1) }
db <- toupper(args[1]); tag <- tolower(db)

SEP  <- paste(rep("=", 70), collapse = "")
SEP2 <- paste(rep("-", 70), collapse = "")
cat(sprintf("\n%s\n01d_hte_mice.R \u2014 HTE Analysis: %s\n%s\n", SEP, db, SEP))

# ── Load matched set + Cr time series ─────────────────────────────────────
matched <- read.csv(file.path(RESULTS, sprintf("did_matched_%s_mice_24h.csv", tag)),
                    stringsAsFactors = FALSE)
cr_all  <- read.csv(file.path(RESULTS, sprintf("did_cr_all_%s.csv", tag)),
                    stringsAsFactors = FALSE)

cr_id <- if ("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
cr_all$pid <- cr_all[[cr_id]]
if (!"labresultoffset" %in% names(cr_all)) cr_all$labresultoffset <- cr_all$offset_min
cr_all$offset_h <- cr_all$labresultoffset / 60

n_trt <- sum(matched$treated == 1); n_ctl <- sum(matched$treated == 0)
cat(sprintf("  Matched set: %d treated + %d control\n", n_trt, n_ctl))

# ── PS covariates available for DR adjustment ─────────────────────────────
ps_all <- c("age","is_female","bmi","surg_cabg","surg_valve","surg_combined",
            "heart_failure","hypertension","diabetes","ckd","copd","pvd",
            "stroke","liver_disease","egfr","first_heartrate",
            "first_potassium","first_calcium","first_lactate",
            "lactate_missing","first_mg_value")
ps_avail <- intersect(ps_all, names(matched))

# ── Compute AKI staging (need Cr_pre / Cr_post for ratio) ────────────────
cat(sprintf("\n%s\n  AKI STAGING\n%s\n", SEP2, SEP2))

cr24 <- build_dcr(matched$pid, cr_all, 24)
matched <- merge(matched, cr24[, c("pid","cr_ratio","labresult_pre","labresult_post")],
                 by = "pid", all.x = TRUE)

matched$aki1 <- with(matched, ifelse(
  is.na(delta_cr), NA_integer_,
  as.integer(delta_cr >= 0.3 | (!is.na(cr_ratio) & cr_ratio >= 1.5))))
matched$aki2 <- with(matched, ifelse(
  is.na(delta_cr), NA_integer_,
  as.integer(!is.na(cr_ratio) & cr_ratio >= 2.0)))
matched$aki3 <- with(matched, ifelse(
  is.na(delta_cr), NA_integer_,
  as.integer((!is.na(cr_ratio) & cr_ratio >= 3.0) |
             (!is.na(labresult_post) & labresult_post >= 4.0))))

for (stg in c("aki1","aki2","aki3")) {
  a <- aki_rd(matched, stg)
  lab <- c(aki1 = "KDIGO \u2265 Stage 1", aki2 = "KDIGO \u2265 Stage 2",
           aki3 = "KDIGO \u2265 Stage 3")[stg]
  if (!is.null(a))
    cat(sprintf("  %s: trt=%.1f%%, ctl=%.1f%%, RD=%+.1f%% [%+.1f,%+.1f], P=%.4f, %s\n",
                lab, 100*a$rate_trt, 100*a$rate_ctl, 100*a$rd,
                100*a$ci_lo, 100*a$ci_hi, a$p, fmt_nnt(a)))
}

# ======================================================================== #
# SECTION 1 — Risk-Based Stratified HTE
# ======================================================================== #
cat(sprintf("\n%s\n  SECTION 1: Risk-Based HTE (Kent et al.)\n%s\n", SEP, SEP))

risk_vars <- intersect(c("age","egfr","diabetes","ckd","heart_failure","hypertension",
                          "surg_cabg","surg_combined","bmi","first_lactate","first_mg_value"),
                       names(matched))
risk_vars <- risk_vars[vapply(risk_vars, function(v) {
  x <- matched[[v]][matched$treated == 0]
  !all(is.na(x)) && var(x, na.rm = T) > 1e-10
}, logical(1))]

controls_aki <- matched[matched$treated == 0 & !is.na(matched$aki1), ]
risk_fml <- as.formula(paste("aki1 ~", paste(risk_vars, collapse = "+")))
risk_fit <- tryCatch(glm(risk_fml, data = controls_aki, family = binomial()),
                     error = function(e) { cat("  Risk model failed\n"); NULL })

if (!is.null(risk_fit)) {
  matched$pred_risk <- predict(risk_fit, newdata = matched, type = "response")
  matched$pred_risk[is.na(matched$pred_risk)] <- median(matched$pred_risk, na.rm = T)

  # — Quintile cuts (for ΔCr) —
  q5_brk <- quantile(matched$pred_risk, probs = seq(0, 1, 0.2), na.rm = T)
  if (length(unique(q5_brk)) == 6) {
    matched$risk_q5 <- cut(matched$pred_risk, breaks = q5_brk,
                           labels = paste0("Q", 1:5), include.lowest = TRUE)
  } else {
    matched$risk_q5 <- cut(matched$pred_risk,
                           breaks = unique(q5_brk), include.lowest = TRUE)
    levels(matched$risk_q5) <- paste0("Q", seq_along(levels(matched$risk_q5)))
  }

  # — Tertile cuts (for AKI binary) —
  t3_brk <- quantile(matched$pred_risk, probs = seq(0, 1, 1/3), na.rm = T)
  if (length(unique(t3_brk)) == 4) {
    matched$risk_t3 <- cut(matched$pred_risk, breaks = t3_brk,
                           labels = paste0("T", 1:3), include.lowest = TRUE)
  } else {
    matched$risk_t3 <- cut(matched$pred_risk,
                           breaks = unique(t3_brk), include.lowest = TRUE)
    levels(matched$risk_t3) <- paste0("T", seq_along(levels(matched$risk_t3)))
  }

  # Report quintiles
  cat("\n  \u2500\u2500 Risk Quintiles (\u0394Cr) \u2500\u2500\n")
  cat("  Quintile  n_trt  n_ctl  DiD        95% CI               P       AKI1_trt AKI1_ctl\n")

  q5_rows <- list()
  for (q in levels(matched$risk_q5)) {
    sub <- matched[!is.na(matched$risk_q5) & matched$risk_q5 == q & !is.na(matched$delta_cr), ]
    r <- did_dr(sub, ps_avail); a <- aki_rd(sub, "aki1")
    if (!is.null(r)) {
      cat(sprintf("  %-8s  %5d  %5d  %+.4f   [%+.4f,%+.4f]  %.4f%s  %5.1f%%   %5.1f%%\n",
                  q, r$n_trt, r$n_ctl, r$did, r$ci_lo, r$ci_hi, r$p, fmt_sig(r$p),
                  if (!is.null(a)) 100*a$rate_trt else NA,
                  if (!is.null(a)) 100*a$rate_ctl else NA))
      q5_rows[[length(q5_rows)+1]] <- data.frame(
        stratum = q, n_trt = r$n_trt, n_ctl = r$n_ctl,
        did = r$did, did_se = r$se, did_p = r$p,
        did_ci_lo = r$ci_lo, did_ci_hi = r$ci_hi,
        aki1_trt = if (!is.null(a)) a$rate_trt else NA,
        aki1_ctl = if (!is.null(a)) a$rate_ctl else NA,
        aki1_rd  = if (!is.null(a)) a$rd else NA,
        aki1_p   = if (!is.null(a)) a$p else NA,
        stringsAsFactors = FALSE)
    }
  }

  ip_cr <- interaction_p(matched, "pred_risk", "delta_cr")
  cat(sprintf("\n  Interaction (treated \u00d7 pred_risk, \u0394Cr): P=%.4f\n", ip_cr))

  # Report tertiles
  cat("\n  \u2500\u2500 Risk Tertiles (AKI \u2265 Stage 1) \u2500\u2500\n")
  cat("  Tertile  n_trt  n_ctl  AKI_trt AKI_ctl  RD        95% CI               P       NNT\n")

  t3_rows <- list()
  for (tt in levels(matched$risk_t3)) {
    sub <- matched[!is.na(matched$risk_t3) & matched$risk_t3 == tt, ]
    a <- aki_rd(sub, "aki1"); r <- did_dr(sub, ps_avail)
    if (!is.null(a)) {
      cat(sprintf("  %-8s  %5d  %5d  %5.1f%%  %5.1f%%  %+.4f   [%+.4f,%+.4f]  %.4f%s  %s\n",
                  tt, a$n_trt, a$n_ctl, 100*a$rate_trt, 100*a$rate_ctl,
                  a$rd, a$ci_lo, a$ci_hi, a$p, fmt_sig(a$p), fmt_nnt(a)))
      t3_rows[[length(t3_rows)+1]] <- data.frame(
        stratum = tt, n_trt = a$n_trt, n_ctl = a$n_ctl,
        aki1_trt = a$rate_trt, aki1_ctl = a$rate_ctl,
        rd = a$rd, rd_se = a$se, rd_p = a$p,
        rd_ci_lo = a$ci_lo, rd_ci_hi = a$ci_hi,
        nnt = a$nnt, nnh = a$nnh,
        did = if (!is.null(r)) r$did else NA,
        did_p = if (!is.null(r)) r$p else NA,
        stringsAsFactors = FALSE)
    }
  }

  ip_aki <- interaction_p(matched, "pred_risk", "aki1")
  cat(sprintf("\n  Interaction (treated \u00d7 pred_risk, AKI): P=%.4f\n", ip_aki))

  if (length(q5_rows) > 0)
    write.csv(do.call(rbind, q5_rows),
              file.path(RESULTS, sprintf("hte_risk_q5_%s.csv", tag)), row.names = FALSE)
  if (length(t3_rows) > 0)
    write.csv(do.call(rbind, t3_rows),
              file.path(RESULTS, sprintf("hte_risk_t3_%s.csv", tag)), row.names = FALSE)
}

# ======================================================================== #
# SECTION 2 — Pre-Specified Clinical Subgroups (Binary Splits)
# ======================================================================== #
cat(sprintf("\n%s\n  SECTION 2: Clinical Subgroups\n%s\n", SEP, SEP))

subgroups <- list(
  list(name = "eGFR < 60",       var = "egfr",            op = "<",  val = 60),
  list(name = "eGFR < 90",       var = "egfr",            op = "<",  val = 90),
  list(name = "Mg < 1.8",        var = "first_mg_value",  op = "<",  val = 1.8),
  list(name = "Mg < 2.0",        var = "first_mg_value",  op = "<",  val = 2.0),
  list(name = "CABG",            var = "surg_cabg",       op = "==", val = 1),
  list(name = "Diabetes",        var = "diabetes",        op = "==", val = 1),
  list(name = "CKD",             var = "ckd",             op = "==", val = 1),
  list(name = "Heart failure",   var = "heart_failure",   op = "==", val = 1),
  list(name = "BMI >= 30",       var = "bmi",             op = ">=", val = 30),
  list(name = "Age >= 65",       var = "age",             op = ">=", val = 65)
)

cat("  Subgroup          Side   n_trt n_ctl  DiD       P_did   AKI_RD    P_aki   NNT    Int_P\n")
cat("  ────────────────  ────   ───── ─────  ────────  ──────  ────────  ──────  ─────  ──────\n")

sg_rows <- list()
for (sg in subgroups) {
  if (!sg$var %in% names(matched)) next

  # Build yes/no index
  v <- matched[[sg$var]]
  yes_idx <- switch(sg$op,
    "<"  = !is.na(v) & v <  sg$val,
    ">=" = !is.na(v) & v >= sg$val,
    "==" = !is.na(v) & v == sg$val)

  # Interaction P (compute once)
  matched$.sg_flag <- as.integer(yes_idx)
  ip <- interaction_p(matched, ".sg_flag", "delta_cr")

  for (side in c("yes", "no")) {
    idx <- if (side == "yes") yes_idx else (!is.na(v) & !yes_idx)
    sub <- matched[idx & !is.na(matched$delta_cr), ]
    r <- did_dr(sub, ps_avail); a <- aki_rd(sub, "aki1")
    if (is.null(r)) next

    label <- if (side == "yes") sg$name else paste0("NOT ", sg$name)
    ip_str <- if (side == "yes" && !is.na(ip)) sprintf("%.4f", ip) else ""

    cat(sprintf("  %-18s %-4s   %5d %5d  %+.4f  %.4f%s %+.4f  %.4f%s %-5s  %s\n",
                label, side, r$n_trt, r$n_ctl,
                r$did, r$p, fmt_sig(r$p),
                if (!is.null(a)) a$rd else NA,
                if (!is.null(a)) a$p  else NA,
                if (!is.null(a)) fmt_sig(a$p) else "  ",
                fmt_nnt(a), ip_str))

    sg_rows[[length(sg_rows)+1]] <- data.frame(
      subgroup = sg$name, side = side,
      n_trt = r$n_trt, n_ctl = r$n_ctl,
      did = r$did, did_se = r$se, did_p = r$p,
      did_ci_lo = r$ci_lo, did_ci_hi = r$ci_hi,
      aki_rd     = if (!is.null(a)) a$rd     else NA,
      aki_rd_se  = if (!is.null(a)) a$se     else NA,
      aki_rd_p   = if (!is.null(a)) a$p      else NA,
      aki_rd_cilo = if (!is.null(a)) a$ci_lo else NA,
      aki_rd_cihi = if (!is.null(a)) a$ci_hi else NA,
      nnt = if (!is.null(a)) a$nnt else NA,
      nnh = if (!is.null(a)) a$nnh else NA,
      interaction_p = ip,
      stringsAsFactors = FALSE)
  }
}
matched$.sg_flag <- NULL

if (length(sg_rows) > 0)
  write.csv(do.call(rbind, sg_rows),
            file.path(RESULTS, sprintf("hte_subgroups_%s.csv", tag)), row.names = FALSE)

# ======================================================================== #
# SECTION 3 — Crossed Subgroups: Benefit + Harm Detection
# ======================================================================== #
cat(sprintf("\n%s\n  SECTION 3: Crossed Subgroups (Benefit + Harm)\n%s\n", SEP, SEP))

crossed <- list(
  # -- Expected benefit phenotypes --
  list(name = "DM + BMI>=30 (metabolic)",
       cond = quote(diabetes == 1 & !is.na(bmi) & bmi >= 30)),
  list(name = "HF + CABG/combined (cardiac burden)",
       cond = quote(heart_failure == 1 & (surg_cabg == 1 | surg_combined == 1))),
  list(name = "CABG + eGFR<90",
       cond = quote(surg_cabg == 1 & !is.na(egfr) & egfr < 90)),
  list(name = "DM + CKD",
       cond = quote(diabetes == 1 & ckd == 1)),
  list(name = "Mg 2.0-2.3 (prophylactic)",
       cond = quote(!is.na(first_mg_value) & first_mg_value >= 2.0 & first_mg_value < 2.3)),
  # -- Expected harm / confounding-by-indication --
  list(name = "Mg<1.8 + CKD",
       cond = quote(!is.na(first_mg_value) & first_mg_value < 1.8 & ckd == 1)),
  list(name = "Mg<1.8 + eGFR<60",
       cond = quote(!is.na(first_mg_value) & first_mg_value < 1.8 & !is.na(egfr) & egfr < 60)),
  list(name = "Mg<1.8 + eGFR<45",
       cond = quote(!is.na(first_mg_value) & first_mg_value < 1.8 & !is.na(egfr) & egfr < 45)),
  list(name = "eGFR < 30",
       cond = quote(!is.na(egfr) & egfr < 30)),
  list(name = "eGFR < 45",
       cond = quote(!is.na(egfr) & egfr < 45)),
  list(name = "Mg>2.3 + CKD",
       cond = quote(!is.na(first_mg_value) & first_mg_value > 2.3 & ckd == 1))
)

cat("  Group                                n_trt n_ctl  DiD       P_did   AKI_RD    P_aki   AKI2_RD\n")
cat("  ──────────────────────────────────── ───── ─────  ──────── ──────  ──────── ──────  ────────\n")

cr_rows <- list()
for (cg in crossed) {
  idx <- tryCatch(eval(cg$cond, matched), error = function(e) rep(FALSE, nrow(matched)))
  idx[is.na(idx)] <- FALSE
  sub <- matched[idx & !is.na(matched$delta_cr), ]
  r  <- did_dr(sub, ps_avail)
  a1 <- aki_rd(sub, "aki1")
  a2 <- aki_rd(sub, "aki2")

  if (!is.null(r)) {
    cat(sprintf("  %-40s %5d %5d  %+.4f  %.4f%s %+.4f  %.4f%s %+.4f\n",
                cg$name, r$n_trt, r$n_ctl,
                r$did, r$p, fmt_sig(r$p),
                if (!is.null(a1)) a1$rd else NA,
                if (!is.null(a1)) a1$p  else NA,
                if (!is.null(a1)) fmt_sig(a1$p) else "  ",
                if (!is.null(a2)) a2$rd else NA))
    cr_rows[[length(cr_rows)+1]] <- data.frame(
      group = cg$name, n_trt = r$n_trt, n_ctl = r$n_ctl,
      did = r$did, did_se = r$se, did_p = r$p,
      did_ci_lo = r$ci_lo, did_ci_hi = r$ci_hi,
      aki1_rd = if (!is.null(a1)) a1$rd else NA,
      aki1_p  = if (!is.null(a1)) a1$p  else NA,
      aki2_rd = if (!is.null(a2)) a2$rd else NA,
      aki2_p  = if (!is.null(a2)) a2$p  else NA,
      stringsAsFactors = FALSE)
  } else {
    n_t <- sum(idx & matched$treated == 1, na.rm = T)
    n_c <- sum(idx & matched$treated == 0, na.rm = T)
    cat(sprintf("  %-40s  (too few: trt=%d, ctl=%d)\n", cg$name, n_t, n_c))
  }
}

if (length(cr_rows) > 0)
  write.csv(do.call(rbind, cr_rows),
            file.path(RESULTS, sprintf("hte_crossed_%s.csv", tag)), row.names = FALSE)

# ======================================================================== #
# SECTION 4 — Time-Course Individual Trajectories
# ======================================================================== #
cat(sprintf("\n%s\n  SECTION 4: Individual Trajectories for Spaghetti Plots\n%s\n", SEP, SEP))

time_pts <- seq(6, 36, by = 6)

# Columns to carry over for plotting / faceting
carry_cols <- intersect(
  c("pid","treated","risk_q5","risk_t3","pred_risk",
    "first_mg_value","egfr","surg_cabg","diabetes","ckd",
    "heart_failure","bmi","aki1"),
  names(matched))
carry <- matched[, carry_cols]

traj_list <- list()
for (th in time_pts) {
  cr_th <- build_dcr(matched$pid, cr_all, th)
  if (nrow(cr_th) == 0) next
  cr_th$target_h <- th
  cr_th <- merge(cr_th[, c("pid","delta_cr","cr_ratio","target_h")], carry, by = "pid")
  traj_list[[length(traj_list)+1]] <- cr_th
}

if (length(traj_list) > 0) {
  traj <- do.call(rbind, traj_list)
  write.csv(traj, file.path(RESULTS, sprintf("hte_trajectories_%s.csv", tag)),
            row.names = FALSE)
  cat(sprintf("  Saved %d patient-timepoint rows (%d unique patients, %d time points)\n",
              nrow(traj), length(unique(traj$pid)), length(time_pts)))

  cat("\n  target_h  n_obs  mean_trt   mean_ctl   diff\n")
  cat("  ────────  ─────  ────────   ────────   ────────\n")
  for (th in time_pts) {
    ts <- traj[traj$target_h == th, ]
    m1 <- mean(ts$delta_cr[ts$treated == 1], na.rm = T)
    m0 <- mean(ts$delta_cr[ts$treated == 0], na.rm = T)
    cat(sprintf("  %4dh     %5d  %+.4f    %+.4f    %+.4f\n", th, nrow(ts), m1, m0, m1 - m0))
  }
}

# ======================================================================== #
# DONE
# ======================================================================== #
cat(sprintf("\n%s\n01d_hte_mice.R \u2014 %s COMPLETE\n%s\n", SEP, db, SEP))

# Output manifest
cat("\n  Files written:\n")
for (f in c(sprintf("hte_risk_q5_%s.csv", tag),
            sprintf("hte_risk_t3_%s.csv", tag),
            sprintf("hte_subgroups_%s.csv", tag),
            sprintf("hte_crossed_%s.csv", tag),
            sprintf("hte_trajectories_%s.csv", tag))) {
  fp <- file.path(RESULTS, f)
  if (file.exists(fp)) cat(sprintf("    %s  (%s)\n", f, format(file.size(fp), big.mark = ",")))
}
