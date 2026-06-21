#!/usr/bin/env Rscript
# ============================================================================
# 02_psm.R — Risk-Set PSM with Temporal Alignment (v5 Design)
#
# Primary endpoint: ΔCr at T₀ + 24h
# Lab timing: FIRST (admission baseline)
# MICE m=20, single run
# Control pool: Mg-free through T₀ + 24h
# Effect monitoring: 6–36h (contamination check for >24h)
#
# Usage: Rscript 02_psm.R eicu
#        Rscript 02_psm.R mimic
# ============================================================================

suppressPackageStartupMessages({ library(sandwich); library(lmtest); library(mice) })

RESULTS    <- path.expand("~/mg_aki/results")
PRIMARY_H  <- 24       # primary: ΔCr at T₀ + 24h
CR_WINDOW  <- 12       # ±12h tolerance for Cr lookup
CALIPER_SD <- 0.2
M_IMP      <- 20
TARGETS    <- c(6, 12, 18, 24, 30, 36)

# PS covariates
PS_FIXED <- c("age","is_female","bmi",
              "surg_cabg","surg_valve","surg_combined",
              "heart_failure","hypertension","diabetes","ckd",
              "copd","pvd","stroke","liver_disease","egfr",
              "ppi_chronic","loop_diuretic_chronic","acei_arb_chronic","nsaid_chronic")
PS_LAB_BASES <- c("magnesium","potassium","calcium","lactate","heartrate")

# ═══════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════

find_cr <- function(cr_pt, target_h, window = CR_WINDOW) {
  if (is.null(cr_pt) || nrow(cr_pt) == 0) return(c(NA, NA))
  cand <- cr_pt[cr_pt$offset_h >= (target_h - window) &
                cr_pt$offset_h <= (target_h + window), ]
  if (nrow(cand) == 0) return(c(NA, NA))
  best <- cand[which.min(abs(cand$offset_h - target_h)), ]
  c(best$labresult, best$offset_h)
}

find_cr_pre <- function(cr_pt, t_h) {
  if (is.null(cr_pt) || nrow(cr_pt) == 0) return(c(NA, NA))
  cand <- cr_pt[cr_pt$offset_h >= 0 & cr_pt$offset_h < t_h, ]
  if (nrow(cand) == 0) return(c(NA, NA))
  best <- cand[which.max(cand$offset_h), ]
  c(best$labresult, best$offset_h)
}

safe_coeftest <- function(fit) {
  ct <- tryCatch(suppressWarnings(coeftest(fit, vcov. = vcovHC(fit, type = "HC1"))),
                 error = function(e) NULL)
  if (!is.null(ct) && is.matrix(ct) && "treated" %in% rownames(ct) &&
      ncol(ct) >= 4 && !any(is.nan(ct["treated", ]))) return(ct)
  tryCatch(coeftest(fit), error = function(e) NULL)
}

# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) { cat("Usage: Rscript 02_psm.R <db>\n"); quit(status = 1) }
db <- toupper(args[1]); tag <- tolower(db)

SEP <- paste(rep("=", 70), collapse = "")
cat(sprintf("\n%s\n02_psm.R — Risk-Set PSM: %s\n  Primary: ΔCr at T₀+%dh | Labs: FIRST | MICE m=%d\n%s\n",
            SEP, db, PRIMARY_H, M_IMP, SEP))

# ── Load data ─────────────────────────────────────────────────────
all_pts <- read.csv(file.path(RESULTS, sprintf("did_all_%s.csv", tag)), stringsAsFactors = FALSE)
cr_all  <- read.csv(file.path(RESULTS, sprintf("did_cr_all_%s.csv", tag)), stringsAsFactors = FALSE)

# Normalize cr_all pid + offset
cr_id <- if ("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
cr_all$pid <- cr_all[[cr_id]]
if (!"offset_h" %in% names(cr_all)) cr_all$offset_h <- cr_all$labresultoffset / 60
cr_all <- cr_all[order(cr_all$pid, cr_all$offset_h), ]

# Load labs (electrolytes + downsampled HR)
cat("  Loading labs...\n")
labs_raw <- read.csv(file.path(RESULTS, sprintf("did_labs_all_%s.csv", tag)), stringsAsFactors = FALSE)
pid_col_lab <- if ("patientunitstayid" %in% names(labs_raw)) "patientunitstayid" else "stay_id"
if (pid_col_lab %in% names(labs_raw)) labs_raw$pid <- labs_raw[[pid_col_lab]]
labs_elec <- labs_raw[labs_raw$lab_name %in% c("magnesium","potassium","calcium","lactate"), ]
labs_hr   <- labs_raw[labs_raw$lab_name == "heartrate", ]
if (nrow(labs_hr) > 500000) {
  labs_hr$hour_bin <- floor(labs_hr$offset_h)
  labs_hr <- labs_hr[order(labs_hr$pid, labs_hr$offset_h), ]
  labs_hr <- labs_hr[!duplicated(paste(labs_hr$pid, labs_hr$hour_bin)), ]
  labs_hr$hour_bin <- NULL
  cat(sprintf("    HR downsampled: %d → %d\n", sum(labs_raw$lab_name=="heartrate"), nrow(labs_hr)))
}
labs <- rbind(labs_elec, labs_hr)
rm(labs_raw, labs_elec, labs_hr); gc()

N <- nrow(all_pts)
cat(sprintf("  Patients: %d (%d treated, %d control)\n", N,
            sum(all_pts$treated == 1, na.rm=TRUE),
            sum(all_pts$treated == 0, na.rm=TRUE)))

# ── Pre-compute FIRST lab values ──────────────────────────────────
cat("  Computing FIRST lab values...\n")
for (ln in PS_LAB_BASES) {
  sub <- labs[labs$lab_name == ln, ]
  if (nrow(sub) == 0) next
  sub$mg_offset_h <- all_pts$mg_offset_h[match(sub$pid, all_pts$pid)]
  sub <- sub[sub$offset_h >= 0 & (is.na(sub$mg_offset_h) | sub$offset_h < sub$mg_offset_h), ]
  if (nrow(sub) == 0) next
  s <- sub[order(sub$offset_h), ]
  s <- s[!duplicated(s$pid), ]
  idx <- match(all_pts$pid, s$pid)
  all_pts[[paste0("first_",ln)]] <- s$value[idx]
  nf <- sum(!is.na(all_pts[[paste0("first_",ln)]]))
  cat(sprintf("    %s: %d (%.0f%%)\n", ln, nf, 100*nf/N))
}
all_pts$first_lactate_missing <- as.integer(is.na(all_pts$first_lactate))

# ── Treatment indices ─────────────────────────────────────────────
trt_idx <- which(all_pts$treated == 1 & !is.na(all_pts$mg_offset_h))
n_trt <- length(trt_idx)
cat(sprintf("  Treated with valid t_mg: %d\n", n_trt))

# ── Risk sets ─────────────────────────────────────────────────────
cat("  Pre-computing risk sets...\n")

cr_list <- split(cr_all[, c("labresult","offset_h")], cr_all$pid)

earliest_cr <- sapply(cr_list, function(x) min(x$offset_h, na.rm = TRUE))
all_pts$earliest_cr_h <- earliest_cr[as.character(all_pts$pid)]
all_pts$earliest_cr_h[is.na(all_pts$earliest_cr_h)] <- Inf

cat("    Computing first AKI time...\n")
first_aki <- sapply(cr_list, function(cr_pt) {
  if (nrow(cr_pt) < 2) return(NA_real_)
  cr_pt <- cr_pt[order(cr_pt$offset_h), ]
  bl <- cr_pt$labresult[1]
  if (is.na(bl) || bl <= 0) return(NA_real_)
  for (k in 2:nrow(cr_pt)) {
    d <- cr_pt$labresult[k] - bl
    r <- cr_pt$labresult[k] / bl
    if (!is.na(d) && (d >= 0.3 || (!is.na(r) && r >= 1.5)))
      return(cr_pt$offset_h[k])
  }
  NA_real_
})
all_pts$first_aki_h <- first_aki[as.character(all_pts$pid)]

cat(sprintf("    Building risk sets (Mg-free through t_mg + %dh)...\n", PRIMARY_H))
trt_pids <- all_pts$pid[trt_idx]
trt_tmg  <- all_pts$mg_offset_h[trt_idx]
risk_sets <- vector("list", n_trt)

for (k in seq_len(n_trt)) {
  t_mg <- trt_tmg[k]
  risk_sets[[k]] <- which(
    all_pts$icu_discharge_h > t_mg &
    (is.na(all_pts$first_aki_h) | all_pts$first_aki_h > t_mg) &
    all_pts$earliest_cr_h <= t_mg &
    (is.na(all_pts$mg_offset_h) | all_pts$mg_offset_h > t_mg + PRIMARY_H) &
    all_pts$pid != trt_pids[k])
}

rs_sizes <- vapply(risk_sets, length, integer(1))
cat(sprintf("    Risk sets: median=%.0f, IQR=[%.0f,%.0f], empty=%d\n",
            median(rs_sizes), quantile(rs_sizes, 0.25), quantile(rs_sizes, 0.75),
            sum(rs_sizes == 0)))

# ── PS covariates ─────────────────────────────────────────────────
lab_cols <- paste0("first_", PS_LAB_BASES)
ps_vars <- c(intersect(PS_FIXED, names(all_pts)),
             intersect(c(lab_cols, "first_lactate_missing"), names(all_pts)))
ps_vars <- ps_vars[vapply(ps_vars, function(v) {
  x <- all_pts[[v]]; !all(is.na(x)) && var(x, na.rm = TRUE) > 1e-10
}, logical(1))]
cat(sprintf("  PS covariates (%d): %s\n", length(ps_vars), paste(ps_vars, collapse=", ")))

# ── MICE ──────────────────────────────────────────────────────────
to_impute <- ps_vars[vapply(ps_vars, function(v) any(is.na(all_pts[[v]])), logical(1))]
cat(sprintf("  MICE m=%d on %d vars: %s\n", M_IMP, length(to_impute), paste(to_impute, collapse=", ")))

if (length(to_impute) > 0) {
  mice_df <- all_pts[, c("treated", ps_vars)]
  meth <- rep("", ncol(mice_df)); names(meth) <- names(mice_df)
  for (v in to_impute) meth[v] <- "pmm"
  imp <- mice(mice_df, m = M_IMP, method = meth, printFlag = FALSE, maxit = 10)
  cat(sprintf("  MICE done. Logged events: %d\n", nrow(imp$loggedEvents)))
} else {
  imp <- NULL
}

# ── Per-imputation matching + outcome ─────────────────────────────
cat("  Matching...\n")
all_dids <- list()
m1_match_trt <- NULL; m1_match_ctl <- NULL; m1_matched <- NULL; m1_d <- NULL

for (m_idx in 1:M_IMP) {
  if (m_idx %% 5 == 1) cat(sprintf("  m=%d...", m_idx))

  d <- all_pts
  if (!is.null(imp)) {
    imp_df <- complete(imp, m_idx)
    for (v in to_impute) d[[v]] <- imp_df[[v]]
  }

  ps_fml <- as.formula(paste("treated ~", paste(ps_vars, collapse = "+")))
  ps_fit <- tryCatch(glm(ps_fml, data = d, family = binomial()),
                     error = function(e) { cat("PS failed\n"); return(NULL) })
  if (is.null(ps_fit)) next
  d$ps <- predict(ps_fit, type = "response")
  caliper <- CALIPER_SD * sd(d$ps, na.rm = TRUE)

  match_trt <- integer(n_trt); match_ctl <- integer(n_trt)
  match_ps  <- numeric(n_trt); matched <- logical(n_trt)

  for (k in seq_len(n_trt)) {
    rs <- risk_sets[[k]]
    if (length(rs) == 0) next
    ps_i <- d$ps[trt_idx[k]]
    ps_dist <- abs(d$ps[rs] - ps_i)
    within_cal <- which(ps_dist <= caliper)
    if (length(within_cal) == 0) next
    best <- within_cal[which.min(ps_dist[within_cal])]
    match_trt[k] <- trt_idx[k]; match_ctl[k] <- rs[best]
    match_ps[k] <- ps_dist[best]; matched[k] <- TRUE
  }

  n_matched <- sum(matched)
  if (n_matched < 50) next
  if (m_idx == 1) {
    m1_match_trt <- match_trt; m1_match_ctl <- match_ctl
    m1_matched <- matched; m1_d <- d
  }

  for (target_h in TARGETS) {
    dcr_trt <- dcr_ctl <- numeric(n_matched)
    valid <- logical(n_matched); idx <- 0
    for (kk in which(matched)) {
      idx <- idx + 1
      tpid <- as.character(all_pts$pid[match_trt[kk]])
      cpid <- as.character(all_pts$pid[match_ctl[kk]])
      t_mg <- trt_tmg[kk]

      # Contamination check: control received Mg during outcome window?
      ctl_mg_h <- all_pts$mg_offset_h[match_ctl[kk]]
      if (!is.na(ctl_mg_h) && ctl_mg_h < t_mg + target_h) {
        valid[idx] <- FALSE; next
      }

      pre_t <- find_cr_pre(cr_list[[tpid]], t_mg)
      pre_c <- find_cr_pre(cr_list[[cpid]], t_mg)
      post_t <- find_cr(cr_list[[tpid]], t_mg + target_h)
      post_c <- find_cr(cr_list[[cpid]], t_mg + target_h)

      if (any(is.na(c(pre_t[1], pre_c[1], post_t[1], post_c[1])))) {
        valid[idx] <- FALSE; next
      }
      dcr_trt[idx] <- post_t[1] - pre_t[1]
      dcr_ctl[idx] <- post_c[1] - pre_c[1]
      valid[idx] <- TRUE
    }

    n_valid <- sum(valid)
    if (n_valid < 30) {
      all_dids[[length(all_dids)+1]] <- data.frame(
        m=m_idx, target_h=target_h, n=n_valid,
        did=NA, se=NA, n_matched=n_matched, stringsAsFactors=FALSE)
      next
    }

    pair_df <- data.frame(
      delta_cr = c(dcr_trt[valid], dcr_ctl[valid]),
      treated = rep(c(1,0), each = n_valid))
    fit <- lm(delta_cr ~ treated, data = pair_df)
    ct <- safe_coeftest(fit)
    est <- if (!is.null(ct) && "treated" %in% rownames(ct)) ct["treated",1] else coef(fit)["treated"]
    se  <- if (!is.null(ct) && "treated" %in% rownames(ct)) ct["treated",2] else NA

    all_dids[[length(all_dids)+1]] <- data.frame(
      m=m_idx, target_h=target_h, n=n_valid,
      did=est, se=se, n_matched=n_matched, stringsAsFactors=FALSE)
  }
}
cat("\n")

res <- do.call(rbind, all_dids)

# ── Rubin's rules ─────────────────────────────────────────────────
cat(sprintf("\n  ── RESULTS (FIRST labs, primary=%dh) ──\n", PRIMARY_H))
cat("  target   DiD        SE       P        95%% CI                n    FMI\n")
cat("  ──────   ────────   ──────   ──────   ────────────────────   ───  ─────\n")

pooled <- list()
for (th in TARGETS) {
  sub <- res[res$target_h == th & !is.na(res$did) & !is.na(res$se), ]
  if (nrow(sub) < 2) { cat(sprintf("  %4dh    (insufficient data)\n", th)); next }
  Q <- mean(sub$did); U <- mean(sub$se^2); B <- var(sub$did)
  m_v <- nrow(sub)
  T_var <- U + (1 + 1/m_v) * B; T_se <- sqrt(T_var)
  lam <- ((1 + 1/m_v) * B) / T_var
  df <- max(3, (m_v - 1) / lam^2)
  p <- 2 * pt(abs(Q/T_se), df = df, lower.tail = FALSE)
  ci_lo <- Q - qt(0.975, df) * T_se; ci_hi <- Q + qt(0.975, df) * T_se
  sig <- if (!is.na(p) && p < 0.05) " *" else "  "
  pri <- if (th == PRIMARY_H) " << PRIMARY" else ""
  n_avg <- mean(sub$n)
  cat(sprintf("  %4dh    %+.4f   %.4f   %.4f%s [%+.4f,%+.4f]   %3.0f  %.3f%s\n",
              th, Q, T_se, p, sig, ci_lo, ci_hi, n_avg, lam, pri))
  pooled[[length(pooled)+1]] <- data.frame(
    target_h=th, did=Q, se=T_se, p=p, ci_lo=ci_lo, ci_hi=ci_hi,
    fmi=lam, n_avg=n_avg, stringsAsFactors=FALSE)
}

# ── Balance ───────────────────────────────────────────────────────
if (!is.null(m1_d) && !is.null(m1_matched) && sum(m1_matched) > 50) {
  cat(sprintf("\n  Matched (m=1): %d/%d (%.1f%%)\n",
              sum(m1_matched), n_trt, 100*sum(m1_matched)/n_trt))
  cat("\n  ── COVARIATE BALANCE (m=1) ──\n")
  cat("  Covariate                raw_SMD  matched_SMD  status\n")
  cat("  ───────────────────────  ───────  ───────────  ──────\n")

  trt_all <- which(m1_d$treated == 1); ctl_all <- which(m1_d$treated == 0)
  trt_m <- m1_match_trt[m1_matched]; ctl_m <- m1_match_ctl[m1_matched]
  n_viol <- 0

  for (v in ps_vars) {
    x1_raw <- m1_d[[v]][trt_all]; x0_raw <- m1_d[[v]][ctl_all]
    x1_mat <- m1_d[[v]][trt_m];   x0_mat <- m1_d[[v]][ctl_m]
    sp_raw <- sqrt((var(x1_raw, na.rm=T) + var(x0_raw, na.rm=T)) / 2)
    sp_mat <- sqrt((var(x1_mat, na.rm=T) + var(x0_mat, na.rm=T)) / 2)
    smd_raw <- if (!is.na(sp_raw) && sp_raw > 1e-10)
                 abs(mean(x1_raw,na.rm=T) - mean(x0_raw,na.rm=T)) / sp_raw else NA
    smd_mat <- if (!is.na(sp_mat) && sp_mat > 1e-10)
                 abs(mean(x1_mat,na.rm=T) - mean(x0_mat,na.rm=T)) / sp_mat else NA
    flag <- if (!is.na(smd_mat) && smd_mat > 0.1) { n_viol <- n_viol+1; "VIOL" } else "ok"
    cat(sprintf("  %-25s  %.3f    %.3f       %s\n", v,
                ifelse(is.na(smd_raw), NA, smd_raw),
                ifelse(is.na(smd_mat), NA, smd_mat), flag))
  }
  cat(sprintf("\n  Max matched SMD: %.3f | Violations (>0.1): %d/%d\n",
              max(vapply(ps_vars, function(v) {
                x1 <- m1_d[[v]][trt_m]; x0 <- m1_d[[v]][ctl_m]
                sp <- sqrt((var(x1,na.rm=T)+var(x0,na.rm=T))/2)
                if (!is.na(sp) && sp > 1e-10) abs(mean(x1,na.rm=T)-mean(x0,na.rm=T))/sp else 0
              }, numeric(1))),
              n_viol, length(ps_vars)))

  write.csv(do.call(rbind, pooled),
            file.path(RESULTS, sprintf("did_riskset_%s.csv", tag)), row.names = FALSE)
}

cat(sprintf("\n%s\n02_psm.R — %s COMPLETE\n  Output: did_riskset_%s.csv\n%s\n",
            SEP, db, tag, SEP))
