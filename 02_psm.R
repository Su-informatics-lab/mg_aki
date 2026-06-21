#!/usr/bin/env Rscript
# ============================================================================
# 02_psm.R — Canonical Risk-Set PSM (v5 Final)
#
# Specs:
#   PRIMARY:     21 vars (base 15 + last K/Mg/Ca/lactate/lactate_miss/HR)
#   SENSITIVITY: 19 vars (drop last_K and last_Mg — treatment-trigger labs)
#
# Pools:    yet-untreated (sequential trial) | never-treated (parallel trial)
# Methods:  PSM plain | PSM + DR (adjust SMD > 0.1 covariates)
# Horizons: 6, 12, 18, 24, 30, 36, 42, 48h
#
# Lab timing: LAST (closest to T₀) — columns named last_* accordingly
# MICE: m=20, averaged → single PS → single match
# Matching: 1:1 with replacement, caliper 0.2 SD, HC1 SE
#
# Usage: Rscript 02_psm.R eicu
#        Rscript 02_psm.R mimic
# ============================================================================

suppressPackageStartupMessages({ library(sandwich); library(lmtest); library(mice) })

RESULTS    <- path.expand("~/mg_aki/results")
PRIMARY_H  <- 24
CR_WINDOW  <- 12
CALIPER_SD <- 0.2
M_IMP      <- 20
TARGETS    <- c(6, 12, 18, 24, 30, 36, 42, 48)

# ── Covariate specs ───────────────────────────────────────────────
PS_BASE <- c("age","is_female","bmi",
             "surg_cabg","surg_valve","surg_combined",
             "heart_failure","hypertension","diabetes","ckd",
             "copd","pvd","stroke","liver_disease","egfr")

PS_LABS_ALL <- c("last_magnesium","last_potassium","last_calcium",
                 "last_lactate","last_lactate_missing","last_heartrate")

PS_LABS_REDUCED <- c("last_calcium",
                     "last_lactate","last_lactate_missing","last_heartrate")

SPECS <- list(
  primary     = list(vars = c(PS_BASE, PS_LABS_ALL),     label = "21var (all labs)"),
  sensitivity = list(vars = c(PS_BASE, PS_LABS_REDUCED), label = "19var (no K/Mg)")
)

LAB_BASES <- c("magnesium","potassium","calcium","lactate","heartrate")

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
# run_spec_pool — one covariate spec x one control pool
# ═══════════════════════════════════════════════════════════════════

run_spec_pool <- function(spec_name, spec_obj, pool_name,
                          all_pts, trt_idx, risk_sets,
                          cr_list, trt_tmg, caliper) {

  ps_vars <- intersect(spec_obj$vars, names(all_pts))
  ps_vars <- ps_vars[vapply(ps_vars, function(v) {
    x <- all_pts[[v]]; !all(is.na(x)) && var(x, na.rm=TRUE) > 1e-10
  }, logical(1))]

  n_trt <- length(trt_idx)
  trt_pids <- all_pts$pid[trt_idx]

  sep <- paste(rep("-", 60), collapse = "")
  cat(sprintf("\n%s\n  SPEC: %s [%s] | POOL: %s\n  Covariates (%d): %s\n%s\n",
              sep, spec_name, spec_obj$label, toupper(pool_name),
              length(ps_vars), paste(ps_vars, collapse=", "), sep))

  # Fit PS
  ps_fml <- as.formula(paste("treated ~", paste(ps_vars, collapse = "+")))
  ps_fit <- suppressWarnings(glm(ps_fml, data = all_pts, family = binomial()))
  all_pts$ps <- predict(ps_fit, type = "response")

  # Match
  match_trt <- integer(n_trt); match_ctl <- integer(n_trt)
  matched <- logical(n_trt)

  for (k in seq_len(n_trt)) {
    rs <- risk_sets[[k]]
    if (length(rs) == 0) next
    ps_i <- all_pts$ps[trt_idx[k]]
    ps_dist <- abs(all_pts$ps[rs] - ps_i)
    within_cal <- which(ps_dist <= caliper)
    if (length(within_cal) == 0) next
    best <- within_cal[which.min(ps_dist[within_cal])]
    match_trt[k] <- trt_idx[k]; match_ctl[k] <- rs[best]
    matched[k] <- TRUE
  }
  n_matched <- sum(matched)
  cat(sprintf("  Matched: %d/%d (%.1f%%)\n", n_matched, n_trt, 100*n_matched/n_trt))
  if (n_matched < 50) { cat("  WARNING: <50 matches\n"); return(NULL) }

  # Temporal diagnostics
  ctl_indices <- match_ctl[matched]
  n_never <- sum(all_pts$treated[ctl_indices] == 0, na.rm=TRUE)
  n_future <- sum(all_pts$treated[ctl_indices] == 1, na.rm=TRUE)
  n_unique_ctl <- length(unique(ctl_indices))
  cat(sprintf("  Controls: %d never (%.0f%%), %d yet-untreated (%.0f%%) | Unique: %d (reuse %.1fx)\n",
              n_never, 100*n_never/n_matched, n_future, 100*n_future/n_matched,
              n_unique_ctl, n_matched/n_unique_ctl))

  # Balance + identify DR adjustment vars
  trt_m <- match_trt[matched]; ctl_m <- match_ctl[matched]
  smds <- sapply(ps_vars, function(v) {
    x1 <- all_pts[[v]][trt_m]; x0 <- all_pts[[v]][ctl_m]
    sp <- sqrt((var(x1,na.rm=T) + var(x0,na.rm=T)) / 2)
    if (is.na(sp) || sp < 1e-10) 0 else abs(mean(x1,na.rm=T) - mean(x0,na.rm=T)) / sp
  })
  max_smd <- max(smds, na.rm=TRUE)
  n_viol <- sum(smds > 0.1, na.rm=TRUE)
  adj_vars <- names(smds[smds > 0.1])

  cat(sprintf("  Balance: max SMD=%.3f, violations=%d/%d",
              max_smd, n_viol, length(ps_vars)))
  if (n_viol > 0) cat(sprintf(" (%s)", paste(adj_vars, collapse=", ")))
  cat("\n")

  # Outcomes at each horizon
  results <- list()

  for (target_h in TARGETS) {
    dcr_trt <- dcr_ctl <- numeric(n_matched)
    valid <- logical(n_matched); idx <- 0

    for (kk in which(matched)) {
      idx <- idx + 1
      tpid <- as.character(all_pts$pid[match_trt[kk]])
      cpid <- as.character(all_pts$pid[match_ctl[kk]])
      t_mg <- trt_tmg[kk]

      if (pool_name == "yet_untreated") {
        ctl_mg_h <- all_pts$mg_offset_h[match_ctl[kk]]
        if (!is.na(ctl_mg_h) && ctl_mg_h < t_mg + target_h) {
          valid[idx] <- FALSE; next
        }
      }

      pre_t  <- find_cr_pre(cr_list[[tpid]], t_mg)
      pre_c  <- find_cr_pre(cr_list[[cpid]], t_mg)
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
      results[[length(results)+1]] <- data.frame(
        spec=spec_name, pool=pool_name, target_h=target_h,
        method="psm", n=n_valid, did=NA, se=NA, p=NA,
        ci_lo=NA, ci_hi=NA, max_smd=max_smd, n_viol=n_viol)
      results[[length(results)+1]] <- data.frame(
        spec=spec_name, pool=pool_name, target_h=target_h,
        method="psm_dr", n=n_valid, did=NA, se=NA, p=NA,
        ci_lo=NA, ci_hi=NA, max_smd=max_smd, n_viol=n_viol)
      next
    }

    pair_df <- data.frame(
      delta_cr = c(dcr_trt[valid], dcr_ctl[valid]),
      treated = rep(c(1, 0), each = n_valid))

    trt_rows <- match_trt[matched][valid]
    ctl_rows <- match_ctl[matched][valid]
    for (av in adj_vars) {
      if (av %in% names(all_pts))
        pair_df[[av]] <- c(all_pts[[av]][trt_rows], all_pts[[av]][ctl_rows])
    }

    # PSM plain
    fit_p <- lm(delta_cr ~ treated, data = pair_df)
    ct_p <- safe_coeftest(fit_p)
    est_p <- ct_p["treated",1]; se_p <- ct_p["treated",2]; pv_p <- ct_p["treated",4]
    results[[length(results)+1]] <- data.frame(
      spec=spec_name, pool=pool_name, target_h=target_h, method="psm",
      n=n_valid, did=est_p, se=se_p, p=pv_p,
      ci_lo=est_p-1.96*se_p, ci_hi=est_p+1.96*se_p,
      max_smd=max_smd, n_viol=n_viol)

    # PSM + DR
    usable_adj <- intersect(adj_vars, names(pair_df))
    usable_adj <- usable_adj[vapply(usable_adj, function(v)
      var(pair_df[[v]], na.rm=TRUE) > 1e-10, logical(1))]
    if (length(usable_adj) > 0) {
      dr_fml <- as.formula(paste("delta_cr ~ treated +", paste(usable_adj, collapse="+")))
      fit_d <- tryCatch(lm(dr_fml, data = pair_df), error = function(e) fit_p)
    } else { fit_d <- fit_p }

    ct_d <- safe_coeftest(fit_d)
    est_d <- ct_d["treated",1]; se_d <- ct_d["treated",2]; pv_d <- ct_d["treated",4]
    results[[length(results)+1]] <- data.frame(
      spec=spec_name, pool=pool_name, target_h=target_h, method="psm_dr",
      n=n_valid, did=est_d, se=se_d, p=pv_d,
      ci_lo=est_d-1.96*se_d, ci_hi=est_d+1.96*se_d,
      max_smd=max_smd, n_viol=n_viol)
  }

  # Print results
  res <- do.call(rbind, results)
  for (mt in c("psm","psm_dr")) {
    cat(sprintf("\n  -- %s [%s | %s] --\n", toupper(mt), spec_name, pool_name))
    cat("  hour    DiD        SE       P        95% CI                n\n")
    sub <- res[res$method == mt, ]
    for (i in seq_len(nrow(sub))) {
      r <- sub[i,]
      if (is.na(r$did)) { cat(sprintf("  %3dh   (n=%d)\n", r$target_h, r$n)); next }
      sig <- if (!is.na(r$p) && r$p < 0.05) " *" else "  "
      pri <- if (r$target_h == PRIMARY_H) " << PRIMARY" else ""
      cat(sprintf("  %3dh   %+.4f   %.4f   %.4f%s [%+.4f,%+.4f]  %4d%s\n",
                  r$target_h, r$did, r$se, r$p, sig,
                  r$ci_lo, r$ci_hi, r$n, pri))
    }
  }

  # Balance table
  cat(sprintf("\n  -- BALANCE [%s | %s] --\n", spec_name, pool_name))
  trt_all <- which(all_pts$treated == 1); ctl_all <- which(all_pts$treated == 0)
  for (v in ps_vars) {
    x1r <- all_pts[[v]][trt_all]; x0r <- all_pts[[v]][ctl_all]
    spr <- sqrt((var(x1r,na.rm=T)+var(x0r,na.rm=T))/2)
    smd_r <- if (!is.na(spr) && spr>1e-10) abs(mean(x1r,na.rm=T)-mean(x0r,na.rm=T))/spr else NA
    smd_m <- smds[v]
    flag <- if (!is.na(smd_m) && smd_m > 0.1) "VIOL" else "ok"
    cat(sprintf("  %-28s  %.3f -> %.3f  %s\n", v,
                ifelse(is.na(smd_r), NA, smd_r), ifelse(is.na(smd_m), NA, smd_m), flag))
  }

  pairs_df <- data.frame(
    trt_pid = all_pts$pid[match_trt[matched]],
    ctl_pid = all_pts$pid[match_ctl[matched]],
    t_mg = trt_tmg[matched])

  list(results = res, pairs = pairs_df, n_matched = n_matched,
       max_smd = max_smd, n_viol = n_viol)
}


# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) { cat("Usage: Rscript 02_psm.R <db>\n"); quit(status = 1) }
db <- toupper(args[1]); tag <- tolower(db)

SEP <- paste(rep("=", 70), collapse = "")
cat(sprintf("\n%s\n02_psm.R — Canonical Risk-Set PSM: %s\n", SEP, db))
cat(sprintf("  Primary: 21 vars (base + all labs)  |  Sensitivity: 19 vars (no K/Mg)\n"))
cat(sprintf("  Pools: yet-untreated, never-treated  |  Methods: PSM, PSM+DR\n"))
cat(sprintf("  Horizons: %s  |  Labs: LAST  |  MICE m=%d\n%s\n",
            paste(TARGETS, "h", sep="", collapse=","), M_IMP, SEP))

# Load data
all_pts <- read.csv(file.path(RESULTS, sprintf("did_all_%s.csv", tag)),
                    stringsAsFactors = FALSE)
cr_all  <- read.csv(file.path(RESULTS, sprintf("did_cr_all_%s.csv", tag)),
                    stringsAsFactors = FALSE)

cr_id <- if ("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
cr_all$pid <- cr_all[[cr_id]]
if (!"offset_h" %in% names(cr_all)) cr_all$offset_h <- cr_all$labresultoffset / 60
cr_all <- cr_all[order(cr_all$pid, cr_all$offset_h), ]

cat("  Loading labs...\n")
labs_raw <- read.csv(file.path(RESULTS, sprintf("did_labs_all_%s.csv", tag)),
                     stringsAsFactors = FALSE)
pid_col_lab <- if ("patientunitstayid" %in% names(labs_raw)) "patientunitstayid" else "stay_id"
if (pid_col_lab %in% names(labs_raw)) labs_raw$pid <- labs_raw[[pid_col_lab]]
labs_elec <- labs_raw[labs_raw$lab_name %in% c("magnesium","potassium","calcium","lactate"), ]
labs_hr   <- labs_raw[labs_raw$lab_name == "heartrate", ]
if (nrow(labs_hr) > 500000) {
  labs_hr$hour_bin <- floor(labs_hr$offset_h)
  labs_hr <- labs_hr[order(labs_hr$pid, labs_hr$offset_h), ]
  labs_hr <- labs_hr[!duplicated(paste(labs_hr$pid, labs_hr$hour_bin)), ]
  labs_hr$hour_bin <- NULL
  cat(sprintf("    HR downsampled: %d -> %d\n",
              sum(labs_raw$lab_name == "heartrate"), nrow(labs_hr)))
}
labs <- rbind(labs_elec, labs_hr)
rm(labs_raw, labs_elec, labs_hr); gc()

N <- nrow(all_pts)
cat(sprintf("  Patients: %d (%d treated, %d control)\n\n", N,
            sum(all_pts$treated == 1), sum(all_pts$treated == 0)))

# LAST lab values → last_* columns
cat("  Computing LAST lab values (closest to T0)...\n")
for (ln in LAB_BASES) {
  sub <- labs[labs$lab_name == ln, ]
  if (nrow(sub) == 0) next
  sub$mg_offset_h <- all_pts$mg_offset_h[match(sub$pid, all_pts$pid)]
  sub <- sub[sub$offset_h >= 0 &
             (is.na(sub$mg_offset_h) | sub$offset_h < sub$mg_offset_h), ]
  if (nrow(sub) == 0) next
  s <- sub[order(-sub$offset_h), ]
  s <- s[!duplicated(s$pid), ]
  idx <- match(all_pts$pid, s$pid)
  col_name <- paste0("last_", ln)
  all_pts[[col_name]] <- s$value[idx]
  nf <- sum(!is.na(all_pts[[col_name]]))
  cat(sprintf("    %s: %d/%d (%.0f%%)\n", col_name, nf, N, 100*nf/N))
}
all_pts$last_lactate_missing <- as.integer(is.na(all_pts$last_lactate))

# Treatment indices
trt_idx <- which(all_pts$treated == 1 & !is.na(all_pts$mg_offset_h))
n_trt <- length(trt_idx)
trt_pids <- all_pts$pid[trt_idx]
trt_tmg  <- all_pts$mg_offset_h[trt_idx]
cat(sprintf("  Treated eligible: %d\n", n_trt))

# Cr lists + AKI
cat("  Pre-computing Cr lists and AKI times...\n")
cr_list <- split(cr_all[, c("labresult","offset_h")], cr_all$pid)

earliest_cr <- sapply(cr_list, function(x) min(x$offset_h, na.rm = TRUE))
all_pts$earliest_cr_h <- earliest_cr[as.character(all_pts$pid)]
all_pts$earliest_cr_h[is.na(all_pts$earliest_cr_h)] <- Inf

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

# Risk sets
cat("  Pre-computing risk sets...\n")
rs_yt <- rs_nt <- vector("list", n_trt)
for (k in seq_len(n_trt)) {
  t_mg <- trt_tmg[k]
  rs_yt[[k]] <- which(
    all_pts$icu_discharge_h > t_mg &
    (is.na(all_pts$first_aki_h) | all_pts$first_aki_h > t_mg) &
    all_pts$earliest_cr_h <= t_mg &
    (is.na(all_pts$mg_offset_h) | all_pts$mg_offset_h > t_mg + PRIMARY_H) &
    all_pts$pid != trt_pids[k])
  rs_nt[[k]] <- which(
    all_pts$treated == 0 &
    all_pts$icu_discharge_h > t_mg &
    (is.na(all_pts$first_aki_h) | all_pts$first_aki_h > t_mg) &
    all_pts$earliest_cr_h <= t_mg &
    all_pts$pid != trt_pids[k])
}
cat(sprintf("    YT: median=%.0f  NT: median=%.0f\n",
            median(sapply(rs_yt, length)), median(sapply(rs_nt, length))))

# MICE (once, on all 21 candidate vars)
all_candidate <- unique(c(PS_BASE, PS_LABS_ALL))
all_candidate <- intersect(all_candidate, names(all_pts))
all_candidate <- all_candidate[vapply(all_candidate, function(v) {
  x <- all_pts[[v]]; !all(is.na(x)) && var(x, na.rm=TRUE) > 1e-10
}, logical(1))]

to_impute <- all_candidate[vapply(all_candidate, function(v)
  any(is.na(all_pts[[v]])), logical(1))]
cat(sprintf("\n  MICE m=%d on %d vars: %s\n", M_IMP, length(to_impute),
            paste(to_impute, collapse=", ")))

if (length(to_impute) > 0) {
  mice_df <- all_pts[, c("treated", all_candidate)]
  meth <- rep("", ncol(mice_df)); names(meth) <- names(mice_df)
  for (v in to_impute) meth[v] <- "pmm"
  imp <- mice(mice_df, m = M_IMP, method = meth, printFlag = FALSE, maxit = 10)
  cat(sprintf("  MICE done. Logged events: %d\n", nrow(imp$loggedEvents)))
  for (v in to_impute) {
    vals <- sapply(1:M_IMP, function(m) complete(imp, m)[[v]])
    all_pts[[v]] <- rowMeans(vals, na.rm = TRUE)
  }
}

# Caliper from primary PS
ps_fml_full <- as.formula(paste("treated ~", paste(
  intersect(SPECS$primary$vars, names(all_pts)), collapse = "+")))
ps_fit_full <- suppressWarnings(glm(ps_fml_full, data = all_pts, family = binomial()))
caliper <- CALIPER_SD * sd(predict(ps_fit_full, type = "response"), na.rm = TRUE)
cat(sprintf("  Caliper (from primary PS): %.4f\n", caliper))

# ═══════════════════════════════════════════════════════════════════
# RUN 2 SPECS x 2 POOLS = 4 COMBINATIONS
# ═══════════════════════════════════════════════════════════════════

all_results <- list()
all_pairs   <- list()

for (spec_name in names(SPECS)) {
  for (pool_name in c("yet_untreated", "never_treated")) {
    rs <- if (pool_name == "yet_untreated") rs_yt else rs_nt
    out <- run_spec_pool(spec_name, SPECS[[spec_name]], pool_name,
                         all_pts, trt_idx, rs, cr_list, trt_tmg, caliper)
    if (!is.null(out)) {
      all_results[[length(all_results)+1]] <- out$results
      tag_pair <- sprintf("%s_%s", spec_name, pool_name)
      all_pairs[[tag_pair]] <- out$pairs
    }
  }
}

# ═══════════════════════════════════════════════════════════════════
# COMPARISON TABLE
# ═══════════════════════════════════════════════════════════════════

res_all <- do.call(rbind, all_results)
res_all$db <- db

cat(sprintf("\n%s\n  HEAD-TO-HEAD: %s (PSM_DR, primary=24h)\n%s\n", SEP, db, SEP))
cat("                yet-untreated              never-treated\n")
cat("  spec          DiD       P      n         DiD       P      n\n")
cat("  ────────────  ────────  ─────  ────      ────────  ─────  ────\n")
for (sn in names(SPECS)) {
  parts <- c()
  for (pl in c("yet_untreated","never_treated")) {
    r <- res_all[res_all$spec==sn & res_all$pool==pl &
                 res_all$method=="psm_dr" & res_all$target_h==PRIMARY_H, ]
    if (nrow(r)==1 && !is.na(r$did)) {
      parts <- c(parts, sprintf("%+.4f  %.3f  %4d", r$did, r$p, r$n))
    } else {
      parts <- c(parts, "   --      --    -- ")
    }
  }
  cat(sprintf("  %-14s%s      %s\n", sn, parts[1], parts[2]))
}

# ═══════════════════════════════════════════════════════════════════
# SAVE
# ═══════════════════════════════════════════════════════════════════

write.csv(res_all,
          file.path(RESULTS, sprintf("did_riskset_%s.csv", tag)),
          row.names = FALSE)

for (nm in names(all_pairs)) {
  write.csv(all_pairs[[nm]],
            file.path(RESULTS, sprintf("did_pairs_%s_%s.csv", nm, tag)),
            row.names = FALSE)
}

cat(sprintf("\n%s\n02_psm.R -- %s COMPLETE\n", SEP, db))
cat(sprintf("  Output: did_riskset_%s.csv (%d rows)\n", tag, nrow(res_all)))
cat(sprintf("  Pairs:  did_pairs_{spec}_{pool}_%s.csv\n", tag))
cat(sprintf("  Next:   python 04_fig_timecourse.py\n%s\n", SEP))
