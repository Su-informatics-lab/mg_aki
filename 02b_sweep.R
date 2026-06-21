#!/usr/bin/env Rscript
# ============================================================================
# 02b_sweep.R — Specification curve for v5 risk-set design
#
# Grid:
#   Covariates: 15 base + 6 toggle groups (D1,L1-L5) = 64 specs
#   Methods:    PSM | PSM+DR | AIPW (risk-set, binned sequential trial)
#   Pools:      yet-untreated | never-treated (all 3 methods)
#   Time:       24h (primary), 36h (persistence)
#
# PSM:    risk-set 1:1, caliper 0.2 SD, HC1 SE
# PSM+DR: same + outcome-adjust covariates with matched SMD > 0.1
# AIPW:   risk-set binned sequential trial (dynamic T₀, no immortal time bias)
#
# Usage: Rscript 02b_sweep.R eicu
#        Rscript 02b_sweep.R mimic
# ============================================================================

suppressPackageStartupMessages({ library(sandwich); library(lmtest); library(mice) })

RESULTS    <- path.expand("~/mg_aki/results")
PRIMARY_H  <- 24
CR_WINDOW  <- 12
CALIPER_SD <- 0.2
M_IMP      <- 20
TARGETS    <- c(24, 36)

# ── Base covariates (always included) ─────────────────────────────
BASE <- c("age","is_female","bmi",
          "surg_cabg","surg_valve","surg_combined",
          "heart_failure","hypertension","diabetes","ckd",
          "copd","pvd","stroke","liver_disease","egfr")

# ── Toggle groups ─────────────────────────────────────────────────
TOGGLES <- list(
  D1_drugs   = c("ppi_chronic","loop_diuretic_chronic",
                  "acei_arb_chronic","nsaid_chronic"),
  L1_potass  = c("first_potassium"),
  L2_calcium = c("first_calcium"),
  L3_lactate = c("first_lactate","first_lactate_missing"),
  L4_mg      = c("first_magnesium"),
  L5_hr      = c("first_heartrate")
)

PS_LAB_BASES <- c("magnesium","potassium","calcium","lactate","heartrate")

# ═══════════════════════════════════════════════════════════════════
# HELPERS (same as 02_psm.R)
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
if (length(args) < 1) { cat("Usage: Rscript 02b_sweep.R <db>\n"); quit(status = 1) }
db <- toupper(args[1]); tag <- tolower(db)

SEP <- paste(rep("=", 70), collapse = "")
n_toggles <- length(TOGGLES)
n_specs <- 2^n_toggles
cat(sprintf("\n%s\n02b_sweep.R — Specification Curve: %s\n", SEP, db))
cat(sprintf("  %d base + %d toggle groups = %d specs\n", length(BASE), n_toggles, n_specs))
cat(sprintf("  Methods: PSM, PSM+DR, AIPW | Pools: yet-untreated, never-treated\n"))
cat(sprintf("  Time: %s\n%s\n", paste(TARGETS, "h", sep="", collapse=", "), SEP))

# ── Load data ─────────────────────────────────────────────────────
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
}
labs <- rbind(labs_elec, labs_hr)
rm(labs_raw, labs_elec, labs_hr); gc()

N <- nrow(all_pts)
cat(sprintf("  Patients: %d (%d treated, %d control)\n", N,
            sum(all_pts$treated == 1), sum(all_pts$treated == 0)))

# ── LAST lab values ───────────────────────────────────────────────
cat("  Computing LAST lab values...\n")
for (ln in PS_LAB_BASES) {
  sub <- labs[labs$lab_name == ln, ]
  if (nrow(sub) == 0) next
  sub$mg_offset_h <- all_pts$mg_offset_h[match(sub$pid, all_pts$pid)]
  sub <- sub[sub$offset_h >= 0 &
             (is.na(sub$mg_offset_h) | sub$offset_h < sub$mg_offset_h), ]
  if (nrow(sub) == 0) next
  s <- sub[order(-sub$offset_h), ]
  s <- s[!duplicated(s$pid), ]
  idx <- match(all_pts$pid, s$pid)
  all_pts[[paste0("first_", ln)]] <- s$value[idx]
}
all_pts$first_lactate_missing <- as.integer(is.na(all_pts$first_lactate))

# ══════════════════════════════════════════════════════════════════
# MISSINGNESS TABLE
# ══════════════════════════════════════════════════════════════════
cat(sprintf("\n  ── MISSINGNESS TABLE ──\n"))
cat(sprintf("  %-28s  %8s  %8s  %8s\n", "Variable", "Overall", "Treated", "Control"))
cat(sprintf("  %-28s  %8s  %8s  %8s\n", "────────────────────────────",
            "────────", "────────", "────────"))

all_candidate <- unique(c(BASE, unlist(TOGGLES)))
for (v in all_candidate) {
  if (!v %in% names(all_pts)) {
    cat(sprintf("  %-28s  %8s\n", v, "NOT FOUND"))
    next
  }
  n_miss <- sum(is.na(all_pts[[v]]))
  n_trt_miss <- sum(is.na(all_pts[[v]][all_pts$treated == 1]))
  n_ctl_miss <- sum(is.na(all_pts[[v]][all_pts$treated == 0]))
  n_trt <- sum(all_pts$treated == 1)
  n_ctl <- sum(all_pts$treated == 0)
  cat(sprintf("  %-28s  %5d (%2.0f%%)  %5d (%2.0f%%)  %5d (%2.0f%%)\n",
              v, n_miss, 100*n_miss/N,
              n_trt_miss, 100*n_trt_miss/n_trt,
              n_ctl_miss, 100*n_ctl_miss/n_ctl))
}

# ── Treatment indices + Cr lists ──────────────────────────────────
trt_idx <- which(all_pts$treated == 1 & !is.na(all_pts$mg_offset_h))
n_trt <- length(trt_idx)
trt_pids <- all_pts$pid[trt_idx]
trt_tmg  <- all_pts$mg_offset_h[trt_idx]

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

# ── Precompute risk sets (once per pool) ──────────────────────────
cat("\n  Pre-computing risk sets...\n")
rs_yt <- rs_nt <- vector("list", n_trt)

for (k in seq_len(n_trt)) {
  t_mg <- trt_tmg[k]

  # Yet-untreated: Mg-free through t_mg + PRIMARY_H
  rs_yt[[k]] <- which(
    all_pts$icu_discharge_h > t_mg &
    (is.na(all_pts$first_aki_h) | all_pts$first_aki_h > t_mg) &
    all_pts$earliest_cr_h <= t_mg &
    (is.na(all_pts$mg_offset_h) | all_pts$mg_offset_h > t_mg + PRIMARY_H) &
    all_pts$pid != trt_pids[k])

  # Never-treated: treated == 0
  rs_nt[[k]] <- which(
    all_pts$treated == 0 &
    all_pts$icu_discharge_h > t_mg &
    (is.na(all_pts$first_aki_h) | all_pts$first_aki_h > t_mg) &
    all_pts$earliest_cr_h <= t_mg &
    all_pts$pid != trt_pids[k])
}
cat(sprintf("    YT: median=%.0f  NT: median=%.0f\n",
            median(sapply(rs_yt, length)),
            median(sapply(rs_nt, length))))

# ── MICE (once, on ALL candidate vars) ────────────────────────────
all_ps_vars <- unique(c(BASE, unlist(TOGGLES)))
all_ps_vars <- intersect(all_ps_vars, names(all_pts))
# Drop zero-variance and all-NA
all_ps_vars <- all_ps_vars[vapply(all_ps_vars, function(v) {
  x <- all_pts[[v]]; !all(is.na(x)) && var(x, na.rm = TRUE) > 1e-10
}, logical(1))]

to_impute <- all_ps_vars[vapply(all_ps_vars, function(v)
  any(is.na(all_pts[[v]])), logical(1))]
cat(sprintf("\n  MICE m=%d on %d vars: %s\n", M_IMP, length(to_impute),
            paste(to_impute, collapse = ", ")))

if (length(to_impute) > 0) {
  mice_df <- all_pts[, c("treated", all_ps_vars)]
  meth <- rep("", ncol(mice_df)); names(meth) <- names(mice_df)
  for (v in to_impute) meth[v] <- "pmm"
  imp <- mice(mice_df, m = M_IMP, method = meth, printFlag = FALSE, maxit = 10)
  cat(sprintf("  MICE done. Logged events: %d\n", nrow(imp$loggedEvents)))
  for (v in to_impute) {
    vals <- sapply(1:M_IMP, function(m) complete(imp, m)[[v]])
    all_pts[[v]] <- rowMeans(vals, na.rm = TRUE)
  }
}

# ── Risk-set AIPW (sequential trial: dynamic T₀, binned) ─────────
#
# For each 2h bin of treated T₀:
#   1. Treated in bin: use their actual T₀ for ΔCr
#   2. Controls: everyone eligible at bin midpoint, ΔCr anchored there
#   3. Run AIPW within bin (bin-specific PS + outcome models)
# Pool across bins: weighted by n_treated per bin
#
run_aipw_riskset <- function(all_pts, cr_list, ps_vars, target_h,
                              pool = "never_treated", bin_w = 2) {

  trt <- all_pts[all_pts$treated == 1 & !is.na(all_pts$mg_offset_h), ]
  bins <- seq(0, max(trt$mg_offset_h, na.rm = TRUE) + bin_w, by = bin_w)

  bin_res <- list()

  for (b in seq_len(length(bins) - 1)) {
    t_lo <- bins[b]; t_hi <- bins[b + 1]; t_mid <- (t_lo + t_hi) / 2

    # Treated in this bin
    trt_bin <- trt[trt$mg_offset_h >= t_lo & trt$mg_offset_h < t_hi, ]
    if (nrow(trt_bin) < 5) next

    # Controls eligible at t_lo (alive, in ICU, no AKI, have Cr)
    if (pool == "never_treated") {
      ctl_idx <- which(
        all_pts$treated == 0 &
        all_pts$icu_discharge_h > t_lo &
        (is.na(all_pts$first_aki_h) | all_pts$first_aki_h > t_lo) &
        all_pts$earliest_cr_h <= t_lo)
    } else {
      ctl_idx <- which(
        all_pts$icu_discharge_h > t_lo &
        (is.na(all_pts$first_aki_h) | all_pts$first_aki_h > t_lo) &
        all_pts$earliest_cr_h <= t_lo &
        (is.na(all_pts$mg_offset_h) | all_pts$mg_offset_h > t_lo + target_h) &
        !all_pts$pid %in% trt_bin$pid)
    }
    if (length(ctl_idx) < 20) next

    # ΔCr for treated: each uses own T₀
    trt_dcr <- vapply(seq_len(nrow(trt_bin)), function(j) {
      pid <- as.character(trt_bin$pid[j])
      t0 <- trt_bin$mg_offset_h[j]
      pre <- find_cr_pre(cr_list[[pid]], t0)
      post <- find_cr(cr_list[[pid]], t0 + target_h)
      if (any(is.na(c(pre[1], post[1])))) NA_real_ else post[1] - pre[1]
    }, numeric(1))

    # ΔCr for controls: use bin midpoint as shared T₀
    ctl_dcr <- vapply(ctl_idx, function(idx) {
      pid <- as.character(all_pts$pid[idx])
      pre <- find_cr_pre(cr_list[[pid]], t_mid)
      post <- find_cr(cr_list[[pid]], t_mid + target_h)
      if (any(is.na(c(pre[1], post[1])))) NA_real_ else post[1] - pre[1]
    }, numeric(1))

    # Combine valid observations
    trt_ok <- which(!is.na(trt_dcr)); ctl_ok <- which(!is.na(ctl_dcr))
    if (length(trt_ok) < 5 || length(ctl_ok) < 20) next

    trt_rows <- match(trt_bin$pid[trt_ok], all_pts$pid)
    d <- rbind(
      cbind(all_pts[trt_rows, , drop = FALSE], delta_cr = trt_dcr[trt_ok]),
      cbind(all_pts[ctl_idx[ctl_ok], , drop = FALSE], delta_cr = ctl_dcr[ctl_ok]))
    d$treated <- c(rep(1L, length(trt_ok)), rep(0L, length(ctl_ok)))

    avail <- intersect(ps_vars, names(d))
    avail <- avail[vapply(avail, function(v)
      var(d[[v]], na.rm = TRUE) > 1e-10, logical(1))]
    d <- d[complete.cases(d[, c(avail, "delta_cr", "treated")]), ]
    nt <- sum(d$treated == 1); nc <- sum(d$treated == 0)
    if (nt < 5 || nc < 20) next

    # Bin-specific PS
    ps_fml <- as.formula(paste("treated ~", paste(avail, collapse = "+")))
    ps_fit <- tryCatch(suppressWarnings(
      glm(ps_fml, data = d, family = binomial())), error = function(e) NULL)
    if (is.null(ps_fit)) next
    d$ps <- pmax(pmin(fitted(ps_fit), 0.99), 0.01)

    # Bin-specific outcome models
    out_fml <- as.formula(paste("delta_cr ~", paste(avail, collapse = "+")))
    m1 <- tryCatch(lm(out_fml, data = d[d$treated == 1, ]), error = function(e) NULL)
    m0 <- tryCatch(lm(out_fml, data = d[d$treated == 0, ]), error = function(e) NULL)
    if (is.null(m1) || is.null(m0)) next

    mu1 <- predict(m1, newdata = d); mu0 <- predict(m0, newdata = d)
    Y <- d$delta_cr; Tr <- d$treated; e <- d$ps

    phi <- Tr * (Y - mu1) / e - (1 - Tr) * (Y - mu0) / (1 - e) + (mu1 - mu0)
    bin_res[[length(bin_res) + 1]] <- list(tau = mean(phi), n_trt = nt, phi = phi)
  }

  if (length(bin_res) == 0) return(NULL)

  # Pool: n_trt-weighted average of bin-specific τ̂
  wts  <- sapply(bin_res, `[[`, "n_trt")
  taus <- sapply(bin_res, `[[`, "tau")
  tau_pool <- sum(wts * taus) / sum(wts)

  # SE from pooled influence function
  all_phi <- unlist(lapply(bin_res, `[[`, "phi"))
  se_pool <- sd(all_phi) / sqrt(length(all_phi))

  data.frame(did = tau_pool, se = se_pool,
             p = 2 * pnorm(-abs(tau_pool / se_pool)),
             n_trt = sum(wts), n_bins = length(bin_res))
}

# ═══════════════════════════════════════════════════════════════════
# SWEEP
# ═══════════════════════════════════════════════════════════════════
toggle_names <- names(TOGGLES)
total <- n_specs * length(TARGETS)
cat(sprintf("\n  Starting sweep: %d specs × %d time points...\n\n", n_specs, length(TARGETS)))

all_rows <- list(); ridx <- 0
t_start <- Sys.time()

for (i in 0:(n_specs - 1)) {
  bits <- as.integer(intToBits(i)[1:n_toggles])
  on_names <- toggle_names[bits == 1]

  ps_vars <- BASE
  for (tn in on_names) ps_vars <- c(ps_vars, TOGGLES[[tn]])
  ps_vars <- unique(ps_vars)
  # Keep only available non-constant vars
  ps_vars <- intersect(ps_vars, names(all_pts))
  ps_vars <- ps_vars[vapply(ps_vars, function(v)
    var(all_pts[[v]], na.rm = TRUE) > 1e-10, logical(1))]

  label <- if (length(on_names) == 0) "base_only" else paste(on_names, collapse = "+")
  spec_id <- i + 1

  # ── Fit PS (once per spec) ──────────────────────────────────
  ps_fml <- as.formula(paste("treated ~", paste(ps_vars, collapse = "+")))
  ps_fit <- tryCatch(glm(ps_fml, data = all_pts, family = binomial()),
                     error = function(e) NULL)
  if (is.null(ps_fit)) {
    cat(sprintf("  [%3d] PS fit failed: %s\n", spec_id, label))
    next
  }
  all_pts$ps <- predict(ps_fit, type = "response")
  caliper <- CALIPER_SD * sd(all_pts$ps, na.rm = TRUE)

  # ── Match per pool (once per spec) ─────────────────────────
  for (pool in c("yet_untreated", "never_treated")) {
    rs <- if (pool == "yet_untreated") rs_yt else rs_nt

    match_trt <- integer(n_trt); match_ctl <- integer(n_trt)
    matched <- logical(n_trt)

    for (k in seq_len(n_trt)) {
      r <- rs[[k]]
      if (length(r) == 0) next
      ps_i <- all_pts$ps[trt_idx[k]]
      ps_dist <- abs(all_pts$ps[r] - ps_i)
      within_cal <- which(ps_dist <= caliper)
      if (length(within_cal) == 0) next
      best <- within_cal[which.min(ps_dist[within_cal])]
      match_trt[k] <- trt_idx[k]; match_ctl[k] <- r[best]
      matched[k] <- TRUE
    }
    n_matched <- sum(matched)
    if (n_matched < 50) next

    # Balance: max SMD + which vars violate
    trt_m <- match_trt[matched]; ctl_m <- match_ctl[matched]
    smds <- sapply(ps_vars, function(v) {
      x1 <- all_pts[[v]][trt_m]; x0 <- all_pts[[v]][ctl_m]
      sp <- sqrt((var(x1, na.rm=T) + var(x0, na.rm=T)) / 2)
      if (is.na(sp) || sp < 1e-10) 0 else abs(mean(x1,na.rm=T) - mean(x0,na.rm=T)) / sp
    })
    max_smd <- max(smds, na.rm = TRUE)
    n_viol <- sum(smds > 0.1, na.rm = TRUE)
    adj_vars <- names(smds[smds > 0.1])

    # ── Outcomes per time point ───────────────────────────────
    for (target_h in TARGETS) {
      dcr_trt <- dcr_ctl <- numeric(n_matched)
      valid <- logical(n_matched); idx <- 0

      for (kk in which(matched)) {
        idx <- idx + 1
        tpid <- as.character(all_pts$pid[match_trt[kk]])
        cpid <- as.character(all_pts$pid[match_ctl[kk]])
        t_mg <- trt_tmg[kk]

        if (pool == "yet_untreated") {
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
        ridx <- ridx + 1
        all_rows[[ridx]] <- data.frame(
          spec_id=spec_id, target_h=target_h, pool=pool,
          method="psm", n_covars=length(ps_vars), label=label,
          n=n_valid, did=NA, se=NA, p=NA, max_smd=max_smd, n_viol=n_viol,
          stringsAsFactors=FALSE)
        ridx <- ridx + 1
        all_rows[[ridx]] <- data.frame(
          spec_id=spec_id, target_h=target_h, pool=pool,
          method="psm_dr", n_covars=length(ps_vars), label=label,
          n=n_valid, did=NA, se=NA, p=NA, max_smd=max_smd, n_viol=n_viol,
          stringsAsFactors=FALSE)
        next
      }

      pair_df <- data.frame(
        delta_cr = c(dcr_trt[valid], dcr_ctl[valid]),
        treated = rep(c(1, 0), each = n_valid))

      # Attach covariates for DR adjustment
      trt_rows <- match_trt[matched][valid]
      ctl_rows <- match_ctl[matched][valid]
      for (av in adj_vars) {
        if (av %in% names(all_pts)) {
          pair_df[[av]] <- c(all_pts[[av]][trt_rows], all_pts[[av]][ctl_rows])
        }
      }

      # --- PSM plain ---
      fit_plain <- lm(delta_cr ~ treated, data = pair_df)
      ct_plain <- safe_coeftest(fit_plain)
      est_p <- if (!is.null(ct_plain)) ct_plain["treated", 1] else NA
      se_p  <- if (!is.null(ct_plain)) ct_plain["treated", 2] else NA
      pv_p  <- if (!is.null(ct_plain)) ct_plain["treated", 4] else NA

      ridx <- ridx + 1
      all_rows[[ridx]] <- data.frame(
        spec_id=spec_id, target_h=target_h, pool=pool,
        method="psm", n_covars=length(ps_vars), label=label,
        n=n_valid, did=est_p, se=se_p, p=pv_p,
        max_smd=max_smd, n_viol=n_viol, stringsAsFactors=FALSE)

      # --- PSM + DR (adjust imbalanced covars) ---
      usable_adj <- intersect(adj_vars, names(pair_df))
      usable_adj <- usable_adj[vapply(usable_adj, function(v)
        var(pair_df[[v]], na.rm = TRUE) > 1e-10, logical(1))]

      if (length(usable_adj) > 0) {
        dr_fml <- as.formula(paste("delta_cr ~ treated +",
                                    paste(usable_adj, collapse = "+")))
        fit_dr <- tryCatch(lm(dr_fml, data = pair_df), error = function(e) NULL)
      } else {
        fit_dr <- fit_plain  # no adjustment needed
      }

      ct_dr <- if (!is.null(fit_dr)) safe_coeftest(fit_dr) else NULL
      est_d <- if (!is.null(ct_dr) && "treated" %in% rownames(ct_dr)) ct_dr["treated",1] else est_p
      se_d  <- if (!is.null(ct_dr) && "treated" %in% rownames(ct_dr)) ct_dr["treated",2] else se_p
      pv_d  <- if (!is.null(ct_dr) && "treated" %in% rownames(ct_dr)) ct_dr["treated",4] else pv_p

      ridx <- ridx + 1
      all_rows[[ridx]] <- data.frame(
        spec_id=spec_id, target_h=target_h, pool=pool,
        method="psm_dr", n_covars=length(ps_vars), label=label,
        n=n_valid, did=est_d, se=se_d, p=pv_d,
        max_smd=max_smd, n_viol=n_viol, stringsAsFactors=FALSE)

    }  # target_h
  }  # pool

  # ── AIPW (risk-set, both pools) ───────────────────────────────
  for (pool in c("yet_untreated", "never_treated")) {
    for (target_h in TARGETS) {
      aipw <- tryCatch(
        run_aipw_riskset(all_pts, cr_list, ps_vars, target_h, pool = pool),
        error = function(e) NULL)
      ridx <- ridx + 1
      all_rows[[ridx]] <- data.frame(
        spec_id=spec_id, target_h=target_h, pool=pool,
        method="aipw", n_covars=length(ps_vars), label=label,
        n=if (!is.null(aipw)) aipw$n_trt else NA,
        did=if (!is.null(aipw)) aipw$did else NA,
        se=if (!is.null(aipw)) aipw$se else NA,
        p=if (!is.null(aipw)) aipw$p else NA,
        max_smd=NA, n_viol=NA, stringsAsFactors=FALSE)
    }
  }

  # Progress
  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  eta <- elapsed / spec_id * (n_specs - spec_id)
  cat(sprintf("  [%3d/%d] %2d covars  yt:%s  nt:%s  %.1fmin  ETA %.0fmin  %s\n",
              spec_id, n_specs, length(ps_vars),
              if (n_matched > 0) sprintf("%d", n_matched) else "—",
              if (n_matched > 0) sprintf("%d", n_matched) else "—",
              elapsed, eta, label))
}

# ═══════════════════════════════════════════════════════════════════
# RESULTS
# ═══════════════════════════════════════════════════════════════════
sweep <- do.call(rbind, all_rows)
write.csv(sweep, file.path(RESULTS, sprintf("did_sweep_%s.csv", tag)),
          row.names = FALSE)

cat(sprintf("\n%s\nSWEEP COMPLETE: %d rows\n%s\n", SEP, nrow(sweep), SEP))

# ── Summary per method × pool × time ─────────────────────────────
for (th in TARGETS) {
  cat(sprintf("\n  === %dh ===\n", th))
  for (mt in c("psm", "psm_dr", "aipw")) {
    for (pl in c("yet_untreated", "never_treated")) {
      s <- sweep[sweep$target_h == th & sweep$method == mt & sweep$pool == pl, ]
      s <- s[!is.na(s$did), ]
      if (nrow(s) == 0) next
      n_neg <- sum(s$did < 0)
      n_sig <- sum(s$p < 0.05, na.rm = TRUE)
      med_did <- median(s$did)
      range_did <- range(s$did)
      cat(sprintf("    %-8s %-16s: %d specs, median DiD=%+.4f [%+.4f,%+.4f], neg=%d/%d (%.0f%%), sig=%d\n",
                  mt, pl, nrow(s), med_did, range_did[1], range_did[2],
                  n_neg, nrow(s), 100*n_neg/nrow(s), n_sig))
    }
  }
}

# ── Toggle impact (24h, PSM, yet-untreated) ───────────────────────
cat(sprintf("\n  --- TOGGLE IMPACT (24h, PSM, yet-untreated) ---\n"))
s24 <- sweep[sweep$target_h == 24 & sweep$method == "psm" &
             sweep$pool == "yet_untreated" & !is.na(sweep$did), ]
if (nrow(s24) > 0) {
  cat(sprintf("  %-14s  ON_median  OFF_median  shift    direction\n", "Toggle"))
  cat(sprintf("  %-14s  ─────────  ──────────  ─────    ─────────\n", "──────────────"))
  for (k in seq_along(toggle_names)) {
    # Which specs have this toggle ON vs OFF?
    on_specs <- which(as.integer(intToBits(0:(n_specs-1))[k + (0:(n_specs-1))*32 +1]) == 1)
    # Simpler: reconstruct bits
    on_ids <- c(); off_ids <- c()
    for (ii in 0:(n_specs-1)) {
      b <- as.integer(intToBits(ii)[1:n_toggles])
      if (b[k] == 1) on_ids <- c(on_ids, ii+1) else off_ids <- c(off_ids, ii+1)
    }
    med_on  <- median(s24$did[s24$spec_id %in% on_ids], na.rm = TRUE)
    med_off <- median(s24$did[s24$spec_id %in% off_ids], na.rm = TRUE)
    shift <- med_on - med_off
    dir <- if (abs(shift) < 0.001) "stable" else if (shift < 0) "more protective" else "less protective"
    cat(sprintf("  %-14s  %+.4f    %+.4f     %+.4f   %s\n",
                toggle_names[k], med_on, med_off, shift, dir))
  }
}

cat(sprintf("\n  Saved: did_sweep_%s.csv (%d rows)\n", tag, nrow(sweep)))
cat(sprintf("  Runtime: %.1f minutes\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat(sprintf("%s\n", SEP))
