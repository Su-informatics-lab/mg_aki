#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────
# Mg → Cardiac Surgery Outcomes (eICU) — v2
# 03_models.R — Self-validating multi-outcome TTE
#
#   A. Prognostic: Mg level → AKI
#   B. Positive control: Mg supplementation → POAF (+ BB stratification)
#   C. Primary: Mg supplementation → AKI (severity-stratified + BB)
#   D. Negative controls (bias calibration)
#   E. Neuro outcomes (exploratory)
#   F. Secondary (RRT, mortality)
# ─────────────────────────────────────────────────────────────────────

local({
  pkgs <- c("tidyverse", "survival", "survey", "sandwich", "lmtest", "broom")
  miss <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(miss) > 0)
    install.packages(miss, repos = "https://cloud.r-project.org",
                     quiet = TRUE, Ncpus = parallel::detectCores())
  tryCatch(install.packages("EValue", repos = "https://cloud.r-project.org",
                            quiet = TRUE), error = function(e) NULL)
})
HAS_EVALUE <- requireNamespace("EValue", quietly = TRUE)

suppressPackageStartupMessages({
  library(tidyverse); library(survival); library(survey)
  library(sandwich); library(lmtest); library(broom); library(splines)
})

RESULTS <- path.expand("~/mg_aki/results")

# ─── Load ───────────────────────────────────────────────────────────
cat("Loading data...\n")
load(file.path(RESULTS, "02b_analysis_a_prepared.RData"))
dat_tte_a <- tryCatch(read_csv(file.path(RESULTS, "02d_tte_iptw.csv"),
                               show_col_types = FALSE), error = function(e) NULL)
dat_tte_b <- tryCatch(read_csv(file.path(RESULTS, "02e_tteb_iptw.csv"),
                               show_col_types = FALSE), error = function(e) NULL)
dat_m_b   <- tryCatch(read_csv(file.path(RESULTS, "02f_tteb_matched.csv"),
                               show_col_types = FALSE), error = function(e) NULL)

add_aki <- function(d) {
  if (!"aki_kdigo1" %in% names(d))
    d$aki_kdigo1 <- as.integer(d$aki_primary == 1 | d$aki_delta03 == 1)
  d
}
dat_a <- add_aki(dat_a)
if (!is.null(dat_tte_a)) dat_tte_a <- add_aki(dat_tte_a)
if (!is.null(dat_tte_b)) dat_tte_b <- add_aki(dat_tte_b)
if (!is.null(dat_m_b))   dat_m_b   <- add_aki(dat_m_b)

cov_rhs <- paste(c(
  "age_num", "is_female", "bmi", "surgery_type",
  "hx_chf", "hx_hypertension", "hx_diabetes", "hx_ckd",
  "hx_copd", "hx_pvd", "hx_stroke",
  "baseline_cr", "baseline_egfr",
  "nephrotox_loop_diuretic", "nephrotox_nsaid",
  "nephrotox_acei_arb", "nephrotox_ppi"
), collapse = " + ")
for (v in c("preop_antiarrhythmic", "first_k_value",
            "has_vasopressor", "first_hr", "first_ca_value")) {
  if (v %in% names(dat_a)) cov_rhs <- paste(cov_rhs, "+", v)
}

results_list <- list()

# ─── Helper: IPTW OR for one outcome ───────────────────────────────
iptw_or <- function(dat, outcome, label = outcome, wt = "iptw") {
  if (!outcome %in% names(dat)) return(NULL)
  d <- dat %>% filter(!is.na(.data[[outcome]]))
  nev <- sum(d[[outcome]], na.rm = TRUE)
  if (nev < 5) return(NULL)
  tryCatch({
    des <- if ("hospitalid" %in% names(d))
      svydesign(ids = ~hospitalid, weights = as.formula(paste0("~", wt)), data = d)
    else svydesign(ids = ~1, weights = as.formula(paste0("~", wt)), data = d)
    m <- svyglm(as.formula(paste(outcome, "~ trt")), design = des,
                family = quasibinomial())
    s <- tidy(m, conf.int = TRUE, exponentiate = TRUE) %>% filter(term == "trt")
    if (nrow(s) > 0) { s$n_events <- nev; s$label <- label; s$n_total <- nrow(d) }
    s
  }, error = function(e) NULL)
}

ptbl <- function(rows, title) {
  cat(sprintf("\n  %s\n", title))
  cat(sprintf("  %-44s %6s %16s %8s\n", "Outcome", "Events", "OR (95% CI)", "p"))
  cat("  ", strrep("-", 78), "\n")
  for (r in rows) if (!is.null(r) && nrow(r) > 0)
    cat(sprintf("  %-44s %6d %5.2f (%.2f-%.2f) %8.4f\n",
                r$label, r$n_events, r$estimate, r$conf.low, r$conf.high, r$p.value))
}

# =====================================================================
# A. PROGNOSTIC
# =====================================================================
cat("\n", strrep("=", 70), "\n")
cat("A. PROGNOSTIC: Mg level -> AKI\n")
cat(strrep("=", 70), "\n")

dat_a$mg_neg <- -dat_a$first_mg_value
aki_defs <- c(aki_kdigo1 = "KDIGO Stage >=1", aki_primary = "Ratio >=1.5x",
              aki_delta03 = "Delta >=0.3", aki_stage2 = "Stage >=2",
              aki_stage3 = "Stage >=3")

cat("\n  A1: OR per 1 mg/dL Mg increase:\n")
cat(sprintf("  %-44s %6s %16s %8s\n", "Outcome", "Events", "OR (95% CI)", "p"))
cat("  ", strrep("-", 78), "\n")
for (ov in names(aki_defs)) {
  if (!ov %in% names(dat_a)) next
  nev <- sum(dat_a[[ov]], na.rm = TRUE)
  tryCatch({
    m <- glm(as.formula(paste(ov, "~ mg_neg +", cov_rhs)),
             data = dat_a, family = binomial())
    s <- tidy(m, conf.int = TRUE, exponentiate = TRUE) %>% filter(term == "mg_neg")
    if (nrow(s) > 0) {
      or <- 1/s$estimate; lo <- 1/s$conf.high; hi <- 1/s$conf.low
      cat(sprintf("  %-44s %6d %5.2f (%.2f-%.2f) %8.4f\n",
                  aki_defs[[ov]], nev, or, lo, hi, s$p.value))
      results_list[[paste0("A1_", ov)]] <<- tibble(
        model = paste0("A1_", ov), term = "mg_per_1_increase",
        estimate = or, conf.low = lo, conf.high = hi,
        p.value = s$p.value, n_events = nev, label = aki_defs[[ov]])
    }
  }, error = function(e) cat(sprintf("  %-44s FAILED\n", aki_defs[[ov]])))
}

# Quartiles
cat("\n  A2: Mg quartiles (Q4 = ref), KDIGO >=1:\n")
tryCatch({
  dat_a$mg_q <- relevel(factor(dat_a$mg_quartile), ref = "Q4")
  m <- glm(as.formula(paste("aki_kdigo1 ~ mg_q +", cov_rhs)),
           data = dat_a, family = binomial())
  s <- tidy(m, conf.int = TRUE, exponentiate = TRUE) %>% filter(grepl("mg_q", term))
  print(s %>% select(term, estimate, conf.low, conf.high, p.value))
  results_list[["A2_quartiles"]] <<- s
}, error = function(e) cat(sprintf("  Failed: %s\n", e$message)))

# A3: Surgery-type interaction (cardioplegia hypothesis test)
cat("\n  A3: Surgery-type × Mg interaction (cardioplegia hypothesis):\n")
tryCatch({
  dat_a$surg_complex <- ifelse(dat_a$surgery_type %in% c("combined", "valve"),
                                "complex", "simple")
  m_int <- glm(as.formula(paste("aki_kdigo1 ~ first_mg_value * surg_complex +", cov_rhs)),
               data = dat_a, family = binomial())
  s_int <- tidy(m_int, conf.int = TRUE, exponentiate = TRUE)
  main <- s_int %>% filter(term == "first_mg_value")
  inter <- s_int %>% filter(grepl(":", term))
  if (nrow(main) > 0)
    cat(sprintf("    Mg effect (simple surgery): OR %.2f (%.2f-%.2f), p=%.4f\n",
                main$estimate, main$conf.low, main$conf.high, main$p.value))
  if (nrow(inter) > 0)
    cat(sprintf("    Interaction (complex × Mg):  OR %.2f (%.2f-%.2f), p=%.4f\n",
                inter$estimate, inter$conf.low, inter$conf.high, inter$p.value))

  for (stype in c("simple", "complex")) {
    dat_sub <- dat_a %>% filter(surg_complex == stype)
    m_sub <- glm(as.formula(paste("aki_kdigo1 ~ first_mg_value +", cov_rhs)),
                 data = dat_sub, family = binomial())
    s_sub <- tidy(m_sub, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(term == "first_mg_value")
    if (nrow(s_sub) > 0)
      cat(sprintf("    Stratified %s (n=%d): OR %.2f (%.2f-%.2f), p=%.4f\n",
                  stype, nrow(dat_sub), s_sub$estimate, s_sub$conf.low,
                  s_sub$conf.high, s_sub$p.value))
  }
}, error = function(e) cat(sprintf("  Interaction test failed: %s\n", e$message)))

# A4: APACHE sensitivity
cat("\n  A4: APACHE sensitivity (mediator bias test):\n")
if ("apachescore" %in% names(dat_a)) {
  tryCatch({
    m_no <- glm(as.formula(paste("aki_kdigo1 ~ first_mg_value +", cov_rhs)),
                data = dat_a, family = binomial())
    s_no <- tidy(m_no, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(term == "first_mg_value")
    m_yes <- glm(as.formula(paste("aki_kdigo1 ~ first_mg_value +", cov_rhs, "+ apachescore")),
                 data = dat_a %>% filter(!is.na(apachescore)), family = binomial())
    s_yes <- tidy(m_yes, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(term == "first_mg_value")
    if (nrow(s_no) > 0)
      cat(sprintf("    Without APACHE (primary): OR %.2f (%.2f-%.2f), p=%.4f\n",
                  s_no$estimate, s_no$conf.low, s_no$conf.high, s_no$p.value))
    if (nrow(s_yes) > 0)
      cat(sprintf("    With APACHE (sensitivity): OR %.2f (%.2f-%.2f), p=%.4f\n",
                  s_yes$estimate, s_yes$conf.low, s_yes$conf.high, s_yes$p.value))
  }, error = function(e) cat(sprintf("  APACHE sensitivity failed: %s\n", e$message)))
}


# =====================================================================
# TTE SUITE (runs for both unrestricted and restricted)
# =====================================================================
run_suite <- function(dat, dat_m, title, pfx) {
  if (is.null(dat) || nrow(dat) < 20) {
    cat(sprintf("\n  Skipping %s\n", title)); return() }
  dat <- add_aki(dat)
  if (!is.null(dat_m)) dat_m <- add_aki(dat_m)
  n_total <- nrow(dat)  # for E-value rare-outcome check

  cat("\n", strrep("=", 70), "\n")
  cat(sprintf("%s  [N=%d, trt=%d]\n", title, nrow(dat), sum(dat$trt)))
  cat(strrep("=", 70), "\n")

  # ── A*. PIPELINE CHECK: Serum Mg elevation ──────────────────────
  cat("\n  A*. PIPELINE CHECK: Serum Mg elevation\n")
  if ("delta_mg" %in% names(dat)) {
    dat_mg <- dat %>% filter(!is.na(delta_mg))
    if (nrow(dat_mg) > 20) {
      tryCatch({
        if ("hospitalid" %in% names(dat_mg)) {
          des <- svydesign(ids = ~hospitalid, weights = ~iptw, data = dat_mg)
        } else {
          des <- svydesign(ids = ~1, weights = ~iptw, data = dat_mg)
        }
        m <- svyglm(delta_mg ~ trt, design = des)
        s <- tidy(m, conf.int = TRUE) %>% filter(term == "trt")
        if (nrow(s) > 0) {
          cat(sprintf("    IPTW mean diff = %.3f mg/dL (%.3f to %.3f), p = %.4f\n",
                      s$estimate, s$conf.low, s$conf.high, s$p.value))
          results_list[[paste0(pfx, "_mg_elevation")]] <<- tibble(
            model = paste0(pfx, "_mg_elevation"), term = "trt",
            estimate = s$estimate, conf.low = s$conf.low, conf.high = s$conf.high,
            p.value = s$p.value, n_events = nrow(dat_mg), label = "Serum Mg elevation (pipeline)")
        }
      }, error = function(e) cat(sprintf("    Failed: %s\n", e$message)))
    }
  }

  # ── A**. POSITIVE CONTROL: Ventricular arrhythmia (RR 0.52) ───
  cat("\n  A**. POSITIVE CONTROL: Ventricular arrhythmia (RCT RR=0.52)\n")
  r_va <- iptw_or(dat, "vent_arrhythmia", "VT/VF (positive control)")
  if (!is.null(r_va)) {
    cat(sprintf("    IPTW OR = %.2f (%.2f-%.2f), p = %.4f  [n_events=%d]\n",
                r_va$estimate, r_va$conf.low, r_va$conf.high, r_va$p.value, r_va$n_events))
    results_list[[paste0(pfx, "_vent_arrhythmia")]] <<- r_va
  } else {
    cat("    Insufficient VT/VF events\n")
  }

  # ── B. POAF (secondary outcome) ──────────────────────────────────
  r_all  <- iptw_or(dat, "poaf", "POAF dx-only (primary)")
  r_cv   <- iptw_or(dat, "poaf_cardioversion", "POAF cardioversion (Yan)")
  r_comp <- iptw_or(dat, "poaf_composite", "POAF composite (sensitivity)")
  r_nobb <- iptw_or(dat %>% filter(has_betablocker == 0), "poaf", "POAF dx (no BB)")
  r_bb   <- iptw_or(dat %>% filter(has_betablocker == 1), "poaf", "POAF dx (with BB)")
  ptbl(list(r_all, r_cv, r_comp, r_nobb, r_bb), "B. POAF (secondary)")
  for (r in list(r_all, r_nobb, r_bb))
    if (!is.null(r)) results_list[[paste0(pfx, "_", r$label)]] <<- r

  # ── C. AKI severity-stratified ─────────────────────────────────
  aki_oc <- c(aki_delta03="Delta >=0.3", aki_kdigo1="KDIGO >=1",
              aki_primary="Ratio >=1.5x", aki_stage2="Stage >=2",
              aki_stage3="Stage >=3")
  aki_r <- lapply(names(aki_oc), function(o) iptw_or(dat, o, aki_oc[[o]]))
  ptbl(aki_r, "C. PRIMARY: AKI severity-stratified")
  for (r in aki_r) if (!is.null(r))
    results_list[[paste0(pfx, "_aki_", r$label)]] <<- r

  # C1b. Time-windowed AKI (ratio >=1.5x within 24h, 48h, 72h)
  tw_oc <- c(aki_primary_24h="AKI 1.5x <=24h", aki_primary_48h="AKI 1.5x <=48h",
             aki_primary_72h="AKI 1.5x <=72h")
  tw_r <- lapply(names(tw_oc), function(o) iptw_or(dat, o, tw_oc[[o]]))
  ptbl(tw_r, "C1b. Time-windowed AKI")
  for (r in tw_r) if (!is.null(r))
    results_list[[paste0(pfx, "_tw_", r$label)]] <<- r

  # C1c. Peak Cr ratio (continuous outcome)
  if ("peak_cr_ratio" %in% names(dat)) {
    dat_pcr <- dat %>% filter(!is.na(peak_cr_ratio))
    if (nrow(dat_pcr) > 20) {
      tryCatch({
        des <- if ("hospitalid" %in% names(dat_pcr))
          svydesign(ids=~hospitalid, weights=~iptw, data=dat_pcr)
        else svydesign(ids=~1, weights=~iptw, data=dat_pcr)
        m <- svyglm(peak_cr_ratio ~ trt, design=des)
        s <- tidy(m, conf.int=TRUE) %>% filter(term=="trt")
        if (nrow(s) > 0) {
          cat(sprintf("\n  C1c. Peak Cr ratio diff = %.3f (%.3f to %.3f), p=%.4f\n",
                      s$estimate, s$conf.low, s$conf.high, s$p.value))
          results_list[[paste0(pfx, "_peak_cr")]] <<- tibble(
            model=paste0(pfx,"_peak_cr"), term="trt", estimate=s$estimate,
            conf.low=s$conf.low, conf.high=s$conf.high, p.value=s$p.value,
            n_events=nrow(dat_pcr), label="Peak Cr ratio diff")
        }
      }, error = function(e) cat(sprintf("  Peak Cr failed: %s\n", e$message)))
    }
  }

  # AKI x BB
  r_nobb <- iptw_or(dat %>% filter(has_betablocker==0), "aki_primary", "AKI 1.5x (no BB)")
  r_bb   <- iptw_or(dat %>% filter(has_betablocker==1), "aki_primary", "AKI 1.5x (with BB)")
  ptbl(list(r_nobb, r_bb), "C2. AKI x beta-blocker")
  for (r in list(r_nobb, r_bb))
    if (!is.null(r)) results_list[[paste0(pfx, "_", r$label)]] <<- r

  # Sensitivity (Cox, OW, matching) for ratio >=1.5x
  cat("\n  C3. Sensitivity (ratio >=1.5x):\n")
  tryCatch({
    ds <- dat %>% filter(!is.na(time_to_event_hours), time_to_event_hours > 0)
    if (nrow(ds) > 10) {
      m <- coxph(Surv(time_to_event_hours, aki_primary) ~ trt,
                 weights = iptw, data = ds, robust = TRUE)
      s <- tidy(m, conf.int = TRUE, exponentiate = TRUE)
      cat(sprintf("    Cox HR = %.2f (%.2f-%.2f), p=%.4f\n",
                  s$estimate[1], s$conf.low[1], s$conf.high[1], s$p.value[1]))
      results_list[[paste0(pfx, "_cox")]] <<- s
    }
  }, error = function(e) cat(sprintf("    Cox: %s\n", e$message)))
  tryCatch({
    des <- svydesign(ids=~1, weights=~ow, data=dat)
    m <- svyglm(aki_primary ~ trt, design=des, family=quasibinomial())
    s <- tidy(m, conf.int=TRUE, exponentiate=TRUE) %>% filter(term=="trt")
    if (nrow(s)>0) cat(sprintf("    OW OR = %.2f (%.2f-%.2f), p=%.4f\n",
                               s$estimate, s$conf.low, s$conf.high, s$p.value))
    results_list[[paste0(pfx, "_ow")]] <<- s
  }, error = function(e) cat(sprintf("    OW: %s\n", e$message)))
  if (!is.null(dat_m) && nrow(dat_m) > 4) tryCatch({
    m <- glm(aki_primary ~ trt, data=dat_m, family=binomial(), weights=weights)
    s <- tidy(m, conf.int=TRUE, exponentiate=TRUE) %>% filter(term=="trt")
    if (nrow(s)>0) cat(sprintf("    Matched OR = %.2f (%.2f-%.2f), p=%.4f\n",
                               s$estimate, s$conf.low, s$conf.high, s$p.value))
    results_list[[paste0(pfx, "_matched")]] <<- s
  }, error = function(e) cat(sprintf("    Match: %s\n", e$message)))

  # ── D. NEGATIVE CONTROLS ───────────────────────────────────────
  nc_oc <- c(nc_fracture="Fracture", nc_skin_infection="Skin infection", nc_uti="UTI")
  nc_r <- lapply(names(nc_oc), function(o) iptw_or(dat, o, nc_oc[[o]]))
  ptbl(nc_r, "D. NEGATIVE CONTROLS")
  for (r in nc_r) if (!is.null(r)) results_list[[paste0(pfx, "_nc_", r$label)]] <<- r

  # ── E. NEURO (exploratory) ─────────────────────────────────────
  neuro_oc <- c(neuro_delirium="Delirium", neuro_seizure="Seizure",
                neuro_stroke_postop="Stroke (postop)", neuro_encephalopathy="Encephalopathy")
  neuro_r <- lapply(names(neuro_oc), function(o) iptw_or(dat, o, neuro_oc[[o]]))
  ptbl(neuro_r, "E. NEURO (exploratory)")
  for (r in neuro_r) if (!is.null(r)) results_list[[paste0(pfx, "_neuro_", r$label)]] <<- r

  # ── F. SECONDARY ───────────────────────────────────────────────
  sec_oc <- c(rrt_7d="RRT 7d", icu_mortality="ICU mortality", hosp_mortality="Hospital mortality")
  sec_r <- lapply(names(sec_oc), function(o) iptw_or(dat, o, sec_oc[[o]]))
  ptbl(sec_r, "F. SECONDARY")
  for (r in sec_r) if (!is.null(r)) results_list[[paste0(pfx, "_sec_", r$label)]] <<- r

  # ── G. E-VALUES (unmeasured confounding sensitivity) ───────────
  # FIX (bug #1): compute directly from OR/CI, not from n_total column
  if (HAS_EVALUE) {
    cat("\n  G. E-VALUES:\n")
    evalue_targets <- c("KDIGO >=1", "Ratio >=1.5x", "AKI 1.5x <=48h",
                        "Hospital mortality")
    for (nm in evalue_targets) {
      for (k in names(results_list)) {
        obj <- results_list[[k]]
        if (!is.null(obj) && is.data.frame(obj) && "label" %in% names(obj) &&
            any(obj$label == nm) && grepl(pfx, k)) {
          row <- obj[obj$label == nm, ]
          if (nrow(row) > 0 && !is.na(row$estimate[1])) {
            tryCatch({
              rare <- (row$n_events[1] / n_total) < 0.15
              ev <- EValue::evalues.OR(row$estimate[1], row$conf.low[1],
                                       row$conf.high[1], rare = rare)
              cat(sprintf("    %s: E-value = %.2f (CI bound: %.2f)\n",
                          nm, ev["E-values", "point"], ev["E-values", "lower"]))
            }, error = function(e) NULL)
          }
          break
        }
      }
    }
  }
}

# ─── Run both TTE designs ───────────────────────────────────────────
run_suite(dat_tte_b, dat_m_b,
          "TTE-B: ALL PATIENTS (unrestricted)", "TTEB")
run_suite(dat_tte_a, NULL,
          "TTE-A: HYPOMAGNESEMIA ONLY (Mg < 2.0)", "TTEA")


# =====================================================================
# SAVE
# =====================================================================
cat("\n", strrep("=", 70), "\n")
cat("RESULTS SUMMARY\n")
cat(strrep("=", 70), "\n")

rows <- list()
for (nm in names(results_list)) {
  obj <- results_list[[nm]]
  if (is.data.frame(obj) && "estimate" %in% names(obj))
    for (i in seq_len(nrow(obj)))
      rows[[length(rows)+1]] <- tibble(
        model = nm,
        term = if ("label" %in% names(obj)) obj$label[i] else obj$term[i],
        estimate = obj$estimate[i],
        conf.low = if ("conf.low" %in% names(obj)) obj$conf.low[i] else NA,
        conf.high = if ("conf.high" %in% names(obj)) obj$conf.high[i] else NA,
        p.value = if ("p.value" %in% names(obj)) obj$p.value[i] else NA,
        n_events = if ("n_events" %in% names(obj)) obj$n_events[i] else NA)
}
if (length(rows) > 0) {
  df <- bind_rows(rows)
  write_csv(df, file.path(RESULTS, "03_results_summary.csv"))
  print(df, n = 100)
  cat(sprintf("\nSaved: %s\n", file.path(RESULTS, "03_results_summary.csv")))
}
cat("\n03_models.R COMPLETE\n")
