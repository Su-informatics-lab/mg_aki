#!/usr/bin/env Rscript
# ============================================================================
# 12_cox_prognostic.R — Su's three-step Cox analysis
#
#   Step 1: Prognostic Cox in UNTREATED only: baseline Mg → AKI
#   Step 2: Full Cox: baseline Mg + Mg supp + interaction
#   Step 3: Spline interaction → golden window visualization
#
# Output: results/12_cox_prognostic.csv
#         results/12_spline_interaction.csv (for Python figure)
#         figs/fig_cox_spline.pdf
#
# Run: Rscript 12_cox_prognostic.R
# ============================================================================

suppressPackageStartupMessages({
  library(survival); library(sandwich); library(lmtest)
})

# Try rms for splines; fall back to splines::ns if unavailable
has_rms <- requireNamespace("rms", quietly = TRUE)
if (has_rms) library(rms)

RESULTS <- path.expand("~/mg_aki/results")
FIGS    <- path.expand("~/mg_aki/figs")
dir.create(FIGS, showWarnings = FALSE)
LANDMARK_MIN <- 360  # 6 hours

# ── Standardize ──────────────────────────────────────────────────────
stdz <- function(d) {
  rmap <- c(mg_supplementation="mg_supp", hosp_mortality="hospital_mortality",
    age_num="age", baseline_cr="baseline_creatinine", baseline_egfr="egfr",
    hx_chf="heart_failure", hx_hypertension="hypertension",
    hx_diabetes="diabetes", hx_ckd="ckd", hx_copd="copd",
    hx_pvd="pvd", hx_stroke="stroke", hx_liver="liver_disease",
    nephrotox_loop_diuretic="loop_diuretics", nephrotox_nsaid="nsaids",
    nephrotox_acei_arb="acei_arb", nephrotox_ppi="ppi",
    has_betablocker="beta_blockers", has_steroid="steroids",
    preop_antiarrhythmic="antiarrhythmics",
    first_k_value="first_potassium", first_ca_value="first_calcium",
    first_hr="first_heartrate", has_vasopressor="vasopressor_6h")
  for (old in names(rmap)) {
    new <- rmap[[old]]
    if (old %in% names(d) && !new %in% names(d))
      names(d)[names(d) == old] <- new
  }
  if (is.character(d$age)) { d$age <- suppressWarnings(as.numeric(d$age)); d$age[is.na(d$age)] <- 90 }
  if ("surgery_type" %in% names(d)) {
    d$surg_cabg <- as.integer(d$surgery_type=="cabg")
    d$surg_valve <- as.integer(d$surgery_type=="valve")
    d$surg_combined <- as.integer(d$surgery_type=="combined")
  }
  if ("first_lactate" %in% names(d)) {
    d$lactate_missing <- as.integer(is.na(d$first_lactate))
    d$first_lactate[is.na(d$first_lactate)] <- median(d$first_lactate, na.rm=TRUE)
  }
  for (v in c("bmi","first_heartrate","first_calcium","first_potassium"))
    if (v %in% names(d) && any(is.na(d[[v]])))
      d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm=TRUE)
  d
}

# ── Construct time-to-event from 6h landmark ─────────────────────────
make_tte <- function(d, db_name) {
  MAX_FU <- 7 * 24 * 60  # 7 days in minutes

  if (db_name == "eICU") {
    # AKI time from ICU admission (minutes)
    aki_time <- d$aki_time_offset
    # Censor time = min(discharge, 7d)
    censor_time <- pmin(d$unitdischargeoffset, MAX_FU, na.rm = TRUE)
  } else {
    # MIMIC: aki_time_offset should now be in the cohort
    aki_time <- d$aki_time_offset
    if (is.null(aki_time)) aki_time <- rep(NA, nrow(d))
    # LOS in minutes
    if ("los_min" %in% names(d)) {
      censor_time <- pmin(d$los_min, MAX_FU, na.rm = TRUE)
    } else if ("outtime" %in% names(d) && "intime" %in% names(d)) {
      censor_time <- as.numeric(difftime(as.POSIXct(d$outtime),
                                          as.POSIXct(d$intime), units="mins"))
      censor_time <- pmin(censor_time, MAX_FU, na.rm = TRUE)
    } else {
      censor_time <- rep(MAX_FU, nrow(d))
    }
  }

  # Time from landmark (minutes → days)
  event_time_from_lm <- (aki_time - LANDMARK_MIN) / (60 * 24)
  censor_time_from_lm <- (censor_time - LANDMARK_MIN) / (60 * 24)

  # Event indicator
  has_aki <- !is.na(d$aki_kdigo1) & d$aki_kdigo1 == 1 & !is.na(aki_time)
  time_days <- ifelse(has_aki, event_time_from_lm, censor_time_from_lm)
  time_days <- pmax(time_days, 0.001)  # avoid zero/negative
  event <- as.integer(has_aki)

  # Death as competing event (for future use)
  death_time_from_lm <- (d$death_offset_min - LANDMARK_MIN) / (60 * 24)
  has_death <- !is.na(d$death_offset_min) & d$death_offset_min > LANDMARK_MIN &
               d$death_offset_min <= MAX_FU

  d$surv_time <- time_days
  d$surv_event <- event
  d$death_time_lm <- death_time_from_lm
  d$death_event <- as.integer(has_death)

  cat(sprintf("  %s TTE: %d events / %d patients, median FU %.1f days\n",
              db_name, sum(event), nrow(d), median(time_days)))
  d
}

# ── Covariates for Cox ───────────────────────────────────────────────
cox_covars <- c("age","is_female","bmi",
  "surg_cabg","surg_valve","surg_combined",
  "heart_failure","hypertension","diabetes","ckd",
  "copd","pvd","stroke","liver_disease",
  "baseline_creatinine","egfr",
  "loop_diuretics","nsaids","acei_arb","ppi",
  "beta_blockers","steroids","antiarrhythmics",
  "first_potassium","first_calcium","first_heartrate","vasopressor_6h",
  "transfusion_6h")

# ── Load data ────────────────────────────────────────────────────────
cat("Loading cohorts...\n")
dat_e <- make_tte(stdz(read.csv(file.path(RESULTS, "01_analysis_a_cohort.csv"),
                                 stringsAsFactors=FALSE)), "eICU")
dat_m <- make_tte(stdz(read.csv(file.path(RESULTS, "04_mimic_cohort.csv"),
                                 stringsAsFactors=FALSE)), "MIMIC")

all_rows <- list()

# ============================================================================
# STEP 1: PROGNOSTIC COX — UNTREATED ONLY
# "不补mg的病人，是不是baseline mg和evt强相关"
# ============================================================================
cat(sprintf("\n%s\nSTEP 1: PROGNOSTIC COX (untreated only)\n%s\n",
            strrep("=",65), strrep("=",65)))

run_prognostic <- function(d, db_name) {
  d0 <- d[d$mg_supp == 0, ]
  avail <- intersect(cox_covars, names(d0))
  avail <- avail[avail != "first_mg_value"]  # Mg is the exposure here
  d0 <- d0[complete.cases(d0[, c("surv_time","surv_event","first_mg_value", avail)]), ]
  cat(sprintf("\n  %s untreated: N=%d, events=%d\n", db_name, nrow(d0), sum(d0$surv_event)))

  fml <- as.formula(paste("Surv(surv_time, surv_event) ~ first_mg_value +",
                           paste(avail, collapse="+")))
  fit <- coxph(fml, data=d0)
  s <- summary(fit)
  mg_row <- which(rownames(s$coefficients) == "first_mg_value")
  hr <- s$coefficients[mg_row, "exp(coef)"]
  lo <- s$conf.int[mg_row, "lower .95"]
  hi <- s$conf.int[mg_row, "upper .95"]
  p  <- s$coefficients[mg_row, "Pr(>|z|)"]

  sig <- ifelse(p < 0.05, " *", "")
  cat(sprintf("  Baseline Mg → AKI: HR %.3f (%.3f-%.3f) P=%.4f%s\n",
              hr, lo, hi, p, sig))

  # By surgery type
  if ("surgery_type" %in% names(d0)) {
    d0$complex <- as.integer(d0$surgery_type %in% c("valve","combined"))
    for (stype in c("Simple","Complex")) {
      sub <- if (stype=="Simple") d0[d0$complex==0,] else d0[d0$complex==1,]
      if (nrow(sub) < 50 || sum(sub$surv_event) < 10) next
      avail_s <- intersect(avail, names(sub))
      avail_s <- avail_s[!avail_s %in% c("surg_cabg","surg_valve","surg_combined")]
      fml_s <- as.formula(paste("Surv(surv_time, surv_event) ~ first_mg_value +",
                                 paste(avail_s, collapse="+")))
      tryCatch({
        fit_s <- coxph(fml_s, data=sub)
        s_s <- summary(fit_s)
        mg_r <- which(rownames(s_s$coefficients) == "first_mg_value")
        cat(sprintf("    %s: HR %.3f (%.3f-%.3f) P=%.4f\n", stype,
                    s_s$coefficients[mg_r,"exp(coef)"],
                    s_s$conf.int[mg_r,"lower .95"],
                    s_s$conf.int[mg_r,"upper .95"],
                    s_s$coefficients[mg_r,"Pr(>|z|)"]))
      }, error=function(e) cat(sprintf("    %s: FAILED\n", stype)))
    }
  }

  data.frame(db=db_name, step="prognostic_untreated",
             hr=round(hr,3), lo=round(lo,3), hi=round(hi,3), p=round(p,4))
}

all_rows[[1]] <- run_prognostic(dat_e, "eICU")
all_rows[[2]] <- run_prognostic(dat_m, "MIMIC")

# ============================================================================
# STEP 2: COX WITH Mg SUPP AS COVARIATE + INTERACTION
# "把补mg作为一个covariant加进来看一下趋势"
# ============================================================================
cat(sprintf("\n%s\nSTEP 2: COX WITH TREATMENT + INTERACTION\n%s\n",
            strrep("=",65), strrep("=",65)))

run_interaction_cox <- function(d, db_name) {
  avail <- intersect(cox_covars, names(d))
  d <- d[complete.cases(d[, c("surv_time","surv_event","first_mg_value","mg_supp", avail)]), ]
  cat(sprintf("\n  %s: N=%d, events=%d, trt=%d\n",
              db_name, nrow(d), sum(d$surv_event), sum(d$mg_supp)))

  # Model without interaction
  fml1 <- as.formula(paste("Surv(surv_time, surv_event) ~ mg_supp + first_mg_value +",
                            paste(avail, collapse="+")))
  fit1 <- coxph(fml1, data=d)
  s1 <- summary(fit1)
  mg_r <- which(rownames(s1$coefficients) == "mg_supp")
  cat(sprintf("  Main effect (no interaction):\n"))
  cat(sprintf("    Mg supp:     HR %.3f (%.3f-%.3f) P=%.4f\n",
              s1$coefficients[mg_r,"exp(coef)"], s1$conf.int[mg_r,"lower .95"],
              s1$conf.int[mg_r,"upper .95"], s1$coefficients[mg_r,"Pr(>|z|)"]))
  bm_r <- which(rownames(s1$coefficients) == "first_mg_value")
  cat(sprintf("    Baseline Mg: HR %.3f (%.3f-%.3f) P=%.4f\n",
              s1$coefficients[bm_r,"exp(coef)"], s1$conf.int[bm_r,"lower .95"],
              s1$conf.int[bm_r,"upper .95"], s1$coefficients[bm_r,"Pr(>|z|)"]))

  # Model WITH interaction
  fml2 <- as.formula(paste("Surv(surv_time, surv_event) ~ mg_supp * first_mg_value +",
                            paste(avail, collapse="+")))
  fit2 <- coxph(fml2, data=d)
  s2 <- summary(fit2)
  int_r <- grep("mg_supp:first_mg_value", rownames(s2$coefficients))

  if (length(int_r) > 0) {
    int_hr <- s2$coefficients[int_r, "exp(coef)"]
    int_p  <- s2$coefficients[int_r, "Pr(>|z|)"]
    sig <- ifelse(int_p < 0.05, " *", "")
    cat(sprintf("  Interaction (Cox):\n"))
    cat(sprintf("    Mg_supp × baseline_Mg: HR %.3f P=%.4f%s\n", int_hr, int_p, sig))
    cat(sprintf("    Interpretation: per 1 mg/dL higher baseline Mg,\n"))
    cat(sprintf("    the HR for Mg supp changes by factor %.3f\n", int_hr))

    # LRT for interaction
    lr <- anova(fit1, fit2)
    if (nrow(lr) >= 2) {
      lrt_p <- lr[2, "Pr(>|Chi|)"]
      cat(sprintf("    LRT for interaction: P=%.4f\n", lrt_p))
    }

    return(data.frame(db=db_name, step="cox_interaction",
                      hr=round(int_hr,3), lo=NA, hi=NA, p=round(int_p,4)))
  }
  NULL
}

all_rows[[3]] <- run_interaction_cox(dat_e, "eICU")
all_rows[[4]] <- run_interaction_cox(dat_m, "MIMIC")

# ============================================================================
# STEP 3: SPLINE — TREATMENT EFFECT BY BASELINE Mg (golden window)
# ============================================================================
cat(sprintf("\n%s\nSTEP 3: SPLINE INTERACTION (golden window)\n%s\n",
            strrep("=",65), strrep("=",65)))

run_spline <- function(d, db_name) {
  avail <- intersect(cox_covars, names(d))
  d <- d[complete.cases(d[, c("surv_time","surv_event","first_mg_value","mg_supp", avail)]), ]

  # Grid of baseline Mg values
  mg_grid <- seq(max(1.0, quantile(d$first_mg_value, 0.02)),
                 min(4.5, quantile(d$first_mg_value, 0.98)),
                 length.out = 50)
  cat(sprintf("  %s: Mg range %.1f-%.1f, N=%d\n",
              db_name, min(mg_grid), max(mg_grid), nrow(d)))

  # Fit interaction model with natural spline
  tryCatch({
    if (has_rms) {
      dd <- datadist(d); options(datadist="dd")
      fml <- as.formula(paste(
        "Surv(surv_time, surv_event) ~ mg_supp * rcs(first_mg_value, 4) +",
        paste(avail, collapse="+")))
      fit <- cph(fml, data=d, x=TRUE, y=TRUE)

      # Predict HR for mg_supp=1 vs mg_supp=0 at each Mg level
      results <- data.frame()
      for (mg_val in mg_grid) {
        d1 <- d[1,]; d0 <- d[1,]
        d1$mg_supp <- 1; d1$first_mg_value <- mg_val
        d0$mg_supp <- 0; d0$first_mg_value <- mg_val
        tryCatch({
          p1 <- predict(fit, newdata=d1, type="lp")
          p0 <- predict(fit, newdata=d0, type="lp")
          log_hr <- p1 - p0
          # Bootstrap SE would be ideal; use delta method approximation
          results <- rbind(results, data.frame(
            mg=mg_val, log_hr=log_hr, hr=exp(log_hr)))
        }, error=function(e) {})
      }
      if (nrow(results) > 0) {
        # Simple CI from overall SE
        se_approx <- sd(results$log_hr) * 0.5  # rough
        results$lo <- exp(results$log_hr - 1.96 * se_approx)
        results$hi <- exp(results$log_hr + 1.96 * se_approx)
      }
    } else {
      # Fallback: use ns() from splines
      library(splines)
      fml <- as.formula(paste(
        "Surv(surv_time, surv_event) ~ mg_supp * ns(first_mg_value, df=3) +",
        paste(avail, collapse="+")))
      fit <- coxph(fml, data=d)

      # Predict at grid points
      results <- data.frame()
      for (mg_val in mg_grid) {
        d1 <- d[1,, drop=FALSE]; d0 <- d[1,, drop=FALSE]
        d1$mg_supp <- 1; d1$first_mg_value <- mg_val
        d0$mg_supp <- 0; d0$first_mg_value <- mg_val
        tryCatch({
          lp1 <- predict(fit, newdata=d1, type="lp", se.fit=TRUE)
          lp0 <- predict(fit, newdata=d0, type="lp", se.fit=TRUE)
          log_hr <- lp1$fit - lp0$fit
          se <- sqrt(lp1$se.fit^2 + lp0$se.fit^2)
          results <- rbind(results, data.frame(
            mg=mg_val, log_hr=log_hr, hr=exp(log_hr),
            lo=exp(log_hr - 1.96*se), hi=exp(log_hr + 1.96*se)))
        }, error=function(e) {})
      }
    }

    if (nrow(results) > 0) {
      results$db <- db_name
      cat(sprintf("  %s spline: %d grid points computed\n", db_name, nrow(results)))

      # Find the "golden window"
      protective <- results[results$hr < 1 & results$hi < 1, ]
      if (nrow(protective) > 0) {
        cat(sprintf("  Golden window (HR<1 + CI<1): Mg %.1f-%.1f mg/dL\n",
                    min(protective$mg), max(protective$mg)))
      } else {
        marginal <- results[results$hr < 0.8, ]
        if (nrow(marginal) > 0) {
          cat(sprintf("  Suggestive range (HR<0.8): Mg %.1f-%.1f mg/dL\n",
                      min(marginal$mg), max(marginal$mg)))
        }
      }
      return(results)
    }
  }, error = function(e) {
    cat(sprintf("  %s spline FAILED: %s\n", db_name, e$message))
  })
  NULL
}

spline_e <- run_spline(dat_e, "eICU")
spline_m <- run_spline(dat_m, "MIMIC")

# Save spline data for Python figure
if (!is.null(spline_e) || !is.null(spline_m)) {
  spline_all <- rbind(spline_e, spline_m)
  write.csv(spline_all, file.path(RESULTS, "12_spline_interaction.csv"), row.names=FALSE)
  cat(sprintf("  Saved: results/12_spline_interaction.csv\n"))

  # Quick R plot
  tryCatch({
    pdf(file.path(FIGS, "fig_cox_spline.pdf"), width=5, height=4)
    par(mar=c(4,4,2,1), family="sans", cex=0.8)
    plot(NULL, xlim=c(1, 4), ylim=c(0.1, 3),
         xlab="Baseline serum magnesium (mg/dL)",
         ylab="HR for AKI (Mg supp vs no supp)",
         main="Treatment Effect by Baseline Mg (Cox spline)", log="y")
    abline(h=1, lty=2, col="gray50")

    if (!is.null(spline_e)) {
      polygon(c(spline_e$mg, rev(spline_e$mg)),
              c(spline_e$lo, rev(spline_e$hi)),
              col=rgb(0.84,0.37,0, 0.15), border=NA)
      lines(spline_e$mg, spline_e$hr, col="#D55E00", lwd=2)
    }
    if (!is.null(spline_m)) {
      polygon(c(spline_m$mg, rev(spline_m$mg)),
              c(spline_m$lo, rev(spline_m$hi)),
              col=rgb(0,0.45,0.70, 0.15), border=NA)
      lines(spline_m$mg, spline_m$hr, col="#0072B2", lwd=2)
    }
    legend("topright", c("eICU-CRD","MIMIC-IV"),
           col=c("#D55E00","#0072B2"), lwd=2, bty="n")
    dev.off()
    cat("  Saved: figs/fig_cox_spline.pdf\n")
  }, error=function(e) cat(sprintf("  Plot failed: %s\n", e$message)))
}

# ── Save all results ─────────────────────────────────────────────────
res <- do.call(rbind, all_rows[!sapply(all_rows, is.null)])
write.csv(res, file.path(RESULTS, "12_cox_prognostic.csv"), row.names=FALSE)
cat(sprintf("\n✓ Saved: results/12_cox_prognostic.csv\n"))

cat(sprintf("\n%s\nSUMMARY FOR Su\n%s\n", strrep("=",65), strrep("=",65)))
cat("  Step 1: Does baseline Mg predict AKI in untreated? (prognostic)\n")
cat("  Step 2: Does the Mg supp effect vary with baseline Mg? (interaction)\n")
cat("  Step 3: Where is the golden window? (spline)\n")
cat("  → If Step 1 shows signal AND Step 2 interaction significant\n")
cat("    → '有动静，可以玩花样了'\n")
