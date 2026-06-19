#!/usr/bin/env Rscript
# ============================================================================
# 03_did_figures.R — All publication + presentation figures
#
# Fig 1: DiD time course (matching-based, primary)
# Fig 2: KM-style cumulative AKI incidence (overall + key subgroups)
# Fig 3: Subgroup forest plot (AKI KDIGO≥1)
# Fig 4: Secondary outcomes forest plot
# eFig 1: IPTW time course comparison
# eFig 2: Mg-stratified DiD
#
# Run: Rscript 03_did_figures.R
# ============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

RESULTS <- path.expand("~/mg_aki/results")

# ── Nature-style theme ───────────────────────────────────────────────────
WONG <- c(black="#000000", orange="#E69F00", skyblue="#56B4E9",
          green="#009E73", yellow="#F0E442", blue="#0072B2",
          vermil="#D55E00", purple="#CC79A7")

theme_nature <- function(base_size=7) {
  theme_classic(base_size=base_size, base_family="sans") %+replace%
    theme(
      axis.line = element_line(linewidth=0.3),
      axis.ticks = element_line(linewidth=0.3),
      axis.text = element_text(size=6, color="black"),
      axis.title = element_text(size=7),
      plot.title = element_text(size=8, face="bold", hjust=0),
      legend.key.size = unit(0.3, "cm"),
      legend.text = element_text(size=6),
      legend.title = element_text(size=7),
      strip.text = element_text(size=7, face="bold"),
      strip.background = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(3,3,3,3,"mm")
    )
}

save_fig <- function(p, name, w=7.2, h=4) {
  ggsave(file.path(RESULTS, paste0(name, ".pdf")), p, width=w, height=h,
         units="in", device=cairo_pdf)
  ggsave(file.path(RESULTS, paste0(name, ".png")), p, width=w, height=h,
         units="in", dpi=300)
  cat(sprintf("  Saved: %s.pdf / .png\n", name))
}

# ============================================================================
# FIG 1: DiD TIME COURSE (matching-based)
# ============================================================================
fig1_timecourse <- function() {
  cat("\n── Fig 1: DiD Time Course ──\n")

  dfs <- list()
  for (db in c("eicu","mimic")) {
    path <- file.path(RESULTS, sprintf("did_timecourse_%s.csv", db))
    if (!file.exists(path)) { cat("  Missing:", path, "\n"); next }
    d <- read.csv(path, stringsAsFactors=FALSE)
    d$db <- toupper(db)
    d$db_label <- ifelse(db=="eicu", "eICU-CRD", "MIMIC-IV")
    dfs[[db]] <- d
  }
  if (length(dfs)==0) { cat("  No data\n"); return(NULL) }
  tc <- do.call(rbind, dfs)
  tc <- tc[!is.na(tc$did_adj), ]
  tc$tol_label <- factor(sprintf("\u00b1%dh", tc$tol_h), levels=c("\u00b12h","\u00b14h","\u00b16h"))
  tc$sig <- tc$p_adj < 0.05

  p <- ggplot(tc, aes(x=target_h, y=did_adj, color=tol_label, shape=tol_label)) +
    geom_hline(yintercept=0, color="grey70", linewidth=0.3) +
    geom_ribbon(data=tc[tc$tol_h==6,], aes(ymin=ci_lo, ymax=ci_hi, fill=tol_label),
                alpha=0.12, color=NA) +
    geom_line(linewidth=0.6, alpha=0.8) +
    geom_point(size=1.5) +
    geom_point(data=tc[tc$sig,], size=3, shape=8, stroke=0.5) +
    facet_wrap(~db_label, nrow=1) +
    scale_color_manual(values=c(WONG[["skyblue"]], WONG[["orange"]], WONG[["blue"]])) +
    scale_fill_manual(values=c(WONG[["skyblue"]], WONG[["orange"]], WONG[["blue"]])) +
    scale_x_continuous(breaks=seq(6,36,6)) +
    labs(x="Hours after IV MgSO\u2084", y="DiD estimate, \u0394\u0394Cr (mg/dL)",
         color="Temporal tolerance", fill="Temporal tolerance",
         shape="Temporal tolerance") +
    theme_nature() +
    theme(legend.position="bottom")

  save_fig(p, "fig1_did_timecourse", w=7.2, h=3.5)
  return(p)
}

# ============================================================================
# FIG 2: KM-STYLE CUMULATIVE AKI INCIDENCE
# ============================================================================
fig2_km_aki <- function() {
  cat("\n── Fig 2: Cumulative AKI Incidence ──\n")

  plots <- list()
  panel_idx <- 0

  for (db in c("eicu","mimic")) {
    path <- file.path(RESULTS, sprintf("did_matched_%s_24h.csv", db))
    if (!file.exists(path)) { cat("  Missing:", path, "\n"); next }
    df <- read.csv(path, stringsAsFactors=FALSE)

    # Compute time-to-AKI from all Cr measurements
    cr_path <- file.path(RESULTS, sprintf("did_cr_all_%s.csv", db))
    if (!file.exists(cr_path)) next
    cr_all <- read.csv(cr_path, stringsAsFactors=FALSE)
    id_col <- if ("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
    cr_all$pid <- cr_all[[id_col]]

    db_label <- ifelse(db=="eicu", "eICU-CRD", "MIMIC-IV")

    # For each patient in matched set, find time of first ΔCr ≥ 0.3
    time_points <- c(6, 12, 18, 24, 36, 48)
    cum_aki <- list()

    for (grp in c(0, 1)) {
      pids <- df$pid[df$treated == grp]
      n_total <- length(pids)
      cr_sub <- cr_all[cr_all$pid %in% pids, ]

      rates <- numeric(length(time_points))
      for (j in seq_along(time_points)) {
        tp <- time_points[j]
        cr_window <- cr_sub[cr_sub$labresultoffset <= tp * 60, ]
        # Merge baseline Cr
        cr_window <- merge(cr_window, df[df$treated==grp, c("pid","cr_pre")], by="pid")
        cr_window$delta <- cr_window$labresult - cr_window$cr_pre
        aki_pids <- unique(cr_window$pid[cr_window$delta >= 0.3])
        rates[j] <- length(aki_pids) / n_total
      }

      grp_label <- ifelse(grp==1, "IV Mg", "Control")
      cum_aki[[grp_label]] <- data.frame(
        time = time_points, rate = rates * 100,
        group = grp_label, db = db_label, stringsAsFactors = FALSE)
    }

    cum_df <- do.call(rbind, cum_aki)

    # Overall
    panel_idx <- panel_idx + 1
    p <- ggplot(cum_df, aes(x=time, y=rate, color=group, linetype=group)) +
      geom_step(linewidth=0.8) +
      geom_point(size=1.5) +
      scale_color_manual(values=c("IV Mg"=WONG[["blue"]], "Control"=WONG[["vermil"]])) +
      scale_linetype_manual(values=c("IV Mg"="solid", "Control"="dashed")) +
      scale_x_continuous(breaks=time_points) +
      labs(x="Hours from ICU admission", y="Cumulative AKI incidence (%)",
           title=sprintf("(%s) %s — Overall", letters[panel_idx], db_label),
           color=NULL, linetype=NULL) +
      theme_nature() +
      theme(legend.position=c(0.2, 0.85))

    plots[[panel_idx]] <- p

    # Key subgroup: eGFR 60-90 for eICU, Male for MIMIC
    if (db == "eicu") {
      sg_df <- df[!is.na(df$egfr) & df$egfr >= 60 & df$egfr < 90, ]
      sg_label <- "eGFR 60\u201390"
    } else {
      sg_df <- df[!is.na(df$is_female) & df$is_female == 0, ]
      sg_label <- "Males"
    }

    if (nrow(sg_df) >= 40) {
      cum_sg <- list()
      for (grp in c(0, 1)) {
        pids <- sg_df$pid[sg_df$treated == grp]
        n_total <- length(pids)
        if (n_total < 10) next
        cr_sub <- cr_all[cr_all$pid %in% pids, ]
        rates <- numeric(length(time_points))
        for (j in seq_along(time_points)) {
          tp <- time_points[j]
          cr_window <- cr_sub[cr_sub$labresultoffset <= tp * 60, ]
          cr_window <- merge(cr_window, sg_df[sg_df$treated==grp, c("pid","cr_pre")], by="pid")
          cr_window$delta <- cr_window$labresult - cr_window$cr_pre
          aki_pids <- unique(cr_window$pid[cr_window$delta >= 0.3])
          rates[j] <- length(aki_pids) / n_total
        }
        grp_label <- ifelse(grp==1, "IV Mg", "Control")
        cum_sg[[grp_label]] <- data.frame(
          time = time_points, rate = rates * 100,
          group = grp_label, db = db_label, stringsAsFactors = FALSE)
      }
      cum_sg_df <- do.call(rbind, cum_sg)

      panel_idx <- panel_idx + 1
      p2 <- ggplot(cum_sg_df, aes(x=time, y=rate, color=group, linetype=group)) +
        geom_step(linewidth=0.8) +
        geom_point(size=1.5) +
        scale_color_manual(values=c("IV Mg"=WONG[["blue"]], "Control"=WONG[["vermil"]])) +
        scale_linetype_manual(values=c("IV Mg"="solid", "Control"="dashed")) +
        scale_x_continuous(breaks=time_points) +
        labs(x="Hours from ICU admission", y="Cumulative AKI incidence (%)",
             title=sprintf("(%s) %s — %s", letters[panel_idx], db_label, sg_label),
             color=NULL, linetype=NULL) +
        theme_nature() +
        theme(legend.position=c(0.2, 0.85))

      plots[[panel_idx]] <- p2
    }
  }

  if (length(plots) > 0) {
    g <- arrangeGrob(grobs=plots, ncol=2)
    ggsave(file.path(RESULTS, "fig2_km_aki.pdf"), g, width=7.2, height=5, device=cairo_pdf)
    ggsave(file.path(RESULTS, "fig2_km_aki.png"), g, width=7.2, height=5, dpi=300)
    cat("  Saved: fig2_km_aki.pdf / .png\n")
  }
}

# ============================================================================
# FIG 3: SUBGROUP FOREST PLOT (AKI KDIGO ≥1)
# ============================================================================
fig3_forest <- function() {
  cat("\n── Fig 3: Subgroup Forest Plot ──\n")

  dfs <- list()
  for (db in c("eicu","mimic")) {
    path <- file.path(RESULTS, sprintf("did_aki_subgroups_%s.csv", db))
    if (!file.exists(path)) next
    d <- read.csv(path, stringsAsFactors=FALSE)
    d$db <- ifelse(db=="eicu", "eICU-CRD", "MIMIC-IV")
    dfs[[db]] <- d
  }
  if (length(dfs)==0) { cat("  No data\n"); return(NULL) }
  all <- do.call(rbind, dfs)
  k1 <- all[all$endpoint == "aki_kdigo1" & !is.na(all$or) & all$or < 100, ]

  # Select key subgroups for forest
  keep <- c("Overall", "Age < 65", "Age >= 65", "Female", "Male",
            "eGFR < 60", "eGFR 60-90", "eGFR >= 90",
            "Diabetes", "No diabetes", "CKD", "No CKD",
            "CABG", "Valve", "Mg < 1.8", "Mg >= 2.0")
  k1 <- k1[k1$subgroup %in% keep, ]
  k1$subgroup <- factor(k1$subgroup, levels=rev(keep))
  k1$sig <- k1$p < 0.05

  p <- ggplot(k1, aes(x=or, y=subgroup, color=db)) +
    geom_vline(xintercept=1, color="grey70", linewidth=0.3) +
    geom_errorbarh(aes(xmin=or_lo, xmax=or_hi), height=0.2, linewidth=0.4,
                   position=position_dodge(width=0.5)) +
    geom_point(aes(shape=sig), size=2, position=position_dodge(width=0.5)) +
    scale_color_manual(values=c("eICU-CRD"=WONG[["blue"]], "MIMIC-IV"=WONG[["vermil"]])) +
    scale_shape_manual(values=c("TRUE"=16, "FALSE"=1), guide="none") +
    scale_x_log10(limits=c(0.1, 10), breaks=c(0.1, 0.25, 0.5, 1, 2, 4)) +
    labs(x="Odds Ratio (AKI KDIGO \u2265 Stage 1)", y=NULL, color=NULL) +
    theme_nature() +
    theme(legend.position="bottom")

  save_fig(p, "fig3_forest_aki", w=5, h=5)
  return(p)
}

# ============================================================================
# FIG 4: SECONDARY OUTCOMES FOREST PLOT
# ============================================================================
fig4_secondary <- function() {
  cat("\n── Fig 4: Secondary Outcomes Forest ──\n")

  dfs <- list()
  for (db in c("eicu","mimic")) {
    path <- file.path(RESULTS, sprintf("did_secondary_%s.csv", db))
    if (!file.exists(path)) next
    d <- read.csv(path, stringsAsFactors=FALSE)
    d$db <- ifelse(db=="eicu", "eICU-CRD", "MIMIC-IV")
    dfs[[db]] <- d
  }
  if (length(dfs)==0) { cat("  No data\n"); return(NULL) }
  all <- do.call(rbind, dfs)
  all <- all[!is.na(all$or) & all$or < 100, ]
  all$label <- c(hosp_mortality="Mortality", poaf="POAF",
                  encephalopathy="Encephalopathy",
                  vent_arrhythmia="Vent. arrhythmia")[all$outcome]
  all$label <- factor(all$label, levels=rev(c("Mortality","POAF","Encephalopathy","Vent. arrhythmia")))

  p <- ggplot(all, aes(x=or, y=label, color=db)) +
    geom_vline(xintercept=1, color="grey70", linewidth=0.3) +
    geom_errorbarh(aes(xmin=or_lo, xmax=or_hi), height=0.2, linewidth=0.5,
                   position=position_dodge(width=0.5)) +
    geom_point(size=2.5, position=position_dodge(width=0.5)) +
    scale_color_manual(values=c("eICU-CRD"=WONG[["blue"]], "MIMIC-IV"=WONG[["vermil"]])) +
    scale_x_log10(limits=c(0.1, 10), breaks=c(0.1, 0.25, 0.5, 1, 2, 4)) +
    labs(x="Odds Ratio (95% CI)", y=NULL, color=NULL,
         title="Secondary Outcomes (matched comparison)") +
    theme_nature() +
    theme(legend.position="bottom")

  save_fig(p, "fig4_secondary_forest", w=4.7, h=2.5)
  return(p)
}

# ============================================================================
# eFIG 1: IPTW TIME COURSE
# ============================================================================
efig1_iptw <- function() {
  cat("\n── eFig 1: IPTW Time Course ──\n")

  dfs <- list()
  for (db in c("eicu","mimic")) {
    path <- file.path(RESULTS, sprintf("did_iptw_%s.csv", db))
    if (!file.exists(path)) next
    d <- read.csv(path, stringsAsFactors=FALSE)
    d$db_label <- ifelse(db=="eicu", "eICU-CRD", "MIMIC-IV")
    dfs[[db]] <- d
  }
  if (length(dfs)==0) { cat("  No data\n"); return(NULL) }
  all <- do.call(rbind, dfs)

  method_colors <- c(sIPTW=WONG[["skyblue"]], sIPTW_t99=WONG[["orange"]],
                      sIPTW_t95=WONG[["green"]], sIPTW_DR=WONG[["blue"]],
                      AIPW=WONG[["vermil"]])

  p <- ggplot(all, aes(x=target_h, y=did, color=method, linetype=method)) +
    geom_hline(yintercept=0, color="grey70", linewidth=0.3) +
    geom_line(linewidth=0.5) +
    geom_point(size=1) +
    facet_wrap(~db_label, nrow=1) +
    scale_color_manual(values=method_colors) +
    scale_linetype_manual(values=c(sIPTW="dotted", sIPTW_t99="dashed",
                                    sIPTW_t95="dashed", sIPTW_DR="solid",
                                    AIPW="solid")) +
    scale_x_continuous(breaks=seq(6,36,6)) +
    labs(x="Hours from ICU admission", y="DiD estimate, \u0394\u0394Cr (mg/dL)",
         color="IPTW variant", linetype="IPTW variant") +
    theme_nature() +
    theme(legend.position="bottom")

  save_fig(p, "efig1_iptw_timecourse", w=7.2, h=3.5)
  return(p)
}

# ============================================================================
# RUN ALL
# ============================================================================
cat("======================================================================\n")
cat("03_did_figures.R — Publication + presentation figures\n")
cat("======================================================================\n")

fig1_timecourse()
fig2_km_aki()
fig3_forest()
fig4_secondary()
efig1_iptw()

cat("\n======================================================================\n")
cat("All figures saved to", RESULTS, "\n")
cat("======================================================================\n")
