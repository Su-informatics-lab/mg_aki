#!/usr/bin/env Rscript
# ============================================================================
# fig_pairwise_cr.R — Pair-level Cr visualization (Dr. Su)
#
# Panel a: Scatterplot — each dot = one matched pair
#          X = control ΔCr, Y = treated ΔCr
#          Points BELOW diagonal = Mg protective for that pair
#
# Panel b: Histogram of pair-level ΔΔCr (= ΔCr_treated - ΔCr_control)
#          Distribution shifted left of 0 = overall protection
#
# Panel c-f: Mean ΔCr trajectory stratified by eGFR
#            (<45, 45-60, 60-90, >=90)
#
# Run: Rscript fig_pairwise_cr.R eicu
#      Rscript fig_pairwise_cr.R mimic
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

RESULTS <- path.expand("~/mg_aki/results")
WONG <- c(blue="#0072B2", vermil="#D55E00", green="#009E73",
          orange="#E69F00", skyblue="#56B4E9", purple="#CC79A7")

theme_nature <- function(base_size=7) {
  theme_classic(base_size=base_size, base_family="sans") %+replace%
    theme(axis.line=element_line(linewidth=0.3),
          axis.ticks=element_line(linewidth=0.3),
          axis.text=element_text(size=6, color="black"),
          axis.title=element_text(size=7),
          plot.title=element_text(size=8, face="bold", hjust=0),
          panel.grid=element_blank(),
          plot.margin=margin(3,3,3,3,"mm"))
}

# ============================================================================
make_pairwise <- function(db) {
  tag <- tolower(db)
  db_label <- ifelse(tag=="eicu", "eICU-CRD", "MIMIC-IV")
  cat(sprintf("\n== %s: Pairwise Cr ==\n", db_label))

  df <- read.csv(file.path(RESULTS, sprintf("did_matched_%s_24h.csv", tag)),
                  stringsAsFactors=FALSE)

  # ── Build pair-level data ────────────────────────────────────────────
  trt <- df[df$treated == 1, c("match_pair_id","delta_cr","cr_pre","egfr",
                                 "age","is_female","diabetes","ckd")]
  ctl <- df[df$treated == 0, c("match_pair_id","delta_cr")]

  # Aggregate controls per pair (some pairs have multiple controls from replacement)
  ctl_agg <- aggregate(delta_cr ~ match_pair_id, data=ctl, FUN=mean)

  pairs <- merge(trt, ctl_agg, by="match_pair_id", suffixes=c("_trt","_ctl"))
  pairs$dd_cr <- pairs$delta_cr_trt - pairs$delta_cr_ctl  # pair-level DiD

  # eGFR strata
  pairs$egfr_group <- cut(pairs$egfr,
                           breaks=c(0, 45, 60, 90, Inf),
                           labels=c("<45","45-60","60-90",">=90"),
                           right=FALSE)

  cat(sprintf("  Pairs: %d\n", nrow(pairs)))
  cat(sprintf("  Mean DDCr: %.4f (median: %.4f)\n",
              mean(pairs$dd_cr, na.rm=TRUE), median(pairs$dd_cr, na.rm=TRUE)))
  cat(sprintf("  %% pairs with DDCr < 0 (Mg protective): %.1f%%\n",
              100 * mean(pairs$dd_cr < 0, na.rm=TRUE)))

  # ── Panel (a): Scatterplot of pair ΔCr ──────────────────────────────
  lim <- c(-1.5, 2)  # symmetric enough to show the diagonal

  pa <- ggplot(pairs, aes(x=delta_cr_ctl, y=delta_cr_trt)) +
    geom_abline(slope=1, intercept=0, color="grey60", linewidth=0.3, linetype="dashed") +
    geom_point(alpha=0.15, size=0.8, color=WONG[["blue"]]) +
    geom_smooth(method="lm", se=TRUE, color=WONG[["vermil"]], linewidth=0.8,
                fill=WONG[["vermil"]], alpha=0.15) +
    coord_cartesian(xlim=lim, ylim=lim) +
    labs(x="Control DCr (mg/dL)",
         y="IV Mg DCr (mg/dL)",
         title=sprintf("(a) %s: Pair-level Cr change", db_label)) +
    annotate("text", x=lim[2]*0.5, y=lim[1]*0.6, size=2.2, color="grey40",
             label="Below diagonal = Mg protective", hjust=0.5) +
    annotate("text", x=0.05, y=1.8, hjust=0, size=2.2, color=WONG[["vermil"]],
             label=sprintf("Mean DDCr = %.3f mg/dL", mean(pairs$dd_cr, na.rm=T))) +
    theme_nature()

  # ── Panel (b): Histogram of pair-level ΔΔCr ────────────────────────
  pb <- ggplot(pairs, aes(x=dd_cr)) +
    geom_histogram(binwidth=0.05, fill=WONG[["blue"]], color="white",
                   linewidth=0.1, alpha=0.7) +
    geom_vline(xintercept=0, color="grey60", linewidth=0.3, linetype="dashed") +
    geom_vline(xintercept=mean(pairs$dd_cr, na.rm=T),
               color=WONG[["vermil"]], linewidth=0.8) +
    coord_cartesian(xlim=c(-1.5, 1.5)) +
    labs(x="Pair-level DDCr (mg/dL)",
         y="Number of pairs",
         title=sprintf("(b) %s: DDCr distribution", db_label)) +
    annotate("text", x=mean(pairs$dd_cr,na.rm=T)-0.05, y=Inf, vjust=2, hjust=1,
             size=2.2, color=WONG[["vermil"]],
             label=sprintf("Mean = %.3f", mean(pairs$dd_cr,na.rm=T))) +
    theme_nature()

  # ── Panels (c-f): eGFR-stratified mean ΔCr trajectory ──────────────
  # Load cr_all for trajectories
  cr_all <- read.csv(file.path(RESULTS, sprintf("did_cr_all_%s.csv", tag)),
                      stringsAsFactors=FALSE)
  id_col <- if ("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
  cr_all$pid <- cr_all[[id_col]]

  cr_sub <- cr_all[cr_all$pid %in% df$pid, ]
  cr_sub <- merge(cr_sub, df[, c("pid","cr_pre","treated","egfr")], by="pid")
  cr_sub$dcr <- cr_sub$labresult - cr_sub$cr_pre
  cr_sub$hour_bin <- floor(cr_sub$labresultoffset / 60 / 6) * 6
  cr_sub <- cr_sub[cr_sub$hour_bin >= 0 & cr_sub$hour_bin <= 48, ]
  cr_sub$group <- ifelse(cr_sub$treated==1, "IV Mg", "Control")
  cr_sub$egfr_group <- cut(cr_sub$egfr,
                            breaks=c(0, 45, 60, 90, Inf),
                            labels=c("<45","45-60","60-90",">=90"),
                            right=FALSE)

  egfr_panels <- list()
  egfr_labels <- c("<45","45-60","60-90",">=90")
  panel_letters <- c("c","d","e","f")

  for (i in seq_along(egfr_labels)) {
    eg <- egfr_labels[i]
    d_eg <- cr_sub[!is.na(cr_sub$egfr_group) & cr_sub$egfr_group == eg, ]
    if (nrow(d_eg) < 50) {
      egfr_panels[[i]] <- ggplot() + theme_void() +
        labs(title=sprintf("(%s) eGFR %s — too few", panel_letters[i], eg))
      next
    }

    agg <- aggregate(dcr ~ group + hour_bin, data=d_eg,
                      FUN=function(x) c(m=mean(x), se=sd(x)/sqrt(length(x))))
    agg <- do.call(data.frame, agg)
    names(agg) <- c("group","hour_bin","mean","se")

    n_trt <- length(unique(d_eg$pid[d_eg$treated==1]))
    n_ctl <- length(unique(d_eg$pid[d_eg$treated==0]))

    # Pair-level DDCr for this eGFR stratum
    pairs_eg <- pairs[!is.na(pairs$egfr_group) & pairs$egfr_group == eg, ]
    dd_mean <- if (nrow(pairs_eg)>5) sprintf("DDCr=%.3f", mean(pairs_eg$dd_cr,na.rm=T)) else ""

    egfr_panels[[i]] <- ggplot(agg, aes(x=hour_bin, y=mean, color=group, fill=group)) +
      geom_hline(yintercept=0, color="grey80", linewidth=0.2, linetype="dashed") +
      geom_ribbon(aes(ymin=mean-1.96*se, ymax=mean+1.96*se), alpha=0.15, color=NA) +
      geom_line(linewidth=0.7) +
      geom_point(size=1.2) +
      scale_color_manual(values=c("IV Mg"=WONG[["blue"]], "Control"=WONG[["vermil"]])) +
      scale_fill_manual(values=c("IV Mg"=WONG[["blue"]], "Control"=WONG[["vermil"]])) +
      scale_x_continuous(breaks=seq(0,48,12)) +
      labs(x="Hours from ICU",
           y="DCr (mg/dL)",
           title=sprintf("(%s) eGFR %s (n=%d/%d) %s",
                         panel_letters[i], eg, n_trt, n_ctl, dd_mean),
           color=NULL, fill=NULL) +
      theme_nature() +
      theme(legend.position=if(i==1) c(0.2,0.85) else "none",
            legend.background=element_blank())
  }

  # ── Arrange all panels ──────────────────────────────────────────────
  top_row <- arrangeGrob(pa, pb, ncol=2)
  bot_row <- arrangeGrob(grobs=egfr_panels, ncol=4)
  g <- arrangeGrob(top_row, bot_row, nrow=2, heights=c(1, 0.8))

  out <- sprintf("fig_pairwise_%s", tag)
  ggsave(file.path(RESULTS, paste0(out,".pdf")), g, width=7.2, height=6, device=cairo_pdf)
  ggsave(file.path(RESULTS, paste0(out,".png")), g, width=7.2, height=6, dpi=300)
  cat(sprintf("  Saved: %s.pdf / .png\n", out))
}

# ============================================================================
cat("======================================================================\n")
cat("fig_pairwise_cr.R — Pair-level + eGFR-stratified Cr (Dr. Su)\n")
cat("======================================================================\n")

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) args <- c("eicu","mimic")
for (a in args) make_pairwise(a)
