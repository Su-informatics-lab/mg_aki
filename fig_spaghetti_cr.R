#!/usr/bin/env Rscript
# ============================================================================
# fig_spaghetti_cr.R — Individual Cr trajectories aligned to IV Mg time
#
# Dr. Su's request: show matched patients' raw Cr changes
# Panel A: IV Mg group (treated) — each line = one patient
# Panel B: Matched controls — aligned to partner's IV Mg time
# Same X/Y range for direct visual comparison
# Thin semi-transparent individual lines + bold mean LOESS
#
# Run: Rscript fig_spaghetti_cr.R eicu
#      Rscript fig_spaghetti_cr.R mimic
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

RESULTS <- path.expand("~/mg_aki/results")

WONG <- c(blue="#0072B2", vermil="#D55E00", grey="#999999")

theme_nature <- function(base_size=7) {
  theme_classic(base_size=base_size, base_family="sans") %+replace%
    theme(
      axis.line=element_line(linewidth=0.3),
      axis.ticks=element_line(linewidth=0.3),
      axis.text=element_text(size=6, color="black"),
      axis.title=element_text(size=7),
      plot.title=element_text(size=8, face="bold", hjust=0),
      legend.key.size=unit(0.3,"cm"),
      legend.text=element_text(size=6),
      panel.grid=element_blank(),
      plot.margin=margin(3,3,3,3,"mm")
    )
}

# ============================================================================
make_spaghetti <- function(db) {
  tag <- tolower(db)
  db_label <- ifelse(tag=="eicu", "eICU-CRD", "MIMIC-IV")
  cat(sprintf("\n── Spaghetti: %s ──\n", db_label))

  # Load matched dataset
  matched <- read.csv(file.path(RESULTS, sprintf("did_matched_%s_24h.csv", tag)),
                       stringsAsFactors=FALSE)

  # Load treated ETL for mg_offset_min
  trt <- read.csv(file.path(RESULTS, sprintf("did_treated_%s.csv", tag)),
                   stringsAsFactors=FALSE)
  id_col <- if ("patientunitstayid" %in% names(trt)) "patientunitstayid" else "stay_id"
  trt$pid <- trt[[id_col]]
  mg_map <- setNames(trt$mg_offset_min, trt$pid)

  # Load all Cr
  cr_all <- read.csv(file.path(RESULTS, sprintf("did_cr_all_%s.csv", tag)),
                      stringsAsFactors=FALSE)
  cr_id <- if ("patientunitstayid" %in% names(cr_all)) "patientunitstayid" else "stay_id"
  cr_all$pid <- cr_all[[cr_id]]

  # ── Build treated spaghetti data ──────────────────────────────────────
  trt_pids <- matched$pid[matched$treated == 1]
  trt_pairs <- setNames(matched$match_pair_id[matched$treated == 1],
                         matched$pid[matched$treated == 1])
  trt_crpre <- setNames(matched$cr_pre[matched$treated == 1],
                          matched$pid[matched$treated == 1])

  cr_trt <- cr_all[cr_all$pid %in% trt_pids, ]
  cr_trt$mg_off <- mg_map[as.character(cr_trt$pid)]
  cr_trt$t_from_mg <- (cr_trt$labresultoffset - cr_trt$mg_off) / 60
  cr_trt$cr_pre <- trt_crpre[as.character(cr_trt$pid)]
  cr_trt$dcr <- cr_trt$labresult - cr_trt$cr_pre
  cr_trt$group <- "IV Mg"

  # ── Build control spaghetti data ─────────────────────────────────────
  ctl_pids <- matched$pid[matched$treated == 0]
  ctl_pair_ids <- matched$match_pair_id[matched$treated == 0]
  ctl_crpre <- setNames(matched$cr_pre[matched$treated == 0],
                          matched$pid[matched$treated == 0])

  # Map control → their matched treated partner's mg_offset
  pair_to_trt <- setNames(matched$pid[matched$treated == 1],
                            matched$match_pair_id[matched$treated == 1])
  ctl_mg_ref <- sapply(ctl_pair_ids, function(pid) {
    trt_pid <- pair_to_trt[as.character(pid)]
    if (!is.na(trt_pid)) mg_map[as.character(trt_pid)] else NA
  })
  names(ctl_mg_ref) <- ctl_pids

  cr_ctl <- cr_all[cr_all$pid %in% ctl_pids, ]
  cr_ctl$mg_off <- ctl_mg_ref[as.character(cr_ctl$pid)]
  cr_ctl$t_from_mg <- (cr_ctl$labresultoffset - cr_ctl$mg_off) / 60
  cr_ctl$cr_pre <- ctl_crpre[as.character(cr_ctl$pid)]
  cr_ctl$dcr <- cr_ctl$labresult - cr_ctl$cr_pre
  cr_ctl$group <- "Control"

  # ── Combine and filter time range ────────────────────────────────────
  both <- rbind(cr_trt[, c("pid","t_from_mg","dcr","group")],
                cr_ctl[, c("pid","t_from_mg","dcr","group")])
  both <- both[!is.na(both$t_from_mg) & !is.na(both$dcr), ]
  both <- both[both$t_from_mg >= -12 & both$t_from_mg <= 48, ]

  # Winsorize extreme dcr for display
  both$dcr <- pmin(pmax(both$dcr, -2), 3)

  cat(sprintf("  Treated lines: %d patients, %d Cr measurements\n",
              length(unique(both$pid[both$group=="IV Mg"])),
              sum(both$group=="IV Mg")))
  cat(sprintf("  Control lines: %d patients, %d Cr measurements\n",
              length(unique(both$pid[both$group=="Control"])),
              sum(both$group=="Control")))

  # ── Shared Y range ───────────────────────────────────────────────────
  ylim_range <- range(both$dcr, na.rm=TRUE)

  # ── Mean trend per group (binned) ────────────────────────────────────
  both$tbin <- round(both$t_from_mg / 3) * 3  # 3h bins
  mean_trend <- aggregate(dcr ~ group + tbin, data=both,
                           FUN=function(x) c(m=mean(x), se=sd(x)/sqrt(length(x))))
  mean_trend <- do.call(data.frame, mean_trend)
  names(mean_trend) <- c("group","tbin","mean","se")

  # ── Plot function for one panel ──────────────────────────────────────
  make_panel <- function(grp, panel_letter, color) {
    d <- both[both$group == grp, ]
    mt <- mean_trend[mean_trend$group == grp, ]

    # Subsample individuals for readability (max 500 lines)
    all_pids <- unique(d$pid)
    if (length(all_pids) > 500) {
      set.seed(42)
      show_pids <- sample(all_pids, 500)
      d_show <- d[d$pid %in% show_pids, ]
    } else {
      d_show <- d
    }

    n_pts <- length(unique(d$pid))

    ggplot() +
      geom_hline(yintercept=0, color="grey80", linewidth=0.3, linetype="dashed") +
      geom_vline(xintercept=0, color="grey60", linewidth=0.3, linetype="dotted") +
      # Individual lines
      geom_line(data=d_show, aes(x=t_from_mg, y=dcr, group=pid),
                color=color, alpha=0.06, linewidth=0.2) +
      # Mean trend
      geom_ribbon(data=mt, aes(x=tbin, ymin=mean-1.96*se, ymax=mean+1.96*se),
                  fill=color, alpha=0.25) +
      geom_line(data=mt, aes(x=tbin, y=mean), color=color, linewidth=1.2) +
      geom_point(data=mt, aes(x=tbin, y=mean), color=color, size=1.5) +
      scale_x_continuous(breaks=seq(-12, 48, 12),
                         labels=seq(-12, 48, 12)) +
      coord_cartesian(xlim=c(-12, 48), ylim=ylim_range) +
      labs(x="Hours from IV MgSO4 (t=0)",
           y="Cr change from baseline (mg/dL)",
           title=sprintf("(%s) %s: %s (n=%d)", panel_letter, db_label, grp, n_pts)) +
      annotate("text", x=0.5, y=ylim_range[2]*0.9, hjust=0, size=2,
               color="grey40", label="IV Mg given", fontface="italic") +
      theme_nature()
  }

  p1 <- make_panel("IV Mg", "a", WONG[["blue"]])
  p2 <- make_panel("Control", "b", WONG[["vermil"]])

  g <- arrangeGrob(p1, p2, ncol=2,
                    top=textGrob(sprintf("%s: Individual Cr Trajectories (matched pairs, aligned to IV Mg time)",
                                         db_label),
                                 gp=gpar(fontsize=8, fontface="bold")))

  out_base <- sprintf("fig_spaghetti_%s", tag)
  ggsave(file.path(RESULTS, paste0(out_base, ".pdf")), g,
         width=7.2, height=4, device=cairo_pdf)
  ggsave(file.path(RESULTS, paste0(out_base, ".png")), g,
         width=7.2, height=4, dpi=300)
  cat(sprintf("  Saved: %s.pdf / .png\n", out_base))
}

# ============================================================================
cat("======================================================================\n")
cat("fig_spaghetti_cr.R — Individual Cr trajectories (Dr. Su)\n")
cat("======================================================================\n")

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) args <- c("eicu","mimic")
for (a in args) make_spaghetti(a)
