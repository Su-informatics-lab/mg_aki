#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────
# 06_figures.R — Publication figures
#   Figure 2: Forest plot (AC primary + sensitivity + controls)
#   Figure 3: Prognostic Mg-AKI × surgery type
#   eFigure: PS overlap density
#
# Reads: results/02_results.csv
# ─────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({ library(ggplot2); library(grid); library(gridExtra) })
RESULTS <- path.expand("~/mg_aki/results")

theme_jno <- function(base_size = 9) {
  theme_minimal(base_size = base_size, base_family = "Helvetica") +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(color="black", linewidth=0.3),
          axis.ticks = element_line(color="black", linewidth=0.3),
          plot.title = element_text(size=base_size+1, face="bold", hjust=0),
          plot.margin = margin(8,12,8,8), legend.position = "bottom")
}

COL_EICU <- "#2166AC"; COL_MIMIC <- "#B2182B"; COL_POOL <- "black"

# =====================================================================
# FIGURE 2: Forest Plot
# =====================================================================
cat("FIGURE 2: Forest plot\n")

res <- read.csv(file.path(RESULTS, "02_results.csv"), stringsAsFactors=FALSE)
res <- res[is.na(res$m) | res$m == min(res$m, na.rm=TRUE), ]

# Define sections (bottom to top in the plot)
sections <- list(
  list(section="Control", outcomes=list(
    list(label="Fracture (negative control)", a="ow_frac"),
    list(label="Encephalopathy (exploratory)", a="ow_enceph")
  )),
  list(section="Exploratory", outcomes=list(
    list(label="Hospital mortality", a="ow_mort")
  )),
  list(section="Sensitivity: All-Patient", outcomes=list(
    list(label="PS matching", a="psm_aki1"),
    list(label="Overlap weighting", a="ow_aki1"),
    list(label="IPTW", a="iptw_aki1")
  )),
  list(section="Primary: Active Comparator", outcomes=list(
    list(label="KDIGO stage ≥1 (Mg+K vs K-only)", a="ac_aki1")
  ))
)

rows <- list(); y <- 0
for (sec in sections) {
  for (oc in sec$outcomes) {
    e <- res[res$db=="eICU" & res$analysis==oc$a,]
    m <- res[res$db=="MIMIC" & res$analysis==oc$a,]
    p <- res[res$db=="Pooled" & res$analysis==oc$a,]
    if (nrow(p)==0) next
    # Pooled
    y <- y+1; rows[[length(rows)+1]] <- data.frame(
      y=y, label="  Pooled", or=p$or, lo=p$lo, hi=p$hi, source="Pooled", section=sec$section)
    # MIMIC
    if (nrow(m)>0) { y<-y+1; rows[[length(rows)+1]] <- data.frame(
      y=y, label="  MIMIC-IV", or=m$or, lo=m$lo, hi=m$hi, source="MIMIC-IV", section=sec$section) }
    # eICU
    if (nrow(e)>0) { y<-y+1; rows[[length(rows)+1]] <- data.frame(
      y=y, label="  eICU", or=e$or, lo=e$lo, hi=e$hi, source="eICU", section=sec$section) }
    # Header
    y<-y+1; rows[[length(rows)+1]] <- data.frame(
      y=y, label=oc$label, or=NA, lo=NA, hi=NA, source="header", section=sec$section)
  }
  # Section header
  y<-y+1; rows[[length(rows)+1]] <- data.frame(
    y=y, label=sec$section, or=NA, lo=NA, hi=NA, source="section", section=sec$section)
}

fd <- do.call(rbind, rows)
fd$or_text <- ifelse(is.na(fd$or), "", sprintf("%.2f (%.2f\u2013%.2f)", fd$or, fd$lo, fd$hi))
fd$face <- ifelse(fd$source %in% c("section","header","Pooled"), "bold", "plain")
pts <- fd[!is.na(fd$or),]

p_forest <- ggplot(fd, aes(y=y)) +
  geom_vline(xintercept=1, linetype="dashed", color="grey50", linewidth=0.3) +
  geom_segment(data=pts, aes(x=lo, xend=hi, yend=y, color=source), linewidth=0.5, show.legend=FALSE) +
  geom_point(data=pts, aes(x=or, shape=source, color=source, size=source), show.legend=FALSE) +
  geom_text(aes(x=0.28, label=label, fontface=face), hjust=0, size=2.8) +
  geom_text(aes(x=2.6, label=or_text), hjust=0, size=2.4) +
  scale_x_log10(breaks=c(0.5,0.75,1,1.5,2), limits=c(0.25,3.5)) +
  scale_color_manual(values=c("eICU"=COL_EICU, "MIMIC-IV"=COL_MIMIC, "Pooled"=COL_POOL)) +
  scale_shape_manual(values=c("eICU"=15, "MIMIC-IV"=15, "Pooled"=18)) +
  scale_size_manual(values=c("eICU"=2.5, "MIMIC-IV"=2.5, "Pooled"=4)) +
  labs(x="Odds Ratio (95% CI)", y=NULL,
       title="Figure 2. Association of Magnesium Supplementation With AKI and Secondary Outcomes") +
  annotate("text", x=0.5, y=0, label="Favors supplementation", size=2.2, color="grey40") +
  annotate("text", x=1.5, y=0, label="Favors no supplementation", size=2.2, color="grey40") +
  theme_jno(9) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
                        panel.grid.major.y=element_blank(), plot.title=element_text(size=10))

ggsave(file.path(RESULTS, "06_fig2_forest.pdf"), p_forest, width=8.5, height=6, dpi=300)
cat("  Saved: 06_fig2_forest.pdf\n")

# =====================================================================
# FIGURE 3: Surgery-Type Interaction (hardcoded prognostic ORs)
# =====================================================================
cat("\nFIGURE 3: Surgery-type interaction\n")

int_data <- data.frame(
  database = factor(rep(c("eICU","MIMIC-IV"), each=2), levels=c("eICU","MIMIC-IV")),
  surgery = factor(rep(c("Simple\n(CABG/other)","Complex\n(valve/combined)"), 2),
                   levels=c("Simple\n(CABG/other)","Complex\n(valve/combined)")),
  or = c(1.36, 1.74, 1.01, 1.09),
  lo = c(1.19, 1.45, 0.85, 0.85),
  hi = c(1.55, 2.11, 1.19, 1.41)
)

p_int <- ggplot(int_data, aes(x=surgery, y=or, color=database, shape=database)) +
  geom_hline(yintercept=1, linetype="dashed", color="grey50", linewidth=0.3) +
  geom_pointrange(aes(ymin=lo, ymax=hi), position=position_dodge(0.4), size=0.6, linewidth=0.5) +
  scale_color_manual(values=c("eICU"=COL_EICU, "MIMIC-IV"=COL_MIMIC), name=NULL) +
  scale_shape_manual(values=c("eICU"=15, "MIMIC-IV"=17), name=NULL) +
  scale_y_continuous(breaks=seq(0.8,2.2,0.2)) + coord_cartesian(ylim=c(0.7,2.3)) +
  labs(x=NULL, y="OR per 1 mg/dL Serum Mg Increase\n(Prognostic Association With AKI)",
       title="Figure 3. Prognostic Magnesium\u2013AKI Association by Surgery Type") +
  theme_jno(10) + theme(legend.position=c(0.15,0.9),
                         legend.background=element_rect(fill="white",color=NA))

ggsave(file.path(RESULTS, "06_fig3_interaction.pdf"), p_int, width=5, height=5, dpi=300)
cat("  Saved: 06_fig3_interaction.pdf\n")

# =====================================================================
# FIGURE 1: CONSORT (CSV for manual flow diagram)
# =====================================================================
cat("\nFIGURE 1: CONSORT flow (CSV)\n")

consort <- data.frame(
  Step = c("Total ICU admissions", "Adults, cardiac surgery, first stay",
           "Mg within 6h", "Baseline Cr available", "Excluded: Cr >=4.0 or ESKD",
           "Eligible cohort", "  Supplemented", "  Not supplemented",
           "Active comparator (K-repleted)", "  Mg+K", "  K-only"),
  eICU = c("200,859", "26,715", "9,379", "8,650", "-541",
           "8,109", "1,128 (13.9%)", "6,981",
           "1,986", "557", "1,429"),
  MIMIC = c("94,458", "13,706", "3,913", "3,746", "—",
            "3,746", "647 (17.3%)", "3,099",
            "844", "280", "564"),
  stringsAsFactors = FALSE
)
write.csv(consort, file.path(RESULTS, "06_fig1_consort.csv"), row.names=FALSE)
cat("  Saved: 06_fig1_consort.csv\n")

cat("\n06_figures.R COMPLETE\n")
