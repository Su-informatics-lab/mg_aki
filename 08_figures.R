#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────
# 08_figures.R — Publication figures (JNO style)
#   Figure 1: CONSORT flow (both databases)
#   Figure 2: Forest plot (pooled AKI + mortality + controls)
#   Figure 3: Prognostic Mg-AKI × surgery type (cardioplegia hypothesis)
# ─────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyverse); library(grid); library(gridExtra)
})

RESULTS <- path.expand("~/mg_aki/results")

# ── JNO theme ───────────────────────────────────────────────────────
theme_jno <- function(base_size = 9) {
  theme_minimal(base_size = base_size, base_family = "Helvetica") +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.3),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.margin = margin(8, 12, 8, 8),
      legend.position = "bottom",
      legend.key.size = unit(0.4, "cm")
    )
}

COL_EICU  <- "#2166AC"  # blue
COL_MIMIC <- "#B2182B"  # red
COL_POOL  <- "black"


# =====================================================================
# FIGURE 2: Forest Plot (primary — this is the money figure)
# =====================================================================
cat("FIGURE 2: Forest plot\n")

meta <- read_csv(file.path(RESULTS, "06_meta_results.csv"), show_col_types = FALSE)
eicu <- read_csv(file.path(RESULTS, "03_results_summary.csv"), show_col_types = FALSE)
mimic <- read_csv(file.path(RESULTS, "05_mimic_results_summary.csv"), show_col_types = FALSE)

# Build forest data: each outcome gets 3 rows (eICU, MIMIC, Pooled)
sections <- list(
  list(section = "AKI Outcomes",
       outcomes = list(
         list("KDIGO stage ≥1 (primary)", "ALL_aki_KDIGO >=1", "KDIGO Stage >=1"),
         list("Creatinine ratio ≥1.5×", "ALL_aki_Ratio >=1.5x", "Ratio >=1.5x"),
         list("KDIGO stage ≥2", "ALL_aki_Stage >=2", "Stage >=2"),
         list("AKI within 48 h", "ALL_tw_AKI 1.5x <=48h", "AKI 1.5x <=48h")
       )),
  list(section = "Mortality",
       outcomes = list(
         list("Hospital mortality (all+all)", "ALL_sec_Hospital mortality", "Mortality (all+all)"),
         list("Hospital mortality (hypo+all)", "HYPO_sec_Hospital mortality", "Mortality (hypo+all)")
       )),
  list(section = "Control Outcomes",
       outcomes = list(
         list("Fracture (negative control)", "ALL_nc_Fracture", "Fracture (neg ctrl)"),
         list("Encephalopathy (exploratory)", "ALL_neuro_Encephalopathy", "Encephalopathy")
       ))
)

forest_rows <- list()
row_idx <- 0

for (sec in rev(sections)) {
  for (oc in rev(sec$outcomes)) {
    label <- oc[[1]]; model_key <- oc[[2]]; meta_key <- oc[[3]]

    # eICU — use HYPO for hypo+all mortality, else ALL
    e_key <- model_key
    if (grepl("HYPO_", model_key)) e_key <- model_key
    e <- eicu %>% filter(model == e_key)
    m <- mimic %>% filter(model == ifelse(grepl("HYPO_", model_key), "ALL_sec_Hospital mortality", model_key))
    p <- meta %>% filter(outcome == meta_key)

    if (nrow(p) == 0) next

    # Pooled
    row_idx <- row_idx + 1
    forest_rows[[length(forest_rows)+1]] <- tibble(
      y = row_idx, label = "  Pooled", outcome = label,
      or = p$pooled_or, lo = p$pooled_lo, hi = p$pooled_hi,
      source = "Pooled", section = sec$section,
      p_val = p$pooled_p, i2 = p$I2)

    # MIMIC
    if (nrow(m) > 0) {
      row_idx <- row_idx + 1
      forest_rows[[length(forest_rows)+1]] <- tibble(
        y = row_idx, label = "  MIMIC-IV", outcome = label,
        or = m$estimate, lo = m$conf.low, hi = m$conf.high,
        source = "MIMIC-IV", section = sec$section,
        p_val = m$p.value, i2 = NA)
    }

    # eICU
    if (nrow(e) > 0) {
      row_idx <- row_idx + 1
      forest_rows[[length(forest_rows)+1]] <- tibble(
        y = row_idx, label = "  eICU", outcome = label,
        or = e$estimate, lo = e$conf.low, hi = e$conf.high,
        source = "eICU", section = sec$section,
        p_val = e$p.value, i2 = NA)
    }

    # Outcome header
    row_idx <- row_idx + 1
    forest_rows[[length(forest_rows)+1]] <- tibble(
      y = row_idx, label = label, outcome = label,
      or = NA, lo = NA, hi = NA,
      source = "header", section = sec$section,
      p_val = NA, i2 = NA)
  }

  # Section header
  row_idx <- row_idx + 1
  forest_rows[[length(forest_rows)+1]] <- tibble(
    y = row_idx, label = sec$section, outcome = NA,
    or = NA, lo = NA, hi = NA,
    source = "section", section = sec$section,
    p_val = NA, i2 = NA)
}

fd <- bind_rows(forest_rows)

# Format OR text
fd$or_text <- ifelse(is.na(fd$or), "",
  sprintf("%.2f (%.2f–%.2f)", fd$or, fd$lo, fd$hi))
fd$p_text <- ifelse(is.na(fd$p_val), "",
  ifelse(fd$p_val < 0.001, "P<.001",
  sprintf("P=%.3f", fd$p_val)))

# Point shapes and colors
fd$pt_color <- case_when(
  fd$source == "eICU" ~ COL_EICU,
  fd$source == "MIMIC-IV" ~ COL_MIMIC,
  fd$source == "Pooled" ~ COL_POOL,
  TRUE ~ NA_character_)
fd$pt_shape <- case_when(
  fd$source == "Pooled" ~ 18,  # diamond
  fd$source %in% c("eICU", "MIMIC-IV") ~ 15,  # square
  TRUE ~ NA_integer_)
fd$pt_size <- ifelse(fd$source == "Pooled", 4, 2.5)
fd$face <- case_when(
  fd$source == "section" ~ "bold",
  fd$source == "header" ~ "bold",
  fd$source == "Pooled" ~ "bold",
  TRUE ~ "plain")

points_data <- fd %>% filter(!is.na(or))

p_forest <- ggplot(fd, aes(y = y)) +
  # Reference line at 1
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  # CI lines
  geom_segment(data = points_data,
               aes(x = lo, xend = hi, yend = y, color = source),
               linewidth = 0.5, show.legend = FALSE) +
  # Points
  geom_point(data = points_data,
             aes(x = or, shape = source, color = source, size = source),
             show.legend = FALSE) +
  # Labels on left
  geom_text(aes(x = 0.28, label = label, fontface = face),
            hjust = 0, size = 2.8, family = "Helvetica") +
  # OR text on right
  geom_text(aes(x = 2.6, label = or_text), hjust = 0, size = 2.4,
            family = "Helvetica") +
  # Scale
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2),
                limits = c(0.25, 3.5)) +
  scale_color_manual(values = c("eICU" = COL_EICU, "MIMIC-IV" = COL_MIMIC,
                                "Pooled" = COL_POOL)) +
  scale_shape_manual(values = c("eICU" = 15, "MIMIC-IV" = 15, "Pooled" = 18)) +
  scale_size_manual(values = c("eICU" = 2.5, "MIMIC-IV" = 2.5, "Pooled" = 4)) +
  labs(x = "Odds Ratio (95% CI)", y = NULL,
       title = "Figure 2. Pooled Odds Ratios for AKI, Mortality, and Control Outcomes") +
  annotate("text", x = 0.5, y = 0, label = "Favors supplementation",
           size = 2.2, hjust = 0.5, color = "grey40") +
  annotate("text", x = 1.5, y = 0, label = "Favors no supplementation",
           size = 2.2, hjust = 0.5, color = "grey40") +
  theme_jno(base_size = 9) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size = 10))

ggsave(file.path(RESULTS, "08_fig2_forest.pdf"), p_forest,
       width = 8.5, height = 7, dpi = 300)
cat("  Saved: 08_fig2_forest.pdf\n")


# =====================================================================
# FIGURE 1: CONSORT Flow
# =====================================================================
cat("\nFIGURE 1: CONSORT flow\n")

# eICU numbers from 00_consort.csv
consort <- read_csv(file.path(RESULTS, "00_consort.csv"), show_col_types = FALSE)
get_n <- function(step) consort$n[consort$step == step]

# MIMIC numbers (from cohort file)
mimic_cohort <- read_csv(file.path(RESULTS, "04_mimic_cohort.csv"), show_col_types = FALSE)
n_mimic <- nrow(mimic_cohort)
n_mimic_trt <- sum(mimic_cohort$mg_supplementation)

# Simple text-based CONSORT as a table figure
consort_text <- tribble(
  ~Step, ~eICU, ~`MIMIC-IV`,
  "Total ICU admissions", "200,859", "94,458",
  "Adults, cardiac surgery, first stay", "26,715", "13,706",
  "Mg within 6 h of admission", "9,379", "3,913",
  "Baseline creatinine available", "8,650", "3,746",
  "Excluded: Cr ≥4.0 or ESKD", "−541", "—",
  "Eligible cohort", "8,109", "3,746",
  "All-patient TTE (complete covariates)", "7,924", "3,746",
  "  Supplemented", "1,104 (13.9%)", "647 (17.3%)",
  "  Not supplemented", "6,820", "3,099",
  "HypoMg TTE (Mg <2.0, complete)", "3,068", "1,780",
  "  Supplemented", "759 (24.7%)", "509 (28.6%)",
  "  Not supplemented", "2,309", "1,271"
)

out_consort <- file.path(RESULTS, "08_fig1_consort.csv")
write_csv(consort_text, out_consort)
cat(sprintf("  Saved: %s (format as flow diagram in Illustrator/PowerPoint)\n", out_consort))

# Also save as a simple grid table figure
p_consort <- gridExtra::tableGrob(consort_text,
  rows = NULL,
  theme = ttheme_minimal(
    base_size = 9,
    base_family = "Helvetica",
    core = list(fg_params = list(hjust = 0, x = 0.05)),
    colhead = list(fg_params = list(fontface = "bold", hjust = 0, x = 0.05))
  ))

pdf(file.path(RESULTS, "08_fig1_consort.pdf"), width = 7, height = 4)
grid.draw(p_consort)
dev.off()
cat("  Saved: 08_fig1_consort.pdf\n")


# =====================================================================
# FIGURE 3: Surgery-Type Interaction (Cardioplegia Hypothesis)
# =====================================================================
cat("\nFIGURE 3: Surgery-type interaction\n")

# Hardcoded from run log (prognostic ORs — these are from the
# multivariable regression, not the TTE)
interaction_data <- tribble(
  ~database, ~surgery, ~or, ~lo, ~hi,
  "eICU",    "Simple\n(CABG/other)",   1.36, 1.19, 1.55,
  "eICU",    "Complex\n(valve/combined)", 1.74, 1.45, 2.11,
  "MIMIC-IV","Simple\n(CABG/other)",   1.01, 0.85, 1.19,
  "MIMIC-IV","Complex\n(valve/combined)", 1.09, 0.85, 1.41
)

interaction_data$database <- factor(interaction_data$database,
                                    levels = c("eICU", "MIMIC-IV"))
interaction_data$surgery <- factor(interaction_data$surgery,
                                   levels = c("Simple\n(CABG/other)",
                                              "Complex\n(valve/combined)"))

p_interaction <- ggplot(interaction_data,
       aes(x = surgery, y = or, color = database, shape = database)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_pointrange(aes(ymin = lo, ymax = hi),
                  position = position_dodge(width = 0.4),
                  size = 0.6, linewidth = 0.5) +
  scale_color_manual(values = c("eICU" = COL_EICU, "MIMIC-IV" = COL_MIMIC),
                     name = NULL) +
  scale_shape_manual(values = c("eICU" = 15, "MIMIC-IV" = 17), name = NULL) +
  scale_y_continuous(breaks = seq(0.8, 2.2, 0.2)) +
  coord_cartesian(ylim = c(0.7, 2.3)) +
  labs(x = NULL, y = "OR per 1 mg/dL Serum Mg Increase\n(Prognostic Association With AKI)",
       title = "Figure 3. Prognostic Magnesium–AKI Association by Surgery Type") +
  annotate("text", x = 2.3, y = 1.74, label = "More cardioplegia →\nhigher Mg + more AKI",
           size = 2.5, hjust = 0, color = "grey40", fontface = "italic") +
  theme_jno(base_size = 10) +
  theme(legend.position = c(0.15, 0.9),
        legend.background = element_rect(fill = "white", color = NA))

ggsave(file.path(RESULTS, "08_fig3_interaction.pdf"), p_interaction,
       width = 5, height = 5, dpi = 300)
cat("  Saved: 08_fig3_interaction.pdf\n")


# =====================================================================
# eFigure: PS Overlap Density (both databases)
# =====================================================================
cat("\neFigure: PS overlap density\n")

tryCatch({
  dat_all <- read_csv(file.path(RESULTS, "02e_all_iptw.csv"), show_col_types = FALSE)
  p_ps_eicu <- ggplot(dat_all, aes(x = ps, fill = factor(trt))) +
    geom_density(alpha = 0.5, color = NA) +
    scale_fill_manual(values = c("0" = "grey60", "1" = COL_EICU),
                      labels = c("No supplementation", "Mg supplementation"),
                      name = NULL) +
    labs(x = "Propensity Score", y = "Density",
         title = "eICU: Propensity Score Overlap (All Patients)") +
    theme_jno() + theme(legend.position = c(0.7, 0.9))

  ggsave(file.path(RESULTS, "08_efig_ps_eicu.pdf"), p_ps_eicu,
         width = 5, height = 3.5, dpi = 300)
  cat("  Saved: 08_efig_ps_eicu.pdf\n")
}, error = function(e) cat(sprintf("  PS plot failed: %s\n", e$message)))


cat("\n08_figures.R COMPLETE\n")
