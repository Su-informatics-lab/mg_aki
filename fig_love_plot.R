#!/usr/bin/env Rscript
# fig_love_plot.R — Covariate balance Love plot (Dr. Su)
# Shows SMD before (raw) and after (matched) for all 28 PS covariates
# Run: Rscript fig_love_plot.R eicu mimic

suppressPackageStartupMessages({
  library(MatchIt)
  library(ggplot2)
  library(grid)
})

RESULTS <- path.expand("~/mg_aki/results")
CALIPER <- 0.2

PS_COVARS <- c(
  "age", "is_female", "bmi",
  "surg_cabg", "surg_valve", "surg_combined",
  "heart_failure", "hypertension", "diabetes", "ckd",
  "copd", "pvd", "stroke", "liver_disease",
  "egfr",
  "loop_diuretics", "nsaids", "acei_arb", "ppi",
  "beta_blockers", "steroids", "antiarrhythmics",
  "first_potassium", "first_calcium", "first_heartrate",
  "first_mg_value", "first_lactate", "lactate_missing"
)

NICE_NAMES <- c(
  age="Age", is_female="Female sex", bmi="BMI",
  surg_cabg="CABG", surg_valve="Valve surgery", surg_combined="Combined surgery",
  heart_failure="Heart failure", hypertension="Hypertension", diabetes="Diabetes",
  ckd="CKD", copd="COPD", pvd="PVD", stroke="Stroke", liver_disease="Liver disease",
  egfr="eGFR", loop_diuretics="Loop diuretics", nsaids="NSAIDs", acei_arb="ACEi/ARB",
  ppi="PPI", beta_blockers="Beta blockers", steroids="Steroids",
  antiarrhythmics="Antiarrhythmics",
  first_potassium="Potassium", first_calcium="Calcium", first_heartrate="Heart rate",
  first_mg_value="Serum Mg", first_lactate="Lactate", lactate_missing="Lactate missing"
)

WONG <- c(blue="#0072B2", vermil="#D55E00", orange="#E69F00")

theme_nature <- function(base_size=7) {
  theme_classic(base_size=base_size, base_family="sans") %+replace%
    theme(axis.line=element_line(linewidth=0.3),
          axis.ticks=element_line(linewidth=0.3),
          axis.text=element_text(size=6, color="black"),
          axis.title=element_text(size=7),
          plot.title=element_text(size=8, face="bold", hjust=0),
          panel.grid=element_blank(),
          plot.margin=margin(3,5,3,3,"mm"))
}

median_impute <- function(d, vars) {
  for (v in vars)
    if (v %in% names(d) && any(is.na(d[[v]])))
      d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm=TRUE)
  d
}

make_love <- function(db) {
  tag <- tolower(db)
  db_label <- ifelse(tag=="eicu", "eICU-CRD", "MIMIC-IV")
  cat(sprintf("\n== %s Love Plot ==\n", db_label))

  trt <- read.csv(file.path(RESULTS, sprintf("did_treated_%s.csv", tag)), stringsAsFactors=F)
  ctl <- read.csv(file.path(RESULTS, sprintf("did_control_%s.csv", tag)), stringsAsFactors=F)
  id_col <- if ("patientunitstayid" %in% names(trt)) "patientunitstayid" else "stay_id"
  trt$pid <- trt[[id_col]]; ctl$pid <- ctl[[id_col]]

  ps_vars <- intersect(PS_COVARS, intersect(names(trt), names(ctl)))
  stack_cols <- intersect(unique(c("pid","treated",ps_vars)), intersect(names(trt),names(ctl)))
  combined <- rbind(trt[,stack_cols], ctl[,stack_cols])
  rownames(combined) <- seq_len(nrow(combined))
  combined <- median_impute(combined, ps_vars)

  # PS matching
  ps_fml <- as.formula(paste("treated ~", paste(ps_vars, collapse="+")))
  m <- suppressWarnings(matchit(ps_fml, data=combined, method="nearest",
                                 distance="glm", ratio=1, caliper=CALIPER, replace=TRUE))
  md <- match.data(m)

  # Compute SMDs
  smds <- data.frame(variable=ps_vars, stringsAsFactors=FALSE)
  smds$nice_name <- NICE_NAMES[smds$variable]
  smds$nice_name[is.na(smds$nice_name)] <- smds$variable[is.na(smds$nice_name)]

  for (i in seq_along(ps_vars)) {
    v <- ps_vars[i]
    # Raw
    x1 <- combined[[v]][combined$treated==1]; x0 <- combined[[v]][combined$treated==0]
    sp <- sqrt((var(x1,na.rm=T)+var(x0,na.rm=T))/2)
    smds$raw[i] <- if(!is.na(sp)&&sp>1e-10) abs(mean(x1,na.rm=T)-mean(x0,na.rm=T))/sp else 0
    # Matched
    x1m <- md[[v]][md$treated==1]; x0m <- md[[v]][md$treated==0]
    spm <- sqrt((var(x1m,na.rm=T)+var(x0m,na.rm=T))/2)
    smds$matched[i] <- if(!is.na(spm)&&spm>1e-10) abs(mean(x1m,na.rm=T)-mean(x0m,na.rm=T))/spm else 0
  }

  # Sort by raw SMD
  smds <- smds[order(smds$raw), ]
  smds$nice_name <- factor(smds$nice_name, levels=smds$nice_name)

  cat(sprintf("  Top 5 hardest to match:\n"))
  top5 <- tail(smds, 5)
  for (i in seq_len(nrow(top5)))
    cat(sprintf("    %s: raw=%.3f -> matched=%.3f\n",
                top5$nice_name[i], top5$raw[i], top5$matched[i]))

  # Save CSV
  write.csv(smds, file.path(RESULTS, sprintf("did_love_%s.csv", tag)), row.names=FALSE)

  # ── Love plot ──────────────────────────────────────────────────────
  p <- ggplot(smds) +
    geom_vline(xintercept=0.1, color="grey70", linewidth=0.3, linetype="dashed") +
    geom_segment(aes(x=matched, xend=raw, y=nice_name, yend=nice_name),
                 color="grey80", linewidth=0.3) +
    geom_point(aes(x=raw, y=nice_name), color=WONG[["vermil"]],
               shape=4, size=1.5, stroke=0.5) +
    geom_point(aes(x=matched, y=nice_name), color=WONG[["blue"]],
               shape=16, size=2) +
    scale_x_continuous(limits=c(0, max(smds$raw)*1.1),
                       breaks=seq(0, 2, 0.2)) +
    labs(x="Absolute Standardized Mean Difference",
         y=NULL,
         title=sprintf("%s: Covariate Balance", db_label)) +
    annotate("text", x=max(smds$raw)*0.7, y=3, size=2.2,
             label="X = Before matching", color=WONG[["vermil"]]) +
    annotate("text", x=max(smds$raw)*0.7, y=1.5, size=2.2,
             label="O = After matching", color=WONG[["blue"]]) +
    theme_nature()

  ggsave(file.path(RESULTS, sprintf("fig_love_%s.pdf", tag)), p,
         width=4.7, height=5, device=cairo_pdf)
  ggsave(file.path(RESULTS, sprintf("fig_love_%s.png", tag)), p,
         width=4.7, height=5, dpi=300)
  cat(sprintf("  Saved: fig_love_%s.pdf / .png\n", tag))
}

cat("======================================================================\n")
cat("fig_love_plot.R — Covariate Balance Love Plot\n")
cat("======================================================================\n")

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) args <- c("eicu","mimic")
for (a in args) make_love(a)
