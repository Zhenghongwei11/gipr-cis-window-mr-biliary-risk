#!/usr/bin/env Rscript
# Figure 1: Study design and methodological framework (workflow schematic)

source("scripts/mr/lib_io.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(ggforce) 
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  args[idx + 1]
}
plot_outdir <- get_arg("--outdir", "plots/causal/drug_target_mr")
strict_primary_path <- get_arg(
  "--strict-primary",
  "docs/source_data/Source_Data_3_primary_results_hba1c_strictGIPR200kb.tsv"
)

fmt_or_p <- function(or, lo, hi, p) {
  or_txt <- sprintf("OR %.2f [%.2f\u2013%.2f]", or, lo, hi)
  p_txt <- if (is.na(p)) "P = NA" else sprintf("P = %.3f", p)
  paste0(or_txt, "\n", p_txt)
}

or_callout <- "OR (see Table 2)"
if (file.exists(strict_primary_path)) {
  df <- read_tsv_any(strict_primary_path)
  df$pvalue <- as.numeric(df$pvalue)
  df$or <- as.numeric(df$or)
  df$or_ci_lower <- as.numeric(df$or_ci_lower)
  df$or_ci_upper <- as.numeric(df$or_ci_upper)
  hit <- df[df$target_symbol == "GIPR" & df$outcome_name == "cholelithiasis_primary", , drop = FALSE]
  if (nrow(hit) >= 1) {
    or_callout <- fmt_or_p(hit$or[1], hit$or_ci_lower[1], hit$or_ci_upper[1], hit$pvalue[1])
  }
}

# --- Journal Standard Palette ---
pal <- list(
  navy    = "#34495E", 
  blue    = "#2980B9", 
  orange  = "#D35400", 
  green   = "#27AE60", 
  slate   = "#7F8C8D", 
  bg      = "#FFFFFF",
  border  = "#2C3E50"
)

p <- ggplot() +
  # --- PANEL A: GENOMIC DATA INFRASTRUCTURE ---
  annotate("text", x = -7.5, y = 10.5, label = "A", size = 6, fontface = "bold") +
  annotate("text", x = -7.2, y = 10.5, label = "  Genomic Data Infrastructure", size = 3.5, fontface = "bold", hjust = 0, colour = pal$blue) +
  geom_ellipse(aes(x0 = -5, y0 = 9.4, a = 1.0, b = 0.2, angle = 0), fill = pal$blue, alpha = 0.4) +
  geom_ellipse(aes(x0 = -5, y0 = 9.1, a = 1.0, b = 0.2, angle = 0), fill = pal$blue, alpha = 0.7) +
  annotate("text", x = -5, y = 8.5, label = "HbA1c GWAS\n(UK Biobank-derived)", size = 2.4) +
  geom_ellipse(aes(x0 = -2, y0 = 9.4, a = 1.0, b = 0.2, angle = 0), fill = pal$orange, alpha = 0.4) +
  geom_ellipse(aes(x0 = -2, y0 = 9.1, a = 1.0, b = 0.2, angle = 0), fill = pal$orange, alpha = 0.7) +
  annotate("text", x = -2, y = 8.5, label = "Primary biliary outcomes\n+ exploratory pancreatitis", size = 2.4) +

  # --- PANEL B: PHENOTYPIC AND CROSS-COHORT TRIANGULATION ---
  annotate("text", x = 0.5, y = 10.5, label = "B", size = 6, fontface = "bold") +
  annotate("text", x = 0.8, y = 10.5, label = "  Phenotypic & Cross-Cohort Triangulation", size = 3.5, fontface = "bold", hjust = 0, colour = pal$orange) +
  annotate("rect", xmin = 1.5, xmax = 3.7, ymin = 9.2, ymax = 9.9, fill = pal$orange, alpha = 0.12, colour = pal$orange, linewidth = 0.4) +
  annotate("text", x = 2.6, y = 9.55, label = "FinnGen R12\ncholelithiasis", size = 2.3) +
  annotate("rect", xmin = 4.3, xmax = 6.7, ymin = 9.2, ymax = 9.9, fill = pal$blue, alpha = 0.10, colour = pal$blue, linewidth = 0.4) +
  annotate("text", x = 5.5, y = 9.55, label = "UK Biobank\ncholecystectomy", size = 2.3) +
  geom_segment(aes(x = 3.8, y = 9.55, xend = 4.2, yend = 9.55), arrow = arrow(length = unit(0.06, "inches"), type = "closed"), colour = pal$slate, linewidth = 0.3) +
  annotate("text", x = 4.1, y = 8.65, label = "Diagnosis-procedure consistency\nwith potential sample overlap noted", size = 2.15, fontface = "italic", colour = pal$slate) +

  # --- PANEL C: REGIONAL PLEIOTROPY SENSITIVITY ANALYSES ---
  annotate("text", x = -7.5, y = 5.5, label = "C", size = 6, fontface = "bold") +
  annotate("text", x = -7.2, y = 5.5, label = "  Regional Pleiotropy Sensitivity Analyses", size = 3.5, fontface = "bold", hjust = 0, colour = pal$green) +
  geom_polygon(aes(x = c(-5.5, -1.5, -2.5, -4.5), y = c(4.8, 4.8, 3.2, 3.2)), fill = pal$green, alpha = 0.1, color = pal$green, linewidth = 0.4) +
  annotate("text", x = -3.5, y = 4.4, label = "Broad cis-window (\u00b11 Mb)", size = 2.4) +
  annotate("text", x = -3.5, y = 4.0, label = "Regional pleiotropy contamination", size = 2.1, colour = pal$orange, fontface = "italic") +
  annotate("text", x = -3.5, y = 3.5, label = "Strict cis-window (\u00b1200 kb; GIPR)", size = 2.4, fontface = "bold") +
  geom_point(aes(x = -3.5, y = 2.8), shape = 8, size = 2.5, colour = pal$green) +
  annotate("text", x = -3.5, y = 2.4, label = "Instrument PheWAS\n(Chr19q13 characterisation)", size = 2.2) +

  # --- PANEL D: SYNTHESISED LOCUS-PROXY ASSOCIATIONS ---
  annotate("text", x = 0.5, y = 5.5, label = "D", size = 6, fontface = "bold") +
  annotate("text", x = 0.8, y = 5.5, label = "  Synthesised Locus-Proxy Associations", size = 3.5, fontface = "bold", hjust = 0, colour = pal$navy) +
  geom_segment(aes(x = 2, y = 3.5, xend = 6, yend = 3.5), colour = pal$slate, linewidth = 0.3) + 
  geom_segment(aes(x = 3.5, y = 3.0, xend = 3.5, yend = 3.9), linetype = "dotted", colour = pal$slate, linewidth = 0.3) +
  geom_point(aes(x = 5.2, y = 3.5), shape = 15, size = 4, colour = pal$orange) + 
  geom_segment(aes(x = 4.2, y = 3.5, xend = 6.2, yend = 3.5), colour = pal$orange, linewidth = 1.2) +
  annotate("text", x = 4, y = 4.2, label = "GIPR locus proxy \u2192 Cholelithiasis (FinnGen R12)", size = 2.7, fontface = "bold", colour = pal$navy) +
  annotate("text", x = 4, y = 2.8, label = or_callout, size = 2.6, fontface = "bold") +

  # --- ACADEMIC CONNECTORS ---
  geom_curve(aes(x = -3.5, y = 8.6, xend = -3.5, yend = 5.0), curvature = -0.15, arrow = arrow(length = unit(0.08, "inches"), type = "closed"), colour = pal$slate) +
  geom_segment(aes(x = -1.2, y = 3.5, xend = 1.4, yend = 3.5), arrow = arrow(length = unit(0.08, "inches"), type = "closed"), colour = pal$slate) +

  # Final Layout Settings
  coord_cartesian(xlim = c(-8, 8), ylim = c(1.5, 12), expand = FALSE) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = pal$bg, colour = NA),
    plot.margin = margin(10, 10, 10, 10, "pt")
  )

# --- Save Publication-Ready Assets ---
dir.create(plot_outdir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(plot_outdir, "workflow_diagram.pdf"), p, width = 8, height = 7, device = cairo_pdf)
ggsave(file.path(plot_outdir, "workflow_diagram.png"), p, width = 8, height = 7, dpi = 350)

# Final Sync to Publication Folder
file.copy(file.path(plot_outdir, "workflow_diagram.pdf"), "plots/publication/pdf/Figure1.pdf", overwrite = TRUE)
file.copy(file.path(plot_outdir, "workflow_diagram.png"), "plots/publication/png/Figure1.png", overwrite = TRUE)

message("Final Publication-Ready Framework Diagram (Shortened Title, No Figure #) generated.")
