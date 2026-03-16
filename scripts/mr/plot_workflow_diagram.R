#!/usr/bin/env Rscript
# Figure 1: Scientific Workflow Schematic (Final Refined Version)
#
# Logic Alignment:
# Step 1: Instrument Discovery
# Step 2: Primary MR Analysis
# Step 3: Methodological Guardrails
# Step 4: Evidence Triangulation

source("scripts/mr/lib_io.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) stop("Missing value for ", flag)
  args[idx + 1]
}

plot_outdir <- get_arg("--outdir", "plots/causal/drug_target_mr")
pheno_defs <- get_arg("--pheno-defs", "docs/PHENOTYPE_DEFS.tsv")

# --- Aesthetic Configuration ---
cols <- list(
  header  = "#2C3E50", # Dark Navy
  step1   = "#34495E", # Slate
  step2   = "#2980B9", # Muted Blue
  step3   = "#E67E22", # Muted Orange (Sensitivity)
  step3b  = "#27AE60", # Muted Green (Pleiotropy)
  step4   = "#2C3E50", # Dark Navy (Synthesis)
  border  = "#2C3E50",
  text    = "#2C3E50"
)

draw_box <- function(x, y, w, h, fill, label, subtitle = "") {
  data.frame(x=x, y=y, w=w, h=h, fill=fill, label=label, subtitle=subtitle, stringsAsFactors=FALSE)
}

fmt_n_total <- function(n_case, n_control) {
  suppressWarnings({
    nc <- as.numeric(n_case)
    nn <- as.numeric(n_control)
  })
  if (is.na(nc) || is.na(nn)) return(NA_character_)
  format(nc + nn, big.mark = ",", scientific = FALSE)
}

phenos <- NULL
if (file.exists(pheno_defs)) {
  phenos <- read_tsv_any(pheno_defs)
}

get_total_n <- function(outcome_name) {
  if (is.null(phenos) || !all(c("outcome_name", "n_case", "n_control") %in% names(phenos))) return(NA_character_)
  hit <- phenos[phenos$outcome_name == outcome_name, , drop = FALSE]
  if (nrow(hit) < 1) return(NA_character_)
  fmt_n_total(hit$n_case[1], hit$n_control[1])
}

# --- Define Schematic Components ---

# 1. Title / Header (Added per audit)
header <- draw_box(0, 11.2, 15, 0.8, "white", "Study design and methodological framework", "")

# Step 1: Instruments
s1 <- draw_box(0, 9.8, 14, 1.2, cols$step1, "STEP 1: GENETIC INSTRUMENT DISCOVERY", 
               "Exposure: HbA1c GWAS (UK Biobank-derived)")

# Step 2: Primary MR
N_finngen_cholelith <- get_total_n("cholelithiasis_primary")
N_finngen_cholecyst <- get_total_n("cholecystectomy_secondary")
N_ukb_cholecyst <- get_total_n("cholecystectomy_ukb")

finngen_primary_line <- if (!is.na(N_finngen_cholelith)) paste0("Primary: Cholelithiasis (N=", N_finngen_cholelith, ")") else "Primary: Cholelithiasis"
finngen_secondary_line <- if (!is.na(N_finngen_cholecyst)) paste0("Secondary: Cholecystectomy (N=", N_finngen_cholecyst, ")") else "Secondary: Cholecystectomy"
ukb_supporting_line <- if (!is.na(N_ukb_cholecyst)) paste0("Supporting: Cholecystectomy (Procedure; N=", N_ukb_cholecyst, ")") else "Supporting: Cholecystectomy (Procedure)"

s2_a <- draw_box(-3.8, 7.2, 6.4, 2.4, cols$step2, "STEP 2A: INDEPENDENT DISCOVERY", 
                 paste0("FinnGen Release 12\n", finngen_primary_line, "\n", finngen_secondary_line))
s2_b <- draw_box(3.8, 7.2, 6.4, 2.4, cols$step2, "STEP 2B: CLINICAL TRIANGULATION", 
                 paste0("UK Biobank (UKB)\n", ukb_supporting_line, "\n(Transparency: potential sample overlap)"))

# Step 3: Guardrails
s3_a <- draw_box(-3.8, 4.0, 6.4, 2.2, cols$step3, "STEP 3A: SENSITIVITY ANALYSIS", 
                 "\u00b11 Mb Broad Window vs.\n\u00b1200 kb Strict Window (GIPR)")
s3_b <- draw_box(3.8, 4.0, 6.4, 2.2, cols$step3b, "STEP 3B: PLEIOTROPY SCREENING", 
                 "SNP-level Wald Diagnostics\nInstrument PheWAS (Chromosome 19q13)")

# Step 4: Synthesis
s4 <- draw_box(0, 1.2, 14, 1.5, cols$step4, "STEP 4: EVIDENCE SYNTHESIS & TRIANGULATION", 
               "Alignment between independent discovery cohorts, diagnostic stability,\nand regional pleiotropy characterization.")

boxes <- rbind(s1, s2_a, s2_b, s3_a, s3_b, s4)

# Arrows
arrows <- data.frame(
  x    = c(0, 0, -3.8, 3.8, -3.8, 3.8),
  y    = c(9.2, 9.2, 6.0, 6.2, 2.9, 2.9), # Adjusted for new spacing
  xend = c(-3.8, 3.8, -3.8, 3.8, -3.8, 3.8),
  yend = c(8.4, 8.4, 5.1, 5.1, 2.0, 2.0)
)

p <- ggplot() +
  # Draw Header Label (Separate styling)
  geom_text(data = header, aes(x = x, y = y, label = label), 
            colour = cols$header, fontface = "bold", size = 4.2, hjust = 0.5) +
  
  # Draw Boxes
  geom_rect(data = boxes, aes(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2, fill = fill),
            colour = cols$border, linewidth = 0.4) +
  
  # Main Labels
  geom_text(data = boxes, aes(x = x, y = y + h/4, label = label), 
            colour = "white", fontface = "bold", size = 3.6) +
  # Subtitles
  geom_text(data = boxes, aes(x = x, y = y - h/6, label = subtitle), 
            colour = "white", size = 3.0, lineheight = 0.95) +
  
  # Arrows
  geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.08, "inches"), type = "closed"),
               colour = cols$border, linewidth = 0.5) +
  
  # Horizontal connection Step 2
  geom_segment(aes(x = -0.6, y = 7.2, xend = 0.6, yend = 7.2),
               linetype = "dashed", colour = "white", linewidth = 0.4) +
  
  scale_fill_identity() +
  coord_cartesian(xlim = c(-8, 8), ylim = c(0, 12), expand = FALSE) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(10, 10, 10, 10, "pt")
  )

# --- Final Save ---
dir.create(plot_outdir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(plot_outdir, "workflow_diagram.pdf"), p, width = 10, height = 7, device = cairo_pdf)
ggsave(file.path(plot_outdir, "workflow_diagram.png"), p, width = 10, height = 7, dpi = 300)

message("Final journal-style workflow diagram (Steps 1-4) saved to ", plot_outdir)
