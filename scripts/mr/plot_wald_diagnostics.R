#!/usr/bin/env Rscript
# Figure 2: SNP-level Wald ratio diagnostics for GIPR → cholelithiasis
# Panel A: broad cis-window (Source Data 7)   — 2 SNPs, pooled IVW overlaid
# Panel B: strict cis-window (Source Data 8)  — 1 SNP (single-instrument Wald)
# Reads primary results for pooled IVW estimate overlay.

source("scripts/mr/lib_io.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) stop("Missing value for ", flag)
  args[idx + 1]
}

wald_broad_path  <- get_arg("--wald-broad",
  "docs/source_data/Source_Data_7_wald_GIPR_cholelithiasis_wide1mb.tsv")
wald_strict_path <- get_arg("--wald-strict",
  "docs/source_data/Source_Data_8_wald_GIPR_cholelithiasis_strict200kb.tsv")
primary_broad    <- get_arg("--primary-broad",
  "docs/source_data/Source_Data_1_primary_results_hba1c_wide1mb.tsv")
primary_strict   <- get_arg("--primary-strict",
  "docs/source_data/Source_Data_3_primary_results_hba1c_strictGIPR200kb.tsv")
plot_outdir      <- get_arg("--outdir", "plots/causal/drug_target_mr")

for (p in c(wald_broad_path, wald_strict_path)) {
  if (!file.exists(p)) stop("Missing: ", p)
}

# --- Read Wald data ---
read_wald <- function(path) {
  df <- read_tsv_any(path)
  num_cols <- c("wald_or", "wald_or_ci_lower", "wald_or_ci_upper", "wald_pvalue")
  for (col in num_cols) df[[col]] <- as.numeric(df[[col]])
  df
}

wald_broad  <- read_wald(wald_broad_path)
wald_strict <- read_wald(wald_strict_path)

# --- Read pooled IVW for overlay ---
get_pooled <- function(path, target = "GIPR", outcome = "cholelithiasis_primary") {
  if (!file.exists(path)) return(NULL)
  df <- read_tsv_any(path)
  row <- df[df$target_symbol == target & df$outcome_name == outcome, ]
  if (nrow(row) == 0) return(NULL)
  list(
    or = as.numeric(row$or[1]),
    ci_lo = as.numeric(row$or_ci_lower[1]),
    ci_hi = as.numeric(row$or_ci_upper[1])
  )
}

pooled_broad  <- get_pooled(primary_broad)
pooled_strict <- get_pooled(primary_strict)

# --- Colour-blind safe palette ---
snp_cols <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#F0E442")

# --- Plotting function ---
make_wald_panel <- function(df, pooled, panel_label, show_legend = FALSE) {
  df$snp_label <- df$snp
  df$snp_label <- factor(df$snp_label, levels = rev(df$snp_label))

  # Annotation text
  df$annot <- sprintf("OR %.2f (%.2f\u2013%.2f)\np = %.1e",
                       df$wald_or, df$wald_or_ci_lower, df$wald_or_ci_upper, df$wald_pvalue)

  p <- ggplot(df, aes(x = wald_or, y = snp_label)) +
    # Null line
    geom_vline(xintercept = 1, linetype = "solid", colour = "grey70", linewidth = 0.4)

  # Pooled IVW estimate as shaded band
  if (!is.null(pooled) && nrow(df) > 1) {
    p <- p +
      annotate("rect",
               xmin = pooled$ci_lo, xmax = pooled$ci_hi,
               ymin = -Inf, ymax = Inf,
               fill = "#E8E8E8", alpha = 0.5) +
      geom_vline(xintercept = pooled$or, linetype = "dashed",
                 colour = "grey40", linewidth = 0.5)
  }

  p <- p +
    geom_pointrange(aes(xmin = wald_or_ci_lower, xmax = wald_or_ci_upper, colour = snp),
                    size = 1.5, linewidth = 0.7) +
    scale_colour_manual(values = setNames(snp_cols[seq_len(nrow(df))], df$snp),
                        name = "Instrument SNP") +
    scale_x_continuous(
      trans = "log10",
      breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10),
      labels = c("0.1", "0.25", "0.5", "1", "2", "4", "10")
    ) +
    labs(x = "Wald ratio OR (95% CI)", y = NULL,
         subtitle = panel_label) +
    theme_bw(base_size = 9) +
    theme(
      text = element_text(family = "sans"),
      plot.subtitle = element_text(face = "bold", size = 10, hjust = 0),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "grey40", linewidth = 0.4),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(size = 9, face = "italic"),
      axis.text.x = element_text(size = 8),
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7.5, face = "italic"),
      plot.margin = margin(6, 12, 4, 4, "pt")
    )
  p
}

# --- Panel A: Broad window ---
p_a <- make_wald_panel(wald_broad, pooled_broad,
  "A. Broad cis-window (\u00b11 Mb) \u2014 GIPR \u2192 Cholelithiasis",
  show_legend = TRUE)

# --- Panel B: Strict window ---
p_b <- make_wald_panel(wald_strict, pooled_strict,
  "B. Strict cis-window (\u00b1200 kb) \u2014 GIPR \u2192 Cholelithiasis",
  show_legend = TRUE)

# --- Combine panels ---
fig2 <- plot_grid(
  p_a, p_b,
  ncol = 1, rel_heights = c(1, 0.65),
  align = "v", axis = "lr"
)

dir.create(plot_outdir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(plot_outdir, "wald_diagnostics.pdf"), fig2,
       width = 6.5, height = 4.5, device = cairo_pdf)
ggsave(file.path(plot_outdir, "wald_diagnostics.png"), fig2,
       width = 6.5, height = 4.5, dpi = 300)
message("Figure 2 saved to ", plot_outdir)
