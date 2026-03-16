#!/usr/bin/env Rscript
# Figure 1: Broad vs strict cis-window MR forest plot (dual-panel)
# Reads Source Data 1 (broad) and Source Data 3 (strict GIPR)
# plus sensitivity files for heterogeneity annotations.

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

broad_path  <- get_arg("--broad",  "docs/source_data/Source_Data_1_primary_results_hba1c_wide1mb.tsv")
strict_path <- get_arg("--strict", "docs/source_data/Source_Data_3_primary_results_hba1c_strictGIPR200kb.tsv")
sens_broad  <- get_arg("--sens-broad",  "docs/source_data/Source_Data_2_sensitivity_results_hba1c_wide1mb.tsv")
sens_strict <- get_arg("--sens-strict", "docs/source_data/Source_Data_4_sensitivity_results_hba1c_strictGIPR200kb.tsv")
plot_outdir <- get_arg("--outdir", "plots/causal/drug_target_mr")

for (p in c(broad_path, strict_path)) {
  if (!file.exists(p)) stop("Missing: ", p)
}

# --- Read and prepare data ---
read_primary <- function(path) {
  df <- read_tsv_any(path)
  for (col in c("or", "or_ci_lower", "or_ci_upper", "pvalue", "fdr", "nsnp")) {
    df[[col]] <- as.numeric(df[[col]])
  }
  df
}

broad  <- read_primary(broad_path)
strict <- read_primary(strict_path)

# Read heterogeneity p-values
read_het <- function(path) {
  if (!file.exists(path)) return(data.frame())
  s <- read_tsv_any(path)
  s <- s[s$check_name == "heterogeneity_q_p", ]
  s$het_p <- as.numeric(s$check_value)
  s[, c("target_symbol", "outcome_name", "het_p")]
}
het_broad  <- read_het(sens_broad)
het_strict <- read_het(sens_strict)

# --- Outcome display names ---
outcome_labels <- c(
  "cholelithiasis_primary"      = "Cholelithiasis",
  "cholecystectomy_secondary"   = "Cholecystectomy (FinnGen)",
  "cholecystectomy_ukb"         = "Cholecystectomy (UKB)",
  "cholecystitis_finngen"       = "Cholecystitis",
  "acute_pancreatitis_finngen"  = "Acute pancreatitis"
)

# --- Target colours (colour-blind safe) ---
target_cols <- c("GLP1R" = "#0072B2", "GIPR" = "#D55E00", "SLC5A2" = "#009E73")
target_shapes <- c("GLP1R" = 16, "GIPR" = 17, "SLC5A2" = 15)

# --- Build panel data ---
build_panel <- function(df, het_df, panel_label) {
  df$outcome_display <- outcome_labels[df$outcome_name]
  df$outcome_display[is.na(df$outcome_display)] <- df$outcome_name[is.na(df$outcome_display)]
  # Merge heterogeneity
  if (nrow(het_df) > 0) {
    df <- merge(df, het_df, by = c("target_symbol", "outcome_name"), all.x = TRUE)
  } else {
    df$het_p <- NA_real_
  }
  # Significance markers
  df$sig_label <- ifelse(df$fdr < 0.05, "*", "")
  # Row ordering: group by target, then outcome
  target_order <- c("GLP1R", "GIPR", "SLC5A2")
  outcome_order <- rev(names(outcome_labels))
  df$target_symbol <- factor(df$target_symbol, levels = target_order)
  df$outcome_display <- factor(df$outcome_display, levels = rev(outcome_labels))
  # Y-axis label with target
  df$row_label <- paste0(df$target_symbol, "  \u2192  ", df$outcome_display)
  # Sort for consistent vertical order
  df <- df[order(df$target_symbol, df$outcome_display), ]
  df$row_label <- factor(df$row_label, levels = rev(df$row_label))
  df$panel <- panel_label
  df
}

panel_a <- build_panel(broad, het_broad, "A. Broad cis-window (\u00b11 Mb)")
panel_b_strict <- strict[strict$target_symbol == "GIPR", ]
panel_b <- build_panel(panel_b_strict, het_strict[het_strict$target_symbol == "GIPR", ], "B. Strict GIPR cis-window (\u00b1200 kb)")

# --- Common plotting function ---
make_forest <- function(df, x_label = "Odds ratio (95% CI)", x_log = TRUE) {
  # Compute nice axis limits
  all_vals <- c(df$or_ci_lower, df$or_ci_upper)
  x_min <- max(0.01, min(all_vals, na.rm = TRUE) * 0.7)
  x_max <- min(1000, max(all_vals, na.rm = TRUE) * 1.4)

  # Build annotation column: OR (CI), het P
  df$annot <- sprintf("%.2f (%.2f\u2013%.2f)", df$or, df$or_ci_lower, df$or_ci_upper)
  df$het_annot <- ifelse(!is.na(df$het_p) & df$nsnp > 1,
                         sprintf("Q p=%.3f", df$het_p), "")
  # Flag significant heterogeneity
  df$het_flag <- !is.na(df$het_p) & df$het_p < 0.05

  p <- ggplot(df, aes(x = or, y = row_label, colour = target_symbol, shape = target_symbol)) +
    geom_vline(xintercept = 1, linetype = "solid", colour = "grey70", linewidth = 0.4) +
    geom_pointrange(aes(xmin = or_ci_lower, xmax = or_ci_upper),
                    size = 1.2, linewidth = 0.5) +
    # Heterogeneity flag: open diamond behind point if Q p<0.05
    {if (any(df$het_flag))
      geom_point(data = df[df$het_flag, ], aes(x = or, y = row_label),
                 shape = 5, size = 4.5, colour = "grey30", stroke = 0.8, inherit.aes = FALSE)
    } +
    # OR annotation to the right
    geom_text(aes(x = x_max * 0.98, label = annot),
              hjust = 1, size = 2.5, show.legend = FALSE, fontface = "plain") +
    scale_colour_manual(values = target_cols, name = "Target locus") +
    scale_shape_manual(values = target_shapes, name = "Target locus") +
    scale_x_continuous(trans = "log10",
                       breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10),
                       labels = c("0.1", "0.25", "0.5", "1", "2", "4", "10"),
                       limits = c(x_min, x_max),
                       oob = squish) +
    labs(x = x_label, y = NULL) +
    facet_wrap(~ panel, scales = "free_y", ncol = 1) +
    theme_bw(base_size = 9) +
    theme(
      text = element_text(family = "sans"),
      strip.text = element_text(face = "bold", size = 10, hjust = 0),
      strip.background = element_rect(fill = "grey96", colour = NA),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "grey40", linewidth = 0.4),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 8),
      legend.position = "bottom",
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7.5),
      legend.key.size = unit(0.35, "cm"),
      plot.margin = margin(4, 8, 4, 4, "pt")
    )
  p
}

# --- Panel A ---
p_a <- make_forest(panel_a) + theme(legend.position = "none")

# --- Panel B ---
p_b <- make_forest(panel_b) +
  theme(legend.position = "bottom")

# --- Combine ---
fig1 <- plot_grid(
  p_a, p_b,
  ncol = 1, rel_heights = c(1, 0.45),
  align = "v", axis = "lr"
)

dir.create(plot_outdir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(plot_outdir, "forest_primary.pdf"), fig1,
       width = 7.5, height = 8, device = cairo_pdf)
ggsave(file.path(plot_outdir, "forest_primary.png"), fig1,
       width = 7.5, height = 8, dpi = 300)
message("Figure 1 saved to ", plot_outdir)
