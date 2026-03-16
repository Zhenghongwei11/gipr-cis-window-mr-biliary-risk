#!/usr/bin/env Rscript
# Figure 4: Regional pleiotropy locus plot for chr19q13 around GIPR (GRCh37).
#
# Visualizes gene density within the broad cis window (±1 Mb) and highlights the
# strict GIPR window (±200 kb), alongside the two broad-window HbA1c instruments
# (rs17561351 excluded by strict window; rs10407429 retained).

source("scripts/mr/lib_io.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

set.seed(1)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) stop("Missing value for ", flag)
  args[idx + 1]
}

genes_path <- get_arg("--genes",
  "docs/source_data/Source_Data_12_chr19q13_gene_coordinates_grch37.tsv")
inst_broad_path <- get_arg("--inst-broad",
  "docs/source_data/Source_Data_5_GIPR_instruments_wide1mb.tsv")
inst_strict_path <- get_arg("--inst-strict",
  "docs/source_data/Source_Data_6_GIPR_instruments_strict200kb.tsv")
plot_outdir <- get_arg("--outdir", "plots/causal/drug_target_mr")

if (!file.exists(genes_path)) stop("Missing: ", genes_path)
if (!file.exists(inst_broad_path)) stop("Missing: ", inst_broad_path)
if (!file.exists(inst_strict_path)) stop("Missing: ", inst_strict_path)

genes <- read_tsv_any(genes_path)
genes$start <- as.integer(genes$start)
genes$end <- as.integer(genes$end)

gipr <- genes[genes$gene_symbol == "GIPR", , drop = FALSE]
if (nrow(gipr) != 1) stop("Expected exactly one GIPR gene row in: ", genes_path)

cis_broad <- 1000000L
cis_strict <- 200000L

broad_start <- max(1L, as.integer(gipr$start) - cis_broad)
broad_end <- as.integer(gipr$end) + cis_broad
strict_start <- max(1L, as.integer(gipr$start) - cis_strict)
strict_end <- as.integer(gipr$end) + cis_strict

# Keep protein-coding genes for readability.
genes <- genes[genes$biotype == "protein_coding", , drop = FALSE]
genes <- genes[genes$end >= broad_start & genes$start <= broad_end, , drop = FALSE]
genes <- genes[order(genes$start, genes$end, genes$gene_symbol), , drop = FALSE]

assign_rows <- function(df, min_gap = 10000L) {
  if (nrow(df) < 1) return(df)
  row_ends <- integer(0)
  rows <- integer(nrow(df))
  for (i in seq_len(nrow(df))) {
    s <- as.integer(df$start[i])
    e <- as.integer(df$end[i])
    placed <- FALSE
    if (length(row_ends) > 0) {
      for (r in seq_along(row_ends)) {
        if (s > (row_ends[r] + min_gap)) {
          rows[i] <- r
          row_ends[r] <- e
          placed <- TRUE
          break
        }
      }
    }
    if (!placed) {
      row_ends <- c(row_ends, e)
      rows[i] <- length(row_ends)
    }
  }
  df$row <- rows
  df
}

genes <- assign_rows(genes, min_gap = 12000L)
max_row <- ifelse(nrow(genes) > 0, max(genes$row), 1L)

inst_broad <- read_tsv_any(inst_broad_path)
inst_strict <- read_tsv_any(inst_strict_path)
inst_broad$pos <- as.integer(inst_broad$pos)
inst_strict$pos <- as.integer(inst_strict$pos)

inst_broad$role <- ifelse(inst_broad$snp %in% inst_strict$snp, "Retained (strict window)", "Excluded (strict window)")
inst_broad$role <- factor(inst_broad$role, levels = c("Excluded (strict window)", "Retained (strict window)"))

inst_labels <- c(
  "rs17561351" = "rs17561351 (excluded)",
  "rs10407429" = "rs10407429 (retained)"
)
inst_broad$label <- inst_labels[inst_broad$snp]

notable <- c("GIPR", "APOE", "TOMM40", "PVRL2", "APOC1", "APOC2", "APOC4", "BCL3", "BCAM")
genes$is_notable <- genes$gene_symbol %in% notable
genes$is_target <- genes$gene_symbol == "GIPR"

genes$gene_mid <- (genes$start + genes$end) / 2
inst_broad$y <- max_row + 1.8

# Nudge gene labels off their segment line (strand-aware).
genes$label_y <- genes$row + ifelse(as.integer(genes$strand) == 1L, -0.28, 0.28)

target_cols <- c("Target (GIPR)" = "#D55E00", "Other genes" = "grey35")
inst_cols <- c("Excluded (strict window)" = "#D55E00", "Retained (strict window)" = "#0072B2")

genes$gene_group <- ifelse(genes$is_target, "Target (GIPR)", "Other genes")
genes$gene_group <- factor(genes$gene_group, levels = c("Target (GIPR)", "Other genes"))

p <- ggplot() +
  annotate("rect",
    xmin = strict_start, xmax = strict_end,
    ymin = -Inf, ymax = Inf,
    fill = "#56B4E9", alpha = 0.12
  ) +
  geom_segment(
    data = genes,
    aes(x = start, xend = end, y = row, yend = row, colour = gene_group),
    linewidth = 0.7, lineend = "round"
  ) +
  geom_point(
    data = inst_broad,
    aes(x = pos, y = y, colour = role),
    size = 2.2, alpha = 0.95
  ) +
  geom_text_repel(
    data = inst_broad,
    aes(x = pos, y = y, label = label, colour = role),
    size = 2.7,
    min.segment.length = 0,
    segment.size = 0.25,
    box.padding = 0.25,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  geom_label_repel(
    data = genes[genes$is_notable, , drop = FALSE],
    aes(x = gene_mid, y = label_y, label = gene_symbol),
    size = 2.5,
    colour = "grey10",
    fill = "white",
    label.size = 0,
    min.segment.length = 0,
    segment.size = 0.2,
    box.padding = 0.2,
    point.padding = 0.15,
    max.overlaps = Inf
  ) +
  geom_segment(
    aes(x = broad_start, xend = broad_end, y = max_row + 2.8, yend = max_row + 2.8),
    linewidth = 0.7, colour = "grey25"
  ) +
  geom_segment(
    aes(x = strict_start, xend = strict_end, y = max_row + 2.35, yend = max_row + 2.35),
    linewidth = 1.2, colour = "#0072B2"
  ) +
  annotate("text", x = broad_start, y = max_row + 3.05, hjust = 0,
    label = "Broad cis window (±1 Mb)", size = 3.0, colour = "grey20"
  ) +
  annotate("text", x = strict_start, y = max_row + 2.6, hjust = 0,
    label = "Strict GIPR window (±200 kb)", size = 3.0, colour = "#0072B2"
  ) +
  scale_colour_manual(values = c(target_cols, inst_cols), guide = "none") +
  coord_cartesian(xlim = c(broad_start, broad_end), clip = "off") +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 7),
    labels = function(x) sprintf("%.1f", x / 1e6)
  ) +
  scale_y_reverse(expand = expansion(mult = c(0.02, 0.18))) +
  labs(
    x = "Genomic position (chr19, GRCh37; Mb)",
    y = NULL,
    title = "Regional pleiotropy context around GIPR (chr19q13)",
    subtitle = "Gene-dense broad window (±1 Mb) vs strict GIPR window (±200 kb), with broad-window instruments overlaid"
  ) +
  theme_bw(base_size = 9) +
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
    plot.subtitle = element_text(size = 8.5, hjust = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_rect(colour = "grey40", linewidth = 0.4),
    plot.margin = margin(6, 20, 6, 6, "pt")
  )

dir.create(plot_outdir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(plot_outdir, "regional_pleiotropy_locus.pdf"), p,
  width = 11, height = 6.2, device = cairo_pdf
)
ggsave(file.path(plot_outdir, "regional_pleiotropy_locus.png"), p,
  width = 11, height = 6.2, dpi = 300
)
message("Regional pleiotropy locus plot saved to ", plot_outdir)
