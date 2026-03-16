#!/usr/bin/env Rscript
# Figure 3 / Supplementary: PheWAS pleiotropy contrast for GIPR instruments
#
# Compares PheWAS burden between the two GIPR broad-window cis instruments:
#   rs17561351 — EXCLUDED by strict ±200 kb window (high pleiotropy, 136 hits)
#   rs10407429 — RETAINED by strict window  (lower pleiotropy, 39 hits)
#
# Left panel:  horizontal dot plot of top trait associations per domain
# Right panel: domain-level hit-count bar comparison
#
# This directly supports the manuscript's Step C:
# "PheWAS screening explains WHY the broad window introduces heterogeneity"
#
# Input: results/causal/instrument_phewas_hits_p1e8.tsv (full instrument PheWAS)

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

phewas_path <- get_arg("--in",
  "results/causal/instrument_phewas_hits_p1e8.tsv")
plot_outdir <- get_arg("--outdir", "plots/causal/drug_target_mr")

if (!file.exists(phewas_path)) stop("Missing: ", phewas_path)

df <- read_tsv_any(phewas_path)
df$p    <- as.numeric(df$p)
df$beta <- as.numeric(df$beta)

# --- Focus on the two GIPR broad-window instruments ---
snps_gipr <- c("rs17561351", "rs10407429")
df <- df[df$rsid %in% snps_gipr, ]
df$neglog10p <- -log10(pmax(df$p, 1e-300))

# SNP role labels
snp_labels <- c(
  "rs17561351" = "rs17561351 (excluded by strict window)",
  "rs10407429" = "rs10407429 (retained by strict window)"
)
df$snp_label <- snp_labels[df$rsid]

# --- Classify trait domains ---
classify_domain <- function(trait) {
  trait_lc <- tolower(trait)
  domain <- rep("Other", length(trait_lc))
  domain[grepl("cholesterol|ldl|hdl|triglyceride|lipid|apolipoprotein|apob|apoa|vldl|idl|cholesteryl|phospholipid",
               trait_lc)] <- "Lipids / Lipoproteins"
  domain[grepl("c-reactive|crp|inflam|interleukin|cytokine",
               trait_lc)] <- "Inflammation"
  domain[grepl("platelet|mean platelet|plt|mpv|thrombocyte",
               trait_lc)] <- "Platelet indices"
  domain[grepl("alzheimer|dementia|cognit|neurodegen|brain|memory",
               trait_lc)] <- "Neurocognitive"
  domain[grepl("vitamin d|hydroxyvitamin",
               trait_lc)] <- "Vitamin D"
  domain[grepl("bmi|body mass|waist|obesity|fat|adipos",
               trait_lc)] <- "Anthropometric"
  domain[grepl("glucose|hba1c|insulin|diabet|glyc|hemoglobin a1c|haemoglobin",
               trait_lc)] <- "Glycaemic"
  domain[grepl("blood pressure|hypertension|systolic|diastolic",
               trait_lc)] <- "Blood pressure"
  domain[grepl("eosinophil|neutrophil|lymphocyte|monocyte|white blood|wbc|basophil",
               trait_lc)] <- "Blood cell counts"
  domain[grepl("eqtl|ensg",
               trait_lc)] <- "eQTL / Gene expression"
  domain
}

df$domain <- classify_domain(df$trait)

# --- Clean trait names ---
clean_trait <- function(s) {
  s <- gsub("\\s*\\(UKB data field \\d+\\)", "", s)
  s <- gsub(" levels$", "", s)
  s <- gsub("^Direct ", "", s)
  s <- gsub("^Serum ", "", s)
  # Remove ENSG IDs — replace with "Gene expression (ENSG...)"
  s <- ifelse(grepl("^ENSG\\d+$", s), paste0("Gene expr. (", s, ")"), s)
  s <- trimws(s)
  ifelse(nchar(s) > 52, paste0(substr(s, 1, 49), "..."), s)
}
df$trait_clean <- clean_trait(df$trait)

# --- De-duplicate: keep most significant per unique trait per SNP ---
df <- df[order(df$p), ]
df <- df[!duplicated(paste(df$rsid, df$trait_clean)), ]

# --- High contrast palette (no yellow on white) ---
domain_cols <- c(
  "Lipids / Lipoproteins"  = "#C44E52",
  "Inflammation"           = "#8C564B",
  "Platelet indices"       = "#2CA02C",
  "Neurocognitive"         = "#1F77B4",
  "Vitamin D"              = "#9467BD",
  "Anthropometric"         = "#17BECF",
  "Glycaemic"              = "#FF7F0E",
  "Blood pressure"         = "#7F7F7F",
  "Blood cell counts"      = "#BCBD22",
  "eQTL / Gene expression" = "#E377C2",
  "Other"                  = "#AAAAAA"
)

# =====================================================================
# PANEL A: Horizontal dot plot — top 4 traits per domain per SNP
# =====================================================================
df_top <- do.call(rbind, lapply(split(df, paste(df$rsid, df$domain)), function(x) {
  head(x[order(x$p), ], 2)
}))
# Cap per SNP for readability
cap_per_snp <- function(d, n_max = 15) {
  d <- d[order(d$p), ]
  head(d, n_max)
}
df_top <- do.call(rbind, lapply(split(df_top, df_top$rsid), cap_per_snp))

# Facet order: excluded SNP on top
df_top$snp_label <- factor(df_top$snp_label,
  levels = c(snp_labels["rs17561351"], snp_labels["rs10407429"]))

# Row ordering within facets
df_top <- df_top[order(df_top$snp_label, df_top$neglog10p), ]
df_top$trait_clean <- factor(df_top$trait_clean,
  levels = unique(df_top$trait_clean))

p_a <- ggplot(df_top, aes(x = neglog10p, y = trait_clean, colour = domain)) +
  geom_vline(xintercept = -log10(5e-8), linetype = "dashed",
             colour = "grey65", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = neglog10p, yend = trait_clean),
               linewidth = 0.3, alpha = 0.45) +
  geom_point(aes(size = abs(beta)), alpha = 0.85) +
  scale_colour_manual(values = domain_cols, name = "Trait domain", drop = TRUE) +
  scale_size_continuous(range = c(1.5, 4.5),
                        name = expression("|" * beta * "|"),
                        breaks = c(0.05, 0.1, 0.5)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  facet_wrap(~ snp_label, ncol = 1, scales = "free_y") +
  labs(x = expression(-log[10](italic(p))), y = NULL) +
  theme_bw(base_size = 9) +
  theme(
    text = element_text(family = "sans"),
    strip.text = element_text(face = "bold", size = 8.5, hjust = 0),
    strip.background = element_rect(fill = "grey96", colour = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "grey40", linewidth = 0.4),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7.5),
    legend.position = "none",
    plot.margin = margin(4, 4, 4, 4, "pt")
  )

# =====================================================================
# PANEL B: Domain-level hit-count comparison (grouped bar)
# =====================================================================
counts <- as.data.frame(table(df$rsid, df$domain), stringsAsFactors = FALSE)
names(counts) <- c("rsid", "domain", "n_hits")
counts <- counts[counts$n_hits > 0, ]
counts$snp_label <- snp_labels[counts$rsid]
counts$snp_label <- factor(counts$snp_label,
  levels = c(snp_labels["rs17561351"], snp_labels["rs10407429"]))

# Order domains by total hits descending
domain_totals <- tapply(counts$n_hits, counts$domain, sum)
counts$domain <- factor(counts$domain,
  levels = names(sort(domain_totals, decreasing = FALSE)))

snp_fill <- c(
  "rs17561351 (excluded by strict window)" = "#D55E00",
  "rs10407429 (retained by strict window)"  = "#0072B2"
)

p_b <- ggplot(counts, aes(x = n_hits, y = domain, fill = snp_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_text(aes(label = n_hits),
            position = position_dodge(width = 0.7),
            hjust = -0.2, size = 2.5, show.legend = FALSE) +
  scale_fill_manual(values = snp_fill, name = "GIPR instrument") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Number of PheWAS hits (p < 1e-8)", y = NULL,
       subtitle = "Domain-level pleiotropy burden") +
  theme_bw(base_size = 9) +
  theme(
    text = element_text(family = "sans"),
    plot.subtitle = element_text(face = "bold", size = 9, hjust = 0),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "grey40", linewidth = 0.4),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 7.5),
    axis.text.x = element_text(size = 7.5),
    legend.position = "bottom",
    legend.title = element_text(size = 7.5, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.3, "cm"),
    plot.margin = margin(4, 8, 4, 4, "pt")
  )

# =====================================================================
# Combine panels
# =====================================================================
fig <- plot_grid(
  p_a, p_b,
  ncol = 2, rel_widths = c(1.1, 0.9),
  labels = c("A", "B"), label_size = 11
)

# Add a shared title
title_grob <- ggdraw() +
  draw_label(
    "Instrument PheWAS: GIPR broad-window variants (rs17561351 vs rs10407429)",
    fontface = "bold", size = 10, x = 0.02, hjust = 0
  )

fig_final <- plot_grid(title_grob, fig, ncol = 1, rel_heights = c(0.05, 1))

dir.create(plot_outdir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(plot_outdir, "phewas_contrast.pdf"), fig_final,
       width = 11, height = 7, device = cairo_pdf)
ggsave(file.path(plot_outdir, "phewas_contrast.png"), fig_final,
       width = 11, height = 7, dpi = 300)
message("PheWAS contrast figure saved to ", plot_outdir)
