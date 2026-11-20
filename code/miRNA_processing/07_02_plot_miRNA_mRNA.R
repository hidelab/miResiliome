rm(list = ls())

# Usage (optional setup number):
# Rscript code/_miRNA_processing/07_02_plot_miRNA_mRNA.R [setup]
#   args[1] = setup number used for output<setup>/ (default: "")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || is.na(args[1])) {
  setup <- ""
  message("No setup provided; defaulting to output/ directory.")
} else {
  setup <- args[1]
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(stringr)
})

options(stringsAsFactors = FALSE)

# Directories (match prior scripts' structure)
output.dir <- paste0("output", setup, "/")
coexp.dir  <- file.path(output.dir, "coexpr", "mirna_mrna")
infile     <- file.path(coexp.dir, "miRNA_mRNA_cor_by_group_tarbase.csv")

# Fallback: if setup-suffixed path missing, try plain output/
if (!file.exists(infile)) {
  alt_infile <- file.path("output", "coexpr", "mirna_mrna", "miRNA_mRNA_cor_by_group_tarbase.csv")
  if (file.exists(alt_infile)) {
    infile <- alt_infile
    coexp.dir <- dirname(infile)
  } else {
    stop("Input file not found: ", infile)
  }
}

plot.dir <- file.path(coexp.dir, "plots")
if (!dir.exists(plot.dir)) dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)

cat("Reading:", infile, "\n")
df <- suppressMessages(readr::read_csv(infile, show_col_types = FALSE))

# Expect columns: mirna_name, gene_id, gene_name, r_AD, p_adj_AD, r_nAD, p_adj_nAD, r_hat, p_adj_consensus, ...
req_cols <- c("mirna_name", "gene_id", "gene_name",
              "r_AD", "p_adj_AD", "r_nAD", "p_adj_nAD",
              "r_hat", "p_adj_consensus")
missing <- setdiff(req_cols, colnames(df))
if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))

# Establish a universal color scale range across all plots
# Compute signed -log10(FDR) for each group globally, then take a symmetric max
signed_vals_all <- c(
  with(df, -log10(p_adj_consensus) * sign(r_hat))
)
if (any(is.finite(signed_vals_all))) {
  global_signed_max <- max(abs(signed_vals_all[is.finite(signed_vals_all)]), na.rm = TRUE)
  # Guard against degenerate zero range
  if (!is.finite(global_signed_max) || global_signed_max <= 0) {
    global_signed_max <- 1
  }
} else {
  global_signed_max <- 1
}
color_limits <- c(-global_signed_max, global_signed_max)

# Helper: for a given miRNA, get its top-N genes by adjusted p-value for a group
prep_top_genes_for_miRNA <- function(data, miRNA, r_col, p_col, top_n = 30) {
  stopifnot(all(c(r_col, p_col) %in% colnames(data)))
  d <- data %>%
    dplyr::filter(mirna_name == miRNA) %>%
    dplyr::filter(is.finite(.data[[r_col]]), is.finite(.data[[p_col]])) %>%
    dplyr::mutate(gene_label = dplyr::if_else(is.na(gene_name) | gene_name == "",
                                              gene_id, gene_name)) %>%
    dplyr::distinct(gene_id, .keep_all = TRUE) %>%
    dplyr::arrange(.data[[p_col]], dplyr::desc(abs(.data[[r_col]]))) %>%
    dplyr::slice_head(n = min(top_n, nrow(.))) %>%
    dplyr::mutate(
      r = .data[[r_col]],
      p_adj = .data[[p_col]],
      signed_log10_fdr = -log10(p_adj) * sign(r),
      gene_label = forcats::fct_reorder(gene_label, r, .fun = identity, .desc = FALSE)
    )
  d
}

plot_miRNA_group <- function(top_df, miRNA, title_suffix, outfile) {
  if (nrow(top_df) == 0) {
    warning("No rows to plot for ", miRNA, " (", title_suffix, ")")
    return(invisible(NULL))
  }
  p <- ggplot(top_df, aes(x = r, y = gene_label, color = signed_log10_fdr)) +
    geom_point(size = 3, alpha = 0.9) +
    scale_x_continuous(limits = c(-.5, .5), breaks = seq(-.5, .5, by = 0.5)) +
    scale_color_gradient2(
      low = "#3B73B9", mid = "#E0E0E0", high = "#C43B3B", midpoint = 0,
      limits = color_limits, oob = scales::squish,
      name = "signed -log10(FDR)"
    ) +
    labs(
      title = paste0("miRNA–mRNA coexpression: ", miRNA, " (", title_suffix, ")"),
      subtitle = "Top genes by adjusted p-value for this miRNA",
      x = "Correlation",
      y = "Genes (ranked by correlation)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.title.y = element_text(margin = margin(r = 8)),
      axis.title.x = element_text(margin = margin(t = 6)),
      plot.title = element_text(face = "bold")
    ) + coord_fixed(ratio = .4)

  ggsave(filename = outfile, plot = p, width = 7.5, height = 8)
  cat("Saved:", outfile, "\n")
}

# Iterate over miRNAs and produce three plots per miRNA
all_miRNAs <- sort(unique(df$mirna_name))

for (mi in all_miRNAs) {
  # AD
  d_ad  <- prep_top_genes_for_miRNA(df, miRNA = mi, r_col = "r_AD",  p_col = "p_adj_AD",        top_n = 30)
  f_ad  <- file.path(plot.dir, paste0(gsub("[/\\\\ ]", "_", mi), "_AD_top30.pdf"))
  plot_miRNA_group(d_ad, miRNA = mi, title_suffix = "AD", outfile = f_ad)

  # nAD
  d_nad <- prep_top_genes_for_miRNA(df, miRNA = mi, r_col = "r_nAD", p_col = "p_adj_nAD",       top_n = 30)
  f_nad <- file.path(plot.dir, paste0(gsub("[/\\\\ ]", "_", mi), "_nAD_top30.pdf"))
  plot_miRNA_group(d_nad, miRNA = mi, title_suffix = "non-AD", outfile = f_nad)

  # Consensus
  d_hat <- prep_top_genes_for_miRNA(df, miRNA = mi, r_col = "r_hat", p_col = "p_adj_consensus", top_n = 30)
  f_hat <- file.path(plot.dir, paste0(gsub("[/\\\\ ]", "_", mi), "_consensus_top30.pdf"))
  plot_miRNA_group(d_hat, miRNA = mi, title_suffix = "Consensus", outfile = f_hat)
}

cat("Done. Per-miRNA plots are in:", plot.dir, "\n")
