#!/usr/bin/env Rscript

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(stringr)
  library(tools)
  library(dplyr)
  library(limma)
  library(psych)
  library(CovariateAnalysis)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(broom)
  library(VGAM)
  library(foreach)
  library(ggrepel)
  library(patchwork)
  library(cowplot)
})

doParallel::registerDoParallel(10)
source("code/00-02-functions2.R")

# ---- Parameters ----
input_dir   <- file.path("output", "pcxn_out_pfc")
file_glob   <- "*clustNo_*.csv"
out_dir     <- file.path(input_dir, "Figures")

id_col      <- "x"
pval_col    <- "aggInv_pval"
fdr_col     <- "aggInv_fdr"
size_col    <- "cluster_hits"

top_n       <- 30
fdr_thresh  <- 0.1
size_base   <- 0.5
size_range  <- c(2.0, 9.0)

col_low     <- "#FFCCCC"
col_high    <- "#B30000"
col_ns      <- "#BFBFBF"

x_label     <- sprintf("-log10(%s)", "PanomiR FDR")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load data for proportional odds plots ----
covs   <- readRDS("data/preprocessed/clinicalTotal.RDS")
EXPRS  <- readRDS("output/DEPaths_bulk_update/ROSMAP_DLPFC_ADvsControl_Residuals.RDS")

all.ad  <- covs[covs$cogdx.x %in% c(1, 2, 4), ]
gene.AD <- EXPRS[, as.character(all.ad$SampleID)]

clean_names <- read.csv("data/MSigDB_display_names_v2024.csv")
clean_names <- clean_names[, c("standard_name", "display_name")]

prepare_prop_odds_table <- function(df, clean_names) {
  df <- dplyr::left_join(df, clean_names[, c("standard_name", "display_name")],
                         by = c("gene" = "standard_name"))
  df <- df %>%
    dplyr::select(gene, display_name, p.adj, Dir, "Pr...z..")
  names(df) <- c("gene", "display_name", "p.adj", "Dir", "pval")
  return(df)
}

cog.paths        <- read.csv("output/Bulk_paths_prop_odds/cog_paths.csv")
cog.paths        <- prepare_prop_odds_table(cog.paths, clean_names)

pathway_clusters <- read.csv("output/pcxn_out_pfc/Clustering_Cogpaths.csv")
pathway_clusters <- pathway_clusters[pathway_clusters$cluster %in% c(1:6), ]

my_stars <- read.csv("output/deMirs/Heatmap_deMirs_sexSpec_inputData.csv")

# ---- Helper: read PanomiR CSV ----
read_prioritization <- function(path) {
  df <- suppressMessages(readr::read_csv(path, show_col_types = FALSE))
  idx_cols <- names(df)[names(df) %in% c("...1", "X1", "X_1", "X")]
  if (length(idx_cols) == 1 && is.numeric(df[[idx_cols]])) {
    df <- df |> dplyr::select(-all_of(idx_cols))
  }
  expected <- c(id_col, size_col, pval_col, fdr_col)
  missing <- setdiff(expected, names(df))
  if (length(missing) > 0) {
    stop(sprintf("Missing expected columns in %s: %s", basename(path), paste(missing, collapse = ", ")))
  }
  df
}

# ---- Dot plot function (from 03-04) ----
make_dot_plot <- function(path, global_hits_max = NULL) {
  df <- read_prioritization(path)

  n_rows <- min(nrow(df), top_n)
  df_top <- df |>
    arrange(.data[[pval_col]]) |>
    slice_head(n = n_rows) |>
    mutate(
      miRNA = .data[[id_col]],
      fdr = .data[[fdr_col]],
      pval = .data[[pval_col]],
      hits = .data[[size_col]],
      xvar = -log10(pmax(fdr, .Machine$double.eps)),
      sig = fdr < fdr_thresh,
      color_value = if_else(sig, -log10(pmax(fdr, .Machine$double.eps)), as.numeric(NA))
    )

  df_top <- mutate(df_top,
                   miRNA = case_when(
                     miRNA %in% my_stars$miRNA[my_stars$star == "***"] ~ paste0("***", miRNA),
                     miRNA %in% my_stars$miRNA[my_stars$star == "^"] ~ paste0("^", miRNA),
                     TRUE ~ miRNA))

  df_top <- df_top |>
    mutate(miRNA = fct_reorder(miRNA, pval, .desc = FALSE)) |>
    arrange(pval) |>
    mutate(miRNA = fct_rev(miRNA))

  clust_num <- sub(".*clustNo_([0-9]+).*", "\\1", basename(path))
  title_txt <- sprintf("Cluster %s top %d predicted miRNAs", clust_num, nrow(df_top))

  max_x <- max(df_top$xvar)

  local_hits_max <- suppressWarnings(max(df[[size_col]], na.rm = TRUE))
  hits_max <- if (is.null(global_hits_max) || !is.finite(global_hits_max)) local_hits_max else global_hits_max
  if (!is.finite(hits_max) || hits_max < 0) hits_max <- 0
  size_limits <- c(size_base, hits_max + size_base)

  if (hits_max > 0) {
    size_breaks_hits <- pretty(c(0, hits_max), n = 4)
    size_breaks_hits <- unique(pmin(pmax(size_breaks_hits, 0), hits_max))
  } else {
    size_breaks_hits <- 0
  }
  size_breaks <- size_breaks_hits + size_base

  p <- ggplot(df_top, aes(x = xvar, y = miRNA)) +
    geom_point(aes(size = hits + size_base, color = color_value), alpha = 0.9) +
    scale_color_gradient(low = col_low, high = col_high, na.value = col_ns,
                         name = sprintf("FDR<%.02f\n(−log10)", fdr_thresh)) +
    xlim(c(0, max_x)) +
    scale_size_continuous(range = size_range, limits = size_limits,
                          breaks = size_breaks, labels = size_breaks_hits,
                          name = "Cluster hits") +
    guides(size = guide_legend(order = 1),
           color = guide_colorbar(order = 2)) +
    labs(title = title_txt, x = x_label, y = "miRNA") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 12),
      axis.text.y = element_text(size = 12),
      legend.position = "right"
    )

  p
}

# ---- Proportional odds plot function (from 03-07) ----
make_prop_odds_plot <- function(cluster_i) {
  sig.mirs <- pathway_clusters[pathway_clusters$cluster == cluster_i, ]$Pathway
  sig.mirs <- cog.paths[cog.paths$gene %in% sig.mirs, ]$gene
  sig.mirs <- sig.mirs[1:pmin(30, length(sig.mirs))]

  tested_gene <- gene.AD[sig.mirs, ]
  tested_gene <- as.data.frame(t(tested_gene))
  tested_gene <- rownames_to_column(tested_gene, var = "SampleID")
  tested_gene <- pivot_longer(tested_gene, cols = sig.mirs)

  my_met <- all.ad[, c("SampleID", "cogdx.x")]
  my_met$cogdx <- factor(my_met$cogdx.x)
  levels(my_met$cogdx) <- c("No CI", "MCI", "AD")

  tested_gene <- dplyr::left_join(tested_gene, my_met, by = "SampleID")
  more_meta   <- cog.paths[, c("gene", "p.adj", "Dir", "pval")]
  tested_gene <- dplyr::left_join(tested_gene, more_meta, by = c("name" = "gene"))
  tested_gene <- dplyr::left_join(tested_gene, clean_names, by = c("name" = "standard_name"))

  tested_gene_norm <- tested_gene %>%
    group_by(name) %>%
    mutate(
      baseline = mean(value[cogdx == "No CI"], na.rm = TRUE),
      adj_expr = value - baseline,
      logp = ifelse(Dir == "NEG", log10(pval), -log10(pval))
    ) %>%
    ungroup()

  summary_df2 <- tested_gene_norm %>%
    group_by(name, cogdx, Dir, logp, display_name) %>%
    summarise(mean_expr = mean(adj_expr, na.rm = TRUE), .groups = "drop")

  lab_df <- summary_df2 %>% filter(cogdx == "AD")
  lab_df$display_name <- ifelse(nchar(lab_df$display_name) > 80,
                                paste0(substr(lab_df$display_name, 1, 77), "..."),
                                lab_df$display_name)

  p <- ggplot(summary_df2,
              aes(x = cogdx, y = mean_expr, group = name, colour = logp)) +
    geom_point(size = 0) +
    geom_line() +
    geom_text_repel(
      data = lab_df,
      aes(label = display_name),
      nudge_x = 0.1,
      direction = "y",
      hjust = 0,
      segment.colour = "black",
      segment.linetype = "dashed",
      box.padding = 0.25
    ) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.6, linetype = "dashed") +
    scale_x_discrete(expand = expansion(mult = c(0.05, 4))) +
    scale_colour_gradient2(
      low = "blue", mid = "white", high = "firebrick",
      midpoint = 0, limits = c(-8, 9), name = "signed log-p"
    ) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    xlab("") +
    ylab("normalized pathway activity") +
    labs(title = paste0("Cluster ", cluster_i, " top pathways")) +
    theme(
      plot.title = element_text(size = 14),
      axis.text = element_text(size = 14),
      axis.text.x = element_text(face = "bold", size = 16, angle = 60, vjust = 1, hjust = 1),
      axis.text.y = element_text(size = 20),
      axis.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 14, angle = 90),
      axis.line.x = element_blank(),
      plot.margin = margin(5.5, 130, 5.5, 5.5)
    )

  p
}

# ---- Main: combine side by side ----
main <- function() {
  files <- Sys.glob(file.path(input_dir, file_glob))
  if (length(files) == 0) stop(sprintf("No files found at %s matching %s", input_dir, file_glob))

  # Compute global hits max for consistent dot-plot sizing

max_vec <- vapply(files, function(f) {
    val <- tryCatch({
      df <- read_prioritization(f)
      suppressWarnings(max(df[[size_col]], na.rm = TRUE))
    }, error = function(e) NA_real_)
    if (is.finite(val)) val else NA_real_
  }, numeric(1))
  global_hits_max <- suppressWarnings(max(max_vec, na.rm = TRUE))
  if (!is.finite(global_hits_max)) global_hits_max <- NA_real_

  # Sort files by cluster number
  clust_nums <- as.integer(sub(".*clustNo_([0-9]+).*", "\\1", basename(files)))
  files <- files[order(clust_nums)]
  clust_nums <- sort(clust_nums)

  combined_plots <- list()

  for (idx in seq_along(files)) {
    cl <- clust_nums[idx]
    message(sprintf("Processing cluster %d ...", cl))

    dot_p  <- make_dot_plot(files[idx], global_hits_max = global_hits_max)
    odds_p <- make_prop_odds_plot(cl)

    pair <- dot_p + odds_p +
      plot_layout(widths = c(.8, 1.8)) +
      plot_annotation(title = sprintf("Cluster %d", cl),
                      theme = theme(plot.title = element_text(face = "bold", size = 16)))

    combined_plots[[idx]] <- pair

    # Save individual pair
    pair_out <- file.path(out_dir, sprintf("cluster_%d_combined.pdf", cl))
    ggsave(pair_out, pair, width = 22, height = 10, device = cairo_pdf, bg = "white")
    message(sprintf("  Saved: %s", pair_out))
  }
}

#   # Save all clusters in a single stacked PDF
#   all_p <- wrap_plots(combined_plots, ncol = 1) &
#     theme(legend.position = "right")
#   all_out <- file.path(out_dir, "all_clusters_combined.pdf")
#   ggsave(all_out, all_p, width = 22, height = 10 * length(combined_plots),
#          device = cairo_pdf, bg = "white")
#   message(sprintf("Saved combined figure: %s", all_out))
# }

if (identical(environment(), globalenv())) {
  main()
}
