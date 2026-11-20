rm(list = ls())

# Usage (optional setup number):
# Rscript code/_miRNA_processing/07_04_miRNA_Pathway_heatmap.R [setup] [r_col] [p_col]
#   args[1] = setup number used for output<setup>/ (default: "")
#   args[2] = correlation column (default: r_hat)
#   args[3] = adjusted p-value column (default: p_adj_consensus)

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
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

options(stringsAsFactors = FALSE)

# Directories (match prior scripts' structure)
output.dir <- paste0("output", setup, "/")
coexp.dir  <- file.path(output.dir, "coexpr", "mirna_path")
infile     <- file.path(coexp.dir, "miRNA_path_cor_by_group_tarbase.csv")

# Fallbacks if setup-suffixed path missing
if (!file.exists(infile)) {
  candidates <- c(
    file.path(coexp.dir, "resilience_miRNA_pathway_cor_by_group_tarbase.csv"),
    file.path("output", "coexpr", "mirna_path", "miRNA_pathway_cor_by_group_resilience.csv"),
    file.path("output", "coexpr", "mirna_path", "resilience_miRNA_pathway_cor_by_group_tarbase.csv")
  )
  infile <- candidates[file.exists(candidates)][1]
  if (is.na(infile) || !file.exists(infile)) {
    stop("Input file not found. Tried: ", paste(c(infile, candidates), collapse = ", "))
  }
}

plot.dir <- file.path(coexp.dir, "plots")
if (!dir.exists(plot.dir)) dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)

cat("Reading:", infile, "\n")
df <- suppressMessages(readr::read_csv(infile, show_col_types = FALSE))

clean_names <- read.csv("data/MSigDB_display_names_v2024.csv")
clean_names <- clean_names[, c("standard_name", "display_name")]
df <- dplyr::left_join(df, clean_names[, c("standard_name", "display_name")],
                       by = c("pathway_name" = "standard_name"))
df$display_name <- ifelse(nchar(df$display_name) > 80,
                          paste0(substr(df$display_name, 1, 77), "..."),
                          df$display_name)



setup_pairs <- list(c("r_hat", "p_adj_consensus"),
                    c("r_AD", "p_adj_AD"),
                    c("r_nAD", "p_adj_nAD"))




r_col <- ifelse(length(args) >= 2 && !is.na(args[2]), args[2], "r_hat")
p_col <- ifelse(length(args) >= 3 && !is.na(args[3]), args[3], "p_adj_consensus")

prep_top_pathways_for_miRNA <- function(data, miRNA, r_col, p_col, top_n = 30) {
  stopifnot(all(c(r_col, p_col) %in% colnames(data)))
  data %>%
    dplyr::filter(mirna_name == miRNA) %>%
    dplyr::filter(is.finite(.data[[r_col]]), is.finite(.data[[p_col]])) %>%
    dplyr::mutate(pathway_label = dplyr::if_else(is.na(display_name) | display_name == "",
                                                 pathway_name, display_name)) %>%
    dplyr::distinct(pathway_name, .keep_all = TRUE) %>%
    dplyr::arrange(.data[[p_col]], dplyr::desc(abs(.data[[r_col]]))) %>%
    dplyr::slice_head(n = min(top_n, nrow(.)))
}



my_df_all <- list()
for(pair in setup_pairs) {
  if (!all(pair %in% colnames(df))) {
    message("Warning: Missing columns for pair: ", paste(pair, collapse = ", "),
            "; available columns: ", paste(colnames(df), collapse = ", "))
  }

  my_df <- df[df[,pair[[2]]] < 0.05,]

  my_df <- my_df %>%
    dplyr::group_by(mirna_name) %>%
    dplyr::summarise(num_sig_pathways = n(),
                     .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(num_sig_pathways))
  # rename num_sig_pathways to include pair[[2]]
  colnames(my_df)[colnames(my_df) == "num_sig_pathways"] <- paste0("num_sig_pathways_", pair[[2]])
  my_df_all[[pair[[2]]]] <- my_df



}

# Merge all counts into a single data frame
df_counts <- Reduce(function(x, y) merge(x, y, by = "mirna_name", all = TRUE), my_df_all)
df_counts[is.na(df_counts)] <- 0
df_counts <- df_counts[order(df_counts$num_sig_pathways_p_adj_consensus,decreasing = T), ]

write.csv(df_counts,
          file = file.path(coexp.dir, "miRNA_pathway_significant_pathway_counts_by_setup.csv"),
          row.names = FALSE)



my_df_all <- list()
for(pair in setup_pairs) {
  if (!all(pair %in% colnames(df))) {
    message("Warning: Missing columns for pair: ", paste(pair, collapse = ", "),
            "; available columns: ", paste(colnames(df), collapse = ", "))
  }

  my_df <- df[abs(df[,pair[[1]]]) > sqrt(.1),]

  my_df <- my_df %>%
    dplyr::group_by(mirna_name) %>%
    dplyr::summarise(num_sig_pathways = n(),
                     .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(num_sig_pathways))
  # rename num_sig_pathways to include pair[[2]]
  colnames(my_df)[colnames(my_df) == "num_sig_pathways"] <- paste0("num_sig_pathways_", pair[[2]])
  my_df_all[[pair[[1]]]] <- my_df

}

# Merge all counts into a single data frame
df_counts <- Reduce(function(x, y) merge(x, y, by = "mirna_name", all = TRUE), my_df_all)
df_counts[is.na(df_counts)] <- 0
df_counts <- df_counts[order(df_counts$num_sig_pathways_p_adj_consensus,decreasing = T), ]

write.csv(df_counts,
          file = file.path(coexp.dir, "miRNA_pathway_significant_pathway_counts_by_R2geq1.csv"),
          row.names = FALSE)








for (pair in setup_pairs) {
  if (all(pair %in% colnames(df))) {
    r_col <- pair[1]
    p_col <- pair[2]

    message("Using columns: r = '", r_col, "', p = '", p_col, "'")
    req_cols <- c("mirna_name", "pathway_name", "display_name", r_col, p_col)

    missing <- setdiff(req_cols, colnames(df))
    if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))

    # Helper: for a given miRNA, get top-N pathways by adjusted p-value

    # Collect top pathways per miRNA
    all_miRNAs <- sort(unique(df$mirna_name))
    top_list <- lapply(all_miRNAs, function(mi) {
      prep_top_pathways_for_miRNA(df, miRNA = mi, r_col = r_col, p_col = p_col, top_n = 50)
    })
    top_df <- dplyr::bind_rows(top_list)
    top_df <- top_df[top_df[[p_col]] <= 0.05, ]
    

    if (nrow(top_df) == 0) stop("No top pathways found; check input file content.")

    top_pathways <- top_df %>%
      dplyr::distinct(pathway_name, display_name) %>%
      dplyr::mutate(pathway_label = dplyr::if_else(is.na(display_name) | display_name == "",
                                                  pathway_name, display_name))

    # Build miRNA x pathway matrix of correlations
    heat_df <- df %>%
      dplyr::filter(pathway_name %in% top_pathways$pathway_name) %>%
      dplyr::filter(is.finite(.data[[r_col]])) %>%
      dplyr::group_by(mirna_name, pathway_name) %>%
      dplyr::summarise(r_val = mean(.data[[r_col]]), .groups = "drop")

    mat <- heat_df %>%
      tidyr::pivot_wider(names_from = mirna_name, values_from = r_val) %>%
      tibble::column_to_rownames("pathway_name") %>%
      as.matrix()

    # Apply display names and order rows by max absolute correlation
    row_labels <- top_pathways$pathway_label[match(rownames(mat), top_pathways$pathway_name)]
    row_labels[is.na(row_labels)] <- rownames(mat)
    rownames(mat) <- make.unique(row_labels)

    row_order <- order(apply(abs(mat), 1, max, na.rm = TRUE), decreasing = TRUE)
    mat <- mat[row_order, , drop = FALSE]

    max_abs <- max(abs(mat), na.rm = TRUE)
    if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 1
    col_fun <- circlize::colorRamp2(c(-max_abs, 0, max_abs),
                                    c("#3B73B9", "#F2F2F2", "#C43B3B"))

    # replace NA with zero 
    mat[is.na(mat)] <- 0
    mat[abs(mat) < min(abs(top_df[,r_col]))] <- 0
    mat <- mat[rowSums(mat) > 0,]
    ht <- Heatmap(
      mat,
      name = r_col,
      col = col_fun,
      na_col = "#E6E6E6",
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_row_names = TRUE,
      clustering_method_rows = "ward",
      clustering_method_columns = "ward",
      show_row_dend = FALSE,
      show_column_dend = FALSE,
      show_column_names = TRUE,
      row_names_gp = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 8),
      row_split = 5,
      column_split = 3,
      width = ncol(mat) * unit(3, "mm"),
      height = nrow(mat) * unit(2, "mm")
    )

    width_in <- max(15, min(20, ncol(mat) * 0.5))
    height_in <- max(10, min(30, nrow(mat) * 0.4))

    out_pdf <- file.path(plot.dir, paste0("miRNA_pathway_heatmap_",r_col,".pdf"))
    out_png <- file.path(plot.dir, paste0("miRNA_pathway_heatmap_",r_col,".png"))

    message("Saving heatmap to:\n  ", out_pdf, "\n  ", out_png)

    pdf(out_pdf, width = width_in, height = height_in)
    draw(ht, heatmap_legend_side = "left")
    invisible(dev.off())

    png(out_png, width = width_in * 300, height = height_in * 300, res = 300)
    draw(ht, heatmap_legend_side = "left")
    invisible(dev.off())

    message("Done.")


   
  }
}



my_tab <- read.csv(file.path(coexp.dir, "miRNA_pathway_significant_pathway_counts_by_setup.csv"))

# Reshape data for grouped bar plot
my_tab_long <- my_tab %>%
  tidyr::pivot_longer(
    cols = starts_with("num_sig_pathways_"),
    names_to = "group",
    values_to = "count"
  ) %>%
  dplyr::mutate(
    group = gsub("num_sig_pathways_p_adj_", "", group),
    group = ifelse(group == "consensus", "consensus",
                   ifelse(group == "AD", "AD", "non-AD")),
    mirna_name = factor(mirna_name, levels = my_tab$mirna_name)
  )

# Create grouped bar plot
p <- ggplot(my_tab_long, aes(x = mirna_name, y = count, fill = group)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = count), position = position_dodge(width = 0.9),
            vjust = -0.3, size = 2.5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  ) +
  labs(
    x = "",
    y = "Number of pathways with significant correlation",
    fill = "Group",
    title = ""
  )

print(p)

# Optionally save the plot
out_bar_plot <- file.path(plot.dir, "miRNA_pathway_counts_barplot.png")
ggsave(out_bar_plot, plot = p, width = 12, height = 6, dpi = 300)
message("Bar plot saved to: ", out_bar_plot)


plot.dir2 <- "output/figures_paper"
out_bar_plot2 <- file.path(plot.dir2, "Figure 4A miRNA_pathway_counts_barplot.jpg")
ggsave(out_bar_plot2, plot = p, width = 12, height = 6, dpi = 600)
message("Bar plot saved to: ", out_bar_plot)



out_bar_plot2 <- file.path(plot.dir2, "Figure 4A miRNA_pathway_counts_barplot.pdf")
ggsave(out_bar_plot2, plot = p, width = 12, height = 6)
message("Bar plot saved to: ", out_bar_plot)


# Expect columns used for correlation/selection


