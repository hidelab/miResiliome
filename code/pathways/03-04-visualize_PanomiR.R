#!/usr/bin/env Rscript

# Visualize PanomiR prioritization outputs
# - Reads CSV tables like output/pcxn_out/*clustNo_*.csv
# - Takes top 30 miRNAs by aggInv_pval
# - Plots dot plot: y = miRNA, x = log(aggInv_fdr)
# - Colors: significant (FDR < 0.05) on a red gradient; others gray
# - Size: proportional to cluster_hits with a positive baseline so zeros are visible

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(stringr)
  library(tools)
  library(dplyr)
})

# ---- Parameters (edit as needed) ----
input_dir <- file.path("output", "pcxn_out_pfc")
file_glob <- "*clustNo_*.csv"
out_dir <- file.path(input_dir, "Figures")

id_col     <- "x"              # miRNA identifier column
pval_col   <- "aggInv_pval"    # column to rank by (lower is better)
fdr_col    <- "aggInv_fdr"     # used for x-axis (log) and significance coloring
size_col   <- "cluster_hits"   # bubble size

top_n      <- 30                # number of miRNAs to display per file
fdr_thresh <- 0.1              # significance threshold
size_base  <- 0.5               # baseline added to cluster_hits so zeros still show
size_range <- c(2.0, 9.0)       # constant visual size range across plots

# Choose color palette for significant points
col_low <- "#FFCCCC"
col_high <- "#B30000"
col_ns <- "#BFBFBF"            # gray for non-significant

# Axis label uses natural log as requested
x_label <- sprintf("-log10(%s)", "PanomiR FDR")

# Ensure output directory exists
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Helper: try to read CSV and normalize columns
read_prioritization <- function(path) {
  df <- suppressMessages(readr::read_csv(path, show_col_types = FALSE))

  # Drop unnamed index column if present (e.g., "...1" or "X1")
  idx_cols <- names(df)[names(df) %in% c("...1", "X1", "X_1", "X")]
  if (length(idx_cols) == 1 && is.numeric(df[[idx_cols]])) {
    df <- df |> dplyr::select(-all_of(idx_cols))
  }

  # Validate expected columns
  expected <- c(id_col, size_col, pval_col, fdr_col)
  missing <- setdiff(expected, names(df))
  if (length(missing) > 0) {
    stop(sprintf("Missing expected columns in %s: %s", basename(path), paste(missing, collapse = ", ")))
  }

  df
}


my_stars <- read.csv("output/deMirs/Heatmap_deMirs_sexSpec_inputData.csv")
# Create a plot for a single file
plot_file <- function(path, global_hits_max = NULL) {
  df <- read_prioritization(path)

  # Select top N by p-value (ascending)
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
                          miRNA %in% my_stars$miRNA[my_stars$star == "***"] ~ paste0("***",miRNA),
                          miRNA %in% my_stars$miRNA[my_stars$star == "^"] ~ paste0("^",miRNA),
                          TRUE ~ miRNA))

  # Order y-axis by p-value, best at top
  df_top <- df_top |>
    mutate(miRNA = fct_reorder(miRNA, pval, .desc = FALSE)) |>
    arrange(pval) |>
    mutate(miRNA = fct_rev(miRNA))

  clust_num <- sub(".*clustNo_([0-9]+).*", "\\1", basename(path))
  title_txt <- sprintf("Cluster %s top %d pathways",
                       clust_num, nrow(df_top))

  print(title_txt)               

  max_x <- max(df_top$xvar)

  # Determine size scale limits to keep sizes comparable across plots
  local_hits_max <- suppressWarnings(max(df[[size_col]], na.rm = TRUE))
  hits_max <- if (is.null(global_hits_max) || !is.finite(global_hits_max)) local_hits_max else global_hits_max
  if (!is.finite(hits_max) || hits_max < 0) hits_max <- 0
  size_limits <- c(size_base, hits_max + size_base)
  # Consistent size breaks and labels (in raw hits) across plots
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
    xlim(c(0,max_x)) +
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

  # Figure size scales with number of rows
  h <- max(4, 0.30 * nrow(df_top) + 1.5)
  w <- 5.5

  # Save both PNG and PDF
  base_name <- tools::file_path_sans_ext(basename(path))
  png_out <- file.path(out_dir, paste0(base_name, ".dotplot.png"))
  pdf_out <- file.path(out_dir, paste0(base_name, ".dotplot.pdf"))

  suppressMessages(ggsave(png_out, p, width = w, height = h, dpi = 300, bg = "white"))
  suppressMessages(ggsave(pdf_out, p, width = w, height = h, device = cairo_pdf, bg = "white"))

  invisible(list(plot = p, data = df_top, files = c(png_out, pdf_out)))
}

# ---- Main ----
main <- function() {
  files <- Sys.glob(file.path(input_dir, file_glob))
  if (length(files) == 0) {
    stop(sprintf("No files found at %s matching %s", input_dir, file_glob))
  }

  # Compute a global maximum for cluster_hits across all files to standardize sizes
  max_vec <- vapply(files, function(f) {
    val <- tryCatch({
      df <- read_prioritization(f)
      suppressWarnings(max(df[[size_col]], na.rm = TRUE))
    }, error = function(e) NA_real_)
    if (is.finite(val)) val else NA_real_
  }, numeric(1))
  global_hits_max <- suppressWarnings(max(max_vec, na.rm = TRUE))
  if (!is.finite(global_hits_max)) global_hits_max <- NA_real_

  message(sprintf("Found %d file(s). Writing plots to %s", length(files), out_dir))
  res <- lapply(files, function(f) {
    message(sprintf("Processing: %s", basename(f)))
    tryCatch(plot_file(f, global_hits_max = global_hits_max), error = function(e) {
      warning(sprintf("Failed on %s: %s", basename(f), e$message))
      NULL
    })
  })

  n_ok <- sum(vapply(res, function(x) !is.null(x), logical(1)))
  message(sprintf("Completed %d/%d plot(s).", n_ok, length(files)))
}

if (identical(environment(), globalenv())) {
  main()
}
