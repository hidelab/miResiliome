#!/usr/bin/env Rscript

# Build ComplexHeatmap for sex-specific de miRNAs
# - Reads all files matching output*/deMirs/deMirs_sexSpec_*.csv
# - Selects union of miRNAs with adj.P.Val < 0.2
# - Heatmap cells show signed -log10(P.Value) (sign from logFC)
# - Stars overlay: * if adj.P.Val < 0.2, ** if < 0.05
# - Column names formatted as requested, e.g., (Male) (AD vs Control)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# Locate all contrast result files
files_sex <- Sys.glob("output*/deMirs/deMirs_sexSpec_*.csv")
# Also include cohort-wide contrasts (All)
file_all_ctrl_candidates <- c(
  Sys.glob("output*/deMirs/deMirsADvsControl.csv"),
  Sys.glob("output*/deMirs/deMirsADvControl.csv")
)
file_all_ctrl <- unique(file_all_ctrl_candidates)
file_all_ctrl <- file_all_ctrl[file.exists(file_all_ctrl)]
file_all_res  <- Sys.glob("output*/deMirs/deMirsADvsResilient.csv")
files_props   <- Sys.glob("output/deMirs/*prop*.csv")  # not used here

files <- c(files_sex, file_all_ctrl, file_all_res)
if (length(files) == 0) {
  stop("No files found. Expected sex-specific 'deMirs_sexSpec_*.csv' and/or 'deMirsADvsControl.csv'/'deMirsADvControl.csv', 'deMirsADvsResilient.csv'. Run 04-02-deMirs.R first.")
}


# Map internal contrast names to display names
pretty_name <- function(x) {
  switch(
    x,
    "AD_vs_Ctrl_in_Males"              = "(Male) (AD vs Ctrl)",
    "AD_vs_Ctrl_in_Females"            = "(Female) (AD vs Ctrl)",
    "AD_vs_Resilient_in_Males"         = "(Male) (AD vs Res)",
    "AD_vs_Resilient_in_Females"       = "(Female) (AD vs Res)",
    "Male_vs_Female_in_AD"             = "(M vs F) (AD)",
    "Male_vs_Female_in_Control"        = "(M vs F) (Ctrl)",
    "Male_vs_Female_in_Resilient"      = "(M vs F) (Res)",
    "Male_vs_Female_in_AD_vs_Control"  = "(M vs F) (AD vs Ctrl)",
    "Male_vs_Female_in_AD_vs_Res"      = "(M vs F) (AD vs Res)",
    "All_AD_vs_Control"                = "(All) (AD vs Ctrl)",
    "All_AD_vs_Resilient"              = "(All) (AD vs Res)",
    "cognition"                        = "(All) (Cognitive loss)",
    "plaque"                           = "(All) (Plaque burden)",
    "tangle"                          = "(All) (Tangle burden)",
    x
  )
}

# Desired column order if available
desired_order <- c(
  "All_AD_vs_Control",
  "All_AD_vs_Resilient",
  "cognition",
  "plaque",
  "tangle",
  "AD_vs_Ctrl_in_Males",
  "AD_vs_Ctrl_in_Females",
  "AD_vs_Resilient_in_Males",
  "AD_vs_Resilient_in_Females",
  "Male_vs_Female_in_Control",
  "Male_vs_Female_in_Resilient",
  "Male_vs_Female_in_AD",
  "Male_vs_Female_in_AD_vs_Control",
  "Male_vs_Female_in_AD_vs_Res"
)

# Read and gather data per contrast
read_contrast <- function(path, name = NULL) {
  # infer internal contrast name from filename
  if (is.null(name)) {
    base <- basename(path)
    # remove prefix and suffix
    name <- base |>
      str_remove("^deMirs_sexSpec_") |>
      str_remove("\\.csv$")
  }
  # read; first column is rownames (miRNA id) in these outputs
  df <- suppressMessages(readr::read_csv(path, show_col_types = FALSE))
  # ensure a 'miRNA' column exists
  if (!("" %in% names(df) || "...1" %in% names(df))) {
    # try reading via base R if readr didn't preserve first unnamed column
    df_base <- read.csv(path, check.names = FALSE)
    if ("" %in% names(df_base)) {
      df <- df_base
    } else {
      # fallback: assume first col is miRNA id
      colnames(df)[1] <- "miRNA"
      df$miRNA <- df[[1]]
    }
  }
  if ("" %in% names(df)) {
    df <- df |> rename(miRNA = "")
  }
  if ("...1" %in% names(df)) {
    df <- df |> rename(miRNA = "...1")
  }
  # keep necessary columns
  keep_cols <- c("miRNA", "logFC", "P.Value", "adj.P.Val")
  missing_cols <- setdiff(keep_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing expected columns in", path, ":", paste(missing_cols, collapse = ", ")))
  }
  df <- df |> select(all_of(keep_cols))
  df$contrast <- name
  df
}

data_list <- list()
for (p in files) {
  nm <- NULL
  bn <- basename(p)
  if (grepl("^deMirs_sexSpec_", bn)) {
    nm <- NULL  # infer from filename
  } else if (bn %in% c("deMirsADvsControl.csv", "deMirsADvControl.csv")) {
    nm <- "All_AD_vs_Control"
  } else if (bn == "deMirsADvsResilient.csv") {
    nm <- "All_AD_vs_Resilient"
  }
  data_list[[length(data_list) + 1]] <- read_contrast(p, name = nm)
}
all_df <- bind_rows(data_list)



prop_list <- list()
for (p in files_props) {
  nm <- NULL
  bn <- basename(p)
  nm <- gsub("_miRNA_prop_odds\\.csv$", "", bn)

  my_tab <- read.csv(p)
  my_tab <- my_tab[,c("gene", "Estimate","pval","p.adj")]
  my_tab$contrast <- nm
  my_tab$Estimate <- -my_tab$Estimate  # invert sign to match logFC directionality
  my_tab <- my_tab |> rename(
    miRNA = gene,
    logFC = Estimate,
    P.Value = pval,
    adj.P.Val = p.adj
  )


  prop_list[[length(prop_list) + 1]] <- my_tab
}
prop_df <- bind_rows(prop_list)

all_df <- bind_rows(all_df, prop_df)




# Determine union of significant miRNAs across all contrasts
sig_miRNAs <- all_df |>
  filter(!is.na(adj.P.Val) & adj.P.Val < 0.2 & P.Value < 0.05) |>
  distinct(miRNA) |>
  pull(miRNA)

if (length(sig_miRNAs) == 0) {
  stop("No miRNAs with adj.P.Val < 0.1 across contrasts. Nothing to plot.")
}

# Subset to significant set for heatmap
all_df_sig <- all_df |>
  filter(miRNA %in% sig_miRNAs)

# Compute signed -log10(P) values and star annotations
small_eps <- 1e-300
all_df_sig <- all_df_sig |>
  mutate(
    signed_logp = sign(logFC) * (-log10(pmax(P.Value, small_eps))),
    star = case_when(
      !is.na(adj.P.Val) & adj.P.Val < 0.05 ~ "***",
      !is.na(adj.P.Val) & adj.P.Val < 0.2  ~ "^",
      TRUE ~ ""
    )
  )

write.csv(all_df_sig, file = "output/deMirs/Heatmap_deMirs_sexSpec_inputData.csv", row.names = FALSE)

# Pivot to matrices
mat_vals <- all_df_sig |>
  select(miRNA, contrast, signed_logp) |>
  tidyr::pivot_wider(names_from = contrast, values_from = signed_logp) |>
  tibble::column_to_rownames("miRNA") |>
  as.matrix()



mat_stars <- all_df_sig |>
  select(miRNA, contrast, star) |>
  tidyr::pivot_wider(names_from = contrast, values_from = star) |>
  tibble::column_to_rownames("miRNA") |>
  as.matrix()

# Reorder columns by desired order where present
present <- colnames(mat_vals)
ordered_cols <- c(desired_order[desired_order %in% present], setdiff(present, desired_order))
mat_vals <- mat_vals[, ordered_cols, drop = FALSE]
mat_stars <- mat_stars[, ordered_cols, drop = FALSE]

# Apply pretty column names
pretty_cols <- vapply(colnames(mat_vals), pretty_name, character(1))
colnames(mat_vals) <- pretty_cols
colnames(mat_stars) <- pretty_cols

# Build color mapping centered at zero
max_abs <- max(abs(mat_vals), na.rm = TRUE)
if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 1
col_fun <- circlize::colorRamp2(c(-max_abs, 0, max_abs), c("#4575B4", "#F7F7F7", "#D73027"))

# Legend for significance symbols
sig_legend <- Legend(
  labels = c("adj.P < 0.2", "adj.P < 0.05"),
  title = "Significance",
  graphics = list(
    function(x, y, w, h) grid.text("^", x, y, gp = gpar(fontsize = 12, fontface = "bold")),
    function(x, y, w, h) grid.text("***", x, y, gp = gpar(fontsize = 12, fontface = "bold"))
  ),
  labels_gp = gpar(fontsize = 10),
  title_gp = gpar(fontsize = 10, fontface = "bold")
)

# Draw heatmap with star overlay
ht <- Heatmap(
  mat_vals,
  name = "signed -log10(p)",
  col = col_fun,
  na_col = "#efefef",
  show_row_dend = FALSE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 12),
  cell_fun = function(j, i, x, y, width, height, fill) {
    lab <- mat_stars[i, j]
    if (!is.na(lab) && nzchar(lab)) {
      grid.text(lab, x = x, y = y, gp = gpar(fontsize = 10, fontface = "bold"))
    }
  },
  ,
   width = ncol(mat_vals)*unit(5, "mm"), 
    height = nrow(mat_vals)*unit(4, "mm")
)

# Ensure an output directory for plots
out_dir <- unique(dirname(files))
if (length(out_dir) == 0) out_dir <- "output/deMirs"
if (!dir.exists(out_dir[1])) dir.create(out_dir[1], recursive = TRUE, showWarnings = FALSE)
out_pdf <- file.path(out_dir[1], "Heatmap_deMirs_sexSpec.pdf")
out_png <- file.path("output/figures_paper", "Figure 2A Heatmap_deMirs_sexSpec.jpg")

message("Saving heatmap to:\n  ", out_pdf, "\n  ", out_png)

pdf(out_pdf, width = 7, height = 13)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", annotation_legend_list = list(sig_legend))
invisible(dev.off())

jpeg(out_png, width = 7, height = 13, unit = "in", res = 600)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", annotation_legend_list = list(sig_legend))
invisible(dev.off())

# Create transposed version with hsa- removed from miRNA names
mat_vals_t <- t(mat_vals)
mat_stars_t <- t(mat_stars)

# Remove "hsa-" prefix from column names (which are miRNA names after transpose)
colnames(mat_vals_t) <- gsub("^hsa-", "", colnames(mat_vals_t))
colnames(mat_stars_t) <- gsub("^hsa-", "", colnames(mat_stars_t))

# Draw transposed heatmap
ht_vertical <- Heatmap(
  mat_vals_t,
  name = "signed -log10(p)",
  col = col_fun,
  na_col = "#efefef",
  show_column_dend = FALSE,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 12),
  cell_fun = function(j, i, x, y, width, height, fill) {
    lab <- mat_stars_t[i, j]
    if (!is.na(lab) && nzchar(lab)) {
      grid.text(lab, x = x, y = y, gp = gpar(fontsize = 9, fontface = "bold"))
    }
  },
  width = ncol(mat_vals_t)*unit(4, "mm"),
  height = nrow(mat_vals_t)*unit(5, "mm")
)

# Save transposed heatmap
out_pdf_v <- file.path(out_dir[1], "Heatmap_deMirs_sexSpec_vertical.pdf")
out_png_v <- file.path(out_dir[1], "Heatmap_deMirs_sexSpec_vertical.png")

message("Saving transposed heatmap to:\n  ", out_pdf_v, "\n  ", out_png_v)

pdf(out_pdf_v, width = 13, height = 7)
draw(ht_vertical, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", annotation_legend_list = list(sig_legend))
invisible(dev.off())

png(out_png_v, width = 1600, height = 1200, res = 300)
draw(ht_vertical, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", annotation_legend_list = list(sig_legend))
invisible(dev.off())

message("Done.")