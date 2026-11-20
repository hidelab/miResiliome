#!/usr/bin/env Rscript

# Boxplots of residual miRNA expression by Sex within Diagnosis groups
# - Inputs:
#   - data/preprocessed/ResidualmiRNA.RDS (matrix: rows=miRNAs, cols=samples)
#   - data/preprocessed/mirClinical.RDS (clinical with msex, cogdx, braaksc, ceradsc, age_death, mirna_id)
# - For each selected miRNA (default: hsa-miR-7-5p), create boxplots across
#   Diagnosis = Control, Resilient, AD; colored by Sex (Male/Female), using residual expression.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
})

miRNAs <-  c("hsa-miR-7-5p","hsa-miR-218-5p","hsa-miR-132-3p")


data_dir <- "data/preprocessed"
expr_path <- file.path(data_dir, "ResidualmiRNA.RDS")
clin_path <- file.path(data_dir, "mirClinical.RDS")

if (!file.exists(expr_path)) stop("Missing ", expr_path)
if (!file.exists(clin_path)) stop("Missing ", clin_path)

EXPRS <- readRDS(expr_path)   # rows: miRNAs, cols: samples (mirna_id)
CLIN  <- readRDS(clin_path)

# Derive Diagnosis groups consistent with 04-02-deMirs.R
cogdx.resils <- c(1)
cogdx.AD     <- c(4)
cogdx.normal <- c(1)

braak.resils <- c(3,4,5,6)
braak.AD     <- c(3,4,5,6)

plq.resils   <- c(1,2)
plq.AD       <- c(1,2)

age.cut      <- 80
age.cut.AD   <- 80
age.cut.norm <- 0

CLIN <- mutate(CLIN, RESIL = "OTHER")

resilients <- CLIN %>%
  filter(cogdx %in% cogdx.resils,
         braaksc %in% braak.resils,
         ceradsc %in% plq.resils,
         age_death >= age.cut) %>%
  mutate(RESIL = "Resilient")

AD <- CLIN %>%
  filter(cogdx %in% cogdx.AD,
         braaksc %in% braak.AD,
         ceradsc %in% plq.AD,
         age_death >= age.cut.AD) %>%
  mutate(RESIL = "AD")

normals <- CLIN %>%
  filter(cogdx %in% cogdx.normal,
         braaksc %in% c(0,1,2),
         ceradsc %in% c(3,4),
         age_death >= age.cut.norm) %>%
  mutate(RESIL = "Control")

CLIN$RESIL[CLIN$mirna_id %in% resilients$mirna_id] <- "Resilient"
CLIN$RESIL[CLIN$mirna_id %in% AD$mirna_id]        <- "AD"
CLIN$RESIL[CLIN$mirna_id %in% normals$mirna_id]   <- "Control"

CLIN$Diagnosis <- factor(CLIN$RESIL, levels = c("Control", "Resilient", "AD"))

# Harmonize Sex labels
if (!"msex" %in% names(CLIN)) {
  stop("Clinical file missing 'msex' column.")
}
CLIN$Sex <- as.character(CLIN$msex)
CLIN$Sex[CLIN$Sex %in% c("male", "Male", "M", 1)]   <- "Male"
CLIN$Sex[CLIN$Sex %in% c("female", "Female", "F", 0)] <- "Female"
CLIN$Sex <- factor(CLIN$Sex, levels = c("Female", "Male"))

# Align samples present in both clinical and expression
common_samples <- intersect(colnames(EXPRS), CLIN$mirna_id)
if (length(common_samples) == 0) stop("No overlapping sample IDs between expression and clinical.")

CLIN_sub <- CLIN[match(common_samples, CLIN$mirna_id), , drop = FALSE]
stopifnot(identical(CLIN_sub$mirna_id, common_samples))

# Output directory alongside deMirs results
out_dir <- "output/deMirs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load significance data from heatmap pipeline
sig_data_path <- "output/deMirs/Heatmap_deMirs_sexSpec_inputData.csv"
if (!file.exists(sig_data_path)) stop("Missing ", sig_data_path, ". Run 04-10-Heatmap_SexSpecific.R first.")
sig_all <- read.csv(sig_data_path, stringsAsFactors = FALSE) %>%
  filter(!is.na(adj.P.Val) & adj.P.Val < 0.2)

# Dodge width matches the boxplot
dw <- 0.8
# x positions: Female=1, Male=2; within each, 3 Diagnosis levels dodged
# Diagnosis factor order: Control=1, Resilient=2, AD=3 → dodge offsets: -dw/3, 0, +dw/3
diag_offset <- c(Control = -dw/3, Resilient = 0, AD = dw/3)
sex_x <- c(Female = 1, Male = 2)

# Map contrast names to (xmin, xmax) coordinates
xf <- unname(sex_x["Female"])
xm <- unname(sex_x["Male"])
dc <- unname(diag_offset["Control"])
dr <- unname(diag_offset["Resilient"])
da <- unname(diag_offset["AD"])

contrast_coords <- list(
  AD_vs_Ctrl_in_Females      = c(xmin = xf + dc, xmax = xf + da),
  AD_vs_Resilient_in_Females = c(xmin = xf + dr, xmax = xf + da),
  AD_vs_Ctrl_in_Males        = c(xmin = xm + dc, xmax = xm + da),
  AD_vs_Resilient_in_Males   = c(xmin = xm + dr, xmax = xm + da),
  Male_vs_Female_in_AD       = c(xmin = xf + da, xmax = xm + da),
  Male_vs_Female_in_Control  = c(xmin = xf + dc, xmax = xm + dc),
  Male_vs_Female_in_Resilient = c(xmin = xf + dr, xmax = xm + dr)
)

plot_one <- function(mir) {
  if (!(mir %in% rownames(EXPRS))) {
    warning("miRNA not found in residual expression: ", mir)
    return(invisible(NULL))
  }
  y <- as.numeric(EXPRS[mir, common_samples])
  df <- CLIN_sub %>%
    transmute(Sample = mirna_id, Diagnosis, Sex, Residual = y)

  # Keep only target diagnosis groups
  df <- df %>% filter(!is.na(Diagnosis), Diagnosis %in% levels(CLIN$Diagnosis))

  # Build significance bars for this miRNA
  mir_sig <- sig_all %>% filter(miRNA == mir, contrast %in% names(contrast_coords))
  sig_bars <- NULL
  if (nrow(mir_sig) > 0) {
    sig_bars <- data.frame(
      label = ifelse(mir_sig$adj.P.Val < 0.05, "***", "^"),
      xmin  = sapply(mir_sig$contrast, function(ct) contrast_coords[[ct]]["xmin"]),
      xmax  = sapply(mir_sig$contrast, function(ct) contrast_coords[[ct]]["xmax"]),
      stringsAsFactors = FALSE
    )
    sig_bars <- sig_bars[order(abs(sig_bars$xmax - sig_bars$xmin)), ]
    y_max <- max(df$Residual, na.rm = TRUE)
    y_range <- diff(range(df$Residual, na.rm = TRUE))
    sig_bars$y_position <- y_max + y_range * seq(0.08, by = 0.08, length.out = nrow(sig_bars))
  }

  p <- ggplot(df, aes(x = Sex, y = Residual, fill = Diagnosis)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(width = 0.8)) +
    geom_jitter(
                #aes(color = "black"),
                #width = 0.15,
                #height = 0,
                size = 1,
                position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), alpha = 0.2) +
    scale_fill_brewer(palette = "Set2") +
    #scale_color_brewer(palette = "Set2", guide = "none") +
    labs(title = paste0(mir, ""),
         x = NULL, y = "Normalized expression") +
    theme_bw(base_size = 12) +
    theme(legend.position = "top",
          plot.title = element_text(hjust = 0.5))

  if (!is.null(sig_bars) && nrow(sig_bars) > 0) {
    tip <- y_range * 0.015
    for (k in seq_len(nrow(sig_bars))) {
      xmn <- sig_bars$xmin[k]
      xmx <- sig_bars$xmax[k]
      yp  <- sig_bars$y_position[k]
      lbl <- sig_bars$label[k]
      p <- p +
        annotate("segment", x = xmn, xend = xmn, y = yp - tip, yend = yp) +
        annotate("segment", x = xmn, xend = xmx, y = yp, yend = yp) +
        annotate("segment", x = xmx, xend = xmx, y = yp, yend = yp - tip) +
        annotate("text", x = (xmn + xmx) / 2, y = yp, label = lbl,
                 vjust = -0.3, size = 4, fontface = "bold")
    }
    p <- p + coord_cartesian(ylim = c(NA, max(sig_bars$y_position) + y_range * 0.08), clip = "off")
  }

  pdf_path <- file.path(out_dir, paste0("Boxplot_residuals_", gsub("/", "-", mir), ".pdf"))
  png_path <- file.path(out_dir, paste0("Boxplot_residuals_", gsub("/", "-", mir), ".jpg"))
  

  ggsave(pdf_path, p, width = 3.5, height = 5, device = cairo_pdf)
  ggsave(png_path, p, width = 3.5, height = 5, dpi = 600)

  data_path <- file.path(out_dir, paste0("Boxplot_residuals_data_", gsub("/", "-", mir), ".csv"))

  write.csv(df,data_path)

  message("Saved: ", pdf_path, " and ", png_path)
}

invisible(lapply(miRNAs, plot_one))

message("Done.")

