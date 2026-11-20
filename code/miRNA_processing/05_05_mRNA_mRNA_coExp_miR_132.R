rm(list = ls())

# Usage:
# Rscript code/_miRNA_processing/05_03_miRNA_coExp.R 11 [pearson|spearman] [0.05] [0.3]
#   args[1] = setup number used for Output<setup>/
#   args[2] = correlation method (default: pearson)
#   args[3] = FDR (BH) cutoff for edge significance (default: 0.05)
#   args[4] = minimum absolute correlation cutoff (default: 0.30)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || is.na(args[1])) {
  setup <- ""
  warning("The experimental design specifications must be set. Defaulting setup=11.")
} else {
  setup <- args[1]
}

method <- ifelse(length(args) >= 2 && !is.na(args[2]), args[2], "spearman")
fdr_cutoff <- ifelse(length(args) >= 3 && !is.na(as.numeric(args[3])), as.numeric(args[3]), 0.05)
min_abs_cor <- ifelse(length(args) >= 4 && !is.na(as.numeric(args[4])), as.numeric(args[4]), sqrt(0.2))  # ~0.316

cat("miRNA co-expression networks by condition\n")
cat("setup:", setup, " method:", method, " FDR:", fdr_cutoff, " min|r|:", min_abs_cor, "\n")

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(Hmisc)      # rcorr for r and p
  library(igraph)
  library(readr)


  library(limma)
  library(psych)
  library(tidyverse)
  library(CovariateAnalysis)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(broom)
  library(data.table)

})

source("code/00-02-functions2.R")

data.dir   <- "data/preprocessed/"
output.dir <- paste0("output", setup, "/")
input.dir  <- paste0(output.dir, "preprocessed/")

if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
coexp.dir <- file.path(output.dir, "coexpr")
if (!dir.exists(coexp.dir)) dir.create(coexp.dir, recursive = TRUE, showWarnings = FALSE)


exp_covs  <- readRDS(file.path(data.dir, "clinicalTotal.RDS"))
exp_covs$SampleID


uni_covs <- exp_covs 





# Partition samples
uni_all.ad  <- uni_covs[uni_covs$cogdx.x %in% c(1, 2, 4), ]            # all relevant
uni_AD_set  <- uni_covs[uni_covs$cogdx.x %in% c(4), ]                   # AD
uni_nAD_set <- uni_covs[uni_covs$cogdx.x %in% c(1, 2), ]                # non-AD


resilience_targets <- read.csv("output/coexpr/mirna_mrna/miRNA_mRNA_cor_by_group_tarbase.csv")
resilience_targets <- resilience_targets %>%
  filter((p_adj_AD < 0.05 & abs(r_AD) >= 0.2), mirna_name =="hsa-miR-132-3p")



resilience_targets$gene_id

resilience_targets <- (unique(resilience_targets[,c("gene_id","gene_name")]))

mRNA_EXPES <- readRDS("output/DEGenes_bulk/RESIDUAL_gene_expression.RDS")

rownames(mRNA_EXPES) 
mRNA_EXPES <- mRNA_EXPES[(resilience_targets$gene_id), as.character(uni_covs$SampleID)]
rownames(mRNA_EXPES) <- resilience_targets$gene_name[match(rownames(mRNA_EXPES),resilience_targets$gene_id)]

AD_exp <- mRNA_EXPES[, as.character(uni_AD_set$SampleID), drop = FALSE]
nAD_exp <- mRNA_EXPES[, as.character(uni_nAD_set$SampleID), drop = FALSE]
all_exp <- mRNA_EXPES[, as.character(uni_all.ad$SampleID), drop = FALSE]



get_partial_cor <- function(exp){


  Z <- t(exp) # exp has miRNAs in rows
  R <- cor(Z)                 # correlation
  K <- solve(R)               # precision
  P <- -cov2cor(K)            # partial correlations
  diag(P) <- 1


  
  n = nrow(Z)
  p = ncol(Z)


  df <- n - p
  if (df <= 0) stop("df = n - p must be > 0 for classical p-values. Consider shrinkage + resampling.")

    ## Extract upper triangle partial correlations
  ut <- which(upper.tri(P), arr.ind = TRUE)
  r  <- P[ut]

  ## Guard against +/-1 due to numerics
  r  <- pmin(pmax(r, -0.999999), 0.999999)

    ## t-test for partial correlation: t = r * sqrt(df / (1 - r^2)),  df = n - p
  tval <- r * sqrt(df / pmax(1e-12, 1 - r^2))
  pval <- 2 * pt(-abs(tval), df = df)
  ## Multiple testing correction (BH/FDR)
  qval <- p.adjust(pval, method = "BH")

   ## Optional: assemble a tidy edge table
  res_ad <- data.frame(
  var1 = colnames(Z)[ut[,1]],
  var2 = colnames(Z)[ut[,2]],
  r = r,
  t = tval,
  pval = pval,
  p_adj = qval,
  df = df
)
}   





#res_ad <- get_partial_cor(AD_exp)
res_nad <- get_partial_cor(nAD_exp)


View(res_nad)





mir_mir_cor <- res_nad


all_res <- bind_cols(mir_mir_cor)

  

# Write outputs
out_mm_dir <- file.path(coexp.dir, "mrna_mrna")
if (!dir.exists(out_mm_dir)) dir.create(out_mm_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(all_res,
                 file.path(out_mm_dir, "mir_132_mRNA_mRNA_cor_by_group_resilience.csv"))

cat("mRNA–mRNA correlations written:",
    file.path(out_mm_dir, "mir_132_mRNA_mRNA_cor_by_group.csv"), "\n")


