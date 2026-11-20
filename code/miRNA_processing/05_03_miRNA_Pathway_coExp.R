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

# Load covariates and expression
mir_covs  <- readRDS(file.path(data.dir, "mirClinical.RDS"))
miR_EXPRS <- readRDS(file.path(data.dir, "ResidualmiRNA.RDS"))

# Scale miRNA expression across samples and drop zero-MAD rows
miR_EXPRS <- t(scale(t(miR_EXPRS), center = TRUE, scale = TRUE))
Y_resid <- miR_EXPRS 





exp_covs  <- readRDS(file.path(data.dir, "clinicalTotal.RDS"))
exp_covs$SampleID

colnames(exp_covs)
colnames(mir_covs)

uni_covs <- (as.character(mir_covs$Sample_ID) %in% exp_covs$SampleID)
uni_covs <- mir_covs[uni_covs, ]






# Partition samples
uni_all.ad  <- uni_covs[uni_covs$cogdx %in% c(1, 2, 4), ]            # all relevant
uni_AD_set  <- uni_covs[uni_covs$cogdx %in% c(4), ]                   # AD
uni_nAD_set <- uni_covs[uni_covs$cogdx %in% c(1, 2), ]                # non-AD

colnames(Y_resid) 

# Align expression to sample IDs in each set
miRNA_AD_exp  <- Y_resid[, as.character(uni_AD_set$mirna_id), drop = FALSE]
miRNA_nAD_exp <- Y_resid[, as.character(uni_nAD_set$mirna_id), drop = FALSE]
miRNA_all_exp <- Y_resid[, as.character(uni_all.ad$mirna_id), drop = FALSE]

cat("miRNA AD matrix dims:", paste(dim(miRNA_AD_exp), collapse = " x "), "\n")
cat("miRNA nAD matrix dims:", paste(dim(miRNA_nAD_exp), collapse = " x "), "\n")
cat("miRNA all matrix dims:", paste(dim(miRNA_all_exp), collapse = " x "), "\n")


path_EXPES <- readRDS("output/DEPaths_bulk_update/ROSMAP_DLPFC_ADvsControl_Residuals.RDS")
path_EXPES <- path_EXPES[, as.character(uni_covs$Sample_ID)]

path_AD_exp <- path_EXPES[, as.character(uni_AD_set$Sample_ID), drop = FALSE]
path_nAD_exp <- path_EXPES[, as.character(uni_nAD_set$Sample_ID), drop = FALSE]
path_all_exp <- path_EXPES[, as.character(uni_all.ad$Sample_ID), drop = FALSE]






resilience_miRNA <- read.csv("output/deMirs/cognition_miRNA_prop_odds.csv")
resilience_miRNA <- resilience_miRNA[resilience_miRNA$p.adj<0.05,]$gene




my_stars <- read.csv("output/deMirs/Heatmap_deMirs_sexSpec_inputData.csv")
my_stars <- unique(my_stars[my_stars$star=="***",]$miRNA)
resilience_miRNA <- union(resilience_miRNA,my_stars)

write.csv(resilience_miRNA,
          file = file.path("output/deMirs/", "significant_miRNAs_for_table1.csv"),
          row.names = FALSE)

resilience_pathways <- read.csv("output/Bulk_paths_prop_odds/cog_paths.csv")
# resilience_pathways <- resilience_pathways[resilience_pathways$p.adj<0.01,]$gene
resilience_pathways <- resilience_pathways$gene


tb_pairs <- expand.grid(mirna_name = resilience_miRNA,pathway_name = resilience_pathways)



###################################
# miRNA–path correlation (AD vs nAD; TarBase) #
###############################################

cat("\nComputing miRNA–path Spearman correlations with p-values ...\n")

# Helper: compute Spearman correlation for specified miRNA–pathway pairs
correlation_pairs <- function(mi_mat, path_mat, pairs_df) {
  mi_idx  <- match(pairs_df$mirna_name, rownames(mi_mat))
  ge_idx  <- match(pairs_df$pathway_name,   rownames(path_mat))
  ok <- !is.na(mi_idx) & !is.na(ge_idx)
  vals <- rep(NA_real_, nrow(pairs_df))
  ns   <- rep(NA_integer_, nrow(pairs_df))
  pvs  <- rep(NA_real_, nrow(pairs_df))
  if (any(ok)) {
    res <- mapply(
      function(i_mi, i_ge) {
        x <- as.numeric(mi_mat[i_mi, ])
        y <- as.numeric(path_mat[i_ge, ])
        cc <- is.finite(x) & is.finite(y)
        n  <- sum(cc)
        if (n < 3) return(c(r = NA_real_, n = n, p = NA_real_))
        r <- suppressWarnings(cor(x[cc], y[cc], method = "pearson"))
        if (is.na(r) || abs(r) >= 1) {
          p <- NA_real_
        } else {
          tstat <- r * sqrt((n - 2) / pmax(1e-12, 1 - r^2))
          p <- 2 * stats::pt(-abs(tstat), df = n - 2)
        }
        c(r = r, n = n, p = p)
      },
      mi_idx[ok], ge_idx[ok]
    )
    vals[ok] <- res["r", ]
    ns[ok]   <- res["n", ]
    pvs[ok]  <- res["p", ]
  }
  tibble(
    mirna_name = pairs_df$mirna_name,
    pathway_name    = pairs_df$pathway_name,
    r          = vals,
    n          = ns,
    p          = pvs
  )
}

# Limit to TarBase miRNAs/genes to reduce compute and ensure name matching
keep_mirs <- intersect(rownames(miRNA_all_exp), unique(tb_pairs$mirna_name))
keep_genes <- intersect(rownames(path_EXPES),     unique(tb_pairs$pathway_name))

# Subset expression to retained features
miRNA_AD_sub  <- miRNA_AD_exp[rownames(miRNA_AD_exp) %in% keep_mirs, , drop = FALSE]
miRNA_nAD_sub <- miRNA_nAD_exp[rownames(miRNA_nAD_exp) %in% keep_mirs, , drop = FALSE]

path_AD_sub   <- path_AD_exp[rownames(path_AD_exp) %in% keep_genes, , drop = FALSE]
path_nAD_sub  <- path_nAD_exp[rownames(path_nAD_exp) %in% keep_genes, , drop = FALSE]

# Sanity: columns (samples) should correspond positionally (same subject order)
stopifnot(ncol(miRNA_AD_sub)  == ncol(path_AD_sub))
stopifnot(ncol(miRNA_nAD_sub) == ncol(path_nAD_sub))

tb_pairs_unique <- tb_pairs %>% distinct(mirna_name, pathway_name,.keep_all = TRUE)

# Compute Spearman correlations only for TarBase-mentioned pairs
ad_pairs  <- correlation_pairs(miRNA_AD_sub,  path_AD_sub,  tb_pairs_unique) %>% mutate(p_adj = p.adjust(p, method = "BH"))
nad_pairs <- correlation_pairs(miRNA_nAD_sub, path_nAD_sub, tb_pairs_unique) %>% mutate(p_adj = p.adjust(p, method = "BH"))



mir_path_cor <- dplyr::inner_join(ad_pairs,nad_pairs, by= c("mirna_name",'pathway_name'),suffix = c("_AD","_nAD"))

# mir_path_cor <- bind_rows(ad_pairs, nad_pairs) %>%
#   group_by(group) %>%
#   mutate(p_adj = p.adjust(p, method = "BH")) %>%
#   ungroup()


consensus_results_df <- make_consensus (
  r1 = mir_path_cor$r_AD,
  m1 = mir_path_cor$n_AD,
  r2 = mir_path_cor$r_nAD,
  m2 = mir_path_cor$n_nAD,
  k = 0
)
# Test for differential correlation
diff_results_df <- make_diff_cor (
  r1 = mir_path_cor$r_AD,
  m1 = mir_path_cor$n_AD,
  r2 = mir_path_cor$r_nAD,
  m2 = mir_path_cor$n_nAD
)


all_res <- bind_cols(mir_path_cor,consensus_results_df,diff_results_df)


# Write outputs
out_mm_dir <- file.path(coexp.dir, "mirna_path")
if (!dir.exists(out_mm_dir)) dir.create(out_mm_dir, recursive = TRUE, showWarnings = FALSE)






readr::write_csv(all_res,
                 file.path(out_mm_dir, "miRNA_path_cor_by_group_tarbase.csv"))

cat("miRNA–path correlations written:",
    file.path(out_mm_dir, "miRNA_path_cor_by_group_tarbase_spearman.csv"), "\n")


all_res0 <- all_res %>% filter((p_adj_AD<0.05 | p_adj_nAD<0.05 | p_adj_consensus < 0.05)) 




readr::write_csv(all_res0,
                 file.path(out_mm_dir, "sig_miRNA_path_cor_by_group_tarbase.csv"))

cat("Significant miRNA–path correlations written:",
    file.path(out_mm_dir, "miRNA_path_cor_by_group_tarbase_spearman.csv"), "\n")



resilience_pathways <- read.csv("output/Bulk_paths_prop_odds/cog_paths.csv")
resilience_pathways <- resilience_pathways[resilience_pathways$p.adj<0.01,]$gene

all_res <- all_res %>% filter((p_adj_AD<0.05 | p_adj_nAD<0.05 | p_adj_consensus < 0.05) & pathway_name %in% resilience_pathways ) 
 


resilience_pathways <- read.csv("output/Bulk_paths_prop_odds/cog_paths.csv")
resilience_pathways  <- resilience_pathways[resilience_pathways$p.adj<0.01,c("gene","Dir")]
names(resilience_pathways) <- c("pathway_name","Dir_pathway")


resilience_miRNA <- read.csv("output/deMirs/cognition_miRNA_prop_odds.csv")
resilience_miRNA <- resilience_miRNA[resilience_miRNA$p.adj<0.05,c("gene","Dir")]
names(resilience_miRNA) <- c("mirna_name","Dir_miRNA")


all_res <- left_join(all_res,resilience_pathways,by="pathway_name")
all_res <- left_join(all_res,resilience_miRNA,by="mirna_name")


readr::write_csv(all_res,
                 file.path(out_mm_dir, "miRNA_pathway_cor_by_group_resilience.csv"))







