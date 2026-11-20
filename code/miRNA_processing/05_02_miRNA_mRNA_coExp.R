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


mRNA_EXPES <- readRDS("output/DEGenes_bulk/RESIDUAL_gene_expression.RDS")
mRNA_EXPES <- mRNA_EXPES[, as.character(uni_covs$Sample_ID)]

mRNA_AD_exp <- mRNA_EXPES[, as.character(uni_AD_set$Sample_ID), drop = FALSE]
mRNA_nAD_exp <- mRNA_EXPES[, as.character(uni_nAD_set$Sample_ID), drop = FALSE]
mRNA_all_exp <- mRNA_EXPES[, as.character(uni_all.ad$Sample_ID), drop = FALSE]






resilience_miRNA <- read.csv("output/deMirs/cognition_miRNA_prop_odds.csv")
resilience_miRNA <- resilience_miRNA[resilience_miRNA$p.adj<0.05,]$gene




resilience_mRNA  <- read.csv("output/DEGenes_bulk/cognition_genes_prop_odds.csv")
resilience_genes  <- resilience_mRNA[resilience_mRNA$p.adj<0.01,]$gene

tarbase_v9 <- fread("data/Homo_sapiens_TarBase-v9.tsv.gz")

tarbase_filter <- tarbase_v9 %>% dplyr::filter(confidence == 1, 
                                     !is.na(microt_score),
                                     interaction_group == "primary",
                                     regulation == "Negative",
                                     microt_score > 0.5,
                                     species == "Homo sapiens",
                                     gene_location == "3UTR"
                                     ) %>%
                                dplyr::select(mirna_name,gene_name,gene_id)
tarbase_filter <- unique(tarbase_filter)
tarbase_filter <- tarbase_filter %>% 
            dplyr::filter(.,mirna_name %in% resilience_miRNA,
            gene_id %in% rownames(mRNA_EXPES))

sort(table(tarbase_filter$gene_name))



# Barplot of genes with >= 5 miRNAs (from tarbase_filter$gene_name)
gene_counts <- tarbase_filter %>%
  dplyr::count(gene_name, name = "miRNA_count") %>%
  dplyr::arrange(desc(miRNA_count))

genes_ge5 <- gene_counts %>% dplyr::filter(miRNA_count >= 5)

if (nrow(genes_ge5) == 0) {
  warning("No genes with >= 5 miRNAs found in tarbase_filter")
} else {
  # Save counts table
  readr::write_csv(genes_ge5, file.path(coexp.dir, "/tarbase_gene_miRNA_counts_ge5.csv"))
  
  # Create horizontal bar plot ordered by count
  p <- ggplot2::ggplot(genes_ge5, ggplot2::aes(y = reorder(gene_name, miRNA_count), x = miRNA_count)) +
    coord_fixed(ratio = 1) +
    ggplot2::geom_col(fill = "gray30") +
    #ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Genes with >= 5 miRNAs (TarBase filter)",
      x = NULL,
      y = "Number of targeting miRNAs"
    ) +
    ggplot2::theme_minimal(base_size = 12) 
  
  # Print to current device and save files
  print(p)
  ggplot2::ggsave(filename = file.path(coexp.dir, "/tarbase_gene_miRNA_counts_ge5.png"), plot = p, width = 8, height = max(4, 0.25 * nrow(genes_ge5)))
  ggplot2::ggsave(filename = file.path(coexp.dir, "/tarbase_gene_miRNA_counts_ge5.pdf"), plot = p, width = 8, height = max(4, 0.25 * nrow(genes_ge5)))
}

###################################
# miRNA–mRNA correlation (AD vs nAD; TarBase) #
###############################################

cat("\nComputing miRNA–mRNA Spearman correlations with p-values ...\n")

# Helper: compute Spearman correlation for specified miRNA–gene pairs
correlation_pairs <- function(mi_mat, mrna_mat, pairs_df) {
  mi_idx  <- match(pairs_df$mirna_name, rownames(mi_mat))
  ge_idx  <- match(pairs_df$gene_id,   rownames(mrna_mat))
  ok <- !is.na(mi_idx) & !is.na(ge_idx)
  vals <- rep(NA_real_, nrow(pairs_df))
  ns   <- rep(NA_integer_, nrow(pairs_df))
  pvs  <- rep(NA_real_, nrow(pairs_df))
  if (any(ok)) {
    res <- mapply(
      function(i_mi, i_ge) {
        x <- as.numeric(mi_mat[i_mi, ])
        y <- as.numeric(mrna_mat[i_ge, ])
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
    gene_id    = pairs_df$gene_id,
    gene_name  = pairs_df$gene_name,
    r          = vals,
    n          = ns,
    p          = pvs
  )
}

# Limit to TarBase miRNAs/genes to reduce compute and ensure name matching
tb_pairs <- tarbase_filter %>%
  dplyr::select(mirna_name, gene_id,gene_name) %>%
  distinct(mirna_name,gene_id,,.keep_all = TRUE)

keep_mirs <- intersect(rownames(miRNA_all_exp), unique(tb_pairs$mirna_name))
keep_genes <- intersect(rownames(mRNA_EXPES),     unique(tb_pairs$gene_id))

# Subset expression to retained features
miRNA_AD_sub  <- miRNA_AD_exp[rownames(miRNA_AD_exp) %in% keep_mirs, , drop = FALSE]
miRNA_nAD_sub <- miRNA_nAD_exp[rownames(miRNA_nAD_exp) %in% keep_mirs, , drop = FALSE]

mRNA_AD_sub   <- mRNA_AD_exp[rownames(mRNA_AD_exp) %in% keep_genes, , drop = FALSE]
mRNA_nAD_sub  <- mRNA_nAD_exp[rownames(mRNA_nAD_exp) %in% keep_genes, , drop = FALSE]

# Sanity: columns (samples) should correspond positionally (same subject order)
stopifnot(ncol(miRNA_AD_sub)  == ncol(mRNA_AD_sub))
stopifnot(ncol(miRNA_nAD_sub) == ncol(mRNA_nAD_sub))

tb_pairs_unique <- tb_pairs %>% distinct(mirna_name, gene_id,.keep_all = TRUE)

# Compute Spearman correlations only for TarBase-mentioned pairs
ad_pairs  <- correlation_pairs(miRNA_AD_sub,  mRNA_AD_sub,  tb_pairs_unique) %>% mutate(p_adj = p.adjust(p, method = "BH"))
nad_pairs <- correlation_pairs(miRNA_nAD_sub, mRNA_nAD_sub, tb_pairs_unique) %>% mutate(p_adj = p.adjust(p, method = "BH"))



mir_mrna_cor <- dplyr::inner_join(ad_pairs,nad_pairs, by= c("mirna_name",'gene_id',"gene_name"),suffix = c("_AD","_nAD"))




# diff_cor_test <- function(r1, m1, r2, m2) {
#   z1 <- atanh(r1); z2 <- atanh(r2)
#   se <- sqrt(1/(m1 - 3) + 1/(m2 - 3))
#   Z  <- (z1 - z2)/se
#   p  <- 2 * pnorm(-abs(Z))
#   c(Z = Z, p = p)
# }

# consensus_cor <- function(r1, m1, r2, m2, k = 0) {
#   z1 <- atanh(r1); z2 <- atanh(r2)
#   w1 <- m1 - k - 3; w2 <- m2 - k - 3
#   stopifnot(w1 > 0, w2 > 0)
#   zhat <- (w1*z1 + w2*z2) / (w1 + w2)
#   se   <- 1 / sqrt(w1 + w2)
#   Z    <- zhat / se
#   p    <- 2*pnorm(-abs(Z))
#   ci_z <- zhat + c(-1,1) * qnorm(0.975) * se
#   list(
#     r_hat = tanh(zhat),
#     ci = tanh(ci_z),
#     Z = Z, p = p,
#     #Cochran's Q statistic for heterogeneity
#     Q = w1*(z1-zhat)^2 + w2*(z2-zhat)^2,
#     w = c(w1=w1, w2=w2)
#   )
# }


# consensus_results <- mapply(
#   consensus_cor,
#   r1 = mir_mrna_cor$r_AD,
#   m1 = mir_mrna_cor$n_AD,
#   r2 = mir_mrna_cor$r_nAD,
#   m2 = mir_mrna_cor$n_nAD,
#   SIMPLIFY = FALSE
# )

# consensus_results_df <- do.call(rbind, lapply(consensus_results, function(res) {
#   data.frame(
#     r_hat = res$r_hat,
#     ci_lower = res$ci[1],
#     ci_upper = res$ci[2],
#     Z_consensus = res$Z,
#     p_consensus = res$p,
#     Q = res$Q,
#     stringsAsFactors = FALSE
#   )
# }))

# consensus_results_df$p_adj_consensus <- p.adjust(consensus_results_df$p_consensus, method = "BH")


consensus_results_df <- make_consensus (
  r1 = mir_mrna_cor$r_AD,
  m1 = mir_mrna_cor$n_AD,
  r2 = mir_mrna_cor$r_nAD,
  m2 = mir_mrna_cor$n_nAD,
  k = 0
)

diff_results_df <- make_diff_cor (
  r1 = mir_mrna_cor$r_AD,
  m1 = mir_mrna_cor$n_AD,
  r2 = mir_mrna_cor$r_nAD,
  m2 = mir_mrna_cor$n_nAD
)

all_res <- bind_cols(mir_mrna_cor,consensus_results_df,diff_results_df)


# Write outputs
out_mm_dir <- file.path(coexp.dir, "mirna_mrna")
if (!dir.exists(out_mm_dir)) dir.create(out_mm_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(all_res,
                 file.path(out_mm_dir, "miRNA_mRNA_cor_by_group_tarbase.csv"))

cat("miRNA–mRNA correlations written:",
    file.path(out_mm_dir, "miRNA_mRNA_cor_by_group_tarbase.csv"), "\n")


all_res0 <- all_res %>% filter((p_adj_AD<0.05 & r_AD < 0) | (p_adj_nAD<0.05 & r_nAD < 0) | (p_adj_consensus < 0.05 & r_hat < 0) ) 
table(all_res0$mirna_name)

readr::write_csv(all_res0,
                 file.path(out_mm_dir, "miRNA_mRNA_negative_cor_by_group_tarbase.csv"))


all_res <- all_res %>% filter((p_adj_AD<0.05 | p_adj_nAD<0.05 | p_adj_consensus < 0.05) & gene_name %in% resilience_genes ) 
 


resilience_mRNA  <- read.csv("output/DEGenes_bulk/cognition_genes_prop_odds.csv")
resilience_mRNA  <- resilience_mRNA[resilience_mRNA$p.adj<0.01,c("gene","Dir")]
names(resilience_mRNA) <- c("gene_name","Dir_gene")


resilience_miRNA <- read.csv("output/deMirs/cognition_miRNA_prop_odds.csv")
resilience_miRNA <- resilience_miRNA[resilience_miRNA$p.adj<0.05,c("gene","Dir")]
names(resilience_miRNA) <- c("mirna_name","Dir_miRNA")

resilience_miRNA <- (resilience_miRNA[!duplicated(resilience_miRNA$mirna_name),])




all_res <- left_join(all_res,resilience_mRNA,by="gene_name")
all_res <- left_join(all_res,resilience_miRNA,by="mirna_name")


readr::write_csv(all_res,
                 file.path(out_mm_dir, "miRNA_mRNA_cor_by_group_resilience.csv"))


