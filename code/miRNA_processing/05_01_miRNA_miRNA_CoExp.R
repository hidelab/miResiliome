rm(list = ls())

args = commandArgs(trailingOnly=TRUE)
if(is.na(args[1])){
    setup <- ""
    warning("The experimental design specifications must be set.")
} else {
    setup <- args[1]
}

print("miRNA-pathology correlation analysis")

options(StringAsFactors = F)

library(limma)
library(dplyr)
library(psych)
library(tidyverse)
library(CovariateAnalysis)
library(clusterProfiler)
library(org.Hs.eg.db)
library(broom)
library(VGAM)
library(foreach)
library(ggrepel)

doParallel::registerDoParallel(10)

source("code/00-02-functions2.R")






data.dir        <- "data/preprocessed/"
output.dir      <- paste0("Output", setup,"/")
input.dir       <- paste0(output.dir,"preprocessed/")


if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
coexp.dir <- file.path(output.dir, "coexpr")
if (!dir.exists(coexp.dir)) dir.create(coexp.dir, recursive = TRUE, showWarnings = FALSE)



# Loading covariates, corrected miRNA, and corrected Pathways
covs   <- readRDS(file = "data/preprocessed/mirClinical.RDS")
EXPRS  <-  readRDS("data/preprocessed/ResidualmiRNA.RDS")


DEmiRs  <- read.csv("output/deMirs/deMirsADvsResilient.csv")
DEmiRs  <- DEmiRs[DEmiRs$adj.P.Val < 0.05, ]$X
cogmiRs <- read.csv("output/deMirs/cognition_miRNA_prop_odds.csv")
cogmiRs <- cogmiRs[cogmiRs$p.adj < 0.05, ]$gene


DEmiRs <- union(DEmiRs, cogmiRs)






EXPRS <- EXPRS[rownames(EXPRS) %in% DEmiRs, ]

# Scale mirna expression
EXPRS <- t(scale(t(EXPRS), center = TRUE, scale = TRUE))  # columns: miRNAs


mad_vals <- apply(EXPRS, 1, mad)
keep_var <- mad_vals > 0
sum(keep_var)
Y_resid <- EXPRS[keep_var,, drop=FALSE]

hist(mad_vals, breaks = 30, main = "miRNA MAD values", xlab = "MAD")
dim(Y_resid)


all.ad  <- covs[covs$cogdx %in% c(1,2,4),]
gene.AD <- EXPRS[,as.character(all.ad$mirna_id)]



AD_set  <- covs[covs$cogdx %in% c(4),]
AD_exp  <- Y_resid[,as.character(AD_set$mirna_id)]

nAD_set <- covs[covs$cogdx %in% c(1,2),]
nAD_exp <- Y_resid[,as.character(nAD_set$mirna_id)]


dim(AD_exp)
dim(nAD_exp)



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





res_ad <- get_partial_cor(AD_exp)
res_nad <- get_partial_cor(nAD_exp)

mir_mir_cor <- dplyr::inner_join(res_ad,res_nad, by= c("var1",'var2'),suffix = c("_AD","_nAD"))

consensus_results_df <- make_consensus (
  r1 = mir_mir_cor$r_AD,
  m1 = mir_mir_cor$df_AD,
  r2 = mir_mir_cor$r_nAD,
  m2 = mir_mir_cor$df_nAD,
  k = 0
)

diff_results_df <- make_diff_cor (
  r1 = mir_mir_cor$r_AD,
  m1 = mir_mir_cor$df_AD,
  r2 = mir_mir_cor$r_nAD,
  m2 = mir_mir_cor$df_nAD
)


all_res <- bind_cols(mir_mir_cor,consensus_results_df,diff_results_df)

  

# Write outputs
out_mm_dir <- file.path(coexp.dir, "mirna_mirna")
if (!dir.exists(out_mm_dir)) dir.create(out_mm_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(all_res,
                 file.path(out_mm_dir, "miRNA_miRNA_cor_by_group_resilience.csv"))

cat("miRNA–miRNA correlations written:",
    file.path(out_mm_dir, "miRNA_moRNA_cor_by_group.csv"), "\n")




all(res_ad$var1 == res_nad$var1 &
    res_ad$var2 == res_nad$var2)




