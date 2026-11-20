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
output.dir      <- paste0("output", setup,"/deMirs/")
if(!dir.exists(output.dir)){
  dir.create(output.dir,recursive = T)
}



input.dir       <- paste0(output.dir,"preprocessed/")



# Loading covariates, corrected miRNA, and corrected Pathways
covs   <- readRDS(file = "data/preprocessed/mirClinical.RDS")
EXPRS  <-  readRDS("data/preprocessed/ResidualmiRNA.RDS")

all.ad  <- covs[covs$cogdx %in% c(1,2,4),]
gene.AD <- EXPRS[,as.character(all.ad$mirna_id)]
#
#write.csv(all.ad,"RUSH_data_request/only_mirs.csv")

####  Generating proportioanl odds model
#*  proportioanl odds model for cognition *# ####

cog.paths <- foreach(j = 1:nrow(gene.AD),.combine = 'rbind') %dopar% {

  tested_gene0    <- rownames(gene.AD)[j]
  tested_gene     <- gene.AD[tested_gene0,]
  this.gene       <- data.frame("tested_gene" = tested_gene,
                               "cog"  = factor(all.ad$cogdx,ordered = T))

  prop.fit        <- proportional.summary("cog","tested_gene ",
                                          this.gene, mir.name = tested_gene0)
}

cog.paths <- clean.write.podds(cog.paths,an.address = paste0(output.dir,"/cognition_miRNA_prop_odds.csv"))

#*  proportioanl odds model for tangles   *# ####

tangle.paths <- foreach(j = 1:nrow(gene.AD),.combine = 'rbind') %dopar% {

  tested_gene0 <- rownames(gene.AD)[j]
  tested_gene <- gene.AD[tested_gene0,]
  this.gene      <- data.frame("tested_gene" = tested_gene,
                               "braak"  = factor(all.ad$braaksc, ordered = T))
  prop.fit       <- proportional.summary("braak","tested_gene",
                                    this.gene, mir.name = tested_gene0)
}

tangle.paths <- clean.write.podds(tangle.paths,paste0(output.dir,"/tangle_miRNA_prop_odds.csv"))

#*  proportioanl odds model for plaques   *# ####


plaque.paths <- foreach(j = 1:nrow(gene.AD),.combine = 'rbind') %dopar% {
  tested_gene0   <- rownames(gene.AD)[j]
  tested_gene    <- gene.AD[tested_gene0,]
  this.gene      <- data.frame("tested_gene" = tested_gene,
                               "cerad"  = factor(4- all.ad$ceradsc,ordered = T))
  prop.fit       <- proportional.summary("cerad","tested_gene",
                                         this.gene, mir.name = tested_gene0)
}

plaque.paths <- clean.write.podds(plaque.paths,paste0(output.dir,"/plaque_miRNA_prop_odds.csv"))




all.paths <- foreach(j = 1:nrow(gene.AD),.combine = 'rbind') %dopar% {
  tested_gene0   <- rownames(gene.AD)[j]
  tested_gene    <- gene.AD[tested_gene0,]
  this.gene      <- data.frame("tested_gene" = tested_gene,
                               "cerad"  = factor(4- all.ad$ceradsc,ordered = T),
                               "braak"  = factor(all.ad$braaksc, ordered = T),
                               "cog"  = factor(all.ad$cogdx, ordered = T))
  prop.fit       <- proportional.summary("cog","tested_gene + cerad + braak",
                                         this.gene, mir.name = tested_gene0)
}

all.paths <- clean.write.podds(all.paths,paste0(output.dir,"/adjusted_cognition_miRNA.csv"))
