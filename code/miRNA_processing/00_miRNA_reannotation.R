
# Library initialization 
rm(list =ls())
# args = commandArgs(trailingOnly=TRUE)
# if(is.na(args[1])){
#     comparison <- "ADvsResilient"  
# } else {
#     comparison <- args[1]
# }
# 

library(data.table)
library(limma)
library(dplyr)
library(synapser) 
library(CovariateAnalysis)
library(stringr)
library(parallel)
library(nanostringr)
library("variancePartition")
library("edgeR")
library("BiocParallel")

num.cores <- detectCores()

options(stringsAsFactors = F)
source("code/00-functions.R")


# synLogin()
# 
# mir_ID <- 'syn5856115';
# mir_raw <- synGet(mir_ID)
# Inspect the zip file above and see what's necessary


file <- read_rcc(path = "data/ROSMAP_arraymiRNA_raw/")

ROSMAP.mirs <- file$raw

# The following file corresponds to Nanostring nCounter Human V1.2 provided by
# the Nanostring staff. 

NANOST.mirs <- read.csv("data/GeneSets/mirRefNanostring.csv")
tail(NANOST.mirs)



NANOST.mirs$miRNA.Probe.ID <- str_sub(NANOST.mirs$miRNA.Probe.ID,end = -3)
dup.inds <- which(table(NANOST.mirs$miRNA.Probe.ID) >1)
dup.inds <- names(dup.inds)
NANOST.mirs[NANOST.mirs$miRNA.Probe.ID %in% dup.inds,]


# Latest miRBase (22.1) as of January 2020

library(Biostrings)

# 5/22/18,8:00PM from: ftp://mirbase.org/pub/mirbase/CURRENT/database_files/
# mirBase 22.1 
# Variable names 
  #`auto_mature`
  #`mature_name` 
  #`previous_mature_id` 
  #`mature_acc` 
  #`evidence`
  #`experiment`  
  #`similarity` 
  #`dead_flag` 
## We filter for miRNAs with dead flag as well

mat.mirs <- fread("data/GeneSets/miRBase_mirna_mature.txt")
mat.mirs <- mat.mirs[mat.mirs$V8 == 0,]
mat.mirs <- mat.mirs[,c(2,4)]
names(mat.mirs) <- c("miRNA","Accession")
mat.mirs <- (inner_join(mat.mirs,NANOST.mirs, by = c("Accession" = "Accession")))
mat.dups <- mat.mirs[,c(1,5)] 
#mat.dups <- mat.mirs[mat.mirs$miRNA.Probe.ID %in% dup.inds,c(1,5)] 

mat.dups <- as.data.table(mat.dups)
mat.dups <- mat.dups[, list(miRNA = paste(miRNA, collapse="+")), by = miRNA.Probe.ID]


mat.dups <- inner_join(mat.dups,ROSMAP.mirs[,2:3],by = c("miRNA.Probe.ID" = "Accession"))
dim(mat.dups)
names(mat.dups)

mat.dups$Name   <- gsub("\\s\\(.*\\)","", mat.dups$Name)
mat.dups$Name   <- gsub("\\+","_", mat.dups$Name)
mat.dups$miRNA  <- gsub("\\+","_", mat.dups$miRNA)

names(mat.dups) <- c("miRNA.Probe.ID", "miRNA_new",  "miRNA_old")

saveRDS(mat.dups,"data/preprocessed/miRNA_old_new_conversion.RDS")


# It's all the test from here. 