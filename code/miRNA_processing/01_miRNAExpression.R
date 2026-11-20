args = commandArgs(trailingOnly=TRUE)
if(is.na(args[1])){
    setup <- ""
    warning("The experimental design specifications must be set.")
} else {
    setup <- args[1]
}

options(StringAsFactors = F)



library(data.table)
library(limma)
library(dplyr)
library(synapser) 
library(CovariateAnalysis)
library(stringr)
library(plyr)
library(psych)
library(tidyverse)



age.update <- T

data.dir      <- "data/preprocessed/"
output.dir    <- paste0("Output", setup,"/")
output.DEG    <- paste0(output.dir,"DEGenes")
if(!dir.exists(output.DEG)) dir.create(output.DEG,recursive = T)

# 
# input.dir      <- paste0(output.dir,"preprocessed/")
# # The first output dir is the current one for pathway summaries. Keep if things fall apart
# 
# output.dir.dep     <- paste0(output.dir,"DEPandSummaries/")





#covs       <- readRDS(paste0(input.dir,"clinicalTotal.RDS"))
mirExp     <- fread("data/ROSMAP_mirna/ROSMAP_arraymiRNA.gct",skip = 2)

mirExp     <- data.frame(mirExp)
rownames(mirExp) <- mirExp[,1]
mirExp     <- mirExp[,-c(1,2)]

# synLogin()
# MIR_META_ID  <- 'syn5857921'
# MIR_META_OBJ <- synGet(MIR_META_ID)
# MIR_META     <- read.table(MIR_META_OBJ$path, header=T, sep=',', check.names = F)

MIR_META <- read.table("data/ROSMAP_mirna/ROSMAP_arraymiRNA_covariates.csv",
                       header = T, sep = ",",check.names = F)

clinical2 <- fread("data/Covariates and Metadata/ROSMAP_Clinical_2019-05_v2.csv")
clinical3 <- fread("data/Covariates and Metadata/ROSMAP_IDkey.csv")
if(age.update){
  temp <- readxl::read_excel("data/ROSMAP Full Metadata/dataset_874_basic_03-17-2020.xlsx")
  temp$projid <- as.numeric(temp$projid)
  
  temp <- temp[,c("projid","age_death")]
  temp <- dplyr::left_join(clinical2[,c("projid","age_death")],
                           temp, by = "projid")
  all(temp$projid == clinical2$projid)
  clinical2$age_death <- temp$age_death.y 
}





# Merging clinical variables for future purposes

inds           <- clinical3 %>% filter(., mirna_id !="", 
                                       !is.na(mirna_id))
inds           <- inds[inds$mirna_id %in% colnames(mirExp),]
inds           <- (inds[!duplicated(inds$projid),])
MIR_META       <- (MIR_META[!duplicated(MIR_META$Sample_ID),])

clinical2      <- inner_join(inds, clinical2, by = "projid")
clinical.tot   <- inner_join(clinical2, MIR_META,
                             by ="projid")

#clinical.tot[clinical.tot$age_death == "90+",]$age_death <- 90

clinical.tot   <- clinical.tot[clinical.tot$RIN >= 5.0,]
clinical.tot   <- clinical.tot[!duplicated(clinical.tot$mirna_id),]
clinical.tot   <- clinical.tot[!is.na(clinical.tot$RIN) & !is.na(clinical.tot$pmi),]

saveRDS(clinical.tot,file = "data/preprocessed/mirClinical.RDS")


keep           <- colnames(mirExp) %in% clinical.tot$mirna_id
mirExp         <- mirExp[,keep]



## We have a consistent miRNA expresison and covariates now. 
#  Let's correct miRNA's for fixed effects.

postAdjustCovars <- c('msex');
primaryVariable  <- c("cogdx","apoe_genotype")
FactorCovariates <- c("msex", "race", "spanish","cogdx","apoe_genotype","plate")
ContCovariates   <- c("RIN", "age_death", "pmi","educ")


clinical.tot <- as.data.frame(clinical.tot)

COVARIATES <- clinical.tot[,c(FactorCovariates,ContCovariates),drop = F]
COVARIATES <- data.frame(lapply(COVARIATES,
                                function(x){
                                  x <- sapply(x,function(y){
                                    str_replace_all(as.character(y),'\\+','')
                                    })
                                  }))
rownames(COVARIATES) <- clinical.tot$mirna_id

# Convert factor covariates to factors
COVARIATES[,FactorCovariates] = lapply(COVARIATES[,FactorCovariates], factor)
COVARIATES[,ContCovariates] = lapply(COVARIATES[,ContCovariates], as.character)
COVARIATES[,ContCovariates] = lapply(COVARIATES[,ContCovariates], as.numeric)


colnames(mirExp) == rownames(COVARIATES)
mirExp <- mirExp[,rownames(COVARIATES)]


limma::plotDensities(mirExp)

pc <- prcomp(t(mirExp), rank.=14,center = T, scale. = T)
summary(pc)
pc_mat <- pc$x



# Make sure the package below is installed. 
library(variancePartition)

form <- ~  msex + RIN + age_death + pmi + educ + apoe_genotype + cogdx + plate
#form <- ~  (1|msex) + RIN + age_death + pmi + educ + (1|apoe_genotype) + (1|cogdx) + (1|plate)



varPart <- fitExtractVarPartModel(t(pc_mat), form, COVARIATES)
res <- fitVarPartModel(t(pc_mat), form, COVARIATES)


# evaluate the collinearity score on the first model fit
# this reports the correlation matrix between coefficient estimates
# for fixed effects
# the collinearity score is the maximum absolute correlation value
# If the collinearity score > .99 then the variance partition
# estimates may be problematic
# In that case, a least one variable should be omitted
colinearityScore(res[[1]])

summary(varPart)

vp <- sortCols(varPart)
plotPercentBars(vp[1:10, ])

ggsave("output/qc/pca_covariates.pdf",width = 4,height = 5)

plotVarPart(vp[1:10,])
ggsave("output/qc/variancepart_covariates.pdf",width = 4,height = 3)


saveRDS(COVARIATES,file = "data/preprocessed/mirCovars.RDS")


corrections <- c("msex", "RIN", "age_death", "pmi","apoe_genotype")


colnames(mirExp) == rownames(COVARIATES)
    
  DESIGN         <- getDesignMatrix(COVARIATES[,(corrections),
                                              drop = F], Intercept = F)
  DESIGN$design  <- DESIGN$design[,linColumnFinder(DESIGN$design)$indepCols]
  FIT            <- lmFit(mirExp, DESIGN$design)
  RESIDUAL.miRNA <-  residuals.MArrayLM(FIT,mirExp)
  new.names      <- rownames(RESIDUAL.miRNA)
  name.converse  <- readRDS("data/preprocessed/miRNA_old_new_conversion.RDS")
  new.names      <- data.frame("miRNA_old" = new.names)
  new.names      <- inner_join(new.names,name.converse, by ="miRNA_old")
  RESIDUAL.miRNA <- RESIDUAL.miRNA[new.names$miRNA_old,]
  
  rownames(RESIDUAL.miRNA) <- new.names$miRNA_new
  saveRDS(rownames(RESIDUAL.miRNA),"data/preprocessed/brainmiRNA.RDS")
  saveRDS(RESIDUAL.miRNA,"data/preprocessed/ResidualmiRNA.RDS")
  
  
  
  mirExp <- mirExp[new.names$miRNA_old,]
  
  rownames(mirExp) <- new.names$miRNA_new
  saveRDS(mirExp,"data/preprocessed/miRNAExpression.RDS")