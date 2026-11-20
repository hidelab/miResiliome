rm(list = ls())
print("Finding DE genes")

options(StringAsFactors = F)

## Load required libraries
#devtools::install_github('th1vairam/CovariateAnalysis@dev')

library(CovariateAnalysis) # get the package from 
library(data.table)
library(plyr)
library(dplyr)
library(tidyverse)
library(psych)

library(limma)
library(biomaRt)
library(RColorBrewer)
library(cqn)
library(knitr)
library(edgeR)


#library(synapser) # new synapse R client from install.packages("synapser", repos=c("https://sage-bionetworks.github.io/staging-ran", 'http://cran.fhcrc.org'))
#library(githubr) # get the package from 
#devtools::install_github('brian-bot/githubr')





data.dir      <- "data/preprocessed/"
output.dir    <- paste0("output","/")
input.dir     <- paste0(output.dir,"preprocessed/")
output.DEG    <- paste0(output.dir,"DEGenes_bulk")
if(!dir.exists(output.DEG)) dir.create(output.DEG,recursive = T)





primaryVariable     <- c("cogdx", "Diagnosis", "APOE4")
FactorCovariates    <- c("Batch", "Sex", "race", "spanish", "cogdx", "Diagnosis", "APOE4")
ContCovariates      <- c("RINcontinuous", "RINcontinuous2", "age_death", "pmi", "educ", 
                         "PCT_PF_READS_ALIGNED", "PCT_CODING_BASES",
                         "PCT_INTERGENIC_BASES", "PCT_INTRONIC_BASES", 
                         "PCT_RIBOSOMAL_BASES")
  

NEW.COUNTS          <- readRDS(paste0(data.dir,"/NEWCOUNTS.RDS"))
COVARIATES          <- readRDS(paste0(data.dir,"/COVARIATES.RDS"))
GENE.PARAM          <- readRDS(paste0(data.dir,"/GENEPARAM.RDS"))
postAdjustCovars    <- readRDS(paste0(data.dir,"/POSTCOVARS.RDS"))
cnts                <- NEW.COUNTS[,rownames(COVARIATES)]


# From Castanho and Naderi et al 2025

clinical.Tot           <- readRDS(paste0(data.dir,"/clinicalTotal.RDS"))
clinical.Tot$Diagnosis <- as.character(clinical.Tot$Diagnosis)



#Defining resilience

cogdx.resils    <- c(1)
cogdx.AD        <- c(4)
cogdx.normal    <- c(1)

braak.resils    <- c(3,4,5,6)
braak.AD        <- c(3,4,5,6)

plq.resils      <- c(1,2)
plq.AD          <- c(1,2)
age.cut         <- 80
age.cut.AD      <- 0
age.cut.norm    <- 0





clinical.Tot   <- mutate(clinical.Tot, RESIL = "OTHER")

resilients     <- clinical.Tot %>% dplyr::filter(., 
                                          cogdx.x %in% cogdx.resils,
                                          braaksc %in% braak.resils,
                                          ceradsc %in% plq.resils,
                                          age_death.x >= age.cut) %>% 
  mutate(.,RESIL = "Resilient")#%>% nrow

ADcontrols     <- clinical.Tot %>% filter(., 
                                          cogdx.x %in% cogdx.AD,
                                          braaksc %in% braak.AD,
                                          ceradsc %in% plq.AD,
                                          age_death.x >= age.cut.AD) %>% 
  mutate(.,RESIL = "AD")

normals     <- clinical.Tot %>% filter(., 
                                       cogdx.x %in% cogdx.normal,
                                       braaksc %in% c(0,1,2),
                                       ceradsc %in% c(3,4),
                                       age_death.x >= age.cut.norm) %>% 
  mutate(.,RESIL = "Control")

presym     <- clinical.Tot %>% filter(., 
                                      cogdx.x %in% c(2),
                                      braaksc %in% braak.resils,
                                      ceradsc %in% plq.resils) %>% 
  mutate(.,RESIL = "Presym")

clinical.resil  <- rbind(resilients,ADcontrols,normals,presym)

clinical.Tot[clinical.Tot$SampleID %in% resilients$SampleID,]$RESIL <- "Resilient"
clinical.Tot[clinical.Tot$SampleID %in% ADcontrols$SampleID,]$RESIL <- "AD"
clinical.Tot[clinical.Tot$SampleID %in% normals$SampleID,   ]$RESIL <- "Control"
clinical.Tot[clinical.Tot$SampleID %in% presym$SampleID,   ]$RESIL <- "Presym"

age.update <- TRUE
if(age.update){
  temp <- readxl::read_excel("data/ROSMAP Full Metadata/dataset_874_basic_03-17-2020.xlsx")
  temp$projid <- as.numeric(temp$projid)
  
  temp <- temp[,c("projid","age_death")]
  temp <- dplyr::left_join(clinical.Tot[,c("projid","age_death.x")],
                           temp, by = "projid")
  all(temp$projid == clinical.Tot$projid)
  clinical.Tot$age_death.x <- temp$age_death 
}




clinical.Tot$Diagnosis <- clinical.Tot$RESIL


clinical.Tot$Diagnosis <- as.factor(clinical.Tot$Diagnosis)




COVARIATES           <- clinical.Tot  
rownames(COVARIATES) <- COVARIATES$SampleID
colnames(COVARIATES) <- gsub("\\.x", "",colnames(COVARIATES))


setdiff(postAdjustCovars,colnames(COVARIATES))
# DE Gene Expression Resilience and AD.


DESIGN = CovariateAnalysis::getDesignMatrix(COVARIATES[, c("Diagnosis", postAdjustCovars), drop = F],
                         Intercept = F)
DESIGN$design = DESIGN$design[,linColumnFinder(DESIGN$design)$indepCols]


TMP.COUNTS <- NEW.COUNTS[, rownames(DESIGN$design)]


DESIGN$design <- DESIGN$design[,linColumnFinder(DESIGN$design)$indepCols]
VOOM.WEIGHTS  <- voom(TMP.COUNTS, design=DESIGN$design, plot=T)
FIT           <- lmFit(log(TMP.COUNTS), DESIGN$design,
                       weights = VOOM.WEIGHTS$weights)




      contrast      <- makeContrasts(contrasts=c("DiagnosisAD-DiagnosisResilient"),
                                     levels = colnames(FIT$coefficients))
      FIT.CONTR     <- contrasts.fit(FIT, contrasts=contrast)
      FIT.CONTR     <- eBayes(FIT.CONTR)
      tT            <- topTable(FIT.CONTR, adjust="fdr", sort.by="p", number=Inf)
          
      
      
      DE = lapply(1:dim(contrast)[2], function(i, FIT){
        topTable(FIT, coef=i, number = 50000, confint = T) %>%
          rownameToFirstColumn('ensembl_gene_id')
      }, FIT.CONTR)
      names(DE) = c('AD-RESILIENT')
      
      DE = DE %>% 
        rbindlist(idcol = 'Comparison') %>%
        left_join(GENE.PARAM %>%
                    dplyr::select(ensembl_gene_id, hgnc_symbol,
                                  percentage_gc_content, gene.length) %>%
                    unique()) %>%
        dplyr::mutate(Direction = logFC/abs(logFC),
                      Direction = factor(Direction, c(-1,1), c('-1' = 'DOWN', '1' = 'UP')),
                      Direction = as.character(Direction))
      DE$Direction[DE$adj.P.Val > 0.1 | abs(DE$logFC) < log2(1.1)] = 'NONE'

      
      
      
      log_con <- file(paste0(output.DEG,"/README.txt"), open = "a")
      sink(log_con)
      cat("\n =================================================",file = log_con,sep = "\n",append = T)
      cat( "\n Comparison: " , "AD-RESILIENT",file = log_con,sep = " ",append = T)
      cat("\n =========== DE Genes disribution ================ \n",file = log_con,sep = "\n",append = T)
      print(table(DE$Direction))
      cat("\n =================================================",file = log_con,sep = "\n",append = T)
      sink(NULL)
      
      
      
      
      
      
      
      

saveRDS(DE,paste0(output.DEG,"/ADvsResilient.RDS"))
write.csv(DE,paste0(output.DEG,"/ADvsResilient.csv"))




#
#
# Comparing AD and Normal
#
#


     
      contrast      <- makeContrasts(contrasts=c("DiagnosisAD-DiagnosisControl"),
                                     levels = colnames(FIT$coefficients))
      FIT.CONTR     <- contrasts.fit(FIT, contrasts=contrast)
      FIT.CONTR     <- eBayes(FIT.CONTR)
      tT            <- topTable(FIT.CONTR, adjust="fdr", sort.by="p", number=Inf)
      
      
      DE = lapply(1:dim(contrast)[2], function(i, FIT){
        topTable(FIT, coef=i, number = 50000, confint = T) %>%
          rownameToFirstColumn('ensembl_gene_id')
      }, FIT.CONTR)
      names(DE) = c('AD-CONTROL')
      
      DE = DE %>% 
        rbindlist(idcol = 'Comparison') %>%
        left_join(GENE.PARAM %>%
                    dplyr::select(ensembl_gene_id, hgnc_symbol,
                                  percentage_gc_content, gene.length) %>%
                    unique()) %>%
        dplyr::mutate(Direction = logFC/abs(logFC),
                      Direction = factor(Direction, c(-1,1), c('-1' = 'DOWN', '1' = 'UP')),
                      Direction = as.character(Direction))
      DE$Direction[DE$adj.P.Val > 0.1 | abs(DE$logFC) < log2(1.1)] = 'NONE'


      log_con <- file(paste0(output.DEG,"/README.txt"), open = "a")
      sink(log_con)
      cat("\n =================================================",file = log_con,sep = "\n",append = T)
      cat( "\n Comparison: " , "AD-CONTROL",file = log_con,sep = " ",append = T)
      cat("\n =========== DE Genes disribution ================ \n",file = log_con,sep = "\n",append = T)
      print(table(DE$Direction))
      cat("\n =================================================",file = log_con,sep = "\n",append = T)
      sink(NULL)
      
          
saveRDS(DE,paste0(output.DEG,"/ADvsControl.RDS"))
write.csv(DE,paste0(output.DEG,"/ADvsControl.csv"))





#
#
# Comparing Normal and Resilience.
#
#


      
      contrast      <- makeContrasts(contrasts=c("DiagnosisResilient-DiagnosisControl"),
                                     levels = colnames(FIT$coefficients))
      FIT.CONTR     <- contrasts.fit(FIT, contrasts=contrast)
      FIT.CONTR     <- eBayes(FIT.CONTR)
      tT            <- topTable(FIT.CONTR, adjust="fdr", sort.by="p", number=Inf)
      
      
      DE = lapply(1:dim(contrast)[2], function(i, FIT){
        topTable(FIT, coef=i, number = 50000, confint = T) %>%
          rownameToFirstColumn('ensembl_gene_id')
      }, FIT.CONTR)
      names(DE) = c('RESILIENT-CONTROL')
      
      DE = DE %>% 
        rbindlist(idcol = 'Comparison') %>%
        left_join(GENE.PARAM %>%
                    dplyr::select(ensembl_gene_id, hgnc_symbol,
                                  percentage_gc_content, gene.length) %>%
                    unique()) %>%
        dplyr::mutate(Direction = logFC/abs(logFC),
                      Direction = factor(Direction, c(-1,1), c('-1' = 'DOWN', '1' = 'UP')),
                      Direction = as.character(Direction))
      DE$Direction[DE$adj.P.Val > 0.1 | abs(DE$logFC) < log2(1.1)] = 'NONE'
      
      
      log_con <- file(paste0(output.DEG,"/README.txt"), open = "a")
      sink(log_con)
      cat("\n =================================================",file = log_con,sep = "\n",append = T)
      cat( "\n Comparison: " , "RESILIENT-CONTROL",file = log_con,sep = " ",append = T)
      cat("\n =========== DE Genes disribution ================ \n",file = log_con,sep = "\n",append = T)
      print(table(DE$Direction))
      cat("\n =================================================",file = log_con,sep = "\n",append = T)
      sink(NULL)
      


saveRDS(DE,paste0(output.DEG,"/ResilientvsControl.RDS"))
write.csv(DE,paste0(output.DEG,"/ResilientvsControl.csv"))






# Presym

contrast      <- makeContrasts(contrasts=c("DiagnosisResilient-DiagnosisPresym"),
                               levels = colnames(FIT$coefficients))
FIT.CONTR     <- contrasts.fit(FIT, contrasts=contrast)
FIT.CONTR     <- eBayes(FIT.CONTR)
tT            <- topTable(FIT.CONTR, adjust="fdr", sort.by="p", number=Inf)


DE = lapply(1:dim(contrast)[2], function(i, FIT){
  topTable(FIT, coef=i, number = 50000, confint = T) %>%
    rownameToFirstColumn('ensembl_gene_id')
}, FIT.CONTR)
names(DE) = c('RESILIENT-PRESYM')

DE = DE %>% 
  rbindlist(idcol = 'Comparison') %>%
  left_join(GENE.PARAM %>%
              dplyr::select(ensembl_gene_id, hgnc_symbol,
                            percentage_gc_content, gene.length) %>%
              unique()) %>%
  dplyr::mutate(Direction = logFC/abs(logFC),
                Direction = factor(Direction, c(-1,1), c('-1' = 'DOWN', '1' = 'UP')),
                Direction = as.character(Direction))
DE$Direction[DE$adj.P.Val > 0.1 | abs(DE$logFC) < log2(1.1)] = 'NONE'


log_con <- file(paste0(output.DEG,"/README.txt"), open = "a")
sink(log_con)
cat("\n =================================================",file = log_con,sep = "\n",append = T)
cat( "\n Comparison: " , "RESILIENT-PRESYM",file = log_con,sep = " ",append = T)
cat("\n =========== DE Genes disribution ================ \n",file = log_con,sep = "\n",append = T)
print(table(DE$Direction))
cat("\n =================================================",file = log_con,sep = "\n",append = T)
sink(NULL)



saveRDS(DE,paste0(output.DEG,"/ResilientvsPresym.RDS"))
write.csv(DE,paste0(output.DEG,"/ResilientvsPresym.csv"))





# Presym + Resilient vs AD

contrast      <- makeContrasts(contrasts=c("DiagnosisAD-(DiagnosisPresym+ DiagnosisResilient)/2"),
                               levels = colnames(FIT$coefficients))
FIT.CONTR     <- contrasts.fit(FIT, contrasts=contrast)
FIT.CONTR     <- eBayes(FIT.CONTR)
tT            <- topTable(FIT.CONTR, adjust="fdr", sort.by="p", number=Inf)


DE = lapply(1:dim(contrast)[2], function(i, FIT){
  topTable(FIT, coef=i, number = 50000, confint = T) %>%
    rownameToFirstColumn('ensembl_gene_id')
}, FIT.CONTR)
names(DE) = c('AD vs RESILIENT+PRESYM')

DE = DE %>% 
  rbindlist(idcol = 'Comparison') %>%
  left_join(GENE.PARAM %>%
              dplyr::select(ensembl_gene_id, hgnc_symbol,
                            percentage_gc_content, gene.length) %>%
              unique()) %>%
  dplyr::mutate(Direction = logFC/abs(logFC),
                Direction = factor(Direction, c(-1,1), c('-1' = 'DOWN', '1' = 'UP')),
                Direction = as.character(Direction))
DE$Direction[DE$adj.P.Val > 0.1 | abs(DE$logFC) < log2(1.1)] = 'NONE'


log_con <- file(paste0(output.DEG,"/README.txt"), open = "a")
sink(log_con)
cat("\n =================================================",file = log_con,sep = "\n",append = T)
cat( "\n Comparison: " , "AD vs RESILIENT+ PRESYM",file = log_con,sep = " ",append = T)
cat("\n =========== DE Genes disribution ================ \n",file = log_con,sep = "\n",append = T)
print(table(DE$Direction))
cat("\n =================================================",file = log_con,sep = "\n",append = T)
sink(NULL)



saveRDS(DE,paste0(output.DEG,"/ADvsResilientPresym.RDS"))
write.csv(DE,paste0(output.DEG,"/ADvsResilientPresym.csv"))







# Presym vs Control

contrast      <- makeContrasts(contrasts=c("DiagnosisPresym-DiagnosisControl"),
                               levels = colnames(FIT$coefficients))
FIT.CONTR     <- contrasts.fit(FIT, contrasts=contrast)
FIT.CONTR     <- eBayes(FIT.CONTR)
tT            <- topTable(FIT.CONTR, adjust="fdr", sort.by="p", number=Inf)


DE = lapply(1:dim(contrast)[2], function(i, FIT){
  topTable(FIT, coef=i, number = 50000, confint = T) %>%
    rownameToFirstColumn('ensembl_gene_id')
}, FIT.CONTR)
names(DE) = c('PRESYM vs Control')

DE = DE %>% 
  rbindlist(idcol = 'Comparison') %>%
  left_join(GENE.PARAM %>%
              dplyr::select(ensembl_gene_id, hgnc_symbol,
                            percentage_gc_content, gene.length) %>%
              unique()) %>%
  dplyr::mutate(Direction = logFC/abs(logFC),
                Direction = factor(Direction, c(-1,1), c('-1' = 'DOWN', '1' = 'UP')),
                Direction = as.character(Direction))
DE$Direction[DE$adj.P.Val > 0.1 | abs(DE$logFC) < log2(1.1)] = 'NONE'


log_con <- file(paste0(output.DEG,"/README.txt"), open = "a")
sink(log_con)
cat("\n =================================================",file = log_con,sep = "\n",append = T)
cat( "\n Comparison: " , "PRESYM vs Control",file = log_con,sep = " ",append = T)
cat("\n =========== DE Genes disribution ================ \n",file = log_con,sep = "\n",append = T)
print(table(DE$Direction))
cat("\n =================================================",file = log_con,sep = "\n",append = T)
sink(NULL)



saveRDS(DE,paste0(output.DEG,"/PresymvsControl.RDS"))
write.csv(DE,paste0(output.DEG,"/PresymvsControl.csv"))



# AD vs Presym

contrast      <- makeContrasts(contrasts=c("DiagnosisAD-DiagnosisPresym"),
                               levels = colnames(FIT$coefficients))
FIT.CONTR     <- contrasts.fit(FIT, contrasts=contrast)
FIT.CONTR     <- eBayes(FIT.CONTR)
tT            <- topTable(FIT.CONTR, adjust="fdr", sort.by="p", number=Inf)


DE = lapply(1:dim(contrast)[2], function(i, FIT){
  topTable(FIT, coef=i, number = 50000, confint = T) %>%
    rownameToFirstColumn('ensembl_gene_id')
}, FIT.CONTR)
names(DE) = c('AD vs PRESYM')

DE = DE %>% 
  rbindlist(idcol = 'Comparison') %>%
  left_join(GENE.PARAM %>%
              dplyr::select(ensembl_gene_id, hgnc_symbol,
                            percentage_gc_content, gene.length) %>%
              unique()) %>%
  dplyr::mutate(Direction = logFC/abs(logFC),
                Direction = factor(Direction, c(-1,1), c('-1' = 'DOWN', '1' = 'UP')),
                Direction = as.character(Direction))
DE$Direction[DE$adj.P.Val > 0.1 | abs(DE$logFC) < log2(1.1)] = 'NONE'


log_con <- file(paste0(output.DEG,"/README.txt"), open = "a")
sink(log_con)
cat("\n =================================================",file = log_con,sep = "\n",append = T)
cat( "\n Comparison: " , "AD vs PRESYM",file = log_con,sep = " ",append = T)
cat("\n =========== DE Genes disribution ================ \n",file = log_con,sep = "\n",append = T)
print(table(DE$Direction))
cat("\n =================================================",file = log_con,sep = "\n",append = T)
sink(NULL)



saveRDS(DE,paste0(output.DEG,"/ADvsPRESYM.RDS"))
write.csv(DE,paste0(output.DEG,"/ADvsPRESYM.csv"))



DESIGN = CovariateAnalysis::getDesignMatrix(COVARIATES[, c(postAdjustCovars), drop = F],
                                            Intercept = F)
DESIGN$design = DESIGN$design[,linColumnFinder(DESIGN$design)$indepCols]


TMP.COUNTS <- NEW.COUNTS[, rownames(DESIGN$design)]


DESIGN$design <- DESIGN$design[,linColumnFinder(DESIGN$design)$indepCols]
VOOM.WEIGHTS  <- voom(TMP.COUNTS, design=DESIGN$design, plot=T)
FIT           <- lmFit(log(TMP.COUNTS), DESIGN$design,
                       weights = VOOM.WEIGHTS$weights)


# Residuals after normalisation
RESIDUAL.GENE_EXPRESSION = residuals.MArrayLM(FIT,
                                              log(TMP.COUNTS))

saveRDS(RESIDUAL.GENE_EXPRESSION,paste0(output.DEG,"/RESIDUAL_gene_expression.RDS"))



COVARIATES$Diagnosis != "OTHER"
temp_resid <- RESIDUAL.GENE_EXPRESSION[,COVARIATES$Diagnosis != "OTHER"]

pca_genes <- prcomp((temp_resid),scale.=T, center = T)


# Plot first 2 PCs
plotdata <- data.frame(SampleID=rownames(pca_genes$rotation), 
                       PC1=pca_genes$rotation[,1], 
                       PC2=pca_genes$rotation[,2])

plotdata <- dplyr::left_join(plotdata, rownameToFirstColumn(COVARIATES[COVARIATES$Diagnosis != "OTHER",], 'SampleID')) %>%
  dplyr::mutate(label = SampleID)

p <- ggplot(plotdata, aes(x=PC1, y=PC2))
p <- p + geom_point(aes(color=Diagnosis))
p <- p + theme_bw() + theme(legend.position="right")
p


