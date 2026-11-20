rm(list = ls())
library(PanomiR)
library(dplyr)
library(qvalue)




# Covariate definition



data.dir      <- "data/preprocessed/"
output.dir    <- paste0("output","/")
input.dir     <- paste0(output.dir,"preprocessed/")
output.DEG    <- paste0(output.dir,"DEPaths_bulk_update/")
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

clinical.Tot$age_death.x <- clinical.Tot$age_death.y




clinical.Tot$Diagnosis <- clinical.Tot$RESIL


clinical.Tot$Diagnosis <- as.factor(clinical.Tot$Diagnosis)

table(clinical.Tot$Diagnosis)
#saveRDS(clinical.Tot,"data/updated2023_bulk_ROSMAP_clinicalTotal.RDS")


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


table(clinical.Tot$Diagnosis)


COVARIATES           <- clinical.Tot  
rownames(COVARIATES) <- COVARIATES$SampleID
colnames(COVARIATES) <- gsub("\\.x", "",colnames(COVARIATES))
#colnames(COVARIATES) <- gsub("\\.y", "",colnames(COVARIATES))







primaryVariable     <- c("cogdx", "Diagnosis", "APOE4")
FactorCovariates    <- c("Batch", "Sex", "race", "spanish", "cogdx", "Diagnosis", "APOE4")
ContCovariates      <- c("RINcontinuous", "RINcontinuous2", "age_death", "pmi", "educ", 
                         "PCT_PF_READS_ALIGNED", "PCT_CODING_BASES",
                         "PCT_INTERGENIC_BASES", "PCT_INTRONIC_BASES", 
                         "PCT_RIBOSOMAL_BASES")



saveRDS(COVARIATES,"data/preprocessed/UPDATED_COVARIATES_ROSMAP_withResilience2025.RDS")


# Directory to the pathway gene membership
pathways2 <- readRDS('data/preprocessed/MSigDBPathGeneTab2024.RDS')

# Directory to the normalized gene expression
genes.counts2 <- readRDS("data/preprocessed/NEWCOUNTS.RDS")
#genes.counts2 <- readRDS("output/DEGenes_bulk/RESIDUAL_gene_expression.RDS")

genes.counts2 <- as.data.frame(genes.counts2)


AD.covariates <- COVARIATES

# Condition for finding DE pathways 
# TO DO: should it be a column in AD.covariates?
condition = 'Diagnosis'

# Where to write the output data, DE pathways

# These are the covariates that you may need to adjust for. 
# Can include both technical such as RIN or clinical such as Age at Death
# A character vector?
adj.covars <- readRDS("data/preprocessed/POSTCOVARS.RDS")

# Not necessary if you are not using ROSMAP data remove this and all other instances
# TO DO: could it be used as analysis name?
tisss <- "ROSMAP_DLPFC_ADvsControl"
# Making sure factor covariates are treated as such
# Any correction should be on either numeric or factor. No character type allowed

# Making show gene expression and covariates data confirm
genes.counts2 <- genes.counts2[,(AD.covariates$SampleID)]


de.paths <- PanomiR::differentialPathwayAnalysis(geneCounts = genes.counts2, 
                                     pathways = pathways2,
                                     covariates = AD.covariates,
                                     condition = condition,
                                     outDir  = output.DEG,
                                     saveOutName =  paste(tisss,"_DEP_AD_vs_Control.RDS"),
                                     adjustCovars =adj.covars,
                                     contrastConds = "DiagnosisAD-DiagnosisControl",
                                     trim = 0)

temp <- de.paths$DEP
temp <- temp[order(temp$P.Value),]
temp$qvalue <- qvalue(temp$P.Value)$lfdr
temp %>% 
    write.csv(file = paste0(output.DEG, tisss, "_diffPathways", ".csv"),row.names = T)

temp <- tibble::rownames_to_column(temp)
temp$qvalue <- qvalue(temp$P.Value)$lfdr


de.paths$PathwayResiduals |> saveRDS(
                                     paste0(output.DEG,
                                            tisss, "_Residuals",
                                            ".RDS")
                                     )



tisss <- "ROSMAP_DLPFC_ADvsResilient"
de.paths <- PanomiR::differentialPathwayAnalysis(geneCounts = genes.counts2,
                                                 pathways = pathways2,
                                                 covariates = AD.covariates,
                                                 condition = condition,
                                                 outDir  = output.DEG,
                                                 saveOutName =  paste(tisss,"_DEP_AD_vs_Control.RDS"),
                                                 adjustCovars =adj.covars,
                                                 contrastConds = "DiagnosisAD-DiagnosisResilient",
                                                 trim = 0)

temp <- de.paths$DEP
temp <- temp[order(temp$P.Value),]
temp$qvalue <- qvalue(temp$P.Value)$lfdr


temp %>% 
    write.csv(file = paste0(output.DEG, tisss, "_diffPathways", ".csv"),row.names = T)

temp <- tibble::rownames_to_column(temp)


de.paths$PathwayResiduals |> saveRDS(
    paste0(output.DEG,
           tisss, "_Residuals",
           ".RDS")
)



tisss <- "ROSMAP_DLPFC_ResilientvsControl"
de.paths <- PanomiR::differentialPathwayAnalysis(geneCounts = genes.counts2,
                                                 pathways = pathways2,
                                                 covariates = AD.covariates,
                                                 condition = condition,
                                                 outDir  = output.DEG,
                                                 saveOutName =  paste(tisss,"_DEP_AD_vs_Control.RDS"),
                                                 adjustCovars =adj.covars,
                                                 contrastConds = "DiagnosisResilient-DiagnosisControl",
                                                 trim = 0)

temp <- de.paths$DEP
temp <- temp[order(temp$P.Value),]
temp$qvalue <- qvalue(temp$P.Value)$lfdr

temp %>% 
    write.csv(file = paste0(output.DEG, tisss, "_diffPathways", ".csv"),row.names = T)

temp <- tibble::rownames_to_column(temp)


de.paths$PathwayResiduals |> saveRDS(
    paste0(output.DEG,
           tisss, "_Residuals",
           ".RDS")
)



tisss <- "ROSMAP_DLPFC_ResilientvsPresym"
de.paths <- PanomiR::differentialPathwayAnalysis(geneCounts = genes.counts2,
                                                 pathways = pathways2,
                                                 covariates = AD.covariates,
                                                 condition = condition,
                                                 outDir  = output.DEG,
                                                 saveOutName =  paste(tisss,"_DEP_AD_vs_Control.RDS"),
                                                 adjustCovars =adj.covars,
                                                 contrastConds = "DiagnosisResilient-DiagnosisPresym",
                                                 trim = 0)

temp <- de.paths$DEP
temp <- temp[order(temp$P.Value),]
temp$qvalue <- qvalue(temp$P.Value)$lfdr
temp %>% 
    write.csv(file = paste0(output.DEG, tisss, "_diffPathways",  ".csv"),row.names = T)

temp <- tibble::rownames_to_column(temp)


de.paths$PathwayResiduals |> saveRDS(
    paste0(output.DEG,
           tisss, "_Residuals",
           ".RDS")
)









tisss <- "ROSMAP_DLPFC_ADvsResilientPresym"
de.paths <- PanomiR::differentialPathwayAnalysis(geneCounts = genes.counts2,
                                                 pathways = pathways2,
                                                 covariates = AD.covariates,
                                                 condition = condition,
                                                 outDir  = output.DEG,
                                                 saveOutName =  paste(tisss,"_DEP_AD_vs_Control.RDS"),
                                                 adjustCovars =adj.covars,
                                                 contrastConds = "DiagnosisAD-(DiagnosisPresym+DiagnosisResilient)/2",
                                                 trim = 0)

temp <- de.paths$DEP
temp <- temp[order(temp$P.Value),]
temp$qvalue <- qvalue(temp$P.Value)$lfdr
temp %>% 
    write.csv(file = paste0(output.DEG, tisss, "_diffPathways", ".csv"),row.names = T)

temp <- tibble::rownames_to_column(temp)


de.paths$PathwayResiduals |> saveRDS(
    paste0(output.DEG,
           tisss, "_Residuals",
           ".RDS")
)





tisss <- "ROSMAP_DLPFC_ADvsPresym"
de.paths <- PanomiR::differentialPathwayAnalysis(geneCounts = genes.counts2,
                                                 pathways = pathways2,
                                                 covariates = AD.covariates,
                                                 condition = condition,
                                                 outDir  = output.DEG,
                                                 saveOutName =  paste(tisss,"_DEP_AD_vs_Control.RDS"),
                                                 adjustCovars =adj.covars,
                                                 contrastConds = "DiagnosisAD-DiagnosisPresym",
                                                 trim = 0)

temp <- de.paths$DEP
temp <- temp[order(temp$P.Value),]
temp$qvalue <- qvalue(temp$P.Value)$lfdr

temp %>% 
    write.csv(file = paste0(output.DEG, tisss, "_diffPathways", ".csv"),row.names = T)

temp <- tibble::rownames_to_column(temp)


de.paths$PathwayResiduals |> saveRDS(
    paste0(output.DEG,
           tisss, "_Residuals",
           ".RDS")
)





tisss <- "ROSMAP_DLPFC_ADvsNonAD"
de.paths <- PanomiR::differentialPathwayAnalysis(geneCounts = genes.counts2,
                                                 pathways = pathways2,
                                                 covariates = AD.covariates,
                                                 condition = condition,
                                                 outDir  = output.DEG,
                                                 saveOutName =  paste(tisss,"_DEP_AD_vs_Control.RDS"),
                                                 adjustCovars =adj.covars,
                                                 contrastConds = "DiagnosisAD-(DiagnosisPresym+DiagnosisResilient+ DiagnosisControl)/3",
                                                 trim = 0)

temp <- de.paths$DEP
temp <- temp[order(temp$P.Value),]
temp$qvalue <- qvalue(temp$P.Value)$lfdr

temp %>% 
  write.csv(file = paste0(output.DEG, tisss, "_diffPathways", ".csv"),row.names = T)

temp <- tibble::rownames_to_column(temp)


de.paths$PathwayResiduals |> saveRDS(
  paste0(output.DEG,
         tisss, "_Residuals",
         ".RDS")
)





