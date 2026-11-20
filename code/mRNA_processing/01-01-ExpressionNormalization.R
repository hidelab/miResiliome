rm(list = ls())


POURYA_VOOM    <- T
PLOT.META      <- F
OUTLIER.DETECT <- F
SVA.ADJUST     <- F
DE.ANALYSIS    <- F
age.update     <- T

# The above controls are meant to be developed into the pipeline later. 
# Please do not change. If need the functionality just run the code blocks 
# within respective logicals. 
## It is assumed your working directory is where this file
## Load required libraries
#devtools::install_github('th1vairam/CovariateAnalysis@dev')

library(CovariateAnalysis) # get the package from 
library(data.table)
library(plyr)
library(tidyverse)
library(psych)
library(limma)
library(biomaRt)
library(RColorBrewer)
library(cqn)
library(knitr)
library(edgeR)
library(synapser) # new synapse R client from install.packages("synapser", repos=c("https://sage-bionetworks.github.io/staging-ran", 'http://cran.fhcrc.org'))
library(githubr) # get the package from 
#devtools::install_github('brian-bot/githubr')

INIT = T
if(INIT){
  synLogin()
  
  library(doParallel)
  library(foreach)
  
  cl = makeCluster(detectCores()-2)
  registerDoParallel(cl)
  
  options(xtable.type="html")
  
  
  # Download expression data
  COUNT_ID <- 'syn8691134';
  ALL_USED_IDs <- COUNT_ID
  COUNT_OBJ <- synGet(COUNT_ID)
  COUNT <- read.table(COUNT_OBJ$path, header=T, sep='\t', check.names = F, row.names = 1)
  COUNT[,grep('150_120419', colnames(COUNT))[2]] = NULL
  
  # Convert rownames of counts from tracking id to ensemble gene id
  tmp = data.frame(Gene.ID = rownames(COUNT)) %>%
    dplyr::mutate(ID = Gene.ID) %>%
    tidyr::separate(ID, c('ensembl_gene_id', 'position'), sep = '\\.')
  rownames(tmp) = tmp$Gene.ID
  rownames(COUNT) = tmp[rownames(COUNT), 'ensembl_gene_id']
  
  # Get clinical metadata
  METADATA.CLINICAL_ID <- 'syn3191087'
  ALL_USED_IDs[length(ALL_USED_IDs)+1] = METADATA.CLINICAL_ID
  METADATA.CLINICAL_OBJ <- synGet(METADATA.CLINICAL_ID)
  METADATA.CLINICAL <- read.table(METADATA.CLINICAL_OBJ$path,sep=',',header=T)
  
  # Updating age variable.
 
   if(age.update){
    temp <- readxl::read_excel("data/ROSMAP Full Metadata/dataset_874_basic_03-17-2020.xlsx")
    temp$projid <- as.numeric(temp$projid)
    
    temp <- temp[,c("projid","age_death")]
    temp <- dplyr::left_join(METADATA.CLINICAL[,c("projid","age_death")],
                             temp, by = "projid")
    all(temp$projid == METADATA.CLINICAL$projid)
    METADATA.CLINICAL$age_death <- temp$age_death.y 
  }

  
    # Get clinical metadata with uncensored ages
  # METADATA.CLINICAL_ID1 <- 'syn7116000'
  # ALL_USED_IDs[length(ALL_USED_IDs)+1] = METADATA.CLINICAL_ID1
  # METADATA.CLINICAL_OBJ1 <- synGet(METADATA.CLINICAL_ID1)
  # METADATA.CLINICAL1 <- read.table(METADATA.CLINICAL_OBJ1$path,sep=',',header=T)
  
  
  # METADATA.CLINICAL1 <- read.table(file = "data/Covariates and Metadata/ROSMAP_Clinical_2019-05_v2.csv",
  #                                  sep = ",",header = T)
  
  # Get technical covariates
  METADATA.TECH_ID <- 'syn4300313'
  ALL_USED_IDs[length(ALL_USED_IDs)+1] = METADATA.TECH_ID
  METADATA.TECH_OBJ <- synGet(METADATA.TECH_ID)
  METADATA.TECH <- read.table(METADATA.TECH_OBJ$path,sep='\t',header=T)
  
  # Get picard metrics from synapse
  METADATA.PICARD_ID <- 'syn8698240';
  ALL_USED_IDs[length(ALL_USED_IDs)+1] = METADATA.PICARD_ID
  METADATA.PICARD <- synGet(METADATA.PICARD_ID)$path %>%
    data.table::fread() %>%
    dplyr::rename(Sampleid = sample)
  colnames(METADATA.PICARD) = gsub('AlignmentSummaryMetrics__','',colnames(METADATA.PICARD))
  colnames(METADATA.PICARD) = gsub('RnaSeqMetrics__','',colnames(METADATA.PICARD))
  
  # Pourya added the line below, remove if things collapse
  colnames(METADATA.PICARD)[42]   <- "PF_ALIGNED_BASES1"
  METADATA.PICARD[,37]
  
  # Fix error in technical covariates data
  KEY_ID <- 'syn3382527'
  ALL_USED_IDs[length(ALL_USED_IDs)+1] = KEY_ID  
  KEY <- synGet(KEY_ID)$path %>%
    read_csv %>% 
    dplyr::filter(!is.na(rnaseq_id)) %>%
    dplyr::select(projid, rnaseq_id) %>%
    tidyr::separate(rnaseq_id, c('a','b','batch'), sep = '_') %>% 
    unite(Sampleid, a, b) %>%
    dplyr::select(-batch) %>%
    unique
  
  # Match technical and clinical covariates
  METADATA <- METADATA.TECH %>%
    dplyr::left_join(METADATA.PICARD) %>%
    dplyr::select(-projid) %>%
    dplyr::left_join(KEY) %>%
    dplyr::left_join(METADATA.CLINICAL) # %>%
  #dplyr::select(-age_first_ad_dx, -age_death, -age_at_visit_max) #%>%
  #dplyr::left_join(METADATA.CLINICAL1)
  
  # Pick higher quality RIN batch for sample 492_120515
  METADATA <- METADATA %>%
    dplyr::group_by(Sampleid) %>%
    dplyr::top_n(1, RINcontinuous)
  
  # Get gene specific parameters from synapse
  GENE.PARAM = synGet('syn8449369')$path %>%
    data.table::fread(data.table = FALSE)
  ALL_USED_IDs = c(ALL_USED_IDs, 'syn8449369')
  
  GENE.LEN = dplyr::select(GENE.PARAM, ensembl_gene_id, gene.length) %>% 
    unique() 
  rownames(GENE.LEN) = GENE.LEN$ensembl_gene_id
  
  GENE.GC = dplyr::select(GENE.PARAM, ensembl_gene_id, percentage_gc_content) %>% 
    unique() 
  rownames(GENE.GC) = GENE.GC$ensembl_gene_id 
  
  
  
  # Remove samples with no cogdx, RIN, PMI scores and age_death
  METADATA <- METADATA %>%
    ungroup %>%
    dplyr::filter(Sampleid %in% colnames(COUNT)) %>%
    dplyr::filter(!is.na(cogdx), !is.na(braaksc), !is.na(ceradsc)) %>%
    dplyr::filter(!is.na(RINcontinuous)) %>%
    dplyr::filter(!is.na(pmi)) %>%
    dplyr::filter(!is.na(PCT_INTRONIC_BASES)) %>%
    dplyr::filter(!is.na(age_death)) %>%
    as.data.frame()
  
  # Add harmonised case-control status
  METADATA$Diagnosis = 'OTHER'
  METADATA$Diagnosis[(METADATA$cogdx == 1 | METADATA$cogdx == 2) &
                      METADATA$braaksc <= 2 & METADATA$ceradsc >= 3] = 'CONTROL'
  METADATA$Diagnosis[METADATA$cogdx == 4 & METADATA$braaksc >= 3 & METADATA$ceradsc <= 2] = 'AD'
  
  # Add sex variable 
  METADATA$Sex = 'FEMALE'
  METADATA$Sex[METADATA$msex == 1] = 'MALE'
  
  # Add apoe4 genotype (0, 1, 2)
  METADATA$APOE4 = 0
  METADATA$APOE4[METADATA$apoe_genotype %in% c(24, 34)] = 1
  METADATA$APOE4[METADATA$apoe_genotype %in% c(44)] = 2
  # METADATA$APOE4[is.na(METADATA$apoe_genotype)] = NA
  
  # Get square of RIN
  METADATA$RINcontinuous2 = METADATA$RINcontinuous^2
  
  
  
  # Match covariates to expression data
  
  indToRetain = intersect(METADATA$Sampleid, colnames(COUNT))
  removedIDs = setdiff(colnames(COUNT), METADATA$Sampleid)
  
  COUNT = COUNT[,indToRetain]
  
  rownames(METADATA) = METADATA$Sampleid
  METADATA = METADATA[indToRetain,]
  
  
  
  
  ## Warning! 
  #
  #

  
  # Dorsolateral prefrontal cortex of `r dim(COUNT)[2]` subjects from the ROS and MAP cohorts are used for the analysis. 
  #Following sample are removed due to missing metadata `r paste(removedIDs, collapse = ',')`
  
  ### Covariate clustering
  ## Determine relationship between covariates
  
  primaryVariable <- c("cogdx", "Diagnosis", "APOE4")
  FactorCovariates <- c("Batch", "Sex", "race", "spanish", "cogdx", "Diagnosis", "APOE4")
  ContCovariates <- c("RINcontinuous", "RINcontinuous2", "age_death", "pmi", "educ", 
                      "PCT_PF_READS_ALIGNED", "PCT_CODING_BASES",
                      "PCT_INTERGENIC_BASES", "PCT_INTRONIC_BASES", 
                      "PCT_RIBOSOMAL_BASES")
  
  # Find inter relation between factor covariates
  COVARIATES = METADATA[,c(FactorCovariates,ContCovariates),drop=F]
  COVARIATES <- data.frame(lapply(COVARIATES,function(x){x <- sapply(x,function(y){str_replace_all(as.character(y),'\\+','')})}))
  rownames(COVARIATES) <- METADATA$Sampleid
  
  # Convert factor covariates to factors
  COVARIATES[,FactorCovariates] = lapply(COVARIATES[,FactorCovariates], factor)
  COVARIATES[,ContCovariates] = lapply(COVARIATES[,ContCovariates], as.character)
  COVARIATES[,ContCovariates] = lapply(COVARIATES[,ContCovariates], as.numeric)
  
}






### Filter genes
# * Remove genes that have less than 1 cpm counts in at least 50% of samples per Diagnosis
# * Remove genes with missing gene length and percentage GC content

genesToAnalyze = COVARIATES %>%
  rownameToFirstColumn('Sampleid') %>%
  dlply(.(Diagnosis), .fun = function(mtd, count){
    processed.counts = getGeneFilteredGeneExprMatrix(count[,mtd$Sampleid],
                                                     MIN_GENE_CPM=1, 
                                                     MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM=0.5,
                                                     verbose = T)
    processed.counts$filteredExprMatrix$genes
  }, COUNT)

genesToAnalyze = unlist(genesToAnalyze) %>% 
  unique() %>%
  intersect(GENE.GC$ensembl_gene_id[!is.na(GENE.GC$percentage_gc_content)]) %>%
  intersect(GENE.LEN$ensembl_gene_id[!is.na(GENE.LEN$gene.length)]) %>%
  setdiff(c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"))

PROCESSED_COUNTS = getGeneFilteredGeneExprMatrix(COUNT[genesToAnalyze, ], 
                                                 MIN_GENE_CPM=0, 
                                                 MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM=0.0,
                                                 verbose = T)




#  Check gene biotype
## Define biomart object
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                host = "dec2016.archive.ensembl.org", 
                dataset = "hsapiens_gene_ensembl")

## Query biomart
Ensemble2HGNC <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                       filters = "ensembl_gene_id", 
                       values = PROCESSED_COUNTS$filteredExprMatrix$genes,
                       mart = mart)

summary(factor(Ensemble2HGNC$gene_biotype)) %>%
  rownameToFirstColumn('Biotype') %>%
  dplyr::rename(fraction = DF) %>%
  dplyr::mutate(fraction = fraction/dim(PROCESSED_COUNTS$filteredExprMatrix$genes)[1]) %>%
  dplyr::filter(fraction >= 0.01) %>% kable()




# Compute offset for gene length and gc content
CQN.GENE_EXPRESSION = cqn(PROCESSED_COUNTS$filteredExprMatrix$counts, 
                          x = GENE.GC[PROCESSED_COUNTS$filteredExprMatrix$genes$genes, 'percentage_gc_content'],
                          lengths = GENE.LEN[PROCESSED_COUNTS$filteredExprMatrix$genes$genes, 'gene.length'],
                          lengthMethod = "smooth", 
                          verbose = FALSE)

CQN.GENE_EXPRESSION$E = CQN.GENE_EXPRESSION$y + CQN.GENE_EXPRESSION$offset


### Outlier Analysis
#### Sample outliers
#    Outlier analysis is performed before library normalisation with raw cpm counts


indToRemove = c('380_120503', '500_120515')

# Find principal components of expression to plot


if(OUTLIER.DETECT){PC <- prcomp(CQN.GENE_EXPRESSION$E, scale.=T, center = T)

# Plot first 2 PCs
plotdata <- data.frame(SampleID=rownames(PC$rotation), 
                       PC1=PC$rotation[,1], 
                       PC2=PC$rotation[,2])

plotdata <- dplyr::left_join(plotdata, rownameToFirstColumn(COVARIATES, 'SampleID')) %>%
  dplyr::mutate(label = SampleID)
plotdata$label[!(plotdata$SampleID %in% indToRemove)] = ''

p <- ggplot(plotdata, aes(x=PC1, y=PC2))
p <- p + geom_point(aes(color=Batch, shape=Diagnosis, size=RINcontinuous))
p <- p + theme_bw() + theme(legend.position="right")
p <- p + geom_text(aes(label= label), size=4, hjust=0)
p


# Plot abberent distribution of logcpm counts
tmp1 = CQN.GENE_EXPRESSION$E %>%
  rownameToFirstColumn('Gene.ID') %>%
  tidyr::gather(SampleID, logCPM, -Gene.ID) %>%
  dplyr::left_join(COVARIATES %>%
                     rownameToFirstColumn('SampleID'))

p = ggplot(tmp1 %>%
             dplyr::filter(SampleID %in% indToRemove),
           aes(x = logCPM, color = SampleID)) + geom_density() 
p = p + theme(legend.position = 'top')
p
}



indToRetain = setdiff(colnames(PROCESSED_COUNTS$filteredExprMatrix$counts), indToRemove)
PROCESSED_COUNTS$filteredExprMatrix$counts = PROCESSED_COUNTS$filteredExprMatrix$counts[,indToRetain]
CQN.GENE_EXPRESSION$E = CQN.GENE_EXPRESSION$E[,indToRetain]
COVARIATES = COVARIATES[indToRetain,]


tmp = COVARIATES %>%
  dplyr::group_by(Diagnosis, cogdx) %>%
  dplyr::summarise(count = n()) %>%
  tidyr::spread(Diagnosis, count) %>% kable()









CQN.GENE_EXPRESSION$E.no.na = CQN.GENE_EXPRESSION$E
CQN.GENE_EXPRESSION$E.no.na[is.na(CQN.GENE_EXPRESSION$E.no.na)] = 0
LIB.SIZE = colSums(PROCESSED_COUNTS$filteredExprMatrix$counts)
NEW.COUNTS = (2^CQN.GENE_EXPRESSION$E) 



# Find correlation between PC's of gene expression with covariates

preAdjustedSigCovars = runPCAandPlotCorrelations(CQN.GENE_EXPRESSION$E.no.na, 
                                                 COVARIATES,
                                                 'NULL design(voom-normalized)', 
                                                 isKeyPlot=TRUE, 
                                                 MIN_PVE_PCT_PC = 1)




### Normalisation (iterative design)
#   Primary variable of interest
postAdjustCovars = c('Batch', 'Sex');

# Assign residual covariates
residualCovars    = setdiff(preAdjustedSigCovars$significantCovars,
                            c(postAdjustCovars, primaryVariable))
residualSigCovars = preAdjustedSigCovars
covariatesEffects = preAdjustedSigCovars$Effects.significantCovars[residualCovars]
postAdjustCovars  = c(postAdjustCovars, names(which.max(covariatesEffects))) %>%
  unique()

loopCount = 0 
while(length(residualSigCovars$significantCovars)!=0 && loopCount <= 20){
  writeLines(paste('Using following covariates in the model:',
                   paste(postAdjustCovars, collapse=', '),
                   'as fixed effects'))
  
  # Post adjusted design matrix
  # I have added the intercept term. Pourya
  DM1 = getDesignMatrix(COVARIATES[,postAdjustCovars,drop=F],Intercept = T)
  DM1$design = DM1$design[,linColumnFinder(DM1$design)$indepCols]
  
  # Estimate voom weights for dispersion control
  cnts = NEW.COUNTS
  cnts[is.na(cnts)] = 0
  VOOM.GENE_EXPRESSION = voom(cnts, 
                              design=DM1$design, 
                              plot=F,
                              na.rm = T)
  
  # Fit linear model using new weights and new design
  VOOM.ADJUSTED.FIT = lmFit(CQN.GENE_EXPRESSION$E,
                            design = DM1$design,
                            weights = VOOM.GENE_EXPRESSION$weights)
  
  # Residuals after normalisation
  RESIDUAL.GENE_EXPRESSION = residuals.MArrayLM(VOOM.ADJUSTED.FIT,
                                                CQN.GENE_EXPRESSION$E)
  
  # Residual covariates to choose from
  residCovars <- setdiff(c(FactorCovariates,ContCovariates), c(postAdjustCovars, primaryVariable))
  
  # Find PC of residual gene expression and significant covariates that are highly correlated with PCs
  expr = RESIDUAL.GENE_EXPRESSION
  expr[is.na(expr)] = 0
  residualSigCovars = runPCAandPlotCorrelations(expr, 
                                                COVARIATES[, residCovars, drop=F], 
                                                'adjusted design(voom-normalized)',
                                                isKeyPlot=TRUE)
  
  # Add postadjusted covariates (if any)
  residCovars = setdiff(residualSigCovars$significantCovars, c(postAdjustCovars, primaryVariable))
  covariatesEffects = residualSigCovars$Effects.significantCovars[residCovars]
  
  postAdjustCovars = c(postAdjustCovars, names(which.max(covariatesEffects)))
  loopCount = loopCount + 1
}

modelStr <- paste(paste(gsub('_','\\\\_',postAdjustCovars), collapse=', '),
                  'as fixed effects')
tmp      <- paste('Using following covariates in the final model:', modelStr)



#saveRDS(DESIGN,"data/preprocessed/design.RDS")
saveRDS(NEW.COUNTS,"data/preprocessed/NEWCOUNTS.RDS")
saveRDS(rownames(NEW.COUNTS),"data/preprocessed/brainGenes.RDS")
saveRDS(expr,"data/preprocessed/expr.RDS")
saveRDS(CQN.GENE_EXPRESSION,"data/preprocessed/CQNNORMALIZED.RDS")
saveRDS(COVARIATES,"data/preprocessed/COVARIATES.RDS")
saveRDS(GENE.PARAM,"data/preprocessed/GENEPARAM.RDS")
saveRDS(postAdjustCovars,"data/preprocessed/POSTCOVAR.RDS")


