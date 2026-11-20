rm(list = ls())

args = commandArgs(trailingOnly=TRUE)
if(is.na(args[1])){
    setup <- ""
    warning("The experimental design specifications must be set.")
} else {
    setup <- args[1]
}

print("Differential miRNA Expression")

options(StringAsFactors = F)



library(limma)
library(dplyr)
library(psych)
library(tidyverse)
library(CovariateAnalysis)

data.dir        <- "data/preprocessed/"
output.dir      <- paste0("output", setup,"/")
output.dir2     <- paste0("output", setup,"/deMirs/")
if(!dir.exists(output.dir2)){
  dir.create(output.dir2)
}


input.dir       <- paste0(output.dir,"preprocessed/")
# The first output dir is the current one for pathway summaries. Keep if things fall apart


clinical.Tot   <- readRDS(file = "data/preprocessed/mirClinical.RDS")
#RESIDUAL.paths <- readRDS(paste0(input.dir,"resils.RDS"))
miRExp <- readRDS("data/preprocessed/miRNAExpression.RDS")


#Defining resilience

cogdx.resils    <- c(1)
cogdx.AD        <- c(4)
cogdx.normal    <- c(1)

braak.resils    <- c(3,4,5,6)
braak.AD        <- c(3,4,5,6)

plq.resils      <- c(1,2)
plq.AD          <- c(1,2)
age.cut         <- 80
age.cut.AD      <- 80
age.cut.norm    <- 0

clinical.Tot   <- mutate(clinical.Tot, RESIL = "OTHER")

resilients     <- clinical.Tot %>% dplyr::filter(., 
                                                 cogdx %in% cogdx.resils,
                                                 braaksc %in% braak.resils,
                                                 ceradsc %in% plq.resils,
                                                 age_death >= age.cut) %>% 
  mutate(.,RESIL = "Resilient")#%>% nrow

AD     <- clinical.Tot %>% filter(., 
                                          cogdx %in% cogdx.AD,
                                          braaksc %in% braak.AD,
                                          ceradsc %in% plq.AD,
                                          age_death >= age.cut.AD) %>% 
  mutate(.,RESIL = "AD")

normals     <- clinical.Tot %>% filter(., 
                                       cogdx %in% cogdx.normal,
                                       braaksc %in% c(0,1,2),
                                       ceradsc %in% c(3,4),
                                       age_death >= age.cut.norm) %>% 
  mutate(.,RESIL = "Control")

presym     <- clinical.Tot %>% filter(., 
                                      cogdx %in% c(2),
                                      braaksc %in% braak.resils,
                                      ceradsc %in% plq.resils) %>% 
  mutate(.,RESIL = "Presym")

clinical.resil  <- rbind(resilients,AD,normals,presym)

clinical.Tot[clinical.Tot$mirna_id %in% resilients$mirna_id,]$RESIL <- "Resilient"
clinical.Tot[clinical.Tot$mirna_id %in% AD$mirna_id,]$RESIL <- "AD"
clinical.Tot[clinical.Tot$mirna_id %in% normals$mirna_id,   ]$RESIL <- "Control"
clinical.Tot[clinical.Tot$mirna_id %in% presym$mirna_id,   ]$RESIL <- "Presym"


clinical.Tot$Diagnosis <- clinical.Tot$RESIL
clinical.Tot$Diagnosis <- as.factor(clinical.Tot$Diagnosis)


COVARIATES           <- clinical.Tot  
rownames(COVARIATES) <- COVARIATES$mirna_id
#colnames(COVARIATES) <- gsub("\\.x", "",colnames(COVARIATES))

table(COVARIATES$Diagnosis)

COVARIATES <- as.data.frame(COVARIATES)
rownames(COVARIATES) <- COVARIATES$mirna_id


FactorCovariates <- c("msex","apoe_genotype", "braaksc","ceradsc",                    
                      "cogdx","plate")

COVARIATES[,FactorCovariates] = lapply(COVARIATES[,FactorCovariates], factor)





corrections <- c("msex", "RIN", "age_death", "pmi","apoe_genotype")


colSums(miRExp)

miRExp <- miRExp[,rownames(COVARIATES)]


library(qvalue)
library(EnhancedVolcano)



    DESIGN         <- getDesignMatrix(COVARIATES[,c("Diagnosis",corrections),
                                           drop = F], Intercept = F)
    DESIGN$design  <- DESIGN$design[,linColumnFinder(DESIGN$design)$indepCols]
    FIT            <- lmFit(miRExp, DESIGN$design)
    contrast  <- makeContrasts(contrasts=c("DiagnosisAD-DiagnosisResilient"),
                               levels = colnames(FIT$coefficients))
    FIT.CONTR <- contrasts.fit(FIT, contrasts=contrast)
    FIT.CONTR <- eBayes(FIT.CONTR)
    tT        <- topTable(FIT.CONTR, adjust="fdr", sort.by="p", number=Inf)
    tT$qvalue <- qvalue(tT$P.Value)$qval
    
    write.csv(tT,paste0(output.dir2,"deMirsADvsResilient.csv"))
    
    
    res1 <- tT
    res1$miRNA <- rownames(res1)
    
    sel_labs <- (res1[   res1$adj.P.Val < 0.05 , ]$miRNA)
    
    sel_labs <- sel_labs[sel_labs!= ""]
    keyvals <- ifelse(
      (res1$adj.P.Val > 0.05), 'grey60',
      ifelse((res1$logFC) < 0, 'royalblue',
             'red3'))
    
    keyvals[is.na(keyvals)] <- 'grey60'
    names(keyvals)[keyvals == 'red3'] <- 'Up in AD'
    names(keyvals)[keyvals == 'grey60'] <- 'n.s.'
    names(keyvals)[keyvals == 'royalblue'] <- 'Down in AD'
    
    


    # p-cut corresponding to FDR 0.05

    pcut <- sum(res1$adj.P.Val < 0.05)
    if(pcut >0)
      {pcut <- mean(c((res1$P.Value[pcut]),(res1$P.Value[pcut+1]) ))
      }else
    {pcut <- 0.05/nrow(res1)}
      

    p_25_cut <- sum(res1$adj.P.Val < 0.2)
    if(p_25_cut >0)
    {p_25_cut <- mean(c(res1$P.Value[p_25_cut],res1$P.Value[p_25_cut+1]))
    }else
    {p_25_cut <- 0.2/nrow(res1)}


    # This only works because AD vs control has 
    # min_p_ad_res <- min(res1$P.Value)
    # min_p_ad_res <- ceiling(-log10(min_p_ad_res))
    min_p_ad_res <- 9
  

    p1 <-   EnhancedVolcano(res1,
                            lab = res1$miRNA, 
                            selectLab = sel_labs,
                            x = "logFC",
                            y = "P.Value",
                            xlim = c(-.7,.7),
                            ylim = c(0, min_p_ad_res+1),
                            caption = NULL,
                            xlab = bquote(~Log[2]~ "fold change"),
                            ylab = bquote(~-Log[10]~""~italic(P)),
                            title = NULL,
                            subtitle = "",
                            pCutoff = pcut,   
                            pointSize = 3.0,
                            labSize = 5,
                            labCol = 'black',
                            labFace = 'plain',
                            drawConnectors = T, arrowheads = F,
                            boxedLabels = FALSE,
                            colCustom = keyvals,
                            widthConnectors = 1.2,
                            colConnectors = 'black',
                            legendPosition = "none",
                            axisLabSize = 30,
                            legendLabSize = 20,
                            legendIconSize = 10
    ) 
  
    ggsave("output/figures/AD_vs_Res_volcano.pdf", p1,width = 6,height = 8)
    ggsave("output/figures_paper/Figure 1B AD_vs_Res_volcano.jpg", p1,width = 6,height = 8,dpi = 600)

    
    
    
    
    

    DESIGN         <- getDesignMatrix(COVARIATES[,c("Diagnosis",corrections),
                                                 drop = F], Intercept = F)
    DESIGN$design  <- DESIGN$design[,linColumnFinder(DESIGN$design)$indepCols]
    FIT            <- lmFit(miRExp, DESIGN$design)
    contrast  <- makeContrasts(contrasts=c("DiagnosisAD-DiagnosisControl"),
                               levels = colnames(FIT$coefficients))
    FIT.CONTR <- contrasts.fit(FIT, contrasts=contrast)
    FIT.CONTR <- eBayes(FIT.CONTR)
    tT        <- topTable(FIT.CONTR, adjust="fdr", sort.by="p", number=Inf)
    tT$qvalue <- qvalue(tT$P.Value)$qval
    write.csv(tT,paste0(output.dir2,"deMirsADvsControl.csv"))

    
    res1 <- tT
    res1$miRNA <- rownames(res1)
    
    sel_labs <- (res1[   res1$adj.P.Val < 0.05 , ]$miRNA)
    
    sel_labs <- sel_labs[sel_labs!= ""]
    keyvals <- ifelse(
      (res1$adj.P.Val > 0.05), 'grey60',
      ifelse((res1$logFC) < 0, 'royalblue',
             'red3'))
    
    keyvals[is.na(keyvals)] <- 'grey60'
    names(keyvals)[keyvals == 'red3'] <- 'Up in AD'
    names(keyvals)[keyvals == 'grey60'] <- 'n.s.'
    names(keyvals)[keyvals == 'royalblue'] <- 'Down in AD'
    
    
    # p-cut corresponding to FDR 0.05

    pcut <- sum(res1$adj.P.Val < 0.05)
    if(pcut >0)
      {pcut <- mean(c((res1$P.Value[pcut]),(res1$P.Value[pcut+1]) ))
      }else
    {pcut <- 0.05/nrow(res1)}
      

    p_25_cut <- sum(res1$adj.P.Val < 0.2)
    if(p_25_cut >0)
    {p_25_cut <- mean(c(res1$P.Value[p_25_cut],res1$P.Value[p_25_cut+1]))
    }else
    {p_25_cut <- 0.2/nrow(res1)}


    # This only works because AD vs control has 
    # min_p_ad_res <- min(res1$P.Value)
    # min_p_ad_res <- ceiling(-log10(min_p_ad_res))
    min_p_ad_res <- 9


  

    p1 <-   EnhancedVolcano(res1,
                            lab = res1$miRNA, 
                            selectLab = sel_labs,
                            x = "logFC",
                            y = "P.Value",
                            xlim = c(-.7,.7),
                            ylim = c(0, min_p_ad_res+1),
                            caption = NULL,
                            xlab = bquote(~Log[2]~ "fold change"),
                            ylab = bquote(~-Log[10]~""~italic(P)),
                            title = NULL,
                            subtitle = "",
                            pCutoff = pcut,   
                            pointSize = 3.0,
                            labSize = 5,
                            labCol = 'black',
                            labFace = 'plain',
                            drawConnectors = T, arrowheads = F,
                            boxedLabels = FALSE,
                            colCustom = keyvals,
                            widthConnectors = 1.2,
                            colConnectors = 'black',
                            legendPosition = "none",
                            axisLabSize = 30,
                            legendLabSize = 20,
                            legendIconSize = 10
    )

    ggsave("output/figures/AD_vs_Ctrl_volcano.pdf", p1,width = 6,height = 8)
    ggsave("output/figures_paper/1B AD_vs_Ctrl_volcano.jpg", p1,width = 6,height = 8,dpi = 600)

    
    DESIGN         <- getDesignMatrix(COVARIATES[,c("Diagnosis",corrections),
                                                 drop = F], Intercept = F)
    DESIGN$design  <- DESIGN$design[,linColumnFinder(DESIGN$design)$indepCols]
    FIT            <- lmFit(miRExp, DESIGN$design)
    contrast  <- makeContrasts(contrasts=c("DiagnosisResilient-DiagnosisControl"),
                               levels = colnames(FIT$coefficients))
    FIT.CONTR <- contrasts.fit(FIT, contrasts=contrast)
    FIT.CONTR <- eBayes(FIT.CONTR)
    tT        <- topTable(FIT.CONTR, adjust="fdr", sort.by="p", number=Inf)
    
    write.csv(tT,paste0(output.dir2,"deMirsResilientvsNormal.csv"))
    
    

    res1 <- tT
    res1$miRNA <- rownames(res1)
    
    sel_labs <- (res1[   res1$adj.P.Val < 0.05 , ]$miRNA)
    
    sel_labs <- sel_labs[sel_labs!= ""]
    keyvals <- ifelse(
      (res1$adj.P.Val > 0.05), 'grey60',
      ifelse((res1$logFC) < 0, 'royalblue',
             'red3'))
    
    keyvals[is.na(keyvals)] <- 'grey60'
    names(keyvals)[keyvals == 'red3'] <- 'Up in Resilient'
    names(keyvals)[keyvals == 'grey60'] <- 'n.s.'
    names(keyvals)[keyvals == 'royalblue'] <- 'Down in Resilient'
    
     # p-cut corresponding to FDR 0.05

    pcut <- sum(res1$adj.P.Val < 0.05)
    if(pcut >0)
      {pcut <- mean(c((res1$P.Value[pcut]),(res1$P.Value[pcut+1]) ))
      }else
    {pcut <- 0.05/nrow(res1)}
      

    p_25_cut <- sum(res1$adj.P.Val < 0.2)
    if(p_25_cut >0)
    {p_25_cut <- mean(c(res1$P.Value[p_25_cut],res1$P.Value[p_25_cut+1]))
    }else
    {p_25_cut <- 0.2/nrow(res1)}


    # This only works because AD vs control has 
    # min_p_ad_res <- min(res1$P.Value)
    # min_p_ad_res <- ceiling(-log10(min_p_ad_res))
    min_p_ad_res <- 9


  

    p1 <-   EnhancedVolcano(res1,
                            lab = res1$miRNA, 
                            selectLab = sel_labs,
                            x = "logFC",
                            y = "P.Value",
                            xlim = c(-.7,.7),
                            ylim = c(0, min_p_ad_res+1),
                            caption = NULL,
                            xlab = bquote(~Log[2]~ "fold change"),
                            ylab = bquote(~-Log[10]~""~italic(P)),
                            title = NULL,
                            subtitle = "",
                            pCutoff = pcut,   
                            pointSize = 3.0,
                            labSize = 5,
                            labCol = 'black',
                            labFace = 'plain',
                            drawConnectors = T, arrowheads = F,
                            boxedLabels = FALSE,
                            colCustom = keyvals,
                            widthConnectors = 1.2,
                            colConnectors = 'black',
                            legendPosition = "none",
                            axisLabSize = 30,
                            legendLabSize = 20,
                            legendIconSize = 10
    )
    ggsave("output/figures/Res_vs_Ctrl_volcano.pdf", p1,width = 6,height = 8)
    ggsave("output/figures_paper/Figure 1B Res_vs_Ctrl_volcano.jpg", p1,width = 6,height = 8, dpi = 600)

    
    



    
    
    DESIGN         <- getDesignMatrix(COVARIATES[,c("Diagnosis",corrections),
                                                 drop = F], Intercept = F)
    DESIGN$design  <- DESIGN$design[,linColumnFinder(DESIGN$design)$indepCols]
    FIT            <- lmFit(miRExp, DESIGN$design)
    contrast  <- makeContrasts(contrasts=c("DiagnosisAD-DiagnosisPresym"),
                               levels = colnames(FIT$coefficients))
    FIT.CONTR <- contrasts.fit(FIT, contrasts=contrast)
    FIT.CONTR <- eBayes(FIT.CONTR)
    tT        <- topTable(FIT.CONTR, adjust="fdr", sort.by="p", number=Inf)
    tT$qvalue <- qvalue(tT$P.Value)$qval
    write.csv(tT,paste0(output.dir2,"deMirsADvsPresym.csv"))
    
    
    
    
    DESIGN         <- getDesignMatrix(COVARIATES[,c("Diagnosis",corrections),
                                                 drop = F], Intercept = F)
    DESIGN$design  <- DESIGN$design[,linColumnFinder(DESIGN$design)$indepCols]
    FIT            <- lmFit(miRExp, DESIGN$design)
    contrast  <- makeContrasts(contrasts=c("DiagnosisResilient-DiagnosisPresym"),
                               levels = colnames(FIT$coefficients))
    FIT.CONTR <- contrasts.fit(FIT, contrasts=contrast)
    FIT.CONTR <- eBayes(FIT.CONTR)
    tT        <- topTable(FIT.CONTR, adjust="fdr", sort.by="p", number=Inf)
    tT$qvalue <- qvalue(tT$P.Value)$qval
    write.csv(tT,paste0(output.dir2,"deMirsResilientvsPresym.csv"))
    
    
    
    DESIGN         <- getDesignMatrix(COVARIATES[,c("Diagnosis",corrections),
                                                 drop = F], Intercept = F)
    DESIGN$design  <- DESIGN$design[,linColumnFinder(DESIGN$design)$indepCols]
    FIT            <- lmFit(miRExp, DESIGN$design)
    contrast  <- makeContrasts(contrasts=c("DiagnosisPresym-DiagnosisControl"),
                               levels = colnames(FIT$coefficients))
    FIT.CONTR <- contrasts.fit(FIT, contrasts=contrast)
    FIT.CONTR <- eBayes(FIT.CONTR)
    tT        <- topTable(FIT.CONTR, adjust="fdr", sort.by="p", number=Inf)
    tT$qvalue <- qvalue(tT$P.Value)$qval
    write.csv(tT,paste0(output.dir2,"deMirsPresymvsControl.csv"))
    
    
    
    DESIGN         <- getDesignMatrix(COVARIATES[,c("Diagnosis",corrections),
                                                 drop = F], Intercept = F)
    DESIGN$design  <- DESIGN$design[,linColumnFinder(DESIGN$design)$indepCols]
    FIT            <- lmFit(miRExp, DESIGN$design)
    contrast  <- makeContrasts(contrasts=c("DiagnosisAD-(DiagnosisControl+DiagnosisPresym+DiagnosisResilient)/3"),
                               levels = colnames(FIT$coefficients))
    FIT.CONTR <- contrasts.fit(FIT, contrasts=contrast)
    FIT.CONTR <- eBayes(FIT.CONTR)
    tT        <- topTable(FIT.CONTR, adjust="fdr", sort.by="p", number=Inf)
    tT$qvalue <- qvalue(tT$P.Value)$qval
    write.csv(tT,paste0(output.dir2,"deMirsADvsNonAD.csv"))





    COVARIATES$Diagnosis <- relevel(COVARIATES$Diagnosis, ref="OTHER")
    COVARIATES$sex <- factor(as.character(COVARIATES$msex))
    levels(COVARIATES$sex) <-  c("female","male")





    design <- model.matrix(~ sex * Diagnosis + age_death + pmi + apoe_genotype + RIN, data=COVARIATES)
        colnames(design)
    colnames(design) <- make.names(colnames(design))


COVARIATES$Diagnosis
    
    fit <- lmFit(miRExp, design)
    
    fit <- eBayes(fit)

    # check coefficient names

    # interaction test (does AD effect differ by sex?)
    topTable(fit, coef="DiagnosisAD", number=Inf)

    # create contrasts for disease effect in each sex
    cont <- makeContrasts(
      AD_vs_Ctrl_in_Females   = DiagnosisAD - DiagnosisControl,                          # Female is reference
      AD_vs_Resilient_in_Females = DiagnosisAD - DiagnosisResilient,
      AD_vs_Ctrl_in_Males     = (DiagnosisAD + sexmale.DiagnosisAD) - (DiagnosisControl + sexmale.DiagnosisControl),
      AD_vs_Resilient_in_Males = (DiagnosisAD + sexmale.DiagnosisAD) - (DiagnosisResilient + sexmale.DiagnosisResilient),
      Male_vs_Female_in_AD = sexmale + sexmale.DiagnosisAD,
      Male_vs_Female_in_Control = sexmale + sexmale.DiagnosisControl,
      Male_vs_Female_in_Resilient = sexmale + sexmale.DiagnosisResilient,
      Male_vs_Female_in_AD_vs_Control = sexmale.DiagnosisAD - sexmale.DiagnosisControl,
      Male_vs_Female_in_AD_vs_Res = sexmale.DiagnosisAD - sexmale.DiagnosisResilient,
      levels = design
    )
    fit2 <- contrasts.fit(fit, cont)
    fit2 <- eBayes(fit2)

    my_constrasts <- colnames(cont)

    for(test in my_constrasts){
      out_file <- paste0(output.dir2,"deMirs_sexSpec_",test,".csv")
      tt <- topTable(fit2, coef=test, number=Inf)
      tt$qvalue <- qvalue(tt$P.Value)$qval
      write.csv(tt, out_file)
    }
