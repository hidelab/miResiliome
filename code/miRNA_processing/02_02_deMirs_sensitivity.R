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
output.dir2     <- paste0("output", setup,"/deMirs_sensitivity/")
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





corrections_with_apoe <- c("msex", "RIN", "age_death", "pmi","apoe_genotype")
corrections_without_apoe <- c("msex", "RIN", "age_death", "pmi")


colSums(miRExp)

miRExp <- miRExp[,rownames(COVARIATES)]


library(qvalue)
library(EnhancedVolcano)



    DESIGN         <- getDesignMatrix(COVARIATES[,c("Diagnosis",corrections_without_apoe),
                                           drop = F], Intercept = F)
    DESIGN$design  <- DESIGN$design[,linColumnFinder(DESIGN$design)$indepCols]
    FIT            <- lmFit(miRExp, DESIGN$design)
    contrast  <- makeContrasts(contrasts=c("DiagnosisAD-DiagnosisResilient"),
                               levels = colnames(FIT$coefficients))
    FIT.CONTR <- contrasts.fit(FIT, contrasts=contrast)
    FIT.CONTR <- eBayes(FIT.CONTR)
    tT        <- topTable(FIT.CONTR, adjust="fdr", sort.by="p", number=Inf)
    tT$qvalue <- qvalue(tT$P.Value)$qval
    
    write.csv(tT,paste0(output.dir2,"deMirsADvsResilient_withoutAPOE.csv"))



    
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
  
    ggsave("output/figures/sensitivity_AD_vs_Res_volcano_NoAPOE.pdf", p1,width = 6,height = 8)
    ggsave("output/figures_paper/sensitivity_Figure 1B AD_vs_Res_volcano_NoAPOE.jpg", p1,width = 6,height = 8,dpi = 600)

    
    










cellprops <- readRDS("./data/theta_final_celltype.rds")

cellprops <- rownames_to_column(as.data.frame(cellprops),"Sample_ID")
colnames(cellprops)
head(cellprops)
COVARIATES_extended <- dplyr::inner_join(COVARIATES,cellprops)




FactorCovariates <- c("msex","apoe_genotype", "braaksc","ceradsc",                    
                      "cogdx","plate")

COVARIATES_extended[,FactorCovariates] = lapply(COVARIATES_extended[,FactorCovariates], factor)
rownames(COVARIATES_extended) <- COVARIATES_extended$mirna_id




corrections_with <- c("msex", "RIN", "age_death", "pmi","apoe_genotype","Ex")
corrections_without <- c("msex", "RIN", "age_death", "pmi","apoe_genotype")


colSums(miRExp)

miRExp_extended <- miRExp[,rownames(COVARIATES_extended)]




    DESIGN         <- getDesignMatrix(COVARIATES_extended[,c("Diagnosis",corrections_with),
                                           drop = F], Intercept = F)
    DESIGN$design  <- DESIGN$design[,linColumnFinder(DESIGN$design)$indepCols]
    FIT            <- lmFit(miRExp_extended, DESIGN$design)
    contrast  <- makeContrasts(contrasts=c("DiagnosisAD-DiagnosisResilient"),
                               levels = colnames(FIT$coefficients))
    FIT.CONTR <- contrasts.fit(FIT, contrasts=contrast)
    FIT.CONTR <- eBayes(FIT.CONTR)
    tT        <- topTable(FIT.CONTR, adjust="fdr", sort.by="p", number=Inf)
    tT$qvalue <- qvalue(tT$P.Value)$qval
    
    write.csv(tT,paste0(output.dir2,"deMirsADvsResilient_ExNeurons.csv"))



    
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
  
    ggsave("output/figures/sensitivity_AD_vs_Res_volcano_ExNeurons.pdf", p1,width = 6,height = 8)
    ggsave("output/figures_paper/sensitivity_Figure 1B AD_vs_Res_volcano_ExNeurons.jpg", p1,width = 6,height = 8,dpi = 600)

    
    




corrections_without <- c("msex", "RIN", "age_death", "pmi","apoe_genotype")

table(COVARIATES_extended$Diagnosis)



    DESIGN         <- getDesignMatrix(COVARIATES_extended[,c("Diagnosis",corrections_without),
                                           drop = F], Intercept = F)
    DESIGN$design  <- DESIGN$design[,linColumnFinder(DESIGN$design)$indepCols]
    FIT            <- lmFit(miRExp_extended, DESIGN$design)
    contrast  <- makeContrasts(contrasts=c("DiagnosisAD-DiagnosisResilient"),
                               levels = colnames(FIT$coefficients))
    FIT.CONTR <- contrasts.fit(FIT, contrasts=contrast)
    FIT.CONTR <- eBayes(FIT.CONTR)
    tT        <- topTable(FIT.CONTR, adjust="fdr", sort.by="p", number=Inf)
    tT$qvalue <- qvalue(tT$P.Value)$qval
    
    write.csv(tT,paste0(output.dir2,"deMirsADvsResilient_NoExNeurons.csv"))



    
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
  
    ggsave("output/figures/sensitivity_AD_vs_Res_volcano_NoExNeurons.pdf", p1,width = 6,height = 8)
    ggsave("output/figures_paper/sensitivity_Figure 1B AD_vs_Res_volcano_NoExNeurons.jpg", p1,width = 6,height = 8,dpi = 600)




# Excitatory neuron proportion (Ex) by Diagnosis

ex_dat <- COVARIATES_extended %>%
  dplyr::filter(Diagnosis %in% c("AD", "Resilient", "Control"))
ex_dat$Diagnosis <- factor(ex_dat$Diagnosis, levels = c("AD", "Resilient", "Control"))

tapply(ex_dat$Ex, ex_dat$Diagnosis, summary)


# Wilcoxon rank-sum test, AD vs each group, with BH correction across the comparisons
ex_groups <- c("Resilient", "Control")
ex_contrasts <- data.frame(
  contrast  = paste0("AD vs ", ex_groups),
  median_AD = median(ex_dat$Ex[ex_dat$Diagnosis == "AD"]),
  median_grp = sapply(ex_groups, function(g) median(ex_dat$Ex[ex_dat$Diagnosis == g])),
  P.Value   = sapply(ex_groups, function(g)
    wilcox.test(ex_dat$Ex[ex_dat$Diagnosis == "AD"],
                ex_dat$Ex[ex_dat$Diagnosis == g])$p.value),
  row.names = NULL
)
ex_contrasts$adj.P.Val <- p.adjust(ex_contrasts$P.Value, method = "BH")
ex_contrasts


# Box plot with significance bars (Wilcoxon p-values, raw and BH-adjusted)
library(ggsignif)

diag_cols <- c(AD = "#D55E00", Resilient = "#009E73", Control = "#0072B2")

ex_comparisons <- lapply(ex_groups, function(g) c("AD", g))
ex_p_labels    <- paste0("Wilcoxon p = ", formatC(ex_contrasts$P.Value,   format = "f", digits = 4),
                         "\nBH adj. p = ", formatC(ex_contrasts$adj.P.Val, format = "f", digits = 4))
ex_ymax        <- max(ex_dat$Ex, na.rm = TRUE)
ex_yrange      <- diff(range(ex_dat$Ex, na.rm = TRUE))
ex_ypos        <- ex_ymax + ex_yrange * c(0.08, 0.30)

    p_ex <- ggplot(ex_dat, aes(x = Diagnosis, y = Ex, fill = Diagnosis)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(aes(color = Diagnosis), width = 0.2, alpha = 0.7, size = 1.8) +
      geom_signif(comparisons = ex_comparisons,
                  annotations = ex_p_labels,
                  y_position  = ex_ypos,
                  tip_length  = 0.01,
                  textsize    = 6) +
      scale_fill_manual(values = diag_cols) +
      scale_color_manual(values = diag_cols) +
      scale_y_continuous(breaks = pretty(range(ex_dat$Ex, na.rm = TRUE)),
                         expand = expansion(mult = c(0.05, 0.20))) +
      ylab("Excitatory neuron proportion (Ex)") +
      xlab(NULL) +
      theme_bw(base_size = 20) +
      theme(legend.position = "none",
            axis.text.x     = element_text(size = 20, color = "black"),
            axis.text.y     = element_text(size = 16, color = "black"),
            axis.title.y    = element_text(size = 20))

    p_ex
    ggsave("output/figures/sensitivity_Ex_by_Diagnosis_boxplot.pdf", p_ex, width = 7, height = 7)

    
    





