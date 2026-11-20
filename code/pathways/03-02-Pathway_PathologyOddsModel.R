rm(list = ls())

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
doParallel::registerDoParallel(10)

source("code/00-02-functions2.R")



# Loading covariates, corrected miRNA, and corrected Pathways
EXPRS <- readRDS("output/DEPaths_bulk_update/ROSMAP_DLPFC_ADvsResilient_Residuals.RDS")
covs  <- readRDS(paste0("data/preprocessed/clinicalTotal.RDS"))




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




DE.genes <- read.csv(paste0("output/DEPaths_bulk_update/ROSMAP_DLPFC_ADvsResilient_diffPathways.csv"), 
                     stringsAsFactors = F, row.names = 1)

net.list <- list()
de.list  <- list()
braak.list  <- list()
cerad.list  <- list()
sum.list    <- list()
inter.list  <- list()
binom.list  <- list()
pom.cog     <- list()
pom.braak   <- list()
pom.cerad   <- list()


proportional.summary <- function(Y = NA, X1 = NA , DAT = NA, mir.name = NA){
                        
                        formz <- paste0(eval(Y), " ~ " , eval(X1))
                        fit   <- vglm(as.formula(formz) ,
                                      family=cumulative(parallel=TRUE),
                                      data=DAT)
                        
                        temp <- summaryvglm(fit)
                        temp <- temp@coef3
                        temp <- as.data.frame(temp)
                        temp <- rownames_to_column(temp,var =  "term")
                        temp$gene <- mir.name
                        temp$R2lat <- R2latvar(fit)
                        return(temp)
                        
                    }


all.ad  <- covs[covs$cogdx.x %in% c(1,2,4),]
#all.ad  <- covs[covs$RESIL %in% c("OTHER") & covs$cogdx.x %in% c(1,2,4),]
gene.AD <- EXPRS[,as.character(all.ad$SampleID)]


gene.names <- rownames(DE.genes)
iter <- length(gene.names)

cog.paths <- foreach(j = 1:iter,.combine = 'rbind') %dopar% {
    
    tested_gene0 <- gene.names[j]
    tested_gene  <- gene.AD[tested_gene0,]
    this.gene    <- data.frame("tested_gene" = tested_gene,
                                 "cog"  = factor(all.ad$cogdx.y, ordered = T))
    
    prop.fit        <- proportional.summary("cog","tested_gene",          
                                            this.gene, mir.name = gene.names[j])
}


cog.paths0 <- clean.write.podds2(cog.paths,"output/Bulk_paths_prop_odds/cog_paths.csv")


tangle.paths <- foreach(j = 1:iter,.combine = 'rbind') %dopar% {
    
  tested_gene0 <- gene.names[j]
  tested_gene  <- gene.AD[tested_gene0,]
    this.gene      <- data.frame("tested_gene" = tested_gene,
                                 "braak"  = factor(all.ad$braaksc,ordered = T))
    prop.fit       <- proportional.summary("braak","tested_gene",
                                           this.gene, mir.name = gene.names[j])
}
tangle.paths0 <- clean.write.podds2(tangle.paths,"output/Bulk_paths_prop_odds/tangle_paths.csv")


plaque.paths <- foreach(j = 1:iter,.combine = 'rbind') %dopar% {
    
    
    tested_gene0 <- gene.names[j]
    tested_gene  <- gene.AD[tested_gene0,]
    this.gene      <- data.frame("tested_gene" = tested_gene,
                                 "cerad"  = factor(4- all.ad$ceradsc, ordered = T))
    prop.fit       <- proportional.summary("cerad","tested_gene",
                                           this.gene, mir.name =gene.names[j])
}
plaque.paths0 <- clean.write.podds2(plaque.paths,"output/Bulk_paths_prop_odds/plaque_paths.csv")



all.paths <- foreach(j = 1:iter,.combine = 'rbind') %dopar% {
    tested_gene0 <- gene.names[j]
    tested_gene  <- gene.AD[tested_gene0,]
    this.gene      <- data.frame("tested_gene" = tested_gene,
                                 "cerad"  = factor(4- all.ad$ceradsc,ordered = T),
                                 "braak"  = factor(all.ad$braaksc, ordered = T),
                                 "cog"  = factor(all.ad$cogdx.y, ordered = T))
    prop.fit       <- proportional.summary("cog","tested_gene + cerad + braak",
                                           this.gene, mir.name = gene.names[j])
}


#all.paths <- clean.write.podds(all.paths,"Output11/NIH_resiliome/all_cog_corrected_all_genes.csv")
all.paths0 <- clean.write.podds2(all.paths,"output/Bulk_paths_prop_odds/cog_corrected_paths.csv")

# 


# sessionInfo()
# > sessionInfo()
# R version 4.4.1 (2024-06-14)
# Platform: aarch64-apple-darwin20
# Running under: macOS 15.7.1

# Matrix products: default
# BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# time zone: America/New_York
# tzcode source: internal

# attached base packages:
# [1] splines   stats4    stats     graphics  grDevices utils     datasets 
# [8] methods   base     

# other attached packages:
#  [1] qvalue_2.38.0           PanomiR_1.10.0          foreach_1.5.2          
#  [4] VGAM_1.1-6              broom_1.0.7             org.Hs.eg.db_3.20.0    
#  [7] AnnotationDbi_1.68.0    IRanges_2.40.1          S4Vectors_0.44.0       
# [10] Biobase_2.66.0          BiocGenerics_0.52.0     clusterProfiler_4.14.4 
# [13] CovariateAnalysis_0.1.0 lubridate_1.9.4         forcats_1.0.0          
# [16] tibble_3.2.1            tidyverse_2.0.0         psych_1.8.12           
# [19] limma_3.62.1            scales_1.3.0            ggplot2_3.5.1          
# [22] stringr_1.5.1           purrr_1.0.2             tidyr_1.3.1            
# [25] readr_2.1.5             dplyr_1.1.4            

# loaded via a namespace (and not attached):
#  [1] DBI_1.2.3               mnormt_2.1.1            gson_0.1.0             
#  [4] readxl_1.4.3            rlang_1.1.4             magrittr_2.0.3         
#  [7] DOSE_4.0.0              compiler_4.4.1          RSQLite_2.3.9          
# [10] png_0.1-8               systemfonts_1.1.0       vctrs_0.6.5            
# [13] reshape2_1.4.4          pkgconfig_2.0.3         crayon_1.5.3           
# [16] fastmap_1.2.0           backports_1.5.0         XVector_0.46.0         
# [19] labeling_0.4.3          utf8_1.2.4              tzdb_0.4.0             
# [22] enrichplot_1.26.3       UCSC.utils_1.2.0        ragg_1.3.3             
# [25] bit_4.5.0.1             zlibbioc_1.52.0         cachem_1.1.0           
# [28] aplot_0.2.3             GenomeInfoDb_1.42.1     jsonlite_1.8.9         
# [31] blob_1.2.4              BiocParallel_1.40.0     parallel_4.4.1         
# [34] R6_2.5.1                stringi_1.8.4           RColorBrewer_1.1-3     
# [37] cellranger_1.1.0        GOSemSim_2.32.0         iterators_1.0.14       
# [40] Rcpp_1.0.13-1           ggtangle_0.0.5          R.utils_2.12.3         
# [43] Matrix_1.7-1            igraph_2.1.2            timechange_0.3.0       
# [46] tidyselect_1.2.1        doParallel_1.0.17       codetools_0.2-20       
# [49] lattice_0.22-6          plyr_1.8.9              treeio_1.30.0          
# [52] withr_3.0.2             KEGGREST_1.46.0         foreign_0.8-87         
# [55] gridGraphics_0.5-1      Biostrings_2.74.0       pillar_1.9.0           
# [58] ggtree_3.14.0           ggfun_0.1.8             generics_0.1.3         
# [61] vroom_1.6.5             hms_1.1.3               munsell_0.5.1          
# [64] tidytree_0.4.6          glue_1.8.0              lazyeval_0.2.2         
# [67] tools_4.4.1             data.table_1.16.4       fgsea_1.32.0           
# [70] fs_1.6.5                fastmatch_1.1-4         cowplot_1.1.3          
# [73] grid_4.4.1              ape_5.8                 colorspace_2.1-1       
# [76] nlme_3.1-166            GenomeInfoDbData_1.2.13 patchwork_1.3.0        
# [79] cli_3.6.3               textshaping_0.4.1       fansi_1.0.6            
# [82] gtable_0.3.6            R.methodsS3_1.8.2       yulab.utils_0.1.8      
# [85] digest_0.6.37           ggrepel_0.9.6           ggplotify_0.1.2        
# [88] farver_2.1.2            memoise_2.0.1           R.oo_1.27.0            
# [91] lifecycle_1.0.4         httr_1.4.7              GO.db_3.20.0           
# [94] statmod_1.5.0           bit64_4.5.2            