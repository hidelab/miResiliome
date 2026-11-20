# MicroRNA correlates of resilience to Alzheimer's disease

This repository contains the code to reproduce the analyses in:

> Naderi P, Matai L, Adewale Q, Castanho I, Rodrigues A, Mavrikaki M, Vlachos IS, Slack FJ, Hide W.
> **MicroRNA correlates of resilience to Alzheimer's disease identify candidate therapeutic targets for neuroprotection.**
> *bioRxiv* [Preprint]. 2026 Aug 10:2026.08.04.737285.
> doi: [10.64898/2026.08.04.737285](https://doi.org/10.64898/2026.08.04.737285).
> PMID: 42845305; PMCID: PMC13640498.

## Repository structure

| Folder | Contents |
| --- | --- |
| `code/` | Analysis scripts (`00-02-functions2.R` holds shared helper functions) |
| `code/miRNA_processing/` | miRNA reannotation, normalization, differential expression, pathology/cognition models, co-expression, and miRNA–mRNA integration |
| `code/mRNA_processing/` | mRNA expression normalization and resilience differential expression |
| `code/pathways/` | Pathway analysis, pathology/cognition models, and PanomiR analysis |
| `setting/` | Parameter files used by the analysis scripts |
| `data/` | Placeholder for input data (not included; see below) |
| `output/` | Analysis outputs (only `Bulk_paths_prop_odds/` is included) |

Scripts are numbered in the order they are meant to run.

## Data

The input data is not distributed with this repository.

Gene and miRNA expression datasets used in this study are publicly available through the AMP-AD data
portal hosted on Synapse under their original ROSMAP accession code
[syn8456629](https://www.synapse.org/#!Synapse:syn8456629) under controlled use conditions subject to
privacy regulations. Bulk RNAseq data and miRNA expression data are accessible under the accession codes
from their original publications [syn8456629](https://www.synapse.org/#!Synapse:syn8456629) and
[syn3387325](https://www.synapse.org/#!Synapse:syn3387325), respectively, under controlled use
conditions due to human privacy regulations.

Put the downloaded files in `data/`.

## Software environment

The analyses were run in R 4.4.1. The full list of packages and versions is in
[`sessionInfo.txt`](sessionInfo.txt). To regenerate it, run this from the project root:

```sh
Rscript code/sessionInfo.R
```

## Citation

If you use this code, please cite the publication above.

## Session info

Output of `code/sessionInfo.R` (also saved in [`sessionInfo.txt`](sessionInfo.txt)):

```
Generated: 2026-11-20 11:50:09 EST
Scripts scanned: 30

R version 4.4.1 (2024-06-14)
Platform: aarch64-apple-darwin20
Running under: macOS 26.7.1

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] C.UTF-8/C.UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
 [1] splines   parallel  grid      stats4    stats     graphics  grDevices
 [8] utils     datasets  methods   base     

other attached packages:
 [1] VGAM_1.1-13              variancePartition_1.36.3 lubridate_1.9.4         
 [4] tidyverse_2.0.0          tidyr_1.3.1              tidygraph_1.3.1         
 [7] tibble_3.2.1             synapser_2.1.1.259       stringr_1.5.1           
[10] scales_1.3.0             readxl_1.4.3             readr_2.1.5             
[13] RCy3_2.26.0              RColorBrewer_1.1-3       qvalue_2.38.0           
[16] purrr_1.0.2              psych_1.8.12             plyr_1.8.9              
[19] patchwork_1.3.0          PanomiR_1.10.0           org.Hs.eg.db_3.20.0     
[22] AnnotationDbi_1.68.0     Biobase_2.66.0           nanostringr_0.5.0       
[25] knitr_1.49               igraph_2.1.2             Hmisc_5.2-1             
[28] gtsummary_2.1.0          gridExtra_2.3            githubr_0.9.1           
[31] ggsignif_0.6.4           ggraph_2.2.2             forcats_1.0.0           
[34] EnhancedVolcano_1.24.0   ggrepel_0.9.6            ggplot2_3.5.1           
[37] edgeR_4.4.1              limma_3.62.1             dplyr_1.1.4             
[40] doParallel_1.0.17        iterators_1.0.14         foreach_1.5.2           
[43] data.table_1.16.4        cqn_1.52.0               mclust_6.1.1            
[46] cowplot_1.1.3            CovariateAnalysis_0.1.0  ComplexHeatmap_2.22.0   
[49] clusterProfiler_4.14.4   circlize_0.4.16          broom_1.0.7             
[52] Biostrings_2.74.0        GenomeInfoDb_1.42.1      XVector_0.46.0          
[55] IRanges_2.40.1           S4Vectors_0.44.0         BiocGenerics_0.52.0     
[58] biomaRt_2.62.0           BiocParallel_1.40.0     

loaded via a namespace (and not attached):
  [1] fs_1.6.5                matrixStats_1.4.1       bitops_1.0-9           
  [4] enrichplot_1.26.3       httr_1.4.7              numDeriv_2016.8-1.1    
  [7] gh_1.4.1                repr_1.1.7              tools_4.4.1            
 [10] backports_1.5.0         utf8_1.2.4              R6_2.5.1               
 [13] lazyeval_0.2.2          GetoptLong_1.0.5        withr_3.0.2            
 [16] prettyunits_1.2.0       quantreg_5.99.1         cli_3.6.5              
 [19] mvtnorm_1.3-2           yulab.utils_0.1.8       pbdZMQ_0.3-14          
 [22] gson_0.1.0              foreign_0.8-87          DOSE_4.0.0             
 [25] R.utils_2.12.3          rstudioapi_0.17.1       RSQLite_2.3.9          
 [28] generics_0.1.3          gridGraphics_0.5-1      shape_1.4.6.1          
 [31] gtools_3.9.5            GO.db_3.20.0            Matrix_1.7-1           
 [34] fansi_1.0.6             R.methodsS3_1.8.2       lifecycle_1.0.4        
 [37] gplots_3.2.0            BiocFileCache_2.14.0    blob_1.2.4             
 [40] crayon_1.5.3            ggtangle_0.0.5          lattice_0.22-6         
 [43] KEGGREST_1.46.0         pillar_1.9.0            fgsea_1.32.0           
 [46] boot_1.3-31             rjson_0.2.21            corpcor_1.6.10         
 [49] codetools_0.2-20        fastmatch_1.1-4         glue_1.8.0             
 [52] ggfun_0.1.8             Rdpack_2.6.2            vctrs_0.6.5            
 [55] png_0.1-8               treeio_1.30.0           cellranger_1.1.0       
 [58] gtable_0.3.6            cachem_1.1.0            xfun_0.51              
 [61] rbibutils_2.3           reformulas_0.4.0        survival_3.7-0         
 [64] statmod_1.5.0           nlme_3.1-166            pbkrtest_0.5.3         
 [67] ggtree_3.14.0           bit64_4.5.2             EnvStats_3.0.0         
 [70] progress_1.2.3          filelock_1.0.3          rprojroot_2.0.4        
 [73] nor1mix_1.3-3           KernSmooth_2.23-24      rpart_4.1.23           
 [76] colorspace_2.1-1        DBI_1.2.3               nnet_7.3-19            
 [79] mnormt_2.1.1            tidyselect_1.2.1        bit_4.5.0.1            
 [82] compiler_4.4.1          curl_6.0.1              httr2_1.0.7            
 [85] graph_1.84.0            htmlTable_2.4.3         SparseM_1.84-2         
 [88] xml2_1.3.6              checkmate_2.3.2         caTools_1.18.3         
 [91] remaCor_0.0.18          rappdirs_0.3.3          digest_0.6.37          
 [94] minqa_1.2.8             rmarkdown_2.29          aod_1.3.3              
 [97] RhpcBLASctl_0.23-42     htmltools_0.5.8.1       pkgconfig_2.0.3        
[100] base64enc_0.1-3         lme4_1.1-36             dbplyr_2.5.0           
[103] fastmap_1.2.0           rlang_1.1.7             GlobalOptions_0.1.2    
[106] htmlwidgets_1.6.4       UCSC.utils_1.2.0        farver_2.1.2           
[109] jsonlite_1.8.9          GOSemSim_2.32.0         R.oo_1.27.0            
[112] RCurl_1.98-1.16         magrittr_2.0.3          Formula_1.2-5          
[115] GenomeInfoDbData_1.2.13 ggplotify_0.1.2         IRkernel_1.3.2         
[118] munsell_0.5.1           Rcpp_1.0.13-1           ape_5.8                
[121] viridis_0.6.5           reticulate_1.40.0       stringi_1.8.4          
[124] RJSONIO_2.0.0           zlibbioc_1.52.0         MASS_7.3-61            
[127] graphlayouts_1.2.2      IRdisplay_1.1           hms_1.1.3              
[130] locfit_1.5-9.10         uuid_1.2-1              base64url_1.4          
[133] reshape2_1.4.4          XML_3.99-0.17           evaluate_1.0.1         
[136] nloptr_2.2.0            tzdb_0.5.0              tweenr_2.0.3           
[139] MatrixModels_0.5-3      polyclip_1.10-7         clue_0.3-66            
[142] ggforce_0.4.2           fANCOVA_0.6-1           tidytree_0.4.6         
[145] viridisLite_0.4.2       lmerTest_3.1-3          aplot_0.2.3            
[148] memoise_2.0.1           cluster_2.1.7           timechange_0.3.0       
[151] here_1.0.1             

Packages referenced in code/ but not installed / failed to load:
  cardx
  flextable
```
