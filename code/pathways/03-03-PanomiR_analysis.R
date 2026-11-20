
rm(list = ls())
library(PanomiR)
library(data.table)
library(dplyr)
library(tibble)




my_paths <- read.csv("output/Bulk_paths_prop_odds/cog_paths.csv", row.names = 1)
rownames(my_paths) <- my_paths$gene
my_paths$adj.P.Val <- my_paths$p.adj
my_paths$AveExpr <- 0
my_paths$logFC <- ifelse(my_paths$Dir == "POS", 1, -1)
my_genes <- read.csv("output/DEGenes_bulk/ADvsResilient.csv",row.names = 1)
my_mirs  <- read.csv("output/deMirs/deMirsADvsResilient.csv", row.names = 1)
#PCXN <- read.csv("data/pcxn_MSigDB2024_pathway_estimates.tsv",sep = "\t")
PCXN <- read.csv("data/pathway_estimates_msigdb2024_refactor.tsv",sep = "\t")


# PCXN <- read.csv("data/pcxn_MSigDB2024_pathway_estimates.tsv",sep = "\t")
# # names(PCXN)


#  PCXN <- read.csv("data/Brain_BrainFrontalCortex-BA9_cpad_pathcor.tsv",sep = "\t")
#   PCXN <- PCXN[,c("Pathway.A","Pathway.B","estimate","p.value")]
 
#   PCXN$p.adj <- p.adjust(PCXN$p.value,method = "BH")
#   names(PCXN) <- c("Pathway.A","Pathway.B","PathCor","p.value","p.adj")

#PCXN0 <- PCXN



out.dir <- "output/pcxn_out_pfc/"
dir.create(out.dir,showWarnings = FALSE,recursive = TRUE)
set.seed(2)
pathwayClusts <- mappingPathwaysClusters(
  pcxn = PCXN,
  dePathways = my_paths[1:500,],
  topClusters = 6,
  topPathways = 500,
  outDir= paste0(out.dir),
  plotOut = TRUE,
  pathwayFDR = 0.01,
  subplot = FALSE,
  prefix='',
  edgeFDR = 0.05,
  clusteringFunction = "cluster_louvain",
  correlationCutOff = sqrt(.2))

temp <- as.data.frame(igraph::as_long_data_frame(pathwayClusts$DE_PCXN))

temp <- pathwayClusts$Clustering

write.csv(temp, file = paste0(out.dir,"/Clustering_Cogpaths.csv"))






new_cluster_plot <- function (subNet, subplot = FALSE, topClusters = 6, prefix = "", 
          outDir = ".", plotSave = TRUE) {
    if (!dir.exists(outDir)) {
        stop("Output directory does not exist.")
    }
    figDir <- paste0(outDir, "/", prefix)
    legendCats <- data.frame(attr = c("Up-regulated", "Down-regulated"), 
                             shape = unique(igraph::V(subNet)$shape))
    cols <- RColorBrewer::brewer.pal(8, "Pastel1")
    clstMems <- igraph::V(subNet)$cluster
    nodeColors <- cols[clstMems]
    smallClust <- which(table(clstMems) <= table(clstMems)[5], 
                        useNames = TRUE)
    nodeColors[clstMems %in% smallClust] <- NA
    clusterPlotHelper(plotSave, figDir, subNet, nodeColors, 
                       legendCats)
    clusterSubPlotHelper(subplot, topClusters, clstMems, subNet, 
                          figDir, cols, legendCats)
    "%notin%" <- Negate("%in%")
    remove <- which(table(clstMems) < 5)
    remove <- which((clstMems %notin% remove))
    subNet2 <- igraph::induced_subgraph(subNet, remove)
    if (subplot == TRUE) {
        grDevices::pdf(paste0(figDir, "ConnectedPathways_PCxNCorGraph.pdf"), 
                       width = 18, height = 11)
        plot(subNet2, edge.width = 1.3, vertex.size = 5, vertex.label = NA, 
             vertex.color = cols[clstMems[remove]], legend = TRUE, 
             layout = igraph::layout_components)
        graphics::legend(x = "bottomright", y = 200, legend = legendCats$attr, 
                         pch = c(0, 1), bty = "n", cex = 1.4)
        graphics::legend(x = "topleft", legend = c("Positive Cor", 
                                                   "Negative Cor"), col = c("#E41A1C", "#377EB8"), lty = 1, 
                         lwd = 2, cex = 1.4, bty = "n")
        grDevices::dev.off()
        saveRDS(subNet2,paste0(figDir,"network_plot_data.RDS"))
    }
}


clusterPlotHelper <- function(plotSave, figDir, subNet, nodeColors,
                               legendCats) {
    if (plotSave) {
        grDevices::pdf(paste0(figDir, "PCxNCorGraph.pdf"),
                       width = 18, height = 11)
    }
    plot(subNet, vertex.size = 5, vertex.label = NA, vertex.color = nodeColors)
    graphics::legend(x = "bottomleft", legend = legendCats$attr, pch = c(0, 1),
                     bty = "n", cex = 1.6)
    graphics::legend(x = "topleft", legend = c("Positive Cor", "Negative Cor"),
                     col = c("#E41A1C", "#377EB8"), lty = 1, lwd = 2,
                     cex = 1.6, bty = "n")
    if (plotSave)
        grDevices::dev.off()
}

clusterSubPlotHelper <- function(subplot, topClusters, clstMems, subNet,
                                  figDir, cols, legendCats) {
    if (subplot == TRUE) {
        for (k in seq_len(topClusters)) {
            keep <- which((clstMems) == k)
            subNet2 <- igraph::induced_subgraph(subNet, keep)
            if (length(igraph::V(subNet2)) < 2) next
            
            grDevices::pdf(paste0(figDir, "PCxNCorGraph_",
                                  "Cluster_", k, ".pdf"))
            plot(subNet2, edge.width = 1.3, vertex.size = 5, vertex.label = NA,
                 vertex.color = cols[clstMems[keep]], legend = TRUE,
                 layout = igraph::layout.fruchterman.reingold)
            graphics::legend(x = "bottomleft", legend = legendCats$attr,
                             pch = c(0, 1), bty = "n", cex = 1.4)
            graphics::legend(x = "topleft",
                             legend = c("Positive Cor", "Negative Cor"),
                             col = c("#E41A1C", "#377EB8"),
                             lty = 1, lwd = 2, cex = 1.4, bty = "n")
            grDevices::dev.off()
        }
    }
}



set.seed(3)
new_cluster_plot(pathwayClusts$DE_PCXN,
            subplot = T, plotSave = T, outDir = out.dir, prefix = "costumized")






table(pathwayClusts$Clustering$cluster)

# instead of the following line, run this line. 
#data("msigdb_c2")


msigdb_c2 <- readRDS("data/MSigDB2024.RDS")
length(msigdb_c2)


genestargets <- fread("data/Homo_sapiens_TarBase-v9.tsv.gz", header = TRUE)

# Preparing  Tarbase
Tar.dir <- "data/Homo_sapiens_TarBase-v9.tsv.gz"


tarbase.v9.prepare <- function(Tar.dir = NULL,
                               species = "Homo sapiens",
                               microt_thresh = 0.5,
                               gene_type = c("ensembl", "symbol"),
                               min.size.thresh = 10,
                               max.size.thresh = 5000){
  
  if(length(gene_type) > 1) {
    message("more than one gene name selected")
    test_me <- match.arg( c("ensembl", "symbol"),
                          gene_type,
                          several.ok = TRUE)
    gene_type <- test_me[[1]]
  }
  
  if(gene_type =="ensembl"){
    gene_var = c("gene_id")
  }else if(gene_type =="symbol"){
    gene_var = c("gene_name")
  }
  tarbase   <- fread(Tar.dir,header = T)
  mir.sets  <- tarbase %>% dplyr::filter(.,species == "Homo sapiens",
                                         regulation == "Negative",
                                         confidence == 1,
                                         microt_score > microt_thresh,
                                         !is.na(microt_score)) 
  
  mir.sets  <- mir.sets %>% dplyr::select(.,mirna_name,!!sym(gene_var)) %>% group_split(.,mirna_name) 
  mir.list  <- sapply(mir.sets, function(X){dplyr::slice(X,1) %>% 
      dplyr::select(.,mirna_name) %>% pull})
  mir.sets  <- sapply(mir.sets,function(X){dplyr::select(X,!!sym(gene_var)) %>% unique})
  
  names(mir.sets) <- mir.list
  mir.sets <- mir.sets[unlist(lapply(mir.sets, length)) >= min.size.thresh]
  mir.sets <- mir.sets[unlist(lapply(mir.sets, length)) <= max.size.thresh]
  
  
  return(mir.sets)
}




miRNA_targets <- tarbase.v9.prepare(Tar.dir = ,Tar.dir, gene_type = "ensembl")



library(org.Hs.eg.db)

miRNA_targets2 <- lapply(
  miRNA_targets,
  function(x) {
    res <- clusterProfiler::bitr(x, fromType = "ENSEMBL", toType = "ENTREZID",
                                 OrgDb = org.Hs.eg.db::org.Hs.eg.db)
    return(unique(res$ENTREZID))
  }
)


#Work it out from here

customeTableEnrich <- miRNAPathwayEnrichment(mirSets = miRNA_targets2,
                                                  pathwaySets = msigdb_c2,
                                                  geneSelection = my_genes$ensembl_gene_id,
                                                  mirSelection = rownames(my_mirs),
                                                  fromID = "ENSEMBL",
                                                  toID = "ENTREZID",
                                                  minPathSize = 10,
                                                  numCores = 8,
                                                  outDir = ".",
                                                  saveOutName = NULL)

write.csv(customeTableEnrich, file = "output/Enrich_msigdb2024_tarbase.csv")
# customeTableEnrich <- read.csv("output/Enrich_msigdb2024_tarbase.csv", row.names = 1)



set.seed(1)
output2_tgct <- prioritizeMicroRNA(enriches0 = customeTableEnrich,
                                   pathClust = pathwayClusts$Clustering,
                                   topClust = 6,
                                   enrichmentFDR = 0.05,
                                   sampRate = 2000,
                                   method = c("aggInv"),
                                   outDir = out.dir,
                                   dataDir = out.dir,
                                   saveSampling = TRUE,
                                   runJackKnife = FALSE,
                                   numCores = 4,
                                   prefix = "ROSMAP_cog_aggInv",
                                   saveCSV = TRUE)

