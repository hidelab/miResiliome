rm(list = ls())

args = commandArgs(trailingOnly=TRUE)
if(is.na(args[1])){
    setup <- ""
    warning("The experimental design specifications must be set.")
} else {
    setup <- args[1]
}

print("Pathway-pathology correlation analysis (cognition update)")

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
library(ggrepel)

doParallel::registerDoParallel(10)

source("code/00-02-functions2.R")



data.dir        <- "data/preprocessed/"
output.dir      <- paste0("output", setup,"/Bulk_paths_prop_odds/")
if(!dir.exists(output.dir)){
  dir.create(output.dir,recursive = T)
}



input.dir       <- paste0(output.dir,"preprocessed/")



# Loading covariates, corrected Pathway, and corrected Pathways
covs   <- readRDS(paste0("data/preprocessed/clinicalTotal.RDS"))
EXPRS  <- readRDS("output/DEPaths_bulk_update/ROSMAP_DLPFC_ADvsControl_Residuals.RDS")


all.ad  <- covs[covs$cogdx.x %in% c(1,2,4),]
gene.AD <- EXPRS[,as.character(all.ad$SampleID)]



clean_names <- read.csv("data/MSigDB_display_names_v2024.csv")
clean_names <- clean_names[,c("standard_name","display_name")]

prepare_prop_odds_table <- function(df,clean_names){
  df <- dplyr::left_join(df,clean_names[,c("standard_name","display_name")],by = c("gene" = "standard_name"))
  df <- df %>%
    dplyr::select(gene,display_name,p.adj,Dir,"Pr...z.." )
  names(df) <- c("gene","display_name","p.adj","Dir","pval")
  return(df)
}

cog.paths    <- read.csv("output/Bulk_paths_prop_odds/cog_paths.csv")
cog.paths    <- prepare_prop_odds_table(cog.paths,clean_names)


pathway_clusters <- read.csv("output/pcxn_out_pfc/Clustering_Cogpaths.csv")
pathway_clusters <- pathway_clusters[pathway_clusters$cluster %in% c(1:6), ]
names(pathway_clusters)
#*  plotting odds model for cognition (updated labels)     *# ####


library(gridExtra)
library(cowplot)


color_limits <- c(-4, 8)

plot_list <- list()
for (i in 1:6)
{

sig.mirs   <- pathway_clusters[pathway_clusters$cluster == i,]$Pathway
sig.mirs   <- cog.paths[cog.paths$gene %in% sig.mirs,]$gene
sig.mirs   <- sig.mirs[1:pmin(30,length(sig.mirs))]



tested_gene    <- gene.AD[sig.mirs,]
tested_gene    <- as.data.frame(t(tested_gene))

tested_gene    <- rownames_to_column(tested_gene,var = "SampleID")
tested_gene    <- pivot_longer(tested_gene,cols = sig.mirs)

my_met  <- all.ad[,c("SampleID","cogdx.x")]
my_met$cogdx <- factor(my_met$cogdx.x)
levels(my_met$cogdx) <- c("No CI", "MCI", "AD")


tested_gene <- dplyr::left_join(tested_gene,my_met)

tested_gene <- dplyr::left_join(tested_gene,my_met)
more_meta <- cog.paths[,c( "gene","p.adj","Dir","pval" )]

tested_gene <- dplyr::left_join(tested_gene,more_meta,by = c("name" = "gene"))
tested_gene <- dplyr::left_join(tested_gene,clean_names,by = c("name" = "standard_name"))

tested_gene_norm <- tested_gene %>%
  group_by(name) %>%
  mutate(
    baseline = mean(value[cogdx == "No CI"], na.rm = TRUE),
    adj_expr = value - baseline,
    logp = ifelse(Dir == "NEG",log10(pval),-log10(pval))
  ) %>%
  ungroup()


summary_df2 <- tested_gene_norm %>%
  group_by(name, cogdx,Dir,logp,display_name) %>%
  summarise(
    mean_expr = mean(adj_expr, na.rm = TRUE),
    .groups = "drop"
  )




lab_df <- summary_df2 %>% 
  filter(cogdx == "AD")

# replace lab_df$display_name with a shorter version if string length is larger than 80 characters
lab_df$display_name <- ifelse(nchar(lab_df$display_name) > 80,
                              paste0(substr(lab_df$display_name,1,77),"..."),
                              lab_df$display_name)



p1 <- ggplot(summary_df2, 
             aes(x=cogdx, y=mean_expr, 
                 group = name,
                 colour = logp)) +
      geom_point(size = 0) + 
      geom_line() +
      geom_text_repel(                    # add Pathway names
        data = lab_df,
        aes(label = display_name),
        nudge_x = 0.1,                # slight push right of the line end
        direction = "y",
        hjust = 0,                    # left-align: text starts near line end
        segment.colour = "black",
        segment.linetype = "dashed",
        box.padding  = 0.25
      ) +
        geom_hline(yintercept = 0,                 # x‑axis line at zero
                   colour = "black",
                   linewidth = 0.6, linetype = "dashed") +
      scale_x_discrete(expand = expansion(mult = c(0.05, 4))) +
      scale_colour_gradient2(             # red = positive, blue = negative
        low  = "blue",
        mid  = "white",
        high = "firebrick",
        midpoint = 0,
        limits = c(-8,9),
        name = "signed log-p"
      ) +
      coord_cartesian(clip = "off") +
      theme_classic() +
      xlab("")+ 
      ylab(paste0("normalized pathway activity")) +
      labs(title = paste0("Cognition cluster ",i))+
      theme(plot.title = element_text(size=14),
            axis.text=element_text(size=14),
            axis.text.x = element_text(face = "bold",size = 16, angle = 60,
                                       vjust = 1,hjust = 1),
            axis.text.y =  element_text(size = 20),
            axis.title=element_text(size=14,face="bold"),
            strip.text = element_text(size = 14,angle = 90),
            #legend.position = "none",
            axis.line.x = element_blank(),
            plot.margin = margin(5.5, 130, 5.5, 5.5) 
            )
    

plot_list[[i]] <- p1


ggsave(paste0(output.dir,"/cluster_",i," prop_odds_cognitive_pathways.pdf"),p1, width = 14, height = 10)


}


library(patchwork)

p0 = plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + plot_list[[4]] + plot_list[[5]] + plot_list[[6]] + 
  plot_layout(ncol = 2, nrow = 3, guides = "collect") & theme(legend.position = 'bottom')

ggsave(paste0(output.dir,"/prop_odds_cognitive_pathways.pdf"),p0, width = 22, height = 34)

