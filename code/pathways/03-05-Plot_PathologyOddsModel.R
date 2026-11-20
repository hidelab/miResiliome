rm(list = ls())

args = commandArgs(trailingOnly=TRUE)
if(is.na(args[1])){
    setup <- ""
    warning("The experimental design specifications must be set.")
} else {
    setup <- args[1]
}

print("Pathway-pathology correlation analysis")

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
tangle.paths <- read.csv("output/Bulk_paths_prop_odds/tangle_paths.csv")
plaque.paths <- read.csv("output/Bulk_paths_prop_odds/plaque_paths.csv")
all.paths    <- read.csv("output/Bulk_paths_prop_odds/cog_corrected_paths.csv")


cog.paths    <- prepare_prop_odds_table(cog.paths,clean_names)
tangle.paths <- prepare_prop_odds_table(tangle.paths,clean_names) 
plaque.paths <- prepare_prop_odds_table(plaque.paths,clean_names)
all.paths    <- prepare_prop_odds_table(all.paths,clean_names)






#*  plotting odds model for cognition     *# ####


library(gridExtra)
library(cowplot)


sig.mirs   <- cog.paths[1:30,]$gene

tested_gene    <- gene.AD[sig.mirs,]
tested_gene    <- as.data.frame(t(tested_gene))

tested_gene    <- rownames_to_column(tested_gene,var = "SampleID")
tested_gene    <- pivot_longer(tested_gene,cols = sig.mirs)

my_met  <- all.ad[,c("SampleID","cogdx.x")]
my_met$cogdx <- factor(my_met$cogdx.x)
levels(my_met$cogdx) <- c("No impairment", "Mild cognitive impairment", "AD dementia")


tested_gene <- dplyr::left_join(tested_gene,my_met)

tested_gene <- dplyr::left_join(tested_gene,my_met)
more_meta <- cog.paths[,c( "gene","p.adj","Dir","pval" )]

tested_gene <- dplyr::left_join(tested_gene,more_meta,by = c("name" = "gene"))
tested_gene <- dplyr::left_join(tested_gene,clean_names,by = c("name" = "standard_name"))

tested_gene_norm <- tested_gene %>%
  group_by(name) %>%
  mutate(
    baseline = mean(value[cogdx == "No impairment"], na.rm = TRUE),
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
  filter(cogdx == "AD dementia")

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
        nudge_x = 3,                   # push slightly to the right
        direction = "y",
        hjust = 0,
        segment.colour = NA,              # no leader lines
        box.padding  = 0.25
      ) +
        
        geom_hline(yintercept = 0,                 # x‑axis line at zero
                   colour = "black",
                   linewidth = 0.6, linetype = "dashed") +
        
        ## --- NEW: loess smoothing per Pathway -------------
      # geom_smooth(
      #   se       = FALSE,              # hide confidence ribbon
      #   method   = "loess",            # or "gam"
      #   formula  = y ~ as.numeric(x),  # treat the factor as numeric
      #   size     = 0.6,
      #   linetype = "dashed"            # makes the smooth easy to spot
      # ) +
      #   ## -------------------------------------------------

      scale_colour_gradient2(             # red = positive, blue = negative
        low  = "blue",
        mid  = "white",
        high = "firebrick",
        midpoint = 0,
        name = "signed log-p"
      ) +
      theme_classic() +
      xlab("")+
      ylab(paste0("normalized pathway activity")) +
      labs(title = paste0("Cognition"))+
      theme(plot.title = element_text(size=14),
            axis.text=element_text(size=14),
            axis.text.x = element_text(face = "bold",size = 16, angle = 60,
                                       vjust = 1,hjust = 1),
            axis.text.y =  element_text(size = 20),
            axis.title=element_text(size=14,face="bold"),
            strip.text = element_text(size = 14,angle = 90),
            #legend.position = "none",
            axis.line.x = element_blank() 
            )
    





# Doing plaques
 

sig.mirs   <- plaque.paths[1:30,]$gene

tested_gene    <- gene.AD[sig.mirs,]
tested_gene    <- as.data.frame(t(tested_gene))
tested_gene    <- rownames_to_column(tested_gene,var = "SampleID")
tested_gene    <- pivot_longer(tested_gene,cols = sig.mirs)



my_met  <- all.ad[,c("SampleID","ceradsc")]
my_met$ceradsc <- factor(4- my_met$ceradsc,ordered = T)
my_met$cogdx <- factor(my_met$ceradsc)
levels(my_met$cogdx) <- paste0("CERAD ",(as.numeric(levels(my_met$cogdx))))


tested_gene <- dplyr::left_join(tested_gene,my_met)
more_meta <- cog.paths[,c( "gene","p.adj","Dir","pval" )]

tested_gene <- dplyr::left_join(tested_gene,more_meta,by = c("name" = "gene"))
tested_gene <- dplyr::left_join(tested_gene,clean_names,by = c("name" = "standard_name"))

tested_gene_norm <- tested_gene %>%
  group_by(name) %>%
  mutate(
    baseline = mean(value[cogdx == "CERAD 0"], na.rm = TRUE),
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
  filter(cogdx == "CERAD 3")



# replace lab_df$display_name with a shorter version if string length is larger than 80 characters
lab_df$display_name <- ifelse(nchar(lab_df$display_name) > 80,
                              paste0(substr(lab_df$display_name,1,77),"..."),
                              lab_df$display_name)






p2 <- ggplot(summary_df2, 
             aes(x=cogdx, y=mean_expr, 
                 group = name,
                 colour = logp)) +
  geom_point(size = 0) + 
  geom_line() +
  geom_text_repel(                    # add Pathway names
    data = lab_df,
    aes(label = display_name),
    nudge_x = 0.65,                   # push slightly to the right
    direction = "y",
    hjust = 0,
    segment.colour = NA,              # no leader lines
    box.padding  = 0.0
  ) +
  
  geom_hline(yintercept = 0,                 # x‑axis line at zero
             colour = "black",
             linewidth = 0.6, linetype = "dashed") +
  
  ## --- NEW: loess smoothing per Pathway -------------
# geom_smooth(
#   se       = FALSE,              # hide confidence ribbon
#   method   = "loess",            # or "gam"
#   formula  = y ~ as.numeric(x),  # treat the factor as numeric
#   size     = 0.6,
#   linetype = "dashed"            # makes the smooth easy to spot
# ) +
#   ## -------------------------------------------------

scale_colour_gradient2(             # red = positive, blue = negative
  low  = "blue",
  mid  = "white",
  high = "firebrick",
  midpoint = 0,
  name = "signed log-p"
) +
  theme_classic() +
  xlab("")+
  ylab(paste0("normalized pathway activity")) +
  labs(title = paste0("Plaques"))+
  theme(plot.title = element_text(size=14),
        axis.text=element_text(size=14),
        axis.text.x = element_text(face = "bold",size = 16, angle = 60,
                                   vjust = 1,hjust = 1),
        axis.text.y =  element_text(size = 20),
        axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size = 14,angle = 90),
        #legend.position = "none",
        axis.line.x = element_blank() 
  )



# Doing tangles


sig.mirs   <- tangle.paths[1:30,]$gene

tested_gene    <- gene.AD[sig.mirs,]
tested_gene    <- as.data.frame(t(tested_gene))
tested_gene    <- rownames_to_column(tested_gene,var = "SampleID")
tested_gene    <- pivot_longer(tested_gene,cols = sig.mirs)



my_met  <- all.ad[,c("SampleID","braaksc")]
my_met$cogdx <- factor(my_met$braaksc)
levels(my_met$cogdx) <- c("No Tangles",paste0("Braak ",as.roman(1:6)))



tested_gene <- dplyr::left_join(tested_gene,my_met)
more_meta <- cog.paths[,c( "gene","p.adj","Dir","pval" )]

tested_gene <- dplyr::left_join(tested_gene,more_meta,by = c("name" = "gene"))
tested_gene <- dplyr::left_join(tested_gene,clean_names,by = c("name" = "standard_name"))


tested_gene_norm <- tested_gene %>%
  group_by(name) %>%
  mutate(
    baseline = mean(value[cogdx == "No Tangles"], na.rm = TRUE),
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
  filter(cogdx == "Braak VI")


# replace lab_df$display_name with a shorter version if string length is larger than 80 characters
lab_df$display_name <- ifelse(nchar(lab_df$display_name) > 80,
                              paste0(substr(lab_df$display_name,1,77),"..."),
                              lab_df$display_name)





p3 <- ggplot(summary_df2, 
             aes(x=cogdx, y=mean_expr, 
                 group = name,
                 colour = logp)) +
  geom_point(size = 0) + 
  geom_line() +
  geom_text_repel(                    # add Pathway names
    data = lab_df,
    aes(label = display_name),
    nudge_x = 0.65,                   # push slightly to the right
    direction = "y",
    hjust = 0,
    segment.colour = NA,              # no leader lines
    box.padding  = 0.25
  ) +
  
  geom_hline(yintercept = 0,                 # x‑axis line at zero
             colour = "black",
             linewidth = 0.6, linetype = "dashed") +
  
  ## --- NEW: loess smoothing per Pathway -------------
# geom_smooth(
#   se       = FALSE,              # hide confidence ribbon
#   method   = "loess",            # or "gam"
#   formula  = y ~ as.numeric(x),  # treat the factor as numeric
#   size     = 0.6,
#   linetype = "dashed"            # makes the smooth easy to spot
# ) +
#   ## -------------------------------------------------

scale_colour_gradient2(             # red = positive, blue = negative
  low  = "blue",
  mid  = "white",
  high = "firebrick",
  midpoint = 0,
  name = "signed log-p"
) +
  theme_classic() +
  xlab("")+
  ylab(paste0("normalized pathway activity")) +
  labs(title = paste0("Tangles"))+
  theme(plot.title = element_text(size=14),
        axis.text=element_text(size=14),
        axis.text.x = element_text(face = "bold",size = 16, angle = 60,
                                   vjust = 1,hjust = 1),
        axis.text.y =  element_text(size = 20),
        axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size = 14,angle = 90),
        #legend.position = "none",
        axis.line.x = element_blank() 
  )
 


library(patchwork)

p0 = p1 + p2 + p3

ggsave(paste0(output.dir,"/prop_odds_top_pathways.pdf"), width = 35, height = 18)





