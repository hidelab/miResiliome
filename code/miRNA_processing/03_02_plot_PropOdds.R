rm(list = ls())

args = commandArgs(trailingOnly=TRUE)
if(is.na(args[1])){
    setup <- ""
    warning("The experimental design specifications must be set.")
} else {
    setup <- args[1]
}

print("Plotting proportional odds models")

options(StringAsFactors = F)

library(dplyr)
library(tidyverse)
library(ggrepel)
library(gridExtra)
library(cowplot)
library(patchwork)

output.dir      <- paste0("output", setup,"/deMirs/")

# Loading covariates and corrected miRNA
covs   <- readRDS(file = "data/preprocessed/mirClinical.RDS")
EXPRS  <-  readRDS("data/preprocessed/ResidualmiRNA.RDS")

all.ad  <- covs[covs$cogdx %in% c(1,2,4),]
gene.AD <- EXPRS[,as.character(all.ad$mirna_id)]

# Loading proportional odds results
cog.paths    <- read.csv(paste0(output.dir,"/cognition_miRNA_prop_odds.csv"))
plaque.paths <- read.csv(paste0(output.dir,"/plaque_miRNA_prop_odds.csv"))
tangle.paths <- read.csv(paste0(output.dir,"/tangle_miRNA_prop_odds.csv"))

# Compute shared color scale across all three plots
all_logp <- c(
  ifelse(cog.paths$Dir == "NEG", log10(cog.paths$pval), -log10(cog.paths$pval)),
  ifelse(plaque.paths$Dir == "NEG", log10(plaque.paths$pval), -log10(plaque.paths$pval)),
  ifelse(tangle.paths$Dir == "NEG", log10(tangle.paths$pval), -log10(tangle.paths$pval))
)
logp_range <- range(all_logp, na.rm = TRUE)

#*  plotting odds model for cognition     *# ####


sig.mirs   <- c("hsa-miR-129-5p", "hsa-miR-132-3p", "hsa-miR-29a-3p",
                "hsa-miR-362-3p")

sig.mirs   <- cog.paths[cog.paths$p.adj < 0.05,]$gene

tested_gene    <- gene.AD[sig.mirs,]
tested_gene    <- as.data.frame(t(tested_gene))

tested_gene    <- rownames_to_column(tested_gene,var = "mirna_id")
tested_gene    <- pivot_longer(tested_gene,cols = sig.mirs)

my_met  <- all.ad[,c("mirna_id","cogdx")]
my_met$cogdx <- factor(my_met$cogdx)
levels(my_met$cogdx) <- c("No impairment", "Mild cognitive impairment", "AD dementia")


tested_gene <- dplyr::left_join(tested_gene,my_met)

names(cog.paths)
tested_gene <- dplyr::left_join(tested_gene,my_met)
more_meta <- cog.paths[,c( "gene","R2lat","p.adj","Dir","pval","qval" )]

tested_gene <- dplyr::left_join(tested_gene,more_meta,by = c("name" = "gene"))


tested_gene_norm <- tested_gene %>%
  group_by(name) %>%
  mutate(
    baseline = mean(value[cogdx == "No impairment"], na.rm = TRUE),
    adj_expr = value - baseline,
    logp = ifelse(Dir == "NEG",log10(pval),-log10(pval))
  ) %>%
  ungroup()


summary_df2 <- tested_gene_norm %>%
  group_by(name, cogdx,Dir,logp) %>%
  summarise(
    mean_expr = mean(adj_expr, na.rm = TRUE),
    .groups = "drop"
  )


lab_df <- summary_df2 %>%
  filter(cogdx == "AD dementia")


p1 <- ggplot(summary_df2,
             aes(x=cogdx, y=mean_expr,
                 group = name,
                 colour = logp)) +
      geom_point(size = 0) +
      geom_line() +
      geom_text_repel(                    # add miRNA names
        data = lab_df,
        aes(label = name),
        nudge_x = 0.65,                   # push slightly to the right
        direction = "y",
        hjust = 0,
        segment.colour = NA,              # no leader lines
        box.padding  = 0.25
      ) +

        geom_hline(yintercept = 0,                 # x‑axis line at zero
                   colour = "black",
                   linewidth = 0.6, linetype = "dashed") +

        ## --- NEW: loess smoothing per miRNA -------------
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
        mid  = "grey70",
        high = "firebrick",
        midpoint = 0,
        limits = logp_range,
        name = "signed log-p"
      ) +
      theme_classic() +
      xlab("")+
      ylab(paste0("normalized miRNA expression")) +
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


sig.mirs   <- plaque.paths[plaque.paths$p.adj < 0.05,]$gene

tested_gene    <- gene.AD[sig.mirs,]
tested_gene    <- as.data.frame(t(tested_gene))

tested_gene    <- rownames_to_column(tested_gene,var = "mirna_id")
tested_gene    <- pivot_longer(tested_gene,cols = sig.mirs)

my_met  <- all.ad[,c("mirna_id","ceradsc")]
my_met$ceradsc <- factor(4- my_met$ceradsc,ordered = T)
my_met$cogdx <- factor(my_met$ceradsc)
levels(my_met$cogdx) <- paste0("CERAD ",(as.numeric(levels(my_met$cogdx))))


tested_gene <- dplyr::left_join(tested_gene,my_met)

names(cog.paths)
tested_gene <- dplyr::left_join(tested_gene,my_met)
more_meta <- plaque.paths[,c( "gene","R2lat","p.adj","Dir","pval","qval" )]

tested_gene <- dplyr::left_join(tested_gene,more_meta,by = c("name" = "gene"))


tested_gene_norm <- tested_gene %>%
  group_by(name) %>%
  mutate(
    baseline = mean(value[cogdx == "CERAD 0"], na.rm = TRUE),
    adj_expr = value - baseline,
    logp = ifelse(Dir == "NEG",log10(pval),-log10(pval))
  ) %>%
  ungroup()


summary_df2 <- tested_gene_norm %>%
  group_by(name, cogdx,Dir,logp) %>%
  summarise(
    mean_expr = mean(adj_expr, na.rm = TRUE),
    .groups = "drop"
  )


lab_df <- summary_df2 %>%
  filter(cogdx == "CERAD 3")


p2 <- ggplot(summary_df2,
             aes(x=cogdx, y=mean_expr,
                 group = name,
                 colour = logp)) +
    geom_point(size = 0) +
    geom_line() +
    geom_text_repel(                    # add miRNA names
      data = lab_df,
      aes(label = name),
      nudge_x = 0.65,                   # push slightly to the right
      direction = "y",
      hjust = 0,
      segment.colour = NA,              # no leader lines
      box.padding  = 0.25
    ) +

    geom_hline(yintercept = 0,                 # x‑axis line at zero
              colour = "black",
              linewidth = 0.6, linetype = "dashed") +
    scale_colour_gradient2(             # red = positive, blue = negative
      low  = "blue",
      mid  = "grey70",
      high = "firebrick",
      midpoint = 0,
      limits = logp_range,
      name = "signed log-p"
    ) +
      theme_classic() +
      xlab("")+
      ylab(paste0("normalized miRNA expression")) +
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


sig.mirs   <- tangle.paths[tangle.paths$p.adj < 0.05,]$gene

tested_gene    <- gene.AD[sig.mirs,]
tested_gene    <- as.data.frame(t(tested_gene))

tested_gene    <- rownames_to_column(tested_gene,var = "mirna_id")
tested_gene    <- pivot_longer(tested_gene,cols = sig.mirs)

my_met  <- all.ad[,c("mirna_id","braaksc")]
my_met$cogdx <- factor(my_met$braaksc)
levels(my_met$cogdx) <- c("No Tangles",paste0("Braak ",as.roman(1:6)))


tested_gene <- dplyr::left_join(tested_gene,my_met)

names(cog.paths)
tested_gene <- dplyr::left_join(tested_gene,my_met)
more_meta <- tangle.pathsa[,c( "gene","R2lat","p.adj","Dir","pval","qval" )]

tested_gene <- dplyr::left_join(tested_gene,more_meta,by = c("name" = "gene"))


tested_gene_norm <- tested_gene %>%
  group_by(name) %>%
  mutate(
    baseline = mean(value[cogdx == "No Tangles"], na.rm = TRUE),
    adj_expr = value - baseline,
    logp = ifelse(Dir == "NEG",log10(pval),-log10(pval))
  ) %>%
  ungroup()


summary_df2 <- tested_gene_norm %>%
  group_by(name, cogdx,Dir,logp) %>%
  summarise(
    mean_expr = mean(adj_expr, na.rm = TRUE),
    .groups = "drop"
  )


lab_df <- summary_df2 %>%
  filter(cogdx == "Braak VI")


p3 <- ggplot(summary_df2,
             aes(x=cogdx, y=mean_expr,
                 group = name,
                 colour = logp)) +
  geom_point(size = 0) +
  geom_line() +
  geom_text_repel(                    # add miRNA names
    data = lab_df,
    aes(label = name),
    nudge_x = 0.65,                   # push slightly to the right
    direction = "y",
    hjust = 0,
    segment.colour = NA,              # no leader lines
    box.padding  = 0.25
  ) +

  geom_hline(yintercept = 0,                 # x‑axis line at zero
             colour = "black",
             linewidth = 0.6, linetype = "dashed") +



scale_colour_gradient2(             # red = positive, blue = negative
  low  = "blue",
  mid  = "grey70",
  high = "firebrick",
  midpoint = 0,
  limits = logp_range,
  name = "signed log-p"
) +
  theme_classic() +
  xlab("")+
  ylab(paste0("normalized miRNA expression")) +
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



p0 = p1 + p2 + p3

ggsave(paste0(output.dir,"/prop_odds_top_miRNAs.pdf"), width = 22, height = 11)
ggsave(paste0("output/figures_paper","/Figure 1C prop_odds_top_miRNAs.jpg"), width = 22, height = 11,
dpi = 600)

