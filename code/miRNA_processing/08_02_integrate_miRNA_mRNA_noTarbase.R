rm(list = ls())

# Integrate miRNA–miRNA, miRNA–mRNA, and mRNA–mRNA correlations
# into a single mixed network (two node types: miRNA and mRNA),
# filtering by p_adj_consensus < 0.05, and plot the network.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(igraph)
  library(stringr)
})

options(stringsAsFactors = FALSE)

# Input files (as requested)
in_mirmir   <- file.path("output", "coexpr", "mirna_mirna", "miRNA_miRNA_cor_by_group_resilience.csv")
in_mirmrna  <- file.path("output", "coexpr", "mirna_mrna",  "miRNA_mRNA_cor_by_group_tarbase.csv")
in_mrnamrna <- file.path("output", "coexpr", "mrna_mrna",   "mRNA_mRNA_cor_by_group_resilience.csv")

# Fallback for mRNA–mRNA path (existing script saved into mrna_mirna)
if (!file.exists(in_mrnamrna)) {
  alt_mrnamrna <- file.path("output", "coexpr", "mrna_mirna", "mRNA_mRNA_cor_by_group_resilience.csv")
  if (file.exists(alt_mrnamrna)) in_mrnamrna <- alt_mrnamrna
}

# Validate inputs
stopifnot(file.exists(in_mirmir))
stopifnot(file.exists(in_mirmrna))
stopifnot(file.exists(in_mrnamrna))

cat("Reading:\n",
    "  ", in_mirmir,  "\n",
    "  ", in_mirmrna, "\n",
    "  ", in_mrnamrna,"\n", sep = "")

df_mirmir   <- suppressMessages(readr::read_csv(in_mirmir,  show_col_types = FALSE))
df_mirmrna  <- suppressMessages(readr::read_csv(in_mirmrna, show_col_types = FALSE))
df_mrnamrna <- suppressMessages(readr::read_csv(in_mrnamrna,show_col_types = FALSE))


# Preparing node color info 

mir_cols <- read.csv("output/deMirs/cognition_miRNA_prop_odds.csv")
mir_cols <- mir_cols[mir_cols$p.adj<0.05,c("gene","Dir")]

mrna_cols  <- read.csv("output/DEGenes_bulk/cognition_genes_prop_odds.csv")
mrna_cols  <- mrna_cols[mrna_cols$p.adj<0.05,c("gene","Dir")]


colnames(mir_cols) <- c("name","color")
colnames(mrna_cols) <- c("name","color")

node_cols <- rbind(mir_cols,mrna_cols)
node_cols$Dir <- node_cols$color
node_cols <- node_cols[!duplicated(node_cols$name),]

node_cols$color <- ifelse(node_cols$color == "POS", "red", ifelse(node_cols$color == "NEG", "blue", "gray"))
# Check required columns
req_mm  <- c("var1", "var2", "r_hat", "p_adj_consensus")
req_mrm <- c("mirna_name", "gene_name", "r_hat", "p_adj_consensus")
for (rc in req_mm)  if (!rc %in% names(df_mirmir))   stop("miRNA–miRNA missing column: ", rc)
for (rc in req_mm)  if (!rc %in% names(df_mrnamrna)) stop("mRNA–mRNA missing column: ", rc)
for (rc in req_mrm) if (!rc %in% names(df_mirmrna))  stop("miRNA–mRNA missing column: ", rc)

# Filter by consensus FDR
thresh <- 0.05
r_thresh <- .2
f_mirmir   <- df_mirmir   %>% filter(is.finite(p_adj_consensus), p_adj_consensus < thresh, abs(r_hat)> r_thresh)
f_mirmrna  <- df_mirmrna  %>% filter(is.finite(p_adj_consensus), p_adj_consensus < thresh, abs(r_hat)> r_thresh)
f_mrnamrna <- df_mrnamrna %>% filter(is.finite(p_adj_consensus), p_adj_consensus < thresh, abs(r_hat)> r_thresh)
f_mrnamrna <- f_mrnamrna %>% filter(var1 %in% f_mirmrna$gene_name, var2 %in% f_mirmrna$gene_name) # remove self-loops
cat("Edges after FDR < ", thresh, ":\n",
    "  miRNA–miRNA: ", nrow(f_mirmir),   "\n",
    "  miRNA–mRNA : ", nrow(f_mirmrna),  "\n",
    "  mRNA–mRNA : ",  nrow(f_mrnamrna), "\n", sep = "")

# Build unified edge list
# Ensure undirected edges use a consistent ordering to deduplicate
order_pair <- function(a, b) {
  ifelse(a <= b, paste(a, b, sep = "|"), paste(b, a, sep = "|"))
}

edges_mirmir <- f_mirmir %>%
  transmute(
    from = var1,
    to   = var2,
    edge_type = "miRNA–miRNA",
    r_hat = r_hat,
    sign = ifelse(r_hat >= 0, "pos", "neg"),
    key = order_pair(var1, var2)
  )

edges_mrnamrna <- f_mrnamrna %>%
  transmute(
    from = var1,
    to   = var2,
    edge_type = "mRNA–mRNA",
    r_hat = r_hat,
    sign = ifelse(r_hat >= 0, "pos", "neg"),
    key = order_pair(var1, var2)
  )

edges_mirmrna <- f_mirmrna %>%
  transmute(
    from = mirna_name,
    to   = gene_name,
    edge_type = "miRNA–mRNA",
    r_hat = r_hat,
    sign = ifelse(r_hat >= 0, "pos", "neg"),
    key = order_pair(mirna_name, gene_name)
  )

edges_all <- bind_rows(edges_mirmir, edges_mirmrna, edges_mrnamrna) %>%
  distinct(key, .keep_all = TRUE) %>%
  dplyr::select(-key)

# Derive node table with types
nodes_miRNA <- unique(c(edges_mirmir$from, edges_mirmir$to, edges_mirmrna$from))
nodes_mRNA  <- unique(c(edges_mrnamrna$from, edges_mrnamrna$to, edges_mirmrna$to))

nodes_df <- tibble(name = unique(c(nodes_miRNA, nodes_mRNA))) %>%
  mutate(node_type = case_when(
    name %in% nodes_miRNA ~ "miRNA",
    name %in% nodes_mRNA  ~ "mRNA",
    str_detect(name, "^hsa-", negate = FALSE) ~ "miRNA", # fallback heuristic
    TRUE ~ "mRNA"
  ))

# Quick sanity check for conflicts (should be none)
conflicts <- intersect(nodes_miRNA, nodes_mRNA)
if (length(conflicts) > 0) {
  warning("Nodes appearing as both miRNA and mRNA: ", paste(head(conflicts, 10), collapse = ", "),
          if (length(conflicts) > 10) " ..." else "")
}

# Create output directory
out_dir <- file.path("output", "coexpr", "integrated")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# If no edges, write empty outputs and exit quietly
if (nrow(edges_all) == 0) {
  warning("No edges passed the FDR threshold; nothing to plot.")
  readr::write_csv(edges_all, file.path(out_dir, "integrated_edges_resilience.csv"))
  readr::write_csv(nodes_df,  file.path(out_dir, "integrated_nodes_resilience.csv"))
  quit(save = "no")
}

# Optionally cap edges for plotting if very large
max_plot_edges <- 2000L
if (nrow(edges_all) > max_plot_edges) {
  message("Large graph (", nrow(edges_all), " edges). Keeping top ", max_plot_edges,
          " by |r_hat| for plotting.")
  edges_plot <- edges_all %>% arrange(desc(abs(r_hat))) %>% slice_head(n = max_plot_edges)
} else {
  edges_plot <- edges_all
}

# Build graph
g <- graph_from_data_frame(d = edges_plot, directed = FALSE, vertices = nodes_df)

# Style
# Color nodes by node type
V(g)$color <- ifelse(V(g)$node_type == "miRNA", "#f39c12", "#3498db")
V(g)$color2 <- node_cols$color[match(V(g)$name, node_cols$name)]
V(g)$size  <- 4 + 1.5 # * log1p(degree(g))
V(g)$Dir    <- node_cols$Dir[match(V(g)$name, node_cols$name)]

# Edge aesthetics: width by |r_hat|, color by r_hat (neg=blue, pos=red)
rw <- abs(E(g)$r_hat)
mx <- if (length(rw) && is.finite(max(rw, na.rm = TRUE))) max(rw, na.rm = TRUE) else 1
E(g)$width <- 1 + 2 * (rw / ifelse(mx > 0, mx, 1))

# Simple sign-based edge colors without scaling r_hat
rh <- E(g)$r_hat
E(g)$color <- ifelse(rh < 0, "#4575b4", ifelse(rh > 0, "#d73027", "#ffffff"))

## Force-directed layout and then separate node types along x-axis
set.seed(12)
lay <- layout_with_fr(g, niter = 1000)

# Create horizontal separation between miRNA (left) and mRNA (right)

separate_x_by_type <- FALSE
if (separate_x_by_type &&
    nrow(lay) == vcount(g) &&
    all(c("node_type") %in% colnames(igraph::as_data_frame(g, what = "vertices")))) {
  # Indices for each type
  is_miRNA <- V(g)$node_type == "miRNA"
  is_mRNA  <- V(g)$node_type == "mRNA"

  # Proceed only if both groups exist
  if (any(is_miRNA) && any(is_mRNA)) {
    x <- lay[, 1]
    y <- lay[, 2]

    # Center each group on its own mean (preserve within-group structure)
    x_mi_c <- x[is_miRNA] - mean(x[is_miRNA])
    x_mr_c <- x[is_mRNA]  - mean(x[is_mRNA])

    # Group half-widths along x
    w_mi <- if (length(x_mi_c)) max(abs(x_mi_c)) else 0
    w_mr <- if (length(x_mr_c)) max(abs(x_mr_c)) else 0

    # Determine a margin based on the overall span
    span <- diff(range(x))
    if (!is.finite(span) || span <= 0) span <- 1
    margin <- 0.1 * span # 10% of current span

    # Target centers to ensure a clear gap = margin
    c_mi <- - (w_mi + margin/2)
    c_mr <- + (w_mr + margin/2)

    # Compose new x with separation; keep y as-is
    x_new <- x
    x_new[is_miRNA] <- x_mi_c + c_mi
    x_new[is_mRNA]  <- x_mr_c + c_mr
    lay[, 1] <- x_new
    lay[, 2] <- y
  }
}


png(file.path(out_dir, "integrated_miRNA_mRNA_network.png"), width = 3000, height = 1500, res = 200)
par(mar = c(1,1,1,1))
plot(g,
     layout = lay,
     vertex.label = V(g)$name,
     vertex.color = V(g)$color2,
     vertex.frame.color = NA,
     label.dist = .1,
     main = "Integrated miRNA–mRNA Co-expression Network (p_adj_consensus < 0.05)")
# legend("topleft",
#        legend = c("miRNA", "mRNA"),
#        col    = c("#f39c12", "#3498db"),
#        pch    = c(16, 16),
#        lty    = c(NA, NA),
#        pt.cex = 1.5, bty = "n")

legend("bottomleft",
       title  = "Edge color (r_hat)",
       legend = c("negative (< -0.2)", "zero (=0)", "positive (> 0.2)"),
       col    = c("#4575b4", "#ffffff", "#d73027"),
       lty    = 1, bty = "n")
dev.off()

# Save edge/node tables
readr::write_csv(edges_all, file.path(out_dir, "integrated_edges_resilience.csv"))
readr::write_csv(nodes_df,  file.path(out_dir, "integrated_nodes_resilience.csv"))


cat("\nSaved:\n  ", file.path(out_dir, "integrated_miRNA_mRNA_network.png"), "\n",
    "  ", file.path(out_dir, "integrated_edges_resilience.csv"),   "\n",
    "  ", file.path(out_dir, "integrated_nodes_resilience.csv"),   "\n", sep = "")



    # Plot using ggraph/ggplot2
    suppressPackageStartupMessages({
      library(ggraph)
      library(ggplot2)
      library(tidygraph)
    })

    # Convert to tbl_graph and create a manual layout from the precomputed coords
    tg <- as_tbl_graph(g)
    layout_g <- create_layout(tg, layout = "manual", x = lay[, 1], y = lay[, 2])

     
    set.seed(12)
    p <- ggraph(g, layout = "stress") +geom_edge_link(aes(colour = r_hat)) +
      geom_node_point(aes(color = Dir), size = 5) +
      scale_edge_colour_gradient2(
        low = ("blue"),
        mid = "white",
        high = ("red"),
        midpoint = 0,
        space = "Lab",
        na.value = "grey50",
        guide = "edge_colourbar",
        aesthetics = "edge_colour"
      ) +
      geom_node_text(aes(label = name), vjust = 1.5, size = 5) +
      scale_color_manual(values = c("POS" = "red", "NEG" = "blue", "gray" = "gray")) +
      theme_void() + coord_fixed(ratio =  1)
      ggtitle("Integrated miRNA–mRNA Co-expression Network (p_adj_consensus < 0.05)")
    # Save to the same output directory (match previous PNG dimensions / resolution)
    out_file <- file.path(out_dir, "integrated_miRNA_mRNA_network_ggraph.png")
    ggsave(out_file, p, width = 15, height = 7.5, dpi = 200)


# library(RCy3)

# cytoscapePing()

# createNetworkFromIgraph(g,"miRNA_mRNA_network")
