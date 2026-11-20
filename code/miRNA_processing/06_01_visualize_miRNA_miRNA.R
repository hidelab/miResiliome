rm(list = ls())

# Visualize miRNA–miRNA consensus co-expression network
# - Reads the output from 05_01_miRNA_miRNA_CoExp.R
# - Keeps only edges with p_adj_consensus < 0.05
# - Nodes are miRNAs; edge color encodes sign of r_hat; width encodes |r_hat|
# - Saves a PNG plot and the filtered edge list

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || is.na(args[1])) {
  setup <- ""
  message("No setup provided; defaulting to output/ directory.")
} else {
  setup <- args[1]
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(igraph)
})

options(stringsAsFactors = FALSE)

# Input file produced by 05_01_miRNA_miRNA_CoExp.R
in_primary <- file.path("output", paste0("", setup), "coexpr", "mirna_mirna", "miRNA_miRNA_cor_by_group_resilience.csv")
in_fallback <- file.path("Output", paste0("", setup), "coexpr", "mirna_mirna", "miRNA_miRNA_cor_by_group_resilience.csv")

infile <- if (file.exists(in_primary)) in_primary else in_fallback
if (!file.exists(infile)) stop("Input not found in either path. Looked for: ", in_primary, " and ", in_fallback)

cat("Reading:", infile, "\n")
df <- suppressMessages(readr::read_csv(infile, show_col_types = FALSE))

# Expect columns present from make_consensus(): var1, var2, r_hat, p_adj_consensus
req_cols <- c("var1", "var2", "r_hat", "p_adj_consensus")
missing <- setdiff(req_cols, names(df))
if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))

# Filter to significant consensus edges
fdr_thresh <- 0.05
edges <- df %>%
  filter(is.finite(p_adj_consensus), p_adj_consensus < fdr_thresh) %>%
  transmute(
    from = var1,
    to   = var2,
    r_hat = r_hat,
    p_adj_consensus = p_adj_consensus,
    sign = ifelse(r_hat >= 0, "pos", "neg")
  )

cat("Edges passing FDR < ", fdr_thresh, ": ", nrow(edges), "\n", sep = "")

# Output directory for plots/exports
out_dir <- file.path(dirname(infile))
plot_path <- file.path(out_dir, "miRNA_miRNA_consensus_network.png")
edges_path <- file.path(out_dir, "miRNA_miRNA_consensus_edges_FDR_lt_0.05.csv")

# If nothing passes, save empty file and exit
if (nrow(edges) == 0) {
  warning("No edges passed the FDR threshold; nothing to plot.")
  readr::write_csv(edges, edges_path)
  quit(save = "no")
}

# Build graph (undirected)
g <- graph_from_data_frame(d = edges, directed = FALSE)

# Vertex styling
V(g)$color <- "#737373"  # neutral gray for miRNAs
V(g)$size  <- 4 + 1.5     # slight base size
V(g)$label.cex <- 0.7

# Edge styling by r_hat
rw <- abs(E(g)$r_hat)
mx <- if (length(rw) && is.finite(max(rw, na.rm = TRUE))) max(rw, na.rm = TRUE) else 1
E(g)$width <- 1 + 2 * (rw / ifelse(mx > 0, mx, 1))
rh <- E(g)$r_hat
E(g)$color <- ifelse(rh < 0, "#4575b4", ifelse(rh > 0, "#d73027", "#ffffff"))

# Layout
set.seed(123)
lay <- layout_with_fr(g, niter = 1000)

# Render
png(plot_path, width = 2400, height = 1800, res = 200)
par(mar = c(1, 1, 2, 1))
plot(g,
     layout = lay,
     vertex.label = V(g)$name,
     main = paste0("miRNA–miRNA Consensus Network (FDR < ", fdr_thresh, ")"))
legend("bottomleft",
       title  = "Edge color (r_hat)",
       legend = c("negative (<0)", "zero (=0)", "positive (>0)"),
       col    = c("#4575b4", "#ffffff", "#d73027"),
       lty    = 1, bty = "n")
dev.off()

# Save filtered edges
readr::write_csv(edges, edges_path)

cat("Saved:\n  ", plot_path, "\n  ", edges_path, "\n", sep = "")

# Optional: export to Cytoscape if RCy3 is available and Cytoscape is running
try({
  suppressPackageStartupMessages(library(RCy3))
  if (RCy3::cytoscapePing()) {
    RCy3::createNetworkFromIgraph(g, "miRNA_miRNA_consensus_network")
  }
}, silent = TRUE)

