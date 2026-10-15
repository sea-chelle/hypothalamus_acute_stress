# load required packages
library(WGCNA)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(igraph)
library(ggraph)
library(ggplot2)

# set working directory
setwd(here("plots"))

# calculate module membership (kME) for grey60 genes
grey60_kME <- cor(vsd_wgcna[, bwnet$colors == "grey60"], bwnet$MEs[, "MEgrey60"], use = "p")

# convert to data frame and order by module membership
grey60_kME <- data.frame(gene = rownames(grey60_kME), kME = grey60_kME[, 1])

grey60_kME$symbol <- mapIds(org.Mm.eg.db,
                            keys = sub("\\..*", "", grey60_kME$gene),
                            column = "SYMBOL",
                            keytype = "ENSEMBL",
                            multiVals = "first")
grey60_kME <- grey60_kME[, c("symbol", "gene", "kME")]
head(grey60_kME, 30)

# calculate adjacency within grey60 module
grey60_adj <- adjacency(vsd_wgcna[, bwnet$colors == "grey60"], power = soft_power, type = "signed")

# calculate intramodular connectivity
grey60_kWithin <- rowSums(grey60_adj) - 1

# add to existing hub gene table
grey60_kME$kWithin <- grey60_kWithin[grey60_kME$gene]

# rank by intramodular connectivity
grey60_kME <- grey60_kME[order(-grey60_kME$kWithin), ]
head(grey60_kME, 20)

# extract unique pairwise adjacency values, excluding self-connections
grey60_edges <- grey60_adj[upper.tri(grey60_adj)]
summary(grey60_edges)
quantile(grey60_edges, c(.90, .95, .975, .99))

# retain the strongest 1% of co-expression connections
grey60_network <- which(grey60_adj >= quantile(grey60_edges, 0.99) & upper.tri(grey60_adj), arr.ind = TRUE)
grey60_network <- data.frame(from = rownames(grey60_adj)[grey60_network[, 1]],
                             to = colnames(grey60_adj)[grey60_network[, 2]],
                             weight = grey60_adj[grey60_network])

# check the number of retained edges and unique genes
nrow(grey60_network)
length(unique(c(grey60_network$from, grey60_network$to)))

# construct an undirected gene co-expression network
grey60_graph <- graph_from_data_frame( grey60_network,directed = FALSE)

# annotate network nodes with gene symbols, module membership and connectivity
V(grey60_graph)$symbol <- grey60_kME$symbol[match(V(grey60_graph)$name, grey60_kME$gene)]
V(grey60_graph)$kME <- grey60_kME$kME[match(V(grey60_graph)$name, grey60_kME$gene)]
V(grey60_graph)$kWithin <- grey60_kME$kWithin[match(V(grey60_graph)$name, grey60_kME$gene)]

# select genes for labelling
label_genes <- c(
  # highly connected genes
  "Nfkbia", "Arl4d", "Tsc22d3", "Slc2a1", "Sgk1", "Cebpd", "Maff", "Ddit4",
  "Gimap6", "Cdkn1a", "Plekhf1", "Trib1", "Pim3", "Csrnp1", "Txnip", "Klf15",
  "Per1", "Fosb",
  
  # informative IEGs
  "Fos", "Junb", "Egr1", "Nr4a1", "Nr4a3", "Fosl2", "Fosb", "Egr4")

V(grey60_graph)$label <- ifelse(V(grey60_graph)$symbol %in% label_genes, V(grey60_graph)$symbol, "")

# plot
set.seed(1234)

ggraph(grey60_graph, layout = "fr") +
  geom_edge_link(aes(width = weight), alpha = 0.6, colour = "grey50") +
  geom_node_point(aes(size = kWithin), colour = "#A995C9", alpha = 0.8) +
  geom_node_text(aes(label = label),
                 repel = TRUE,
                 size = 2.5,
                 fontface = "italic") +
  scale_edge_width(range = c(0.2, 1.2), guide = "none") +
  scale_size_continuous(range = c(1, 3), guide = "none") +
  theme_void(base_size = 8) +
  theme(legend.position = "right")

ggsave("WGCNA_grey60_network.tiff", width = 3.14, height = 4.5, dpi = 600, compression = "lzw")
