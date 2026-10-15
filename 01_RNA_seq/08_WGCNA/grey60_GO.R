# load required packages
library(here)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)
library(grid)

# set working directory
setwd(here("plots"))

# annotate background genes included in wgcna cluster
wgcna_annot <- bitr(colnames(vsd_wgcna),
                    fromType = "ENSEMBL",
                    toType = c("SYMBOL", "ENTREZID"),
                    OrgDb = org.Mm.eg.db)

# convert grey 60 genes to SYMBOL and ENTREZ ID
grey60_genes <- bitr(grey60_genes, 
                    fromType = "ENSEMBL",
                    toType = c("SYMBOL", "ENTREZID"),
                    OrgDb = org.Mm.eg.db)

# perform GO on grey60 module genes vs all WGCNA genes
grey60_GO <- clusterProfiler::enrichGO(gene = grey60_genes$ENTREZID,
                                       universe = wgcna_annot$ENTREZID,
                                       OrgDb = org.Mm.eg.db,
                                       keyType = "ENTREZID",
                                       ont = "BP",
                                       pAdjustMethod = "BH",
                                       pvalueCutoff = 0.05,
                                       qvalueCutoff = 0.05,
                                       readable = TRUE)

# look at significant terms
grey60_GO@result %>%
  filter(p.adjust <= 0.05, Count >=2) %>%
  arrange(p.adjust) %>%
  pull(Description)

# perform KEGG on grey60 module genes vs all WGCNA genes
grey60_KEGG <- enrichKEGG(gene = grey60_genes$ENTREZID,
                          universe = wgcna_annot$ENTREZID,
                          organism = "mmu",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)

# look at significant terms
grey60_KEGG@result %>%
  filter(p.adjust <= 0.05, Count >=2) %>%
  arrange(p.adjust) %>%
  pull(Description)

# select terms
grey60_GO_select <- c("p38MAPK cascade", "response to corticosteroid", "response to progesterone", 
                      "response to reactive oxygen species", "regulation of inflammatory response")

# filter for selected terms, calculate GeneRatio, FDR and shortern plot labels
grey60_GO_plot <- grey60_GO@result %>%
  filter(Description %in% grey60_GO_select) %>%
  mutate(GeneRatio = sapply(GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  }),
  neglog10FDR = -log10(p.adjust),
  PlotLabel = case_when(Description == "response to corticosteroid" ~ "Corticosteroid\nresponse",
                        Description == "response to reactive oxygen species" ~ "ROS response",
                        Description == "p38MAPK cascade" ~ "p38 MAPK\ncascade",
                        Description == "response to progesterone" ~ "Progesterone\nresponse",
                        Description == "regulation of inflammatory response" ~ "Inflammatory\nresponse",
                        TRUE ~ Description))

# dot plot
ggplot(grey60_GO_plot, aes(x = GeneRatio, y = reorder(PlotLabel, GeneRatio), size = Count,  colour = neglog10FDR)) +
  geom_point() +
  scale_colour_gradient(low = "#E9DFFC", high = "#7656A8") +
  scale_size_continuous(range = c(1.5, 5), breaks = c(4, 6, 8, 10)) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.12))) +
  guides(size = guide_legend(override.aes = list(size = c(1, 2, 3, 4)))) +
  labs(x = "GeneRatio",  y = NULL,size = "Count", colour = expression(-log[10]("FDR"))) +
  theme_bw() +
  theme(panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold"),
        axis.text = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(3, "mm"),
        legend.key.height = unit(3, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "pt"))

ggsave("grey60_GO_dotplot.tiff", height = 2.2, width = 3.2, dpi = 600, compression = "lzw")
