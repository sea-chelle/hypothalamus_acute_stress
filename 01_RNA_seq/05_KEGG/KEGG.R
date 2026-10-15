# load required packages
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(dplyr)

# define background genes with valid adjusted p-values
background_basal_vs_10 <- rownames(res_basal_vs_10)[!is.na(res_basal_vs_10$padj)]

# convert Ensembl IDs to Entrez IDs
entrez_background_basal_vs_10 <- mapIds(org.Mm.eg.db,
                                        keys = background_basal_vs_10,
                                        column = "ENTREZID",
                                        keytype = "ENSEMBL",
                                        multiVals = "first")

# remove missing IDs and duplicates
entrez_background_basal_vs_10 <- unique(na.omit(unname(entrez_background_basal_vs_10)))

# convert to Ensembl IDs to entrez IDs
entrez_basal_vs_10 <- bitr(rownames(sig_DEG_basal_vs_10),
                           fromType = "ENSEMBL",
                           toType = "ENTREZID",
                           OrgDb = "org.Mm.eg.db") %>%
  pull(ENTREZID)

# run KEGG pathway analysis on all DEGs
KEGG_basal_vs_10 <- enrichKEGG(gene = entrez_basal_vs_10,
                               universe = entrez_background_basal_vs_10,
                               organism = "mmu", 
                               pvalueCutoff = 0.05)

# create column for significant terms only (p.adj < 0.05 and more than 1 gene driving term)
KEGG_basal_vs_10@result <- KEGG_basal_vs_10@result %>%
  mutate(significant = ifelse(p.adjust <= 0.05 & Count >= 2, "Significant", "Not significant"))

# define background genes with valid adjusted p-values
background_basal_vs_60 <- rownames(res_basal_vs_60)[!is.na(res_basal_vs_60$padj)]

# convert Ensembl IDs to Entrez IDs
entrez_background_basal_vs_60 <- mapIds(org.Mm.eg.db,
                                        keys = background_basal_vs_60,
                                        column = "ENTREZID",
                                        keytype = "ENSEMBL",
                                        multiVals = "first")

# remove missing IDs and duplicates
entrez_background_basal_vs_60 <- unique(na.omit(unname(entrez_background_basal_vs_60)))

# convert to Ensembl IDs to entrez IDs
entrez_upregulated_basal_vs_60 <- bitr(upregulated_basal_vs_60,
                                       fromType = "ENSEMBL",
                                       toType = "ENTREZID",
                                       OrgDb = "org.Mm.eg.db") %>%
  pull(ENTREZID)

entrez_downregulated_basal_vs_60 <- bitr(downregulated_basal_vs_60,
                                         fromType = "ENSEMBL",
                                         toType = "ENTREZID",
                                         OrgDb = "org.Mm.eg.db") %>%
  pull(ENTREZID)


# run KEGG pathway analysis on upregulated DEGs
KEGG_upregulated_basal_vs_60 <- enrichKEGG(gene = entrez_upregulated_basal_vs_60,
                                           universe = entrez_background_basal_vs_60,
                                           organism = "mmu", 
                                           pvalueCutoff = 0.05)

# create column for significant terms only (p.adj < 0.05 and more than 1 gene driving term)
KEGG_upregulated_basal_vs_60@result <- KEGG_upregulated_basal_vs_60@result %>%
  mutate(significant = ifelse(p.adjust <= 0.05 & Count >= 2, "Significant", "Not significant"))

# run KEGG pathway analysis on downregulated DEGs
KEGG_downregulated_basal_vs_60 <- enrichKEGG(gene = entrez_downregulated_basal_vs_60,
                                             universe = entrez_background_basal_vs_60,
                                             organism = "mmu", 
                                             pvalueCutoff = 0.05)

# create column for significant terms only (p.adj < 0.05 and more than 1 gene driving term)
KEGG_downregulated_basal_vs_60@result <- KEGG_downregulated_basal_vs_60@result %>%
  mutate(significant = ifelse(p.adjust <= 0.05 & Count >= 2, "Significant", "Not significant"))

# define background genes with valid adjusted p-values
background_10_vs_60 <- rownames(res_10_vs_60)[!is.na(res_10_vs_60$padj)]

# convert Ensembl IDs to Entrez IDs
entrez_background_10_vs_60 <- mapIds(org.Mm.eg.db,
                                     keys = background_10_vs_60,
                                     column = "ENTREZID",
                                     keytype = "ENSEMBL",
                                     multiVals = "first")

# remove missing IDs and duplicates
entrez_background_10_vs_60 <- unique(na.omit(unname(entrez_background_10_vs_60)))

# convert to Ensembl IDs to entrez IDs 
entrez_upregulated_10_vs_60 <- bitr(upregulated_10_vs_60,
                                    fromType = "ENSEMBL",
                                    toType = "ENTREZID",
                                    OrgDb = "org.Mm.eg.db") %>%
  pull(ENTREZID)

entrez_downregulated_10_vs_60 <- bitr(downregulated_10_vs_60,
                                      fromType = "ENSEMBL",
                                      toType = "ENTREZID",
                                      OrgDb = "org.Mm.eg.db") %>%
  pull(ENTREZID)


# run KEGG pathway analysis on upregulated DEGs
KEGG_upregulated_10_vs_60 <- enrichKEGG(gene = entrez_upregulated_10_vs_60,
                                        universe = entrez_background_10_vs_60,
                                        organism = "mmu", 
                                        pvalueCutoff = 0.05)

# create column for significant terms only (p.adj < 0.05 and more than 1 gene driving term)
KEGG_upregulated_10_vs_60@result <- KEGG_upregulated_10_vs_60@result %>%
  mutate(significant = ifelse(p.adjust <= 0.05 & Count >= 2, "Significant", "Not significant"))

# run KEGG pathway analysis on downregulated DEGs
KEGG_downregulated_10_vs_60 <- enrichKEGG(gene = entrez_downregulated_10_vs_60,
                                          universe = entrez_background_10_vs_60,
                                          organism = "mmu", 
                                          pvalueCutoff = 0.05)

# create column for significant terms only (p.adj < 0.05 and more than 1 gene driving term)
KEGG_downregulated_10_vs_60@result <- KEGG_downregulated_10_vs_60@result %>%
  mutate(significant = ifelse(p.adjust <= 0.05 & Count >= 2, "Significant", "Not significant"))

# define background genes with valid adjusted p-values
background_male_basal_vs_10 <- rownames(res_male_basal_vs_10)[!is.na(res_male_basal_vs_10$padj)]

# convert Ensembl IDs to Entrez IDs
entrez_background_male_basal_vs_10 <- mapIds(org.Mm.eg.db,
                                        keys = background_male_basal_vs_10,
                                        column = "ENTREZID",
                                        keytype = "ENSEMBL",
                                        multiVals = "first")

# remove missing IDs and duplicates
entrez_background_male_basal_vs_10 <- unique(na.omit(unname(entrez_background_male_basal_vs_10)))

# convert to Ensembl IDs to entrez IDs 
entrez_male_basal_vs_10 <- bitr(rownames(sig_DEG_male_basal_vs_10),
                                fromType = "ENSEMBL",
                                toType = "ENTREZID",
                                OrgDb = "org.Mm.eg.db") %>%
  pull(ENTREZID)

# run KEGG pathway analysis on all DEGs
KEGG_male_basal_vs_10 <- enrichKEGG(gene = entrez_male_basal_vs_10,
                                    universe = entrez_background_male_basal_vs_10,
                                    organism = "mmu", 
                                    pvalueCutoff = 0.05) # no downregulated genes


# define background genes with valid adjusted p-values
background_male_basal_vs_60 <- rownames(res_male_basal_vs_60)[!is.na(res_male_basal_vs_60$padj)]

# convert Ensembl IDs to Entrez IDs
entrez_background_male_basal_vs_60 <- mapIds(org.Mm.eg.db,
                                             keys = background_male_basal_vs_60,
                                             column = "ENTREZID",
                                             keytype = "ENSEMBL",
                                             multiVals = "first")

# remove missing IDs and duplicates
entrez_background_male_basal_vs_60 <- unique(na.omit(unname(entrez_background_male_basal_vs_60)))

# convert to Ensembl IDs to entrez IDs
entrez_upregulated_male_basal_vs_60 <- bitr(upregulated_male_basal_vs_60,
                                            fromType = "ENSEMBL",
                                            toType = "ENTREZID",
                                            OrgDb = "org.Mm.eg.db") %>%
  pull(ENTREZID)

entrez_downregulated_male_basal_vs_60 <- bitr(downregulated_male_basal_vs_60,
                                              fromType = "ENSEMBL",
                                              toType = "ENTREZID",
                                              OrgDb = "org.Mm.eg.db") %>%
  pull(ENTREZID)

# run KEGG pathway analysis on upregulated DEGs
KEGG_upregulated_male_basal_vs_60 <- enrichKEGG(gene = entrez_upregulated_male_basal_vs_60,
                                                universe = entrez_background_male_basal_vs_60,
                                                organism = "mmu", 
                                                pvalueCutoff = 0.05)

# run KEGG pathway analysis on downregulated DEGs
KEGG_downregulated_male_basal_vs_60 <- enrichKEGG(gene = entrez_downregulated_male_basal_vs_60,
                                                  universe = entrez_background_male_basal_vs_60,
                                                  organism = "mmu", 
                                                  pvalueCutoff = 0.05)


# define background genes with valid adjusted p-values
background_male_10_vs_60 <- rownames(res_male_10_vs_60)[!is.na(res_male_10_vs_60$padj)]

# convert Ensembl IDs to Entrez IDs
entrez_background_male_10_vs_60 <- mapIds(org.Mm.eg.db,
                                             keys = background_male_10_vs_60,
                                             column = "ENTREZID",
                                             keytype = "ENSEMBL",
                                             multiVals = "first")

# remove missing IDs and duplicates
entrez_background_male_10_vs_60 <- unique(na.omit(unname(entrez_background_male_10_vs_60)))

# convert to Ensembl IDs to entrez IDs
entrez_upregulated_male_10_vs_60 <- bitr(upregulated_male_10_vs_60,
                                         fromType = "ENSEMBL",
                                         toType = "ENTREZID",
                                         OrgDb = "org.Mm.eg.db") %>%
  pull(ENTREZID)

entrez_downregulated_male_10_vs_60 <- bitr(downregulated_male_10_vs_60,
                                           fromType = "ENSEMBL",
                                           toType = "ENTREZID",
                                           OrgDb = "org.Mm.eg.db") %>%
  pull(ENTREZID)

# run KEGG pathway analysis on upregulated DEGs
KEGG_upregulated_male_10_vs_60 <- enrichKEGG(gene = entrez_upregulated_male_10_vs_60,
                                             universe = entrez_background_male_10_vs_60,
                                             organism = "mmu", 
                                             pvalueCutoff = 0.05)

# run KEGG pathway analysis on downregulated DEGs
KEGG_downregulated_male_10_vs_60 <- enrichKEGG(gene = entrez_downregulated_male_10_vs_60,
                                               universe = entrez_background_male_10_vs_60,
                                               organism = "mmu", 
                                               pvalueCutoff = 0.05)


# define background genes with valid adjusted p-values
background_female_basal_vs_10 <- rownames(res_female_basal_vs_10)[!is.na(res_female_basal_vs_10$padj)]

# convert Ensembl IDs to Entrez IDs
entrez_background_female_basal_vs_10 <- mapIds(org.Mm.eg.db,
                                             keys = background_female_basal_vs_10,
                                             column = "ENTREZID",
                                             keytype = "ENSEMBL",
                                             multiVals = "first")

# remove missing IDs and duplicates
entrez_background_female_basal_vs_10 <- unique(na.omit(unname(entrez_background_female_basal_vs_10)))

# convert to Ensembl IDs to entrez IDs 
entrez_female_basal_vs_10 <- bitr(rownames(sig_DEG_female_basal_vs_10),
                                  fromType = "ENSEMBL",
                                  toType = "ENTREZID",
                                  OrgDb = "org.Mm.eg.db") %>%
  pull(ENTREZID)

# run KEGG pathway analysis on all DEGs
KEGG_female_basal_vs_10 <- enrichKEGG(gene = entrez_female_basal_vs_10,
                                      universe = entrez_background_female_basal_vs_10,
                                      organism = "mmu", 
                                      pvalueCutoff = 0.05) # no downregulated genes


# define background genes with valid adjusted p-values
background_female_basal_vs_60 <- rownames(res_female_basal_vs_60)[!is.na(res_female_basal_vs_60$padj)]

# convert Ensembl IDs to Entrez IDs
entrez_background_female_basal_vs_60 <- mapIds(org.Mm.eg.db,
                                               keys = background_female_basal_vs_60,
                                               column = "ENTREZID",
                                               keytype = "ENSEMBL",
                                               multiVals = "first")

# remove missing IDs and duplicates
entrez_background_female_basal_vs_60 <- unique(na.omit(unname(entrez_background_female_basal_vs_60)))

# convert to Ensembl IDs to entrez IDs
entrez_upregulated_female_basal_vs_60 <- bitr(upregulated_female_basal_vs_60,
                                              fromType = "ENSEMBL",
                                              toType = "ENTREZID",
                                              OrgDb = "org.Mm.eg.db") %>%
  pull(ENTREZID)

entrez_downregulated_female_basal_vs_60 <- bitr(downregulated_female_basal_vs_60,
                                              fromType = "ENSEMBL",
                                              toType = "ENTREZID",
                                              OrgDb = "org.Mm.eg.db") %>%
  pull(ENTREZID)


# run KEGG pathway analysis on upregulated DEGs
KEGG_upregulated_female_basal_vs_60 <- enrichKEGG(gene = entrez_upregulated_female_basal_vs_60,
                                                  universe = entrez_background_female_basal_vs_60,
                                                  organism = "mmu", 
                                                  pvalueCutoff = 0.05)

# run KEGG pathway analysis on downregulated DEGs
KEGG_downregulated_female_basal_vs_60 <- enrichKEGG(gene = entrez_downregulated_female_basal_vs_60,
                                                    universe = entrez_background_female_basal_vs_60,
                                                    organism = "mmu", 
                                                    pvalueCutoff = 0.05)

# define background genes with valid adjusted p-values
background_female_10_vs_60 <- rownames(res_female_10_vs_60)[!is.na(res_female_10_vs_60$padj)]

# convert Ensembl IDs to Entrez IDs
entrez_background_female_10_vs_60 <- mapIds(org.Mm.eg.db,
                                            keys = background_female_10_vs_60,
                                            column = "ENTREZID",
                                            keytype = "ENSEMBL",
                                            multiVals = "first")

# remove missing IDs and duplicates
entrez_background_female_10_vs_60 <- unique(na.omit(unname(entrez_background_female_10_vs_60)))

# convert to Ensembl IDs to entrez IDs
entrez_upregulated_female_10_vs_60 <- bitr(upregulated_female_10_vs_60,
                                           fromType = "ENSEMBL",
                                           toType = "ENTREZID",
                                           OrgDb = "org.Mm.eg.db") %>%
  pull(ENTREZID)

entrez_downregulated_female_10_vs_60 <- bitr(downregulated_female_10_vs_60,
                                             fromType = "ENSEMBL",
                                             toType = "ENTREZID",
                                             OrgDb = "org.Mm.eg.db") %>%
  pull(ENTREZID)

# run KEGG pathway analysis on upregulated DEGs
KEGG_upregulated_female_10_vs_60 <- enrichKEGG(gene = entrez_upregulated_female_10_vs_60,
                                               universe = entrez_background_female_10_vs_60,
                                               organism = "mmu", 
                                               pvalueCutoff = 0.05)

# run KEGG pathway analysis on downregulated DEGs
KEGG_downregulated_female_10_vs_60 <- enrichKEGG(gene = entrez_downregulated_female_10_vs_60,
                                                 universe = entrez_background_female_10_vs_60,
                                                 organism = "mmu", 
                                                 pvalueCutoff = 0.05)
