# load required packages
library(clusterProfiler)
library(org.Mm.eg.db)

# subset dataframe to only include significant DEGs 
sig_DEG_basal_vs_10 <- sig_DEG_basal_vs_10[sig_DEG_basal_vs_10$diffexpressed != "NO",] # remove non-significant DEGs

# run GO on all significant DEGs, irrespective of direction
GO_basal_vs_10 <- enrichGO(gene = rownames(sig_DEG_basal_vs_10),
                           universe = rownames(res_basal_vs_10),
                           OrgDb = "org.Mm.eg.db",
                           keyType = "ENSEMBL",
                           ont = "BP")

# subset dataframe to only include significant DEGs 
sig_DEG_basal_vs_60 <- sig_DEG_basal_vs_60[sig_DEG_basal_vs_60$diffexpressed != "NO",] 

# get list of upregulated and downregulated genes to test separately 
genes_to_test <- split(sig_DEG_basal_vs_60, sig_DEG_basal_vs_60$diffexpressed) 
upregulated_basal_vs_60 <- rownames(genes_to_test$UP) 
downregulated_basal_vs_60 <- rownames(genes_to_test$DOWN) 

# run GO on unregulated genes only
GO_upregulated_basal_vs_60 <- enrichGO(gene = upregulated_basal_vs_60, 
                                       universe = rownames(res_basal_vs_60),
                                       OrgDb = "org.Mm.eg.db",
                                       keyType = "ENSEMBL",
                                       ont = "BP")

# run gene ontology on downregulated genes only
GO_downregulated_basal_vs_60 <- enrichGO(gene = downregulated_basal_vs_60, 
                                         universe = rownames(res_basal_vs_60),
                                         OrgDb = "org.Mm.eg.db",
                                         keyType = "ENSEMBL",
                                         ont = "BP")

# subset dataframe to only include significant DEGs 
sig_DEG_10_vs_60 <- sig_DEG_10_vs_60[sig_DEG_10_vs_60$diffexpressed != "NO",]


# get list of upregulated and downregulated genes to test separately 
genes_to_test <- split(sig_DEG_10_vs_60, sig_DEG_10_vs_60$diffexpressed) 
downregulated_10_vs_60 <- rownames(genes_to_test$DOWN)

# run GO on unregulated genes only
GO_upregulated_10_vs_60 <- enrichGO(gene = upregulated_10_vs_60,
                                    universe = rownames(res_10_vs_60),
                                    OrgDb = "org.Mm.eg.db",
                                    keyType = "ENSEMBL", 
                                    ont = "BP")

# run gene ontology on downregulated genes only
GO_downregulated_10_vs_60 <- enrichGO(gene = downregulated_10_vs_60, 
                                      universe = rownames(res_10_vs_60),
                                      OrgDb = "org.Mm.eg.db",
                                      keyType = "ENSEMBL",
                                      ont = "BP")

                                       mutate(Comparison = "Basal vs 60 min ARS"),
                                     GO_upregulated_10_vs_60_plot %>%
                                       mutate(Comparison = "10 vs 60 min ARS"))

# subset dataframe to only include significant DEGs 
sig_DEG_male_basal_vs_10 <- sig_DEG_male_basal_vs_10[sig_DEG_male_basal_vs_10$diffexpressed != "NO",]

# run GO on all significant DEGs, irrespective of direction
GO_male_basal_vs_10 <- enrichGO(gene = rownames(sig_DEG_male_basal_vs_10),
                                universe = rownames(res_male_basal_vs_10),
                                OrgDb = "org.Mm.eg.db",
                                keyType = "ENSEMBL",
                                ont = "BP")

# subset dataframe to only include significant DEGs 
sig_DEG_male_basal_vs_60 <- sig_DEG_male_basal_vs_60[sig_DEG_male_basal_vs_60$diffexpressed != "NO",] 

# get list of upregulated and downregulated genes to test separately 
genes_to_test <- split(sig_DEG_male_basal_vs_60, sig_DEG_male_basal_vs_60$diffexpressed) 
upregulated_male_basal_vs_60 <- rownames(genes_to_test$UP)
downregulated_male_basal_vs_60 <- rownames(genes_to_test$DOWN)

# run gene ontology on upregulated genes only
GO_upregulated_male_basal_vs_60 <- enrichGO(gene = upregulated_male_basal_vs_60, 
                                            universe = rownames(res_male_basal_vs_60),
                                            OrgDb = "org.Mm.eg.db",
                                            keyType = "ENSEMBL",
                                            ont = "BP")

# run gene ontology on downregulated genes only
GO_downregulated_male_basal_vs_60 <- enrichGO(gene = downregulated_male_basal_vs_60,
                                              universe = rownames(res_male_basal_vs_60),
                                              OrgDb = "org.Mm.eg.db",
                                              keyType = "ENSEMBL",
                                              ont = "BP")

# subset dataframe to only include significant DEGs 
sig_DEG_male_10_vs_60 <- sig_DEG_male_10_vs_60[sig_DEG_male_10_vs_60$diffexpressed != "NO",]

# get list of upregulated and downregulated genes to test separately 
genes_to_test <- split(sig_DEG_male_10_vs_60, sig_DEG_male_10_vs_60$diffexpressed)
upregulated_male_10_vs_60 <- rownames(genes_to_test$UP) 
downregulated_male_10_vs_60 <- rownames(genes_to_test$DOWN) 

# run GO on unregulated genes only
GO_upregulated_male_10_vs_60 <- enrichGO(gene = upregulated_male_10_vs_60, 
                                         universe = rownames(res_male_10_vs_10),
                                         OrgDb = "org.Mm.eg.db",
                                         keyType = "ENSEMBL",
                                         ont = "BP")

# run GO on unregulated genes only
GO_downregulated_male_10_vs_60 <- enrichGO(gene = downregulated_male_10_vs_60, 
                                           universe = rownames(res_male_10_vs_10),
                                           OrgDb = "org.Mm.eg.db",
                                           keyType = "ENSEMBL",
                                           ont = "BP")

# subset dataframe to only include significant DEGs 
sig_DEG_female_basal_vs_10 <- sig_DEG_female_basal_vs_10[sig_DEG_female_basal_vs_10$diffexpressed != "NO",] 

# run GO on all significant DEGs, irrespective of direction
GO_female_basal_vs_10 <- enrichGO(gene = rownames(sig_DEG_female_basal_vs_10),
                                  universe = rownames(res_female_basal_vs_10),
                                  OrgDb = "org.Mm.eg.db",
                                  keyType = "ENSEMBL",
                                  ont = "BP")

# subset dataframe to only include significant DEGs 
sig_DEG_female_basal_vs_60 <- sig_DEG_female_basal_vs_60[sig_DEG_female_basal_vs_60$diffexpressed != "NO",]


# get list of upregulated and downregulated genes to test separately 
genes_to_test <- split(sig_DEG_female_basal_vs_60, sig_DEG_female_basal_vs_60$diffexpressed) 
upregulated_female_basal_vs_60 <- rownames(genes_to_test$UP)
downregulated_female_basal_vs_60 <- rownames(genes_to_test$DOWN)

# run GO on unregulated genes only
GO_upregulated_female_basal_vs_60 <- enrichGO(gene = upregulated_female_basal_vs_60, 
                                              universe = rownames(res_female_basal_vs_60),
                                              OrgDb = "org.Mm.eg.db",
                                              keyType = "ENSEMBL",
                                              ont = "BP")

# run gene ontology on downregulated genes only
GO_downregulated_female_basal_vs_60 <- enrichGO(gene = downregulated_female_basal_vs_60, 
                                                universe = rownames(res_female_basal_vs_60),
                                                OrgDb = "org.Mm.eg.db",
                                                keyType = "ENSEMBL",
                                                ont = "BP")


# subset dataframe to only include significant DEGs 
sig_DEG_female_10_vs_60 <- sig_DEG_female_10_vs_60[sig_DEG_female_10_vs_60$diffexpressed != "NO",] 

# get list of upregulated and downregulated genes to test separately 
genes_to_test <- split(sig_DEG_female_10_vs_60, sig_DEG_female_10_vs_60$diffexpressed) 
upregulated_female_10_vs_60 <- rownames(genes_to_test$UP) 
downregulated_female_10_vs_60 <- rownames(genes_to_test$DOWN)

# run GO on unregulated genes only
GO_upregulated_female_10_vs_60 <- enrichGO(gene = upregulated_female_10_vs_60, 
                                           universe = rownames(res_female_10_vs_60),
                                           OrgDb = "org.Mm.eg.db",
                                           keyType = "ENSEMBL",
                                           ont = "BP")

# run GO on unregulated genes only
GO_downregulated_female_10_vs_60 <- enrichGO(gene = downregulated_female_10_vs_60, 
                                             universe = rownames(res_female_10_vs_60),
                                             OrgDb = "org.Mm.eg.db",
                                             keyType = "ENSEMBL",
                                             ont = "BP")
