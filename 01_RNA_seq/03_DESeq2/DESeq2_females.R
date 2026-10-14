# load required packages
library(DESeq2)
library(dplyr)
library(AnnotationDbi)
library(org.Mm.eg.db)

# subset data by female
dds_female <- dds[, dds$sex == "f"]
dds_female$sex <- droplevels(dds_female$sex) # drop unused levels / "f"

# change DESeq2 design to condition only
design(dds_female) <- ~ condition

# relevel so basal is the reference
dds$condition <- relevel(dds$condition, ref = "basal")

# perform DESEq2 on feamle subsetted dataset
dds_female <- DESeq(dds_female)

# check coefficient for correct comparison
resultsNames(dds_female) 

# generate results object for basal vs 10 min
res_female_basal_vs_10 <- results(dds_female, 
                                contrast = c("condition", "10", "basal"),
                                alpha = 0.05)
summary(res_female_basal_vs_10) 

# convert to DESeq2 results to dataframe
res_female_basal_vs_10_df <- as.data.frame(res_female_basal_vs_10)

# create column to label down or upregulated DEGs
res_female_basal_vs_10_df <- res_female_basal_vs_10_df %>%
  mutate(diffexpressed = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "UP",
                                   log2FoldChange <= -0.5 & padj <= 0.05 ~ "DOWN",
                                   TRUE ~ "NO"))

# create column with ensembl IDs
res_female_basal_vs_10_df$ensembl <- rownames(res_female_basal_vs_10_df)

# map ENSEMBL IDs to gene symbols
res_female_basal_vs_10_df$symbol <- mapIds(org.Mm.eg.db,
                                         keys = res_female_basal_vs_10_df$ensembl,
                                         column = "SYMBOL",
                                         keytype = "ENSEMBL",
                                         multiVals = "first")

# generate results object for basal vs 10 min
res_female_basal_vs_60 <- results(dds_female, 
                                contrast = c("condition", "60", "basal"),
                                alpha = 0.05)
summary(res_female_basal_vs_60) 

# convert to DESeq2 results to dataframe
res_female_basal_vs_60_df <- as.data.frame(res_female_basal_vs_60)

# create column to label down or upregulated DEGs
res_female_basal_vs_60_df <- res_female_basal_vs_60_df %>%
  mutate(diffexpressed = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "UP",
                                   log2FoldChange <= -0.5 & padj <= 0.05 ~ "DOWN",
                                   TRUE ~ "NO"))

# create column with ensembl IDs
res_female_basal_vs_60_df$ensembl <- rownames(res_female_basal_vs_60_df)

# map ENSEMBL IDs to gene symbols
res_female_basal_vs_60_df$symbol <- mapIds(org.Mm.eg.db,
                                         keys = res_female_basal_vs_60_df$ensembl,
                                         column = "SYMBOL",
                                         keytype = "ENSEMBL",
                                         multiVals = "first")

# generate results object for basal vs 10 min
res_female_10_vs_60 <- results(dds_female, 
                             contrast = c("condition", "60", "10"),
                             alpha = 0.05)
summary(res_female_10_vs_60) 

# convert to DESeq2 results to dataframe
res_female_10_vs_60_df <- as.data.frame(res_female_10_vs_60)

# create column to label down or upregulated DEGs
res_female_10_vs_60_df <- res_female_10_vs_60_df %>%
  mutate(diffexpressed = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "UP",
                                   log2FoldChange <= -0.5 & padj <= 0.05 ~ "DOWN",
                                   TRUE ~ "NO"))

# create column with ensembl IDs
res_female_10_vs_60_df$ensembl <- rownames(res_female_10_vs_60_df)

# map ENSEMBL IDs to gene symbols
res_female_10_vs_60_df$symbol <- mapIds(org.Mm.eg.db,
                                      keys = res_female_10_vs_60_df$ensembl,
                                      column = "SYMBOL",
                                      keytype = "ENSEMBL",
                                      multiVals = "first")
