# load required packages
library(DESeq2)
library(dplyr)
library(AnnotationDbi)
library(org.Mm.eg.db)

# subset data by male
dds_male <- dds[, dds$sex == "m"]
dds_male$sex <- droplevels(dds_male$sex) # drop unused levels / "f"

# change DESeq2 design to condition only
design(dds_male) <- ~ condition

# relevel so basal is the reference
dds$condition <- relevel(dds$condition, ref = "basal")

# perform DESEq2 on male subsetted dataset
dds_male <- DESeq(dds_male)

# check coefficient for correct comparison
resultsNames(dds_male) 

# generate results object for basal vs 10 min
res_male_basal_vs_10 <- results(dds_male, 
                                contrast = c("condition", "10", "basal"),
                                alpha = 0.05)
summary(res_male_basal_vs_10) 

# convert to DESeq2 results to dataframe
res_male_basal_vs_10_df <- as.data.frame(res_male_basal_vs_10)

# create column to label down or upregulated DEGs
res_male_basal_vs_10_df <- res_male_basal_vs_10_df %>%
  mutate(diffexpressed = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "UP",
                                   log2FoldChange <= -0.5 & padj <= 0.05 ~ "DOWN",
                                   TRUE ~ "NO"))
# create column with ensembl IDs
res_male_basal_vs_10_df$ensembl <- rownames(res_male_basal_vs_10_df)

# map ENSEMBL IDs to gene symbols
res_male_basal_vs_10_df$symbol <- mapIds(org.Mm.eg.db,
                                         keys = res_male_basal_vs_10_df$ensembl,
                                         column = "SYMBOL",
                                         keytype = "ENSEMBL",
                                         multiVals = "first")

# generate results object for basal vs 10 min
res_male_basal_vs_60 <- results(dds_male, 
                                contrast = c("condition", "60", "basal"),
                                alpha = 0.05)
summary(res_male_basal_vs_60) 

# convert to DESeq2 results to dataframe
res_male_basal_vs_60_df <- as.data.frame(res_male_basal_vs_60)

# create column to label down or upregulated DEGs
res_male_basal_vs_60_df <- res_male_basal_vs_60_df %>%
  mutate(diffexpressed = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "UP",
                                   log2FoldChange <= -0.5 & padj <= 0.05 ~ "DOWN",
                                   TRUE ~ "NO"))

# create column with ensembl IDs
res_male_basal_vs_60_df$ensembl <- rownames(res_male_basal_vs_60_df)

# map ENSEMBL IDs to gene symbols
res_male_basal_vs_60_df$symbol <- mapIds(org.Mm.eg.db,
                                         keys = res_male_basal_vs_60_df$ensembl,
                                         column = "SYMBOL",
                                         keytype = "ENSEMBL",
                                         multiVals = "first")

# generate results object for basal vs 10 min
res_male_10_vs_60 <- results(dds_male, 
                             contrast = c("condition", "60", "10"),
                             alpha = 0.05)
summary(res_male_10_vs_60) 

# convert to DESeq2 results to dataframe
res_male_10_vs_60_df <- as.data.frame(res_male_10_vs_60)

# create column to label down or upregulated DEGs
res_male_10_vs_60_df <- res_male_10_vs_60_df %>%
  mutate(diffexpressed = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "UP",
                                   log2FoldChange <= -0.5 & padj <= 0.05 ~ "DOWN",
                                   TRUE ~ "NO"))

# create column with ensembl IDs
res_male_10_vs_60_df$ensembl <- rownames(res_male_10_vs_60_df)

# map ENSEMBL IDs to gene symbols
res_male_10_vs_60_df$symbol <- mapIds(org.Mm.eg.db,
                                      keys = res_male_10_vs_60_df$ensembl,
                                      column = "SYMBOL",
                                      keytype = "ENSEMBL",
                                      multiVals = "first")
