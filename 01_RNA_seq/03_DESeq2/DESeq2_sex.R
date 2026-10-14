# load required packages
library(DESeq2)
library(dplyr)
library(AnnotationDbi)
library(org.Mm.eg.db)

# subset data by basal
dds_basal <- dds[, dds$condition == "basal"]
dds_basal$condition <- droplevels(dds_basal$condition)

# change DESeq2 design to condition only
design(dds_basal) <- ~ sex

# perform DESEq2 on basal subsetted dataset
dds_basal <- DESeq(dds_basal)

# check coefficient for correct comparison
resultsNames(dds_basal) 

# generate results object for basal vs 10 min
res_basal_male_vs_female <- results(dds_basal, 
                                    contrast = c("sex", "m", "f"),
                                    alpha = 0.05)
summary(res_basal_male_vs_female) 

# convert to DESeq2 results to dataframe
res_basal_male_vs_female_df <- as.data.frame(res_basal_male_vs_female)

# create column to label down or upregulated DEGs
res_basal_male_vs_female_df <- res_basal_male_vs_female_df %>%
  mutate(diffexpressed = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "UP",
                                   log2FoldChange <= -0.5 & padj <= 0.05 ~ "DOWN",
                                   TRUE ~ "NO"))
# create column with ensembl IDs
res_basal_male_vs_female_df$ensembl <- rownames(res_basal_male_vs_female_df)

# map ENSEMBL IDs to gene symbols
res_basal_male_vs_female_df$symbol <- mapIds(org.Mm.eg.db,
                                             keys = res_basal_male_vs_female_df$ensembl,
                                             column = "SYMBOL",
                                             keytype = "ENSEMBL",
                                             multiVals = "first")

# subset data by 10 min 
dds_10min <- dds[, dds$condition == "10"]
dds_10min$condition <- droplevels(dds_10min$condition) 

# change DESeq2 design to condition only
design(dds_10min) <- ~ sex

# perform DESEq2 on 10 min subsetted dataset
dds_10min <- DESeq(dds_10min)

# check coefficient for correct comparison
resultsNames(dds_10min) 

# generate results object for male vs female at 10 min
res_10min_male_vs_female <- results(dds_10min, 
                                    contrast = c("sex", "m", "f"),
                                    alpha = 0.05)
summary(res_10min_male_vs_female) 

# convert to DESeq2 results to dataframe
res_10min_male_vs_female_df <- as.data.frame(res_10min_male_vs_female)

# create column to label down or upregulated DEGs
res_10min_male_vs_female_df <- res_10min_male_vs_female_df %>%
  mutate(diffexpressed = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "UP",
                                   log2FoldChange <= -0.5 & padj <= 0.05 ~ "DOWN",
                                   TRUE ~ "NO"))

# create column with ensembl IDs
res_10min_male_vs_female_df$ensembl <- rownames(res_10min_male_vs_female_df)

# map ENSEMBL IDs to gene symbols
res_10min_male_vs_female_df$symbol <- mapIds(org.Mm.eg.db,
                                             keys = res_10min_male_vs_female_df$ensembl,
                                             column = "SYMBOL",
                                             keytype = "ENSEMBL",
                                             multiVals = "first")

# subset data by 60 min
dds_60min <- dds[, dds$condition == "60"]
dds_60min$condition <- droplevels(dds_60min$condition) 

# change DESeq2 design to condition only
design(dds_60min) <- ~ sex

# perform DESEq2 on 60 min subsetted dataset
dds_60min <- DESeq(dds_60min)

# check coefficient for correct comparison
resultsNames(dds_60min) 

# generate results object for 60min vs 10 min
res_60min_male_vs_female <- results(dds_60min, 
                                    contrast = c("sex", "m", "f"),
                                    alpha = 0.05)
summary(res_60min_male_vs_female) 

# convert to DESeq2 results to dataframe
res_60min_male_vs_female_df <- as.data.frame(res_60min_male_vs_female)

# create column to label down or upregulated DEGs
res_60min_male_vs_female_df <- res_60min_male_vs_female_df %>%
  mutate(diffexpressed = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "UP",
                                   log2FoldChange <= -0.5 & padj <= 0.05 ~ "DOWN",
                                   TRUE ~ "NO"))

# create column with ensembl IDs
res_60min_male_vs_female_df$ensembl <- rownames(res_60min_male_vs_female_df)

# map ENSEMBL IDs to gene symbols
res_60min_male_vs_female_df$symbol <- mapIds(org.Mm.eg.db,
                                             keys = res_60min_male_vs_female_df$ensembl,
                                             column = "SYMBOL",
                                             keytype = "ENSEMBL",
                                             multiVals = "first")

# set interaction model for design
design(dds) <- ~ sex * condition
dds <- DESeq(dds) # this asks the question, does the effect of restraint on gene expression depend on sex?

# check coefficient for correct comparison
resultsNames(dds) 

# extract interaction terms for basal vs 10 min
res_interaction_10 <- results(dds,
                              contrast = list("sexm.condition10"),
                              alpha = 0.05) # no significant genes

# extract interaction terms for basal vs 60 min
res_interaction_60 <- results(dds, 
                              name = "sexm.condition60",
                              alpha = 0.05)
summary(res_interaction_60) 

# convert to DESeq2 results to dataframe
res_interaction_60_df <- as.data.frame(res_interaction_60)

# create column to label down or upregulated DEGs
res_interaction_60_df <- res_interaction_60_df %>%
  mutate(diffexpressed = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "UP",
                                   log2FoldChange <= -0.5 & padj <= 0.05 ~ "DOWN",
                                   TRUE ~ "NO"))

# create column with ensembl IDs
res_interaction_60_df$ensembl <- rownames(res_interaction_60_df)

# map ENSEMBL IDs to gene symbols
res_interaction_60_df$symbol <- mapIds(org.Mm.eg.db,
                                       keys = res_interaction_60_df$ensembl,
                                       column = "SYMBOL",
                                       keytype = "ENSEMBL",
                                       multiVals = "first")

# get names of all DEGs
res_interaction_60_df$symbol[res_interaction_60_df$diffexpressed %in% c("UP", "DOWN")]
