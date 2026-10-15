# load required packages
library(DESeq2)
library(dplyr)

# change DESeq2 design to condition only
design(dds) <- ~ condition

# relevel so basal is the reference
dds$condition <- relevel(dds$condition, ref = "basal")

# run DESeq2
dds_condition <- DESeq(dds)

# check coefficient for correct comparison
resultsNames(dds_condition) 

# generate results object for basal vs 60 min
res_basal_vs_60 <- results(dds_condition, 
                           contrast = c("condition", "60", "basal"),
                           alpha = 0.05)

summary(res_basal_vs_60) 

# convert to DESeq2 results to dataframe
res_basal_vs_60_df <- as.data.frame(res_basal_vs_60)

# create column to label down or upregulated DEGs
res_basal_vs_60_df <- res_basal_vs_60_df %>%
  mutate(diffaccessible = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "UP",
                                    log2FoldChange <= -0.5 & padj <= 0.05 ~ "DOWN",
                                    TRUE ~ "NO"))

# see how many DARs there are
table(res_basal_vs_60_df$diffaccessible)

# get significant DARs
sig_basal_vs_60 <- res_basal_vs_60_df %>% 
  filter(diffaccessible != "NO")

# subset data by basal
dds_basal <- dds[, dds$condition == "basal"]
dds_basal$condition <- droplevels(dds_basal$condition)

# change DESeq2 design to condition only
design(dds_basal) <- ~ sex

# perform DESEq2 on basal subsetted dataset
dds_basal <- DESeq(dds_basal)

# check coefficient for correct comparison
resultsNames(dds_basal) 

# generate results object for basal male vs female
res_basal_male_vs_female <- results(dds_basal, 
                                    contrast = c("sex", "male", "female"),
                                    alpha = 0.05)
summary(res_basal_male_vs_female) 

# convert to DESeq2 results to dataframe
res_basal_male_vs_female_df <- as.data.frame(res_basal_male_vs_female)

# create column to label down or upregulated DEGs
res_basal_male_vs_female_df <- res_basal_male_vs_female_df %>%
  mutate(diffaccessible = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "UP",
                                    log2FoldChange <= -0.5 & padj <= 0.05 ~ "DOWN",
                                    TRUE ~ "NO"))

# see how many DARs there are
table(res_basal_male_vs_female_df$diffaccessible)

# get significant DARs
sig_basal_male_vs_female <- res_basal_male_vs_female_df %>% 
  filter(diffaccessible != "NO")

# subset data by 60 min ARS samples
dds_60min <- dds[, dds$condition == "60"]
dds_60min$condition <- droplevels(dds_60min$condition)

# change DESeq2 design to condition only
design(dds_60min) <- ~ sex

# perform DESEq2 on 60 min subsetted dataset
dds_60min <- DESeq(dds_60min)

# check coefficient for correct comparison
resultsNames(dds_60min) 

# generate results object for 60 min ARS male vs female 
res_60min_male_vs_female <- results(dds_60min, 
                                    contrast = c("sex", "male", "female"),
                                    alpha = 0.05)
summary(res_60min_male_vs_female) 

# convert to DESeq2 results to dataframe
res_60min_male_vs_female_df <- as.data.frame(res_60min_male_vs_female)

# create column to label down or upregulated DEGs
res_60min_male_vs_female_df <- res_60min_male_vs_female_df %>%
  mutate(diffaccessible = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "UP",
                                    log2FoldChange <= -0.5 & padj <= 0.05 ~ "DOWN",
                                    TRUE ~ "NO"))

# see how many DARs there are
table(res_60min_male_vs_female_df$diffaccessible)

# get significant DARs
sig_60min_male_vs_female <- res_60min_male_vs_female_df %>% 
  filter(diffaccessible != "NO")

# subset data by male samples
dds_male <- dds[, dds$sex == "male"]
dds_male$condition <- droplevels(dds_male$condition)

# change DESeq2 design to condition only
design(dds_male) <- ~ condition

# perform DESEq2 on male subsetted dataset
dds_male <- DESeq(dds_male)

# check coefficient for correct comparison
resultsNames(dds_male) 

# generate results object for male basal vs 60 min
res_male_basal_vs_60 <- results(dds_male, 
                                    contrast = c("condition", "60", "basal"),
                                    alpha = 0.05)
summary(res_male_basal_vs_60) 

# convert to DESeq2 results to dataframe
res_male_basal_vs_60_df <- as.data.frame(res_male_basal_vs_60)

# create column to label down or upregulated DEGs
res_male_basal_vs_60_df <- res_male_basal_vs_60_df %>%
  mutate(diffaccessible = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "UP",
                                    log2FoldChange <= -0.5 & padj <= 0.05 ~ "DOWN",
                                    TRUE ~ "NO"))

# see how many DARs there are
table(res_male_basal_vs_60_df$diffaccessible)

# get significant DARs
sig_male_basal_vs_60 <- res_male_basal_vs_60_df %>% 
  filter(diffaccessible != "NO")
                           
# subset data by female samples
dds_female <- dds[, dds$sex == "female"]
dds_female$condition <- droplevels(dds_female$condition)

# change DESeq2 design to condition only
design(dds_female) <- ~ condition

# perform DESEq2 on female subsetted dataset
dds_female <- DESeq(dds_female)

# check coefficient for correct comparison
resultsNames(dds_female) 

# generate results object for female basal vs 60 min
res_female_basal_vs_60 <- results(dds_female, 
                                contrast = c("condition", "60", "basal"),
                                alpha = 0.05)
summary(res_female_basal_vs_60) 

# convert to DESeq2 results to dataframe
res_female_basal_vs_60_df <- as.data.frame(res_female_basal_vs_60)

# create column to label down or upregulated DEGs
res_female_basal_vs_60_df <- res_female_basal_vs_60_df %>%
  mutate(diffaccessible = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "UP",
                                    log2FoldChange <= -0.5 & padj <= 0.05 ~ "DOWN",
                                    TRUE ~ "NO"))

# see how many DARs there are
table(res_female_basal_vs_60_df$diffaccessible)

# get significant DARs
sig_female_basal_vs_60 <- res_female_basal_vs_60_df %>% 
  filter(diffaccessible != "NO")

# set interaction model for design
design(dds) <- ~ sex * condition
dds <- DESeq(dds) 

# check coefficient for correct comparison
resultsNames(dds) 

# extract interaction terms for basal vs 10 min
res_interaction_60 <- results(dds,
                              contrast = list("sexfemale.condition60"),
                              alpha = 0.05)
summary(res_interaction_60) # no significant genes
