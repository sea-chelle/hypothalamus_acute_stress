# load required packages
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(tibble)
library(AnnotationDbi)

# create vector for names of each GSEA comparison
GSEA_comparisons <- c("male_basal_vs_10", "male_basal_vs_60", "male_10_vs_60", 
                      "female_basal_vs_10", "female_basal_vs_60", "female_10_vs_60")

# run for loop for GSEA GO 
for (comparison in GSEA_comparisons) {
  
  # grab the DESeq2 object
  res <- get(paste0("res_", comparison))
  
  # convert DESeq2 object into df and map Ensembl -> ENTERZ ID
  res_df <- as.data.frame(res) %>%
    rownames_to_column("ENSEMBL") %>%
    left_join(bitr(rownames(res),
                   fromType = "ENSEMBL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db),
              by = "ENSEMBL")
  
  # make the ranked gene list 
  gene_list <- res_df %>%
    filter(!is.na(stat), !is.na(ENTREZID)) %>% # keep genes with both a DESew2 Wald stat and an ENTERZ ID
    arrange(desc(abs(stat))) %>% # temporarily sort genes by the absolute magnitude of the Wald stat
    distinct(ENTREZID, .keep_all = TRUE) %>% # one row per Enterz ID
    select(ENTREZID, stat) %>% # if multiple rows correspond to the same ENTEREZ ID, retain the one with the strongest Wald stat
    deframe() %>%
    sort(decreasing = TRUE) # signed ranking for GSEA
  
  # run GSEA GO
  GSEA_result <- gseGO(geneList = gene_list,
                       OrgDb = org.Mm.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       minGSSize = 10,
                       maxGSSize = 500,
                       pvalueCutoff = 1,
                       pAdjustMethod = "BH",
                       eps = 0,
                       verbose = FALSE)
  
  # store GSEA results in a new variable
  assign(paste0("GSEA_GO_", comparison),
         GSEA_result,
         envir = .GlobalEnv)
  
  # remove temporary objects
  rm(res, res_df, gene_list, GSEA_result, comparison)
  }

# create one table containing every significant GO term from all analyses
GSEA_GO_all <- bind_rows(GSEA_GO_male_basal_vs_10@result %>%
                           mutate(Sex = "Male", Comparison = "Basal vs 10"),
                         
                         GSEA_GO_male_basal_vs_60@result %>%
                           mutate(Sex = "Male", Comparison = "Basal vs 60"),
                         
                         GSEA_GO_male_10_vs_60@result %>%
                           mutate(Sex = "Male", Comparison = "10 vs 60"),
                         
                         GSEA_GO_female_basal_vs_10@result %>%
                           mutate(Sex = "Female", Comparison = "Basal vs 10"),
                         
                         GSEA_GO_female_basal_vs_60@result %>%
                           mutate(Sex = "Female", Comparison = "Basal vs 60"),
                         
                         GSEA_GO_female_10_vs_60@result %>%
                           mutate(Sex = "Female", Comparison = "10 vs 60")) %>%
  filter(p.adjust <= 0.05)

# see how many significant terms for each comparison
GSEA_GO_all %>% 
  dplyr::count(Sex, Comparison)

# see how many positive and negative terms for each comparison
GSEA_GO_all %>%
  mutate(Direction = ifelse(NES > 0, "Positive NES", "Negative NES")) %>%
  dplyr::count(Sex, Comparison, Direction)

# run for loop to get significant terms only
for (comparison in GSEA_comparisons) {
  
  # grab GSEA result
  GSEA_result <- get(paste0("GSEA_GO_", comparison))
  
  # keep significant GSEA terms only
  GSEA_sig <- GSEA_result
  GSEA_sig@result <- GSEA_sig@result %>%
    filter(p.adjust <= 0.05)
  
  # remove redundant GO terms
  GSEA_simplified <- simplify(GSEA_sig,
                              cutoff = 0.7,
                              by = "p.adjust",
                              select_fun = min)
  
  # store simplified GSEA results in a new object
  assign(paste0("GSEA_GO_", comparison, "_simplified"),
         GSEA_simplified,
         envir = .GlobalEnv)
  
  # remove temporary objects
  rm(GSEA_result, GSEA_sig, GSEA_simplified, comparison)
  }

# check how many postivie and negative terms are left
bind_rows(GSEA_GO_male_basal_vs_10_simplified@result %>%
            mutate(Sex = "Male", Comparison = "Basal vs 10"),
          GSEA_GO_male_basal_vs_60_simplified@result %>%
            mutate(Sex = "Male", Comparison = "Basal vs 60"),
          GSEA_GO_male_10_vs_60_simplified@result %>%
            mutate(Sex = "Male", Comparison = "10 vs 60"),
          GSEA_GO_female_basal_vs_10_simplified@result %>%
            mutate(Sex = "Female", Comparison = "Basal vs 10"),
          GSEA_GO_female_basal_vs_60_simplified@result %>%
            mutate(Sex = "Female", Comparison = "Basal vs 60"),
          GSEA_GO_female_10_vs_60_simplified@result %>%
            mutate(Sex = "Female", Comparison = "10 vs 60")) %>%
  mutate(Direction = ifelse(NES > 0, "Positive NES", "Negative NES")) %>%
  dplyr::count(Sex, Comparison, Direction)

# create vector for names of each GSEA comparison
GSEA_comparisons <- c("basal_male_vs_female", "10min_male_vs_female", "60min_male_vs_female")

# run for loop for GSEA GO 
for (comparison in GSEA_comparisons) {
  
  # grab the DESeq2 object
  res <- get(paste0("res_", comparison))
  
  # convert DESeq2 object into df and map Ensembl -> ENTERZ ID
  res_df <- as.data.frame(res) %>%
    rownames_to_column("ENSEMBL") %>%
    left_join(bitr(rownames(res),
                   fromType = "ENSEMBL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db),
              by = "ENSEMBL")
  
  # make the ranked gene list 
  gene_list <- res_df %>%
    filter(!is.na(stat), !is.na(ENTREZID)) %>% # keep genes with both a DESew2 Wald stat and an ENTERZ ID
    arrange(desc(abs(stat))) %>% # temporarily sort genes by the absolute magnitude of the Wald stat
    distinct(ENTREZID, .keep_all = TRUE) %>% # one row per Enterz ID
    select(ENTREZID, stat) %>% # if multiple rows correspond to the same ENTEREZ ID, retain the one with the strongest Wald stat
    deframe() %>%
    sort(decreasing = TRUE) # signed ranking for GSEA
  
  # run GSEA GO
  GSEA_result <- gseGO(geneList = gene_list,
                       OrgDb = org.Mm.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       minGSSize = 10,
                       maxGSSize = 500,
                       pvalueCutoff = 1,
                       pAdjustMethod = "BH",
                       eps = 0,
                       verbose = FALSE)
  
  # store GSEA results in a new variable
  assign(paste0("GSEA_GO_", comparison),
         GSEA_result,
         envir = .GlobalEnv)
  
  # remove temporary objects
  rm(res, res_df, gene_list, GSEA_result, comparison)
}

# create one table containing every significant GO term from all analyses
GSEA_GO_all <- bind_rows(GSEA_GO_basal_male_vs_female@result %>%
                           mutate(Comparison = "Basal"),
                         
                         GSEA_GO_10min_male_vs_female@result %>%
                           mutate(Comparison = "10 min"),
                         
                         GSEA_GO_60min_male_vs_female@result %>%
                           mutate(Comparison = "60 min")) %>%
  filter(p.adjust <= 0.05)

# see how many significant terms for each comparison
GSEA_GO_all %>% 
  dplyr::count(Comparison)

# see how many positive and negative terms for each comparison
GSEA_GO_all %>%
  mutate(Direction = ifelse(NES > 0, "Positive NES", "Negative NES")) %>%
  dplyr::count(Comparison, Direction)

= 0.05, NES < 0) %>%
    pull(Description),

# create vector for names of each GSEA comparison
GSEA_comparisons <- c("interaction_10", "interaction_60")

# run for loop for GSEA GO 
for (comparison in GSEA_comparisons) {
  
  # grab the DESeq2 object
  res <- get(paste0("res_", comparison))
  
  # convert DESeq2 object into df and map Ensembl -> ENTERZ ID
  res_df <- as.data.frame(res) %>%
    rownames_to_column("ENSEMBL") %>%
    left_join(bitr(rownames(res),
                   fromType = "ENSEMBL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db),
              by = "ENSEMBL")
  
  # make the ranked gene list 
  gene_list <- res_df %>%
    filter(!is.na(stat), !is.na(ENTREZID)) %>% # keep genes with both a DESew2 Wald stat and an ENTERZ ID
    arrange(desc(abs(stat))) %>% # temporarily sort genes by the absolute magnitude of the Wald stat
    distinct(ENTREZID, .keep_all = TRUE) %>% # one row per Enterz ID
    select(ENTREZID, stat) %>% # if multiple rows correspond to the same ENTEREZ ID, retain the one with the strongest Wald stat
    deframe() %>%
    sort(decreasing = TRUE) # signed ranking for GSEA
  
  # run GSEA GO
  GSEA_result <- gseGO(geneList = gene_list,
                       OrgDb = org.Mm.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       minGSSize = 10,
                       maxGSSize = 500,
                       pvalueCutoff = 1,
                       pAdjustMethod = "BH",
                       eps = 0,
                       verbose = FALSE)
  
  # store GSEA results in a new variable
  assign(paste0("GSEA_GO_", comparison),
         GSEA_result,
         envir = .GlobalEnv)
  
  # remove temporary objects
  rm(res, res_df, gene_list, GSEA_result, comparison)
  }

# create one table containing every significant GO term from all analyses
GSEA_GO_all <- bind_rows(GSEA_GO_interaction_10@result %>%
                           mutate(Comparison = "10 min"),
                         
                         GSEA_GO_interaction_60@result %>%
                           mutate(Comparison = "60 min")) %>%
  
  filter(p.adjust <= 0.05)

# see how many significant terms for each comparison
GSEA_GO_all %>% 
  dplyr::count(Comparison)

# see how many positive and negative terms for each comparison
GSEA_GO_all %>%
  mutate(Direction = ifelse(NES > 0, "Positive NES", "Negative NES")) %>%
  dplyr::count(Comparison, Direction)
  pull(Description),
