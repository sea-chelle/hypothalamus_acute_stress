# load required packages
library(here)
library(tidyverse)

# grab RNA-seq file names
rna_files <- list.files(file.path(here("RNA_seq")),
                        pattern = "^rna_.*\\.rds$",
                        full.names = TRUE)

# import RNA-seq DESeq2 df
rna_results <- lapply(rna_files, readRDS)

# set names for RNA-seq results
names(rna_results) <- gsub("^rna_|\\.rds$", "",
                           basename(rna_files))
names(rna_results)

# see how many ENSEMBL ID genes are in both RNA-seq and DAR annotated genes
lapply(names(rna_results), function(comparison) {
  
  atac_subset <- master_dar %>%
    filter(.data$comparison == .env$comparison)
  
  data.frame(comparison = comparison,
             n_DARs = nrow(atac_subset),
             n_ENSEMBL_present = sum(!is.na(atac_subset$ENSEMBL) &
                                       atac_subset$ENSEMBL != ""),
    n_RNA_matched = sum(atac_subset$ENSEMBL %in% rna_results[[comparison]]$ensembl,
                        na.rm = TRUE))
}) %>%
  bind_rows()

# determine whether ChIPseeker nearest annotated TTS is a diff expressed gene from RNA-seq data
nearest_rna <- list()

for (comparison in names(rna_results)) {
  
  nearest_rna[[comparison]] <- master_dar %>%
    filter(.data$comparison == .env$comparison) %>%
    left_join(rna_results[[comparison]] %>%
        dplyr::select(ensembl,
                      rna_symbol = symbol,
                      rna_log2FoldChange = log2FoldChange,
                      rna_padj = padj,
                      rna_diffexpressed = diffexpressed),
        by = c("ENSEMBL" = "ensembl"))
}

# check how many ChIPseeker annotated target genes for each comparison
sapply(nearest_rna, nrow)

# see how many DAR target from ChIPseeker nearest annotated TSS are contained in RNA-seq genes
nearest_rna_summary <- lapply(nearest_rna, function(x) {
  
  x %>%
    summarise(n_DAR_target_pairs = n(),
              n_DARs = n_distinct(peak),
              n_unique_targets = dplyr::n_distinct(ENSEMBL),
              n_RNA_matched_pairs = sum(!is.na(rna_log2FoldChange)),
              n_RNA_matched_targets = dplyr::n_distinct(ENSEMBL[!is.na(rna_log2FoldChange)]))
})

nearest_rna_summary

# retain DARs whose nearest gene is also differentially expressed
nearest_deg <- lapply(nearest_rna, function(x) {
  
  x %>%
    filter(!is.na(rna_diffexpressed),
           rna_diffexpressed != "NO") %>%
    dplyr::select(peak, ENSEMBL, SYMBOL, region, log2FoldChange, padj,
                  rna_symbol, rna_log2FoldChange, rna_padj, rna_diffexpressed,
                  H3K27ac, PLAC, n_PLAC_targets, starts_with("CRE_"))
})

nearest_deg_summary <- lapply(nearest_deg, function(x) {
  
  x %>%
    summarise(
      n_DARs = n_distinct(peak),
      n_unique_DEGs = n_distinct(rna_symbol),
      n_concordant = sum(sign(log2FoldChange) == sign(rna_log2FoldChange)))
})

nearest_deg_summary

# check which DEGs have DAR annotated to corresponding nearest TTS
nearest_deg$basal_vs_60
nearest_deg$male_basal_vs_60
nearest_deg$female_basal_vs_60

# integrate PLAC-linked target genes with RNA-seq
plac_rna <- list()

for (comparison in names(rna_results)) {
  
  plac_rna[[comparison]] <- comparisons[[comparison]]$plac_peak_gene %>%
    
    # add ATAC / regulatory annotation for each DAR
    left_join(
      master_dar %>%
        filter(.data$comparison == .env$comparison) %>%
        dplyr::select(peak, region, log2FoldChange, padj,
                      H3K27ac, starts_with("CRE_")), 
      by = "peak") %>%
    
    # add RNA result for each PLAC-linked target
    left_join(rna_results[[comparison]] %>%
        dplyr::select(ensembl,
                      rna_symbol = symbol,
                      rna_log2FoldChange = log2FoldChange,
                      rna_padj = padj,
                      rna_diffexpressed = diffexpressed),
      by = c("ENSEMBL_ID" = "ensembl"))
}

# check how many PLAC-seq annotated target genes for each comparison
sapply(plac_rna, nrow)

# see how many DAR target from PLAC-seq are contained in RNA-seq genes
plac_rna_summary <- lapply(plac_rna, function(x) {
  
  x %>%
    summarise(n_DAR_target_pairs = n(),
              n_DARs = n_distinct(peak),
              n_unique_targets = dplyr::n_distinct(ENSEMBL_ID),
              n_RNA_matched_pairs = sum(!is.na(rna_log2FoldChange)),
              n_RNA_matched_targets = dplyr::n_distinct(ENSEMBL_ID[!is.na(rna_log2FoldChange)]))
})

plac_rna_summary

# retain PLAC-linked targets that are also differentially expressed
plac_deg <- lapply(plac_rna, function(x) {
  
  x %>%
    filter(!is.na(rna_diffexpressed), rna_diffexpressed != "NO") %>%
    dplyr::select(peak, GeneName, ENSEMBL_ID, region, log2FoldChange, padj,
                  rna_symbol, rna_log2FoldChange, rna_padj, rna_diffexpressed,
                  H3K27ac, starts_with("CRE_"))
})

plac_deg_summary <- lapply(plac_deg, function(x) {
  
  x %>%
    summarise(n_DAR_target_pairs = n(),
              n_DARs = n_distinct(peak),
              n_unique_DEGs = n_distinct(ENSEMBL_ID), 
              n_concordant = sum(sign(log2FoldChange) == sign(rna_log2FoldChange)))
})

plac_deg_summary

# check which DEGs have DAR annotated to PLAC-seq associated contact
plac_deg$basal_vs_60
plac_deg$male_basal_vs_60
plac_deg$female_basal_vs_60

# combine nearest gene, PLAC-seq contact genes and RNA-seq data together 
nearest_vs_plac <- lapply(names(rna_results), function(comparison) {
  
  nearest_genes <- unique(nearest_deg[[comparison]]$rna_symbol)
  plac_genes <- unique(plac_deg[[comparison]]$rna_symbol)
  
  tibble(comparison = comparison,
         nearest_n = length(nearest_genes),
         plac_n = length(plac_genes),
         shared_n = length(intersect(nearest_genes, plac_genes)),
         nearest_only_n = length(setdiff(nearest_genes, plac_genes)),
         plac_only_n = length(setdiff(plac_genes, nearest_genes)),
         shared = paste(intersect(nearest_genes, plac_genes), collapse = ", "),
         nearest_only = paste(setdiff(nearest_genes, plac_genes), collapse = ", "),
         plac_only = paste(setdiff(plac_genes, nearest_genes), collapse = ", "))
}) %>%
  bind_rows()

nearest_vs_plac
