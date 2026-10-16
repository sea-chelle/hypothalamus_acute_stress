# load required packages
library(tidyverse)

# create function to make master DAR table
make_master_dar_table <- function(comparison_name, comparisons) {
  
  anno <- comparisons[[comparison_name]]$anno %>%
    mutate(ENSEMBL = sub("\\..*$", "", geneId))
  
  plac_targets <- comparisons[[comparison_name]]$plac_peak_gene %>%
    group_by(peak) %>%
    summarise(PLAC_genes = paste(sort(unique(na.omit(GeneName))), collapse = "; "),
              PLAC_ensembl = paste(sort(unique(na.omit(ENSEMBL_ID))), collapse = "; "),
              .groups = "drop")
  anno %>%
    left_join(plac_targets, by = "peak") %>%
    mutate(comparison = comparison_name)
}

# make master DAR tables
master_dar_tables <- lapply(c("basal_vs_60", "male_basal_vs_60", "female_basal_vs_60"),
                            make_master_dar_table, comparisons = comparisons)

# rename each comparison
names(master_dar_tables) <- c("Overall", "Male", "Female")

# combine all three tables
master_dar <- bind_rows(master_dar_tables, .id = "group")

# keep one row per DAR
master_dar <- master_dar %>%
  distinct(group, peak, .keep_all = TRUE)

# check dimensions and number of peaks
dim(master_dar)
table(master_dar$group)

# inspect
master_dar %>%
  dplyr::select(group, peak, geneId, ENSEMBL, SYMBOL, 
                log2FoldChange, padj, annotation,
                H3K27ac, PLAC, n_PLAC_targets,
                PLAC_genes, starts_with("CRE_")) %>%
  head()
