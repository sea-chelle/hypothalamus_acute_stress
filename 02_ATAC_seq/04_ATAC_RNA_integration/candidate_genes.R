# load required packages
library(tidyverse)
library(GenomicRanges)
library(IRanges)

# combine full male and female ATAC results
atac_full <- bind_rows(res_male_basal_vs_60_df %>%
                         rownames_to_column("peak") %>%
                         mutate(sex = "Male"),
                       res_female_basal_vs_60_df %>%
                         rownames_to_column("peak") %>%
                         mutate(sex = "Female"))

# add peak name
peak_info <- peak_info %>%
  mutate(peak = paste0(Chr, ":", Start, "-", End))

# create genomic ranges
peak_gr <- GRanges(seqnames = peak_info$Chr,
                   ranges = IRanges(start = peak_info$Start, end = peak_info$End),
                   peak = peak_info$peak)

# annotate H3K27ac and PLAC overlap
peak_info <- peak_info %>%
  mutate(H3K27ac = overlapsAny(peak_gr, ferris_h3k27ac_mm39),
         PLAC = overlapsAny(peak_gr, ferris_plac_mm39, ignore.strand = TRUE))

# identify PLAC contacts
plac_hits <- findOverlaps(peak_gr, ferris_plac_mm39, ignore.strand = TRUE)

plac_map_full <- tibble(peak = peak_gr$peak[queryHits(plac_hits)],
                        bin_id = ferris_plac_mm39$bin_id[subjectHits(plac_hits)]) %>%
  distinct() %>%
  left_join(ferris_plac_contacts_mm39 %>%
      dplyr::select(bin_id, GeneName, ENSEMBL_ID) %>%
      distinct(),
      by = "bin_id",
      relationship = "many-to-many")

# summarise PLAC target genes per peak
plac_targets_full <- plac_map_full %>%
  group_by(peak) %>%
  summarise(PLAC_genes = paste(sort(unique(na.omit(GeneName))), collapse = "; "),
            PLAC_ensembl = paste(sort(unique(na.omit(ENSEMBL_ID))), collapse = "; "),
            n_PLAC_targets = n_distinct(ENSEMBL_ID, na.rm = TRUE),
            .groups = "drop")

# add H3K27ac and PLAC 
atac_full <- atac_full %>%
  left_join(peak_info %>%
      dplyr::select(peak, PeakID, Chr, Start, End, H3K27ac, PLAC),
      by = "peak") %>%
  left_join(plac_targets_full, by = "peak") %>%
  mutate(PLAC_genes = na_if(PLAC_genes, ""),
         PLAC_ensembl = na_if(PLAC_ensembl, "")) %>%
  dplyr::select(sex, peak, PeakID, Chr, Start, End,
                log2FoldChange, padj, H3K27ac, PLAC, n_PLAC_targets, PLAC_genes, PLAC_ensembl)

# inspect
atac_full %>%
  group_by(sex) %>%
  summarise(n_peaks = n(),
            n_H3K27ac = sum(H3K27ac, na.rm = TRUE),
            n_PLAC = sum(PLAC, na.rm = TRUE),
            n_PLAC_targets = sum(!is.na(PLAC_ensembl)),
            .groups = "drop")

# combine full male and female RNA results
rna_full <- bind_rows(rna_results$male_basal_vs_60 %>%
                        mutate(sex = "Male"),
                      rna_results$female_basal_vs_60 %>%
                        mutate(sex = "Female")) %>%
  dplyr::select(sex, ensembl, symbol, log2FoldChange, padj, diffexpressed)

# combine male and female ATAC and RNA results together
rna_atac_full <- atac_full %>%
  filter(!is.na(PLAC_ensembl)) %>%
  separate_rows(PLAC_ensembl, sep = ";\\s*") %>%
  mutate(PLAC_ensembl = sub("\\..*$", "", PLAC_ensembl)) %>%
  left_join(rna_full %>%
              dplyr::select(sex, ensembl, symbol,
                            rna_log2FC = log2FoldChange,
                            rna_padj = padj,
                            diffexpressed),
    by = c("sex", "PLAC_ensembl" = "ensembl")) %>%
  dplyr::rename(ensembl = PLAC_ensembl,
                atac_log2FC = log2FoldChange,
                atac_padj = padj)

# extract significant PLAC-linked DAR–gene relationships
plac_dar_genes <- rna_atac_full %>%
  filter(sex %in% c("Male", "Female"),
         !is.na(ensembl),
         !is.na(atac_padj),
         atac_padj < 0.05) %>%
  transmute(sex, peak, ensembl, atac_log2FC, atac_padj, H3K27ac, PLAC) %>%
  distinct(sex, peak, ensembl, .keep_all = TRUE)

# check number of relationships and unique DARs
plac_dar_genes %>%
  group_by(sex) %>%
  summarise(n_DARs = n_distinct(peak),
            n_DAR_gene_pairs = n(),
            n_genes = n_distinct(ensembl),
            .groups = "drop")

# classify PLAC-linked DAR–gene relationships
plac_candidates <- rna_atac_full %>%
  semi_join(plac_dar_genes %>% 
              dplyr::select(sex, peak, ensembl),
            by = c("sex", "peak", "ensembl")) %>%
  mutate(DEG = diffexpressed %in% c("UP", "DOWN"),
         category = if_else(DEG, "A", "C"))

# identify Category B candidates
candidate_B <- rna_full %>%
  filter(diffexpressed %in% c("UP", "DOWN")) %>%
  anti_join(plac_candidates %>%
              filter(category == "A") %>%
              distinct(sex, ensembl),
            by = c("sex", "ensembl")) %>%
  transmute(sex,
            category = "B",
            gene = symbol,
            ensembl,
            rna_log2FC = log2FoldChange,
            rna_padj = padj,
            diffexpressed)

# combine Categories A, B and C
plac_candidates <- bind_rows(plac_candidates,
                             candidate_B %>%
                               mutate(symbol = gene,
                                      DEG = TRUE) %>%
                               dplyr::select(-gene))

# number of unique genes per category
# category A: DEGs with at least 1 significant DAR that has PLAC-seq supported contact of that gene
# category B: DEGs with no signfiicant PLAC-seq supported DAR
# category C: Genes with PLAC-seq supported DAR but no significant DEG
plac_candidates %>%
  group_by(sex, category) %>%
  summarise(n_genes = n_distinct(ensembl),
            n_DARs = n_distinct(peak, na.rm = TRUE),
            n_pairs = sum(!is.na(peak)),
            .groups = "drop")

# extract unique candidate genes and their categories
candidate_genes <- plac_candidates %>%
  distinct(sex, category, ensembl)

# retrieve significant PLAC-linked DARs for A/C and all PLAC-linked peaks for B
candidate_full <- bind_rows(plac_candidates %>%
                              filter(category %in% c("A", "C")) %>%
                              distinct(sex, category, ensembl, peak, H3K27ac),
                            
                            candidate_genes %>%
                              filter(category == "B") %>%
                              left_join(rna_atac_full %>%
                                          distinct(ensembl, peak, H3K27ac),
                                        by = "ensembl", 
                                        relationship = "many-to-many")) %>%
  
  group_by(category, ensembl, peak) %>%
  mutate(sex = if_else(category %in% c("A", "C") & n_distinct(sex) == 2,
                       "Both", sex)) %>%
  distinct(category, ensembl, peak, sex, .keep_all = TRUE) %>%
  ungroup()

# add full RNA-seq and ATAC-seq results for both sexes
candidate_full <- candidate_full %>%
  left_join(rna_full %>%
      dplyr::select(ensembl, sex, log2FoldChange, padj) %>%
      pivot_wider(names_from = sex,
                  values_from = c(log2FoldChange, padj),
                  names_glue = "rna_{.value}_{sex}"),
    by = "ensembl") %>%
  left_join(
    atac_full %>%
      dplyr::select(peak, sex, log2FoldChange, padj) %>%
      pivot_wider(names_from = sex,
                  values_from = c(log2FoldChange, padj),
                  names_glue = "atac_{.value}_{sex}"),
    by = "peak") %>%
  left_join(rna_full %>% distinct(ensembl, symbol),
            by = "ensembl") %>%
  dplyr::select(sex, category, symbol, ensembl, peak, H3K27ac,
                starts_with("rna_"), starts_with("atac_")) %>%
  arrange(category, sex, symbol, peak)

# split into the 3 categories
candidate_A_full <- candidate_full %>% filter(category == "A")
candidate_B_full <- candidate_full %>%
  filter(category == "B") %>%
  group_by(ensembl, peak) %>%
  mutate(sex = if (n_distinct(sex) == 2) "Both" else sex) %>%
  distinct(ensembl, peak, sex, .keep_all = TRUE) %>%
  ungroup()
candidate_C_full <- candidate_full %>% filter(category == "C")

# inspect candidate A gene-peak pairs
view(candidate_A_full %>%
  dplyr::select(sex, symbol, peak, H3K27ac,
                starts_with("rna_"), starts_with("atac_")) %>%
  arrange(symbol, peak, sex) %>%
  print(n = Inf, width = Inf))

# interesting candidates: Sgk1, Dio2, Egr2, Errfi1, Gpr3, Nfkbia, Nr4a1, tSC22D3, Zfp36

# inspect Category B DEGs and their associated accessibility changes
View(candidate_B_full %>%
       group_by(sex, symbol, ensembl) %>%
       summarise(n_peaks = n_distinct(peak, na.rm = TRUE),
                 n_H3K27ac = sum(H3K27ac %in% TRUE),
                 
                 rna_log2FC_Male = dplyr::first(rna_log2FoldChange_Male),
                 rna_log2FC_Female = dplyr::first(rna_log2FoldChange_Female),
                 rna_padj_Male = dplyr::first(rna_padj_Male),
                 rna_padj_Female = dplyr::first(rna_padj_Female),
                 
                 # number of peaks for each DEG with concordant direction
                 n_concordant_Male = sum(rna_log2FoldChange_Male * atac_log2FoldChange_Male > 0,
                                         na.rm = TRUE),
                 n_concordant_Female = sum(rna_log2FoldChange_Female * atac_log2FoldChange_Female > 0,
                                           na.rm = TRUE),
                 
                 mean_atac_log2FC_Male = mean(atac_log2FoldChange_Male, na.rm = TRUE),
                 mean_atac_log2FC_Female = mean(atac_log2FoldChange_Female, na.rm = TRUE),
                 .groups = "drop") %>%
       mutate(pct_concordant_Male = 100 * n_concordant_Male / n_peaks,
              pct_concordant_Female = 100 * n_concordant_Female / n_peaks) %>%
       arrange(sex, symbol))

# interesting candidates: DDIT4, EGR1, FOS, FOSB, PER1, TSC22D3,NFKBIA, CTLA2A

# inspect Category C genes with multiple significant linked DARs
View(candidate_C_full %>%
       group_by(sex, symbol, ensembl) %>%
       summarise(
         n_DARs = n_distinct(peak),
         n_H3K27ac = sum(H3K27ac %in% TRUE),
         rna_log2FC_Male = dplyr::first(rna_log2FoldChange_Male),
         rna_log2FC_Female = dplyr::first(rna_log2FoldChange_Female),
         rna_padj_Male = dplyr::first(rna_padj_Male),
         rna_padj_Female = dplyr::first(rna_padj_Female),
         atac_log2FC_Male = dplyr::first(atac_log2FoldChange_Male),
         atac_log2FC_Female = dplyr::first(atac_log2FoldChange_Female),
         atac_padj_Male = dplyr::first(atac_padj_Male),
         atac_padj_Female = dplyr::first(atac_padj_Female),
         .groups = "drop") %>%
       filter(n_DARs > 1) %>%
       arrange(desc(n_DARs)))

# interesting candidates: SKA1, TEDDMIB, CCDC182, ZBTB16, SGK3, HTR3B, FOXO1, FKBP5
