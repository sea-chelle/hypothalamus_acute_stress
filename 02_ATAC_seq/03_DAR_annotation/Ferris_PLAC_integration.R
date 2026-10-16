# load required packages
library(here)
library(readxl)
library(tidyverse)
library(GenomicRanges)
library(IRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(org.Mm.eg.db)

# set working directory
setwd(here("ferris_2025"))

# import Ferris et al. PLAC-seq data
ferris_plac <- read_excel("PLAC/TableS6_HyrbidPLACseqDistalBinsWithContactedGenes.xlsx", sheet = 2)

# get unique PLAC-seq distal bins
ferris_plac_bins <- ferris_plac %>%
  distinct(DistalChr, DistalStart, DistalEnd) %>%
  mutate(bin_id = paste0(DistalChr, ":", DistalStart, "-", DistalEnd))

dim(ferris_plac_bins)

ferris_plac_gr <- GRanges(seqnames = ferris_plac_bins$DistalChr,
                          ranges = IRanges(start = ferris_plac_bins$DistalStart + 1,
                                           end = ferris_plac_bins$DistalEnd),
                          bin_id = ferris_plac_bins$bin_id)

# lift over to mm39
ferris_plac_mm39 <- liftOver(ferris_plac_gr, mm10_to_mm39)
table(lengths(ferris_plac_mm39))

# check chromosome naming after liftOver
seqlevelsStyle(ferris_plac_mm39)
table(seqnames(ferris_plac_mm39))

# number of original vs uniquely mapped PLAC bins
nrow(ferris_plac_bins)
length(ferris_plac_mm39)

# keep PLAC bins with a unique mm39 mapping
ferris_plac_mm39 <- unlist(ferris_plac_mm39[lengths(ferris_plac_mm39) == 1])
length(ferris_plac_mm39)

# create mm39 coordinate table for retained PLAC bins
ferris_plac_mm39_df <- data.frame(bin_id = ferris_plac_mm39$bin_id,
                                  chr_mm39 = as.character(seqnames(ferris_plac_mm39)),
                                  start_mm39 = start(ferris_plac_mm39),
                                  end_mm39 = end(ferris_plac_mm39))

# add bin_id to original PLAC contact table
ferris_plac <- ferris_plac %>%
  mutate(bin_id = paste0(DistalChr, ":", DistalStart, "-", DistalEnd))

# reconnect each lifted bin to all contacted promoters
ferris_plac_contacts_mm39 <- ferris_plac %>%
  inner_join(ferris_plac_mm39_df, by = "bin_id")

dim(ferris_plac_contacts_mm39)
head(ferris_plac_contacts_mm39)

# check how many unique genes are repsented
n_distinct(ferris_plac_contacts_mm39$GeneName)

# annotate PLAC overlap for all comparison peak sets
for (comparison in names(comparisons)) {
  
  peak_gr <- GRanges(seqnames = comparisons[[comparison]]$anno$seqnames,
    ranges = IRanges(start = comparisons[[comparison]]$anno$start, 
                     end = comparisons[[comparison]]$anno$end),
    peak = comparisons[[comparison]]$anno$peak)
  
  plac_hits <- findOverlaps(peak_gr, ferris_plac_mm39, ignore.strand = TRUE)
  
  # TRUE/FALSE PLAC overlap annotation
  comparisons[[comparison]]$anno$PLAC <- overlapsAny( peak_gr, ferris_plac_mm39, ignore.strand = TRUE)
  
  # only build promoter-target mapping for DAR sets
  if (comparisons[[comparison]]$type == "DAR") {
    
    plac_map <- data.frame(peak = peak_gr$peak[queryHits(plac_hits)],
                           bin_id = ferris_plac_mm39$bin_id[subjectHits(plac_hits)]) %>%
      dplyr::distinct() %>%
      dplyr::left_join(ferris_plac_contacts_mm39 %>%
          dplyr::select(bin_id, GeneName, ENSEMBL_ID) %>%
          dplyr::distinct(),
          by = "bin_id",
          relationship = "many-to-many")
    
    comparisons[[comparison]]$plac_map <- plac_map
  }
}

# proportion of ATAC peaks that have PLAC contacts for each comparison
for (comparison in names(comparisons)) {
  cat("\n", comparison, "\n")
  print(table(comparisons[[comparison]]$anno$PLAC))
  print(prop.table(table(comparisons[[comparison]]$anno$PLAC)))
}

# summarise unique PLAC-linked promoter targets for each DAR
for (comparison in names(comparisons)) {
  
  if (comparisons[[comparison]]$type == "DAR") {
    
    plac_peak_gene <- comparisons[[comparison]]$plac_map %>%
      dplyr::select(peak, GeneName, ENSEMBL_ID) %>%
      dplyr::distinct()
    
    plac_targets_per_peak <- plac_peak_gene %>%
      dplyr::group_by(peak) %>%
      dplyr::summarise(n_PLAC_targets = dplyr::n_distinct(GeneName), .groups = "drop")
    
    comparisons[[comparison]]$plac_peak_gene <- plac_peak_gene
    comparisons[[comparison]]$plac_targets_per_peak <- plac_targets_per_peak
  }
}

# empty vector to store results
plac_enrichment <- list()

# for each comparison run Fishers exact test to determine if ATAC peaks with PLAC contacts are significantly enriched
for (comparison in names(comparisons)) {
  
  if (comparisons[[comparison]]$type == "DAR") {
    
    # DARs
    dar <- comparisons[[comparison]]$anno
    
    # remaining non-DAR consensus peaks
    background <- comparisons$consensus$anno %>%
      filter(!peak %in% dar$peak)
    
    # contingency table
    tab <- matrix(c(sum(dar$PLAC), sum(!dar$PLAC), sum(background$PLAC), sum(!background$PLAC)),
                  nrow = 2, byrow = TRUE,
                  dimnames = list(set = c("DAR", "Background"),
                                  PLAC = c("Yes", "No")))
    
    fisher <- fisher.test(tab)
    
    plac_enrichment[[comparison]] <- data.frame(comparison = comparison,
                                                odds_ratio = unname(fisher$estimate),
                                                CI_low = fisher$conf.int[1],
                                                CI_high = fisher$conf.int[2],
                                                pvalue = fisher$p.value)
  }
}

plac_enrichment <- bind_rows(plac_enrichment) %>% 
  mutate(padj = p.adjust(pvalue, method = "BH"))

plac_enrichment

# compare ChIPseeker nearest genes with PLAC-linked promoter targets
for (comparison in names(comparisons)) {
  
  if (comparisons[[comparison]]$type == "DAR") {
    
    # map nearest-gene Entrez IDs to Ensembl IDs
    entrez_to_ensembl <- AnnotationDbi::select(org.Mm.eg.db,
                                               keys = unique(na.omit(comparisons[[comparison]]$anno$ENTREZID)),
                                               keytype = "ENTREZID",
                                               columns = "ENSEMBL")
    
    # collapse multiple Ensembl IDs per Entrez ID
    entrez_to_ensembl_collapsed <- entrez_to_ensembl %>%
      group_by(ENTREZID) %>%
      summarise(ENSEMBL = paste(sort(unique(na.omit(ENSEMBL))), collapse = ";"),
                .groups = "drop")
    
    # add Ensembl ID to annotation table
    comparisons[[comparison]]$anno <- comparisons[[comparison]]$anno %>%
      dplyr::select(-dplyr::any_of("ENSEMBL")) %>%
      left_join(entrez_to_ensembl_collapsed, by = "ENTREZID")
    
    # retain peak-to-Ensembl mappings for nearest-vs-PLAC comparison
    nearest_ensembl <- comparisons[[comparison]]$anno %>%
      dplyr::select(peak, ENSEMBL) %>%
      separate_rows(ENSEMBL, sep = ";") %>%
      rename(nearest_ensembl = ENSEMBL)
    
    # compare nearest-gene IDs with PLAC-linked promoter targets
    nearest_vs_plac <- comparisons[[comparison]]$plac_peak_gene %>%
      left_join(nearest_ensembl, by = "peak", relationship = "many-to-many") %>%
      mutate(nearest_gene_is_PLAC_target = nearest_ensembl == ENSEMBL_ID)
    
    # collapse to one result per DAR
    nearest_vs_plac_summary <- nearest_vs_plac %>%
      group_by(peak) %>%
      summarise(nearest_gene_supported = any(nearest_gene_is_PLAC_target, na.rm = TRUE),
                .groups = "drop")
    
    comparisons[[comparison]]$nearest_vs_plac <- nearest_vs_plac
    comparisons[[comparison]]$nearest_vs_plac_summary <- nearest_vs_plac_summary
  }
}

# number and proportion of nearest genes supported by PLAC-seq
for (comparison in names(comparisons)) {
  
  if (comparisons[[comparison]]$type == "DAR") {
    
    cat("\n", comparison, "\n")
    print(table(comparisons[[comparison]]$nearest_vs_plac_summary$nearest_gene_supported))
    print(prop.table(table(comparisons[[comparison]]$nearest_vs_plac_summary$nearest_gene_supported)))
  }
}

# add PLAC summary information to DAR annotation tables
for (comparison in names(comparisons)) {
  
  if (comparisons[[comparison]]$type == "DAR") {
    
    comparisons[[comparison]]$anno <- comparisons[[comparison]]$anno %>%
      left_join(comparisons[[comparison]]$plac_targets_per_peak, by = "peak") %>%
      left_join(comparisons[[comparison]]$nearest_vs_plac_summary, by = "peak") %>%
      mutate(n_PLAC_targets = dplyr::coalesce(n_PLAC_targets, 0L))
  }
}
