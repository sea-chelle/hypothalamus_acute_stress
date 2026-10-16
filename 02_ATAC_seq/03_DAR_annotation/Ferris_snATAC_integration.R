# load required packages
library(here)
library(tidyverse)
library(GenomicRanges)
library(IRanges)
library(rtracklayer)

# set working directory
setwd(here("ferris_2025"))

# import fed cell-type scATAC promoter-assigned files
snATAC_files <- list.files("snATAC/", pattern = "^Fed_.*_atac_promo_assigned_hibar_hibdel.*\\.txt$", full.names = TRUE)
ferris_snATAC <- lapply(snATAC_files, read.delim)

# add cell type name to each file
names(ferris_snATAC) <- sub("^Fed_(.*?)_atac.*", "\\1", basename(snATAC_files))

# cell types
names(ferris_snATAC)

# number of cells
sapply(ferris_snATAC, nrow)

# retain unique scATAC peaks for each cell type
ferris_snATAC_peaks <- lapply(ferris_snATAC, \(x) x %>%
                                distinct(chrom, chromStart, chromEnd, name))
sapply(ferris_snATAC_peaks, nrow)

# convert unique scATAC peaks to GRanges
ferris_snATAC_gr <- lapply(ferris_snATAC_peaks,
                           \(x) GRanges(seqnames = x$chrom,
                                        ranges = IRanges(start = x$chromStart + 1,
                                                         end = x$chromEnd),
                                        peak = x$name))

# lift over mm10 to mm39 and retain 1:1 mappings
ferris_snATAC_mm39 <- lapply(ferris_snATAC_gr, \(x) {
  lifted <- liftOver(x, mm10_to_mm39)
  lifted <- lifted[lengths(lifted) == 1]
  unlist(lifted)
  })

# check how many peaks have been retained
sapply(ferris_snATAC_mm39, length)

# annotate peaks with Ferris fed hypothalamic cell-type CREs
for (comparison in names(comparisons)) {
  
  peak_gr <- GRanges(seqnames = comparisons[[comparison]]$anno$seqnames,
                     ranges = IRanges(start = comparisons[[comparison]]$anno$start,
                                      end = comparisons[[comparison]]$anno$end))
  
  for (celltype in names(ferris_snATAC_mm39)) {
    
    comparisons[[comparison]]$anno[[paste0("CRE_", celltype)]] <- overlapsAny(peak_gr,
                                                                              ferris_snATAC_mm39[[celltype]],
                                                                              ignore.strand = TRUE)
  }
}

# proportion of peaks overlapping each cell-type CRE set
cre_overlap <- lapply(names(comparisons), function(comp) {
  
  cre_cols <- paste0("CRE_", names(ferris_snATAC_mm39))
  comparisons[[comp]]$anno %>%
    summarise(across(all_of(cre_cols), ~ mean(.x) * 100)) %>%
    mutate(comparison = comp)
  }) %>%
  bind_rows() %>%
  pivot_longer(cols = starts_with("CRE_"), names_to = "celltype", values_to = "percent") %>%
  mutate(celltype = sub("^CRE_", "", celltype),
    comparison = factor(comparison, 
                        levels = c("consensus", "basal_vs_60", "male_basal_vs_60", "female_basal_vs_60"),
      labels = c("Consensus", "Overall", "Male", "Female")))
cre_overlap

# test enrichment of DARs in cell-type CREs relative to non-DAR consensus peaks
cre_enrichment <- list()

for (comparison in names(comparisons)) {
  
  if (comparisons[[comparison]]$type == "DAR") {
    dar_peaks <- comparisons[[comparison]]$anno$peak
    background <- comparisons$consensus$anno %>%
      filter(!peak %in% dar_peaks)
    
    for (celltype in names(ferris_snATAC_mm39)) {
      cre_col <- paste0("CRE_", celltype)
      tab <- matrix(c(sum(comparisons[[comparison]]$anno[[cre_col]]),
                      sum(!comparisons[[comparison]]$anno[[cre_col]]),
                      sum(background[[cre_col]]),
                      sum(!background[[cre_col]])),
                    nrow = 2,
                    byrow = TRUE)
      test <- fisher.test(tab)
      
      cre_enrichment[[paste(comparison, celltype, sep = "_")]] <-
        data.frame(comparison = comparison,
                   celltype = celltype,
                   odds_ratio = unname(test$estimate),
                   CI_low = test$conf.int[1],
                   CI_high = test$conf.int[2],
                   pvalue = test$p.value)
    }
  }
}

# combine results and correct for multiple testing within each comparison
cre_enrichment <- bind_rows(cre_enrichment) %>%
  group_by(comparison) %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  ungroup()
