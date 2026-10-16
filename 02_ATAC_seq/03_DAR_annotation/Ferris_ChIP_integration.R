# load required packages
library(here)
library(rtracklayer)
library(GenomicRanges)
library(IRanges)
library(GenomeInfoDb)
library(dplyr)

# set working directory 
setwd(here("ferris_2025"))

# import ferris et al. H3K27ac consensus peak bed file
ferris_h3k27ac <- import("H3K27ac_ChIP/h3k27ac_consensus_2of3.bed", format = "BED")

# add chr prefix
seqlevels(ferris_h3k27ac) <- paste0("chr", seqlevels(ferris_h3k27ac))

# mm10 -> mm39 command line download: wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz
# unzip: gunzip -k mm10ToMm39.over.chain.gz
# move to H3K27ac_ChIP folder

# load the chain file for mm10 -> mm39
mm10_to_mm39 <- import.chain("H3K27ac_ChIP/mm10ToMm39.over.chain")

# lift over ferris h3k27ac conensus peak to mm39
ferris_h3k27ac_mm39 <- liftOver(ferris_h3k27ac, mm10_to_mm39)

# retain H3K27ac peaks with a unique (1:1) mm39 mapping
ferris_h3k27ac_mm39 <- ferris_h3k27ac_mm39[lengths(ferris_h3k27ac_mm39) == 1]

# convert GRangesList to GRanges
ferris_h3k27ac_mm39 <- unlist(ferris_h3k27ac_mm39)

# check how many peaks there are after liftOver
length(ferris_h3k27ac_mm39)

# overlap ferris h3k27ac peaks with atac-seq peaks
for (comparison in names(comparisons)) {
  
  peak_gr <- GRanges(seqnames = comparisons[[comparison]]$anno$seqnames,
                     ranges = IRanges(start = comparisons[[comparison]]$anno$start,
                                      end = comparisons[[comparison]]$anno$end))
  
  comparisons[[comparison]]$anno$H3K27ac <- overlapsAny(peak_gr, ferris_h3k27ac_mm39)
}

# number of peaks that overlap with h3k27ac peaks for each comparison
for (comparison in names(comparisons)) {
  cat("\n", comparison, "\n")
  print(table(comparisons[[comparison]]$anno$H3K27ac))
}

# percentage of peaks that overlap with h3k27ac peaks for each comparison
for (comparison in names(comparisons)) {
  cat(comparison, round(mean(comparisons[[comparison]]$anno$H3K27ac) * 100, 1),
      "% H3K27ac+\n")
}

# use fishers exact test to determine whether stress-responsive DARs are significantly enriched for h3k27ac peaks
h3k27ac_enrichment <- lapply(c("basal_vs_60", "male_basal_vs_60", "female_basal_vs_60"),
  function(comparison) {
    
    dar <- comparisons[[comparison]]$anno
    
    background <- comparisons$consensus$anno %>%
      dplyr::filter(!peak %in% dar$peak)
    
    test <- fisher.test(matrix(c(sum(dar$H3K27ac), sum(!dar$H3K27ac), sum(background$H3K27ac), sum(!background$H3K27ac)), nrow = 2, byrow = TRUE))
    
    data.frame(comparison = comparison,
               odds_ratio = unname(test$estimate),
               conf_low = test$conf.int[1],
               conf_high = test$conf.int[2],
               pvalue = test$p.value)
  }
  ) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(padj = p.adjust(pvalue, method = "BH"))

h3k27ac_enrichment
