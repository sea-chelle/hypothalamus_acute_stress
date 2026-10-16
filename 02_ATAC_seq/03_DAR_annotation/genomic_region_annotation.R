# load required packages
library(GenomicRanges)
library(IRanges)
library(GenomeInfoDb)
library(ChIPseeker)
library(org.Mm.eg.db)
library(dplyr)
library(GenomicFeatures)

# keep peaks on canonical chromosomes only
canonical_chr <- c(paste0("chr", 1:19), "chrX", "chrY", "chrM")

# all consensus ATAC peaks as GRanges
consensus_gr <- GRanges(seqnames = peak_info$Chr,
                        ranges = IRanges(start = peak_info$Start,
                                         end = peak_info$End))

consensus_gr <- keepSeqlevels(consensus_gr,
                              canonical_chr[canonical_chr %in% seqlevels(consensus_gr)],
                              pruning.mode = "coarse")

# collapse ChIPseeker categories into simpler labels
collapse_annotation <- function(x) {
  case_when(grepl("^Promoter", x) ~ "Promoter",
            grepl("Intron", x) ~ "Intron",
            grepl("Exon", x) ~ "Exon",
            grepl("UTR", x) ~ "UTR",
            grepl("Downstream", x) ~ "Downstream",
            grepl("Distal Intergenic", x) ~ "Distal intergenic",
            TRUE ~ "Other")
}

# annotate full accessible-peak background
peakAnno_all <- annotatePeak(consensus_gr, 
                             TxDb = txdb_mm39,
                             annoDb = "org.Mm.eg.db",
                             tssRegion = c(-3000, 3000))
peakAnno_all_df <- as.data.frame(peakAnno_all)

# change peakAnno files to new labels
peakAnno_all_df$region <- collapse_annotation(peakAnno_all_df$annotation)
peakAnno_basal_vs_60_df$region <- collapse_annotation(peakAnno_basal_vs_60_df$annotation)
peakAnno_male_basal_vs_60_df$region <- collapse_annotation(peakAnno_male_basal_vs_60_df$annotation)
peakAnno_female_basal_vs_60_df$region <- collapse_annotation(peakAnno_female_basal_vs_60_df$annotation)

# comparisons as a vector
comparisons <- list(consensus = list(anno = peakAnno_all_df,
                                     type = "background"),
                    
                    basal_vs_60 = list(sig = sig_basal_vs_60,
                                       anno = peakAnno_basal_vs_60_df,
                                       type = "DAR"),
                    
                    male_basal_vs_60 = list(sig = sig_male_basal_vs_60,
                                            anno = peakAnno_male_basal_vs_60_df,
                                            type = "DAR"),
                    
                    female_basal_vs_60 = list(sig = sig_female_basal_vs_60,
                                              anno = peakAnno_female_basal_vs_60_df,
                                              type = "DAR"))

# calculate proportion of peaks in each category for each comparison
for (comparison in names(comparisons)) {
  
  comparisons[[comparison]]$region_prop <-
    comparisons[[comparison]]$anno %>%
    dplyr::count(region) %>%
    mutate(proportion = n / sum(n),
           comparison = comparison)
}

# combine results
region_prop <- bind_rows(lapply(comparisons, `[[`, "region_prop"))

# add peak ID to all comparisons
for (comparison in names(comparisons)) {
  comparisons[[comparison]]$anno$peak <- paste0(
    comparisons[[comparison]]$anno$seqnames, ":",
    comparisons[[comparison]]$anno$start, "-",
    comparisons[[comparison]]$anno$end)
}

# get DESeq2 comparisons 
dar_comps <- c("basal_vs_60", "male_basal_vs_60", "female_basal_vs_60")

# grab peak regions
regions <- levels(region_prop$region)

# create empty vector for results
enrichment_results <- list()

# run Fisher exact test to see if DARs are enriched in certain regions compared to conensus background
for (comp in dar_comps) {
  
  dar <- comparisons[[comp]]$anno
  background <- comparisons$consensus$anno %>%
    filter(!peak %in% dar$peak)
  
  enrichment_results[[comp]] <- lapply(regions, function(reg) {
    
    test <- fisher.test(
      matrix(c( sum(dar$region == reg),
                sum(dar$region != reg),
                sum(background$region == reg),
                sum(background$region != reg)), nrow = 2, byrow = TRUE))
    
    data.frame(comparison = comp,
               region = reg,
               odds_ratio = unname(test$estimate),
               conf_low = test$conf.int[1],
               conf_high = test$conf.int[2],
               pvalue = test$p.value)
  }) %>%
    bind_rows() %>%
    mutate(padj = p.adjust(pvalue, method = "BH"))
}

enrichment_results <- bind_rows(enrichment_results)
