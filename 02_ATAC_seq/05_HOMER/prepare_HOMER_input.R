# set working directory
setwd(here("HOMER"))

# make consensus background bed file
consensus_bed <- comparisons$consensus$anno %>%
  transmute(chr = seqnames, 
            start = start - 1, 
            end = end)
write.table(consensus_bed, "bed/consensus_background.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# make basal vs 60 min bed file with upregulated DARs
overall_up_bed <- comparisons$basal_vs_60$anno %>%
  filter(log2FoldChange >= 0.5, padj <= 0.05) %>%
  transmute(chr = seqnames,
            start = start - 1,
            end = end)
write.table(overall_up_bed, "./bed/overall_UP_DARs.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# make basal vs 60 min bed file with downregulated DARs
overall_down_bed <- comparisons$basal_vs_60$anno %>%
  filter(log2FoldChange <= -0.5, padj <= 0.05) %>%
  transmute(chr = seqnames,
            start = start - 1,
            end = end)
write.table(overall_down_bed, "./bed/overall_DOWN_DARs.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# make basal vs 60 min bed file with upregulated DARs
male_up_bed <- comparisons$male_basal_vs_60$anno %>%
  filter(log2FoldChange >= 0.5, padj <= 0.05) %>%
  transmute(chr = seqnames,
            start = start - 1,
            end = end)
write.table(male_up_bed, "./bed/male_UP_DARs.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# make basal vs 60 min bed file with upregulated DARs
female_up_bed <- comparisons$female_basal_vs_60$anno %>%
  filter(log2FoldChange >= 0.5, padj <= 0.05) %>%
  transmute(chr = seqnames,
            start = start - 1,
            end = end)
write.table(female_up_bed, "./bed/female_UP_DARs.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# skip downregulated for male and female comparisons as only 1 and 2 DOWN DARs, respectively
