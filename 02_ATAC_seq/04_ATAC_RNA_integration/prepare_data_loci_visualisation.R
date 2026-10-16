# load required packages
library(here)
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(IRanges)
library(rtracklayer)

# ATAC-seq bigWigs (mm39)
atac_bw_files <- list.files(here("bigwig"),
                            pattern = "\\.(bw|bigWig)$",
                            full.names = TRUE,
                            ignore.case = TRUE)

# add metadata
atac_bw_metadata <- tibble(file = atac_bw_files,
                           sample = sub("-.*", "", basename(atac_bw_files),
                           sex = if_else(grepl("-m-", basename(atac_bw_files)), "Male", "Female"),
                           condition = if_else(grepl("-60_", basename(atac_bw_files)), "60 min", "Basal")) %>%
  arrange(factor(sex, levels = c("Male", "Female")),
          factor(condition, levels = c("Basal", "60 min")), sample)

# exclude samples removed during ATAC-seq QC
excluded_samples <- c("2095", "2124", "2127", "2131")

atac_bw_metadata <- atac_bw_metadata %>%
  filter(!sample %in% excluded_samples)

# H3K27ac ChIP-seq bigWigs (mm10)
h3k27ac_bw_files <- list.files(here("ferris_2025", "H3K27ac_ChIP"),
                               pattern = "\\.(bw|bigWig)$",
                               full.names = TRUE,
                               ignore.case = TRUE)

# convert reference h3k27ac peaks to bed file
h3k27ac_bed <- data.frame(chr = as.character(seqnames(ferris_h3k27ac_mm39)),
                          start = start(ferris_h3k27ac_mm39) - 1L,
                          end = end(ferris_h3k27ac_mm39))

write.table(h3k27ac_bed, file.path(here("ferris_2025, "H3k27ac_ChIP", "h3k27ac_reference_mm39.bed")),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# locate PLAC-seq interaction file
plac_file <- list.files(here("ferris_2025", "plac")),
                        pattern = "sig3Dinteractions.*\\.bedpe$",
                        full.names = TRUE)

# check filenames
basename(atac_bw_files)
basename(h3k27ac_bw_files)
basename(plac_file)

# get plac loops
plac_loops <- fread(plac_file)
dim(plac_loops)
colnames(plac_loops)
head(plac_loops)

# convert BEDPE coordinates (0-based starts) to GRanges
anchor1 <- GRanges(seqnames = plac_loops$chr1,
                   ranges = IRanges(plac_loops$start1 + 1, plac_loops$end1))

anchor2 <- GRanges(seqnames = plac_loops$chr2,
                   ranges = IRanges(plac_loops$start2 + 1, plac_loops$end2))

# lift both anchors from mm10 to mm39
anchor1_mm39 <- liftOver(anchor1, mm10_to_mm39)
anchor2_mm39 <- liftOver(anchor2, mm10_to_mm39)

# keep loops which map uniquely only
keep <- lengths(anchor1_mm39) == 1 &
  lengths(anchor2_mm39) == 1
table(keep)

# reconstruct plac loops in mm39
a1 <- unlist(anchor1_mm39[keep])
a2 <- unlist(anchor2_mm39[keep])

plac_loops_mm39 <- plac_loops[keep, ] %>%
  mutate(chr1 = as.character(seqnames(a1)),
         start1 = start(a1),
         end1 = end(a1),
         chr2 = as.character(seqnames(a2)),
         start2 = start(a2),
         end2 = end(a2))

# check
dim(plac_loops_mm39)
colnames(plac_loops_mm39)
head(plac_loops_mm39)

# check converted anchor widths
summary(plac_loops_mm39$end1 - plac_loops_mm39$start1 + 1)
summary(plac_loops_mm39$end2 - plac_loops_mm39$start2 + 1)

# check chromosome consistency
table(plac_loops_mm39$chr1 == plac_loops_mm39$chr2)

# create genomic ranges for both PLAC-seq anchors
plac_anchor1 <- GRanges(plac_loops_mm39$chr1,
                        IRanges(plac_loops_mm39$start1, plac_loops_mm39$end1))

plac_anchor2 <- GRanges(plac_loops_mm39$chr2,
                        IRanges(plac_loops_mm39$start2, plac_loops_mm39$end2))

# convert to bedpe file
plac_bedpe <- plac_loops_mm39 %>%
  mutate(start1 = start1 - 1L,
         start2 = start2 - 1L) %>%
  dplyr::select(chr1, start1, end1, chr2, start2, end2)

write.table(plac_bedpe,
            here("ferris_2025, "PLAC", "plac_interactions_mm39.bedpe"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# convert lifted PLAC distal bins to a data frame
plac_bin_mapping <- tibble(bin_id = mcols(ferris_plac_mm39)$bin_id,
                           chr = as.character(seqnames(ferris_plac_mm39)),
                           start = start(ferris_plac_mm39) - 1L,
                           end = end(ferris_plac_mm39)) %>%
  distinct()

# export mapping
write.csv(plac_bin_mapping, here("ferris_2025", "PLAC", "plac_bin_mapping_mm39.csv"), row.names = FALSE)

# combine RNA-seq, ATAC-seq and PLAC-seq associations
rna_atac_plac_all <- rna_atac_full %>%
  filter(!is.na(ensembl), !is.na(peak)) %>%
  inner_join(plac_map_full, by = c("peak", "ensembl" = "ENSEMBL_ID"), relationship = "many-to-many") %>%
  transmute(sex, symbol, ensembl, peak, chr = Chr, start = Start, end = End, bin_id,
            H3K27ac, atac_log2FC, atac_padj, rna_log2FC, rna_padj, diffexpressed) %>%
  distinct()

write.csv(rna_atac_plac_all, here("csv", "rna_atac_plac_all.csv"), row.names = FALSE)
