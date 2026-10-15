# load required packages
library(GenomicFeatures)
library(GenomicRanges)
library(IRanges)
library(ChIPseeker)
library(org.Mm.eg.db)

# build TxDb from mm39 gtf file
txdb_mm39 <- makeTxDbFromGFF("/Volumes/ritd-ag-project-rd0227-nhhan20/Michelle/genomes/mm39/gencode.vM36.basic.annotation.gtf")

# get coordinates of each peak from row names
coords <- strcapture("([^:]+):(\\d+)-(\\d+)",
                     rownames(sig_basal_vs_60),
                     proto = list(Chr = character(),
                                  Start = integer(),
                                  End = integer()))

sig_basal_vs_60_gr <- GRanges(seqnames = coords$Chr,
                              ranges = IRanges(start = coords$Start,
                                               end = coords$End))

# peak annotation with ChIPseekerpeak
peakAnno_basal_vs_60 <- annotatePeak(sig_basal_vs_60_gr,
                                     TxDb = txdb_mm39,
                                     annoDb = "org.Mm.eg.db",
                                     tssRegion = c(-3000, 3000))
peakAnno_basal_vs_60_df <- as.data.frame(peakAnno_basal_vs_60)

# add differential accessibility information
peakAnno_basal_vs_60_df$log2FoldChange <- sig_basal_vs_60$log2FoldChange
peakAnno_basal_vs_60_df$padj <- sig_basal_vs_60$padj
peakAnno_basal_vs_60_df$diffaccessible <- sig_basal_vs_60$diffaccessible
peakAnno_basal_vs_60_df$baseMean <- sig_basal_vs_60$baseMean
peakAnno_basal_vs_60_df$pvalue <- sig_basal_vs_60$pvalue

# create coordinate ID in annotated peaks
peakAnno_basal_vs_60_df$peak <- paste0(peakAnno_basal_vs_60_df$seqnames, ":",
                                       peakAnno_basal_vs_60_df$start, "-",
                                       peakAnno_basal_vs_60_df$end)

# get coordinates of each peak from row names
coords <- strcapture("([^:]+):(\\d+)-(\\d+)",
                     rownames(sig_basal_male_vs_female),
                     proto = list(Chr = character(),
                                  Start = integer(),
                                  End = integer()))

sig_basal_male_vs_female_gr <- GRanges(seqnames = coords$Chr,
                              ranges = IRanges(start = coords$Start,
                                               end = coords$End))

# only include those in TxDb
sig_basal_male_vs_female_gr <- keepSeqlevels(sig_basal_male_vs_female_gr,
  intersect(seqlevels(sig_basal_male_vs_female_gr), seqlevels(txdb_mm39)),
  pruning.mode = "coarse")

# peak annotation with ChIPseekerpeak
peakAnno_basal_male_vs_female <- annotatePeak(sig_basal_male_vs_female_gr,
                                     TxDb = txdb_mm39,
                                     annoDb = "org.Mm.eg.db",
                                     tssRegion = c(-3000, 3000))
peakAnno_basal_male_vs_female_df <- as.data.frame(peakAnno_basal_male_vs_female)

# subset sig DAR file to only include annotated peaks
sig_basal_male_vs_female <- sig_basal_male_vs_female[grepl("^chr", rownames(sig_basal_male_vs_female)), ]

# add differential accessibility information
peakAnno_basal_male_vs_female_df$log2FoldChange <- sig_basal_male_vs_female$log2FoldChange
peakAnno_basal_male_vs_female_df$padj <- sig_basal_male_vs_female$padj
peakAnno_basal_male_vs_female_df$diffaccessible <- sig_basal_male_vs_female$diffaccessible
peakAnno_basal_male_vs_female_df$baseMean <- sig_basal_male_vs_female$baseMean
peakAnno_basal_male_vs_female_df$pvalue <- sig_basal_male_vs_female$pvalue

# create coordinate ID in annotated peaks
peakAnno_basal_male_vs_female_df$peak <- paste0(peakAnno_basal_male_vs_female_df$seqnames, ":",
                                       peakAnno_basal_male_vs_female_df$start, "-",
                                       peakAnno_basal_male_vs_female_df$end)

# get coordinates of each peak from row names
coords <- strcapture("([^:]+):(\\d+)-(\\d+)",
                     rownames(sig_60min_male_vs_female),
                     proto = list(Chr = character(),
                                  Start = integer(),
                                  End = integer()))

sig_60min_male_vs_female_gr <- GRanges(seqnames = coords$Chr,
                                       ranges = IRanges(start = coords$Start,
                                                        end = coords$End))

# peak annotation with ChIPseekerpeak
peakAnno_60min_male_vs_female <- annotatePeak(sig_60min_male_vs_female_gr,
                                              TxDb = txdb_mm39,
                                              annoDb = "org.Mm.eg.db",
                                              tssRegion = c(-3000, 3000))
peakAnno_60min_male_vs_female_df <- as.data.frame(peakAnno_60min_male_vs_female)

# add differential accessibility information
peakAnno_60min_male_vs_female_df$log2FoldChange <- sig_60min_male_vs_female$log2FoldChange
peakAnno_60min_male_vs_female_df$padj <- sig_60min_male_vs_female$padj
peakAnno_60min_male_vs_female_df$diffaccessible <- sig_60min_male_vs_female$diffaccessible
peakAnno_60min_male_vs_female_df$baseMean <- sig_60min_male_vs_female$baseMean
peakAnno_60min_male_vs_female_df$pvalue <- sig_60min_male_vs_female$pvalue

# create coordinate ID in annotated peaks
peakAnno_60min_male_vs_female_df$peak <- paste0(peakAnno_60min_male_vs_female_df$seqnames, ":",
                                                peakAnno_60min_male_vs_female_df$start, "-",
                                                peakAnno_60min_male_vs_female_df$end)

# get coordinates of each peak from row names
coords <- strcapture("([^:]+):(\\d+)-(\\d+)",
                     rownames(sig_male_basal_vs_60),
                     proto = list(Chr = character(),
                                  Start = integer(),
                                  End = integer()))

sig_male_basal_vs_60_gr <- GRanges(seqnames = coords$Chr,
                                       ranges = IRanges(start = coords$Start,
                                                        end = coords$End))

# only include those in TxDb
sig_male_basal_vs_60_gr <- keepSeqlevels(sig_male_basal_vs_60_gr,
                                             intersect(seqlevels(sig_male_basal_vs_60_gr), seqlevels(txdb_mm39)),
                                             pruning.mode = "coarse")

# peak annotation with ChIPseekerpeak
peakAnno_male_basal_vs_60 <- annotatePeak(sig_male_basal_vs_60_gr,
                                              TxDb = txdb_mm39,
                                              annoDb = "org.Mm.eg.db",
                                              tssRegion = c(-3000, 3000))
peakAnno_male_basal_vs_60_df <- as.data.frame(peakAnno_male_basal_vs_60)

# add differential accessibility information
peakAnno_male_basal_vs_60_df$log2FoldChange <- sig_male_basal_vs_60$log2FoldChange
peakAnno_male_basal_vs_60_df$padj <- sig_male_basal_vs_60$padj
peakAnno_male_basal_vs_60_df$diffaccessible <- sig_male_basal_vs_60$diffaccessible
peakAnno_male_basal_vs_60_df$baseMean <- sig_male_basal_vs_60$baseMean
peakAnno_male_basal_vs_60_df$pvalue <- sig_male_basal_vs_60$pvalue

# create coordinate ID in annotated peaks
peakAnno_male_basal_vs_60_df$peak <- paste0(peakAnno_male_basal_vs_60_df$seqnames, ":",
                                                peakAnno_male_basal_vs_60_df$start, "-",
                                                peakAnno_male_basal_vs_60_df$end)

# get coordinates of each peak from row names
coords <- strcapture("([^:]+):(\\d+)-(\\d+)",
                     rownames(sig_female_basal_vs_60),
                     proto = list(Chr = character(),
                                  Start = integer(),
                                  End = integer()))

sig_female_basal_vs_60_gr <- GRanges(seqnames = coords$Chr,
                                   ranges = IRanges(start = coords$Start,
                                                    end = coords$End))

# peak annotation with ChIPseekerpeak
peakAnno_female_basal_vs_60 <- annotatePeak(sig_female_basal_vs_60_gr,
                                          TxDb = txdb_mm39,
                                          annoDb = "org.Mm.eg.db",
                                          tssRegion = c(-3000, 3000))
peakAnno_female_basal_vs_60_df <- as.data.frame(peakAnno_female_basal_vs_60)

# add differential accessibility information
peakAnno_female_basal_vs_60_df$log2FoldChange <- sig_female_basal_vs_60$log2FoldChange
peakAnno_female_basal_vs_60_df$padj <- sig_female_basal_vs_60$padj
peakAnno_female_basal_vs_60_df$diffaccessible <- sig_female_basal_vs_60$diffaccessible
peakAnno_female_basal_vs_60_df$baseMean <- sig_female_basal_vs_60$baseMean
peakAnno_female_basal_vs_60_df$pvalue <- sig_female_basal_vs_60$pvalue

# create coordinate ID in annotated peaks
peakAnno_female_basal_vs_60_df$peak <- paste0(peakAnno_female_basal_vs_60_df$seqnames, ":",
                                            peakAnno_female_basal_vs_60_df$start, "-",
                                            peakAnno_female_basal_vs_60_df$end)
