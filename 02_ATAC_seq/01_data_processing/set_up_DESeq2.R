# load required packages
library(DESeq2)

# create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design = ~ sex * condition)

# filter out peaks with less than 10 total reads across all samples
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep,]

# estimate size factors
dds <- estimateSizeFactors(dds)

# normalise counts by variance-stablising transformation
vsd <- vst(dds, blind = T)
