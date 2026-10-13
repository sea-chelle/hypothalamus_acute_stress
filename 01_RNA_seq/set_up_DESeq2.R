# load required packages
library(here)
library(DESeq2)
library(AnnotationDbi)
library(org.Mm.eg.db)

# extract column names from counts matrix
sample_ids <- colnames(cts)

# create metadata
colData <- do.call(rbind, strsplit(sample_ids, "\\."))
colnames(colData) <- c("ID", "sex", "condition")
colData <- as.data.frame(colData)
rownames(colData) <- sample_ids

# add batch column
colData$batch <- c(rep("nov_2025", 19),
                   rep("may_2025", 18))

# create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design = ~ batch + sex + condition)

# filter out genes with less than 10 total counts across all samples
keep <- rowSums(counts(dds))  >= 10
dds <- dds[keep,]

# normalise counts
vsd <- vst(dds, blind = T)

# estimate size factors
dds <- estimateSizeFactors(dds)

# create cts file with symbol column
cts_symbol <- cts
cts_symbol$symbol <-  mapIds(org.Mm.eg.db,
                             keys = rownames(cts_symbol),
                             column = "SYMBOL",
                             keytype = "ENSEMBL",
                             multiVals = "first")
