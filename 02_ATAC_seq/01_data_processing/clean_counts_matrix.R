# load required libraries
library(here)

# eet working directory
setwd(here("featureCounts"))

# import count matrix
cts <- read.delim("consensus_peak_count_matrix.tsv", row.names = 1, check.names = FALSE)

# simplify sample IDs
colnames(cts) <- sub("_S[0-9]+_L[0-9]+$", "", colnames(cts))

# create metadata for count matrix
colData <- data.frame(sex = ifelse(grepl("-f-", colnames(cts)), "female", "male"),
                      condition = ifelse(grepl("-basal", colnames(cts)), "basal", "60"))
rownames(colData) <- colnames(cts)

# ensure metadata columns are factors
colData$sex <- factor(colData$sex, levels = c("male", "female"))
colData$condition <- factor(colData$condition, levels = c("basal", "60"))

# import full count matrix with peak location
peak_info <- read.delim("consensus_peak_counts.txt", comment.char = "#", check.names = FALSE)

# keep chromosome location of each peak only
peak_info <- peak_info[, c("Geneid", "Chr", "Start", "End")]
colnames(peak_info) <- c("PeakID", "Chr", "Start", "End")

# replace peak names with chromosome info
rownames(cts) <- paste0(peak_info$Chr, ":",
                        peak_info$Start, "-",
                        peak_info$End)
              full.names = TRUE)
