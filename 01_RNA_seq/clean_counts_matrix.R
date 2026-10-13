# load required packages
library(here)
library(DESeq2)
library(org.Mm.eg.db)

# set working directory
setwd(here("featureCounts"))

# import featureCounts count matrix
cts <- read.delim("combined_gene_counts.txt", comment.char = "#",  check.names = F, header = T)

# remove gene info columns
cts <- cts[, -c(2:6)]

# set gene IDs as row names
rownames(cts) <- cts$Geneid

# remove gene ID column
cts$Geneid <- NULL

# create function to clean up col names - only include sample ID, sex and condition
clean_names <- function(x) {
  x <- gsub(".*/", "", x)                                   # drop directories
  x <- gsub("^\\.+", "", x)                                 # drop leading dots
  x <- gsub("^X\\.", "", x)                                 # drop leading X.
  x <- gsub("Aligned.sortedByCoord.out.bam", "", x)         # remove STAR suffix
  x <- gsub("_S[0-9]+_?$", "", x)                           # remove _SX or _SX_
  x <- gsub("\\.STAR\\.", ".", x)                           # collapse .STAR.
  x <- gsub("-", ".", x)                                    # '-' → '.'
  x
}

# apply function to clean col names
colnames(cts) <- clean_names(colnames(cts))

# remove ensembl ID version
rownames(cts) <- sub("\\..*", "", rownames(cts))
