# load required packages
library(WGCNA)
library(DESeq2)

# detect outlier genes
gsg <- goodSamplesGenes(t(vsd_mat))
gsg$allOK # if TRUE then no outliers, if FALSE then outliers in genes or samples

# check if outliers exist in samples or genes
table(gsg$goodGenes) # checking genes
table(gsg$goodSamples) # checking samples

# retain genes with at least 10 counts in at least four samples
keep_wgcna <- rowSums(counts(dds) >= 10) >= 4

# see how many genes were filtered out (FALSE) and how many were retained (TRUE)
table(keep_wgcna) # kept 24826 genes, removed 19594

# apply filtering to vsd object
vsd_wgcna <- vsd_mat[keep_wgcna, ]
dim(vsd_wgcna)

# create sample metadata for WGCNA module-trait analysis 
traits <- data.frame(sample = sub("\\..*", "", colnames(vsd_wgcna)),
                     sex = ifelse(grepl("\\.m\\.", colnames(vsd_wgcna)), "m", "f"),
                     condition = sub(".*\\.", "", colnames(vsd_wgcna)))

rownames(traits) <- colnames(vsd_wgcna)

# check traits object
head(traits)
table(traits$sex, traits$condition)

# transpose (genes x sample)
vsd_wgcna <- t(vsd_wgcna)

# choose a set of soft-thresholding powers
powers <- c(seq(4, 10, by = 1), seq(12, 20, by = 2))

# assess scale-free topology across candidate powers
sft <- pickSoftThreshold(vsd_wgcna,
                         powerVector = powers,
                         networkType = "signed",
                         verbose = 5)

# inpsect output
soft_power <- sft$powerEstimate
soft_power
sft$fitIndices

# extract soft-thresholding results
sft.data <- sft$fitIndices

# construct signed co-expression network and identify modules
bwnet <- blockwiseModules(vsd_wgcna,
                          power = soft_power,
                          networkType = "signed",
                          TOMType = "signed",
                          maxBlockSize = 5000,
                          minModuleSize = 30,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

# check output
table(bwnet$colors)
length(unique(bwnet$colors))
dim(bwnet$MEs)

# binarise categorical traits
traits$male <- ifelse(traits$sex == "m", 1, 0)
traits$ARS_10 <- ifelse(traits$condition == "10", 1, 0)
traits$ARS_60 <- ifelse(traits$condition == "60", 1, 0)

# number of samples
nSamples_wgcna <- nrow(vsd_wgcna)

# correlate module eigengenes with traits
module.trait.corr <- cor(bwnet$MEs[, colnames(bwnet$MEs) != "MEgrey"], # exclude grey module as these genes were not assigned to a co-expression module
                         traits[, c("male", "ARS_10", "ARS_60")],
                         use = "p")

# calculate p-values for module-trait correlations
module.trait.corr.pval <- corPvalueStudent(module.trait.corr,
                                           nSamples_wgcna)

# adjust p-values for multiple testing
module.trait.corr.FDR <- matrix(p.adjust(module.trait.corr.pval, method = "BH"),
                                nrow = nrow(module.trait.corr.pval),
                                ncol = ncol(module.trait.corr.pval),
                                dimnames = dimnames(module.trait.corr.pval))
