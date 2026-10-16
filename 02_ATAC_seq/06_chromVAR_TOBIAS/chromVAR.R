# load required libraries
library(tidyverse)
library(SummarizedExperiment)
library(GenomicRanges)
library(chromVAR)
library(motifmatchr)
library(JASPAR2024)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm39)
library(limma)

# create peak IDs for canonical consensus peaks
consensus_peak_ids <- paste0(as.character(seqnames(consensus_gr)), ":",
                             start(consensus_gr), "-", end(consensus_gr))

# subset raw featureCounts matrix to canonical consensus peaks
cts_chromvar <- cts[consensus_peak_ids, ]

# create chromVAR input object
atac_se <- SummarizedExperiment(assays = list(counts = as.matrix(cts_chromvar)),
                                rowRanges = consensus_gr,
                                colData = colData(dds))

# add GC-content bias
atac_se <- addGCBias(atac_se, genome = BSgenome.Mmusculus.UCSC.mm39)

# get the path to the JASPAR 2024 transcription factor motif database
jaspar_db <- JASPAR2024::db(JASPAR2024::JASPAR2024())
opts <- list(species = 10090, collection = "CORE")

# retrieve mouse CORE transcription factor position frequency matrices (PFMs) from JASPAR 2024
motifs <- getMatrixSet(x = jaspar_db, opts = opts)

# identify which consensus ATAC-seq peaks contain each JASPAR TF motif
motif_ix <- matchMotifs(motifs, atac_se, genome = BSgenome.Mmusculus.UCSC.mm39)

# calculate chromVAR motif accessibility deviations for each sample
dev <- computeDeviations(object = atac_se, annotations = motif_ix)

# extract standardised motif deviation Z-scores
chromvar_z <- deviationScores(dev)

# create annotation table linking JASPAR motif IDs to TF names
motif_annotation <- tibble(motif_id = names(motifs),
                           tf = vapply(motifs, name, character(1)))

# create metadata for chromVAR samples
chromvar_meta <- as.data.frame(colData(dds)) %>%
  rownames_to_column("sample") %>%
  mutate(group = factor(paste(sex, condition, sep = "_"),
                        levels = c("female_basal", "female_60", "male_basal", "male_60")))

# design matrix for female/male × basal/60 min groups
design_chromvar <- model.matrix(~ 0 + group, data = chromvar_meta)
colnames(design_chromvar) <- levels(chromvar_meta$group)
design_chromvar

# define contrasts
chromvar_contrasts <- makeContrasts(
  # average restraint response across both sexes
  Overall_60vsBasal = ((female_60 - female_basal) + (male_60 - male_basal)) / 2,
  
  # restraint response in males
  Male_60vsBasal = male_60 - male_basal,
  
  # restraint response in females
  Female_60vsBasal = female_60 - female_basal,
  
  # formal test of whether the restraint response differs by sex
  SexInteraction = (male_60 - male_basal) - (female_60 - female_basal),
  levels = design_chromvar)

# fit linear models to chromVAR motif deviation Z-scores
fit_chromvar <- lmFit(chromvar_z, design_chromvar)
fit_chromvar <- contrasts.fit(fit_chromvar, chromvar_contrasts)
fit_chromvar <- eBayes(fit_chromvar)

# extract results tables
chromvar_results <- list(Overall = topTable(fit_chromvar, coef = "Overall_60vsBasal", number = Inf, sort.by = "P"),
                         Male = topTable(fit_chromvar, coef = "Male_60vsBasal", number = Inf, sort.by = "P"),
                         Female = topTable(fit_chromvar, coef = "Female_60vsBasal", number = Inf, sort.by = "P"),
                         Interaction = topTable(fit_chromvar, coef = "SexInteraction", number = Inf, sort.by = "P"))

# rename logFC to delta_z
chromvar_results <- lapply(chromvar_results,
  function(x) {x %>%
      rownames_to_column("motif_id") %>%
      dplyr::rename(delta_z = logFC)
    })

# add TF names
chromvar_results <- lapply(chromvar_results,
  function(x) {x %>%
      left_join(motif_annotation, by = "motif_id")
  })

# check number of significant motifs in each comparison
sapply(chromvar_results, function(x) sum(x$adj.P.Val < 0.05))

# look at significant per comparison 
chromvar_results$Overall %>%
  filter(adj.P.Val < 0.05) %>%
  dplyr::select(motif_id, tf, delta_z, P.Value, adj.P.Val)

chromvar_results$Male %>%
  filter(adj.P.Val < 0.05) %>%
  dplyr::select(motif_id, tf, delta_z, P.Value, adj.P.Val)

chromvar_results$Female %>%
  filter(adj.P.Val < 0.05) %>%
  dplyr::select(motif_id, tf, delta_z, P.Value, adj.P.Val)
