# load required libraries
library(here)
library(JASPAR2024)
library(TFBSTools)

# set working directory
setwd(here("tobias"))

# write the same JASPAR2024 mouse CORE motifs used for chromVAR
jaspar_file <- "JASPAR2024_mouse_CORE.jaspar"

con <- file(jaspar_file, open = "w")

for (i in seq_along(motifs)) {
  
  motif <- motifs[[i]]
  motif_id <- TFBSTools::ID(motif)
  motif_name <- TFBSTools::name(motif)
  mat <- TFBSTools::Matrix(motif)
  
  writeLines(paste0(">", motif_id, "\t", motif_name), con)
  
  for (base in c("A", "C", "G", "T")) {
    writeLines(paste0(base, " [ ", paste(mat[base, ], collapse = " "), " ]"), con)
  }
}

close(con)
