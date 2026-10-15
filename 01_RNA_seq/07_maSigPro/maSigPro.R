# load required packages
library(here)
library(tidyverse)
library(DESeq2)
library(maSigPro)
library(org.Mm.eg.db)

# set working directory
setwd(here("csv"))

# create maSigPro sample metadata
masig_df <- as.data.frame(colData(vsd))
masig_df <- masig_df[colnames(vsd_mat), , drop = FALSE]

# convert condition to actual time in minutes
masig_df$condition <- as.character(masig_df$condition)
masig_df$condition[masig_df$condition == "basal"] <- "0"
masig_df$condition <- as.numeric(masig_df$condition)

# create experimental design
masig_design <- data.frame(Time = masig_df$condition,
                           Replicate = seq_len(nrow(masig_df)),
                           Male = as.integer(masig_df$sex == "m"),
                           Female = as.integer(masig_df$sex == "f"),
                           row.names = rownames(masig_df))

# build maSigPro design matrix
design <- make.design.matrix(masig_design, degree = 2)

# fit time-course model
fit <- p.vector(vsd_mat, design, Q = 0.05, MT.adjust = "BH")

tstep <- T.fit(fit,  step.method = "backward", alfa = 0.05)

# extract significant trajectory genes
sigs <- get.siggenes(tstep, rsq = 0.5, vars = "groups")


# decide how many clusters are needed
for (k in c(4, 5, 6, 7, 8)) {
  see.genes(sigs$sig.genes$FemalevsMale,
            k = k,
            show.fit = TRUE,
            dis = design$dis,
            edesign = masig_design,
            cluster.method = "hclust",
            cluster.data = 1)
  }

# generate final clustering and save result
see <- see.genes(sigs$sig.genes$FemalevsMale,
                 k = 5,
                 show.fit = TRUE,
                 dis = design$dis,
                 edesign = masig_design,
                 cluster.method = "hclust",
                 cluster.data = 1)

# extract gene names per cluster
cluster_genes <- split(names(see$cut), see$cut)

# check number of genes in each cluster
sapply(cluster_genes, length)

# get all unique Ensembl IDs across clusters
all_ens <- unique(unlist(cluster_genes))

# map to symbol IDs
map <- AnnotationDbi::select(org.Mm.eg.db,
                             keys = all_ens,
                             keytype = "ENSEMBL",
                             columns = "SYMBOL")

# replace ensembl with symbol IDs inside each cluster 
ens2sym <- setNames(map$SYMBOL, map$ENSEMBL)

cluster_genes <- lapply(cluster_genes, function(x) {
  ifelse(!is.na(ens2sym[x]), ens2sym[x], x)
})

# make table with maSigPro results
masig_table <- data.frame(Cluster = as.integer(see$cut),
                          Gene = unname(ens2sym[names(see$cut)]),
                          Ensembl_ID = names(see$cut),
                          P_value = sigs$sig.genes$FemalevsMale$sig.pvalues[names(see$cut), "p-value"],
                          R_squared = sigs$sig.genes$FemalevsMale$sig.pvalues[names(see$cut), "R-squared"]) %>%
  mutate(P_value = format(P_value, scientific = TRUE, digits = 3),
         R_squared = round(R_squared, 3)) %>%
  arrange(Cluster, Gene)

write.csv(masig_table, "maSigPro_cluster_genes.csv", row.names = FALSE)
