# load required packages
library(tidyverse)
library(here)
library(ComplexHeatmap)
library(circlize)
library(grid)

# set working directory
setwd(here("plots"))

# integrate RNA-seq results across comparisons
rna_integrated <- imap_dfr(rna_results, ~ .x %>%
                             dplyr::select(symbol, ensembl, baseMean, log2FoldChange, pvalue, padj) %>%
                             mutate(comparison = dplyr::recode(.y,
                                                               "basal_vs_60" = "Overall",
                                                               "male_basal_vs_60" = "Male",
                                                               "female_basal_vs_60" = "Female"))) %>%
  pivot_wider(id_cols = c(symbol, ensembl),
              names_from = comparison,
              values_from = c(baseMean, log2FoldChange, pvalue, padj),
              names_glue = "rna_{comparison}_{.value}")

# select motifs to plot
motifs_plot <- c("GRE", "PGR", "ERb", "THRb", "FOXP1", "Foxo3", "Foxo1", "HRE", "RORa", 
                 "Atf2", "Atf7", "Jun-AP1", "NFkB-p65", "STAT6", "Klf15", "SF1")

# combine HOMER motif enrichment results across comparisons
homer_selected <- bind_rows(motif_results[c("Overall", "Male", "Female")],
                            .id = "comparison") %>%
  filter(motif_name %in% motifs_plot) %>%
  dplyr::select(motif_name, motif, comparison, fold_enrichment, log2_fold_enrichment, qvalue) %>%
  pivot_wider(id_cols = c(motif_name, motif),
              names_from = comparison,
              values_from = c(fold_enrichment, log2_fold_enrichment, qvalue),
              names_glue = "homer_{comparison}_{.value}")

# select one representative HOMER motif per group
homer_heatmap <- homer_selected %>%
  group_by(motif_name) %>%
  slice_max(order_by = homer_Overall_fold_enrichment, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(motif_name = factor(motif_name, levels = motifs_plot)) %>%
  arrange(motif_name)

# check selected motifs
homer_heatmap %>%
  dplyr::select(motif_name, motif, homer_Overall_fold_enrichment, homer_Male_fold_enrichment, homer_Female_fold_enrichment) %>%
  print(n = Inf)

# standardise motif labels and representative TF genes
homer_heatmap_mapping <- tibble::tribble(~motif_name, ~label,       ~gene,
                                         "GRE",       "GRE",        "Nr3c1",
                                         "PGR",       "PGR",        "Pgr",
                                         "ERb",       "ERβ",        "Esr2",
                                         "THRb",      "THRβ",       "Thrb",
                                         "FOXP1",     "FOXP1",      "Foxp1",
                                         "Foxo3",     "FOXO3",      "Foxo3",
                                         "Foxo1",     "FOXO1",      "Foxo1",
                                         "HRE",       "HSE (HSF)",  "Hsf1",
                                         "RORa",      "RORα",       "Rora",
                                         "Atf2",      "ATF2",       "Atf2",
                                         "Atf7",      "ATF7",       "Atf7",
                                         "Jun-AP1",   "JUN/AP-1",   "Jun",
                                         "NFkB-p65",  "NF-κB p65",  "Rela",
                                         "STAT6",     "STAT6",     "Stat6",
                                         "Klf15",     "KLF15",     "Klf15",
                                         "SF1",       "SF-1",       "Nr5a1")

# join mapping and RNA-seq results
homer_heatmap <- homer_heatmap %>%
  mutate(motif_name = as.character(motif_name)) %>%
  left_join(homer_heatmap_mapping, by = "motif_name") %>%
  left_join(rna_integrated %>%
              filter(!is.na(symbol)),
            by = c("gene" = "symbol"),
            na_matches = "never")

# construct HOMER and RNA-seq matrices for plotting
homer_comparisons <- c("Overall", "Male", "Female")
homer_mat <- as.matrix(homer_heatmap[paste0("homer_", homer_comparisons, "_log2_fold_enrichment")])
rna_mat <- as.matrix(homer_heatmap[paste0("rna_", homer_comparisons, "_log2FoldChange")])

# use motif name + representative TF gene as row labels
motif_plotmath <- c('"GRE"', '"PGR"', '"ER"*beta', '"THR"*beta',
                    '"FOXP1"', '"FOXO3"', '"FOXO1"', '"HSE (HSF)"',
                    '"ROR"*alpha', '"ATF2"', '"ATF7"', '"JUN/AP-1"',
                    '"NF-"*kappa*"B p65"', '"STAT6"', '"KLF15"', '"SF-1"')

row_labels <- parse(text = paste0(motif_plotmath, ' ~ "(" * italic("', homer_heatmap$gene, '") * ")"'))
rownames(homer_mat) <- rownames(rna_mat) <- row_labels
colnames(homer_mat) <- colnames(rna_mat) <- homer_comparisons

# construct HOMER significance matrix (q < 0.05)
homer_sig <- as.matrix(homer_heatmap[paste0("homer_", homer_comparisons, "_qvalue")]) < 0.05

# construct RNA-sew significance matrix (padj < 0.05 AND |log2FC| > 0.5)
rna_sig <- as.matrix(homer_heatmap[paste0("rna_", homer_comparisons, "_padj")]) < 0.05 & abs(rna_mat) > 0.5

# remove any N/A values
homer_sig[is.na(homer_sig)] <- FALSE
rna_sig[is.na(rna_sig)] <- FALSE

# define colour scales
homer_col <- colorRamp2(c(0, max(homer_mat, na.rm = TRUE)), c("white", "#7656A8"))

rna_limit <- max(abs(rna_mat), na.rm = TRUE)
rna_col <- colorRamp2(c(-rna_limit, 0, rna_limit), c("#2674B8", "white", "#D9364A"))

# create HOMER heatmap
ht_homer <- Heatmap(homer_mat,
                    name = "HOMER",
                    col = homer_col,
                    cluster_rows = FALSE,
                    cluster_columns = FALSE,
                    column_title = "HOMER enrichment",
                    column_title_gp = gpar(fontsize = 9, fontface = "bold"),
                    column_names_rot = 0,
                    column_names_centered = TRUE,
                    column_names_gp = gpar(fontsize = 8),
                    row_names_side = "left",
                    row_names_gp = gpar(fontsize = 8),
                    row_labels = row_labels,
                    heatmap_legend_param = list(title = "HOMER\nlog2 enrichment",
                                                title_gp = gpar(fontsize = 7),
                                                labels_gp = gpar(fontsize = 7)),
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      if (homer_sig[i, j]) {
                        grid.text("*", x, y, gp = gpar(fontsize = 10))}})

# create RNA-seq heatmap
ht_rna <- Heatmap(rna_mat,
                  name = "RNA-seq",
                  col = rna_col,
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  column_title = "TF gene expression",
                  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
                  column_names_rot = 0,
                  column_names_centered = TRUE,
                  column_names_gp = gpar(fontsize = 8),
                  show_row_names = FALSE,
                  heatmap_legend_param = list(title = "RNA-seq\nlog2FC",
                                              title_gp = gpar(fontsize = 7),
                                              labels_gp = gpar(fontsize = 7)),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (rna_sig[i, j]) {
      grid.text("*", x, y, gp = gpar(fontsize = 10))}})

# combine and export
pdf("HOMER_RNA_integrated_heatmap.pdf", width = 4.6, height = 4)
draw(ht_homer + ht_rna,
     heatmap_legend_side = "right",
     gap = unit(3, "mm"))
dev.off()
