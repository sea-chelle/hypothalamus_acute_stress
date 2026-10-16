# load required libraries
library(here)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)

# set working directory
setwd(here("plots"))

# integrate chromVAR results across comparisons
chromvar_integrated <-imap_dfr(chromvar_results, ~ .x %>%
                                 dplyr::select(motif_id, tf, delta_z, P.Value, adj.P.Val) %>%
                                 mutate(comparison = .y)) %>%
  pivot_wider(id_cols = c(motif_id, tf),
              names_from = comparison,
              values_from = c(delta_z, P.Value, adj.P.Val),
              names_glue = "chromvar_{comparison}_{.value}")

# integrate chromVAR and TOBIAS using exact JASPAR motif IDs
tf_integrated <- chromvar_integrated %>%
  full_join(tobias_results, by = "motif_id") %>%
  mutate(tf = coalesce(tf, name),
         chromvar_significant = coalesce(chromvar_Male_adj.P.Val < 0.05, FALSE) |
           coalesce(chromvar_Female_adj.P.Val < 0.05, FALSE),
         tobias_highlighted = coalesce(male_highlighted, FALSE) |
           coalesce(female_highlighted, FALSE))

# retain JASPAR motifs identified by either analysis
tf_candidates <- tf_integrated %>%
  filter(chromvar_significant | tobias_highlighted)

# check candidate motifs
tf_candidates %>%
  dplyr::select(motif_id, tf, chromvar_significant, tobias_highlighted) %>%
  arrange(tf) %>%
  print(n = Inf)

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

# map JASPAR motifs to their corresponding TF-encoding genes heterodimer motifs retain both component genes
tf_gene_mapping <- tf_candidates %>%
  dplyr::select(motif_id, tf) %>%
  mutate(gene = if_else(tf == "Hnf1A", "Hnf1a", tf)) %>%
  separate_longer_delim(gene, delim = "::")

# join RNA-seq expression using case-insensitive gene symbols
tf_final <- tf_candidates %>%
  left_join(tf_gene_mapping, by = c("motif_id", "tf")) %>%
  mutate(gene_key = tolower(gene)) %>%
  left_join(rna_integrated %>%
              filter(!is.na(symbol)) %>%
              mutate(gene_key = tolower(symbol)),
            by = "gene_key", na_matches = "never") %>%
  dplyr::select(-gene_key)

# create descriptive labels for the heatmap
tf_final <- tf_final %>%
  mutate(label = if_else(grepl("::", tf), paste0(tf, " (", gene, ")"), tf))

# check integrated results
tf_final %>%
  dplyr::select(motif_id, tf, gene, label,
                chromvar_significant, tobias_highlighted, rna_Male_log2FoldChange, rna_Female_log2FoldChange) %>% 
  arrange(tf, gene) %>%
  print(n = Inf)

# check dimensions and missing RNA-seq mappings
tf_final %>%
  summarise(n_JASPAR_motifs = n_distinct(motif_id),
            n_genes = n_distinct(gene),
            n_rows = n(),
            n_RNA_matched = sum(!is.na(symbol)))

tf_final %>%
  filter(is.na(symbol)) %>%
  dplyr::select(motif_id, tf, gene)

# retain the existing candidate order
tf_heatmap <- tf_final %>%
  mutate(label = make.unique(label))

# create matrices for each analysis
chromvar_mat <- as.matrix(tf_heatmap %>%
                            dplyr::select(chromvar_Male_delta_z, chromvar_Female_delta_z))

tobias_mat <- -as.matrix(tf_heatmap %>%
                           dplyr::select(male_change, female_change))

rownames(tobias_mat) <- tf_heatmap$label
colnames(tobias_mat) <- c("Male", "Female")

rna_mat <- as.matrix(tf_heatmap %>%
                       dplyr::select(rna_Male_log2FoldChange, rna_Female_log2FoldChange))

# assign matching row and column names
for (mat_name in c("chromvar_mat", "tobias_mat", "rna_mat")) {
  mat <- get(mat_name)
  rownames(mat) <- tf_heatmap$label
  colnames(mat) <- c("Male", "Female")
  assign(mat_name, mat)
}

# chromVAR: adjusted p < 0.05
chromvar_sig <- cbind(Male = tf_heatmap$chromvar_Male_adj.P.Val < 0.05,
                      Female = tf_heatmap$chromvar_Female_adj.P.Val < 0.05)

# TOBIAS: highlighted motifs
tobias_sig <- cbind(Male = tf_heatmap$male_highlighted,
                    Female = tf_heatmap$female_highlighted)

# RNA-seq: adjusted p < 0.05 and |log2FC| > 0.5
rna_sig <- cbind(Male = tf_heatmap$rna_Male_padj < 0.05 & abs(tf_heatmap$rna_Male_log2FoldChange) > 0.5,
                 Female = tf_heatmap$rna_Female_padj < 0.05 & abs(tf_heatmap$rna_Female_log2FoldChange) > 0.5)

# treat missing values as non-significant for plotting
chromvar_sig[is.na(chromvar_sig)] <- FALSE
tobias_sig[is.na(tobias_sig)] <- FALSE
rna_sig[is.na(rna_sig)] <- FALSE

# use robust limits to avoid extreme values dominating the colours
colour_limit <- function(mat) {
  max(as.numeric(quantile(abs(mat), 0.95, na.rm = TRUE)), 0.01)
}

chromvar_limit <- colour_limit(chromvar_mat)
tobias_limit <- colour_limit(tobias_mat)
rna_limit <- colour_limit(rna_mat)

# set colour scales
chromvar_col <- colorRamp2(c(-chromvar_limit, 0, chromvar_limit),
                           c("#2674B8", "white", "#D9364A"))

tobias_col <- colorRamp2(c(-tobias_limit, 0, tobias_limit),
                         c("#2674B8", "white", "#D9364A"))

rna_col <- colorRamp2(c(-rna_limit, 0, rna_limit),
                      c("#2674B8", "white", "#D9364A"))

# add asterisks to selected cells
add_stars <- function(sig) {
  function(j, i, x, y, width, height, fill) {
    if (isTRUE(sig[i, j])) {
      grid.text("*", x, y, gp = gpar(fontsize = 10))
    }
  }
}

# create chromVAR heatmap
ht_chromvar <- Heatmap(chromvar_mat,
                       name = chromVAR,
                       col = chromvar_col,
                       na_col = "grey90",
                       cluster_rows = FALSE,
                       cluster_columns = FALSE,
                       show_row_names = TRUE,
                       row_names_side = "left",
                       row_names_gp = gpar(fontsize = 8),
                       column_names_gp = gpar(fontsize = 8),
                       column_title_gp = gpar(fontsize = 9, fontface = "bold"),
                       column_names_rot = 90,
                       column_names_centered = TRUE,
                       column_title = "chromVAR",
                       cell_fun = add_stars(chromvar_sig),
                       width = unit(13, "mm"),
                       heatmap_legend_param = list(title = expression(atop("", Delta * "Z-score")),
                                                   direction = "horizontal",
                                                   title_position = "topcenter",
                                                   legend_width = unit(10, "mm"),
                                                   title_gp = gpar(fontsize = 7),
                                                   labels_gp = gpar(fontsize = 7),
                                                   at = c(-20, 0, 20)))
# create TOBIAS heatmap
ht_tobias <- Heatmap(tobias_mat,
                     name = "TOBIAS\nBinding\nchange",
                     col = tobias_col,
                     na_col = "grey90",
                     cluster_rows = FALSE,
                     cluster_columns = FALSE,
                     show_row_names = FALSE,
                     column_names_gp = gpar(fontsize = 8),
                     column_names_rot = 90,
                     column_names_centered = TRUE,
                     column_title = "TOBIAS",
                     column_title_gp = gpar(fontsize = 9, fontface = "bold"),
                     cell_fun = add_stars(tobias_sig),
                     width = unit(13, "mm"),
                     heatmap_legend_param = list(title = "Binding\nchange",
                                                 direction = "horizontal",
                                                 title_position = "topcenter",
                                                 legend_width = unit(10, "mm"),
                                                 title_gp = gpar(fontsize = 7),
                                                 labels_gp = gpar(fontsize = 7),
                                                 at = c(-0.4, 0, 0.4)))

# create RNA-seq heatmap
ht_rna <- Heatmap(rna_mat,
                  name = "RNA-seq",
                  col = rna_col,
                  na_col = "grey90",
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  show_row_names = FALSE,
                  column_names_gp = gpar(fontsize = 8),
                  column_names_rot = 90,
                  column_names_centered = TRUE,
                  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
                  column_title = "RNA-seq",
                  cell_fun = add_stars(rna_sig),
                  width = unit(13, "mm"),
                  heatmap_legend_param = list(title = "\nlog2FC",
                                              direction = "horizontal",
                                              title_position = "topcenter",
                                              legend_width = unit(10, "mm"),
                                              title_gp = gpar(fontsize = 7),
                                              labels_gp = gpar(fontsize = 7),
                                              at = c(-1, 0, 1)))

# combine heatmaps and export
ht <- ht_chromvar + ht_tobias + ht_rna

pdf("JASPAR_RNA_integrated_heatmap.pdf", width = 2.6, height = 6)
draw(ht,
     heatmap_legend_side = "bottom",
     show_heatmap_legend = TRUE,
     show_annotation_legend = FALSE)
dev.off()
