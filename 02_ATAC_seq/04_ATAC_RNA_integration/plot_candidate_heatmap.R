
# candidates to plot
heatmap_genes <- tibble(symbol = c("Sgk1", "Dio2", "Gpr3", "Egr2", "Nr4a1", "Nfkbia", "Tsc22d3", "Zfp36",
                                   "Ddit4", "Egr1", "Fos", "Fosb", "Fosl2", "Per1", "Ctla2a", 
                                   "Ska1", "Teddm1b", "Ccdc182", "Zbtb16", "Sgk3", "Htr3a", "Htr3b", "Foxo1", "Fkbp5")

# extract RNA log2FC, adjusted p-values and DEG classifications
heatmap_rna <- rna_full %>%
  semi_join(heatmap_genes, by = "symbol") %>%
  dplyr::select(symbol, sex, log2FoldChange, padj, diffexpressed) %>%
  pivot_wider(names_from = sex,
              values_from = c(log2FoldChange, padj, diffexpressed))

# select the most significant PLAC-linked ATAC peak per gene
heatmap_peaks <- candidate_full %>%
  semi_join(heatmap_genes, by = "symbol") %>%
  filter(!is.na(peak)) %>%
  distinct(symbol, peak, .keep_all = TRUE) %>%
  mutate(min_padj = pmin(atac_padj_Male,
                         atac_padj_Female,
                         na.rm = TRUE),
         min_padj = if_else(is.infinite(min_padj), NA_real_, min_padj)) %>%
  arrange(symbol, is.na(min_padj), min_padj, peak) %>%
  group_by(symbol) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  dplyr::select(symbol, peak, H3K27ac,
                starts_with("atac_log2FoldChange"),
                starts_with("atac_padj"))

# combine RNA, ATAC and H3K27ac data
heatmap_data <- heatmap_genes %>%
  left_join(heatmap_rna, by = "symbol") %>%
  left_join(heatmap_peaks, by = "symbol") %>%
  mutate(H3K27ac = case_when(is.na(peak) ~ "?",
                             H3K27ac %in% TRUE ~ "+",
                             TRUE ~ "-"))
# create RNA and ATAC matrices
rna_mat <- heatmap_data %>%
  dplyr::select(log2FoldChange_Male, log2FoldChange_Female) %>%
  as.matrix()

atac_mat <- heatmap_data %>%
  dplyr::select(atac_log2FoldChange_Male, atac_log2FoldChange_Female) %>%
  as.matrix()

rownames(rna_mat) <- rownames(atac_mat) <- heatmap_data$symbol
colnames(rna_mat) <- colnames(atac_mat) <- c("Male", "Female")

# create H3K27ac matrix
h3_mat <- matrix(heatmap_data$H3K27ac, nrow = 1,
                 dimnames = list("H3K27ac", heatmap_data$symbol))

# create significance matrices
rna_sig <- heatmap_data %>%
  transmute(Male = diffexpressed_Male %in% c("UP", "DOWN"),
            Female = diffexpressed_Female %in% c("UP", "DOWN")) %>%
  as.matrix()

atac_sig <- heatmap_data %>%
  transmute(Male = !is.na(atac_padj_Male) & atac_padj_Male < 0.05,
            Female = !is.na(atac_padj_Female) & atac_padj_Female < 0.05) %>%
  as.matrix()

# transpose matrices so genes are columns
rna_mat <- t(rna_mat)
atac_mat <- t(atac_mat)
rna_sig <- t(rna_sig)
atac_sig <- t(atac_sig)

# match the colours to your previous RNA heatmap
lfc_colours <- colorRamp2(c(-1, 0, 3), c("#2674B8", "white", "#D9364A"))
sex_colours <- c("Male" = "#78A8D1", "Female" = "#D78FA9")

# create sex annotations
row_ha <- rowAnnotation(Sex = c("Male", "Female"),
                        col = list(Sex = sex_colours),
                        show_annotation_name = FALSE,
                        simple_anno_size = unit(3, "mm"),
                        annotation_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                                                       labels_gp = gpar(fontsize = 8),
                                                       grid_width = unit(2.5, "mm"),
                                                       grid_height = unit(2.5, "mm")))

# RNA-seq heatmap
ht_rna <- Heatmap(rna_mat, 
                  name = "log2FC",
                  col = lfc_colours,
                  left_annotation = row_ha,
                  cluster_rows = FALSE,
                  column_title_gp = gpar(fontsize = 8, fontface = "bold"),
                  cluster_columns = FALSE,
                  show_row_names = TRUE,
                  row_names_side = "left",
                  row_names_gp = gpar(fontsize = 8),
                  show_column_names = TRUE,
                  column_names_side = "top",
                  column_names_rot = 90,
                  column_names_gp = gpar(fontsize = 8),
                  row_title = "RNA",
                  row_title_side = "left",
                  row_title_gp = gpar(fontsize = 8, fontface = "bold"),
                  height = unit(10.5, "mm"),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    if (rna_sig[i, j]) {
                      grid.text("*", x, y, gp = gpar(fontsize = 10, fontface = "bold"))
                      }
                    },
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 7),
                                              labels_gp = gpar(fontsize = 7),
                                              grid_width = unit(2.5, "mm"),
                                              grid_height = unit(2.5, "mm"),
                                              legend_height = unit(10, "mm")))

# ATAC-seq heatmap
ht_atac <- Heatmap(atac_mat,
                   name = "log2FC",
                   col = lfc_colours,
                   left_annotation = row_ha,
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   show_row_names = TRUE,
                   row_names_side = "left",
                   row_names_gp = gpar(fontsize = 8),
                   row_title = "ATAC",
                   row_title_side = "left",
                   row_title_gp = gpar(fontsize = 8, fontface = "bold"),
                   height = unit(10.5, "mm"),
                   show_column_names = FALSE,
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     if (atac_sig[i, j]) {
                       grid.text("*", x, y, gp = gpar(fontsize = 10, fontface = "bold"))
                       }
                     },
                   heatmap_legend_param = list(title_gp = gpar(fontsize = 7),
                                               labels_gp = gpar(fontsize = 7),
                                               grid_width = unit(2.5, "mm"),
                                               grid_height = unit(2.5, "mm"),
                                               legend_height = unit(10, "mm")))

# H3K27ac heatmap
ht_h3 <- Heatmap(h3_mat,
                 name = "H3K27ac",
                 col = c("+" = "white", "-" = "white", "?" = "white"),
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 show_row_names = TRUE,
                 row_names_side = "left",
                 row_names_gp = gpar(fontsize = 8, fontface = "bold"),
                 show_column_names = FALSE,
                 show_heatmap_legend = FALSE,
                 height = unit(1, "mm"),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(h3_mat[i, j], x, y, gp = gpar(fontsize = 9))
                   })

# combine heatmaps vertically
ht <- ht_rna %v% ht_atac %v% ht_h3

# plot and export
pdf("candidate_heatmap_landscape.pdf", width = 6.27, height = 1.8)
draw(ht,
     main_heatmap = "log2FC",
     cluster_columns = FALSE,
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()
