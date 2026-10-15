# load required packages
library(here)
library(ComplexHeatmap)
library(circlize)
library(grid)

# set working directory
setwd(here("plots"))

# endocrine genes to plot on heatmap
genes_endocrine <- c(
  
  # neuroendocrine 
  "Crh", "Avp", "Oxt", "Trh", "Gnrh1", "Ghrh", 
  
  # neuropeptide
  "Agrp", "Pomc", "Sst", "Npy", "Tac1", "Tac2",
  
  # urocortins
  "Ucn", "Ucn2", "Ucn3")

# endocrine receptor genes to plot on heatmap
genes_receptors <- c(
  # neuroendocrine receptors
  "Crhr1", "Crhr2", "Avpr1a", "Oxtr", "Trhr","Gnrhr", "Ghr", 
  
  # hypothalamic / peptide receptors
  "Tacr2", "Mc3r", "Mc4r", "Sstr1", "Sstr4",
  
  # gonadal steroid receptors
  "Esr1", "Esr2", "Pgr", "Ar", 
  
  # other receptors
  "Prlr", "Nr3c1", "Nr3c2")

# labels for endocrine heatmap
endocrine_blocks <- c(rep("Neuroendocrine", 6),
                      rep("Neuropeptide", 6),
                      rep("Urocortins", 3))
names(endocrine_blocks) <- genes_endocrine

# labels for receptors heatmap
receptor_blocks <- c(rep("Neuroendocrine", 7),
                     rep("Neuropeptide", 5),
                    rep("Steroid\nhormone", 4),
                    rep("Other", 3))
names(receptor_blocks) <- genes_receptors

# create function for plotting heat map
plot_supervised_heatmap <- function(genes, gene_blocks) {
  
  # keep genes present
  genes_present <- genes[genes %in% rownames(grouped_expr)]
  
  # report missing genes
  missing_genes <- setdiff(genes, rownames(grouped_expr))
  if(length(missing_genes) > 0) {
    message("Missing genes: ", paste(missing_genes, collapse = ", "))
  }
  
  # subset expression
  expr_subset <- grouped_expr[genes_present, , drop = FALSE]
  
  # column order
  desired_cols <- c("Male basal", "Male 10", "Male 60",
                    "Female basal", "Female 10", "Female 60")
  expr_subset <- expr_subset[, desired_cols, drop = FALSE]
  
  # labels
  time_labels <- gsub(".* ", "", colnames(expr_subset))
  time_labels[time_labels == "basal"] <- "0"
  sex_labels <- sub(" .*", "", colnames(expr_subset))
  
  # row Z-score
  expr_subset <- t(scale(t(expr_subset)))
  
  # transpose so genes are columns and groups are rows
  expr_subset <- t(expr_subset)
  
  # column blocks
  column_blocks <- factor(gene_blocks[colnames(expr_subset)],
                          levels = unique(gene_blocks[genes_present]))
  
  # sex annotation
  row_ha <- rowAnnotation(Sex = sex_labels,
                          col = list(Sex = c("Male" = "#78A8D1", "Female" = "#D78FA9")),
                          show_annotation_name = FALSE,
                          simple_anno_size = unit(3, "mm"),
                          annotation_legend_param = list(Sex = list(nrow = 1, by_row = TRUE, 
                                                                    title_position = "topcenter",
                                                                    title_gp = gpar(fontsize = 7),
                                                                    labels_gp = gpar(fontsize = 7),
                                                                    grid_width = unit(3, "mm"),
                                                                    grid_height = unit(3, "mm"))))
  
  # make cell size the same across all heatmaps
  cell_size <- 4
  
  # plot heatmap
  heatmap <- Heatmap(expr_subset,
                     name = "Z-score",
                     left_annotation = row_ha,
                     col = colorRamp2(c(-2, 0, 2), c("#2674B8", "white", "#D9364A")),
                     column_split = column_blocks,
                     column_gap = unit(1.5, "mm"),
                     cluster_rows = FALSE,
                     cluster_columns = FALSE,
    
    # number of genes × 4 mm wide
    width = unit(ncol(expr_subset) * cell_size, "mm"),
    
    # 6 groups × 4 mm tall
    height = unit(nrow(expr_subset) * cell_size, "mm"),
    
    # gene names
    column_names_side = "bottom",
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 8),
    
    # gene block labels
    column_title_side = "top",
    column_title_rot = 0,
    column_title_gp = gpar(fontsize = 8, col = "grey30"),
    
    # time labels
    row_labels = time_labels,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8),
    
    row_title = "Time (min)",
    row_title_side = "left",
    row_title_gp = gpar(fontsize = 8),
    
    heatmap_legend_param = list(direction = "horizontal",
                                title_position = "topcenter",
                                title_gp = gpar(fontsize = 7),
                                labels_gp = gpar(fontsize = 7),
                                legend_height = unit(15, "mm"),
                                grid_width = unit(3, "mm"),
                                grid_height = unit(3, "mm"),
                                at = c(-2, 0, 2)))
  
  return(heatmap)
}

# plot endocrine heatmap
endocrine_hm <- plot_supervised_heatmap(genes_endocrine,
                                        endocrine_blocks)

# plot hormone receptor heatmap
receptors_hm <- plot_supervised_heatmap(genes_receptors,
                                        receptor_blocks)
receptors_hm@left_annotation <- NULL

# combine and plot
endocrine_receptors_hm <- endocrine_hm + receptors_hm

png("fig8_supervised_heatmap.png", height = 2.2, width = 6.27, units = "in", res = 300)
draw(endocrine_receptors_hm,
     ht_gap = unit(5, "mm"),
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     align_heatmap_legend = "global_center",
     align_annotation_legend = "global_center",
     merge_legends = TRUE)
dev.off()
