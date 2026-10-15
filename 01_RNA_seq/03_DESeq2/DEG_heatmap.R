# load required packages
library(here)
library(DESeq2)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(tidyverse)
library(tibble)
library(ComplexHeatmap)
library(circlize)

# set working directory
setwd(here("plots"))

# get expression matrix and match to gene symbols
expr_mat <- assay(vsd)

# map Ensembl to gene symbols
gene_map <- AnnotationDbi::select(org.Mm.eg.db,
                                  keys = gsub("\\.\\d+$", "", rownames(vsd)),
                                  columns = c("SYMBOL"),
                                  keytype = "ENSEMBL")
gene_map <- gene_map[!duplicated(gene_map$ENSEMBL), ]
rownames(gene_map) <- gene_map$ENSEMBL

# add gene symbols to VSD
gene_symbols <- gene_map[gsub("\\.\\d+$", "", rownames(vsd)), "SYMBOL"]
vsd_symbol <- vsd
rownames(vsd_symbol) <- gene_symbols
vsd_symbol <- vsd_symbol[!is.na(rownames(vsd_symbol)), ]

# extract group info
group_labels <- colnames(vsd_symbol)
group_info <- do.call(rbind, strsplit(group_labels, "\\.")) # splits into ID, sex, time
group_df <- data.frame(Sample = group_labels,
                       Sex = ifelse(group_info[,2] == "f", "Female", "Male"),
                       Time = gsub("min", " min", group_info[,3]),
                       Group = paste(ifelse(group_info[,2] == "f", "Female", "Male"),
                                     gsub("min", " min", group_info[,3])))

# average expression
vsd_mat <- assay(vsd_symbol)
grouped_expr <- as.data.frame(vsd_mat) %>%
  t() %>%
  as.data.frame() %>%
  mutate(Group = group_df$Group) %>%
  group_by(Group) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("Group") %>%
  t()

# reorder genes by category
gene_ordered <- gene_categories %>% 
  arrange(factor(category, levels = unique(category))) %>%
  pull(gene)

# subset and reorder expression matrix
expr_subset <- grouped_expr[gene_ordered, ]

#  column order
desired_cols <- c("Male basal", "Male 10", "Male 60",
                  "Female basal", "Female 10", "Female 60")

expr_subset <- expr_subset[, desired_cols, drop = FALSE]

# make column labels = time only (basal -> 0)
time_labels <- gsub(".* ", "", colnames(expr_subset))  # basal / 10 / 60
time_labels[time_labels == "basal"] <- "0"

# keep only Sex as the top annotation (no time annotation -> no time legend)
sex_labels <- sub(" .*", "", colnames(expr_subset))  # male / female
expr_subset <- t(scale(t(expr_subset)))

# column annotation
col_ha <- HeatmapAnnotation(Sex = sex_labels,
                            col = list(Sex = c("Male" = "#78A8D1", "Female" = "#D78FA9")),  
                            annotation_name_side = "left",  
                            annotation_name_gp = gpar(fontsize = 8),
                            show_annotation_name = FALSE,
                            simple_anno_size = unit(3, "mm"),
                            annotation_legend_param = list(title_gp = gpar(fontsize = 8),  
                                                           labels_gp = gpar(fontsize = 8),
                                                           annotation_height = unit(3, "cm"),
                                                           grid_width = unit(3, "mm"),
                                                           grid_height = unit(3, "mm")))

# category colours
category_colours <- c("IEG" = "#5FA9A2",  
                      "Stress" = "#A995C9", 
                      "Neuronal" = "#DD805E",
                      "Endocrine" = "#A3A64C",
                      "Immune" = "#D1B654",
                      "Developmental" = "#858B91")

# row annotation
row_ha <- rowAnnotation(Category = gene_categories$category[match(rownames(expr_subset), gene_categories$gene)],
                        col = list(Category = category_colours),
                        show_annotation_name = FALSE,
  annotation_name_side = NULL,  
  annotation_name_gp = gpar(fontsize = 0),
  simple_anno_size = unit(3, "mm"),
  annotation_legend_param = list(title_gp = gpar(fontsize = 8),  
                                 labels_gp = gpar(fontsize = 8),
                                 grid_width = unit(3, "mm"),
                                 grid_height = unit(3, "mm")))


heatmap <- Heatmap(expr_subset,
                   name = "Z-score",
                   col = colorRamp2(c(-2, 0, 2), c("#2674B8", "white", "#D9364A")),
                   top_annotation = col_ha,
                   left_annotation = row_ha,     
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   row_names_side = "left",
                   column_labels = time_labels,
                   column_names_gp = gpar(fontsize = 8),
                   column_names_side = "bottom",
                   column_names_rot = 0,
                   column_title = "Time (min)",
                   column_title_side = "bottom",
                   column_title_gp = gpar(fontsize = 8),
                   row_names_gp = gpar(fontsize = 8),
                   heatmap_legend_param = list(title_gp = gpar(fontsize = 7),
                                               labels_gp = gpar(fontsize = 7),
                                               grid_width = unit(3, "mm"),
                                               grid_height = unit(3, "mm"),
                                               at = c(-2, 0, 2)))


png("heatmap.png", width = 3.2, height = 5, units = "in", res = 600)  
draw(heatmap, 
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off() 

