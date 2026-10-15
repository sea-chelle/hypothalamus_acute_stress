# load required libraries
library(here)
library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(grid)

# set working directory

# run PCA on all samples
pcaData <- plotPCA(vsd,
                   intgroup = c("sex", "condition"),
                   returnData = TRUE)

# calculate percentage variance
percentVar <- round(100 * attr(pcaData, "percentVar"))

# simplify sample name
pcaData$name <- sub("-.*", "", pcaData$name)

# plot PCA by condition
ggplot(pcaData, aes(PC1, PC2, colour = group, label = name)) +
  geom_point(size = 1.5) +
  geom_text(vjust = -0.7, size = 2.5, show.legend = FALSE) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(values = c("male:basal" = "#5FA9A2",
                                "male:60" = "#78A8D1",
                                "female:basal" = "#A995C9",
                                "female:60" = "#D78FA9"),
                     labels = c("male:basal" = "Male basal",
                                "male:60" = "Male 60 min ARS",
                                "female:basal" = "Female basal",
                                "female:60" = "Female 60 min ARS")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  coord_fixed() +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        legend.box.spacing = unit(0, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("fig4.12_atac_pca.tiff", height = 3.5, width = 5, dpi = 600, compression = "lzw")

# turn VST-normalised counts into a matrix
vsd_mat <- assay(vsd)

# calculate Euclidean distance between samples
sample_dists <- dist(t(vsd_mat))
sample_dists <- as.matrix(sample_dists)

# use animal number only for heatmap labels
sample_ids <- sub("-.*", "", colnames(sample_dists))

colnames(sample_dists) <- sample_ids
rownames(sample_dists) <- sample_ids

# create annotation dataframe from biological metadata
annotation_df <- data.frame(Sex = colData$sex,
                            Condition = colData$condition)

# match annotation row names to animal-number labels
rownames(annotation_df) <- sub("-.*", "", rownames(colData))

# plot heatmap
png("sample_heatmap.png", height = 1500, width = 2000, res = 300)
Heatmap(sample_dists,
        name = "Distance",
        col = colorRamp2(c(min(sample_dists), max(sample_dists)), c("#fff7ec", "#D9364A")),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        top_annotation = HeatmapAnnotation(df = annotation_df,
                                           col = list(Sex = c("female" = "#D78FA9", "male" = "#78A8D1"), 
                                                      Condition = c("basal" = "#A995C9", "60" = "#5FA9A2")), 
                                           annotation_name_side = "left"),
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 11),
        column_names_gp = gpar(fontsize = 11),
        column_title = "Sample-to-Sample Distance",
        heatmap_legend_param = list(title = "Euclidean Distance"))
dev.off()
