# load required packages
library(here)
library(EnhancedVolcano)
library(ggplot2)
library(patchwork)
library(cowplot)
library(grid)

# set working directory 
setwd(here("plots"))

# genes to label
genes_to_label <- c("Cebpd", "Fkbp5", "Notch1", "Slc1a2", "Ddit4", "Kdm3a",  "Hcrtr1")

# filter based on genes to label
labels <- peakAnno_basal_vs_60_df$SYMBOL[match(rownames(res_basal_vs_60_df), 
                                               peakAnno_basal_vs_60_df$peak)]
labels <- ifelse(labels %in% genes_to_label, labels, "")

# for duplicate gene labels, priortise promoter peak first, otherwise largest absolute log2FC
for (gene in genes_to_label) {
  idx <- which(labels == gene)
  
  if (length(idx) > 1) {
    anno <- peakAnno_basal_vs_60_df$annotation[
      match(rownames(res_basal_vs_60_df)[idx],
            peakAnno_basal_vs_60_df$peak)
    ]
    
    promoter_idx <- idx[grepl("^Promoter", anno)]
    
    if (length(promoter_idx) > 0) {
      keep <- promoter_idx[which.max(abs(res_basal_vs_60_df$log2FoldChange[promoter_idx]))]
    } else {
      keep <- idx[
        which.max(abs(res_basal_vs_60_df$log2FoldChange[idx]))]
    }
    
    labels[setdiff(idx, keep)] <- ""
  }
}

# plot volcano plot
volcano_basal_vs_60 <- EnhancedVolcano(res_basal_vs_60_df, x = "log2FoldChange", y = "padj",
                                       lab = labels,
                                       pCutoff = 0.05,
                                       FCcutoff = 0.5,
                                       cutoffLineWidth = 0.3,
                                       labSize = 2,
                                       pointSize = 1,
                                       shape = 16,
                                       colAlpha = 1,
                                       legendIconSize = 2,
                                       col = c("#b3b3b3", "#A995C9", "#DD805E", "#5FA9A2"),
                                       subtitle = NULL,
                                       drawConnectors = T,
                                       max.overlaps = Inf) +
  ggtitle("Overall") +
  theme_bw() +
  theme(plot.caption = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "none")

# genes to label
genes_to_label <- c("Cebpd","Ddit4", "Fkbp5", "Kdm3a", "Slc6a1", "Mobp")

# filter based on genes to label
labels <- peakAnno_male_basal_vs_60_df$SYMBOL[match(rownames(res_male_basal_vs_60_df), 
                                                    peakAnno_male_basal_vs_60_df$peak)]
labels <- ifelse(labels %in% genes_to_label, labels, "")

# for duplicate gene labels, priortise promoter peak first, otherwise largest absolute log2FC
for (gene in genes_to_label) {
  idx <- which(labels == gene)
  
  if (length(idx) > 1) {
    anno <- peakAnno_male_basal_vs_60_df$annotation[
      match(rownames(res_male_basal_vs_60_df)[idx],
            peakAnno_male_basal_vs_60_df$peak)
    ]
    
    promoter_idx <- idx[grepl("^Promoter", anno)]
    
    if (length(promoter_idx) > 0) {
      keep <- promoter_idx[which.max(abs(res_male_basal_vs_60_df$log2FoldChange[promoter_idx]))]
    } else {
      keep <- idx[
        which.max(abs(res_male_basal_vs_60_df$log2FoldChange[idx]))]
    }
    
    labels[setdiff(idx, keep)] <- ""
  }
}

# plot volcano plot
volcano_male_basal_vs_60 <- EnhancedVolcano(res_male_basal_vs_60_df, x = "log2FoldChange", y = "padj",
                                            lab = labels,
                                            pCutoff = 0.05,
                                            FCcutoff = 0.5,
                                            cutoffLineWidth = 0.3,
                                            labSize = 2,
                                            pointSize = 1,
                                            shape = 16,
                                            colAlpha = 1,
                                            legendIconSize = 2,
                                            col = c("#b3b3b3", "#A995C9", "#DD805E", "#5FA9A2"),
                                            subtitle = NULL,
                                            drawConnectors = T,
                                            max.overlaps = Inf) +
  ggtitle("Male") +
  theme_bw() +
  theme(plot.caption = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "none")

# genes to label
genes_to_label <- c("Ddit4", "Fkbp5", "Mobp")

# filter based on genes to label
labels <- peakAnno_female_basal_vs_60_df$SYMBOL[match(rownames(res_female_basal_vs_60_df), 
                                                      peakAnno_female_basal_vs_60_df$peak)]
labels <- ifelse(labels %in% genes_to_label, labels, "")

# for duplicate gene labels, priortise promoter peak first, otherwise largest absolute log2FC
for (gene in genes_to_label) {
  idx <- which(labels == gene)
  
  if (length(idx) > 1) {
    anno <- peakAnno_female_basal_vs_60_df$annotation[
      match(rownames(res_female_basal_vs_60_df)[idx],
            peakAnno_female_basal_vs_60_df$peak)
    ]
    
    promoter_idx <- idx[grepl("^Promoter", anno)]
    
    if (length(promoter_idx) > 0) {
      keep <- promoter_idx[which.max(abs(res_female_basal_vs_60_df$log2FoldChange[promoter_idx]))]
    } else {
      keep <- idx[
        which.max(abs(res_female_basal_vs_60_df$log2FoldChange[idx]))]
    }
    
    labels[setdiff(idx, keep)] <- ""
  }
}

# plot volcano
volcano_female_basal_vs_60 <- EnhancedVolcano(res_female_basal_vs_60_df, x = "log2FoldChange", y = "padj",
                                              lab = labels,
                                              pCutoff = 0.05,
                                              FCcutoff = 0.5,
                                              cutoffLineWidth = 0.3,
                                              labSize = 2,
                                              pointSize = 1,
                                              shape = 16,
                                              colAlpha = 1,
                                              legendIconSize = 2,
                                              col = c("#b3b3b3", "#A995C9", "#DD805E", "#5FA9A2"),
                                              subtitle = NULL,
                                              drawConnectors = T,
                                              max.overlaps = Inf) +
  ggtitle("Female") +
  theme_bw() +
  theme(plot.caption = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "none")

# plot volcanoes all togehter
(plot_spacer() | volcano_basal_vs_60) /
  (volcano_male_basal_vs_60 | volcano_female_basal_vs_60)
ggsave("fig4.13b-d_DAR_volcano.tiff", height = 5.5, width = 6.27, dpi = 600, compression = "lzw")

# extract vertical volcano legend
volcano_legend <- get_legend(volcano_basal_vs_60 +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 8),
          legend.spacing.x = unit(0.1, "cm"),
          legend.spacing.y = unit(0, "cm"),
          legend.key.width = unit(0.3, "cm"),
          legend.key.height = unit(0.25, "cm")) +
    guides(colour = guide_legend(nrow = 2,
                                 byrow = TRUE)))

# convert legend grob to drawable object
volcano_legend <- ggdraw() + draw_grob(volcano_legend)

# export with transparent background
ggsave("fig4.13_volcano_legend.tiff", volcano_legend,
       width = 2.8, height = 0.6, units = "in", dpi = 600, bg = "transparent", compression = "lzw")
