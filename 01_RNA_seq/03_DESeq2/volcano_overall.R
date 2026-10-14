# load required packages
library(EnhancedVolcano)
library(ggplot2)

# set working directory
setwd(here("plots"))

# create a vector of DEGs to label
labels <- c("Fos", "Arc", "Nr4a1", "Egr1", "Npas4")

# in labels, show if label is in the list, if not leave blank
labels <- ifelse(res_basal_vs_10_df$symbol %in% labels, 
                 res_basal_vs_10_df$symbol, "")

# plot volcano plot
volcano_basal_vs_10 <- EnhancedVolcano(res_basal_vs_10_df, x = "log2FoldChange", y = "padj",
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
  ggtitle("Basal vs 10 min") +
  theme_bw() +
  theme(plot.caption = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "none")
ggsave("fig4.2a_basal_vs_10_volcano.tiff", height = 3, width = 3, dpi = 300)

# create a vector of DEGs to label
labels <- c("Fos", "Dusp1", "Sgk1", "Tsc22d3", "Ddit4", "Per1", "Dio2", "Nfkbia")

# in labels, show if label is in the list, if not leave blank
labels <- ifelse(res_basal_vs_60_df$symbol %in% labels, 
                 res_basal_vs_60_df$symbol, "")

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
  ggtitle("Basal vs 60 min") +
  theme_bw() +
  theme(plot.caption = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "none")
ggsave("fig4.2b_basal_vs_60_volcano.tiff", height = 3, width = 3, dpi = 300)


# create a vector of DEGs to label
labels <- c("Npas4", "Sgk1", "Tsc22d3", "Ddit4", "Per1", "Nfkbia")

# in labels, show if label is in the list, if not leave blank
labels <- ifelse(res_10_vs_60_df$symbol %in% labels, 
                 res_10_vs_60_df$symbol, "")

# plot volcano plot
volcano_10_vs_60 <- EnhancedVolcano(res_10_vs_60_df, x = "log2FoldChange", y = "padj",
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
  ggtitle("10 vs 60 min") +
  theme_bw() +
  theme(plot.caption = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "none")
ggsave("fig4.2c_10_vs_60_volcano.tiff", height = 6, width = 6, dpi = 300)
