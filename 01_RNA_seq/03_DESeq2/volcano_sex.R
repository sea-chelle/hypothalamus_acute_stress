# load required packages
library(EnhancedVolcano)
library(ggplot2)

# set working directory
setwd(here("plots"))

# create a vector of DEGs to label
labels <- c("Iapp", "Xist", "Uty", "Ddx3y", "Kdm5d", "Sytl4") 

# in labels, show if label is in the list, if not leave blank
labels <- ifelse(res_basal_male_vs_female_df$symbol %in% labels, 
                 res_basal_male_vs_female_df$symbol, "")

# plot volcano plot
volcano_basal_male_vs_female <- EnhancedVolcano(res_basal_male_vs_female_df, x = "log2FoldChange", y = "padj",
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
  ggtitle("Basal") +
  theme_bw() +
  theme(plot.caption = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "none")

ggsave(plot = volcano_basal_male_vs_female, "fig4.5a_basal_male_vs_female_volcano.tiff", height = 3, width = 3, dpi = 600, compression = "lzw")

# create a vector of DEGs to label
labels <- c("Iapp", "Sytl4", "Xist", "Uty", "Ddx3y", "Kdm5d", "Kiss1", "Greb1") 

# in labels, show if label is in the list, if not leave blank
labels <- ifelse(res_10min_male_vs_female_df$symbol %in% labels, 
                 res_10min_male_vs_female_df$symbol, "")

# plot volcano plot
volcano_10min_male_vs_female <- EnhancedVolcano(res_10min_male_vs_female_df, x = "log2FoldChange", y = "padj",
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
  ggtitle("10 min ARS") +
  theme_bw() +
  theme(plot.caption = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 9),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "none")

ggsave(plot = volcano_10min_male_vs_female, "fig4.5b_10min_male_vs_female_volcano.tiff", height = 3, width = 3, dpi = 600, compression = "lzw")

