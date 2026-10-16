# load required libraries
library(here)
library(tidyverse)
library(ggrepel)
library(grid)

# set working directory
setwd(here("plots"))

# combine nearest-TSS DEG pairs
nearest_plot <- bind_rows(nearest_deg, .id = "comparison") %>%
  mutate(comparison = dplyr::recode(comparison, 
                                    basal_vs_60 = "Overall", 
                                    male_basal_vs_60 = "Male",
                                    female_basal_vs_60 = "Female"),
         comparison = factor(comparison, levels = c("Overall", "Male", "Female")),
         concordance = if_else(sign(log2FoldChange) == sign(rna_log2FoldChange),
                               "Concordant", "Discordant"),
         label = if_else(concordance == "Discordant", SYMBOL, NA_character_)) 

# plot
nearest_concordance_plot <- ggplot(nearest_plot, aes(x = log2FoldChange, y = rna_log2FoldChange, colour = comparison)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_text_repel(aes(label = label), 
                  size = 2.5, 
                  color = "black",
                  show.legend = FALSE,
                  na.rm = TRUE,
                  box.padding = 0.4,
                  point.padding = 0.3,
                  min.segment.length = 0) +
  scale_colour_manual(values = c("Overall" = "#A995C9", "Male" = "#78A8D1", "Female" = "#D78FA9")) +
  labs(x = "ATAC-seq log2 fold change", y = "RNA-seq log2 fold change", colour = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.position = "top",
        legend.text = element_text(size = 7),
        legend.box.spacing = unit(0.2, "cm"),
        legend.key.size =  unit(0.3, "cm"),
        legend.margin = margin(l = 0))

nearest_concordance_plot

ggsave("fig4.16a_RNA_ChIPseeker.tiff", height = 4, width = 4.5, dpi = 600, compression = "lzw")

# combine PLAC-seq DEG pairs
plac_plot <- bind_rows(plac_deg, .id = "comparison") %>%
  mutate(comparison = dplyr::recode(comparison, 
                                    basal_vs_60 = "Overall", 
                                    male_basal_vs_60 = "Male",
                                    female_basal_vs_60 = "Female"),
         comparison = factor(comparison, levels = c("Overall", "Male", "Female")),
         concordance = if_else(sign(log2FoldChange) == sign(rna_log2FoldChange),
                               "Concordant", "Discordant"),
         label = if_else(concordance == "Discordant", GeneName, NA_character_)) 

# plot
plac_concordance_plot <- ggplot(plac_plot, aes(x = log2FoldChange, y = rna_log2FoldChange, colour = comparison)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_text_repel(aes(label = label), 
                  size = 2.5, 
                  color = "black",
                  show.legend = FALSE,
                  na.rm = TRUE,
                  box.padding = 0.4,
                  point.padding = 0.3,
                  min.segment.length = 0) +
  scale_colour_manual(values = c("Overall" = "#A995C9", "Male" = "#78A8D1", "Female" = "#D78FA9")) +
  labs(x = "ATAC-seq log2 fold change", y = "RNA-seq log2 fold change", colour = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.position = "top",
        legend.text = element_text(size = 7),
        legend.box.spacing = unit(0.2, "cm"),
        legend.key.size =  unit(0.3, "cm"),
        legend.margin = margin(l = 0))

plac_concordance_plot

ggsave("fig4.16b_RNA_PLAC.tiff", height = 4, width = 4.5, dpi = 600, compression = "lzw")
