# load required packages
library(here)
library(tidyverse)
library(grid)

# set working directory
setwd(here("plots"))

# plot
ggplot(cre_overlap, aes(x = celltype, y = percent, fill = comparison)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("Consensus" = "#b3b3b3",  "Overall" = "#A995C9", "Male" = "#78A8D1", "Female" = "#D78FA9")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Peaks overlapping CREs (%)", fill = NULL) +
  theme_classic() +
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.position = "right",
        legend.box.spacing = unit(0.1, "cm"),
        legend.key.size =  unit(0.3, "cm"),
        legend.margin = margin(l = 0))
ggsave("fig4.14c_snATACseq_proportion.tiff", height = 2, width = 6.3, dpi = 600, compression = "lzw")

# labels for plotting
cre_enrichment <- cre_enrichment %>%
  mutate(comparison = factor(comparison,
                             levels = c("basal_vs_60", "male_basal_vs_60", "female_basal_vs_60"),
                             labels = c("Overall", "Male", "Female")),
         celltype = factor(celltype, levels = c("astrocyte", 
                                                "endothelial", 
                                                "gaba", 
                                                "glut", 
                                                "microglia", 
                                                "neuron", 
                                                "oligodendrocyte", 
                                                "OPC"),
                           labels = c("Astrocyte", 
                                      "Endothelial", 
                                      "GABAergic", 
                                      "Glutamatergic", 
                                      "Microglia", 
                                      "Neuron", 
                                      "Oligodendrocyte", 
                                      "OPC")))

# plot
ggplot(cre_enrichment, aes(x = odds_ratio, y = celltype)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.15, linewidth = 0.3) +
  geom_point(aes(fill = padj < 0.05), shape = 21, size = 2.5, colour = "black", stroke = 0.5) +
  facet_wrap(~comparison, nrow = 1) +
  scale_x_log10() +
  scale_fill_manual(values = c(`TRUE` = "#A995C9", `FALSE` = "white"),
                    labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "NS")) +
  labs(x = "Odds ratio", y = NULL, fill = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.heights = unit(8 * 0.5, "cm"),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.box.spacing = unit(0.1, "cm"),
        legend.key.size =  unit(0.3, "cm"),
        legend.margin = margin(l = 0),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave("fig4.15d_snATAC_enrichment_plot.tiff", height = 2.2, width = 6.27, dpi = 600, compression = "lzw")
