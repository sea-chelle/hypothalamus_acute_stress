# load required packages
library(here)
library(tidyverse)
library(grid)

# set working directory
setwd(here("plots"))

# calculate proportion of peaks overlapping H3K27ac and PLAC-seq regions
regulatory_overlap <- lapply(names(comparisons), function(comp) {
  
  comparisons[[comp]]$anno %>%
    summarise(H3K27ac = mean(H3K27ac) * 100,
              `PLAC-seq` = mean(PLAC) * 100) %>%
    mutate(comparison = comp)
  }) %>%
  bind_rows() %>%
  pivot_longer(cols = c(H3K27ac, `PLAC-seq`), names_to = "annotation", values_to = "percent")

# clean comparison labels and order
regulatory_overlap <- regulatory_overlap %>%
  mutate(comparison = factor(comparison, levels = c("consensus", "basal_vs_60", "male_basal_vs_60", "female_basal_vs_60"),
                             labels = c("Consensus", "Overall", "Male", "Female")),
         annotation = factor(annotation, levels = c("H3K27ac", "PLAC-seq")))
regulatory_overlap

# plot proportion of peaks in H3K27ac+ and PLAC-seq contact regions
ggplot(regulatory_overlap,aes(x = annotation, y = percent, fill = comparison)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("Consensus" = "#b3b3b3", "Overall" = "#A995C9", "Male" = "#78A8D1", "Female" = "#D78FA9")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Peaks overlapping\nannotation (%)", fill = NULL) +
  theme_classic() +
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.position = "top",
        legend.key.size = unit(0.3, "cm"),
        legend.box.spacing = unit(0.1, "cm"),
        legend.margin = margin(l = 0)) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  coord_flip()
ggsave("fig4.15a_DAR_H3K27ac_plac_proportion.tiff", height = 2, width = 2.1, dpi = 600, compression = "lzw")

# combine H3K27ac and PLAC-seq enrichment results
regulatory_enrichment <- bind_rows(h3k27ac_enrichment %>%
                                     transmute(comparison, 
                                               annotation = "H3K27ac",
                                               odds_ratio,
                                               CI_low = conf_low,
                                               CI_high = conf_high,
                                               pvalue,
                                               padj),
                                   
                                   plac_enrichment %>%
                                     transmute(comparison,
                                               annotation = "PLAC-seq",
                                               odds_ratio,
                                               CI_low,
                                               CI_high,
                                               pvalue,
                                               padj)) %>%
  mutate(comparison = factor(comparison, levels = c("basal_vs_60", "male_basal_vs_60", "female_basal_vs_60"),
                             labels = c("Overall", "Male", "Female")),
         annotation = factor(annotation, levels = c("H3K27ac", "PLAC-seq")))
regulatory_enrichment

# plot fisher's exact results for H3K27ac and PLAC-seq enrichment
ggplot(regulatory_enrichment, aes(x = odds_ratio, y = annotation)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0.15, colour = "black", linewidth = 0.3) +
  geom_point(aes(fill = padj < 0.05), shape = 21, size = 2.5, colour = "black", stroke = 0.5) +
  facet_wrap(~comparison, nrow = 1) +
  scale_x_log10() +
  scale_fill_manual(values = c(`TRUE` = "#A995C9", `FALSE` = "white"),
                    labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "NS")) +
  labs(x = "Odds ratio", y = NULL, fill = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.position = "top",
        legend.box.spacing = unit(0.1, "cm"),
        legend.key.size =  unit(0.3, "cm"),
        legend.margin = margin(l = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave("fig4.15b_DAR_H3K27ac_plac_enrichment.tiff", height = 2, width = 4.2, dpi = 600, compression = "lzw")
