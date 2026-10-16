# load required packages
library(here)
library(tidyverse)
library(grid)

# set working directory
setwd(here("plots"))

# reorder region proportion df
region_prop <- region_prop %>%
  mutate(comparison = factor(comparison,
                             levels = c("consensus", "basal_vs_60", "male_basal_vs_60", "female_basal_vs_60")),
         
         region = factor(region,
                         levels = c("Promoter", "UTR", "Exon", "Intron", "Downstream", "Distal intergenic")))


# plot proportion of peaks in each region
ggplot(region_prop, aes(x = region, y = proportion * 100, fill = comparison)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7 ) +
  labs(x = NULL, y = "Peaks (%)",fill = NULL) +
  scale_x_discrete(labels = c("Promoter", "UTR", "Exon", "Intron", "Downstream", "Distal\nintergenic")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("consensus" = "#b3b3b3",
                               "basal_vs_60" = "#A995C9",
                               "male_basal_vs_60" = "#78A8D1",
                               "female_basal_vs_60" = "#D78FA9"),
                    labels = c("consensus" = "Consensus",
                               "basal_vs_60" = "Overall",
                               "male_basal_vs_60" = "Male",
                               "female_basal_vs_60" = "Female")) +
  theme_classic() +
  theme(axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, "cm"))
ggsave("fig4.15a_region_proportion.tiff", height = 2.3, width = 6.3, dpi = 600, compression = "lzw")

# reorder enrichment results
enrichment_results$comparison <- factor(
  rep(c("basal_vs_60", "male_basal_vs_60", "female_basal_vs_60"), each = 6),
  levels = c("basal_vs_60", "male_basal_vs_60", "female_basal_vs_60"))

enrichment_results <- enrichment_results %>%
  mutate(comparison = factor(comparison, levels = c("basal_vs_60", "male_basal_vs_60", "female_basal_vs_60"),
                             labels = c("Overall", "Male", "Female")),
         region = factor(region, levels = c("Distal intergenic", "Downstream", "Intron", "Exon", "UTR", "Promoter")))

# plot enrichment results
ggplot(enrichment_results %>%
         filter(region != "Downstream"), aes(x = odds_ratio, y = region)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = 0.15, colour = "black", linewidth = 0.3) +
  geom_point(aes(fill = padj < 0.05), shape = 21, size = 2.5, colour = "black", stroke = 0.5) +
  facet_wrap(~comparison, nrow = 1) +
  scale_x_log10() +
  scale_fill_manual(values = c(`TRUE` = "#A995C9", `FALSE` = "white"),
                    labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "NS")) +
  labs(x = "Odds ratio", y = NULL, fill = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.heights = unit(5 * 0.6, "cm"),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.position = "right",
        legend.key.size = unit(0.3, "cm"),
        legend.box.spacing = unit(0, "pt"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave("fig4.15b_region_enrichment_plot.tiff", height = 1.8, width = 6.27, dpi = 600, compression = "lzw")


