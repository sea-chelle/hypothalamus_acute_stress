# load required packages
library(here)
library(dplyr)
library(ggplot2)

# set working directory
setwd(here("plot"))

# create vector of selected terms to plot 
KEGG_basal_vs_10_select <- c("IL-17 signaling pathway", "MAPK signaling pathway", "TNF signaling pathway",
                             "Growth hormone synthesis, secretion and action")

# filter for selected terms to plot and calculate GeneRatio and -log10 FDR
KEGG_basal_vs_10_plot <- KEGG_basal_vs_10@result %>%
  filter(Description %in% KEGG_basal_vs_10_select) %>%
  mutate(GeneRatio = sapply(GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  }),
  neglog10FDR = -log10(p.adjust),
  PlotLabel = case_when(Description == "IL-17 signaling pathway" ~ "IL-17 signalling",
                        Description == "MAPK signaling pathway" ~ "MAPK signalling",
                        Description == "TNF signaling pathway" ~ "TNF signalling",
                        Description == "Growth hormone synthesis, secretion and action" ~ "GH signalling",
                        Description == "FoxO signaling pathway" ~ "FoxO signalling",
                        Description == "NF-kappa B signalling pathway" ~ "NF-κB signalling",
                        TRUE ~ Description))

# create vector of selected terms to plot 
KEGG_upregulated_basal_vs_60_select <- c("IL-17 signaling pathway", "MAPK signaling pathway", "TNF signaling pathway", "FoxO signaling pathway",
                                         "NF-kappa B signaling pathway", "Apoptosis")

# filter for selected terms to plot and calculate GeneRatio and -log10 FDR
KEGG_upregulated_basal_vs_60_plot <- KEGG_upregulated_basal_vs_60@result %>%
  filter(Description %in% KEGG_upregulated_basal_vs_60_select) %>%
  mutate(GeneRatio = sapply(GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  }),
  neglog10FDR = -log10(p.adjust),
  PlotLabel = case_when(Description == "IL-17 signaling pathway" ~ "IL-17 signalling",
                        Description == "MAPK signaling pathway" ~ "MAPK signalling",
                        Description == "TNF signaling pathway" ~ "TNF signalling",
                        Description == "Growth hormone synthesis, secretion and action" ~ "GH signalling",
                        Description == "FoxO signaling pathway" ~ "FoxO signalling",
                        Description == "NF-kappa B signaling pathway" ~ "NF-κB signalling",
                        TRUE ~ Description))

# create vector of selected terms to plot 
KEGG_upregulated_10_vs_60_select <- c("IL-17 signaling pathway", "FoxO signaling pathway")

# filter for selected terms to plot and calculate GeneRatio and -log10 FDR
KEGG_upregulated_10_vs_60_plot <- KEGG_upregulated_10_vs_60@result %>%
  filter(Description %in% KEGG_upregulated_10_vs_60_select) %>%
  mutate(GeneRatio = sapply(GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  }),
  neglog10FDR = -log10(p.adjust),
  PlotLabel = case_when(Description == "IL-17 signaling pathway" ~ "IL-17 signalling",
                        Description == "MAPK signaling pathway" ~ "MAPK signalling",
                        Description == "TNF signaling pathway" ~ "TNF signalling",
                        Description == "Growth hormone synthesis, secretion and action" ~ "GH signalling",
                        Description == "FoxO signaling pathway" ~ "FoxO signalling",
                        Description == "NF-kappa B signalling pathway" ~ "NF-κB signalling",
                        TRUE ~ Description))


# make FDR, gene count and gene ratio scales the same across all upregulated plots
fdr_limits <- range(c(KEGG_basal_vs_10_plot$neglog10FDR,
                      KEGG_upregulated_basal_vs_60_plot$neglog10FDR,
                      KEGG_upregulated_10_vs_60_plot$neglog10FDR))

count_limits <- range(c(KEGG_basal_vs_10_plot$Count,
                        KEGG_upregulated_basal_vs_60_plot$Count,
                        KEGG_upregulated_10_vs_60_plot$Count))

gene_ratio_max <- max(c(KEGG_basal_vs_10_plot$GeneRatio,
                        KEGG_upregulated_basal_vs_60_plot$GeneRatio,
                        KEGG_upregulated_10_vs_60_plot$GeneRatio))

# combine all upregulated comparisons
KEGG_upregulated_combined <- bind_rows(KEGG_basal_vs_10_plot %>%
                                       mutate(Comparison = "Basal vs 10 min ARS"),
                                       KEGG_upregulated_basal_vs_60_plot %>%
                                       mutate(Comparison = "Basal vs 60 min ARS"),
                                       KEGG_upregulated_10_vs_60_plot %>%
                                       mutate(Comparison = "10 vs 60 min ARS"))

# order terms
KEGG_upregulated_combined$PlotLabel <- factor(KEGG_upregulated_combined$PlotLabel,
                                              levels = rev(c("IL-17 signalling", "MAPK signalling", "TNF signalling", "GH signalling", 
                                                             "FoxO signalling", "NF-κB signalling", "Apoptosis")))

# order comparisons
KEGG_upregulated_combined$Comparison <- factor(KEGG_upregulated_combined$Comparison,
                                             levels = c("Basal vs 10 min ARS", "Basal vs 60 min ARS", "10 vs 60 min ARS"),
                                             labels = c("Basal vs 10 min", "Basal vs 60 min", "10 vs 60 min"))

# plot all three upregulated dotplots into one plot
ggplot(KEGG_upregulated_combined,
       aes(x = GeneRatio, y = PlotLabel, 
           size = Count,
           colour = neglog10FDR)) +
  geom_point() +
  facet_grid(. ~ Comparison, scales = "fixed") +
  scale_colour_gradient(low = "#FADADA", high = "#D9364A") +
  scale_size_continuous(breaks = scales::breaks_width(1)) +
  scale_x_continuous(breaks = seq(0, 0.3, 0.1), limits = c(0, 0.3)) +
  labs(x = "GeneRatio", y = NULL, 
       size = "Count", 
       colour = expression(-log[10]("FDR"))) +
  theme_bw() +
  theme(panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold"),
        axis.text = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(3, "mm"),
        legend.key.height = unit(3, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "pt"),
        panel.spacing.x = unit(0.5, "lines"))

ggsave("fig4.4._KEGG_upregulated_combined.tiff", width = 6.27, height = 2.8, dpi = 300)

# create vector of selected terms to plot 
KEGG_downregulated_basal_vs_60_select <- c("Necroptosis")

# filter for selected terms to plot and calculate GeneRatio and -log10 FDR
KEGG_downregulated_basal_vs_60_plot <- KEGG_downregulated_basal_vs_60@result %>%
  filter(Description %in% KEGG_downregulated_basal_vs_60_select) %>%
  mutate(GeneRatio = sapply(GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  }),
  neglog10FDR = -log10(p.adjust))

# create vector of selected terms to plot 
KEGG_downregulated_10_vs_60_select <- c("Necroptosis")

# filter for selected terms to plot and calculate GeneRatio and -log10 FDR
KEGG_downregulated_10_vs_60_plot <- KEGG_downregulated_10_vs_60@result %>%
  filter(Description %in% KEGG_downregulated_10_vs_60_select) %>%
  mutate(GeneRatio = sapply(GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  }),
  neglog10FDR = -log10(p.adjust))

# make FDR, gene count and gene ratio scales the same across all downregulated plots
fdr_limits <- range(c(KEGG_downregulated_basal_vs_60_plot$neglog10FDR,
                      KEGG_downregulated_10_vs_60_plot$neglog10FDR))

count_limits <- range(c(KEGG_downregulated_basal_vs_60_plot$Count,
                        KEGG_downregulated_10_vs_60_plot$Count))

gene_ratio_max <- max(c(KEGG_downregulated_basal_vs_60_plot$GeneRatio,
                        KEGG_downregulated_10_vs_60_plot$GeneRatio))

# combine all downregulated comparisons
KEGG_downregulated_combined <- bind_rows(KEGG_downregulated_basal_vs_60_plot %>%
                                         mutate(Comparison = "Basal vs 60 min ARS"),
                                         KEGG_downregulated_10_vs_60_plot %>%
                                         mutate(Comparison = "10 vs 60 min ARS"))

# order comparisons
KEGG_downregulated_combined$Comparison <- factor(KEGG_downregulated_combined$Comparison,
                                               levels = c("Basal vs 60 min ARS", "10 vs 60 min ARS"),
                                               labels = c("Basal vs 60 min", "10 vs 60 min"))

# plot all three downregulated dotplots into one plot
ggplot(KEGG_downregulated_combined,
       aes(x = GeneRatio, y = Description, 
           size = Count,
           colour = neglog10FDR)) +
  geom_point() +
  facet_grid(. ~ Comparison, scales = "fixed") +
  scale_colour_gradient(low = "#C9E5F5", high = "#2674B8") +
  scale_size_continuous(breaks = scales::breaks_width(1)) +
  scale_x_continuous(breaks = seq(0, 0.5, 0.1), limits = c(0, 0.5)) +
  labs(x = "GeneRatio", y = NULL, 
       size = "Count", 
       colour = expression(-log[10]("FDR"))) +
  theme_bw() +
  theme(panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold"),
        axis.text = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(3, "mm"),
        legend.key.height = unit(3, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "pt"),
        panel.spacing.x = unit(0.5, "lines"))

ggsave("AppendixC_KEGG_downregulated_combined.tiff", width = 6.27, height = 1.5, dpi = 300)
