# load required packages
library(here)
library(dplyr)
library(ggplot2)

# set working directory
setwd(here("plots"))

# create vector of selected terms to plot 
GO_basal_vs_10_select <- c("response to corticosteroid", "response to interleukin-1", 
                           "response to reactive oxygen species", "regulation of synaptic plasticity",
                           "p38MAPK cascade")

# filter for selected terms to plot and calculate GeneRatio and -log10 FDR
GO_basal_vs_10_plot <- GO_basal_vs_10@result %>%
  filter(Description %in% GO_basal_vs_10_select) %>%
  mutate(GeneRatio = sapply(GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
    }),
    neglog10FDR = -log10(p.adjust),
    PlotLabel = case_when(Description == "response to corticosteroid" ~ "Corticosteroid response",
                          Description == "response to interleukin-1" ~ "Interleukin-1 response",
                          Description == "response to reactive oxygen species" ~ "ROS response",
                          Description == "regulation of synaptic plasticity" ~ "Synaptic plasticity",
                          Description == "p38MAPK cascade" ~ "MAPK cascade",
                          Description == "response to hypoxia" ~ "Hypoxia response",
                          Description == "response to estradiol" ~ "Estradiol response",
                          Description == "Notch signaling pathway" ~ "Notch signalling pathway",
                          Description == "stem cell differentiation" ~ "Stem cell differentiation",
                          Description == "intrinsic apoptotic signaling pathway" ~ "Apoptotic signalling",
                          Description == "lymphocyte differentiation" ~ "Lymphocyte differentiation",
                          Description == "canonical NF-kappaB signal transduction" ~ "NF-κB signalling",
                          TRUE ~ Description))

# create vector of selected terms to plot 
GO_upregulated_basal_vs_60_select <- c("response to corticosteroid", "p38MAPK cascade",
                                       "response to reactive oxygen species", "response to interleukin-1",
                                       "response to hypoxia", "response to estradiol")

# filter for selected terms to plot and calculate GeneRatio and -log10 FDR
GO_upregulated_basal_vs_60_plot <- GO_upregulated_basal_vs_60@result %>%
  filter(Description %in% GO_upregulated_basal_vs_60_select) %>%
  mutate(GeneRatio = sapply(GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  }),
  neglog10FDR = -log10(p.adjust),
  PlotLabel = case_when(Description == "response to corticosteroid" ~ "Corticosteroid response",
                        Description == "response to interleukin-1" ~ "Interleukin-1 response",
                        Description == "response to reactive oxygen species" ~ "ROS response",
                        Description == "regulation of synaptic plasticity" ~ "Synaptic plasticity",
                        Description == "p38MAPK cascade" ~ "MAPK cascade",
                        Description == "response to hypoxia" ~ "Hypoxia response",
                        Description == "response to estradiol" ~ "Estradiol response",
                        Description == "Notch signaling pathway" ~ "Notch signalling pathway",
                        Description == "stem cell differentiation" ~ "Stem cell differentiation",
                        Description == "intrinsic apoptotic signaling pathway" ~ "Apoptotic signalling",
                        Description == "lymphocyte differentiation" ~ "Lymphocyte differentiation",
                        Description == "canonical NF-kappaB signal transduction" ~ "NF-κB signalling",
                        TRUE ~ Description))

# create vector of selected terms to plot 
GO_upregulated_10_vs_60_select <- c("response to corticosteroid", "p38MAPK cascade", "response to hypoxia",
                                    "response to reactive oxygen species", "response to interleukin-1",
                                    "intrinsic apoptotic signaling pathway")

# filter for selected terms to plot and calculate GeneRatio and -log10 FDR
GO_upregulated_10_vs_60_plot <- GO_upregulated_10_vs_60@result %>%
  filter(Description %in% GO_upregulated_10_vs_60_select) %>%
  mutate(GeneRatio = sapply(GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  }),
  neglog10FDR = -log10(p.adjust),
  PlotLabel = case_when(Description == "response to corticosteroid" ~ "Corticosteroid response",
                        Description == "response to interleukin-1" ~ "Interleukin-1 response",
                        Description == "response to reactive oxygen species" ~ "ROS response",
                        Description == "regulation of synaptic plasticity" ~ "Synaptic plasticity",
                        Description == "p38MAPK cascade" ~ "MAPK cascade",
                        Description == "response to hypoxia" ~ "Hypoxia response",
                        Description == "response to estradiol" ~ "Estradiol response",
                        Description == "Notch signaling pathway" ~ "Notch signalling pathway",
                        Description == "stem cell differentiation" ~ "Stem cell differentiation",
                        Description == "intrinsic apoptotic signaling pathway" ~ "Apoptotic signalling",
                        Description == "lymphocyte differentiation" ~ "Lymphocyte differentiation",
                        Description == "canonical NF-kappaB signal transduction" ~ "NF-κB signalling",
                        TRUE ~ Description))

# make FDR, gene count and gene ratio scales the same across all upregulated plots
fdr_limits <- range(c(GO_basal_vs_10_plot$neglog10FDR,
                      GO_upregulated_basal_vs_60_plot$neglog10FDR,
                      GO_upregulated_10_vs_60_plot$neglog10FDR))

count_limits <- range(c(GO_basal_vs_10_plot$Count,
                        GO_upregulated_basal_vs_60_plot$Count,
                        GO_upregulated_10_vs_60_plot$Count))

gene_ratio_max <- max(c(GO_basal_vs_10_plot$GeneRatio,
                        GO_upregulated_basal_vs_60_plot$GeneRatio,
                        GO_upregulated_10_vs_60_plot$GeneRatio))

# combine all upregulated comparisons
GO_upregulated_combined <- bind_rows(GO_basal_vs_10_plot %>%
                                       mutate(Comparison = "Basal vs 10 min ARS"),
                                     GO_upregulated_basal_vs_60_plot %>%
                                       mutate(Comparison = "Basal vs 60 min ARS"),
                                     GO_upregulated_10_vs_60_plot %>%
                                       mutate(Comparison = "10 vs 60 min ARS"))

# order terms
GO_upregulated_combined$PlotLabel <- factor(
  GO_upregulated_combined$PlotLabel,
  levels = rev(c("Corticosteroid response", "MAPK cascade", "ROS response", "Interleukin-1 response",
                 "Synaptic plasticity", "Hypoxia response", "Estradiol response", "Apoptotic signalling")))

# order comparisons
GO_upregulated_combined$Comparison <- factor(GO_upregulated_combined$Comparison,
                                             levels = c("Basal vs 10 min ARS", "Basal vs 60 min ARS", "10 vs 60 min ARS"),
                                             labels = c("Basal vs 10 min", "Basal vs 60 min", "10 vs 60 min"))

# plot all three upregulated dotplots into one plot
ggplot(GO_upregulated_combined, aes(x = GeneRatio, y = PlotLabel, size = Count, colour = neglog10FDR)) +
  geom_point() + facet_grid(. ~ Comparison, scales = "fixed") +
  scale_colour_gradient(low = "#FADADA", high = "#D9364A") +
  scale_size_continuous(breaks = scales::breaks_width(1)) +
  scale_x_continuous(breaks = seq(0, 0.3, 0.1), limits = c(0, 0.35)) +
  labs(x = "GeneRatio", y = NULL, size = "Count", colour = expression(-log[10]("FDR"))) +
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

ggsave("fig4.3._GO_upregulated_combined.tiff", width = 6.27, height = 2.8, dpi = 300)

# create vector of selected terms to plot 
GO_downregulated_basal_vs_60_select <- c("Notch signaling pathway", "stem cell differentiation")

# filter for selected terms to plot and calculate GeneRatio and -log10 FDR
GO_downregulated_basal_vs_60_plot <- GO_downregulated_basal_vs_60@result %>%
  filter(Description %in% GO_downregulated_basal_vs_60_select) %>%
  mutate(GeneRatio = sapply(GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  }),
  neglog10FDR = -log10(p.adjust),
  PlotLabel = case_when(Description == "response to corticosteroid" ~ "Corticosteroid response",
                        Description == "response to interleukin-1" ~ "Interleukin-1 response",
                        Description == "response to reactive oxygen species" ~ "ROS response",
                        Description == "regulation of synaptic plasticity" ~ "Synaptic plasticity",
                        Description == "p38MAPK cascade" ~ "MAPK cascade",
                        Description == "response to hypoxia" ~ "Hypoxia response",
                        Description == "response to estradiol" ~ "Estradiol response",
                        Description == "Notch signaling pathway" ~ "Notch signalling pathway",
                        Description == "stem cell differentiation" ~ "Stem cell differentiation",
                        Description == "intrinsic apoptotic signaling pathway" ~ "Apoptotic signalling",
                        Description == "lymphocyte differentiation" ~ "Lymphocyte differentiation",
                        Description == "canonical NF-kappaB signal transduction" ~ "NF-κB signalling",
                        TRUE ~ Description))

# create vector of selected terms to plot 
GO_downregulated_10_vs_60_select <- c("stem cell differentiation", "lymphocyte differentiation",
                                      "Notch signaling pathway", "canonical NF-kappaB signal transduction")


# filter for selected terms to plot and calculate GeneRatio and -log10 FDR
GO_downregulated_10_vs_60_plot <- GO_downregulated_10_vs_60@result %>%
  filter(Description %in% GO_downregulated_10_vs_60_select) %>%
  mutate(GeneRatio = sapply(GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  }),
  neglog10FDR = -log10(p.adjust),
  PlotLabel = case_when(Description == "response to corticosteroid" ~ "Corticosteroid response",
                        Description == "response to interleukin-1" ~ "Interleukin-1 response",
                        Description == "response to reactive oxygen species" ~ "ROS response",
                        Description == "regulation of synaptic plasticity" ~ "Synaptic plasticity",
                        Description == "p38MAPK cascade" ~ "MAPK cascade",
                        Description == "response to hypoxia" ~ "Hypoxia response",
                        Description == "response to estradiol" ~ "Estradiol response",
                        Description == "Notch signaling pathway" ~ "Notch signalling pathway",
                        Description == "stem cell differentiation" ~ "Stem cell differentiation",
                        Description == "intrinsic apoptotic signaling pathway" ~ "Apoptotic signalling",
                        Description == "lymphocyte differentiation" ~ "Lymphocyte differentiation",
                        Description == "canonical NF-kappaB signal transduction" ~ "NF-κB signalling",
                        TRUE ~ Description))


# make FDR, gene count and gene ratio scales the same across all downregulated plots
fdr_limits <- range(c(GO_downregulated_basal_vs_60_plot$neglog10FDR,
                      GO_downregulated_10_vs_60_plot$neglog10FDR))

count_limits <- range(c(GO_downregulated_basal_vs_60_plot$Count,
                        GO_downregulated_10_vs_60_plot$Count))

gene_ratio_max <- max(c(GO_downregulated_basal_vs_60_plot$GeneRatio,
                        GO_downregulated_10_vs_60_plot$GeneRatio))

# combine all downregulated comparisons
GO_downregulated_combined <- bind_rows(GO_downregulated_basal_vs_60_plot %>%
                                       mutate(Comparison = "Basal vs 60 min ARS"),
                                       GO_downregulated_10_vs_60_plot %>%
                                       mutate(Comparison = "10 vs 60 min ARS"))

# order terms
GO_downregulated_combined$PlotLabel <- factor(
  GO_downregulated_combined$PlotLabel,
  levels = rev(c("Stem cell differentiation", "Notch signalling pathway",
                 "Lymphocyte differentiation", "NF-κB signalling")))

# order comparisons
GO_downregulated_combined$Comparison <- factor(GO_downregulated_combined$Comparison, 
                                               levels = c("Basal vs 60 min ARS", "10 vs 60 min ARS"),
                                               labels = c("Basal vs 60 min", "10 vs 60 min"))

# plot all three downregulated dotplots into one plot
ggplot(GO_downregulated_combined, aes(x = GeneRatio, y = PlotLabel,size = Count, colour = neglog10FDR)) +
  geom_point() +
  facet_grid(. ~ Comparison, scales = "fixed") +
  scale_colour_gradient(low = "#C9E5F5", high = "#2674B8") +
  scale_size_continuous(breaks = scales::breaks_width(1)) +
  scale_x_continuous(breaks = seq(0, 0.5, 0.1), limits = c(0, 0.4)) +
  labs(x = "GeneRatio", y = NULL, size = "Count", colour = expression(-log[10]("FDR"))) +
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

ggsave("fig.C.1_GO_downregulated_combined.tiff", width = 6.27, height = 2, dpi = 300)
