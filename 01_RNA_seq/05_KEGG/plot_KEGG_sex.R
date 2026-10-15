# load required packages
library(here)
library(dplyr)
library(ggplot2)

# set working directory
setwd(here("plots"))

# selected terms to plot 
KEGG_male_basal_vs_10_select <- c("Growth hormone synthesis, secretion and action", "TNF signaling pathway", "IL-17 signaling pathway")
KEGG_male_basal_vs_10_plot <- KEGG_male_basal_vs_10@result
KEGG_male_basal_vs_10_plot <- KEGG_male_basal_vs_10_plot[KEGG_male_basal_vs_10_plot$Description %in% KEGG_male_basal_vs_10_select, ]

# selected terms to plot 
KEGG_female_basal_vs_10_select <- c("Growth hormone synthesis, secretion and action", "TNF signaling pathway", "MAPK signaling pathway")
KEGG_female_basal_vs_10_plot <- KEGG_female_basal_vs_10@result
KEGG_female_basal_vs_10_plot <- KEGG_female_basal_vs_10_plot[KEGG_female_basal_vs_10_plot$Description %in% KEGG_female_basal_vs_10_select, ]

# selected terms to plot 
KEGG_upregulated_male_basal_vs_60_select <- c("IL-17 signaling pathway", "TNF signaling pathway", "Apoptosis", "FoxO signaling pathway", "MAPK signaling pathway",
                                              "NF-kappa B signaling pathway")
KEGG_upregulated_male_basal_vs_60_plot <- KEGG_upregulated_male_basal_vs_60@result
KEGG_upregulated_male_basal_vs_60_plot <- KEGG_upregulated_male_basal_vs_60_plot[KEGG_upregulated_male_basal_vs_60_plot$Description %in% KEGG_upregulated_male_basal_vs_60_select, ]

# selected terms to plot 
KEGG_upregulated_female_basal_vs_60_select <- c("IL-17 signaling pathway", "TNF signaling pathway", "Apoptosis", "FoxO signaling pathway", "Circadian entrainment")
KEGG_upregulated_female_basal_vs_60_plot <- KEGG_upregulated_female_basal_vs_60@result
KEGG_upregulated_female_basal_vs_60_plot <- KEGG_upregulated_female_basal_vs_60_plot[KEGG_upregulated_female_basal_vs_60_plot$Description %in% KEGG_upregulated_female_basal_vs_60_select, ]

# selected terms to plot 
KEGG_upregulated_male_10_vs_60_select <- c("FoxO signaling pathway", "NF-kappa B signaling pathway", "Apoptosis",
                                           "Integrated stress response (ISR) signaling pathway", "Synaptic vesicle cycle")
KEGG_upregulated_male_10_vs_60_plot <- KEGG_upregulated_male_10_vs_60@result %>% 
  filter(Description %in% KEGG_upregulated_male_10_vs_60_select)

# selected terms to plot 
KEGG_upregulated_female_10_vs_60_select <- c()
KEGG_upregulated_female_10_vs_60_plot <- KEGG_upregulated_female_10_vs_60@result
KEGG_upregulated_female_10_vs_60_plot <- KEGG_upregulated_female_10_vs_60_plot[KEGG_upregulated_female_10_vs_60_plot$Description %in% KEGG_upregulated_female_10_vs_60_select, ]

# combine all upregulated KEGG analyses into 1 df
KEGG_sex_up <- bind_rows(KEGG_male_basal_vs_10_plot %>%
                           mutate(Sex = "Male", Comparison = "Basal vs 10 min"),
                         KEGG_female_basal_vs_10_plot %>%
                           mutate(Sex = "Female", Comparison = "Basal vs 10 min"),
                         KEGG_upregulated_male_basal_vs_60_plot %>%
                           mutate(Sex = "Male", Comparison = "Basal vs 60 min"),
                         KEGG_upregulated_female_basal_vs_60_plot %>%
                           mutate(Sex = "Female", Comparison = "Basal vs 60 min"),
                         KEGG_upregulated_male_10_vs_60_plot %>%
                           mutate(Sex = "Male", Comparison = "10 vs 60 min"),
                         KEGG_upregulated_female_10_vs_60_plot %>%
                           mutate(Sex = "Female", Comparison = "10 vs 60 min")) %>%
  
  # calculate gene ratio and FDR 
  mutate(GeneRatio = sapply(strsplit(GeneRatio, "/"),
                            function(x) as.numeric(x[1]) / as.numeric(x[2])),
         FDR = -log10(p.adjust),
         
         # shortern terms
         PlotLabels = case_when(Description == "Growth hormone synthesis, secretion and action" ~ "GH signalling",
                                Description == "TNF signaling pathway" ~ "TNF signalling",
                                Description == "IL-17 signaling pathway" ~ "IL-17 signalling",
                                Description == "MAPK signaling pathway" ~ "MAPK signalling",
                                Description == "Apoptosis" ~ "Apoptosis",
                                Description == "FoxO signaling pathway" ~ "FoxO signalling",
                                Description == "Synaptic vesicle cycle" ~ "Synaptic vesicle",
                                Description == "Integrated stress response (ISR) signaling pathway" ~ "ISR signalling",
                                Description == "NF-kappa B signaling pathway" ~ "NF-κB signalling",
                                Description == "Circadian entrainment" ~ "Circadian entrainment",
                                TRUE ~ Description),
         
         # set order for sex, comparisons and plot labels
         Sex = factor(Sex, levels = c("Male", 
                                      "Female")),
         
         Comparison = factor(Comparison, levels = c("Basal vs 10 min", 
                                                    "Basal vs 60 min", 
                                                    "10 vs 60 min")),
         
         PlotLabels = factor(PlotLabels, levels = rev(c("GH signalling",
                                                        "TNF signalling",
                                                        "IL-17 signalling",
                                                        "NF-κB signalling",
                                                        "MAPK signalling",
                                                        "FoxO signalling",
                                                        "ISR signalling",
                                                        "Apoptosis",
                                                        "Synaptic vesicle",
                                                        "Circadian entrainment"))))
# plot all upregulated comparisons as a dotplot
ggplot(KEGG_sex_up, aes(x = Count, y = PlotLabels, fill = FDR, shape = Sex)) +
  geom_point(size = 4, colour = "black", stroke = 0.3) +
  
  # make comparison headers span over 2 rows
  facet_wrap(~ Comparison, nrow = 1, 
             labeller = as_labeller(c("Basal vs 10 min" = "Basal vs\n10 min ARS",
                                      "Basal vs 60 min" = "Basal vs\n60 min ARS",
                                      "10 vs 60 min" = "10 vs 60 min\nARS"))) +
  

  scale_shape_manual(values = c("Male" = 21, "Female" = 24)) +
  scale_fill_gradient(low = "#FADADA", high = "#D9364A") +
  scale_x_continuous(breaks = scales::breaks_width(1), expand = expansion(add = c(0.8, 0.8))) +
  labs(x = "Gene count", y = NULL, shape = "Sex", fill = expression(-log[10](FDR))) +
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

ggsave("4.8_KEGG_male_and_female_dotplot.tiff", width = 6.27, height = 3.2, dpi = 300)
