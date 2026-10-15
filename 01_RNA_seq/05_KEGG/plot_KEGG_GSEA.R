# load required packages
library(here)
library(dplyr)
library(ggplot2)

# set working dierectory
setwd(here("plots"))

# select terms
GSEA_KEGG_select <- c("Oxytocin signaling pathway", "Growth hormone synthesis, secretion and action", "FoxO signaling pathway",
                      "Notch signaling pathway", "MAPK signaling pathway", "mTOR signaling pathway",
                      "Apoptosis", "IL-17 signaling pathway", "Integrated stress response (ISR) signaling pathway",
                      "Oxidative phosphorylation", "Dopaminergic synapse", "Serotonergic synapse",
                      "Circadian entrainment", "Glutamatergic synapse", "Synaptic vesicle cycle")

# run for loop
for (comparison in GSEA_comparisons) {
  
  # grab GSEA result
  GSEA_object <- get(paste0("GSEA_KEGG_", comparison))
  
  # filter for significant terms in candidate list for each comparison
  GSEA_plot <- GSEA_object@result %>%
    filter(Description %in% GSEA_KEGG_select, p.adjust <= 0.05)
  
  # add to a new object
  assign(paste0("GSEA_KEGG_", comparison, "_plot"), GSEA_plot)
}

# check for loop worked
for (comparison in GSEA_comparisons) {
  cat("\n", comparison, "\n")
  print(get(paste0("GSEA_KEGG_", comparison, "_plot")) %>%
          select(Description, NES, p.adjust))
}

# add all 6 comparisons into one object
GSEA_KEGG_plot <- bind_rows(lapply(GSEA_comparisons, function(comparison) {
  
  # grab subsetted GSEA object and add comparison column
  get(paste0("GSEA_KEGG_", comparison, "_plot")) %>%
    mutate(Sex = ifelse(grepl("^male", comparison), "Male", "Female"),
           Comparison = case_when(grepl("basal_vs_10", comparison) ~ "Basal vs 10 min",
                                  grepl("basal_vs_60", comparison) ~ "Basal vs 60 min",
                                  grepl("10_vs_60", comparison) ~ "10 vs 60 min"))
  })) %>%
  
  # shortern KEGG terms
  mutate(PlotLabel = recode(Description,
                            "Oxytocin signaling pathway" = "Oxytocin signalling",
                            "Growth hormone synthesis, secretion and action" = "GH signalling",
                            "FoxO signaling pathway" = "FoxO signalling",
                            "Notch signaling pathway" = "Notch signalling",
                            "MAPK signaling pathway" = "MAPK signalling",
                            "mTOR signaling pathway" = "mTOR signalling",
                            "Apoptosis" = "Apoptosis",
                            "IL-17 signaling pathway" = "IL-17 signalling",
                            "Integrated stress response (ISR) signaling pathway" = "ISR signalling",
                            "Oxidative phosphorylation" = "Oxidative phosphorylation",
                            "Dopaminergic synapse" = "Dopaminergic synapse",
                            "Serotonergic synapse" = "Serotonergic synapse",
                            "Circadian entrainment" = "Circadian entrainment",
                            "Glutamatergic synapse" = "Glutamatergic synapse",
                            "Synaptic vesicle cycle" = "Synaptic vesicle cycle"),
         
         # order comparison, sex and plot labels
         Comparison = factor(Comparison, levels = c("Basal vs 10 min", 
                                                    "Basal vs 60 min", 
                                                    "10 vs 60 min")),
         
         Sex = factor(Sex, levels = c("Male",
                                      "Female")),
         
         PlotLabel = factor(PlotLabel,
                            levels = rev(c("Oxytocin signalling",
                                          "GH signalling",
                                          "Notch signalling",
                                          "IL-17 signalling",
                                          "mTOR signalling",
                                          "MAPK signalling",
                                          "FoxO signalling",
                                          "ISR signalling",
                                          "Apoptosis",
                                          "Oxidative phosphorylation",
                                          "Glutamatergic synapse",
                                          "Dopaminergic synapse",
                                          "Serotonergic synapse",
                                          "Synaptic vesicle cycle",
                                          "Circadian entrainment"))))

ggplot(GSEA_KEGG_plot, aes(x = NES, y = PlotLabel, fill = -log10(p.adjust), shape = Sex)) +
  
  # zero reference line
  geom_vline(xintercept = 0,
             linetype = "dashed",
             linewidth = 0.4) +
  
  geom_point(size = 3, colour = "black", stroke = 0.3) +
  
  # make comparison headers span over 2 rows
  facet_wrap(~ Comparison, nrow = 1,) +
  scale_shape_manual(values = c("Male" = 21, "Female" = 24)) +
  scale_fill_gradient(low = "#EEEAF7", high = "#7656A8") +
  scale_x_continuous(breaks = scales::breaks_width(1),
                     expand = expansion(add = c(0.5, 0.5))) +
  labs(x = "Normalized enrichment score (NES)", y = NULL, shape = "Sex", fill = expression(-log[10](FDR))) +
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
        panel.spacing.x = unit(0.5, "lines"),
        panel.widths = unit(rep(3, 3), "cm"),
        panel.heights = unit(15 * 0.5, "cm"))

ggsave("fig4.7a_GSEA_KEGG_male_and_female_dotplot.tiff", width = 6.27, height = 3.5, dpi = 600, compression = "lzw")


# select terms
GSEA_KEGG_select <- c("Oxidative phosphorylation", "Estrogen signaling pathway", "Glutamatergic synapse",
                      "Growth hormone synthesis, secretion and action", "Cushing syndrome", "Cholinergic synapse",
                      "Dopaminergic synapse", "mTOR signaling pathway")

# run for loop
for (comparison in GSEA_comparisons) {
  
  # grab GSEA result
  GSEA_object <- get(paste0("GSEA_KEGG_", comparison))
  
  # filter for significant terms in candidate list for each comparison
  GSEA_plot <- GSEA_object@result %>%
    filter(Description %in% GSEA_KEGG_select, p.adjust <= 0.05)
  
  # add to a new object
  assign(paste0("GSEA_KEGG_", comparison, "_plot"), GSEA_plot)
}

# check for loop worked
for (comparison in GSEA_comparisons) {
  cat("\n", comparison, "\n")
  print(get(paste0("GSEA_KEGG_", comparison, "_plot")) %>%
          select(Description, NES, p.adjust))
}

# add all 6 comparisons into one object
GSEA_KEGG_plot <- bind_rows(lapply(GSEA_comparisons, function(comparison) {
  
  # grab subsetted GSEA object and add comparison column
  get(paste0("GSEA_KEGG_", comparison, "_plot")) %>%
    mutate(Comparison = case_when(grepl("basal", comparison) ~ "Basal",
                                  grepl("10min", comparison) ~ "10 min",
                                  grepl("60min", comparison) ~ "60 min"))
})) %>%
  
  # hshortern KEGG terms
  mutate(PlotLabel = recode(Description,
                            "Oxidative phosphorylation" = "Oxidative phosphorylation",
                            "Glutamatergic synapse" = "Glutamatergic synapse",
                            "Cholinergic synapse" = "Cholinergic synapse",
                            "mTOR signaling pathway" = "mTOR signalling",
                            "Estrogen signaling pathway" = "Estrogen signalling",
                            "Growth hormone synthesis, secretion and action" = "GH signalling",
                            "Cushing syndrome" = "Cushing syndrome"),
         
         # order comparison, sex and plot labels
         Comparison = factor(Comparison, levels = c("Basal",
                                                    "10 min",
                                                    "60 min")),
         PlotLabel = factor(PlotLabel,
                            levels = rev(c("GH signalling",
                                           "Estrogen signalling",
                                           "Cushing syndrome",
                                           "mTOR signalling",
                                           "Oxidative phosphorylation",
                                           "Glutamatergic synapse",
                                           "Cholinergic synapse"))))

ggplot(GSEA_KEGG_plot, aes(x = NES, y = PlotLabel, colour = -log10(p.adjust))) +
  
  # zero reference line
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(size = 3) +
  
  # comparison headers
  facet_wrap(~ Comparison, nrow = 1) +
  scale_colour_gradient(low = "#E9DFFC", high = "#7656A8") +
  scale_x_continuous(breaks = scales::breaks_width(1), expand = expansion(add = c(0.5, 0.5))) +
  labs(x = "Normalized enrichment score (NES)", y = NULL, colour = expression(-log[10](FDR))) +
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
        panel.spacing.x = unit(0.5, "lines"),
        panel.widths = unit(rep(3, 3), "cm"),
        panel.heights = unit(6 * 0.5, "cm"))

ggsave("C.6_GSEA_KEGG_male_vs_female_dotplot.tiff", width = 6.27, height = 1.8, dpi = 600, compression = "lzw")


# select terms
GSEA_KEGG_select <- c("Circadian entrainment", "Glutamatergic synapse", 
                      "GnRH secretion", "mTOR signaling pathway",
                      "Growth hormone synthesis, secretion and action", "Dopaminergic synapse", "Oxytocin signaling pathway",
                      "Thyroid hormone signaling pathway", "MAPK signaling pathway", "Ribosome",
                      "Neurotrophin signaling pathway")
# run for loop
for (comparison in GSEA_comparisons) {
  
  # grab GSEA result
  GSEA_object <- get(paste0("GSEA_KEGG_", comparison))
  
  # filter for significant terms in candidate list for each comparison
  GSEA_plot <- GSEA_object@result %>%
    filter(Description %in% GSEA_KEGG_select, p.adjust <= 0.05)
  
  # add to a new object
  assign(paste0("GSEA_KEGG_", comparison, "_plot"), GSEA_plot)
}

# check for loop worked
for (comparison in GSEA_comparisons) {
  cat("\n", comparison, "\n")
  print(get(paste0("GSEA_KEGG_", comparison, "_plot")) %>%
          select(Description, NES, p.adjust))
}

# add all 6 comparisons into one object
GSEA_KEGG_interaction_plot <- bind_rows(lapply(GSEA_comparisons, function(comparison) {
  
  # grab subsetted GSEA object and add comparison column
  get(paste0("GSEA_KEGG_", comparison, "_plot")) %>%
    mutate(Comparison = ifelse(grepl("_10", comparison),
                               "10 min",
                               "60 min"))
})) %>%
  
  # hshortern KEGG terms
  mutate(NES = -NES, # make males positive and females negative 
         
         PlotLabel = recode(Description,
                            "Oxytocin signaling pathway" = "Oxytocin signalling",
                            "Growth hormone synthesis, secretion and action" = "GH signalling",
                            "GnRH secretion" = "GnRH secretion",
                            "Thyroid hormone signaling pathway" = "TH signalling",
                            "mTOR signaling pathway" = "mTOR signalling",
                            "MAPK signaling pathway" = "MAPK signalling",
                            "Neurotrophin signaling pathway" = "Neurotrophin signalling",
                            "Dopaminergic synapse" = "Dopaminergic synapse",
                            "Circadian entrainment" = "Circadian entrainment",
                            "Glutamatergic synapse" = "Glutamatergic synapse",
                            "Ribosome" = "Ribosome"),
         
         # order comparison and plot labels
         Comparison = factor(Comparison, levels = c("10 min",
                                                    "60 min")),
         
         
         PlotLabel = factor(PlotLabel,
                            levels = rev(c("Oxytocin signalling",
                                           "TH signalling",
                                          "GH signalling",
                                          "GnRH secretion",
                                          "mTOR signalling",
                                          "MAPK signalling",
                                          "Neurotrophin signalling",
                                          "Glutamatergic synapse",
                                          "Dopaminergic synapse",
                                          "Circadian entrainment",
                                          "Ribosome"))))

ggplot(GSEA_KEGG_interaction_plot, aes(x = NES, y = PlotLabel, colour = -log10(p.adjust))) +
  
  # zero reference line
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(size = 3) +
  
  # interaction comparison headers
  facet_wrap(~ Comparison, nrow = 1, labeller = as_labeller(c("10 min" = "10 min interaction",
                                                              "60 min" = "60 min interaction"))) +
  scale_colour_gradient(low = "#E9DFFC", high = "#7656A8") +
  scale_x_continuous(breaks = scales::breaks_width(1), expand = expansion(add = c(0.5, 0.5))) +
  labs(x = "Normalized enrichment score (NES)\nFemale-shifted response \u2190        \u2192 Male-shifted response",
       y = NULL,
       colour = expression(-log[10](FDR))) +
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
        panel.spacing.x = unit(0.5, "lines"),
        panel.widths = unit(rep(3, 3), "cm"),
        panel.heights = unit(11 * 0.5, "cm"))

ggsave("fig4.7b_GSEA_KEGG_interaction_dotplot.tiff", width = 4.5, height = 2.9, dpi = 600, compression = "lzw")
