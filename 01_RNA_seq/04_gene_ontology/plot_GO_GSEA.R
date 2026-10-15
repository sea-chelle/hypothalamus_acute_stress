# set working directory
setwd(here("plots"))

# select terms
GSEA_GO_select <- c("response to calcium ion", "response to hypoxia", "circadian regulation of gene expression",
                    "forebrain development", "regulation of synapse organization", "response to corticosteroid",
                    "glutamate receptor signaling pathway", "oxidative phosphorylation", "TOR signaling",
                    "intrinsic apoptotic signaling pathway", "neurotransmitter reuptake", "p38MAPK cascade",
                    "response to progesterone", "gonadotropin secretion", "heterochromatin organization")

# run for loop
for (comparison in GSEA_comparisons) {
  
  # grab GSEA result
  GSEA_object <- get(paste0("GSEA_GO_", comparison))
  
  # filter for significant terms in candidate list for each comparison
  GSEA_plot <- GSEA_object@result %>%
    filter(Description %in% GSEA_GO_select, p.adjust <= 0.05)
  
  # add to a new object
  assign(paste0("GSEA_GO_", comparison, "_plot"), GSEA_plot)
}

# check for loop worked
for (comparison in GSEA_comparisons) {
  cat("\n", comparison, "\n")
  print(get(paste0("GSEA_GO_", comparison, "_plot")) %>%
          select(Description, NES, p.adjust))
}

# add all 6 comparisons into one object
GSEA_GO_plot <- bind_rows(lapply(GSEA_comparisons, function(comparison) {
  
  # grab subsetted GSEA object and add comparison column
  get(paste0("GSEA_GO_", comparison, "_plot")) %>%
      mutate(Sex = ifelse(grepl("^male", comparison), "Male", "Female"),
        
        Comparison = case_when(grepl("basal_vs_10", comparison) ~ "Basal vs 10 min",
                               grepl("basal_vs_60", comparison) ~ "Basal vs 60 min",
                               grepl("10_vs_60", comparison) ~ "10 vs 60 min"))
  })) %>%
  
  # hshortern GO terms
  mutate(PlotLabel = recode(Description,
                            "response to calcium ion" = "Calcium ion response",
                            "response to hypoxia" = "Hypoxia response",
                            "circadian regulation of gene expression" = "Circadian regulation",
                            "forebrain development" = "Forebrain development",
                            "regulation of synapse organization" = "Synapse organisation",
                            "response to corticosteroid" = "Corticosteroid response",
                            "glutamate receptor signaling pathway" = "Glutamate receptor signalling",
                            "oxidative phosphorylation" = "Oxidative phosphorylation",
                            "TOR signaling" = "TOR signalling",
                            "intrinsic apoptotic signaling pathway" = "Apoptotic signalling",
                            "neurotransmitter reuptake" = "Neurotransmitter reuptake",
                            "p38MAPK cascade" = "p38 MAPK cascade",
                            "response to progesterone" = "Progesterone response",
                            "gonadotropin secretion" = "Gonadotropin secretion",
                            "heterochromatin organization" = "Heterochromatin organisation",
                            "Forebrain development" = "Forebrain development"),
    
    # order comparison, sex and plot labels
    Comparison = factor(Comparison, levels = c("Basal vs 10 min", "Basal vs 60 min", "10 vs 60 min")),
    Sex = factor(Sex, levels = c("Male","Female")),
    PlotLabel = factor(PlotLabel,
                       levels = rev(c("Corticosteroid response",
                                     "Progesterone response",
                                     "Gonadotropin secretion",
                                     "p38 MAPK cascade",
                                     "TOR signalling",
                                     "Oxidative phosphorylation",
                                     "Hypoxia response",
                                     "Apoptotic signalling",
                                     "Calcium ion response",
                                     "Neurotransmitter reuptake",
                                     "Glutamate receptor signalling",
                                     "Synapse organisation",
                                     "Circadian regulation",
                                     "Heterochromatin organisation",
                                     "Forebrain development"))))

ggplot(GSEA_GO_plot, aes(x = NES, y = PlotLabel, fill = -log10(p.adjust), shape = Sex)) +
  
  # zero reference line
  geom_vline(xintercept = 0,
             linetype = "dashed",
             linewidth = 0.4) +
  
  geom_point(size = 3, colour = "black", stroke = 0.3) +
  
  # make comparison headers span over 2 rows
  facet_wrap(~ Comparison, nrow = 1) +
  scale_shape_manual(values = c("Male" = 21, "Female" = 24)) +
  scale_fill_gradient(low = "#EEEAF7", high = "#7656A8") +
  scale_x_continuous(breaks = scales::breaks_width(1),
                     expand = expansion(add = c(0.5, 0.5))) +
  labs(x = "Normalised enrichment score (NES)", y = NULL, shape = "Sex", fill = expression(-log[10](FDR))) +
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

ggsave("fig4.6a_GSEA_GO_sex_dotplot.tiff", width = 6.27, height = 3.5 ,dpi = 600, compression = "lzw")

# select terms
GSEA_GO_select <- c("oxidative phosphorylation", "synapse assembly", "gonadotropin secretion",
                    "synaptic transmission, glutamatergic", "synaptic vesicle cycle", "circadian rhythm", "response to hypoxia",
                    "neurotransmitter transport")

# run for loop
for (comparison in GSEA_comparisons) {
  
  # grab GSEA result
  GSEA_object <- get(paste0("GSEA_GO_", comparison))
  
  # filter for significant terms in candidate list for each comparison
  GSEA_plot <- GSEA_object@result %>%
    filter(Description %in% GSEA_GO_select, p.adjust <= 0.05)
  
  # add to a new object
  assign(paste0("GSEA_GO_", comparison, "_plot"), GSEA_plot)
}
  
# check for loop worked
for (comparison in GSEA_comparisons) {
  cat("\n", comparison, "\n")
  print(get(paste0("GSEA_GO_", comparison, "_plot")) %>%
          select(Description, NES, p.adjust))
}

# add all 6 comparisons into one object
GSEA_GO_plot <- bind_rows(lapply(GSEA_comparisons, function(comparison) {
  
  # grab subsetted GSEA object and add comparison column
  get(paste0("GSEA_GO_", comparison, "_plot")) %>%
    mutate(Comparison = case_when(grepl("basal", comparison) ~ "Basal",
                                  grepl("10min", comparison) ~ "10 min",
                                  grepl("60min", comparison) ~ "60 min"))
})) %>%
  
  # hshortern GO terms
  mutate(PlotLabel = recode(Description,
                            "oxidative phosphorylation" = "Oxidative phosphorylation",
                            "synaptic transmission, glutamatergic" = "Glutamatergic transmission",
                            "synaptic vesicle cycle" = "Synaptic vesicle",
                            "circadian rhythm" = "Circadian rhythm",
                            "response to hypoxia" = "Hypoxia response",
                            "neurotransmitter transport" = "Neurotransmitter transport",
                            "synapse assembly" = "Synapse assembly",
                            "gonadotropin secretion" = "Gonadotropin secretion"),
         
         # order comparison, sex and plot labels
         Comparison = factor(Comparison, levels = c("Basal",
                                                    "10 min",
                                                    "60 min")),
         
         
         PlotLabel = factor(PlotLabel,
                            levels = rev(c("Gonadotropin secretion",
                                          "Oxidative phosphorylation",
                                          "Hypoxia response",
                                          "Neurotransmitter transport",
                                          "Glutamatergic transmission",
                                          "Synaptic vesicle",
                                          "Synapse assembly",
                                          "Circadian rhythm"))))

ggplot(GSEA_GO_plot, aes(x = NES, y = PlotLabel, colour = -log10(p.adjust))) +
  
  # zero reference line
  geom_vline(xintercept = 0,
             linetype = "dashed",
             linewidth = 0.4) +
  geom_point(size = 3) +
  
  # comparison headers
  facet_wrap(~ Comparison, nrow = 1) +
  scale_colour_gradient(low = "#E9DFFC", high = "#7656A8") +
  scale_x_continuous(breaks = scales::breaks_width(1), expand = expansion(add = c(0.5, 0.5))) +
  labs(x = "Normalised enrichment score (NES)", y = NULL, colour = expression(-log[10](FDR))) +
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
        panel.heights = unit(8 * 0.5, "cm"))

ggsave("GSEA_GO_intersex_dotplot.tiff", width = 6.27, height = 4, dpi = 300)

# select terms
GSEA_GO_select <- c("regulation of ERK1 and ERK2 cascade", "TOR signaling", "neurotransmitter secretion", "neurotransmitter reuptake",
                    "regulation of dopamine secretion", "regulation of synapse organization", "synaptic transmission, glutamatergic",
                    "translation at presynapse", "regulation of circadian rhythm", "response to hypoxia", "transforming growth factor beta receptor signaling pathway")
# run for loop
for (comparison in GSEA_comparisons) {
  
  # grab GSEA result
  GSEA_object <- get(paste0("GSEA_GO_", comparison))
  
  # filter for significant terms in candidate list for each comparison
  GSEA_plot <- GSEA_object@result %>%
    filter(Description %in% GSEA_GO_select, p.adjust <= 0.05)
  
  # add to a new object
  assign(paste0("GSEA_GO_", comparison, "_plot"), GSEA_plot)
}

# check for loop worked
for (comparison in GSEA_comparisons) {
  cat("\n", comparison, "\n")
  print(get(paste0("GSEA_GO_", comparison, "_plot")) %>%
          select(Description, NES, p.adjust))
}

# add all 6 comparisons into one object
GSEA_GO_interaction_plot <- bind_rows(lapply(GSEA_comparisons, function(comparison) {
  
  # grab subsetted GSEA object and add comparison column
  get(paste0("GSEA_GO_", comparison, "_plot")) %>%
    mutate(Comparison = case_when(grepl("_10", comparison) ~ "10 min",
                                  grepl("_60", comparison) ~ "60 min"))
})) %>%
  
  # hshortern GO terms
  mutate(PlotLabel = recode(Description,
                            "translation at presynapse" = "Translation at presynapse",
                            "synaptic transmission, glutamatergic" = "Glutamatergic transmission",
                            "regulation of synapse organization" = "Synapse organisation",
                            "regulation of circadian rhythm" = "Circadian rhythm",
                            "response to hypoxia" = "Hypoxia response",
                            "regulation of dopamine secretion" = "Dopamine secretion",
                            "TOR signaling" = "TOR signalling",
                            "regulation of ERK1 and ERK2 cascade" = "ERK1/2 cascade",
                            "neurotransmitter reuptake" = "Neurotransmitter reuptake",
                            "neurotransmitter secretion" = "Neurotransmitter secretion",
                            "transforming growth factor beta receptor signaling pathway" = "TGFβ receptor signalling"),
         
         # order comparison and plot labels
         Comparison = factor(Comparison, levels = c("10 min",
                                                    "60 min")),
         
         
         PlotLabel = factor(PlotLabel,
                            levels = rev(c("ERK1/2 cascade",
                                           "TOR signalling",
                                           "TGFβ receptor signalling",
                                           "Hypoxia response",
                                           "Translation at presynapse",
                                           "Synapse organisation",
                                           "Neurotransmitter secretion",
                                           "Neurotransmitter reuptake",
                                           "Glutamatergic transmission",
                                           "Dopamine secretion",
                                           "Circadian rhythm"))))

# plot
ggplot(GSEA_GO_interaction_plot, aes(x = NES, y = PlotLabel, colour = -log10(p.adjust))) +
  
  # zero reference line
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(size = 3.0) +
  
  # interaction comparison headers
  facet_wrap(~ Comparison, nrow = 1,
             labeller = as_labeller(c("10 min" = "10 min interaction",
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

ggsave("GSEA_GO_interaction_dotplot.tiff", width = 6.27, height = 2.9, dpi = 600, compression = "lzw")
