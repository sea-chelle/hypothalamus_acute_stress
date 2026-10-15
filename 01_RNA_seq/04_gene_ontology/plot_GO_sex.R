# load required packages
library(here)
library(dplyr)
library(ggplot2)

# set working directory
setwd(here("plots"))

# create vector of selected terms to plot 
GO_male_basal_vs_10_select <- c("response to corticosteroid", "response to reactive oxygen species", "response to calcium ion",
                                "regulation of synaptic plasticity", "response to interleukin-1", "response to progesterone")

# create vector of selected terms to plot 
GO_female_basal_vs_10_select <- c("response to corticosteroid", "response to reactive oxygen species", "response to calcium ion")

# copy enrichResult object to new plot variable
GO_male_basal_vs_10_plot <- GO_male_basal_vs_10

# subset the result slot based on selected terms
GO_male_basal_vs_10_plot@result <- GO_male_basal_vs_10@result[GO_male_basal_vs_10@result$Description %in% GO_male_basal_vs_10_select,]

# copy enrichResult object to new plot variable
GO_female_basal_vs_10_plot <- GO_female_basal_vs_10

# subset the result slot based on selected terms
GO_female_basal_vs_10_plot@result <- GO_female_basal_vs_10@result[GO_female_basal_vs_10@result$Description %in% GO_female_basal_vs_10_select,]

# create vector of selected terms to plot 
GO_upregulated_male_basal_vs_60_select <- c("p38MAPK cascade", "response to corticosteroid", "response to reactive oxygen species",
                                            "response to progesterone", "response to hypoxia", "response to calcium ion",
                                            "intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator",
                                            "circadian regulation of gene expression", "interleukin-1-mediated signaling pathway", "response to estradiol",
                                            "neurotransmitter reuptake")

# create vector of selected terms to plot 
GO_upregulated_female_basal_vs_60_select <- c("p38MAPK cascade", "response to corticosteroid", "response to reactive oxygen species",
                                              "response to progesterone", "response to hypoxia", "response to calcium ion",
                                              "intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator",
                                              "circadian regulation of gene expression", "interleukin-1-mediated signaling pathway", "regulation of TOR signaling")

# copy enrichResult object to new plot variable
GO_upregulated_male_basal_vs_60_plot <- GO_upregulated_male_basal_vs_60

# subset the result slot based on selected terms
GO_upregulated_male_basal_vs_60_plot@result <- GO_upregulated_male_basal_vs_60@result[GO_upregulated_male_basal_vs_60@result$Description %in% GO_upregulated_male_basal_vs_60_select,]

# copy enrichResult object to new plot variable
GO_upregulated_female_basal_vs_60_plot <- GO_upregulated_female_basal_vs_60

# subset the result slot based on selected terms
GO_upregulated_female_basal_vs_60_plot@result <- GO_upregulated_female_basal_vs_60@result[GO_upregulated_female_basal_vs_60@result$Description %in% GO_upregulated_female_basal_vs_60_select,]

# create vector of selected terms to plot 
GO_upregulated_male_10_vs_60_select <- c("response to corticosteroid", "response to estradiol", "p38MAPK cascade", 
                                         "intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator",
                                         "neurotransmitter reuptake", "response to progesterone", "response to hypoxia")

# create vector of selected terms to plot 
GO_upregulated_female_10_vs_60_select <- c("response to corticosteroid", "response to estradiol", "p38MAPK cascade", 
                                           "intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator")

# copy enrichResult object to new plot variable
GO_upregulated_male_10_vs_60_plot <- GO_upregulated_male_10_vs_60

# subset the result slot based on selected terms
GO_upregulated_male_10_vs_60_plot@result <- GO_upregulated_male_10_vs_60@result[GO_upregulated_male_10_vs_60@result$Description %in% GO_upregulated_male_10_vs_60_select,]

# copy enrichResult object to new plot variable
GO_upregulated_female_10_vs_60_plot <- GO_upregulated_female_10_vs_60

# subset the result slot based on selected terms
GO_upregulated_female_10_vs_60_plot@result <- GO_upregulated_female_10_vs_60@result[GO_upregulated_female_10_vs_60@result$Description %in% GO_upregulated_female_10_vs_60_select,]

# combine all upregulated GO analyses into 1 df
GO_sex_up <- bind_rows(GO_male_basal_vs_10_plot@result %>%
                         mutate(Sex = "Male", Comparison = "Basal vs 10 min"),
                       GO_female_basal_vs_10_plot@result %>%
                         mutate(Sex = "Female", Comparison = "Basal vs 10 min"),
                       GO_upregulated_male_basal_vs_60_plot@result %>%
                         mutate(Sex = "Male", Comparison = "Basal vs 60 min"),
                       GO_upregulated_female_basal_vs_60_plot@result %>%
                         mutate(Sex = "Female", Comparison = "Basal vs 60 min"),
                       GO_upregulated_male_10_vs_60_plot@result %>%
                         mutate(Sex = "Male", Comparison = "10 vs 60 min"),
                       GO_upregulated_female_10_vs_60_plot@result %>%
                         mutate(Sex = "Female", Comparison = "10 vs 60 min")) %>%
  
  # calculate GeneRatio, FDR and set plotting ordet
  mutate(GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
         FDR = -log10(p.adjust),
         Sex = factor(Sex, levels = c("Male", "Female")),
         Comparison = factor(Comparison, levels = c("Basal vs 10 min", "Basal vs 60 min", "10 vs 60 min")),
         
         # simplify plot labels
         PlotLabels = case_when(Description == "response to corticosteroid" ~ "Corticosteroid response",
                                Description == "response to interleukin-1" ~ "Interleukin-1 response",
                                Description == "interleukin-1-mediated signaling pathway" ~ "Interleukin-1 response",
                                Description == "response to reactive oxygen species" ~ "ROS response",
                                Description == "response to calcium ion" ~ "Calcium ion response",
                                Description == "response to progesterone" ~ "Progesterone response",
                                Description == "regulation of synaptic plasticity" ~ "Synaptic plasticity",
                                Description == "p38MAPK cascade" ~ "p38 MAPK cascade",
                                Description == "response to hypoxia" ~ "Hypoxia response",
                                Description == "response to estradiol" ~ "Estradiol response",
                                Description == "intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator" ~ "Apoptotic signalling",
                                Description == "circadian regulation of gene expression" ~ "Circadian regulation",
                                Description == "neurotransmitter reuptake" ~ "Neurotransmitter reuptake",
                                Description == "regulation of TOR signaling" ~ "TOR signalling regulation",
                                TRUE ~ Description),
         PlotLabels = factor(PlotLabels, levels = rev(c("Corticosteroid response",
                                                        "Progesterone response",
                                                        "Estradiol response",
                                                        "p38 MAPK cascade",
                                                        "TOR signalling regulation",
                                                        "ROS response",
                                                        "Hypoxia response",
                                                        "Interleukin-1 response",
                                                        "Apoptotic signalling",
                                                        "Calcium ion response",
                                                        "Synaptic plasticity",
                                                        "Neurotransmitter reuptake",
                                                        "Circadian regulation"))))

# plot
ggplot(GO_sex_up, aes(x = Count, y = PlotLabels, fill = FDR, shape = Sex)) +
  geom_point(size = 4, colour = "black", stroke = 0.3) +
  facet_wrap(~ Comparison, nrow = 1) +
  scale_shape_manual(values = c("Male" = 21, "Female" = 24)) +
  scale_fill_gradient(low = "#FADADA", high = "#D9364A")  +
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

ggsave("4.7_GO_male_and_female_dotplot.tiff", width = 6.27, height = 3.5, dpi = 600, compress = "lzw")

# create vector of selected terms to plot 
GO_downregulated_male_basal_vs_60_select <- c("Notch signaling pathway", "stem cell differentiation",
                                              "forebrain neuron differentiation", "lymphocyte differentiation")

# create vector of selected terms to plot 
GO_downregulated_female_basal_vs_60_select <- c()

# copy enrichResult object to new plot variable
GO_downregulated_male_basal_vs_60_plot <- GO_downregulated_male_basal_vs_60

# subset the result slot based on selected terms
GO_downregulated_male_basal_vs_60_plot@result <- GO_downregulated_male_basal_vs_60@result[GO_downregulated_male_basal_vs_60@result$Description %in% GO_downregulated_male_basal_vs_60_select,]

# copy enrichResult object to new plot variable
GO_downregulated_female_basal_vs_60_plot <- GO_downregulated_female_basal_vs_60

# subset the result slot based on selected terms
GO_downregulated_female_basal_vs_60_plot@result <- GO_downregulated_female_basal_vs_60@result[GO_downregulated_female_basal_vs_60@result$Description %in% GO_downregulated_female_basal_vs_60_select,]

# create vector of selected terms to plot 
GO_downregulated_male_10_vs_60_select <- c("Notch signaling pathway", "lymphocyte differentiation",
                                           "stem cell differentiation", "forebrain neuron differentiation", "interleukin-18-mediated signaling pathway")

# create vector of selected terms to plot 
GO_downregulated_female_10_vs_60_select <- c("lymphocyte differentiation", "mRNA destabilization", "response to corticosteroid")

# copy enrichResult object to new plot variable
GO_downregulated_male_10_vs_60_plot <- GO_downregulated_male_10_vs_60

# subset the result slot based on selected terms
GO_downregulated_male_10_vs_60_plot@result <- GO_downregulated_male_10_vs_60@result[GO_downregulated_male_10_vs_60@result$Description %in% GO_downregulated_male_10_vs_60_select,]

# copy enrichResult object to new plot variable
GO_downregulated_female_10_vs_60_plot <- GO_downregulated_female_10_vs_60

# subset the result slot based on selected terms
GO_downregulated_female_10_vs_60_plot@result <- GO_downregulated_female_10_vs_60@result[GO_downregulated_female_10_vs_60@result$Description %in% GO_downregulated_female_10_vs_60_select,]

# combine all downregulated GO analyses into 1 df
GO_sex_down <- bind_rows(GO_downregulated_male_basal_vs_60_plot@result %>%
                           mutate(Sex = "Male", Comparison = "Basal vs 60 min"),
                         GO_downregulated_female_basal_vs_60_plot@result %>%
                           mutate(Sex = "Female", Comparison = "Basal vs 60 min"),
                         GO_downregulated_male_10_vs_60_plot@result %>%
                           mutate(Sex = "Male", Comparison = "10 vs 60 min"),
                         GO_downregulated_female_10_vs_60_plot@result %>%
                           mutate(Sex = "Female", Comparison = "10 vs 60 min")) %>%
  mutate(GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
         FDR = -log10(p.adjust),
         Sex = factor(Sex, levels = c("Male", "Female")),
         Comparison = factor(Comparison, levels = c("Basal vs 10 min", "Basal vs 60 min", "10 vs 60 min")),
         PlotLabels = case_when(Description == "response to corticosteroid" ~ "Corticosteroid response",
                                Description == "interleukin-18-mediated signaling pathway" ~ "Interleukin-18 signalling",
                                Description == "Notch signaling pathway" ~ "Notch signalling",
                                Description == "lymphocyte differentiation" ~ "Lymphocyte differentiation",
                                Description == "stem cell differentiation" ~ "Stem cell differentiation",
                                Description == "forebrain neuron differentiation" ~ "Forebrain neuron differentiation",
                                Description == "mRNA destabilization" ~ "mRNA destabilisation",
                                TRUE ~ Description),
         PlotLabels = factor(PlotLabels, levels = rev(c("Corticosteroid response",
                                                        "Stem cell differentiation",
                                                        "Forebrain neuron differentiation",
                                                        "Lymphocyte differentiation",
                                                        "Notch signalling",
                                                        "Interleukin-18 signalling",
                                                        "mRNA destabilisation"))))

# plot
ggplot(GO_sex_down, aes(x = Count, y = PlotLabels, fill = FDR, shape = Sex)) +
  geom_point(size = 4, colour = "black", stroke = 0.3) +
  facet_wrap(~ Comparison, nrow = 1) +
  scale_shape_manual(values = c("Male" = 21, "Female" = 24)) +
  scale_fill_gradient(low = "#C9E5F5", high = "#4A90C2") +
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

ggsave("AppendixC.4_GO_male_and_female_dotplot.tiff", width = 6, height = 2.5, dpi = 600, compression = "lzw")
