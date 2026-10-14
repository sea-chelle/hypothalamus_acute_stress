# load required packages
library(tidyverse)
library(ggplot2)

# set working directory
setwd(here("plots"))

# get DEG count
deg_counts_sex <- bind_rows(res_male_basal_vs_10_df %>%
                              filter(diffexpressed %in% c("UP", "DOWN")) %>%
                              summarise(DEGs = n()) %>%
                              mutate(Comparison = "Basal vs 10 min", Sex = "Male"),
                            
                            res_female_basal_vs_10_df %>%
                              filter(diffexpressed %in% c("UP", "DOWN")) %>%
                              summarise(DEGs = n()) %>%
                              mutate(Comparison = "Basal vs 10 min", Sex = "Female"),
                            
                            res_male_basal_vs_60_df %>%
                              filter(diffexpressed %in% c("UP", "DOWN")) %>%
                              summarise(DEGs = n()) %>%
                              mutate(Comparison = "Basal vs 60 min", Sex = "Male"),
                            
                            res_female_basal_vs_60_df %>%
                              filter(diffexpressed %in% c("UP", "DOWN")) %>%
                              summarise(DEGs = n()) %>%
                              mutate(Comparison = "Basal vs 60 min", Sex = "Female"),
                            
                            res_male_10_vs_60_df %>%
                              filter(diffexpressed %in% c("UP", "DOWN")) %>%
                              summarise(DEGs = n()) %>%
                              mutate(Comparison = "10 vs 60 min", Sex = "Male"),
                            
                            res_female_10_vs_60_df %>%
                              filter(diffexpressed %in% c("UP", "DOWN")) %>%
                              summarise(DEGs = n()) %>%
                              mutate(Comparison = "10 vs 60 min", Sex = "Female"))

# set order for graph
deg_counts_sex$Comparison <- factor(deg_counts_sex$Comparison,
                                    levels = c("Basal vs 10 min", "Basal vs 60 min","10 vs 60 min"))

deg_counts_sex$Sex <- factor(deg_counts_sex$Sex, levels = c("Male", "Female"))

# plot bar plot of DEGs
ggplot(deg_counts_sex, aes(x = Comparison, y = DEGs, fill = Sex)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_text(aes(label = DEGs),
            position = position_dodge(width = 0.75),
            vjust = -0.4,
            size = 2.5) +
  scale_fill_manual(values = c("Male" = "#78A8D1", "Female" = "#D78FA9")) +
  scale_x_discrete(labels = c("Basal vs\n10 min", 
                              "Basal vs\n60 min", 
                              "10 vs 60\nmin")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Number of DEGs", fill = NULL) +
  theme_classic() +
  theme(axis.title = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        legend.box.spacing = unit(0, "pt"))

ggsave("fig4.6a_sex_DEGs_barplot.tiff", height = 2.5, width = 3.5, dpi = 300)
