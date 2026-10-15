# load required packages
library(here)
library(tidyverse)
library(ggplot2)

# set working directory
setwd(here("plots"))

# join both basal vs 10 comparisons for male and female and label based on male and/or female DEG
fc_basal_vs_10 <- inner_join(res_male_basal_vs_10_df %>%
                               mutate(ENSEMBL = rownames(res_male_basal_vs_10_df)) %>%
                               dplyr::select(ENSEMBL, male_log2FC = log2FoldChange, male_DEG = diffexpressed),
                             
                             res_female_basal_vs_10_df %>% 
                               mutate(ENSEMBL = rownames(res_female_basal_vs_10_df)) %>%
                               dplyr::select(ENSEMBL, female_log2FC = log2FoldChange, female_DEG = diffexpressed),
                             
                             by = "ENSEMBL") %>%
  
  filter(!is.na(male_log2FC), !is.na(female_log2FC)) %>%
  mutate(DEG_group = case_when(male_DEG %in% c("UP", "DOWN") & female_DEG %in% c("UP", "DOWN") ~ "Both",
                               male_DEG %in% c("UP", "DOWN") ~ "Male only",
                               female_DEG %in% c("UP", "DOWN") ~ "Female only",
                               TRUE ~ "Neither"))

# get DESeq2 info for all DEGs at basal vs 10 in males or females
fc_basal_vs_10_DEGs <- fc_basal_vs_10 %>%
  filter(DEG_group != "Neither")

# test correlation of expression between male and females
cor_basal_vs_10 <- cor.test(fc_basal_vs_10_DEGs$female_log2FC,
                            fc_basal_vs_10_DEGs$male_log2FC,
                            method = "spearman",
                            exact = FALSE)

# join both basal vs 60 comparisons for male and female and label based on male and/or female DEG
fc_basal_vs_60 <- inner_join(res_male_basal_vs_60_df %>%
                               mutate(ENSEMBL = rownames(res_male_basal_vs_60_df)) %>%
                               dplyr::select(ENSEMBL, male_log2FC = log2FoldChange, male_DEG = diffexpressed),
                             
                             res_female_basal_vs_60_df %>% 
                               mutate(ENSEMBL = rownames(res_female_basal_vs_60_df)) %>%
                               dplyr::select(ENSEMBL, female_log2FC = log2FoldChange, female_DEG = diffexpressed),
                             
                             by = "ENSEMBL") %>%
  
  filter(!is.na(male_log2FC), !is.na(female_log2FC)) %>%
  mutate(DEG_group = case_when(male_DEG %in% c("UP", "DOWN") & female_DEG %in% c("UP", "DOWN") ~ "Both",
                               male_DEG %in% c("UP", "DOWN") ~ "Male only",
                               female_DEG %in% c("UP", "DOWN") ~ "Female only",
                               TRUE ~ "Neither"))

# get DESeq2 info for all DEGs at basal vs 60 in males or females
fc_basal_vs_60_DEGs <- fc_basal_vs_60 %>%
  filter(DEG_group != "Neither")

# test correlation of expression between male and females
cor_basal_vs_60 <- cor.test(fc_basal_vs_60_DEGs$female_log2FC,
                            fc_basal_vs_60_DEGs$male_log2FC,
                            method = "spearman",
                            exact = FALSE)

# join both 10 vs 60 comparisons for male and female and label based on male and/or female DEG
fc_10_vs_60 <- inner_join(res_male_10_vs_60_df %>%
                               mutate(ENSEMBL = rownames(res_male_10_vs_60_df)) %>%
                               dplyr::select(ENSEMBL, male_log2FC = log2FoldChange, male_DEG = diffexpressed),
                             
                             res_female_10_vs_60_df %>% 
                               mutate(ENSEMBL = rownames(res_female_10_vs_60_df)) %>%
                               dplyr::select(ENSEMBL, female_log2FC = log2FoldChange, female_DEG = diffexpressed),
                             
                             by = "ENSEMBL") %>%
  
  filter(!is.na(male_log2FC), !is.na(female_log2FC)) %>%
  mutate(DEG_group = case_when(male_DEG %in% c("UP", "DOWN") & female_DEG %in% c("UP", "DOWN") ~ "Both",
                               male_DEG %in% c("UP", "DOWN") ~ "Male only",
                               female_DEG %in% c("UP", "DOWN") ~ "Female only",
                               TRUE ~ "Neither"))

# get DESeq2 info for all DEGs at 10 vs 60 in males or females
fc_10_vs_60_DEGs <- fc_10_vs_60 %>%
  filter(DEG_group != "Neither")

# test correlation of expression between male and females
cor_10_vs_60 <- cor.test(fc_10_vs_60_DEGs$female_log2FC,
                            fc_10_vs_60_DEGs$male_log2FC,
                            method = "spearman",
                            exact = FALSE)

# scatter plot
fc_basal_vs_10_plot <- ggplot(fc_basal_vs_10_DEGs, aes(x = female_log2FC, y = male_log2FC, colour = DEG_group)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(size = 1, alpha = 0.8) +
  scale_colour_manual(values = c("Both" = "#A995C9", "Male only" = "#78A8D1", "Female only" = "#D78FA9")) +
  coord_cartesian(xlim = c(-3.5, 3.5), ylim = c(-6.5, 7.5)) +
  annotate("label", x = -3, y = 6,
           label = paste0("ρ = ", round(cor_basal_vs_10$estimate, 2)),
           hjust = 0,
           size = 2.5,
           fill = "white",
           colour = "black",
           label.padding = unit(0.25, "lines")) +
  labs(x = NULL,
       y = expression("Male log"[2]*" fold change"),
       colour = NULL) +
  ggtitle("Basal vs 10 min") + 
  theme_classic() +
  theme(plot.title = element_text(size = 9, hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "none")

# scatter plot
fc_basal_vs_60_plot <- ggplot(fc_basal_vs_60_DEGs, aes(x = female_log2FC, y = male_log2FC, colour = DEG_group)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(size = 1, alpha = 0.8) +
  scale_colour_manual(values = c("Both" = "#A995C9", "Male only" = "#78A8D1", "Female only" = "#D78FA9")) +
  coord_cartesian(xlim = c(-3.5, 3.5), ylim = c(-6.5, 7.5)) +
  annotate("label", x = -3, y = 6,
           label = paste0("ρ = ", round(cor_basal_vs_60$estimate, 2)),
           hjust = 0,
           size = 2.5,
           fill = "white",
           colour = "black",
           label.padding = unit(0.25, "lines")) +
  labs(x = NULL,
       y = expression("Male log"[2]*" fold change"),
       colour = NULL) +
  ggtitle("Basal vs 60 min") + 
  theme_classic() +
  theme(plot.title = element_text(size = 9, hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "none")

# scatter plot
fc_10_vs_60_plot <- ggplot(fc_10_vs_60_DEGs, aes(x = female_log2FC, y = male_log2FC, colour = DEG_group)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(size = 1, alpha = 0.8) +
  scale_colour_manual(values = c("Both" = "#A995C9", "Male only" = "#78A8D1", "Female only" = "#D78FA9")) +
  coord_cartesian(xlim = c(-3.5, 3.5), ylim = c(-6.5, 7.5)) +
  annotate("label", x = -3, y = 6,
           label = paste0("ρ = ", round(cor_10_vs_60$estimate, 2)),
           hjust = 0,
           size = 2.5,
           fill = "white",
           colour = "black",
           label.padding = unit(0.25, "lines")) +
  labs(x = NULL,
       y = expression("Male log"[2]*" fold change"),
       colour = NULL) +
  ggtitle("10 vs 60 min") + 
  theme_classic() +
  theme(plot.title = element_text(size = 9, hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "none")
