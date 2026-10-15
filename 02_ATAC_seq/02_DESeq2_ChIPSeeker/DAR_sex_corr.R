# load required libraries
library(here)
library(dplyr)
library(tibble)
library(ggplot2)
library(grid)

# set working directory 
setwd(here("plots"))

# combine all DARs significant in at least one sex
sex_DARs <- union(rownames(sig_male_basal_vs_60), rownames(sig_female_basal_vs_60))

fc_basal_vs_60_DARs <- tibble(peak = sex_DARs, 
                              female_log2FC = res_female_basal_vs_60_df[sex_DARs, "log2FoldChange"],
                              male_log2FC = res_male_basal_vs_60_df[sex_DARs, "log2FoldChange"]) %>%
  mutate(DAR_group = case_when(peak %in% shared_DARs ~ "Both", 
                               peak %in% male_only_DARs ~ "Male only",
                               peak %in% female_only_DARs ~ "Female only"))

# test correlation of accessibility changes between males and females
cor_basal_vs_60_DARs <- cor.test(fc_basal_vs_60_DARs$female_log2FC, fc_basal_vs_60_DARs$male_log2FC,
                                 method = "spearman",
                                 exact = FALSE)
cor_basal_vs_60_DARs
rho_basal_vs_60_DARs <- unname(cor_basal_vs_60_DARs$estimate)

# calculate directional concordance of sex-specific DARs
fc_basal_vs_60_DARs %>%
  filter(DAR_group != "Both") %>%
  mutate(same_direction = sign(female_log2FC) == sign(male_log2FC)) %>%
  group_by(DAR_group) %>%
  summarise(n_DARs = n(), n_same_direction = sum(same_direction),
            percent_same_direction = 100 * mean(same_direction),
            .groups = "drop")

# scatter plot
ggplot(fc_basal_vs_60_DARs, aes(x = female_log2FC, y = male_log2FC, colour = DAR_group)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_colour_manual(values = c("Both" = "#A995C9", "Male only" = "#78A8D1", "Female only" = "#D78FA9"),
                      breaks = c("Both", "Male only", "Female only")) +
  annotate("label", x = -Inf, y = Inf,
           label = paste0("\u03c1 = ", round(rho_basal_vs_60_DARs, 2)),
           hjust = -0.2, vjust = 1.3, size = 3, fill = "white", colour = "black", label.padding = unit(0.25, "lines")) +
  labs(x = expression("Female log"[2]*" fold change"),
       y = expression("Male log"[2]*" fold change"),
       colour = NULL) +
  ggtitle("Basal vs 60 min ARS") +
  theme_classic() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3),
        legend.text = element_text(size = 7),
        legend.position = "right",
        legend.box.spacing = unit(0, "pt"),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(size = 9, hjust = 0.5))
ggsave("fig4.13f_DAR_sex_concordance.tiff", height = 1.9, width = 4, dpi = 600, compression = "lzw")
