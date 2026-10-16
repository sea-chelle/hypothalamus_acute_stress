# load required libraries
library(tidyverse)
library(here)

# set working directory
setwd(here("tobias"))

# import sex-specific BINDetect results
tobias_male <- read.delim("BINDetect/male/male_basal_vs_60_results.txt", check.names = FALSE)
tobias_female <- read.delim("BINDetect/female/female_basal_vs_60_results.txt", check.names = FALSE)

# retain relevant columns and standardise names
tobias_male <- tobias_male %>%
  dplyr::select(motif_id, name, cluster, total_tfbs,
                male_basal_mean_score, male_60_mean_score,
                male_basal_bound, male_60_bound,
                male_change = male_basal_male_60_change,
                male_pvalue = male_basal_male_60_pvalue,
                male_highlighted = male_basal_male_60_highlighted)

tobias_female <- tobias_female %>%
  dplyr::select(motif_id, name, cluster, total_tfbs,
                female_basal_mean_score, female_60_mean_score,
                female_basal_bound, female_60_bound,
                female_change = female_basal_female_60_change,
                female_pvalue = female_basal_female_60_pvalue,
                female_highlighted = female_basal_female_60_highlighted)


# combine male and female results and derive summary variables
tobias_results <- left_join(tobias_male, tobias_female %>%
                              dplyr::select(-name, -cluster, -total_tfbs),
                            by = "motif_id") %>%
  mutate(male_highlighted = male_highlighted == "True",
         female_highlighted = female_highlighted == "True",
         highlight_group = case_when(male_highlighted & female_highlighted ~ "Both",
                                     male_highlighted & !female_highlighted ~ "Male only",
                                     female_highlighted & !male_highlighted ~ "Female only",
                                     TRUE ~ "Neither"),
         sex_difference = male_change - female_change,
         male_bound_change = male_60_bound - male_basal_bound,
         female_bound_change = female_60_bound - female_basal_bound,
         male_bound_percent_change =100 * male_bound_change / male_basal_bound,
         female_bound_percent_change = 100 * female_bound_change / female_basal_bound,
         male_bound_direction = case_when(male_bound_change > 0 ~ "Gain",
                                          male_bound_change < 0 ~ "Loss",
                                          TRUE ~ "No change"),
         female_bound_direction = case_when(female_bound_change > 0 ~ "Gain",
                                            female_bound_change < 0 ~ "Loss",
                                            TRUE ~ "No change"))

# summarise BINDetect-highlighted motifs
table(tobias_results$highlight_group)

tobias_results %>%
  filter(highlight_group != "Neither") %>%
  dplyr::select(name, cluster, female_change, male_change, sex_difference, highlight_group) %>%
  arrange(factor(highlight_group, levels = c("Both", "Female only", "Male only")),
          female_change)
