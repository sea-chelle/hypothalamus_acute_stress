# load required packages
library(here)
library(tidyverse)
library(ggplot2)
library(grid)

# set working directory 
setwd(here("plots"))

# summarise no. of upregulated and downregulated DARs for each comparison
dar_counts <- data.frame(Comparison = c("Overall", "Male", "Female"),
                         Up = c(sum(res_basal_vs_60_df$diffaccessible == "UP", na.rm = TRUE),
                                sum(res_male_basal_vs_60_df$diffaccessible == "UP", na.rm = TRUE),
                                sum(res_female_basal_vs_60_df$diffaccessible == "UP", na.rm = TRUE)),
                         Down = c(sum(res_basal_vs_60_df$diffaccessible == "DOWN", na.rm = TRUE),
                                sum(res_male_basal_vs_60_df$diffaccessible == "DOWN", na.rm = TRUE),
                                sum(res_female_basal_vs_60_df$diffaccessible == "DOWN", na.rm = TRUE)))

# convert to long format for ggplot
dar_counts <- dar_counts %>%
  pivot_longer(cols = c(Up, Down),
               names_to = "Direction",
               values_to = "Count")

dar_counts$Comparison <- factor(dar_counts$Comparison,
                                levels = c("Overall", "Male", "Female"))

dar_counts$Direction <- factor(dar_counts$Direction, levels = c("Down", "Up"))

# calculate total number of DARs for each comparison
dar_totals <- dar_counts %>%
  group_by(Comparison) %>%
  summarise(Total = sum(Count), .groups = "drop")

# plot number of DARs
ggplot(dar_counts, aes(x = Comparison, y = Count, fill = Direction)) +
  geom_col(width = 0.7) +
  geom_text(data = filter(dar_counts, Direction == "Up"), aes(label = Count), 
            position = position_stack(vjust = 0.5), size = 2) +
  geom_text(data = filter(dar_counts, Direction == "Down") %>%
              left_join(dar_totals, by = "Comparison"),
            aes(x = Comparison, y = Total, label = Count),
            inherit.aes = FALSE,
            nudge_y = 70, size = 2) +
  labs(x = NULL, y = "Number of DARs", fill = "Direction") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("Up" = "#A995C9", "Down" = "#5FA9A2")) +
  theme_classic() +
  theme(axis.text = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "right",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.box.spacing = unit(0, "pt"),
        legend.key.size = unit(0.3, "cm")) +
  coord_flip()
ggsave("fig4.13a_barplot_DAR_counts.tiff", height = 1.8, width = 3, dpi = 600, compression = "lzw")
