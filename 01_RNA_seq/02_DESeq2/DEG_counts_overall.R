# load required packages
library(tidyverse)
library(ggplot2)
library(ggvenn

# set working directory
setwd(here("plots"))


# summarise no. of upregulated and downregulated gene for each comparison
deg_counts <- data.frame(Comparison = c("Basal vs 10 min ARS", "Basal vs 60 min ARS", "10 vs 60 min ARS"),
  Up = c(sum(res_basal_vs_10_df$diffexpressed == "UP", na.rm = TRUE),
         sum(res_basal_vs_60_df$diffexpressed == "UP", na.rm = TRUE),
         sum(res_10_vs_60_df$diffexpressed == "UP", na.rm = TRUE)),
  Down = c(sum(res_basal_vs_10_df$diffexpressed == "DOWN", na.rm = TRUE),
           sum(res_basal_vs_60_df$diffexpressed == "DOWN", na.rm = TRUE),
           sum(res_10_vs_60_df$diffexpressed == "DOWN", na.rm = TRUE)))

# convert to long format for ggplot
deg_counts <- deg_counts %>%
  pivot_longer(cols = c(Up, Down),
               names_to = "Direction",
               values_to = "Count")

deg_counts$Comparison <- factor(deg_counts$Comparison,
                                 levels = c("Basal vs 10 min ARS",
                                            "Basal vs 60 min ARS",
                                            "10 vs 60 min ARS"))

# plot no. of DEGs as barplot
ggplot(deg_counts, aes(x = Comparison, y = Count, fill = Direction)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = Count),
            position = position_dodge(width = 0.7),
            vjust = 1.4,
            size = 2.5) +
  labs(x = NULL,
       y = "Number of DEGs") +
  scale_x_discrete(labels = c("Basal vs\n10 min", 
                              "Basal vs\n60 min", 
                              "10 vs 60\nmin")) +
  scale_y_continuous(limits = c(0, NA), 
                     expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(
    values = c("Up" = "#A995C9", "Down" = "#5FA9A2")) +
  theme_classic() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "right",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        legend.box.spacing = unit(0, "pt"))

ggsave("fig4.2a_barplot_DEG_counts.tiff", height = 1.8, width = 3, compression = "lzw")

# plot venn diagram of overlapping DEGs between basal vs 10 and basal vs 60 
deg_10 <- rownames(res_basal_vs_10_df)[res_basal_vs_10_df$diffexpressed %in% c("UP", "DOWN")]
deg_60 <- rownames(res_basal_vs_60_df)[res_basal_vs_60_df$diffexpressed %in% c("UP", "DOWN")]
deg_overlap <- list("10 min" = deg_10,
                    "60 min" = deg_60)

ggvenn(deg_overlap,
  fill_color = c("#5FA9A2", "#A995C9"),
  fill_alpha = 0.5,
  stroke_size = 0,
  set_name_size = 0,
  text_size = 2.5,
  show_percentage = F) +
  ggtitle("DEG overlap") +
  annotate("text", x = -0.7, y = 1.2,
           label = "10 min", size = 2.5) +
  annotate("text", x = 0.7, y = 1.2,
           label = "60 min", size = 2.5) +
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        plot.margin = margin(10, 30, 10, 30)) 

ggsave("fig4.2b_DEG_venn.tiff", height = 2.5, width = 3, dpi = 600, compression = "lzw")
