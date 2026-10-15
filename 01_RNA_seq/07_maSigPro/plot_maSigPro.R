# load required packages
library(here)
library(tidyverse)
library(ggplot2)

# set working directory
setwd(here("plots"))

# format data for plotting
masig_plot <- as.data.frame(vsd_mat[names(see$cut), ]) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expression") %>%
  mutate(cluster = factor(see$cut[gene], levels = 1:5)) %>%
  left_join(masig_df %>%
            rownames_to_column("sample") %>%
            mutate(sex = ifelse(sex == "m", "Male", "Female"),
                   time = as.numeric(condition)),
            by = "sample") %>%
  group_by(cluster, sex, time) %>%
  summarise(expression = median(expression), .groups = "drop")

# plot
ggplot(masig_plot, aes(x = time, y = expression, group = sex, colour = sex)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2.5) +
  facet_wrap(~cluster, ncol = 2, labeller = labeller(cluster = c(`1` = "Cluster 1 (n = 19)",
                                                                 `2` = "Cluster 2 (n = 7)",
                                                                 `3` = "Cluster 3 (n = 8)",
                                                                 `4` = "Cluster 4 (n = 7)",
                                                                 `5` = "Cluster 5 (n = 13)")),
             scales = "free_y") +
  scale_colour_manual(values = c("Male" = "#78A8D1", "Female" = "#D78FA9")) +
  scale_x_continuous(breaks = c(0, 10, 60)) +
  labs(x = "ARS duration (min)", y = "Median VST expression", colour = NULL) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 9),
        legend.position = "top",
        legend.text = element_text(size = 7),
        legend.key.size = unit(3, "mm"),
        legend.key.height = unit(3, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.box.spacing = unit(0, "pt"),
        panel.spacing = unit(1, "lines"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3))

ggsave("4.9_masigpro.tiff", height = 5, width = 6.27, dpi = 600, compression = "lzw")
