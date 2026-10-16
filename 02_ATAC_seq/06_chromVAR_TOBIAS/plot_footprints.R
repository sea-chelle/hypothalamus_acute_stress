# load required libraries
library(tidyverse)
library(here)
library(zoo)

# set working directory
setwd(here("plots"))

# get IDs of all significant motifs and write as text file
tf_final %>%
  distinct(motif_id) %>%
  pull(motif_id) %>%
  writeLines(file.path(proj_dir, "TOBIAS", "aggregate", "tobias_motif_ids.txt"))

# get aggregate directory
aggregate_dir <- file.path(proj_dir, "TOBIAS", "aggregate")

# find all aggregate files
aggregate_files <- list.files(aggregate_dir, pattern = "_aggregate\\.txt$", full.names = TRUE)

# function to import TOBIAS aggregate footprint profiles
read_aggregate <- function(file, tf) {
  
  data <- read.delim(file, skip = 2, header = FALSE, 
                     col.names = c("condition", "regions", "aggregate"), 
                     stringsAsFactors = FALSE)
  
  data %>%
    separate_wider_delim(condition, delim = "_", names = c("sex", "treatment")) %>%
    mutate(TF = tf, aggregate = str_split(aggregate, ",")) %>%
    unnest_longer(aggregate, indices_to = "position") %>%
    mutate(position = as.integer(position), signal = as.numeric(aggregate)) %>%
    dplyr::select(TF, sex, treatment, position, signal)
}

# import and combine
footprint_data <- purrr::map_dfr(aggregate_files, function(file) {
  motif_id <- stringr::str_extract(basename(file), "MA[0-9]+\\.[0-9]+")
  tf <- tf_final$tf[match(motif_id, tf_final$motif_id)][1]
  read_aggregate(file, tf)
})

# add groups
footprint_data <- footprint_data %>%
  mutate(condition_group = case_when(sex == "Male" & treatment == "Basal" ~ "Male Basal",
                                      sex == "Male" & treatment == "60" ~ "Male 60 min",
                                      sex == "Female" & treatment == "Basal" ~ "Female Basal",
                                      sex == "Female" & treatment == "60" ~ "Female 60 min"))

# centre footprint positions on motif
footprint_data <- footprint_data %>%
  group_by(TF, sex, treatment) %>%
  arrange(position, .by_group = TRUE) %>%
  mutate(position = row_number() - (n() + 1) / 2) %>%
  ungroup()

# check imported motifs
dplyr::count(footprint_data, TF) %>% 
  print(n = Inf)

# select and order representative TFs
footprint_selected <- footprint_data %>%
  filter(TF %in% c("Ar", "Atf3", "Irf1", "Lef1")) %>%
  group_by(TF, sex, treatment) %>%
  arrange(position, .by_group = TRUE) %>%
  mutate(TF = factor(TF, levels = c("Ar", "Atf3", "Irf1", "Lef1")),
         signal_smooth = zoo::rollmean(signal, k = 10, fill = NA, align = "center"),
         condition_group = dplyr::recode(condition_group,
                                         "Male 60" = "Male 60 min",
                                         "Female 60" = "Female 60 min")) %>%
  ungroup()

# plot vertically stacked footprints
ggplot(footprint_selected, aes(x = position, y = signal_smooth, colour = condition_group)) +
  geom_line(linewidth = 0.5) +
  facet_grid(TF ~ ., scales = "free_y") +
  scale_colour_manual(values = c("Male Basal" = "#78A8D1",
                                 "Male 60 min" = "#47749B",
                                 "Female Basal" = "#D78FA9",
                                 "Female 60 min" = "#B65679"),
                      breaks = c("Male Basal", "Male 60 min", "Female Basal", "Female 60 min")) +
  scale_x_continuous(breaks = seq(-60, 60, 30)) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(x = "Distance from motif centre (bp)", y = "Mean bias-corrected ATAC-seq signal", colour = "Condition") +
  theme_classic() +
  theme(strip.background = element_blank(),
        
        strip.text.y = element_text(face = "bold", angle = 0, size = 9),
        strip.placement = "outside",
        panel.spacing.y = unit(0.8, "lines"),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.margin = margin(b = 0),
        legend.box.spacing = unit(0, "mm"),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3))
ggsave("TOBIAS_footprinting.tiff", height = 6, width = 3.7, dpi = 300, compression = "lzw")
