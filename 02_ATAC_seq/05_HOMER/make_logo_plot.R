# load required packages
library(tidyverse)
library(here)
library(magick)
library(cowplot)

# set working directory
setwd(here("plots"))

# add manually made Klf15 avg logo file for female comparison
motif_results$Female <- motif_results$Female %>%
  mutate(logo_file = if_else(motif_name == "Klf15",
                             file.path(homer_dir, homer_folders["Female"], "knownResults", "klf15.logo.svg"),
                             logo_file))

# define motifs to show for each panel
shared_logo_motifs <- c("GRE", "ERb", "THRb")
male_logo_motifs <- "Foxo1"
female_logo_motifs <- "Klf15"

# select one representative HOMER motif per simplified name - lowest q-val then highest fold enrichment if tied
select_logo_motifs <- function(data, motifs) {
  
  data %>%
    filter(motif_name %in% motifs) %>%
    group_by(motif_name) %>%
    arrange(qvalue, desc(fold_enrichment), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    ungroup()
}

# get data for each panel
shared_logos <- select_logo_motifs(motif_results$Overall, shared_logo_motifs)
male_logos <- select_logo_motifs(motif_results$Male, male_logo_motifs)
female_logos <- select_logo_motifs(motif_results$Female, female_logo_motifs)

# display names
display_names <- c("GRE" = "GRE",
                   "ERb" = "ERβ",
                   "THRb" = "THRβ",
                   "Foxo1" = "FOXO1",
                   "Klf15" = "KLF15")

shared_logos <- shared_logos %>%
  mutate(display_name = dplyr::recode(motif_name, !!!display_names),
         motif_name = factor(motif_name, levels = shared_logo_motifs)) %>%
  arrange(motif_name)

male_logos <- male_logos %>%
  mutate(display_name = dplyr::recode(motif_name, !!!display_names),
         motif_name = factor(motif_name, levels = male_logo_motifs)) %>%
  arrange(motif_name)

female_logos <- female_logos %>%
  mutate(display_name = dplyr::recode(motif_name, !!!display_names),
         motif_name = factor(motif_name, levels = female_logo_motifs)) %>%
  arrange(motif_name)

# read logo images
shared_logos$logo_img <- map(shared_logos$logo_file, magick::image_read_svg)
male_logos$logo_img <- map(male_logos$logo_file, magick::image_read_svg)
female_logos$logo_img <- map(female_logos$logo_file, magick::image_read_svg)

# function to make motif logo panel
make_motif_panel <- function(df, title, title_x) {
  
  rows <- purrr::map2(df$logo_img, df$display_name, ~ {
    label_plot <- ggdraw() +
      draw_label(.y, x = 0.85, y = 0.5, hjust = 0, vjust = 0.5, fontface = "plain", size = 8)
    
    logo_plot <- ggdraw() + 
      draw_image(.x, x = 0, y = 0.10, width = 0.45, height = 0.80)
    
    plot_grid(label_plot, logo_plot, ncol = 2, rel_widths = c(0.32, 0.68))
    })
  
  body <- plot_grid(plotlist = rows, ncol = 1)
  
  title_plot <- ggdraw() +
    draw_label(title, x = title_x, y = 0.15, hjust = 0.5, vjust = 0.5, fontface = "bold", size = 9)
  
  plot_grid(title_plot, body, ncol = 1, rel_heights = c(0.55, nrow(df)))
}

# combine all plots together
plot_grid(make_motif_panel(shared_logos, "Shared", title_x = 0.475),
  ggdraw(),
  make_motif_panel(male_logos, "Male", title_x = 0.47),
  ggdraw(),
  make_motif_panel(female_logos, "Female", title_x = 0.475),
  ncol = 1,
  rel_heights = c(1.5,   
                  0.1,  
                  1,
                  0.1,
                  1))
ggsave("HOMER_motif_logos.tiff", height = 4, width = 7, dpi = 300, compression = "lzw")
