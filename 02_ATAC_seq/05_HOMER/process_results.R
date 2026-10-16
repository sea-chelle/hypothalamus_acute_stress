# load required packages
library(tidyverse)
library(here)

# set working directory
setwd(here("HOMER"))

# import HOMER known motif enrichment results
motif_results <- list(Overall = read_tsv(file.path("overall_UP", "knownResults.txt"), 
                                         show_col_types = FALSE),
                      
                      Male = read_tsv(file.path("male_UP", "knownResults.txt"), 
                                      show_col_types = FALSE),
                      
                      Female = read_tsv(file.path("female_UP", "knownResults.txt"), 
                                        show_col_types = FALSE),
                      
                      Overall_DOWN = read_tsv(file.path("overall_DOWN", "knownResults.txt"), 
                                              show_col_types = FALSE))

# get homer .svg logo files
homer_dir <- here("HOMER")
homer_folders <- c(Overall = "overall_UP", Male = "male_UP", Female = "female_UP", Overall_DOWN = "overall_DOWN")

# function for grabbing file path for .svg files
get_homer_logo_lookup <- function(folder) {
  
  known_dir <- file.path(homer_dir, folder, "knownResults")
  motif_files <- list.files(known_dir, pattern = "^known[0-9]+\\.motif$", full.names = TRUE)
  
  tibble(motif_file = motif_files) %>%
    mutate(homer_file_rank = as.numeric(str_extract(basename(motif_file), "[0-9]+")),
      header = map_chr(motif_file, ~ readLines(.x, n = 1, warn = FALSE)),
      consensus = map_chr(strsplit(header, "\t"), ~ str_remove(.x[1], "^>")),
      motif = map_chr(strsplit(header, "\t"), ~ .x[2]),
      logo_file = sub("\\.motif$", ".logo.svg", motif_file)) %>%
    arrange(homer_file_rank) %>%
    group_by(motif, consensus) %>% 
    mutate(occurrence = row_number()) %>%
    ungroup() %>%
    dplyr::select(motif, consensus, occurrence, homer_file_rank, motif_file, logo_file)
}

# combine HOMER results into df
motif_results <- lapply(homer_folders, function(folder) {
  
  results <- read_tsv(file.path(homer_dir, folder, "knownResults.txt"),
                      show_col_types = FALSE) %>%
    
    # simplify HOMER column names
    dplyr::rename(motif = `Motif Name`,
                  consensus = Consensus,
                  pvalue = `P-value`,
                  log_pvalue = `Log P-value`,
                  qvalue = `q-value (Benjamini)`,
                  target_n = starts_with("# of Target Sequences"),
                  target_percent = `% of Target Sequences with Motif`,
                  background_n = starts_with("# of Background Sequences"),
                  background_percent = `% of Background Sequences with Motif`) %>%
    
    # add and calculate relevant info
    mutate(target_percent = parse_number(target_percent),
           background_percent = parse_number(background_percent),
           fold_enrichment = target_percent / background_percent,
           log2_fold_enrichment = log2(fold_enrichment),
           motif_short = str_extract(motif, "^[^/]+"),
           motif_name = str_remove(motif_short, "\\(.*$")) %>%
    
    group_by(motif, consensus) %>%
    mutate(occurrence = row_number()) %>%
    ungroup()
  
  # add .svg file using created function
  logo_lookup <- get_homer_logo_lookup(folder)
  
  results %>%
    left_join(logo_lookup, 
              by = c("motif", "consensus", "occurrence"),
              relationship = "one-to-one")
})

#combine all UP DAR comparisons into one df
motif_up <- bind_rows(motif_results[c("Overall", "Male", "Female")],
                      .id = "comparison")

# see how many significantly enriched motifs there are for each comparison
motif_up %>%
  group_by(comparison) %>%
  summarise(n_significant = sum(qvalue < 0.05),
            n_total = n(),
            .groups = "drop")
