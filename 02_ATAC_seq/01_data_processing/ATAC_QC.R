# load required packages
library(here)
library(ggplot2)

# set working directory
setwd(here("plots"))

# get final peak counts
fragment_files <- list.files("qc",
                             pattern = "\\.final_fragments\\.tsv$",
                             full.names = TRUE)

# create df with final number of peak per sample
fragment_df <- do.call(rbind,
                       lapply(fragment_files, read.delim))

# get mitochondrial fraction
mito_files <- list.files("qc",
                         pattern = "\\.mito\\.tsv$",
                         full.names = TRUE)


# calculate mitochondrial fraction as a percentage and create df per sample
mito_df <- do.call(rbind,
                   lapply(mito_files, read.delim))
mito_df$Mito_percent <- mito_df$Mito_fraction * 100

# get MACS2 called peak files
peak_files <- list.files("peaks",
                         pattern = "_peaks\\.narrowPeak$",
                         recursive = TRUE,
                         full.names = TRUE)

# create df with number of called peaks per sample 
peak_df <- data.frame(Sample = basename(dirname(peak_files)),
                      Peaks = sapply(peak_files, function(x) length(readLines(x))))

# get Picard duplication metrics
dedup_files <- list.files("qc",
                          pattern = "\\.dedup\\.metrics$",
                          full.names = TRUE)

# get duplication rate
read_duplication <- function(file) {
  x <- readLines(file)
  
  # find Picard metrics header
  header <- grep("^LIBRARY\\t", x)
  
  # read header + first metrics row
  dat <- read.delim(text = paste(x[c(header, header + 1)], collapse = "\n"),
                    check.names = FALSE)
  
  data.frame(Sample = dat$LIBRARY,
             Duplication_percent = dat$PERCENT_DUPLICATION * 100)
}

# create df with duplication rate per sample 
dup_df <- do.call(rbind,
                  lapply(dedup_files, read_duplication))

# get FRiP scores
frip_files <- list.files("qc/frip",
                         pattern = "\\.frip\\.tsv$",
                         full.names = TRUE)

# create df with FRiP scores per sample
frip_df <- do.call(rbind,
                   lapply(frip_files, read.delim))

# convert FRiP scores to %
frip_df$FRiP_percent <- frip_df$FRiP * 100

# get TSS enrichment score
tss_files <- list.files("qc/tss",
                        pattern = "\\.TSSEscore\\.tsv$",
                        full.names = TRUE)

# create df with TSS enrichment scores per sample
tss_df <- do.call(rbind,
                  lapply(tss_files, read.delim))
colnames(tss_df)[colnames(tss_df) == "sample"] <- "Sample"

# merge all QC metrics
qc_df <- Reduce(function(x, y) merge(x, y, by = "Sample", all = TRUE),
  list(# FRiP
    frip_df[, c("Sample", "Fragments_in_peaks", "FRiP_percent")],
    
    # TSS enrichment
    tss_df[, c("Sample", "TSS_enrichment")],
    
    # MACS2 peaks
    peak_df,
    
    # final fragment count
    fragment_df[, c("Sample", "Final_fragments"
    )],
    
    # mitochondrial QC
    mito_df[, c("Sample", "Total_fragments", "Mito_fragments", "Mito_percent")],
    
    # duplication
    dup_df))

# simplify sample names to just animal ID
rownames(qc_df) <- sub("-.*", "", qc_df$Sample)
qc_df$Sample <- NULL

# add QC status
qc_df$QC <- ifelse(rownames(qc_df) %in% c("2095", "2124", "2127", "2131"), 
                   "Excluded", "Included")

# order by animal number
qc_df <- qc_df[order(as.numeric(rownames(qc_df))), ]

# export QC metrics df
write.csv(qc_df, "ATACseq_QC_metrics.csv", row.names = TRUE)

# plot number of peaks called per sample
peaks_plot <- ggplot(qc_df, aes(x = rownames(qc_df), y = Peaks, colour = QC)) +
  geom_point(size = 1.5) +
  labs(x = "Sample", y = "Number of peaks called") +
  scale_color_manual(values = c("Excluded" = "#A995C9", "Included" = "#5FA9A2")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.position = "none")

# plot FRiP percentage
frip_plot <- ggplot(qc_df, aes(x = rownames(qc_df), y = FRiP_percent, colour = QC)) +
  geom_point(size = 1.5) +
  labs(x = "Sample", y = "FRiP (%)") +
  scale_color_manual(values = c("Excluded" = "#A995C9", "Included" = "#5FA9A2")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.position = "none")

# plot TSS enrichment per sample
tss_plot <- ggplot(qc_df, aes(x = rownames(qc_df), y = TSS_enrichment, colour = QC)) +
  geom_point(size = 1.5) +
  labs(x = "Sample", y = "TSS enrichment score") +
  scale_color_manual(values = c("Excluded" = "#A995C9", "Included" = "#5FA9A2")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.position = "right")
