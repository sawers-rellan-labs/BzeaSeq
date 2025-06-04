# plot_introgressions_by_donor.R
# Updated script to use BzeaSeq package functions

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(cowplot)
library(BzeaSeq)  # Use package functions instead of sourcing

# Read in data from file
data <- read.table("/Volumes/BZea/bzeaseq/ancestry/all_samples_ancestry_segments.bed",
                   sep = "\t", na.strings = c("NA",""), skip = 1)
colnames(data) <- c("sample", "chrom", "start", "end", "genotype", "numeric", "freq")

# Read metadata (update paths as needed)
bzea <- read.csv("~/Desktop/J2Teo_Final_DB.csv", na.strings = c("NA",""))
sample_sheet <- read.csv("~/Desktop/BZea-Sample-List.csv", na.strings = c("NA",""))
sample_sheet$row <- gsub("PV23-","",sample_sheet$Sample_Origin) %>% as.integer()
sample_sheet$is_check <- grepl("B73|Purple",sample_sheet$Line.ID, perl=TRUE)

# Process metadata
in_bzea <- sample_sheet$row > 7369
in_lanteo <- sample_sheet$row < 7288
in_landb <- sample_sheet$row > 7288 & sample_sheet$row <= 7369

sample_sheet$project[in_bzea] <- "bzea"
sample_sheet$project[in_lanteo] <- "lanteo"
sample_sheet$project[in_landb] <- "landb"

metadata <- sample_sheet %>% 
  select(
    field_row = Sample_Origin,
    project,
    sample = Seq_Full_ID,
    is_check,
    pedigree_sample_sheet = Line.ID
  ) %>%
  left_join(bzea %>% select(
    field_row = seed_origin,
    sample = sequencing.id,
    sample_group = taxa_code,
    founder_group = taxa_code,
    donor_id = accession_id,
    pedigree = line_id
  )) 

metadata$sample_group[grepl("Purple",metadata$pedigree_sample_sheet)] <- "Purple"
metadata$sample_group[grepl("B73",metadata$pedigree_sample_sheet)] <- "B73"
metadata$founder_group[grepl("B73",metadata$pedigree_sample_sheet)] <- "B73"
metadata$pedigree_sample_sheet <- NULL

NILs <- metadata$sample[in_bzea & !metadata$is_check]

# Add metadata to data and prepare for processing
data_full <- data %>%
  inner_join(metadata) %>%
  mutate(sample = pedigree) %>%
  select(sample, chrom, start, end, genotype, numeric, freq, donor_id)

cat("Found", length(unique(data_full$donor_id)), "unique donor accessions\n")

# Process introgression data using package function
processed_data <- process_introgression_data(data_full)

# Function to split samples into chunks
chunk_samples <- function(samples, chunk_size = 20) {
  split(samples, ceiling(seq_along(samples) / chunk_size))
}

# Split data by donor and process each
donor_plots <- data_full %>%
  split(.$donor_id) %>%
  lapply(function(donor_data) {
    donor <- unique(donor_data$donor_id)
    unique_samples <- unique(donor_data$sample)
    n_samples <- length(unique_samples)
    
    cat("Processing donor:", donor, "with", n_samples, "NILs\n")
    
    if (n_samples <= 20) {
      # Single plot for this donor using package functions
      plot_data <- donor_data %>% select(-donor_id)
      processed_plot_data <- process_introgression_data(plot_data)
      chromosome_plots <- plot_introgression_stacked(processed_plot_data, order = "sample")
      
      # Theme for combining chromosomes
      chr_theme2 <- theme(
        axis.line.y = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()
      )
      
      # Combine all chromosomes using package plot functions
      all_chromosomes_plot <- ggpubr::ggarrange(
        plotlist = lapply(paste0("chr", 1:10), function(chrom_name) {
          if (!is.null(chromosome_plots[[chrom_name]])) {
            chromosome_plots[[chrom_name]] + 
              xlab(chrom_name) +
              coord_cartesian(expand = FALSE) + 
              coord_flip(expand = FALSE) +
              ggpubr::theme_classic2() + 
              labs(title = NULL) + 
              chr_theme2
          } else {
            ggplot() + theme_void() + labs(title = paste(chrom_name, "(no data)"))
          }
        }), 
        common.legend = TRUE, 
        legend = "bottom", 
        nrow = 5, 
        ncol = 2
      )
      
      # Return single plot with donor title
      list(ggpubr::annotate_figure(
        all_chromosomes_plot,
        top = ggpubr::text_grob(
          paste("Teosinte Donor Accession:", donor), 
          face = "bold", 
          size = 16
        )
      ))
    } else {
      # Multiple plots for this donor
      sample_chunks <- chunk_samples(sort(unique_samples), 20)
      
      lapply(seq_along(sample_chunks), function(chunk_idx) {
        chunk_samples_list <- sample_chunks[[chunk_idx]]
        
        plot_data <- donor_data %>%
          filter(sample %in% chunk_samples_list) %>%
          select(-donor_id)
        
        # Use package functions for processing and plotting
        processed_plot_data <- process_introgression_data(plot_data)
        chromosome_plots <- plot_introgression_stacked(processed_plot_data, order = "sample")
        
        # Theme for combining chromosomes
        chr_theme2 <- theme(
          axis.line.y = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()
        )
        
        # Combine all chromosomes
        all_chromosomes_plot <- ggpubr::ggarrange(
          plotlist = lapply(paste0("chr", 1:10), function(chrom_name) {
            if (!is.null(chromosome_plots[[chrom_name]])) {
              chromosome_plots[[chrom_name]] + 
                xlab(chrom_name) +
                coord_cartesian(expand = FALSE) + 
                coord_flip(expand = FALSE) +
                ggpubr::theme_classic2() + 
                labs(title = NULL) + 
                chr_theme2
            } else {
              ggplot() + theme_void() + labs(title = paste(chrom_name, "(no data)"))
            }
          }), 
          common.legend = TRUE, 
          legend = "bottom", 
          nrow = 5, 
          ncol = 2
        )
        
        # Add donor title with chunk info
        ggpubr::annotate_figure(
          all_chromosomes_plot,
          top = ggpubr::text_grob(
            paste("Teosinte Donor Accession:", donor, "- Part", chunk_idx, "of", length(sample_chunks)), 
            face = "bold", 
            size = 16
          )
        )
      })
    }
  }) %>%
  unlist(recursive = FALSE)  # Flatten the nested list structure

# Create multipage PDF
pdf("~/Desktop/introgressions_by_donor_accession.pdf", width = 14, height = 10)
lapply(donor_plots, print)
dev.off()

cat("PDF saved to: ~/Desktop/introgressions_by_donor_accession.pdf\n")