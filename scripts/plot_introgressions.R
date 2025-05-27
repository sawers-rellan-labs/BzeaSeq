# plot_introgressions.R
# Main script for introgression visualization analysis

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(cowplot)

# Source the introgression plotting functions
source("R/introgression_plots.R")

# Read in data from file
# Replace this with your file path
# data <- read.table("~/Desktop/all_samples_bzea_ALT.bed", sep = "\t", na.strings = na.strings = c("NA",""))
data <- read.table("/Volumes/BZea/bzeaseq/ancestry/all_samples_ancestry_segments.bed",
                   sep = "\t", na.strings = c("NA",""))
colnames(data) <- c("sample", "chrom", "start", "end", "genotype", "numeric", "freq")

# Uncomment and modify these sections as needed for your sample metadata
# bzea <- read.csv("~/Desktop/J2Teo_Final_DB.csv", na.strings = c("NA",""))
# sample_sheet <- read.csv("~/Desktop/BZea-Sample-List.csv", na.strings = c("NA","")) %>%
#   filter(!grepl("LANTEO|LANDB",Line.ID,perl=TRUE)) %>%
#   filter(grepl("B73|Purple",Line.ID,perl=TRUE))

# metadata <- sample_sheet %>% select(
#   field_row = Sample_Origin,
#   bam_preffix = Seq_Full_ID,
#   pedigree_sample_sheet = Line.ID) %>%
#     left_join(bzea %>% select(field_row =  seed_origin,
#     bam_preffix = sequencing.id,
#                              founder_group =taxa_code,
#                              donor_id = accession_id,
#                              pedigree_J2TEO = line_id
#                              )
#              ) %>% filter(!is.na(founder_group))

# Process the introgression data using the sourced functions
processed_data <- process_introgression_data(
  data   # Adjust filtering as needed
)

processed_data


# Alternative processing if you have sequenced samples to include:
# processed_data <- process_introgression_data(
#   data %>% filter(sample %in% sequenced)
# )

# Create the stacked bar plots for each chromosome
chromosome_plots <- plot_introgression_stacked(processed_data, order = "position")


# Display specific chromosome plots
quartz()
print(chromosome_plots$chr1)

# Save individual chromosome plots
pdf("~/Desktop/all_chromosome_introgressions.pdf", width = 10, height = 12)
for(chrom in paste0("chr", 1:10)) {
  print(chromosome_plots[[chrom]])
}
dev.off()

# Create theme for combined plot
chr_theme2 <- theme(
  axis.line.y = element_blank(),
  axis.text.y = element_blank(),
  axis.line.x = element_blank(),
  axis.text.x = element_blank(),
  axis.title.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank()
) 

# Create combined plot for all chromosomes
all_chromosomes_plot <- ggpubr::ggarrange(
  plotlist = lapply(1:10, function(x){
    chromosome_plots[[x]] + xlab(paste0("chr", x)) +
      coord_cartesian(expand = FALSE) + 
      coord_flip(expand = FALSE) +
      ggpubr::theme_classic2() + 
      labs(title = NULL) + 
      chr_theme2
  }), 
  common.legend = TRUE, 
  legend = "bottom", 
  nrow = 5, 
  ncol = 2
)

# Save combined chromosomes plot
ggsave("~/Desktop/introgression_all_chromosomes.png", all_chromosomes_plot, 
       width = 14, height = 10)

# Calculate and display summary statistics
summary_by_chrom <- processed_data %>%
  group_by(chrom, genotype) %>% filter(genotype !="REF") %>%
  summarize(
    segments = n(),
    samples = n_distinct(sample),
    total_mb = sum(span)/1e6,
    avg_size_mb = mean(span)/1e6,
    largest_mb = max(span)/1e6,
    .groups = "drop"
  )

cat("Introgression summary by chromosome:\n")
print(summary_by_chrom)

# Example of converting to matrix format for ComplexHeatmap (optional)
# Uncomment to test matrix conversion
# cat("\nConverting to matrix format for ComplexHeatmap...\n")
# matrix_data <- convert_to_matrix(processed_data, bin_size = 1000000)
# cat("Matrix dimensions:", dim(matrix_data$matrix), "\n")
# cat("Number of samples:", ncol(matrix_data$matrix), "\n")
# cat("Number of genomic bins:", nrow(matrix_data$matrix), "\n")