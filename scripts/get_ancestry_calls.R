#!/usr/bin/env Rscript
# BinAncestry.R
# Script to process GATK AlleleCountsOutput and generate binned ancestry calls
# Usage: Rscript BinAncestry.R input.allelicCounts.tsv output_prefix [bin_size]
# Example: Rscript BinAncestry.R sample1.allelicCounts.tsv sample1 100000

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(Ckmeans.1d.dp)
})

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript BinAncestry.R input.allelicCounts.tsv output_prefix [bin_size]\n")
  cat("Example: Rscript BinAncestry.R sample1.allelicCounts.tsv sample1 1000000\n")
  quit(status = 1)
}

input_file <- args[1]
output_prefix <- args[2]
bin_size <- ifelse(length(args) >= 3, as.numeric(args[3]), 1000000)

# Output messages to both console and log file
log_file <- paste0(output_prefix, ".log")
log_con <- file(log_file, "w")
sink(log_con, type = "output", append = FALSE)
sink(log_con, type = "message", append = FALSE)

cat(paste("Processing file:", input_file, "\n"))
cat(paste("Output prefix:", output_prefix, "\n"))
cat(paste("Bin size:", bin_size, "bp\n"))

# Start timing
start_time <- Sys.time()
cat(paste("Started at:", start_time, "\n"))

# Read the allelic counts data
cat("Reading allelic counts data...\n")
read_count <- try(read.table(input_file, comment.char = "@", header = TRUE))

if (inherits(read_count, "try-error")) {
  cat("Error reading input file. Please check the file format.\n")
  quit(status = 1)
}

# Print basic stats
cat(paste("Total variants:", nrow(read_count), "\n"))
cat(paste("Chromosomes found:", paste(unique(read_count$CONTIG), collapse=", "), "\n"))

# Create bins and calculate statistics
cat("Creating bins and calculating statistics...\n")
read_freq <- read_count %>%
  group_by(CONTIG) %>%
  # Create position bins by dividing positions into bins of size bin_size
  mutate(BIN_POS = ceiling(POSITION/bin_size) %>% as.integer()) %>%
  # Group by these bins
  group_by(CONTIG, BIN_POS) %>%
  # Calculate statistics for each bin
  summarise(
    VARIANT_COUNT = n(),                                   # Number of variants in bin
    READ_COUNT = (sum(ALT_COUNT) + sum(REF_COUNT)),        # Total read count
    ALT_COUNT = sum(ALT_COUNT),                            # Sum of alternative allele counts
    ALT_FREQ = sum(ALT_COUNT)/READ_COUNT,                  # Alternative allele frequency
    BIN_START = min(POSITION),                             # Bin start position
    BIN_END = max(POSITION),                               # Bin end position
    .groups = 'drop'
  )

# Display bin summary
total_bins <- nrow(read_freq)
cat(paste("Total number of bins:", total_bins, "\n"))
cat(paste("Average number of variants per bin:", round(mean(read_freq$VARIANT_COUNT), 1), "\n"))
cat(paste("Median number of variants per bin:", median(read_freq$VARIANT_COUNT), "\n"))
cat(paste("Range of variants per bin:", min(read_freq$VARIANT_COUNT), "to", max(read_freq$VARIANT_COUNT), "\n"))

# Normalize by variant count
cat("Normalizing allele frequencies...\n")
mean_VARIANT_COUNT <- mean(read_freq$VARIANT_COUNT)
read_freq <- read_freq %>%
  mutate(
    # Normalize by variant count
    ALT_FREQ_NORM = ALT_FREQ * VARIANT_COUNT/mean_VARIANT_COUNT,
  )

# Identify zero-ALT bins
zero_alt_bins <- sum(read_freq$ALT_FREQ_NORM == 0)
cat(paste("Bins with ALT_FREQ_NORM = 0:", 
          zero_alt_bins, "out of", total_bins, 
          "(", round(100 * zero_alt_bins / total_bins, 1), "%)\n"))

# Create a separate dataset for non-zero bins (for clustering)
non_zero_bins <- read_freq %>%
  filter(ALT_FREQ_NORM > 0)

cat(paste("Bins with ALT_FREQ_NORM > 0 (to be clustered):", 
          nrow(non_zero_bins), "out of", total_bins,
          "(", round(100 * nrow(non_zero_bins) / total_bins, 1), "%)\n"))

# Apply K-means clustering with k=3 to the non-zero bins
cat("Applying K-means clustering to non-zero bins...\n")
if (nrow(non_zero_bins) > 0) {
  non_zero_bins$K <- as.factor(Ckmeans.1d.dp(non_zero_bins$ALT_FREQ_NORM, 3)$cluster)
  
  # Function to relabel clusters as REF, HET, ALT based on ALT_FREQ_NORM
  relabel_clusters <- function(clusters, data) {
    # Calculate mean ALT_FREQ_NORM for each cluster
    cluster_means <- tapply(data$ALT_FREQ_NORM, clusters, mean)
    
    # Order clusters by mean ALT_FREQ_NORM
    ordered_clusters <- order(cluster_means)
    
    # Create mapping from original cluster to REF, HET, ALT
    cluster_map <- rep(NA, length(unique(as.numeric(clusters))))
    cluster_map[ordered_clusters[1]] <- "REF"
    cluster_map[ordered_clusters[2]] <- "HET"
    cluster_map[ordered_clusters[3]] <- "ALT"
    
    # Apply mapping to original clusters
    return(factor(cluster_map[as.numeric(clusters)], levels = c("REF", "HET", "ALT")))
  }
  
  # Apply relabeling to the clustering
  non_zero_bins$K <- relabel_clusters(non_zero_bins$K, non_zero_bins)
  
  # Check cluster proportions
  cluster_props <- prop.table(table(non_zero_bins$K)) * 100
  cat("K-means cluster proportions (non-zero bins):\n")
  print(round(cluster_props, 1))
} else {
  cat("No non-zero bins found for clustering. All bins will be assigned as REF.\n")
}

# Merge the clustering results back to the full dataset
cat("Merging results back to full dataset...\n")
# Create a new column in read_freq for the clustering
read_freq$GENOTYPE <- factor(NA, levels = c("REF", "HET", "ALT"))

# If there are non-zero bins to cluster
if (nrow(non_zero_bins) > 0) {
  # First, assign the non-zero bin clustering results to the main dataset
  # Create an index for non-zero bins in the original data
  non_zero_idx <- which(read_freq$ALT_FREQ_NORM > 0)
  
  # Assign cluster classifications to non-zero bins
  # Match by CONTIG and BIN_POS to ensure correct assignment
  for (i in 1:nrow(non_zero_bins)) {
    chr <- non_zero_bins$CONTIG[i]
    bin <- non_zero_bins$BIN_POS[i]
    match_idx <- which(read_freq$CONTIG == chr & read_freq$BIN_POS == bin)
    if (length(match_idx) > 0) {
      read_freq$GENOTYPE[match_idx] <- non_zero_bins$K[i]
    }
  }
}

# Then, force bins with ALT_FREQ_NORM = 0 to be "REF"
zero_idx <- which(read_freq$ALT_FREQ_NORM == 0)
read_freq$GENOTYPE[zero_idx] <- "REF"

# Check final genotype proportions
final_props <- prop.table(table(read_freq$GENOTYPE)) * 100
cat("Final genotype proportions:\n")
print(round(final_props, 1))

# Prepare output for each chromosome
cat("Writing results to disk...\n")

# Save detailed binned results
output_file <- paste0(output_prefix, "_bin_genotypes.tsv")
write.table(read_freq, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste("Detailed bin results saved to:", output_file, "\n"))

# Create a simplified BED format for visualization (with genotype as score)
# This is useful for loading into genome browsers
bed_file <- paste0(output_prefix, "_genotypes.bed")
read_freq_bed <- read_freq %>%
  mutate(
    SAMPLE=output_prefix,
    SCORE = case_when(
      GENOTYPE == "REF" ~ 0,
      GENOTYPE == "HET" ~ 1,
      GENOTYPE == "ALT" ~ 2,
      TRUE ~ NA_real_
    )
  ) %>%
  select(SAMPLE, CONTIG, BIN_START, BIN_END, GENOTYPE, SCORE, ALT_FREQ_NORM)

write.table(read_freq_bed, bed_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
cat(paste("BED format results saved to:", bed_file, "\n"))

# Create a chromosome-level summary
chrom_summary <- read_freq %>%
  group_by(CONTIG) %>%
  summarize(
    total_bins = n(),
    REF_bins = sum(GENOTYPE == "REF", na.rm = TRUE),
    HET_bins = sum(GENOTYPE == "HET", na.rm = TRUE),
    ALT_bins = sum(GENOTYPE == "ALT", na.rm = TRUE),
    REF_pct = 100 * sum(GENOTYPE == "REF", na.rm = TRUE) / n(),
    HET_pct = 100 * sum(GENOTYPE == "HET", na.rm = TRUE) / n(),
    ALT_pct = 100 * sum(GENOTYPE == "ALT", na.rm = TRUE) / n(),
    non_ref_pct = 100 * sum(GENOTYPE %in% c("HET", "ALT"), na.rm = TRUE) / n()
  )

summary_file <- paste0(output_prefix, "_chromosome_summary.tsv")
write.table(chrom_summary, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste("Chromosome-level summary saved to:", summary_file, "\n"))

# Calculate overall statistics
overall_summary <- data.frame(
  sample = output_prefix,
  total_bins = nrow(read_freq),
  REF_bins = sum(read_freq$GENOTYPE == "REF", na.rm = TRUE),
  HET_bins = sum(read_freq$GENOTYPE == "HET", na.rm = TRUE),
  ALT_bins = sum(read_freq$GENOTYPE == "ALT", na.rm = TRUE),
  REF_pct = 100 * sum(read_freq$GENOTYPE == "REF", na.rm = TRUE) / nrow(read_freq),
  HET_pct = 100 * sum(read_freq$GENOTYPE == "HET", na.rm = TRUE) / nrow(read_freq),
  ALT_pct = 100 * sum(read_freq$GENOTYPE == "ALT", na.rm = TRUE) / nrow(read_freq),
  non_ref_pct = 100 * sum(read_freq$GENOTYPE %in% c("HET", "ALT"), na.rm = TRUE) / nrow(read_freq)
)

overall_file <- paste0(output_prefix, "_summary.tsv")
write.table(overall_summary, overall_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste("Overall summary saved to:", overall_file, "\n"))

# End timing
end_time <- Sys.time()
time_taken <- end_time - start_time
cat(paste("Completed at:", end_time, "\n"))
cat(paste("Time taken:", round(time_taken, 2), attr(time_taken, "units"), "\n"))

# Close log
sink(type = "output")
sink(type = "message")
close(log_con)