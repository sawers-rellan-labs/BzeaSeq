#!/usr/bin/env Rscript
# get_genotypes_on_all_samples.R
#
# Description: Processes bin-level allelic count data from multiple samples,
# performs clustering across all samples, and applies HMM smoothing.
# Also analyzes segment lengths and transition frequencies.
#
# Usage: Rscript get_genotypes_on_all_samples.R input_file output_prefix [bin_size]
# Example: Rscript get_genotypes_on_all_samples.R ./ancestry/bzea_bin_genotypes.tsv all_samples 1000000

suppressPackageStartupMessages({
  library(data.table) # big table management
  library(tidyr)   # For data manipulation
  library(dplyr)   # For data manipulation
  library(ggplot2) # For visualization
  library(Ckmeans.1d.dp) # unidimensional K means clustering
  library(rebmix)  # For Gaussian mixture modeling
})

source("./R/hmm.R")
source("./R/ancestry_segments.R")

# Process command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript get_genotypes_on_all_samples.R input_file output_prefix [bin_size]\n")
  cat("Example: Rscript get_genotypes_on_all_samples.R ./ancestry/bzea_bin_genotypes.tsv  all_samples 1000000\n")
  quit(status = 1)
}

input_file <- args[1]
output_prefix <- args[2]
bin_size <- ifelse(length(args) >= 3, as.numeric(args[3]), 1000000)

# Create output directory if needed
output_dir <- "ancestry"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set up logging
log_file <- file.path(output_dir, paste0(output_prefix, ".log"))
file.create("log_file")
log_con <- file(log_file, "w")
sink(log_con, type = "output", append = FALSE)
sink(log_con, type = "message", append = FALSE)

cat(paste("Processing file:", input_file, "\n"))
cat(paste("Output prefix:", output_prefix, "\n"))
cat(paste("Bin size:", bin_size, "bp\n"))

# Define correct chromosome order
chrom_order <- paste0("chr", 1:10)

# Start timing
start_time <- Sys.time()
cat(paste("Started at:", start_time, "\n"))

# Read the binned allelic counts data from all samples
cat("Reading allelic counts data...\n")
read_freq <- try(read.table(input_file, comment.char = "@", header = TRUE))

if (inherits(read_freq, "try-error")) {
  cat("Error reading input file. Please check the file format.\n")
  quit(status = 1)
}

# Ensure CONTIG stays as a factor with proper ordering throughout the analysis
read_freq$CONTIG <- factor(read_freq$CONTIG, levels = chrom_order)

# Basic statistics
cat(paste("Total number of bins:", nrow(read_freq), "\n"))
cat(paste("Number of samples:", length(unique(read_freq$SAMPLE)), "\n"))
cat(paste("Average number of informative variants per bin:",
          round(mean(read_freq$INFORMATIVE_VARIANT_COUNT), 1), "\n"))
cat(paste("Median number of informative variants per bin:", 
          median(read_freq$INFORMATIVE_VARIANT_COUNT), "\n"))
cat(paste("Range of informative variants per bin:", 
          min(read_freq$INFORMATIVE_VARIANT_COUNT), "to", 
          max(read_freq$INFORMATIVE_VARIANT_COUNT), "\n"))

# Identify zero-ALT bins
zero_alt_bins <- sum(read_freq$ALT_FREQ == 0)
total_bins <- nrow(read_freq)

cat(paste("Bins with ALT_FREQ = 0:",
          zero_alt_bins, "out of", total_bins,
          "(", round(100 * zero_alt_bins / total_bins, 1), "%)\n"))

# Create a separate dataset for non-zero bins (for clustering)
non_zero_bins <- read_freq %>%
  filter(ALT_FREQ > 0) %>%
  # Maintain chromosome ordering
  mutate(CONTIG = factor(CONTIG, levels = chrom_order))

cat(paste("Bins with ALT_FREQ > 0 (to be clustered):",
          nrow(non_zero_bins), "out of", total_bins,
          "(", round(100 * nrow(non_zero_bins) / total_bins, 1), "%)\n"))

# Create arcsin transformation for better clustering
non_zero_bins$asin <- asin(sqrt(non_zero_bins$ALT_FREQ))

# Function to relabel clusters as REF, HET, ALT based on ALT_FREQ
relabel_clusters <- function(clusters, data) {
  # Calculate mean ALT_FREQ for each cluster
  cluster_means <- tapply(data$ALT_FREQ, clusters, mean)
  
  # Order clusters by mean ALT_FREQ
  ordered_clusters <- order(cluster_means)
  
  # Create mapping from original cluster to REF, HET, ALT
  cluster_map <- rep(NA, length(unique(as.numeric(clusters))))
  cluster_map[ordered_clusters[1]] <- "REF"
  cluster_map[ordered_clusters[2]] <- "HET"
  cluster_map[ordered_clusters[3]] <- "ALT"
  
  # Apply mapping to original clusters
  return(factor(cluster_map[as.numeric(clusters)], levels = c("REF", "HET", "ALT")))
}

# Apply K-means clustering with k=3 to the non-zero bins
K <- 3

cat("Applying K-means clustering to non-zero bins...\n")
if (nrow(non_zero_bins) > 0) {
  # Handle error if too few unique values for clustering
  tryCatch({
    # Apply K-means directly on ALT_FREQ
    non_zero_bins$Kraw <- as.factor(Ckmeans.1d.dp(non_zero_bins$ALT_FREQ, K)$cluster)
    
    # Apply K-means on arcsin-transformed data
    non_zero_bins$Kasin <- as.factor(Ckmeans.1d.dp(non_zero_bins$asin, K)$cluster)
    
    # Apply relabeling to the raw clustering
    non_zero_bins$Kraw <- relabel_clusters(non_zero_bins$Kraw, non_zero_bins)
    
    # Apply relabeling to the arcsin clustering
    non_zero_bins$Kasin <- relabel_clusters(non_zero_bins$Kasin, non_zero_bins)
    
    # Check cluster proportions for raw K-means
    raw_props <- prop.table(table(non_zero_bins$Kraw)) * 100
    cat("K-means on raw ALT_FREQ cluster proportions (non-zero bins):\n")
    print(round(raw_props, 3))
    
    # Check cluster proportions for arcsin K-means
    asin_props <- prop.table(table(non_zero_bins$Kasin)) * 100
    cat("K-means on arcsin-transformed ALT_FREQ cluster proportions (non-zero bins):\n")
    print(round(asin_props, 3))
    
  }, error = function(e) {
    cat("Error in clustering:", e$message, "\n")
    cat("Using simple thresholding instead.\n")
    
    # Apply simple thresholding if clustering fails
    threshold_result <- factor(
      ifelse(non_zero_bins$ALT_FREQ < 0.3, "REF",
             ifelse(non_zero_bins$ALT_FREQ > 0.7, "ALT", "HET")),
      levels = c("REF", "HET", "ALT")
    )
    
    non_zero_bins$Kraw <- threshold_result
    non_zero_bins$Kasin <- threshold_result
  })
} else {
  cat("No non-zero bins found for clustering. All bins will be assigned as REF.\n")
}

# Apply Gaussian Mixture Modeling
cat("Applying Gaussian Mixture Model to non-zero bins...\n")
if (nrow(non_zero_bins) > 0) {
  tryCatch({
    # Fit Gaussian mixture model with 3 components
    normalest <- REBMIX(
      # Dataset = list(data.frame(Value = non_zero_bins$ALT_FREQ)),
      Dataset = list(data.frame(Value = non_zero_bins$asin)),
      Preprocessing = "histogram",
      cmin = 3,  # Force 3 components
      cmax = 3,
      Criterion = "BIC",
      pdf = "normal"
    )
    
    # Assign clusters
    normclu <- RCLRMIX(x = normalest)
    non_zero_bins$Kgmm <- as.factor(normclu@Zp)
    
    # Relabel Kgmm clusters
    non_zero_bins$Kgmm <- relabel_clusters(non_zero_bins$Kgmm, non_zero_bins)
    
    # Check Kgmm cluster proportions
    cat("Gaussian Mixture Model cluster proportions (non-zero bins):\n")
    table(non_zero_bins$Kgmm) %>% prop.table() %>% round(3) %>% print()
  }, error = function(e) {
    cat("Error in Gaussian Mixture Model:", e$message, "\n")
    cat("Using simple thresholding instead.\n")
    
    # Apply simple thresholding if GMM fails
    non_zero_bins$Kgmm <- factor(
      ifelse(non_zero_bins$ALT_FREQ < 0.3, "REF",
             ifelse(non_zero_bins$ALT_FREQ > 0.7, "ALT", "HET")),
      levels = c("REF", "HET", "ALT")
    )
  })
} else {
  cat("No non-zero bins found for clustering. All bins will be assigned as REF.\n")
}

# Merge the clustering results back to the full dataset
cat("Merging results back to full dataset...\n")
# Create new columns in read_freq for each clustering method
read_freq$Kraw <- factor(NA, levels = c("REF", "HET", "ALT"))
read_freq$Kasin <- factor(NA, levels = c("REF", "HET", "ALT"))
read_freq$Kgmm <- factor(NA, levels = c("REF", "HET", "ALT"))

# If there are non-zero bins to cluster
if (nrow(non_zero_bins) > 0) {
  # Use data.table for faster merging
  read_freq <- as.data.table(read_freq)
  non_zero_bins <- as.data.table(non_zero_bins)
  
  # Ensure CONTIG is still a factor with proper ordering
  read_freq$CONTIG <- factor(read_freq$CONTIG, levels = chrom_order)
  non_zero_bins$CONTIG <- factor(non_zero_bins$CONTIG, levels = chrom_order)
  
  # Set keys for merging - INCLUDE SAMPLE HERE
  setkey(read_freq, SAMPLE, CONTIG, BIN_POS)
  setkey(non_zero_bins, SAMPLE, CONTIG, BIN_POS)
  
  # Update read_freq with genotype from non_zero_bins
  read_freq[non_zero_bins, Kraw := i.Kraw]
  read_freq[non_zero_bins, Kasin := i.Kasin]
  read_freq[non_zero_bins, Kgmm := i.Kgmm]
}

# Force bins with ALT_FREQ = 0 to be "REF"
read_freq[ALT_FREQ == 0, Kraw := "REF"]
read_freq[ALT_FREQ == 0, Kasin := "REF"]
read_freq[ALT_FREQ == 0, Kgmm := "REF"]

# Check final genotype proportions for each method
cat("Final cluster proportions for K-means on raw ALT_FREQ:\n")
table(read_freq$Kraw) %>% prop.table() %>% round(3) %>% print()

cat("Final cluster proportions for K-means on arcsin-transformed ALT_FREQ:\n")
table(read_freq$Kasin) %>% prop.table() %>% round(3) %>% print()

cat("Final cluster proportions for Gaussian Mixture Model:\n")
table(read_freq$Kgmm) %>% prop.table() %>% round(3) %>% print()

# Add HMM smoothing based on BC2S3 expectations
cat("Implementing HMM smoothing based on BC2S3 expectations...\n")

# Calculate expected frequencies for BC2S3
bc2s3 <- nil_frequencies_for_hmm(bc=2, s=3)
cat("Expected genotype frequencies for BC2S3:\n")
print(round(bc2s3,2))

# Apply HMM smoothing to each sample and chromosome separately
cat("Applying HMM smoothing to each sample and chromosome...\n")
read_freq <- as.data.frame(read_freq) %>%
  group_by(SAMPLE, CONTIG) %>%
  mutate(
    Kraw_HMM = smooth_ancestry_with_hmm(Kraw, transitions=c(0.995, 0.005)),
    Kasin_HMM = smooth_ancestry_with_hmm(Kasin, transitions=c(0.995, 0.005)),
    Kgmm_HMM = smooth_ancestry_with_hmm(Kgmm, transitions=c(0.995, 0.005))
  ) %>%
  ungroup()

# Make Kraw_HMM the primary genotype call (the best method)
read_freq$GENOTYPE <- read_freq$Kraw_HMM

# Ensure CONTIG is still ordered correctly
read_freq$CONTIG <- factor(read_freq$CONTIG, levels = chrom_order)

# Calculate and display genotype proportions before and after HMM
genotype_counts <- data.frame(
  Method = c("K-means (raw)", "K-means (raw) HMM",
             "K-means (arcsin)", "K-means (arcsin) HMM", 
             "Gaussian Mixture", "Gaussian Mixture HMM",
             "Expected BC2S3"),
  REF = c(
    sum(read_freq$Kraw == "REF", na.rm = TRUE) / nrow(read_freq),
    sum(read_freq$Kraw_HMM == "REF", na.rm = TRUE) / nrow(read_freq),
    sum(read_freq$Kasin == "REF", na.rm = TRUE) / nrow(read_freq),
    sum(read_freq$Kasin_HMM == "REF", na.rm = TRUE) / nrow(read_freq),
    sum(read_freq$Kgmm == "REF", na.rm = TRUE) / nrow(read_freq),
    sum(read_freq$Kgmm_HMM == "REF", na.rm = TRUE) / nrow(read_freq),
    bc2s3["REF"]
  ),
  HET = c(
    sum(read_freq$Kraw == "HET", na.rm = TRUE) / nrow(read_freq),
    sum(read_freq$Kraw_HMM == "HET", na.rm = TRUE) / nrow(read_freq),
    sum(read_freq$Kasin == "HET", na.rm = TRUE) / nrow(read_freq),
    sum(read_freq$Kasin_HMM == "HET", na.rm = TRUE) / nrow(read_freq),
    sum(read_freq$Kgmm == "HET", na.rm = TRUE) / nrow(read_freq),
    sum(read_freq$Kgmm_HMM == "HET", na.rm = TRUE) / nrow(read_freq),
    bc2s3["HET"]
  ),
  ALT = c(
    sum(read_freq$Kraw == "ALT", na.rm = TRUE) / nrow(read_freq),
    sum(read_freq$Kraw_HMM == "ALT", na.rm = TRUE) / nrow(read_freq),
    sum(read_freq$Kasin == "ALT", na.rm = TRUE) / nrow(read_freq),
    sum(read_freq$Kasin_HMM == "ALT", na.rm = TRUE) / nrow(read_freq),
    sum(read_freq$Kgmm == "ALT", na.rm = TRUE) / nrow(read_freq),
    sum(read_freq$Kgmm_HMM == "ALT", na.rm = TRUE) / nrow(read_freq),
    bc2s3["ALT"]
  )
) %>%
  mutate(
    NON_REF = HET + ALT,
    REF_pct = REF * 100,
    HET_pct = HET * 100,
    ALT_pct = ALT * 100,
    NON_REF_pct = NON_REF * 100
  )

genotype_counts[7,] <- c(,)
cat("Genotype proportions across methods:\n")
genotype_counts %>%
  select(Method, REF_pct, HET_pct, ALT_pct, NON_REF_pct) %>%
  print(digits = 3)

# Prepare output for each chromosome
cat("Writing results to disk...\n")

# Sort by sample, chromosome (in numeric order) and bin position
read_freq <- read_freq %>%
  # Ensure chromosome ordering once more
  mutate(CONTIG = factor(CONTIG, levels = chrom_order)) %>%
  arrange(SAMPLE, CONTIG, BIN_POS)

# Save detailed binned results
output_file <- file.path(output_dir,paste0(output_prefix, "_bin_genotypes.tsv"))
write.table(read_freq, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste("Detailed bin results saved to:", output_file, "\n"))

# Find connected components and extract non-REF segments

# Prepare bin data as BED
bin_bed <- read_freq %>%
  mutate(
    sample = SAMPLE,
    chrom = CONTIG,  # Make sure this uses the factor with proper ordering
    start = BIN_START,
    end = BIN_END,
    genotype = GENOTYPE,
    score = as.numeric(factor(GENOTYPE, levels = c("REF", "HET", "ALT"))) - 1,
    alt_freq = ALT_FREQ
  ) %>%
  # Ensure chromosome ordering again
  mutate(chrom = factor(chrom, levels = chrom_order)) %>%
  select(sample, chrom, start, end, genotype, score, alt_freq) %>%
  as.data.table()



# Save non-REF segments to BED file
bed_file <- file.path(output_dir,paste0(output_prefix, "_non_REF.bed"))

if (!is.null(non_ref_segments) && nrow(non_ref_segments) > 0) {
  # Make sure chromosomes are in the right order for output
  non_ref_segments$chrom <- factor(non_ref_segments$chrom, levels = chrom_order)
  non_ref_segments <- non_ref_segments[order(non_ref_segments$sample, non_ref_segments$chrom, non_ref_segments$start)]
  
  # Save the full non_ref_segments data for analysis
  full_segments_file <- file.path(output_dir, paste0(output_prefix, "_segments_data.tsv"))
  write.table(non_ref_segments, file = full_segments_file,
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(paste("Full segment data saved to:", full_segments_file, "\n"))
  
  # Save BED file (without header, simplified for genome browsers)
  write.table(non_ref_segments[, .(sample, chrom, start, end, genotype, score, alt_freq)], 
              file = bed_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat(paste("BED format results with non-REF segments saved to:", bed_file, "\n"))
  
  # Analyze segment length distribution
  cat("\n==== Segment Length Analysis ====\n")
  # Convert length to Mb for better readability
  non_ref_segments[, length_mb := length / 1e6]
  
  # Overall summary statistics
  cat("Overall segment length statistics (in Mb):\n")
  length_stats <- data.frame(
    Statistic = c("Min", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max", "Total Segments"),
    Value = c(
      min(non_ref_segments$length_mb),
      quantile(non_ref_segments$length_mb, 0.25),
      median(non_ref_segments$length_mb),
      mean(non_ref_segments$length_mb),
      quantile(non_ref_segments$length_mb, 0.75),
      max(non_ref_segments$length_mb),
      nrow(non_ref_segments)
    )
  )
  print(length_stats, row.names = FALSE)
  
  # Segment length by genotype
  cat("\nSegment length statistics by genotype (in Mb):\n")
  non_ref_segments[, .(
    count = .N,
    min_length = min(length_mb),
    median_length = median(length_mb),
    mean_length = mean(length_mb),
    max_length = max(length_mb),
    total_length = sum(length_mb)
  ), by = genotype] %>% print()
  
  # Generate segment length histogram
  length_hist_file <- file.path(output_dir,paste0(output_prefix, "_segment_length_hist.pdf"))
  pdf(length_hist_file, width = 8, height = 6)
  hist(non_ref_segments$length_mb, 
       breaks = 30,
       main = "Distribution of Non-REF Segment Lengths",
       xlab = "Segment Length (Mb)",
       col = "steelblue")
  abline(v = mean(non_ref_segments$length_mb), col = "red", lwd = 2)
  text(x = mean(non_ref_segments$length_mb) + 1, 
       y = max(hist(non_ref_segments$length_mb, breaks = 30, plot = FALSE)$counts) * 0.9,
       labels = paste("Mean =", round(mean(non_ref_segments$length_mb), 2), "Mb"),
       col = "red")
  dev.off()
  cat(paste("Segment length histogram saved to:", length_hist_file, "\n"))
  
  # Segment length by genotype (HET vs ALT)
  length_by_genotype_file <- file.path(output_dir,paste0(output_prefix, "_length_by_genotype.pdf"))
  pdf(length_by_genotype_file, width = 10, height = 6)
  boxplot(length_mb ~ genotype, data = non_ref_segments,
          main = "Segment Length by Genotype",
          xlab = "Genotype",
          ylab = "Segment Length (Mb)",
          col = c("springgreen4", "purple4"))
  dev.off()
  cat(paste("Segment length by genotype plot saved to:", length_by_genotype_file, "\n"))
}

# End timing
end_time <- Sys.time()
time_taken <- end_time - start_time
cat(paste("Completed at:", end_time, "\n"))
cat(paste("Time taken:", round(time_taken, 2), attr(time_taken, "units"), "\n"))

# Close log
sink(type = "output")
sink(type = "message")
close(log_con)