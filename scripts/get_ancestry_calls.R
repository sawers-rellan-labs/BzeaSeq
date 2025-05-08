#!/usr/bin/env Rscript
# BinAncestry.R
# Script to process GATK AlleleCountsOutput and generate binned ancestry calls
# Usage: Rscript BinAncestry.R input.allelicCounts.tsv output_prefix [bin_size]
# Example: Rscript BinAncestry.R sample1.allelicCounts.tsv sample1 100000

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
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

informative_variant_count <- nrow(read_count %>% filter(ALT_COUNT + REF_COUNT > 0))

# Print basic stats
cat(paste("Total variant positions in file:", nrow(read_count), "\n"))
cat(paste("Total informative variants:", informative_variant_count, "\n"))
cat(paste("Chromosomes found:", paste(unique(read_count$CONTIG), collapse=", "), "\n"))

# Create bins and calculate statistics
cat("Creating bins and calculating statistics...\n")
read_freq <- read_count %>%
  group_by(CONTIG) %>%
  # Create position bins by dividing positions into bins of size bin_size
  mutate(BIN_POS = ceiling(POSITION/bin_size) %>% as.integer(),
         READ_DEPTH = ALT_COUNT + REF_COUNT) %>% # Total read depth at this position, GATK calls it allelecount
  # Group by these bins
  group_by(CONTIG, BIN_POS) %>%
  # Calculate statistics for each bin
  summarise(
    VARIANT_COUNT = n(),                                      # Number of variants in bin  
    INFORMATIVE_VARIANT_COUNT = sum(ALT_COUNT + REF_COUNT > 0), # Number of variants with actual data for this sample
    DEPTH_SUM = sum(READ_DEPTH),                               # Sum of read depths over all variant positions in the bin
    ALT_COUNT = sum(ALT_COUNT),                               # Sum of alternative allele counts
    ALT_FREQ = ifelse(sum(READ_DEPTH) > 0, sum(ALT_COUNT)/sum(READ_DEPTH), 0), # Alternative allele frequency
    BIN_START = min(POSITION),                                # Bin start position
    BIN_END = max(POSITION),                                  # Bin end position
    .groups = 'drop'
  ) %>%
  # Convert to data.table for faster operations later
  as.data.table()
read_freq$CONTIG <- factor(read_freq$CONTIG, levels = paste0("chr",1:10))

# Display bin summary
total_bins <- nrow(read_freq)
cat(paste("Total number of bins:", total_bins, "\n"))
cat(paste("Average number of informative variants per bin:", 
          round(mean(read_freq$INFORMATIVE_VARIANT_COUNT), 1), "\n"))
cat(paste("Median number of informative variants per bin:", median(read_freq$INFORMATIVE_VARIANT_COUNT), "\n"))
cat(paste("Range of informative variants per bin:", min(read_freq$INFORMATIVE_VARIANT_COUNT), "to", 
          max(read_freq$INFORMATIVE_VARIANT_COUNT), "\n"))

# Identify zero-ALT bins
zero_alt_bins <- sum(read_freq$ALT_FREQ == 0)
cat(paste("Bins with ALT_FREQ = 0:", 
          zero_alt_bins, "out of", total_bins, 
          "(", round(100 * zero_alt_bins / total_bins, 1), "%)\n"))

# Create a separate dataset for non-zero bins (for clustering)
non_zero_bins <- read_freq %>%
  filter(ALT_FREQ > 0) %>%
  as.data.table()

cat(paste("Bins with ALT_FREQ > 0 (to be clustered):", 
          nrow(non_zero_bins), "out of", total_bins,
          "(", round(100 * nrow(non_zero_bins) / total_bins, 1), "%)\n"))

# Apply K-means clustering with k=3 to the non-zero bins
K <- 3

cat("Applying K-means clustering to non-zero bins...\n")
if (nrow(non_zero_bins) > 0) {
  # Handle error if too few unique values for clustering
  tryCatch({
    non_zero_bins$K <- as.factor(Ckmeans.1d.dp(non_zero_bins$ALT_FREQ, K)$cluster)
    non_zero_bins$Klog <- as.factor(Ckmeans.1d.dp(non_zero_bins$log, K)$cluster)
    
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
    
    # Apply relabeling to the clustering
    non_zero_bins$K <- relabel_clusters(non_zero_bins$K, non_zero_bins)
    non_zero_bins$Klog <- relabel_clusters(non_zero_bins$Klog, non_zero_bins)
    
    # Check cluster proportions
    cluster_props <- prop.table(table(non_zero_bins$Klog)) * 100
    cat("K-means cluster proportions (non-zero bins):\n")
    print(round(cluster_props, 1))
  }, error = function(e) {
    cat("Error in clustering:", e$message, "\n")
    cat("Using simple thresholding instead.\n")
    
    # Apply simple thresholding if clustering fails
    non_zero_bins$K <- factor(
      ifelse(non_zero_bins$ALT_FREQ < 0.3, "REF",
             ifelse(non_zero_bins$ALT_FREQ > 0.7, "ALT", "HET")),
      levels = c("REF", "HET", "ALT")
    )
    non_zero_bins$Klog <- non_zero_bins$K
  })
} else {
  cat("No non-zero bins found for clustering. All bins will be assigned as REF.\n")
}

# Merge the clustering results back to the full dataset
cat("Merging results back to full dataset...\n")
# Create a new column in read_freq for the clustering
read_freq$GENOTYPE <- factor(NA, levels = c("REF", "HET", "ALT"))

# If there are non-zero bins to cluster
if (nrow(non_zero_bins) > 0) {
  # Use data.table for faster merging
  setkey(non_zero_bins, CONTIG, BIN_POS)
  setkey(read_freq, CONTIG, BIN_POS)
  
  # Update read_freq with genotype from non_zero_bins
  read_freq[non_zero_bins, GENOTYPE := i.Klog]
}

# Then, force bins with ALT_FREQ = 0 to be "REF"
read_freq[ALT_FREQ == 0, GENOTYPE := "REF"]

# The log transform has nice normal data but high sensitivity to HET
# so I decide to change those hets to REF
read_freq[GENOTYPE == "HET", GENOTYPE := "REF"]

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

# Function to find connected components and extract ALT segments
find_alt_segments <- function(bed_data) {
  # Ensure we have a data.table
  if (!is.data.table(bed_data)) {
    bed_data <- as.data.table(bed_data)
  }
  
  # Verify column names
  required_cols <- c("sample", "chrom", "start", "end", "genotype")
  missing_cols <- setdiff(required_cols, colnames(bed_data))
  
  if (length(missing_cols) > 0) {
    cat("Error: Missing required columns in bed_data:", paste(missing_cols, collapse=", "), "\n")
    return(NULL)
  }
  
  bed_data$chrom <- factor(bed_data$chrom, levels = paste0("chr",1:10))
  
  # Create empty data.table for results
  all_segments <- data.table()
  
  # Process each chromosome separately
  cat("Processing chromosomes for ALT segments...\n")
  
  for (chr in levels(bed_data$chrom)) {
    cat("  Processing", chr, "...\n")
    chr_data <- bed_data[chrom == chr]
    
    # Sort by start position
    setorder(chr_data, start)
    
    # Create a state column (1 for ALT, 0 for others)
    chr_data[, state := ifelse(genotype == "ALT", 1, 0)]
    
    # Skip if no ALT segments
    if (sum(chr_data$state) == 0) {
      cat("    No ALT segments found in", chr, "\n")
      next
    }
    
    # Run length encoding to find consecutive segments with same state
    rle_result <- rle(chr_data$state)
    
    # Calculate segment boundaries
    segment_ends <- cumsum(rle_result$lengths)
    segment_starts <- c(1, segment_ends[-length(segment_ends)] + 1)
    
    # Create segments data frame
    segments <- data.table(
      state = rle_result$values,
      start_idx = segment_starts,
      end_idx = segment_ends
    )
    
    # Extract only ALT segments
    alt_segments <- segments[state == 1]
    
    # If there are ALT segments, process them
    if (nrow(alt_segments) > 0) {
      # Map indices back to positions
      for (i in 1:nrow(alt_segments)) {
        start_idx <- alt_segments$start_idx[i]
        end_idx <- alt_segments$end_idx[i]
        
        # Get the minimum start and maximum end positions for this segment
        min_start <- min(chr_data$start[start_idx:end_idx])
        max_end <- max(chr_data$end[start_idx:end_idx])
        
        # Calculate mean score for the segment if it exists
        mean_score <- 0
        if ("score" %in% colnames(chr_data)) {
          mean_score <- mean(chr_data$score[start_idx:end_idx])
        }
        
        # Calculate mean alt_freq for the segment if it exists
        mean_alt_freq <- 0
        if ("alt_freq" %in% colnames(chr_data)) {
          mean_alt_freq <- mean(chr_data$alt_freq[start_idx:end_idx])
        }
        
        # Add to results
        all_segments <- rbind(all_segments, data.table(
          sample = chr_data$sample[1],
          chrom = chr,
          start = min_start,
          end = max_end,
          genotype = "ALT",
          score = mean_score,
          alt_freq = mean_alt_freq
        ))
      }
    }
  }
  
  if (nrow(all_segments) > 0) {
    cat("Found", nrow(all_segments), "ALT segments across", 
        length(unique(all_segments$chrom)), "chromosomes.\n")
    return(all_segments)
  } else {
    cat("No ALT segments found in any chromosome.\n")
    return(NULL)
  }
}

# Prepare BED data
bed_data <- read_freq %>%
  mutate(
    sample = output_prefix,
    chrom = CONTIG,
    start = BIN_START,
    end = BIN_END,
    genotype = GENOTYPE,
    alt_freq = ALT_FREQ
  ) %>%
  select(sample, chrom, start, end, genotype, alt_freq) %>%
  as.data.table()

# Find ALT segments
alt_segments <- find_alt_segments(bed_data)

# Save ALT segments to BED file
bed_file <- paste0(output_prefix, "_ALT.bed")
if (!is.null(alt_segments) && nrow(alt_segments) > 0) {
  write.table(alt_segments, file = bed_file, 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat(paste("BED format results with ALT segments saved to:", bed_file, "\n"))
} else {
  cat("Warning: No ALT segments found. BED file not created.\n")
}

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
