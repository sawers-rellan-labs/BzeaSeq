#!/usr/bin/env Rscript
# get_genotypes_on_all_samples.R
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(Ckmeans.1d.dp)
})

input_file <- "./ancestry/bzea_bin_genotypes.tsv"
output_prefix <-"./ancestry/all_samples"
cat("Reading allelic counts data...\n")
read_freq <- try(read.table(input_file, comment.char = "@", header = TRUE))

if (inherits(read_freq, "try-error")) {
  cat("Error reading input file. Please check the file format.\n")
  quit(status = 1)
}

cat(paste("Processing file:", input_file, "\n"))
cat(paste("Output prefix:", output_prefix, "\n"))

# Define correct chromosome order
chrom_order <- paste0("chr", 1:10)

# Start timing
start_time <- Sys.time()
cat(paste("Started at:", start_time, "\n"))


# Ensure CONTIG stays as a factor with proper ordering throughout the analysis
read_freq$CONTIG <- factor(read_freq$CONTIG, levels = chrom_order)

# Display bin summary
total_bins <- nrow(read_freq)
cat(paste("Total number of bins:", total_bins, "\n"))
cat(paste("Average number of informative variants per bin:",
          round(mean(read_freq$INFORMATIVE_VARIANT_COUNT), 1), "\n"))
cat(paste("Median number of informative variants per bin:", median(read_freq$INFORMATIVE_VARIANT_COUNT), "\n"))
cat(paste("Range of informative variants per bin:", min(read_freq$INFORMATIVE_VARIANT_COUNT), "to", max(read_freq$INFORMATIVE_VARIANT_COUNT), "\n"))

# Identify zero-ALT bins
zero_alt_bins <- sum(read_freq$ALT_FREQ == 0)
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
  read_freq <- as.data.table(read_freq)
  non_zero_bins <- as.data.table(non_zero_bins)
  
  # Ensure CONTIG is still a factor with proper ordering
  read_freq$CONTIG <- factor(read_freq$CONTIG, levels = chrom_order)
  non_zero_bins$CONTIG <- factor(non_zero_bins$CONTIG, levels = chrom_order)
  
  # Set keys for merging - INCLUDE SAMPLE HERE
  setkey(read_freq, SAMPLE, CONTIG, BIN_POS)
  setkey(non_zero_bins, SAMPLE, CONTIG, BIN_POS)
  
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

# Make sure CONTIG is still ordered correctly
read_freq$CONTIG <- factor(read_freq$CONTIG, levels = chrom_order)

# Prepare output for each chromosome
cat("Writing results to disk...\n")

# Sort by sample, chromosome (in numeric order) and bin position
read_freq <- read_freq %>%
  # Ensure chromosome ordering once more
  mutate(CONTIG = factor(CONTIG, levels = chrom_order)) %>%
  arrange(SAMPLE,CONTIG, BIN_POS)

# Save detailed binned results
output_file <- paste0(output_prefix,"_bzea_bin_genotypes.tsv")
write.table(read_freq, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste("Detailed bin results saved to:", output_file, "\n"))

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
  
  # Ensure chromosome is properly ordered
  bed_data$chrom <- factor(bed_data$chrom, levels = chrom_order)
  
  # Create empty data.table for results
  all_segments <- data.table()
  
  # Process each sample separately
  for (samp in unique(bed_data$sample)) {
    cat("Processing sample:", samp, "\n")
    
    # Get data for this sample only
    sample_data <- bed_data[sample == samp]
    
    # Process each chromosome separately within this sample
    for (chr in chrom_order) {
      # Skip chromosomes not in the data for this sample
      if (!(chr %in% sample_data$chrom)) {
        next
      }
      
      cat("  Processing", chr, "for sample", samp, "...\n")
      chr_data <- sample_data[chrom == chr]
      
      # Sort by start position within chromosome
      setorder(chr_data, start)
      
      # Create a state column (1 for ALT, 0 for others)
      chr_data[, state := ifelse(genotype == "ALT", 1, 0)]
      
      # Skip if no ALT segments
      if (sum(chr_data$state) == 0) {
        cat("    No ALT segments found in", chr, "for sample", samp, "\n")
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
            sample = samp,
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
  }
  
  # Make sure chromosomes are still properly ordered in results
  if (nrow(all_segments) > 0) {
    all_segments$chrom <- factor(all_segments$chrom, levels = chrom_order)
    # Order by sample, chromosome, and start position
    setorder(all_segments, sample, chrom, start)
    
    # Count segments per sample
    sample_counts <- all_segments[, .N, by = sample]
    cat("\nSummary of ALT segments found:\n")
    for (i in 1:nrow(sample_counts)) {
      cat("Sample", sample_counts$sample[i], ":", sample_counts$N[i], "segments\n")
    }
    
    cat("Total:", nrow(all_segments), "ALT segments across", 
        length(unique(all_segments$sample)), "samples and",
        length(unique(all_segments$chrom)), "chromosomes.\n")
    
    return(all_segments)
  } else {
    cat("No ALT segments found in any sample or chromosome.\n")
    return(NULL)
  }
}


# Find ALT segments
alt_segments <- find_alt_segments(bed_data)

# Save ALT segments to BED file
bed_file <- paste0(output_prefix,"_bzea_ALT.bed")
if (!is.null(alt_segments) && nrow(alt_segments) > 0) {
  # Make sure chromosomes are in the right order for output
  alt_segments$chrom <- factor(alt_segments$chrom, levels = chrom_order)
  alt_segments <- alt_segments[order(alt_segments$chrom, alt_segments$start)]
  
  write.table(alt_segments, file = bed_file,
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat(paste("BED format results with ALT segments saved to:", bed_file, "\n"))
} else {
  cat("Warning: No ALT segments found. BED file not created.\n")
}

# Create a chromosome-level summary
chrom_summary <- read_freq %>%
  # Ensure chromosome ordering
  mutate(CONTIG = factor(CONTIG, levels = chrom_order)) %>%
  group_by(SAMPLE,CONTIG) %>%
  summarize(
    total_bins = n(),
    REF_bins = sum(GENOTYPE == "REF", na.rm = TRUE),
    HET_bins = sum(GENOTYPE == "HET", na.rm = TRUE),
    ALT_bins = sum(GENOTYPE == "ALT", na.rm = TRUE),
    REF_pct = 100 * sum(GENOTYPE == "REF", na.rm = TRUE) / n(),
    HET_pct = 100 * sum(GENOTYPE == "HET", na.rm = TRUE) / n(),
    ALT_pct = 100 * sum(GENOTYPE == "ALT", na.rm = TRUE) / n(),
    non_ref_pct = 100 * sum(GENOTYPE %in% c("HET", "ALT"), na.rm = TRUE) / n()
  ) %>%
  # Sort by chromosome number
  arrange(CONTIG)

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
