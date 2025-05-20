#!/usr/bin/env Rscript
# get_genotypes_on_all_samples.R
#
# Description: Processes bin-level allelic count data from multiple samples,
# performs clustering across all samples, and applies HMM smoothing.
# Also analyzes segment lengths and transition frequencies.
#
# Usage: Rscript get_genotypes_on_all_samples.R input_file output_prefix [bin_size]
# Example: Rscript get_genotypes_on_all_samples.R ./ancestry/bzea_bin_genotypes.tsv ./ancestry/all_samples 1000000

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(Ckmeans.1d.dp)
  library(rebmix)  # For Gaussian mixture modeling
  library(HMM)     # For HMM smoothing
  library(ggplot2) # For visualization
  library(tidyr)   # For data manipulation
})

# Process command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript get_genotypes_on_all_samples.R input_file output_prefix [bin_size]\n")
  cat("Example: Rscript get_genotypes_on_all_samples.R ./ancestry/bzea_bin_genotypes.tsv ./ancestry/all_samples 1000000\n")
  quit(status = 1)
}

input_file <- args[1]
output_prefix <- args[2]
bin_size <- ifelse(length(args) >= 3, as.numeric(args[3]), 1000000)

# Create output directory if needed
output_dir <- dirname(output_prefix)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set up logging
log_file <- paste0(output_prefix, ".log")
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

# Apply K-means clustering with k=3 to the non-zero bins
K <- 3

cat("Applying K-means clustering to non-zero bins...\n")
if (nrow(non_zero_bins) > 0) {
  # Handle error if too few unique values for clustering
  tryCatch({
    non_zero_bins$Kasin <- as.factor(Ckmeans.1d.dp(non_zero_bins$asin, K)$cluster)
    
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
    non_zero_bins$Kasin <- relabel_clusters(non_zero_bins$Kasin, non_zero_bins)
    
    # Check cluster proportions
    cluster_props <- prop.table(table(non_zero_bins$Kasin)) * 100
    cat("K-means cluster proportions (non-zero bins):\n")
    print(round(cluster_props, 3))
  }, error = function(e) {
    cat("Error in clustering:", e$message, "\n")
    cat("Using simple thresholding instead.\n")
    
    # Apply simple thresholding if clustering fails
    non_zero_bins$Kasin <- factor(
      ifelse(non_zero_bins$ALT_FREQ < 0.3, "REF",
             ifelse(non_zero_bins$ALT_FREQ > 0.7, "ALT", "HET")),
      levels = c("REF", "HET", "ALT")
    )
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
      Dataset = list(data.frame(Value = non_zero_bins$ALT_FREQ)),
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
  read_freq[non_zero_bins, Kasin := i.Kasin]
  read_freq[non_zero_bins, Kgmm := i.Kgmm]
}

# Force bins with ALT_FREQ = 0 to be "REF"
read_freq[ALT_FREQ == 0, Kasin := "REF"]
read_freq[ALT_FREQ == 0, Kgmm := "REF"]

# Check final genotype proportions
cat("Final cluster proportions for K-means (arcsin):\n")
table(read_freq$Kasin) %>% prop.table() %>% round(3) %>% print()

cat("Final cluster proportions for Gaussian Mixture Model:\n")
table(read_freq$Kgmm) %>% prop.table() %>% round(3) %>% print()

# Add HMM smoothing based on BC2S3 expectations
cat("Implementing HMM smoothing based on BC2S3 expectations...\n")

# Function to create breeding matrices
create_mating_matrices <- function() {
  # Backcross matrix for crossing population with AA donor
  backcross_AA <- matrix(c(
    1, 1/2, 0,    # AA × AA → [1, 0, 0]; Aa × AA → [1/2, 1/2, 0]; aa × AA → [0, 1, 0]
    0, 1/2, 1,    
    0,   0, 0     
  ), nrow = 3, byrow = TRUE)
  
  # Backcross matrix for crossing population with aa donor
  backcross_aa <- matrix(c(
    0,   0, 0,    # AA × aa → [0, 1, 0]; Aa × aa → [0, 1/2, 1/2]; aa × aa → [0, 0, 1]
    1, 1/2, 0,    
    0, 1/2, 1     
  ), nrow = 3, byrow = TRUE)
  
  # Selfing matrix
  selfing <- matrix(c(
    1, 1/4, 0,    # AA selfed → [1, 0, 0]; Aa selfed → [1/4, 1/2, 1/4]; aa selfed → [0, 0, 1]
    0, 1/2, 0,    
    0, 1/4, 1     
  ), nrow = 3, byrow = TRUE)
  
  return(list(
    backcross_AA = backcross_AA,
    backcross_aa = backcross_aa,
    selfing = selfing
  ))
}

# Function to calculate genotype frequencies after breeding scheme
calculate_nil_frequencies <- function(bc=2, s=3, donor_type="aa", from=c(0,1,0)) {
  if (bc < 0 || s < 0) {
    stop("Number of backcrosses and self generations must be non-negative")
  }
  if (!donor_type %in% c("AA", "aa")) {
    stop("Donor type must be 'AA' or 'aa'")
  }
  
  # Get mating matrices
  matrices <- create_mating_matrices()
  
  # Initial F1 generation - all heterozygous after crossing pure lines
  current_population <- matrix(from, ncol = 1)
  
  # Apply backcrosses sequentially
  if (bc > 0) {
    backcross_matrix <- if (donor_type == "aa") {
      matrices$backcross_AA
    } else {
      matrices$backcross_aa
    }
    
    for (i in 1:bc) {
      current_population <- backcross_matrix %*% current_population
    }
  }
  
  # Apply selfing generations sequentially  
  if (s > 0) {
    for (i in 1:s) {
      current_population <- matrices$selfing %*% current_population
    }
  }
  
  # Extract final frequencies as named vector
  final_frequencies <- as.vector(current_population)
  names(final_frequencies) <- c("AA", "Aa", "aa")
  
  return(final_frequencies)
}

# Function to return frequencies in HMM format (REF, HET, ALT)
nil_frequencies_for_hmm <- function(bc, s, donor_type="aa", from=c(0,1,0)) {
  # Calculate raw genotype frequencies
  genotype_freqs <- calculate_nil_frequencies(bc, s, donor_type, from)
  
  # Convert to HMM format based on which allele is being introgressed
  if (donor_type == "AA") {
    hmm_frequencies <- c(
      REF = genotype_freqs["aa"],  # Recurrent parent background
      HET = genotype_freqs["Aa"],  # Heterozygous regions  
      ALT = genotype_freqs["AA"]   # Donor introgression
    )
  } else {
    hmm_frequencies <- c(
      REF = genotype_freqs["AA"],  # Recurrent parent background
      HET = genotype_freqs["Aa"],  # Heterozygous regions
      ALT = genotype_freqs["aa"]   # Donor introgression  
    )
  }
  
  return(hmm_frequencies)
}

# Calculate expected frequencies for BC2S3
bc2s3 <- nil_frequencies_for_hmm(bc=2, s=3)
cat("Expected genotype frequencies for BC2S3:\n")
print(round(bc2s3,2))

# Function to apply HMM smoothing to ancestry calls
smooth_ancestry_with_hmm <- function(genotypes, transitions = c(0.995, 0.005)) {
  # Convert genotypes to numeric (0=REF, 1=HET, 2=ALT)
  geno_numeric <- as.numeric(factor(genotypes, levels=c("REF", "HET", "ALT"))) - 1
  
  # Set up transition probabilities (high probability of staying in same state)
  trans_prob <- matrix(c(
    transitions[1], transitions[2]/2, transitions[2]/2,  # From REF
    transitions[2]/2, transitions[1], transitions[2]/2,  # From HET
    transitions[2]/2, transitions[2]/2, transitions[1]   # From ALT
  ), nrow=3, byrow=TRUE)
  
  # Set up emission probabilities (how likely each state produces observed genotypes)
  emiss_prob <- matrix(c(
    0.9, 0.08, 0.02,  # REF state
    0.1, 0.8, 0.1,    # HET state
    0.02, 0.08, 0.9   # ALT state
  ), nrow=3, byrow=TRUE)
  
  # Initialize HMM with BC2S3 priors
  hmm <- initHMM(c("REF", "HET", "ALT"), c("0", "1", "2"), 
                 startProbs = bc2s3,  # Prior probabilities from breeding design
                 transProbs = trans_prob, 
                 emissionProbs = emiss_prob)
  
  # Run Viterbi algorithm to find most likely state sequence
  viterbi_path <- viterbi(hmm, as.character(geno_numeric))
  
  # Convert back to genotype calls
  smoothed_genotypes <- factor(viterbi_path, levels=c("REF", "HET", "ALT"))
  
  return(smoothed_genotypes)
}

# Apply HMM smoothing to each sample and chromosome separately
cat("Applying HMM smoothing to each sample and chromosome...\n")
read_freq <- as.data.frame(read_freq) %>%
  group_by(SAMPLE, CONTIG) %>%
  mutate(
    Kasin_HMM = smooth_ancestry_with_hmm(Kasin, transitions=c(0.995, 0.005)),
    Kgmm_HMM = smooth_ancestry_with_hmm(Kgmm, transitions=c(0.995, 0.005))
  ) %>%
  ungroup()

# Make Kgmm_HMM the primary genotype call (the best method)
read_freq$GENOTYPE <- read_freq$Kgmm_HMM

# Ensure CONTIG is still ordered correctly
read_freq$CONTIG <- factor(read_freq$CONTIG, levels = chrom_order)

# Calculate and display genotype proportions before and after HMM
genotype_counts <- data.frame(
  Method = c("K-means (arcsin)", "K-means (arcsin) HMM", 
             "Gaussian Mixture", "Gaussian Mixture HMM",
             "Expected BC2S3"),
  REF = c(
    sum(read_freq$Kasin == "REF", na.rm = TRUE) / nrow(read_freq),
    sum(read_freq$Kasin_HMM == "REF", na.rm = TRUE) / nrow(read_freq),
    sum(read_freq$Kgmm == "REF", na.rm = TRUE) / nrow(read_freq),
    sum(read_freq$Kgmm_HMM == "REF", na.rm = TRUE) / nrow(read_freq),
    bc2s3["REF"]
  ),
  HET = c(
    sum(read_freq$Kasin == "HET", na.rm = TRUE) / nrow(read_freq),
    sum(read_freq$Kasin_HMM == "HET", na.rm = TRUE) / nrow(read_freq),
    sum(read_freq$Kgmm == "HET", na.rm = TRUE) / nrow(read_freq),
    sum(read_freq$Kgmm_HMM == "HET", na.rm = TRUE) / nrow(read_freq),
    bc2s3["HET"]
  ),
  ALT = c(
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
output_file <- paste0(output_prefix, "_bin_genotypes.tsv")
write.table(read_freq, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste("Detailed bin results saved to:", output_file, "\n"))

# Function to find connected components and extract non-REF segments
find_non_ref_segments <- function(bed_data) {
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
      
      # Create a state column (1 for non-REF, 0 for REF)
      chr_data[, state := ifelse(genotype != "REF", 1, 0)]
      
      # Skip if no non-REF segments
      if (sum(chr_data$state) == 0) {
        cat("    No non-REF segments found in", chr, "for sample", samp, "\n")
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
      
      # Extract only non-REF segments
      non_ref_segments <- segments[state == 1]
      
      # If there are non-REF segments, process them
      if (nrow(non_ref_segments) > 0) {
        # Map indices back to positions
        for (i in 1:nrow(non_ref_segments)) {
          start_idx <- non_ref_segments$start_idx[i]
          end_idx <- non_ref_segments$end_idx[i]
          
          # Get the minimum start and maximum end positions for this segment
          min_start <- min(chr_data$start[start_idx:end_idx])
          max_end <- max(chr_data$end[start_idx:end_idx])
          
          # Get the most common genotype in this segment (HET or ALT)
          genotypes <- chr_data$genotype[start_idx:end_idx]
          genotype_counts <- table(genotypes)
          dominant_genotype <- names(genotype_counts)[which.max(genotype_counts)]
          
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
            length = max_end - min_start + 1,
            genotype = dominant_genotype,
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
    cat("\nSummary of non-REF segments found:\n")
    for (i in 1:nrow(sample_counts)) {
      cat("Sample", sample_counts$sample[i], ":", sample_counts$N[i], "segments\n")
    }
    
    cat("Total:", nrow(all_segments), "non-REF segments across", 
        length(unique(all_segments$sample)), "samples and",
        length(unique(all_segments$chrom)), "chromosomes.\n")
    
    return(all_segments)
  } else {
    cat("No non-REF segments found in any sample or chromosome.\n")
    return(NULL)
  }
}

# Prepare BED data
bed_data <- read_freq %>%
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

# Find non-REF segments
non_ref_segments <- find_non_ref_segments(bed_data)

# Save non-REF segments to BED file
bed_file <- paste0(output_prefix, "_non_REF.bed")
if (!is.null(non_ref_segments) && nrow(non_ref_segments) > 0) {
  # Make sure chromosomes are in the right order for output
  non_ref_segments$chrom <- factor(non_ref_segments$chrom, levels = chrom_order)
  non_ref_segments <- non_ref_segments[order(non_ref_segments$sample, non_ref_segments$chrom, non_ref_segments$start)]
  
  # Save the full non_ref_segments data for analysis
  full_segments_file <- paste0(output_prefix, "_segments_data.tsv")
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
  length_hist_file <- paste0(output_prefix, "_segment_length_hist.pdf")
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
  length_by_genotype_file <- paste0(output_prefix, "_length_by_genotype.pdf")
  pdf(length_by_genotype_file, width = 10, height = 6)
  boxplot(length_mb ~ genotype, data = non_ref_segments,
          main = "Segment Length by Genotype",
          xlab = "Genotype",
          ylab = "Segment Length (Mb)",
          col = c("springgreen3", "purple4"))
  dev.off()
  cat(paste("Segment length by genotype plot saved to:", length_by_genotype_file, "\n"))
  
  # Analyze state transitions
  cat("\n==== State Transition Analysis ====\n")
  
  # Pre-allocate a data frame for all transitions
  transition_data <- data.frame()
  
  # Process each sample separately
  for (samp in unique(read_freq$SAMPLE)) {
    # Process each chromosome separately
    for (chr in chrom_order) {
      # Get states for this chromosome and sample
      chr_states <- read_freq %>% 
        filter(SAMPLE == samp, CONTIG == chr) %>% 
        arrange(BIN_POS) %>% 
        pull(GENOTYPE)
      
      # Skip if empty
      if (length(chr_states) <= 1) next
      
      # Find transitions - positions where state changes
      for (i in 1:(length(chr_states)-1)) {
        if (chr_states[i] != chr_states[i+1]) {
          transition_data <- rbind(transition_data, data.frame(
            sample = samp,
            chromosome = chr,
            position = i,  # Position in bin sequence
            from_state = chr_states[i],
            to_state = chr_states[i+1]
          ))
        }
      }
    }
  }
  
  # Analyze transitions if any found
  if (nrow(transition_data) > 0) {
    # Count transition types
    transition_counts <- transition_data %>%
      group_by(from_state, to_state) %>%
      summarize(
        count = n(),
        percentage = n() / nrow(transition_data) * 100,
        .groups = "drop"
      ) %>%
      mutate(transition = paste(from_state, "→", to_state)) %>%
      arrange(desc(count))
    
    cat("State transition counts:\n")
    print(transition_counts, row.names = FALSE)
    
    # Calculate transition rate per 100Mb
    # Assuming BIN_SIZE is in bp, convert to Mb and calculate per 100Mb
    total_mb <- nrow(read_freq) * bin_size / 1e6
    
    cat("\nTransition rates per 100Mb:\n")
    transition_rates <- transition_counts %>%
      mutate(
        rate_per_100mb = count / total_mb * 100
      ) %>%
      select(transition, count, rate_per_100mb)
    
    print(transition_rates, row.names = FALSE)
    
    # Calculate segment length statistics 
    # A segment is defined as consecutive bins with the same state
    
    segment_summary <- data.frame()
    
    for (samp in unique(read_freq$SAMPLE)) {
      # Process each chromosome separately
      for (chr in chrom_order) {
        # Get states for this chromosome and sample
        chr_data <- read_freq %>% 
          filter(SAMPLE == samp, CONTIG == chr) %>% 
          arrange(BIN_POS)
        
        # Skip if empty
        if (nrow(chr_data) == 0) next
        
        # Get the state sequence
        chr_states <- chr_data$GENOTYPE
        
        # Use run-length encoding to find segments
        rle_result <- rle(as.character(chr_states))
        
        # Get the values (states) and lengths
        states <- rle_result$values
        lengths <- rle_result$lengths
        
        # Record all segments
        for (i in 1:length(states)) {
          segment_summary <- rbind(segment_summary, data.frame(
            sample = samp,
            chromosome = chr,
            state = states[i],
            length_bins = lengths[i],
            length_mb = lengths[i] * (bin_size / 1e6) # Convert bin count to Mb
          ))
        }
      }
    }
    
    # Calculate summary statistics for each state
    cat("\nSegment length statistics by state (in Mb):\n")
    segment_summary %>%
      group_by(state) %>%
      summarize(
        count = n(),
        min_length = min(length_mb),
        q1_length = quantile(length_mb, 0.25),
        median_length = median(length_mb),
        mean_length = mean(length_mb),
        q3_length = quantile(length_mb, 0.75),
        max_length = max(length_mb),
        total_mb = sum(length_mb)
      ) %>%
      arrange(state) %>%
      print()
    
    # Plot transition count barplot
    transition_plot_file <- paste0(output_prefix, "_transitions.pdf")
    pdf(transition_plot_file, width = 8, height = 6)
    # Create more readable labels
    transition_counts$transition <- factor(transition_counts$transition, 
                                           levels = transition_counts$transition[order(transition_counts$count, decreasing = TRUE)])
    # Color coding
    transition_colors <- c("REF → HET" = "orange", "HET → REF" = "orange", 
                           "HET → ALT" = "purple", "ALT → HET" = "purple",
                           "REF → ALT" = "red", "ALT → REF" = "red")
    barplot(transition_counts$count, names.arg = transition_counts$transition,
            main = "State Transition Counts",
            xlab = "Transition Type",
            ylab = "Count",
            col = transition_colors[as.character(transition_counts$transition)],
            las = 2)  # Rotate labels
    dev.off()
    cat(paste("Transition count plot saved to:", transition_plot_file, "\n"))
  } else {
    cat("No state transitions found in the data.\n")
  }
} else {
  cat("Warning: No non-REF segments found. BED file not created.\n")
}

# Create a chromosome-level summary by sample
chrom_summary <- read_freq %>%
  # Ensure chromosome ordering
  mutate(CONTIG = factor(CONTIG, levels = chrom_order)) %>%
  group_by(SAMPLE, CONTIG) %>%
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
  arrange(SAMPLE, CONTIG)

summary_file <- paste0(output_prefix, "_chromosome_summary.tsv")
write.table(chrom_summary, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste("Chromosome-level summary saved to:", summary_file, "\n"))

# Calculate overall statistics per sample
overall_summary <- read_freq %>%
  group_by(SAMPLE) %>%
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
  # Sort by non-REF percentage descending
  arrange(desc(non_ref_pct))

overall_file <- paste0(output_prefix, "_summary.tsv")
write.table(overall_summary, overall_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste("Overall summary saved to:", overall_file, "\n"))

# Generate a summary plot
summary_plot_file <- paste0(output_prefix, "_summary_plot.pdf")
pdf(summary_plot_file, width = 12, height = 8)
# Create a stacked barplot of REF/HET/ALT percentages by sample
overall_summary_plot <- overall_summary %>%
  select(SAMPLE, REF_pct, HET_pct, ALT_pct) %>%
  # Reshape for ggplot
  pivot_longer(cols = c(REF_pct, HET_pct, ALT_pct), names_to = "Genotype", values_to = "Percentage") %>%
  # Clean up genotype names
  mutate(Genotype = gsub("_pct", "", Genotype)) %>%
  # Order samples by non-REF percentage
  mutate(SAMPLE = factor(SAMPLE, levels = overall_summary$SAMPLE))

# Create the plot
ggplot(overall_summary_plot, aes(x = SAMPLE, y = Percentage, fill = Genotype)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("REF" = "gold", "HET" = "springgreen4", "ALT" = "purple4")) +
  labs(
    title = "Ancestry Proportions by Sample",
    x = "Sample",
    y = "Percentage (%)",
    fill = "Genotype"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )
dev.off()
cat(paste("Summary plot saved to:", summary_plot_file, "\n"))

# End timing
end_time <- Sys.time()
time_taken <- end_time - start_time
cat(paste("Completed at:", end_time, "\n"))
cat(paste("Time taken:", round(time_taken, 2), attr(time_taken, "units"), "\n"))

# Close log
sink(type = "output")
sink(type = "message")
close(log_con)