#!/usr/bin/env Rscript
# get_ancestry_calls.R
# Script to process GATK AlleleCountsOutput and generate binned ancestry calls
# Usage: Rscript get_ancestry_calls.R input.allelicCounts.tsv sample_name [bin_size]
# Example: Rscript get_ancestry_calls.R sample1.allelicCounts.tsv sample1 1000000
 
# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(Ckmeans.1d.dp)
  library(rebmix)  # For Gaussian mixture modeling
  library(HMM)     # For HMM smoothing
})

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript get_ancestry_calls.R input.allelicCounts.tsv sample_name [bin_size]\n")
  cat("Example: Rscript get_ancestry_calls.R sample1.allelicCounts.tsv sample1 1000000\n")
  quit(status = 1)
}

input_file <- args[1]
sample_name <- args[2]
bin_size <- ifelse(length(args) >= 3, as.numeric(args[3]), 1000000)

# Create output directory
output_dir <- "ancestry"
dir.create(file.path(getwd(), output_dir), showWarnings = FALSE, recursive = TRUE)

# Construct output paths
output_prefix <- file.path(output_dir, sample_name)
log_file <- paste0(output_prefix, ".log")

# Output messages to both console and log file
log_con <- file(log_file, "w")
sink(log_con, type = "output", append = FALSE)
sink(log_con, type = "message", append = FALSE)

cat(paste("Processing file:", input_file, "\n"))
cat(paste("Sample name:", sample_name, "\n"))
cat(paste("Output prefix:", output_prefix, "\n"))
cat(paste("Bin size:", bin_size, "bp\n"))

# Define correct chromosome order
chrom_order <- paste0("chr", 1:10)

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

# Ensure chromosome is a factor with correct ordering
read_count$CONTIG <- factor(read_count$CONTIG, levels = chrom_order)

informative_variant_count <- nrow(read_count %>% filter(ALT_COUNT + REF_COUNT > 0))

# Print basic stats
cat(paste("Total variant positions in file:", nrow(read_count), "\n"))
cat(paste("Total informative variants:", informative_variant_count, "\n"))
cat(paste("Chromosomes found:", paste(unique(read_count$CONTIG), collapse=", "), "\n"))


# Create bins and calculate statistics
cat("Creating bins and calculating statistics...\n")
read_freq <- read_count %>%
  # Ensure chromosome ordering
  mutate(CONTIG = factor(CONTIG, levels = chrom_order)) %>%
  group_by(CONTIG) %>%
  # Create position bins by dividing positions into bins of size bin_size
  mutate(BIN_POS = ceiling(POSITION/bin_size) %>% as.integer(),
         READ_DEPTH = ALT_COUNT + REF_COUNT) %>% # Total read depth at this position
  # Group by these bins
  group_by(CONTIG, BIN_POS) %>%
  # Calculate statistics for each bin
  summarise(
    SAMPLE = sample_name,
    VARIANT_COUNT = n(),                                     # Number of variants in bin
    INFORMATIVE_VARIANT_COUNT = sum(ALT_COUNT + REF_COUNT > 0), # Number of variants with actual data for this sample
    DEPTH_SUM = sum(READ_DEPTH),                             # Sum of read depths over all variant positions in the bin
    ALT_COUNT = sum(ALT_COUNT),                              # Sum of alternative allele counts
    ALT_FREQ = ifelse(sum(READ_DEPTH) > 0, sum(ALT_COUNT)/sum(READ_DEPTH), 0), # Alternative allele frequency
    BIN_START = min(POSITION),                               # Bin start position
    BIN_END = max(POSITION),                                 # Bin end position
    .groups = 'drop'
  ) %>% 
  select(SAMPLE, everything())

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

# Create transformations for better clustering - only arcsin transform
non_zero_bins$asin <- asin(sqrt(non_zero_bins$ALT_FREQ))

# Apply K-means clustering with k=3 for arcsin-transformed data
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
    
    # Apply relabeling to the arcsin-transformed clustering
    non_zero_bins$Kasin <- relabel_clusters(non_zero_bins$Kasin, non_zero_bins)
    
    # Check cluster proportions
    cluster_props <- prop.table(table(non_zero_bins$Kasin)) * 100
    cat("K-means (arcsin) cluster proportions (non-zero bins):\n")
    print(round(cluster_props, 3))
  }, error = function(e) {
    cat("Error in K-means clustering:", e$message, "\n")
    cat("Using simple thresholding instead.\n")
    
    # Apply simple thresholding if clustering fails
    non_zero_bins$Kasin <- factor(
      ifelse(non_zero_bins$ALT_FREQ < 0.3, "REF",
             ifelse(non_zero_bins$ALT_FREQ > 0.7, "ALT", "HET")),
      levels = c("REF", "HET", "ALT")
    )
  })
  
  # Apply Gaussian Mixture Modeling with 3 components
  tryCatch({
    cat("Applying Gaussian Mixture Model to non-zero bins...\n")
    
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
    gmm_props <- prop.table(table(non_zero_bins$Kgmm)) * 100
    cat("Gaussian Mixture Model cluster proportions (non-zero bins):\n")
    print(round(gmm_props, 3))
    
  }, error = function(e) {
    cat("Error in Gaussian Mixture Model:", e$message, "\n")
    cat("Using arcsin clustering instead for GMM.\n")
    
    # If GMM fails, use the arcsin clustering results
    non_zero_bins$Kgmm <- non_zero_bins$Kasin
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
  
  # Set keys for merging
  setkey(read_freq, CONTIG, BIN_POS)
  setkey(non_zero_bins, CONTIG, BIN_POS)
  
  # Update read_freq with genotype from non_zero_bins
  read_freq[non_zero_bins, Kasin := i.Kasin]
  read_freq[non_zero_bins, Kgmm := i.Kgmm]
}

# Force bins with ALT_FREQ = 0 to be "REF"
read_freq[ALT_FREQ == 0, Kasin := "REF"]
read_freq[ALT_FREQ == 0, Kgmm := "REF"]

# Check final genotype proportions
cat("Final Kasin genotype proportions:\n")
print(table(read_freq$Kasin) %>% prop.table() %>% round(3))

cat("Final Kgmm genotype proportions:\n")
print(table(read_freq$Kgmm) %>% prop.table() %>% round(3))

# Make sure CONTIG is still ordered correctly
read_freq$CONTIG <- factor(read_freq$CONTIG, levels = chrom_order)

# Add HMM smoothing based on BC2S3 expectations
cat("Applying HMM smoothing...\n")

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
cat("REF:", round(bc2s3["REF"], 4), "\n")
cat("HET:", round(bc2s3["HET"], 4), "\n") 
cat("ALT:", round(bc2s3["ALT"], 4), "\n")

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

# Apply HMM smoothing to each chromosome separately to respect chromosome boundaries
read_freq <- read_freq %>%
  as.data.frame() %>%
  group_by(CONTIG) %>%
  mutate(
    Kasin_HMM = smooth_ancestry_with_hmm(Kasin, transitions=c(0.995, 0.005)),
    Kgmm_HMM = smooth_ancestry_with_hmm(Kgmm, transitions=c(0.995, 0.005))
  ) %>%
  ungroup()

# Make Kgmm_HMM the primary genotype call (since it's the best method)
read_freq$GENOTYPE <- read_freq$Kgmm_HMM

# Ensure chromosome factor ordering
read_freq$CONTIG <- factor(read_freq$CONTIG, levels = chrom_order)

# Prepare output for each chromosome
cat("Writing results to disk...\n")

# Sort by chromosome (in numeric order) and bin position
read_freq <- read_freq %>%
  # Ensure chromosome ordering once more
  mutate(CONTIG = factor(CONTIG, levels = chrom_order)) %>%
  arrange(CONTIG, BIN_POS)

# Save detailed binned results
output_file <- paste0(output_prefix, "_bin_genotypes.tsv")
write.table(read_freq, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste("Detailed bin results saved to:", output_file, "\n"))

# Prepare BED data for non-REF regions (both HET and ALT)
bed_data <- read_freq %>%
  mutate(
    sample = sample_name,
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
  arrange(chrom, start) %>%  # Sort by chromosome and position
  as.data.table()

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
  
  # Process each chromosome separately in the correct order
  cat("Processing chromosomes for non-REF segments...\n")
  
  for (chr in chrom_order) {
    # Skip chromosomes not in the data
    if (!(chr %in% bed_data$chrom)) {
      next
    }
    
    cat("  Processing", chr, "...\n")
    chr_data <- bed_data[chrom == chr]
    
    # Sort by start position within chromosome
    setorder(chr_data, start)
    
    # Create a state column (1 for non-REF, 0 for REF)
    chr_data[, state := ifelse(genotype != "REF", 1, 0)]
    
    # Skip if no non-REF segments
    if (sum(chr_data$state) == 0) {
      cat("    No non-REF segments found in", chr, "\n")
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
        
        # Calculate mean score for the segment
        mean_score <- mean(chr_data$score[start_idx:end_idx])
        
        # Calculate mean alt_freq for the segment
        mean_alt_freq <- mean(chr_data$alt_freq[start_idx:end_idx])
        
        # Add to results
        all_segments <- rbind(all_segments, data.table(
          sample = chr_data$sample[1],
          chrom = chr,  # Use the actual chromosome value, not index
          start = min_start,
          end = max_end,
          genotype = dominant_genotype,  # Use the most common genotype (HET or ALT)
          score = mean_score,
          alt_freq = mean_alt_freq
        ))
      }
    }
  }
  
  # Make sure chromosomes are still properly ordered in results
  if (nrow(all_segments) > 0) {
    all_segments$chrom <- factor(all_segments$chrom, levels = chrom_order)
    setorder(all_segments, chrom, start)
    
    cat("Found", nrow(all_segments), "non-REF segments across",
        length(unique(all_segments$chrom)), "chromosomes.\n")
    return(all_segments)
  } else {
    cat("No non-REF segments found in any chromosome.\n")
    return(NULL)
  }
}

# Find non-REF segments
non_ref_segments <- find_non_ref_segments(bed_data)

# Save non-REF segments to BED file
bed_file <- paste0(output_prefix, "_non_REF.bed")
if (!is.null(non_ref_segments) && nrow(non_ref_segments) > 0) {
  # Make sure chromosomes are in the right order for output
  non_ref_segments$chrom <- factor(non_ref_segments$chrom, levels = chrom_order)
  non_ref_segments <- non_ref_segments[order(non_ref_segments$chrom, non_ref_segments$start)]
  
  write.table(non_ref_segments, file = bed_file,
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat(paste("BED format results with non-REF segments saved to:", bed_file, "\n"))
} else {
  cat("Warning: No non-REF segments found. BED file not created.\n")
}

# Create a chromosome-level summary
chrom_summary <- read_freq %>%
  # Ensure chromosome ordering
  mutate(CONTIG = factor(CONTIG, levels = chrom_order)) %>%
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
  ) %>%
  # Sort by chromosome number
  arrange(CONTIG)

summary_file <- paste0(output_prefix, "_chromosome_summary.tsv")
write.table(chrom_summary, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste("Chromosome-level summary saved to:", summary_file, "\n"))

# Calculate overall statistics
overall_summary <- data.frame(
  sample = sample_name,
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