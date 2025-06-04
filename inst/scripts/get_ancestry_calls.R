#!/usr/bin/env Rscript
# get_joint_ancestry_calls.R
#
# Updated script using BzeaSeq package functions for joint ancestry analysis
# Usage: Rscript get_joint_ancestry_calls.R input_file metadata_file output_prefix [bin_size]

suppressPackageStartupMessages({
  library(BzeaSeq)
  library(data.table)
  library(tidyr)
  library(dplyr)
})

# Process command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("Usage: Rscript get_joint_ancestry_calls.R input_file metadata_file output_prefix [bin_size]\n")
  cat("Example: Rscript get_joint_ancestry_calls.R ./ancestry/bzea_bin_genotypes.tsv sample_metadata.csv all_samples 1000000\n")
  quit(status = 1)
}

input_file <- args[1]
metadata_file <- args[2]
output_prefix <- args[3]
bin_size <- ifelse(length(args) >= 4, as.numeric(args[4]), 1000000)

# Create output directory
output_dir <- "ancestry"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set up logging
log_file <- file.path(output_dir, paste0(output_prefix, ".log"))
log_con <- file(log_file, "w")
sink(log_con, type = "output", append = FALSE)
sink(log_con, type = "message", append = FALSE)

cat(paste("Processing file:", input_file, "\n"))
cat(paste("Output prefix:", output_prefix, "\n"))
cat(paste("Bin size:", bin_size, "bp\n"))

# Define chromosome order
chrom_order <- paste0("chr", 1:10)

# Start timing
start_time <- Sys.time()
cat(paste("Started at:", start_time, "\n"))

# Read metadata and filter valid samples
metadata <- read.csv(metadata_file)
valid_samples <- filter_valid_samples(metadata)
cat(paste("Found", length(valid_samples), "valid samples\n"))

# Read bin-level data
cat("Reading allelic counts data...\n")
read_freq <- try(read.table(input_file, comment.char = "@", header = TRUE))

if (inherits(read_freq, "try-error")) {
  cat("Error reading input file. Please check the file format.\n")
  quit(status = 1)
}

# Filter to valid samples
read_freq <- read_freq %>%
  filter(SAMPLE %in% valid_samples) %>%
  mutate(CONTIG = factor(CONTIG, levels = chrom_order))

# Apply joint clustering across all samples
cat("Applying joint clustering across all samples...\n")
read_freq <- apply_joint_ancestry_clustering(read_freq)

# Calculate BC2S3 expectations
bc2s3 <- nil_frequencies_for_hmm(bc = 2, s = 3)
cat("Expected genotype frequencies for BC2S3:\n")
print(round(bc2s3, 2))

# Apply HMM smoothing per sample and chromosome
cat("Applying HMM smoothing to each sample and chromosome...\n")
read_freq <- read_freq %>%
  group_by(SAMPLE, CONTIG) %>%
  mutate(
    GENOTYPE = smooth_ancestry_with_hmm(Kgmm, transitions = c(0.995, 0.005))
  ) %>%
  ungroup()

# Convert to BED format for segment extraction
bin_bed <- read_freq %>%
  mutate(
    sample = SAMPLE,
    chrom = CONTIG,
    start = BIN_START,
    end = BIN_END,
    genotype = GENOTYPE,
    score = as.numeric(factor(GENOTYPE, levels = c("REF", "HET", "ALT"))) - 1,
    alt_freq = ALT_FREQ
  ) %>%
  select(sample, chrom, start, end, genotype, score, alt_freq)

# Extract ancestry segments
cat("Extracting ancestry segments...\n")
all_segments <- get_ancestry_segments(
  bed_data = bin_bed,
  chrom_order = chrom_order,
  report = TRUE
)

# Filter to non-REF segments
non_ref_segments <- all_segments %>%
  filter(genotype %in% c("HET", "ALT"))

# Save joint analysis results
save_joint_ancestry_results(
  bin_data = read_freq,
  segments_data = non_ref_segments,
  all_segments = all_segments,
  output_prefix = output_prefix,
  output_dir = output_dir,
  bin_size = bin_size
)

# Generate summary statistics and plots
generate_ancestry_summary_plots(
  bin_data = read_freq,
  segments_data = non_ref_segments,
  output_prefix = output_prefix,
  output_dir = output_dir
)

# End timing
end_time <- Sys.time()
time_taken <- end_time - start_time
cat(paste("Completed at:", end_time, "\n"))
cat(paste("Time taken:", round(time_taken, 2), attr(time_taken, "units"), "\n"))

# Close log
sink(type = "output")
sink(type = "message")
close(log_con)