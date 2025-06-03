#' Ancestry Segment Extraction Functions
#'
#' This module provides functions for extracting connected ancestry segments
#' from bin-level genotype data across multiple samples and chromosomes.
#'
#' @name ancestry-segment-functions
NULL

library(data.table)
library(dplyr)

#' Extract Ancestry Segments from Bin Data
#'
#' Main function to find connected components and extract ancestry segments
#' for all genotype states (REF, HET, ALT) across samples and chromosomes.
#'
#' @param bed_data Data frame or data.table with columns: sample, chrom, start, end, genotype
#' @param chrom_order Character vector. Chromosome names in desired order 
#'   (default: paste0("chr", 1:10))
#' @param report Logical. Whether to print progress reports (default: FALSE)
#'
#' @return Data.table with columns:
#'   \item{sample}{Sample identifier}
#'   \item{chrom}{Chromosome name}
#'   \item{start}{Segment start position}
#'   \item{end}{Segment end position}
#'   \item{genotype}{Segment genotype (REF, HET, or ALT)}
#'   \item{score}{Mean score if available}
#'   \item{alt_freq}{Mean alternative allele frequency if available}
#'
#' @details
#' This function processes genomic bins to identify contiguous segments
#' of the same ancestry state. Uses run-length encoding to efficiently
#' find connected components within each sample-chromosome combination.
#'
#' @examples
#' # Create example bin data
#' bin_data <- data.frame(
#'   sample = rep("sample1", 5),
#'   chrom = rep("chr1", 5),
#'   start = c(1, 1000000, 2000000, 3000000, 4000000),
#'   end = c(999999, 1999999, 2999999, 3999999, 4999999),
#'   genotype = c("REF", "ALT", "ALT", "REF", "REF")
#' )
#'
#' segments <- get_ancestry_segments(bin_data)
#' print(segments)
#'
#' @export
get_ancestry_segments <- function(bed_data, chrom_order = paste0("chr", 1:10), 
                                  report = FALSE) {
  # Validate input data
  validated_data <- validate_bed_data(bed_data, chrom_order)
  if (is.null(validated_data)) {
    return(NULL)
  }
  
  # Process data and extract segments
  all_segments <- extract_all_segments(validated_data, chrom_order, report)
  
  all_segments
}

#' Validate BED Format Data
#'
#' Internal function to validate and prepare BED format data for segment extraction.
#'
#' @param bed_data Data frame with genomic interval data
#' @param chrom_order Character vector of chromosome names in desired order
#'
#' @return Data.table with validated and formatted data, or NULL if validation fails
#'
#' @details
#' Checks for required columns, converts to data.table format, and ensures
#' proper chromosome factor ordering for consistent processing.
#'
#' @keywords internal
validate_bed_data <- function(bed_data, chrom_order) {
  # Input validation
  if (!is.data.frame(bed_data)) {
    stop("bed_data must be a data frame")
  }
  if (!is.character(chrom_order) || length(chrom_order) == 0) {
    stop("chrom_order must be a non-empty character vector")
  }
  
  # Ensure we have a data.table
  if (!is.data.table(bed_data)) {
    bed_data <- as.data.table(bed_data)
  }
  
  # Verify required columns
  required_cols <- c("sample", "chrom", "start", "end", "genotype")
  missing_cols <- setdiff(required_cols, colnames(bed_data))
  
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns in bed_data:", 
               paste(missing_cols, collapse = ", ")))
  }
  
  # Validate data types and ranges
  if (!is.numeric(bed_data$start) || !is.numeric(bed_data$end)) {
    stop("start and end columns must be numeric")
  }
  if (any(bed_data$start > bed_data$end, na.rm = TRUE)) {
    stop("start positions must be <= end positions")
  }
  if (any(bed_data$start < 0, na.rm = TRUE) || any(bed_data$end < 0, na.rm = TRUE)) {
    stop("start and end positions must be non-negative")
  }
  
  # Ensure chromosome is properly ordered
  bed_data$chrom <- factor(bed_data$chrom, levels = chrom_order)
  
  # Check for unknown chromosomes
  unknown_chroms <- unique(bed_data$chrom[is.na(bed_data$chrom)])
  if (length(unknown_chroms) > 0) {
    warning(paste("Unknown chromosomes found and will be ignored:", 
                  paste(unknown_chroms, collapse = ", ")))
    bed_data <- bed_data[!is.na(chrom)]
  }
  
  bed_data
}

#' Extract Segments for All Sample-Chromosome Combinations
#'
#' Internal function that processes each sample-chromosome combination to
#' extract ancestry segments for all genotype states.
#'
#' @param validated_data Data.table with validated BED format data
#' @param chrom_order Character vector of chromosome names in order
#' @param report Logical. Whether to print progress messages
#'
#' @return Data.table with all extracted segments
#'
#' @keywords internal
extract_all_segments <- function(validated_data, chrom_order, report = FALSE) {
  # Use dplyr to get unique sample-chromosome combinations
  sample_chroms <- validated_data %>% 
    as_tibble() %>%
    distinct(sample, chrom) %>%
    arrange(sample, chrom)
  
  if (nrow(sample_chroms) == 0) {
    warning("No valid sample-chromosome combinations found")
    return(data.table())
  }
  
  # Initialize results
  all_segments <- data.table()
  
  # Process each sample-chromosome combination
  for (i in seq_len(nrow(sample_chroms))) {
    samp <- sample_chroms$sample[i]
    chr <- sample_chroms$chrom[i]
    
    # Get data for this sample and chromosome
    chr_data <- validated_data[sample == samp & chrom == chr]
    
    if (report && (i %% 10 == 1 || i <= 10)) {
      cat("Processing sample", samp, "chromosome", chr, 
          "(", i, "of", nrow(sample_chroms), ")\n")
    }
    
    # Extract segments for all three genotypes
    for (geno in c("REF", "HET", "ALT")) {
      geno_segments <- extract_genotype_segments(chr_data, samp, chr, geno)
      if (nrow(geno_segments) > 0) {
        all_segments <- rbind(all_segments, geno_segments)
      }
    }
  }
  
  # Ensure chromosomes are properly ordered in results
  if (nrow(all_segments) > 0) {
    all_segments$chrom <- factor(all_segments$chrom, levels = chrom_order)
    setorder(all_segments, sample, chrom, start, genotype)
  }
  
  all_segments
}

#' Extract Segments for Specific Genotype in Chromosome
#'
#' Internal function to extract contiguous segments of a specific genotype
#' within a single sample-chromosome combination.
#'
#' @param chr_data Data.table with data for one sample-chromosome
#' @param samp Character. Sample identifier
#' @param chr Character. Chromosome identifier  
#' @param geno Character. Genotype to extract ("REF", "HET", or "ALT")
#'
#' @return Data.table with segments of the specified genotype
#'
#' @details
#' Uses run-length encoding to efficiently identify contiguous regions
#' of the same genotype state. Calculates segment boundaries and
#' summary statistics (mean score and alt_freq if available).
#'
#' @keywords internal
extract_genotype_segments <- function(chr_data, samp, chr, geno) {
  # Input validation
  if (!geno %in% c("REF", "HET", "ALT")) {
    stop("geno must be one of: REF, HET, ALT")
  }
  if (nrow(chr_data) == 0) {
    return(data.table())
  }
  
  # Use dplyr to create state column
  chr_tbl <- chr_data %>%
    as_tibble() %>%
    arrange(start) %>%
    mutate(state = ifelse(genotype == geno, 1, 0))
  
  # Skip if no segments of this genotype
  if (sum(chr_tbl$state) == 0) {
    return(data.table())
  }
  
  # Run length encoding to find consecutive segments with same state
  rle_result <- rle(chr_tbl$state)
  
  # Calculate segment boundaries
  segment_ends <- cumsum(rle_result$lengths)
  segment_starts <- c(1, segment_ends[-length(segment_ends)] + 1)
  
  # Create segments data frame
  segments <- data.table(
    state = rle_result$values,
    start_idx = segment_starts,
    end_idx = segment_ends
  )
  
  # Extract only segments of the current genotype
  geno_segments <- segments[state == 1]
  
  if (nrow(geno_segments) == 0) {
    return(data.table())
  }
  
  # Process segments
  result_segments <- data.table()
  
  for (i in seq_len(nrow(geno_segments))) {
    start_idx <- geno_segments$start_idx[i]
    end_idx <- geno_segments$end_idx[i]
    
    # Get the minimum start and maximum end positions for this segment
    segment_bins <- chr_tbl[start_idx:end_idx, ]
    min_start <- min(segment_bins$start)
    max_end <- max(segment_bins$end)
    
    # Calculate mean statistics if available
    mean_score <- 0
    if ("score" %in% colnames(segment_bins)) {
      mean_score <- mean(segment_bins$score, na.rm = TRUE)
      if (is.na(mean_score)) mean_score <- 0
    }
    
    mean_alt_freq <- 0
    if ("alt_freq" %in% colnames(segment_bins)) {
      mean_alt_freq <- mean(segment_bins$alt_freq, na.rm = TRUE)
      if (is.na(mean_alt_freq)) mean_alt_freq <- 0
    }
    
    # Add to results
    result_segments <- rbind(result_segments, data.table(
      sample = samp,
      chrom = chr,
      start = min_start,
      end = max_end,
      genotype = geno,
      score = mean_score,
      alt_freq = mean_alt_freq
    ))
  }
  
  result_segments
}

# Example usage (commented out to avoid execution during sourcing)
# data <- read.table("/path/to/all_samples_bin_genotypes.tsv", 
#                    sep = "\t", na.strings = c("NA",""), header = TRUE) %>%
#   select(
#     sample = SAMPLE,
#     chrom = CONTIG,
#     start = BIN_START,
#     end = BIN_END,
#     genotype = GENOTYPE
#   )