library(data.table)
library(dplyr)
get_ancestry_segments <- function(bed_data, chrom_order, report=FALSE) {
  # Validate input data
  validated_data <- validate_bed_data(bed_data, chrom_order)
  if (is.null(validated_data)) {
    return(NULL)
  }
  
  # Process data and extract segments
  all_segments <- extract_all_segments(validated_data, chrom_order)
  
  return(all_segments)
}

# Helper function to validate bed data
validate_bed_data <- function(bed_data, chrom_order) {
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
  
  return(bed_data)
}

# Main function to find connected components and extract ancestry segments

# Extract segments for all samples and chromosomes
extract_all_segments <- function(validated_data, chrom_order) {
  # Use dplyr to get unique sample-chromosome combinations
  sample_chroms <- validated_data %>% 
    as_tibble() %>%
    distinct(sample, chrom) %>%
    arrange(sample, chrom)
  
  # Initialize results
  all_segments <- data.table()
  #    
  
  # Process each sample-chromosome combination
  for (i in 1:nrow(sample_chroms)) {
    samp <- sample_chroms$sample[i]
    chr <- sample_chroms$chrom[i]

    # Get data for this sample and chromosome
    chr_data <- validated_data[sample == samp & chrom == chr]
    
    if(i %% 10 ==1){
      cat("Processing sample", samp, "...\n")
    }
    
    # Extract segments for all three genotypes
    for (geno in c("REF", "HET", "ALT")) {
      geno_segments <- extract_genotype_segments(chr_data, samp, chr, geno)
      all_segments <- rbind(all_segments, geno_segments)
    }
  }
  
  # Ensure chromosomes are properly ordered in results
  if (nrow(all_segments) > 0) {
    all_segments$chrom <- factor(all_segments$chrom, levels = chrom_order)
    setorder(all_segments, sample, chrom, start, genotype)
  }
  
  return(all_segments)
}

# Extract segments for a specific genotype in a chromosome
extract_genotype_segments <- function(chr_data, samp, chr, geno) {
  # Use dplyr to create state column
  chr_tbl <- as_tibble(chr_data) %>%
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
  
  # Process segments
  result_segments <- data.table()
  
  for (i in 1:nrow(geno_segments)) {
    start_idx <- geno_segments$start_idx[i]
    end_idx <- geno_segments$end_idx[i]
    
    # Get the minimum start and maximum end positions for this segment
    min_start <- min(chr_tbl$start[start_idx:end_idx])
    max_end <- max(chr_tbl$end[start_idx:end_idx])
    
    # Calculate mean statistics if available
    mean_score <- 0
    if ("score" %in% colnames(chr_tbl)) {
      mean_score <- mean(chr_tbl$score[start_idx:end_idx])
    }
    
    mean_alt_freq <- 0
    if ("alt_freq" %in% colnames(chr_tbl)) {
      mean_alt_freq <- mean(chr_tbl$alt_freq[start_idx:end_idx])
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
  
  return(result_segments)
}


data <- read.table("/Volumes/BZea/bzeaseq/ancestry/all_samples_bin_genotypes.tsv", sep = "\t", na.strings = c("NA",""), header=TRUE) %>%
  select(
    sample=SAMPLE,
    chrom=CONTIG,
    start=BIN_START,
    end=BIN_END,
    genotype=GENOTYPE
  )




