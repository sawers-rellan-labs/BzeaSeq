# introgression_plots.R
# Functions for processing and plotting introgression data
#
# Main functions:
# - process_introgression_data(): Process raw introgression data
# - prepare_stacked_data(): Prepare data for stacked visualization 
# - plot_introgression_stacked(): Create stacked bar plots by chromosome
# - convert_to_matrix(): Convert run-length encoded data to matrix for ComplexHeatmap

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# Define chromosome order as a constant
CHROMOSOME_ORDER <- paste0("chr", 1:10)

# Generate chromosome length data from assembly index
generate_chromosome_lengths <- function() {
  # Exact chromosome lengths from Zm-B73-REFERENCE-NAM-5.0.fa.fai
  chrom_lengths <- data.frame(
    chrom = paste0("chr", 1:10),
    length = c(
      308452471,  # chr1
      243675191,  # chr2
      238017767,  # chr3
      250330460,  # chr4
      226353449,  # chr5
      181357234,  # chr6
      185808916,  # chr7
      182411202,  # chr8
      163004744,  # chr9
      152435371   # chr10
    )
  )
  return(chrom_lengths)
}

# Re-encode HET and ALT segments as continuous introgression blocks
reencode_introgressions <- function(data) {
  data %>%
    # Ensure chromosome factor ordering
    mutate(chrom = factor(chrom, levels = CHROMOSOME_ORDER)) %>%
    arrange(sample, chrom, bin_start) %>%
    group_by(sample, chrom) %>%
    mutate(
      # Mark non-REF segments as introgressions
      is_introgression = genotype != "REF",
      # Create run-length encoding groups for continuous introgression blocks
      introgression_group = cumsum(is_introgression != lag(is_introgression, default = FALSE))
    ) %>%
    # Calculate introgression block statistics for each group
    group_by(sample, chrom, introgression_group) %>%
    mutate(
      # Check if this group represents an introgression block
      group_is_introgression = any(is_introgression),
      # Calculate block stats only for introgression groups
      block_start = if_else(group_is_introgression, min(bin_start), NA_real_),
      block_end = if_else(group_is_introgression, max(bin_end), NA_real_),
      block_span = if_else(group_is_introgression, block_end - block_start + 1, NA_real_),
      block_midpoint = if_else(group_is_introgression, 
                               floor((block_start + block_end) / 2), NA_real_)
    ) %>%
    # Find largest introgression block per sample-chromosome
    group_by(sample, chrom) %>%
    mutate(
      largest_introgression_midpoint = {
        # Get unique introgression blocks for this sample-chromosome
        introgression_blocks <- cur_data() %>% 
          filter(group_is_introgression) %>%
          distinct(introgression_group, block_span, block_midpoint) %>%
          filter(!is.na(block_span))
        
        if(nrow(introgression_blocks) > 0) {
          # Find the largest block
          largest_midpoint <- introgression_blocks %>% 
            slice_max(block_span, n = 1, with_ties = FALSE) %>%
            pull(block_midpoint)
          largest_midpoint[1]
        } else {
          999999999  # Large value for samples with no introgressions
        }
      }
    ) %>%
    ungroup() %>%
    # Clean up temporary columns
    select(-introgression_group, -is_introgression, -group_is_introgression,
           -block_start, -block_end, -block_span, -block_midpoint)
}

# Add segment IDs and snap to bin boundaries
add_segment_metadata <- function(data, bin_size = 1000000) {
  # Get chromosome lengths
  chrom_lengths <- generate_chromosome_lengths()
  
  # Calculate segment span and snap to bin boundaries
  result <- data %>%
    # Ensure chromosome factor ordering
    mutate(chrom = factor(chrom, levels = CHROMOSOME_ORDER)) %>%
    left_join(chrom_lengths, by = "chrom") %>%
    arrange(sample, chrom, start) %>%
    group_by(sample, chrom) %>%
    mutate(
      segment_id = row_number(),
      # Calculate bin boundaries
      bin_start = (floor((start - 1) / bin_size) * bin_size) + 1,
      bin_end = pmin(ceiling(end / bin_size) * bin_size, length),  # Don't exceed chromosome length
      # Update span based on bin boundaries
      span = bin_end - bin_start + 1
    ) %>%
    ungroup() %>%
    select(-length)  # Remove temporary chromosome length column
  
  # Re-encode introgressions and find largest introgression midpoint
  result <- reencode_introgressions(result)
  
  return(result)
}

# Main function to process introgression data
process_introgression_data <- function(data, bin_size = 1000000) {
  result <- add_segment_metadata(data, bin_size)
  return(result)
}

# Main function to prepare stacked data for visualization
prepare_stacked_data <- function(data, chrom_to_plot) {
  # Filter data for the specified chromosome and process
  result_df <- data %>% 
    # Ensure chromosome factor ordering
    mutate(chrom = factor(chrom, levels = CHROMOSOME_ORDER)) %>%
    filter(chrom == chrom_to_plot) %>%
    mutate(
      span = bin_end - bin_start + 1,  # Use bin span
      bin = row_number()  # Simple bin numbering
    ) %>%
    arrange(sample, bin_start)
  
  return(result_df)
}

# Create sample order based on ordering preference
create_sample_order <- function(stacked_data, order) {
  # If order is a character vector with length > 1, it's a custom ordering
  if (is.character(order) && length(order) > 1) {
    # Validate custom order
    sample_levels <- levels(factor(stacked_data$sample))
    
    if (length(order) != length(sample_levels)) {
      stop(paste("Custom order length (", length(order), 
                 ") must match number of unique samples (", length(sample_levels), ")", sep = ""))
    }
    
    if (!setequal(order, sample_levels)) {
      missing_in_order <- setdiff(sample_levels, order)
      extra_in_order <- setdiff(order, sample_levels)
      
      error_msg <- "Custom order must contain exactly the same samples as in the data.\n"
      if (length(missing_in_order) > 0) {
        error_msg <- paste0(error_msg, "Missing from order: ", paste(missing_in_order, collapse = ", "), "\n")
      }
      if (length(extra_in_order) > 0) {
        error_msg <- paste0(error_msg, "Extra in order: ", paste(extra_in_order, collapse = ", "))
      }
      
      stop(error_msg)
    }
    
    return(order)
  }
  
  # Handle string-based ordering
  if (is.character(order) && length(order) == 1) {
    if (order == "sample") {
      sample_order <- levels(factor(stacked_data$sample)) %>% rev()
    } else if (order == "position") {
      sample_order <- stacked_data %>%
        select(sample, largest_introgression_midpoint) %>%
        distinct() %>%
        arrange(largest_introgression_midpoint) %>%
        pull(sample)
    } else {
      stop(paste("Unknown order type:", order, ". Use 'sample', 'position', or provide a custom character vector."))
    }
    
    return(sample_order)
  }
  
  # If not character, invalid input
  stop("order must be either 'sample', 'position', or a character vector of sample names")
}

# Remove the separate validation function as it's now integrated

# Create single chromosome plot
create_chromosome_plot <- function(stacked_data, chrom, genotype_colors) {
  if (nrow(stacked_data) == 0) {
    return(ggplot() + 
             ggtitle(paste("No introgression data for", chrom)) + 
             theme_minimal())
  }
  
  # Create the plot
  p <- ggplot(stacked_data, aes(x = sample, y = span, fill = genotype, group = bin)) +
    geom_bar(stat = "identity", position = "stack", width = 1) +
    scale_fill_manual(values = genotype_colors) +
    labs(title = paste("Introgression segments in", chrom),
         x = "Sample", 
         y = "Chromosome Length") +
    scale_y_continuous(labels = function(x) paste0(round(x/1e6, 1))) +
    coord_flip() + 
    ggpubr::theme_classic2() +
    theme(
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.text.x = element_text(hjust = 0),
      strip.text.y = element_text(angle = 0, face = "bold"),
      panel.spacing.x = unit(0, "null"),
      panel.spacing.y = unit(0, "null"),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_blank(),
      axis.title.x = element_blank()
    )
  
  return(p)
}

# Main function to plot introgression segments using stacked bars
plot_introgression_stacked <- function(data, chrom_lengths = NULL, order = "sample") {
  # Generate chromosome lengths if not provided
  if (is.null(chrom_lengths)) {
    chrom_lengths <- generate_chromosome_lengths()
  }
  
  # Create color palette
  genotype_colors <- c("REF" = "gold", "HET" = "springgreen4", "ALT" = "purple4")
  
  # Create list to hold plots
  chromosome_plots <- list()
  
  # Process each chromosome in order
  for (chrom in CHROMOSOME_ORDER) {
    # Prepare data for this chromosome
    stacked_data <- prepare_stacked_data(data, chrom)
    
    if (nrow(stacked_data) > 0) {
      # Create sample order and set factor levels
      sample_order <- create_sample_order(stacked_data, order)
      stacked_data$sample <- factor(stacked_data$sample, levels = sample_order)
    }
    
    # Create plot
    chromosome_plots[[chrom]] <- create_chromosome_plot(stacked_data, chrom, genotype_colors)
  }
  
  return(chromosome_plots)
}

# Convert run-length encoded introgression data to matrix format for ComplexHeatmap
convert_to_matrix <- function(processed_data, bin_size = 1000000, 
                              chromosomes = paste0("chr", 1:10),
                              fill_missing = TRUE) {
  
  # Get chromosome lengths for creating full genome bins
  chrom_lengths <- generate_chromosome_lengths()
  
  # Create complete bin positions across genome
  all_bins <- data.frame()
  for (chrom in chromosomes) {
    chr_length <- chrom_lengths[chrom_lengths$chrom == chrom, "length"]
    max_bin <- ceiling(chr_length / bin_size)
    
    chr_bins <- data.frame(
      chrom = chrom,
      bin_pos = 1:max_bin,
      genome_pos = paste0(chrom, "_", sprintf("%06d", 1:max_bin))
    )
    all_bins <- rbind(all_bins, chr_bins)
  }
  
  # Convert segments to bins
  segment_bins <- data.frame()
  
  for (i in 1:nrow(processed_data)) {
    segment <- processed_data[i, ]
    
    # Calculate which bins this segment covers
    start_bin <- ceiling(segment$start / bin_size)
    end_bin <- ceiling(segment$end / bin_size)
    
    # Create records for each bin covered by this segment
    if (start_bin <= end_bin) {
      bins_covered <- data.frame(
        sample = segment$sample,
        chrom = segment$chrom,
        bin_pos = start_bin:end_bin,
        genotype = segment$genotype
      )
      segment_bins <- rbind(segment_bins, bins_covered)
    }
  }
  
  # Add genome position identifier
  segment_bins$genome_pos <- paste0(segment_bins$chrom, "_", 
                                    sprintf("%06d", segment_bins$bin_pos))
  
  # Convert to wide format matrix
  genotype_matrix <- segment_bins %>%
    select(sample, genome_pos, genotype) %>%
    pivot_wider(names_from = sample, values_from = genotype, 
                values_fill = if(fill_missing) "REF" else NA) %>%
    column_to_rownames("genome_pos") %>%
    as.matrix()
  
  # Ensure all bins are represented if fill_missing is TRUE
  if (fill_missing) {
    missing_bins <- setdiff(all_bins$genome_pos, rownames(genotype_matrix))
    if (length(missing_bins) > 0) {
      # Create matrix for missing bins filled with REF
      missing_matrix <- matrix("REF", 
                               nrow = length(missing_bins), 
                               ncol = ncol(genotype_matrix))
      rownames(missing_matrix) <- missing_bins
      colnames(missing_matrix) <- colnames(genotype_matrix)
      
      # Combine matrices
      genotype_matrix <- rbind(genotype_matrix, missing_matrix)
    }
    
    # Sort by genomic position
    genotype_matrix <- genotype_matrix[all_bins$genome_pos[all_bins$genome_pos %in% rownames(genotype_matrix)], ]
  }
  
  # Create metadata for annotations
  bin_metadata <- all_bins %>%
    filter(genome_pos %in% rownames(genotype_matrix)) %>%
    mutate(
      chromosome = factor(chrom, levels = CHROMOSOME_ORDER),
      position_mb = bin_pos * bin_size / 1e6
    )
  
  return(list(
    matrix = genotype_matrix,
    bin_metadata = bin_metadata,
    chromosome_lengths = chrom_lengths
  ))
}