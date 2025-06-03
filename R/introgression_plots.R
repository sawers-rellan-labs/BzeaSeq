#' Introgression Visualization and Analysis Functions
#'
#' This module provides comprehensive functions for processing, analyzing, and
#' visualizing introgression data from maize ancestry analysis pipelines.

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# Define chromosome order as a constant
CHROMOSOME_ORDER <- paste0("chr", 1:10)

#' Generate Chromosome Length Data
#'
#' Returns exact chromosome lengths for the B73 reference genome assembly.
#'
#' @param assembly Character. Genome assembly version (currently only "B73v5" supported)
#'
#' @return Data frame with columns:
#'   \item{chrom}{Chromosome name (chr1-chr10)}
#'   \item{length}{Chromosome length in base pairs}
#'
#' @details
#' Chromosome lengths are from Zm-B73-REFERENCE-NAM-5.0.fa.fai.
#' Values represent the exact nuclear chromosome lengths used for
#' bin calculations and genome-wide visualizations.
#'
#' @examples
#' chrom_lengths <- generate_chromosome_lengths()
#' print(chrom_lengths)
#'
#' @export
generate_chromosome_lengths <- function(assembly = "B73v5") {
  if (!identical(assembly, "B73v5")) {
    stop("Currently only B73v5 assembly is supported")
  }
  
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
    ),
    stringsAsFactors = FALSE
  )
  
  chrom_lengths
}

#' Re-encode Introgression Blocks
#'
#' Identifies continuous introgression blocks by combining HET and ALT segments
#' and calculates block-level statistics for downstream analysis.
#'
#' @param data Data frame with columns: sample, chrom, bin_start, bin_end, genotype
#'
#' @return Data frame with additional column:
#'   \item{largest_introgression_midpoint}{Position of largest introgression block midpoint}
#'
#' @details
#' Combines adjacent HET and ALT bins into continuous introgression blocks.
#' Identifies the largest introgression block per sample-chromosome for
#' ordering and visualization purposes. Uses run-length encoding for efficiency.
#'
#' @examples
#' # Example introgression data
#' introg_data <- data.frame(
#'   sample = rep("sample1", 6),
#'   chrom = rep("chr1", 6), 
#'   bin_start = seq(1, 6e6, 1e6),
#'   bin_end = seq(1e6, 6e6, 1e6),
#'   genotype = c("REF", "HET", "ALT", "HET", "REF", "REF")
#' )
#'
#' result <- reencode_introgressions(introg_data)
#' print(result)
#'
#' @export
reencode_introgressions <- function(data) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  
  required_cols <- c("sample", "chrom", "bin_start", "bin_end", "genotype")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
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

#' Add Segment Metadata and Snap to Bin Boundaries
#'
#' Adds segment identifiers and snaps coordinates to bin boundaries for
#' consistent visualization and analysis.
#'
#' @param data Data frame with segment data
#' @param bin_size Integer. Size of genomic bins in base pairs (default: 1000000)
#'
#' @return Data frame with additional columns:
#'   \item{segment_id}{Sequential segment identifier within each sample-chromosome}
#'   \item{bin_start}{Segment start snapped to bin boundary}
#'   \item{bin_end}{Segment end snapped to bin boundary (not exceeding chromosome length)}
#'   \item{span}{Updated span based on bin boundaries}
#'   \item{largest_introgression_midpoint}{Position of largest introgression block}
#'
#' @details
#' Ensures all segments align to standard bin boundaries for consistent
#' analysis and visualization. Prevents bin coordinates from exceeding
#' known chromosome lengths.
#'
#' @examples
#' segment_data <- data.frame(
#'   sample = "sample1",
#'   chrom = "chr1", 
#'   start = 1500000,
#'   end = 2300000,
#'   genotype = "ALT"
#' )
#' 
#' result <- add_segment_metadata(segment_data, bin_size = 1000000)
#' print(result)
#'
#' @export
add_segment_metadata <- function(data, bin_size = 1000000) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  if (!is.numeric(bin_size) || bin_size <= 0) {
    stop("bin_size must be a positive number")
  }
  
  required_cols <- c("sample", "chrom", "start", "end", "genotype")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
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
  
  result
}

#' Process Introgression Data for Analysis
#'
#' Main preprocessing function that adds metadata and prepares introgression
#' data for downstream visualization and analysis.
#'
#' @param data Data frame with raw introgression segment data
#' @param bin_size Integer. Genomic bin size in base pairs (default: 1000000)
#'
#' @return Processed data frame ready for visualization functions
#'
#' @details
#' This is the main entry point for processing raw introgression data.
#' Combines metadata addition, bin snapping, and introgression re-encoding
#' in a single convenient function call.
#'
#' @examples
#' raw_data <- data.frame(
#'   sample = c("sample1", "sample1"),
#'   chrom = c("chr1", "chr1"),
#'   start = c(1000000, 3000000),
#'   end = c(2000000, 4000000),
#'   genotype = c("HET", "ALT")
#' )
#' 
#' processed <- process_introgression_data(raw_data)
#' print(processed)
#'
#' @export
process_introgression_data <- function(data, bin_size = 1000000) {
  result <- add_segment_metadata(data, bin_size)
  result
}

#' Prepare Stacked Data for Single Chromosome Visualization
#'
#' Filters and processes introgression data for visualization of a single
#' chromosome across multiple samples.
#'
#' @param data Processed introgression data from `process_introgression_data()`
#' @param chrom_to_plot Character. Chromosome name to visualize (e.g., "chr1")
#'
#' @return Data frame formatted for stacked visualization with columns:
#'   \item{span}{Segment span in base pairs}
#'   \item{bin}{Sequential bin number for visualization}
#'   Plus all original columns from input data
#'
#' @details
#' Prepares data for stacked bar plots by calculating spans and adding
#' sequential bin numbers for proper visualization ordering.
#'
#' @examples
#' # Assume processed_data exists from process_introgression_data()
#' chr1_data <- prepare_stacked_data(processed_data, "chr1")
#' print(head(chr1_data))
#'
#' @export
prepare_stacked_data <- function(data, chrom_to_plot) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  if (!is.character(chrom_to_plot) || length(chrom_to_plot) != 1) {
    stop("chrom_to_plot must be a single character string")
  }
  if (!chrom_to_plot %in% CHROMOSOME_ORDER) {
    warning(paste("chrom_to_plot", chrom_to_plot, "not in standard chromosome order"))
  }
  
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
  
  result_df
}

#' Create Sample Order for Visualization
#'
#' Generates sample ordering based on different criteria for consistent
#' visualization across plots.
#'
#' @param stacked_data Data frame prepared by `prepare_stacked_data()`
#' @param order Character or character vector specifying ordering method:
#'   \itemize{
#'     \item{"sample"}{Alphabetical order by sample name}
#'     \item{"position"}{Order by largest introgression midpoint position}
#'     \item{Custom vector}{User-specified order of sample names}
#'   }
#'
#' @return Character vector of sample names in desired order
#'
#' @details
#' Supports multiple ordering strategies:
#' - "sample": Simple alphabetical ordering
#' - "position": Orders samples by the position of their largest introgression
#' - Custom vector: User provides exact sample order
#'
#' @examples
#' # Alphabetical ordering
#' sample_order <- create_sample_order(stacked_data, "sample")
#' 
#' # Order by introgression position
#' position_order <- create_sample_order(stacked_data, "position")
#' 
#' # Custom ordering
#' custom_order <- create_sample_order(stacked_data, c("sample3", "sample1", "sample2"))
#'
#' @export
create_sample_order <- function(stacked_data, order) {
  # Input validation
  if (!is.data.frame(stacked_data)) {
    stop("stacked_data must be a data frame")
  }
  if (!"sample" %in% colnames(stacked_data)) {
    stop("stacked_data must contain a 'sample' column")
  }
  
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
      if (!"largest_introgression_midpoint" %in% colnames(stacked_data)) {
        stop("stacked_data must contain 'largest_introgression_midpoint' column for position ordering")
      }
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

#' Create Single Chromosome Plot
#'
#' Internal function to create a stacked bar plot for a single chromosome.
#'
#' @param stacked_data Data frame prepared for single chromosome
#' @param chrom Character. Chromosome name for plot title
#' @param genotype_colors Named vector of colors for genotypes
#'
#' @return ggplot object with stacked bar visualization
#'
#' @keywords internal
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
  
  p
}

#' Plot Introgression Segments Using Stacked Bars
#'
#' Main visualization function that creates stacked bar plots showing
#' introgression patterns across all chromosomes.
#'
#' @param data Processed introgression data from `process_introgression_data()`
#' @param chrom_lengths Data frame with chromosome lengths (optional, will be generated if NULL)
#' @param order Character or character vector specifying sample ordering method
#'   (default: "sample")
#'
#' @return Named list of ggplot objects, one for each chromosome
#'
#' @details
#' Creates comprehensive visualization showing introgression patterns across
#' the genome. Each chromosome gets its own stacked bar plot with samples
#' ordered consistently according to the specified method.
#'
#' @examples
#' # Create basic stacked plots
#' plots <- plot_introgression_stacked(processed_data)
#' 
#' # Order by introgression position
#' plots_ordered <- plot_introgression_stacked(processed_data, order = "position")
#' 
#' # Display chromosome 1 plot
#' print(plots[["chr1"]])
#'
#' @export
plot_introgression_stacked <- function(data, chrom_lengths = NULL, order = "sample") {
  # Input validation
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  
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
  
  chromosome_plots
}

#' Convert Run-Length Encoded Data to Matrix Format
#'
#' Converts segment-based introgression data to a matrix format suitable
#' for ComplexHeatmap visualization and other matrix-based analyses.
#'
#' @param processed_data Data frame with processed introgression segments
#' @param bin_size Integer. Genomic bin size in base pairs (default: 1000000)
#' @param chromosomes Character vector. Chromosomes to include (default: chr1-chr10)
#' @param fill_missing Logical. Whether to fill missing bins with "REF" (default: TRUE)
#'
#' @return List containing:
#'   \item{matrix}{Matrix with samples as columns, genomic bins as rows}
#'   \item{bin_metadata}{Data frame with bin position information}
#'   \item{chromosome_lengths}{Data frame with chromosome length information}
#'
#' @details
#' Converts variable-length segments to fixed-size genomic bins for
#' matrix-based visualization. Each bin is assigned the genotype of
#' overlapping segments. Missing bins are filled with "REF" by default.
#'
#' @examples
#' # Convert to matrix format
#' matrix_data <- convert_to_matrix(processed_data, bin_size = 1000000)
#' 
#' # Access the genotype matrix
#' genotype_matrix <- matrix_data$matrix
#' 
#' # Get bin metadata for annotations
#' bin_info <- matrix_data$bin_metadata
#'
#' @export
convert_to_matrix <- function(processed_data, bin_size = 1000000, 
                              chromosomes = paste0("chr", 1:10),
                              fill_missing = TRUE) {
  # Input validation
  if (!is.data.frame(processed_data)) {
    stop("processed_data must be a data frame")
  }
  if (!is.numeric(bin_size) || bin_size <= 0) {
    stop("bin_size must be a positive number")
  }
  if (!is.character(chromosomes) || length(chromosomes) == 0) {
    stop("chromosomes must be a non-empty character vector")
  }
  if (!is.logical(fill_missing)) {
    stop("fill_missing must be logical")
  }
  
  required_cols <- c("sample", "chrom", "start", "end", "genotype")
  missing_cols <- setdiff(required_cols, colnames(processed_data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Get chromosome lengths for creating full genome bins
  chrom_lengths <- generate_chromosome_lengths()
  
  # Create complete bin positions across genome
  all_bins <- data.frame()
  for (chrom in chromosomes) {
    chr_length <- chrom_lengths[chrom_lengths$chrom == chrom, "length"]
    if (length(chr_length) == 0) {
      warning(paste("Unknown chromosome:", chrom))
      next
    }
    max_bin <- ceiling(chr_length / bin_size)
    
    chr_bins <- data.frame(
      chrom = chrom,
      bin_pos = 1:max_bin,
      genome_pos = paste0(chrom, "_", sprintf("%06d", 1:max_bin)),
      stringsAsFactors = FALSE
    )
    all_bins <- rbind(all_bins, chr_bins)
  }
  
  # Convert segments to bins
  segment_bins <- data.frame()
  
  for (i in seq_len(nrow(processed_data))) {
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
        genotype = segment$genotype,
        stringsAsFactors = FALSE
      )
      segment_bins <- rbind(segment_bins, bins_covered)
    }
  }
  
  # Add genome position identifier
  segment_bins$genome_pos <- paste0(segment_bins$chrom, "_", 
                                    sprintf("%06d", segment_bins$bin_pos))
  
  # Convert to wide format matrix
  if (nrow(segment_bins) > 0) {
    genotype_matrix <- segment_bins %>%
      select(sample, genome_pos, genotype) %>%
      pivot_wider(names_from = sample, values_from = genotype, 
                  values_fill = if(fill_missing) "REF" else NA) %>%
      column_to_rownames("genome_pos") %>%
      as.matrix()
  } else {
    # Create empty matrix if no segments
    genotype_matrix <- matrix(character(0), nrow = 0, ncol = 0)
  }
  
  # Ensure all bins are represented if fill_missing is TRUE
  if (fill_missing && nrow(genotype_matrix) > 0) {
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
    available_positions <- all_bins$genome_pos[all_bins$genome_pos %in% rownames(genotype_matrix)]
    genotype_matrix <- genotype_matrix[available_positions, , drop = FALSE]
  }
  
  # Create metadata for annotations
  bin_metadata <- all_bins %>%
    filter(genome_pos %in% rownames(genotype_matrix)) %>%
    mutate(
      chromosome = factor(chrom, levels = CHROMOSOME_ORDER),
      position_mb = bin_pos * bin_size / 1e6
    )
  
  list(
    matrix = genotype_matrix,
    bin_metadata = bin_metadata,
    chromosome_lengths = chrom_lengths
  )
}
