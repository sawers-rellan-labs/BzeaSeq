#' Global Variable Declarations
#'
#' This file declares global variables used in NSE (non-standard evaluation)
#' contexts to avoid R CMD CHECK warnings.
#'
#' @name globals
NULL

# Declare global variables to avoid R CMD CHECK warnings
globalVariables(c(
  # Column names used in dplyr operations
  "sample", "chrom", "start", "end", "genotype", "bin_start", "bin_end",
  "bin_pos", "genome_pos", "span", "bin", "state",
  
  # Introgression analysis variables
  "is_introgression", "introgression_group", "group_is_introgression",
  "block_start", "block_end", "block_span", "block_midpoint",
  "largest_introgression_midpoint",
  
  # Statistical variables
  "length", "segment_id"
))