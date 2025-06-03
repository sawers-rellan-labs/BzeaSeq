#' Breeding Genetics and HMM Functions for Ancestry Analysis
#'
#' This module provides functions for calculating expected genotype frequencies
#' in breeding populations and applying Hidden Markov Model smoothing to
#' ancestry calls in maize introgression lines.

library(HMM)

#' Create Mating Matrices for Breeding Schemes
#'
#' Creates transition matrices for backcrossing and selfing operations
#' used in Near Isogenic Line (NIL) development.
#'
#' @return A list containing three 3x3 transition matrices:
#'   \item{backcross_AA}{Matrix for backcrossing to AA donor}
#'   \item{backcross_aa}{Matrix for backcrossing to aa donor}  
#'   \item{selfing}{Matrix for self-fertilization}
#'
#' @details
#' Each matrix represents genotype transitions where:
#' - Rows represent offspring genotypes (AA, Aa, aa)
#' - Columns represent parent genotypes (AA, Aa, aa)
#' - Values are transition probabilities
#'
#' @examples
#' matrices <- create_mating_matrices()
#' matrices$selfing  # Shows selfing transition probabilities
#'
#' @export
create_mating_matrices <- function() {
  if (!requireNamespace("HMM", quietly = TRUE)) {
    stop("Package 'HMM' is required but not available")
  }
  
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
  
  list(
    backcross_AA = backcross_AA,
    backcross_aa = backcross_aa,
    selfing = selfing
  )
}

#' Calculate Genotype Frequencies After Breeding Scheme
#'
#' Calculates expected genotype frequencies following a specific breeding
#' scheme involving backcrossing and selfing generations.
#'
#' @param bc Integer. Number of backcross generations (default: 2)
#' @param s Integer. Number of selfing generations (default: 3)
#' @param donor_type Character. Type of donor parent, either "AA" or "aa" (default: "aa")
#' @param pop Numeric vector of length 3. Initial population frequencies for
#'   AA, Aa, aa genotypes (default: c(0,1,0) representing F1)
#'
#' @return Named numeric vector with final genotype frequencies for AA, Aa, aa
#'
#' @details
#' This function simulates breeding schemes commonly used to develop Near
#' Isogenic Lines (NILs). The default BC2S3 scheme represents:
#' - 2 backcross generations to the recurrent parent
#' - 3 generations of self-fertilization
#'
#' @examples
#' # Calculate BC2S3 frequencies
#' bc2s3_freq <- calculate_nil_frequencies(bc = 2, s = 3, donor_type = "aa")
#' print(bc2s3_freq)
#'
#' # Calculate BC1S2 frequencies  
#' bc1s2_freq <- calculate_nil_frequencies(bc = 1, s = 2, donor_type = "aa")
#' print(bc1s2_freq)
#'
#' @export
calculate_nil_frequencies <- function(bc = 2, s = 3, donor_type = "aa", 
                                      pop = c(0, 1, 0)) {
  # Input validation
  if (!is.numeric(bc) || !is.numeric(s) || bc < 0 || s < 0) {
    stop("Number of backcrosses and self generations must be non-negative integers")
  }
  if (!donor_type %in% c("AA", "aa")) {
    stop("Donor type must be 'AA' or 'aa'")
  }
  if (!is.numeric(pop) || length(pop) != 3 || abs(sum(pop) - 1) > 1e-6) {
    stop("Population frequencies must be numeric vector of length 3 summing to 1")
  }
  
  # Get mating matrices
  matrices <- create_mating_matrices()
  
  # Initial population vector: [AA frequency, Aa frequency, aa frequency]
  current_population <- matrix(pop, ncol = 1)
  
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
  
  final_frequencies
}

#' Convert NIL Frequencies to HMM Format
#'
#' Converts genotype frequencies from breeding scheme calculations to the
#' REF/HET/ALT format used in HMM ancestry analysis.
#'
#' @param bc Integer. Number of backcross generations (default: 2)
#' @param s Integer. Number of selfing generations (default: 3)
#' @param donor_type Character. Type of donor parent, either "AA" or "aa" (default: "aa")
#' @param pop Numeric vector of length 3. Initial population frequencies (default: c(0,1,0))
#' @param ... Additional arguments passed to `calculate_nil_frequencies()`
#'
#' @return Named numeric vector with frequencies for REF, HET, ALT states
#'
#' @details
#' This function maps biological genotypes to ancestry states:
#' - REF: Recurrent parent background
#' - HET: Heterozygous regions (admixed)
#' - ALT: Donor introgression regions
#'
#' The mapping depends on donor_type:
#' - If donor_type = "aa": REF=AA, HET=Aa, ALT=aa
#' - If donor_type = "AA": REF=aa, HET=Aa, ALT=AA
#'
#' @examples
#' # BC2S3 with teosinte (aa) donor
#' bc2s3_hmm <- nil_frequencies_for_hmm(bc = 2, s = 3, donor_type = "aa")
#' print(bc2s3_hmm)
#'
#' @export
nil_frequencies_for_hmm <- function(bc = 2, s = 3, donor_type = "aa", 
                                    pop = c(0, 1, 0), ...) {
  # Calculate raw genotype frequencies
  genotype_freqs <- calculate_nil_frequencies(bc, s, donor_type, pop, ...)
  
  # Convert to HMM format based on which allele is being introgressed
  if (donor_type == "AA") {
    hmm_frequencies <- c(
      ALT = genotype_freqs["AA"],  # Donor introgression
      HET = genotype_freqs["Aa"],  # Heterozygous regions  
      REF = genotype_freqs["aa"]   # Recurrent parent background
    )
  } else {
    hmm_frequencies <- c(
      REF = genotype_freqs["AA"],  # Recurrent parent background
      HET = genotype_freqs["Aa"],  # Heterozygous regions
      ALT = genotype_freqs["aa"]   # Donor introgression  
    )
  }
  
  # Ensure proper ordering
  hmm_frequencies[c("REF", "HET", "ALT")]
}

#' Apply HMM Smoothing to Ancestry Calls
#'
#' Applies Hidden Markov Model smoothing to noisy ancestry classifications
#' to enforce biological constraints and reduce classification errors.
#'
#' @param genotypes Character vector or factor. Ancestry calls with levels
#'   c("REF", "HET", "ALT")
#' @param transitions Numeric vector of length 2. Transition probabilities
#'   where transitions[1] = P(stay in same state) and transitions[2] = P(switch states)
#'   (default: c(0.995, 0.005))
#'
#' @return Factor with smoothed ancestry calls (levels: REF, HET, ALT)
#'
#' @details
#' The HMM uses:
#' - Start probabilities from BC2S3 breeding expectations
#' - High self-transition probabilities to model genomic linkage
#' - Emission probabilities that account for classification uncertainty
#'
#' The default transition probabilities (0.995 stay, 0.005 switch) assume
#' approximately 1 Mb bins and 1 cM/Mb recombination rate.
#'
#' @examples
#' # Simulate noisy ancestry calls
#' noisy_calls <- factor(c("REF", "ALT", "REF", "REF", "HET", "REF"), 
#'                       levels = c("REF", "HET", "ALT"))
#'
#' # Apply smoothing
#' smoothed <- smooth_ancestry_with_hmm(noisy_calls)
#' print(data.frame(original = noisy_calls, smoothed = smoothed))
#'
#' @export
smooth_ancestry_with_hmm <- function(genotypes, transitions = c(0.995, 0.005)) {
  # Input validation
  if (!is.factor(genotypes) && !is.character(genotypes)) {
    stop("genotypes must be a character vector or factor")
  }
  if (!is.numeric(transitions) || length(transitions) != 2) {
    stop("transitions must be numeric vector of length 2")
  }
  if (any(transitions < 0) || any(transitions > 1) || sum(transitions) > 1) {
    stop("transition probabilities must be between 0 and 1, and sum <= 1")
  }
  if (!requireNamespace("HMM", quietly = TRUE)) {
    stop("Package 'HMM' is required but not available")
  }
  
  # Convert genotypes to numeric (0=REF, 1=HET, 2=ALT)
  genotype_levels <- c("REF", "HET", "ALT")
  if (is.factor(genotypes)) {
    if (!all(levels(genotypes) %in% genotype_levels)) {
      stop("genotypes factor must have levels from: REF, HET, ALT")
    }
  }
  
  geno_factor <- factor(genotypes, levels = genotype_levels)
  geno_numeric <- as.numeric(geno_factor) - 1
  
  # Handle missing values
  if (any(is.na(geno_numeric))) {
    warning("Missing values in genotypes will be treated as REF")
    geno_numeric[is.na(geno_numeric)] <- 0
  }
  
  # Set up transition probabilities (high probability of staying in same state)
  trans_prob <- matrix(c(
    transitions[1], transitions[2]/2, transitions[2]/2,  # From REF
    transitions[2]/2, transitions[1], transitions[2]/2,  # From HET
    transitions[2]/2, transitions[2]/2, transitions[1]   # From ALT
  ), nrow = 3, byrow = TRUE)
  
  # Set up emission probabilities (how likely each state produces observed genotypes)
  emiss_prob <- matrix(c(
    0.9, 0.08, 0.02,  # REF state
    0.1, 0.8, 0.1,    # HET state
    0.02, 0.08, 0.9   # ALT state
  ), nrow = 3, byrow = TRUE)
  
  # Get BC2S3 priors
  bc2s3 <- nil_frequencies_for_hmm(bc = 2, s = 3)
  
  # Initialize HMM with BC2S3 priors
  hmm <- HMM::initHMM(
    c("REF", "HET", "ALT"), 
    c("0", "1", "2"), 
    startProbs = bc2s3,  # Prior probabilities from breeding design
    transProbs = trans_prob, 
    emissionProbs = emiss_prob
  )
  
  # Run Viterbi algorithm to find most likely state sequence
  viterbi_path <- HMM::viterbi(hmm, as.character(geno_numeric))
  
  # Convert back to genotype calls
  smoothed_genotypes <- factor(viterbi_path, levels = c("REF", "HET", "ALT"))
  
  smoothed_genotypes
}