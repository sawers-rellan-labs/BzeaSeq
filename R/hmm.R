library(HMM)     # For HMM smoothing

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
calculate_nil_frequencies <- function(bc=2, s=3, donor_type="aa", pop=c(0,1,0)) {
  if (bc < 0 || s < 0) {
    stop("Number of backcrosses and self generations must be non-negative")
  }
  if (!donor_type %in% c("AA", "aa")) {
    stop("Donor type must be 'AA' or 'aa'")
  }
  
  # Get mating matrices
  matrices <- create_mating_matrices()
  
  # Initial F1 generation - all heterozygous after crossing pure lines
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
  
  return(final_frequencies)
}

# Function to return frequencies in HMM format (REF, HET, ALT)
nil_frequencies_for_hmm <- function(bc, s, donor_type="aa", pop=c(0,1,0)) {
  # Calculate raw genotype frequencies
  genotype_freqs <- calculate_nil_frequencies(bc, s, donor_type,pop)
  
  # Convert to HMM format based on which allele is being introgressed
  if (donor_type == "AA") {
    hmm_frequencies <- c(
      REF = genotype_freqs["aa"],  # Recurrent parent background
      HET = genotype_freqs["Aa"],  # Heterozygous regions  
      ALT = genotype_freqs["AA"]   # Donor introgression
    )
    names(hmm_frequencies) <- c("ALT", "HET", "REF")
  } else {
    hmm_frequencies <- c(
      REF = genotype_freqs["AA"],  # Recurrent parent background
      HET = genotype_freqs["Aa"],  # Heterozygous regions
      ALT = genotype_freqs["aa"]   # Donor introgression  
    )
    names(hmm_frequencies) <- c("REF", "HET", "ALT")
  }

  return(hmm_frequencies)
}


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
