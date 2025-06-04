# ============================================================================
# Helper functions to add to your R package
# File: R/ancestry_processing.R
# ============================================================================

#' Process Allelic Counts to Genomic Bins
#'
#' @param read_count Data frame with allelic count data from GATK
#' @param sample_name Character. Sample identifier
#' @param bin_size Integer. Bin size in base pairs (default: 1000000)
#' @param chrom_order Character vector. Chromosome order
#'
#' @return Data frame with bin-level statistics
#' @export
process_allelic_counts_to_bins <- function(read_count, sample_name, 
                                           bin_size = 1000000, 
                                           chrom_order = paste0("chr", 1:10)) {
  # Input validation
  if (!is.data.frame(read_count)) {
    stop("read_count must be a data frame")
  }
  
  # Ensure chromosome factor ordering
  read_count$CONTIG <- factor(read_count$CONTIG, levels = chrom_order)
  
  # Calculate bin-level statistics
  read_freq <- read_count %>%
    group_by(CONTIG) %>%
    mutate(
      BIN_POS = ceiling(POSITION / bin_size),
      READ_DEPTH = ALT_COUNT + REF_COUNT
    ) %>%
    group_by(CONTIG, BIN_POS) %>%
    summarise(
      SAMPLE = sample_name,
      VARIANT_COUNT = n(),
      INFORMATIVE_VARIANT_COUNT = sum(ALT_COUNT + REF_COUNT > 0),
      DEPTH_SUM = sum(READ_DEPTH),
      ALT_COUNT = sum(ALT_COUNT),
      ALT_FREQ = ifelse(sum(READ_DEPTH) > 0, sum(ALT_COUNT) / sum(READ_DEPTH), 0),
      BIN_START = min(POSITION),
      BIN_END = max(POSITION),
      .groups = 'drop'
    ) %>%
    select(SAMPLE, everything())
  
  read_freq$CONTIG <- factor(read_freq$CONTIG, levels = chrom_order)
  read_freq
}

#' Apply Ancestry Clustering to Bin Data
#'
#' @param read_freq Data frame with bin-level allelic frequency data
#'
#' @return Data frame with clustering results added
#' @export
apply_ancestry_clustering <- function(read_freq) {
  # This function applies the clustering logic from the original script
  # but in a modular way
  
  if (!requireNamespace("Ckmeans.1d.dp", quietly = TRUE)) {
    stop("Package 'Ckmeans.1d.dp' is required but not available")
  }
  if (!requireNamespace("rebmix", quietly = TRUE)) {
    stop("Package 'rebmix' is required but not available")
  }
  
  # Separate non-zero bins for clustering
  non_zero_bins <- read_freq %>%
    filter(ALT_FREQ > 0) %>%
    mutate(asin = asin(sqrt(ALT_FREQ)))
  
  # Initialize clustering columns
  read_freq$Kasin <- factor(NA, levels = c("REF", "HET", "ALT"))
  read_freq$Kgmm <- factor(NA, levels = c("REF", "HET", "ALT"))
  
  if (nrow(non_zero_bins) > 0) {
    # Apply clustering (simplified version)
    K <- 3
    
    # K-means on arcsin-transformed data
    non_zero_bins$Kasin <- as.factor(
      Ckmeans.1d.dp::Ckmeans.1d.dp(non_zero_bins$asin, K)$cluster
    )
    
    # Gaussian mixture model (simplified)
    tryCatch({
      normalest <- rebmix::REBMIX(
        Dataset = list(data.frame(Value = non_zero_bins$asin)),
        Preprocessing = "histogram",
        cmin = 3, cmax = 3,
        Criterion = "BIC",
        pdf = "normal"
      )
      normclu <- rebmix::RCLRMIX(x = normalest)
      non_zero_bins$Kgmm <- as.factor(normclu@Zp)
    }, error = function(e) {
      non_zero_bins$Kgmm <- non_zero_bins$Kasin
    })
    
    # Relabel clusters based on ALT_FREQ
    non_zero_bins$Kasin <- relabel_clusters(non_zero_bins$Kasin, non_zero_bins)
    non_zero_bins$Kgmm <- relabel_clusters(non_zero_bins$Kgmm, non_zero_bins)
    
    # Merge back to main dataset
    read_freq <- read_freq %>%
      left_join(
        non_zero_bins %>% select(CONTIG, BIN_POS, Kasin, Kgmm),
        by = c("CONTIG", "BIN_POS")
      ) %>%
      mutate(
        Kasin = ifelse(ALT_FREQ == 0, "REF", Kasin),
        Kgmm = ifelse(ALT_FREQ == 0, "REF", Kgmm)
      )
  }
  
  # Convert to factors
  read_freq$Kasin <- factor(read_freq$Kasin, levels = c("REF", "HET", "ALT"))
  read_freq$Kgmm <- factor(read_freq$Kgmm, levels = c("REF", "HET", "ALT"))
  
  read_freq
}

#' Relabel Clusters Based on ALT_FREQ
#'
#' @param clusters Factor. Cluster assignments
#' @param data Data frame containing ALT_FREQ column
#'
#' @return Factor with relabeled clusters (REF, HET, ALT)
relabel_clusters <- function(clusters, data) {
  cluster_means <- tapply(data$ALT_FREQ, clusters, mean)
  ordered_clusters <- order(cluster_means)
  
  cluster_map <- rep(NA, length(unique(as.numeric(clusters))))
  cluster_map[ordered_clusters[1]] <- "REF"
  cluster_map[ordered_clusters[2]] <- "HET"
  cluster_map[ordered_clusters[3]] <- "ALT"
  
  factor(cluster_map[as.numeric(clusters)], levels = c("REF", "HET", "ALT"))
}

#' Filter Valid Samples from Metadata
#'
#' @param metadata Data frame with sample metadata
#'
#' @return Character vector of valid sample names
#' @export
filter_valid_samples <- function(metadata) {
  # This logic should match your metadata filtering criteria
  in_bzea <- !is.na(metadata$project) & metadata$project == "bzea"
  is_valid_sample <- in_bzea | metadata$is_check == TRUE
  metadata$sample[is_valid_sample]
}

#' Save Ancestry Analysis Results
#'
#' @param bin_data Data frame with bin-level results
#' @param segments_data Data frame with segment-level results  
#' @param output_prefix Character. Output file prefix
#' @param sample_name Character. Sample identifier
#' @param bin_size Integer. Bin size used
#'
#' @export
save_ancestry_results <- function(bin_data, segments_data, output_prefix,
                                  sample_name, bin_size) {
  # Save detailed bin results
  bin_file <- paste0(output_prefix, "_bin_genotypes.tsv")
  write.table(bin_data, bin_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Save segment results if available
  if (!is.null(segments_data) && nrow(segments_data) > 0) {
    segments_file <- paste0(output_prefix, "_non_REF.bed")
    write.table(segments_data, segments_file, sep = "\t", 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  # Save summary statistics
  summary_stats <- create_sample_summary(bin_data, sample_name)
  summary_file <- paste0(output_prefix, "_summary.tsv")
  write.table(summary_stats, summary_file, sep = "\t", 
              row.names = FALSE, quote = FALSE)
}

#' Create Sample Summary Statistics
#'
#' @param bin_data Data frame with bin-level genotype data
#' @param sample_name Character. Sample identifier
#'
#' @return Data frame with summary statistics
create_sample_summary <- function(bin_data, sample_name) {
  data.frame(
    sample = sample_name,
    total_bins = nrow(bin_data),
    REF_bins = sum(bin_data$GENOTYPE == "REF", na.rm = TRUE),
    HET_bins = sum(bin_data$GENOTYPE == "HET", na.rm = TRUE),
    ALT_bins = sum(bin_data$GENOTYPE == "ALT", na.rm = TRUE),
    REF_pct = 100 * sum(bin_data$GENOTYPE == "REF", na.rm = TRUE) / nrow(bin_data),
    HET_pct = 100 * sum(bin_data$GENOTYPE == "HET", na.rm = TRUE) / nrow(bin_data),
    ALT_pct = 100 * sum(bin_data$GENOTYPE == "ALT", na.rm = TRUE) / nrow(bin_data),
    non_ref_pct = 100 * sum(bin_data$GENOTYPE %in% c("HET", "ALT"), na.rm = TRUE) / nrow(bin_data)
  )
}