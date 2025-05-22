
# Analyze state transitions
cat("\n==== State Transition Analysis ====\n")

# Pre-allocate a data frame for all transitions
transition_data <- data.frame()

# Process each sample separately
for (samp in unique(read_freq$SAMPLE)) {
  # Process each chromosome separately
  for (chr in chrom_order) {
    # Get states for this chromosome and sample
    chr_states <- read_freq %>% 
      filter(SAMPLE == samp, CONTIG == chr) %>% 
      arrange(BIN_POS) %>% 
      pull(GENOTYPE)
    
    # Skip if empty
    if (length(chr_states) <= 1) next
    
    # Find transitions - positions where state changes
    for (i in 1:(length(chr_states)-1)) {
      if (chr_states[i] != chr_states[i+1]) {
        transition_data <- rbind(transition_data, data.frame(
          sample = samp,
          chromosome = chr,
          position = i,  # Position in bin sequence
          from_state = chr_states[i],
          to_state = chr_states[i+1]
        ))
      }
    }
  }
}

# Analyze transitions if any found
if (nrow(transition_data) > 0) {
  # Count transition types
  transition_counts <- transition_data %>%
    group_by(from_state, to_state) %>%
    summarize(
      count = n(),
      percentage = n() / nrow(transition_data) * 100,
      .groups = "drop"
    ) %>%
    mutate(transition = paste(from_state, "→", to_state)) %>%
    arrange(desc(count))
  
  cat("State transition counts:\n")
  print(transition_counts, row.names = FALSE)
  
  # Calculate transition rate per 100Mb
  # Assuming BIN_SIZE is in bp, convert to Mb and calculate per 100Mb
  total_mb <- nrow(read_freq) * bin_size / 1e6
  
  cat("\nTransition rates per 100Mb:\n")
  transition_rates <- transition_counts %>%
    mutate(
      rate_per_100mb = count / total_mb * 100
    ) %>%
    select(transition, count, rate_per_100mb)
  
  print(transition_rates, row.names = FALSE)
  
  # Calculate segment length statistics 
  # A segment is defined as consecutive bins with the same state
  
  segment_summary <- data.frame()
  
  for (samp in unique(read_freq$SAMPLE)) {
    # Process each chromosome separately
    for (chr in chrom_order) {
      # Get states for this chromosome and sample
      chr_data <- read_freq %>% 
        filter(SAMPLE == samp, CONTIG == chr) %>% 
        arrange(BIN_POS)
      
      # Skip if empty
      if (nrow(chr_data) == 0) next
      
      # Get the state sequence
      chr_states <- chr_data$GENOTYPE
      
      # Use run-length encoding to find segments
      rle_result <- rle(as.character(chr_states))
      
      # Get the values (states) and lengths
      states <- rle_result$values
      lengths <- rle_result$lengths
      
      # Record all segments
      for (i in 1:length(states)) {
        segment_summary <- rbind(segment_summary, data.frame(
          sample = samp,
          chromosome = chr,
          state = states[i],
          length_bins = lengths[i],
          length_mb = lengths[i] * (bin_size / 1e6) # Convert bin count to Mb
        ))
      }
    }
  }
  
  # Calculate summary statistics for each state
  cat("\nSegment length statistics by state (in Mb):\n")
  segment_summary %>%
    group_by(state) %>%
    summarize(
      count = n(),
      min_length = min(length_mb),
      q1_length = quantile(length_mb, 0.25),
      median_length = median(length_mb),
      mean_length = mean(length_mb),
      q3_length = quantile(length_mb, 0.75),
      max_length = max(length_mb),
      total_mb = sum(length_mb)
    ) %>%
    arrange(state) %>%
    print()
  
  # Plot transition count barplot
  transition_plot_file <- file.path(output_dir,paste0(output_prefix, "_transitions.pdf"))
  pdf(transition_plot_file, width = 8, height = 6)
  # Create more readable labels
  transition_counts$transition <- factor(transition_counts$transition, 
                                         levels = transition_counts$transition[order(transition_counts$count, decreasing = TRUE)])
  # Color coding
  transition_colors <- c("REF → HET" = "orange", "HET → REF" = "orange", 
                         "HET → ALT" = "purple", "ALT → HET" = "purple",
                         "REF → ALT" = "red", "ALT → REF" = "red")
  barplot(transition_counts$count, names.arg = transition_counts$transition,
          main = "State Transition Counts",
          xlab = "Transition Type",
          ylab = "Count",
          col = transition_colors[as.character(transition_counts$transition)],
          las = 2)  # Rotate labels
  dev.off()
  cat(paste("Transition count plot saved to:", transition_plot_file, "\n"))
} else {
  cat("No state transitions found in the data.\n")
}
} else {
  cat("Warning: No non-REF segments found. BED file not created.\n")
}

# Create a chromosome-level summary by sample
chrom_summary <- read_freq %>%
  # Ensure chromosome ordering
  mutate(CONTIG = factor(CONTIG, levels = chrom_order)) %>%
  group_by(SAMPLE, CONTIG) %>%
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
  arrange(SAMPLE, CONTIG)

summary_file <- file.path(output_dir,paste0(output_prefix, "_chromosome_summary.tsv"))
write.table(chrom_summary, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste("Chromosome-level summary saved to:", summary_file, "\n"))

# Calculate overall statistics per sample
overall_summary <- read_freq %>%
  group_by(SAMPLE) %>%
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
  # Sort by non-REF percentage descending
  arrange(desc(non_ref_pct))

overall_file <- file.path(output_dir,paste0(output_prefix, "_summary.tsv"))
write.table(overall_summary, overall_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste("Overall summary saved to:", overall_file, "\n"))

# Generate a summary plot
summary_plot_file <- file.path(output_dir,paste0(output_prefix, "_summary_plot.pdf"))
pdf(summary_plot_file, width = 12, height = 8)
# Create a stacked barplot of REF/HET/ALT percentages by sample
overall_summary_plot <- overall_summary %>%
  select(SAMPLE, REF_pct, HET_pct, ALT_pct) %>%
  # Reshape for ggplot
  pivot_longer(cols = c(REF_pct, HET_pct, ALT_pct), names_to = "Genotype", values_to = "Percentage") %>%
  # Clean up genotype names
  mutate(Genotype = gsub("_pct", "", Genotype)) %>%
  # Order samples by non-REF percentage
  mutate(SAMPLE = factor(SAMPLE, levels = overall_summary$SAMPLE))

# Create the plot
ggplot(overall_summary_plot, aes(x = SAMPLE, y = Percentage, fill = Genotype)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("REF" = "gold", "HET" = "springgreen4", "ALT" = "purple4")) +
  labs(
    title = "Ancestry Proportions by Sample",
    x = "Sample",
    y = "Percentage (%)",
    fill = "Genotype"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )
dev.off()
cat(paste("Summary plot saved to:", summary_plot_file, "\n"))
