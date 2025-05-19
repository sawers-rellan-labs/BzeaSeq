library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(cowplot)

# Read in data from file
# Replace this with your file path
data <- read.table("~/Desktop/all_samples_bzea_ALT.bed", sep = "\t", na.strings = na.strings = c("NA",""))
colnames(data) <- c("sample", "chrom", "start", "end", "genotype", "numeric", "freq")

bzea <- read.csv("~/Desktop/J2Teo_Final_DB.csv", na.strings = c("NA",""))
colnames(bzea)
head(bzea)

sample_sheet <- read.csv("~/Desktop/BZea-Sample-List.csv", na.strings = c("NA",""))
colnames(sample_sheet)
head(sample_sheet)

sum(sample_sheet$Seq_Full_ID %in% bzea$sequencing.id)
sum(bzea$line_id=="B73-Bulk")
bzea$line_id[grepl("B73",bzea$line_id)]
B73-bulk

metadata <- sample_sheet %>% select(
  field_row = Sample_Origin,
  bam_preffix = Seq_Full_ID,
  pedigree_sample_sheet = Line.ID) %>%
    left_join(bzea %>% select(field_row =  seed_origin,
    bam_preffix = sequencing.id,
                             founder_group =taxa_code,
                             donor_id = accession_id,
                             pedigree_J2TEO = line_id
                             )
             ) %>% filter(!is.na(founder_group))


metadata
# Clean up sample names if needed
# data$sample <- gsub("./ancestry/", "", data$sample, fixed = TRUE)

# Function to calculate segment span and assign IDs to segments
process_introgression_data <- function(data) {
  # Calculate segment span
  data$span <- data$end - data$start
  
  # Sort by chromosome, then descending span size, then start position
  # This ensures the largest segment in each chromosome gets ID=1
  data <- data %>%     
    arrange(sample,chrom, desc(span), start) %>%
    group_by(sample,chrom) %>%
    mutate(segment_id = row_number()) %>%     
  # Get first (largest) segment info for each sample and chromosome
    mutate(
      largest_segment_start = first(start),
      largest_span = first(span)
    ) %>%
     ungroup()
  return(data)
}


# Generate chromosome length data from the assembly index file
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

# Function to prepare chromosome data in a format suitable for stacked bar plots
prepare_stacked_data <- function(data, chrom_to_plot) {
  # Filter data for the specified chromosome
  chrom_data <- data %>% filter(chrom == chrom_to_plot)
  
  # Get chromosome length
  chr_length <- chrom_lengths %>% 
    filter(chrom == chrom_to_plot) %>% 
    pull(length)
  
  # Get all unique samples from the entire dataset
  all_samples <- unique(data$sample)
  
  # Get samples that have segments in this chromosome
  samples_with_segments <- unique(chrom_data$sample)
  
  # Identify missing samples for this chromosome
  missing_samples <- setdiff(all_samples, samples_with_segments)
  
  # Create a result dataframe to hold all segment data for stacked bars
  result_df <- data.frame()
  
  # Process each sample with segments
  for(current_sample in samples_with_segments) {
    # Get segments for this sample
    sample_segments <- chrom_data %>% 
      filter(sample == current_sample) %>%
      arrange(start)
    
    # If no segments for this sample (should not happen), skip
    if(nrow(sample_segments) == 0) next
    
    # Extract the largest segment start position for sorting later
    largest_start <- sample_segments$largest_segment_start[1]
    
    # Create a list to build up all segments including gaps
    all_segments <- list()
    segment_index <- 1
    
    # Add a segment for the start of chromosome if needed
    if(sample_segments$start[1] > 1) {
      all_segments[[segment_index]] <- data.frame(
        sample = current_sample,
        chrom = chrom_to_plot,
        start = 1,
        end = sample_segments$start[1] - 1,
        genotype = "REF",
        segment_id = paste0("gap_", segment_index),
        span = sample_segments$start[1] - 1,
        is_gap = TRUE,
        largest_start = largest_start
      )
      segment_index <- segment_index + 1
    }
    
    # Process each introgression segment and gaps between them
    for(i in 1:nrow(sample_segments)) {
      # Add the current segment
      all_segments[[segment_index]] <- data.frame(
        sample = current_sample,
        chrom = chrom_to_plot,
        start = sample_segments$start[i],
        end = sample_segments$end[i],
        genotype = "ALT",  # We're only plotting ALT segments
        segment_id = sample_segments$segment_id[i],
        span = sample_segments$span[i],
        is_gap = FALSE,
        largest_start = largest_start
      )
      segment_index <- segment_index + 1
      
      # Add gap after this segment if it's not the last segment
      if(i < nrow(sample_segments)) {
        gap_start <- sample_segments$end[i] + 1
        gap_end <- sample_segments$start[i+1] - 1
        
        if(gap_end >= gap_start) {  # Ensure valid gap
          all_segments[[segment_index]] <- data.frame(
            sample = current_sample,
            chrom = chrom_to_plot,
            start = gap_start,
            end = gap_end,
            genotype = "REF",
            segment_id = paste0("gap_", segment_index),
            span = gap_end - gap_start + 1,
            is_gap = TRUE,
            largest_start = largest_start
          )
          segment_index <- segment_index + 1
        }
      }
    }
    
    # Add a segment for the end of chromosome if needed
    if(sample_segments$end[nrow(sample_segments)] < chr_length) {
      all_segments[[segment_index]] <- data.frame(
        sample = current_sample,
        chrom = chrom_to_plot,
        start = sample_segments$end[nrow(sample_segments)] + 1,
        end = chr_length,
        genotype = "REF",
        segment_id = paste0("gap_", segment_index),
        span = chr_length - sample_segments$end[nrow(sample_segments)],
        is_gap = TRUE,
        largest_start = largest_start
      )
    }
    
    # Combine all segments for this sample
    sample_df <- do.call(rbind, all_segments)
    
    # Add bin for proper stacking in ggplot
    sample_df$bin <- 1:nrow(sample_df)
    
    # Add to result
    result_df <- rbind(result_df, sample_df)
  }
  
  # Add full REF chromosomes for missing samples
  if(length(missing_samples) > 0) {
    for(missing_sample in missing_samples) {
      # Use a very large value for largest_start to ensure these show up at the end
      largest_start <- chr_length * 2  # Twice the chromosome length ensures it's larger than any real position
      
      # Create a single full REF segment for this sample
      missing_df <- data.frame(
        sample = missing_sample,
        chrom = chrom_to_plot,
        start = 1,
        end = chr_length,
        genotype = "REF",
        segment_id = "full_ref",
        span = chr_length,
        is_gap = TRUE,
        largest_start = largest_start,
        bin = 1
      )
      
      # Add to result
      result_df <- rbind(result_df, missing_df)
    }
  }
  
  # Return prepared data
  return(result_df)
}

chr_theme <- theme(
  legend.position = "bottom",
  strip.background = element_blank(),
  strip.text.x = element_text(hjust = 0),
  strip.text.y = element_text(angle = 0, face = "bold"),
  panel.spacing.x=unit(0, "null"),
  panel.spacing.y=unit(0, "null") ,
  axis.line.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.line.x = element_blank(),
  axis.title.x = element_blank(),
  axis.text.y = element_blank()
) 

# Plot introgression segments using geom_bar for a stacked visualization
plot_introgression_stacked <- function(data, chrom_lengths) {
  # Get all unique chromosomes
  chromosomes <- paste0("chr",1:10)
  
  # Create a list to hold plots for each chromosome
  chromosome_plots <- list()
  
  # Create a color palette for genotypes
  genotype_colors <- c("REF" = "gold", "ALT" = "purple4")
  
  # Process each chromosome
  for(chrom in chromosomes) {
    # Prepare data for this chromosome
    stacked_data <- prepare_stacked_data(data, chrom)
    
    if(nrow(stacked_data) == 0) {
      # Create empty plot if no data
      chromosome_plots[[chrom]] <- ggplot() + 
        ggtitle(paste("No introgression data for", chrom)) + 
        theme_minimal()
      next
    }
    
    # Sort samples by largest segment start position
    # Full REF chromosomes will be at the end due to their large largest_start value
    sample_order <- stacked_data %>%
      select(sample, largest_start) %>%
      distinct() %>%
      arrange(largest_start) %>%
      pull(sample)
    
    # Set sample as a factor with the determined order
    stacked_data$sample <- factor(stacked_data$sample, levels = sample_order)
    
    # Create the plot with horizontal bars (no labels)
    p <- ggplot(stacked_data, aes(x = sample, y = span, fill = genotype, group = bin)) +
      geom_bar(stat = "identity", position = "stack", width = 1) +
      scale_fill_manual(values = genotype_colors) +
      labs(title = paste("Introgression segments in", chrom),
           x = "Sample", 
           y = "Chromosome Length") +
      scale_y_continuous(labels = function(x) paste0(round(x/1e6, 1))) +
      # Flip coordinates to make horizontal
      coord_flip() + 
      ggpubr::theme_classic2()+
      chr_theme 
    
    # Add to list of plots
    chromosome_plots[[chrom]] <- p
  }
  
  return(chromosome_plots)
}



# Process the data
processed_data <- process_introgression_data(data)

# Generate chromosome lengths
chrom_lengths <- generate_chromosome_lengths()

# Create the stacked bar plots for each chromosome
chromosome_plots <- plot_introgression_stacked(processed_data, chrom_lengths)

# Arrange all chromosome plots in a grid
# all_chromosomes_plot <- plot_all_chromosomes(chromosome_plots)

# Display a specific chromosome plot (e.g., chromosome 1)
# quartz()
# print(chromosome_plots$chr1)



# Save specific chromosome plots
pdf("~/Desktop/all_chromosome_introgressions.pdf", width = 10, height = 12)
for(chrom in paste0("chr",1:10)) {
  print(chromosome_plots[[chrom]])
}
dev.off()

chr_theme2 <- theme(
  axis.line.y = element_blank(),
  axis.text.y = element_blank(),
  axis.line.x = element_blank(),
  axis.text.x = element_blank(),
  axis.title.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank()
) 



all_chromosomes_plot <- ggpubr::ggarrange(
  plotlist = lapply(1:10, function(x){
    chromosome_plots[[x]] + xlab(paste0("chr",x)) +
      coord_cartesian(expand = FALSE) + 
      coord_flip(expand = FALSE) +
      ggpubr::theme_classic2() + labs(title = NULL) + chr_theme2
  }), 
common.legend = TRUE, legend ="bottom", nrow=5,ncol=2)

# Save all chromosomes plot
ggsave("~/Desktop/introgression_all_chromosomes.png", all_chromosomes_plot, 
       width = 14, height = 10)

# Calculate chromosome summary stats
summary_by_chrom <- processed_data %>%
  group_by(chrom) %>%
  summarize(
    segments = n(),
    samples = n_distinct(sample),
    total_mb = sum(span)/1e6,
    avg_size_mb = mean(span)/1e6,
    largest_mb = max(span)/1e6
  )

# Print summary information
cat("Introgression summary by chromosome:\n")
print(summary_by_chrom)
