#!/bin/bash

#BSUB -J wideseq_analysis
#BSUB -n 4
#BSUB -q normal
#BSUB -R "rusage[mem=16GB]"
#BSUB -o wideseq_%J.out
#BSUB -e wideseq_%J.err

# Set variables
FASTQ_DIR="/path/to/wideseq_fastq"
REFERENCE="/path/to/B73_v4_reference.fasta"
OUTPUT_DIR="/path/to/wideseq_output"
EXTENDED_HAPMAP="/path/to/expanded_hapmap3.vcf.gz"  # From previous pipeline
HAPMAP_PANEL="/path/to/hapmap3_panel.vcf.gz"
THREADS=4
BIN_SIZE=1000000

# Create output directories
mkdir -p ${OUTPUT_DIR}/{trimmed,mapping,variants,hapmap_filtered,bins,haplotypes,merged,logs}

# List of chromosomes in maize
CHROMOSOMES=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10")

# 1. Prepare chromosome-specific reference files (once)
bsub -J "prep_ref" \
     -n 2 \
     -R "rusage[mem=8GB]" \
     -o ${OUTPUT_DIR}/logs/prep_ref.out \
     -e ${OUTPUT_DIR}/logs/prep_ref.err \
     "
     # Index main reference
     bwa index ${REFERENCE}
     samtools faidx ${REFERENCE}
     
     # Create chromosome-specific reference files
     for chrom in ${CHROMOSOMES[@]}; do
         samtools faidx ${REFERENCE} \${chrom} > ${OUTPUT_DIR}/mapping/\${chrom}.fasta
         bwa index ${OUTPUT_DIR}/mapping/\${chrom}.fasta
         samtools faidx ${OUTPUT_DIR}/mapping/\${chrom}.fasta
     done
     
     # Create genome file for bedtools
     cut -f1,2 ${REFERENCE}.fai > ${OUTPUT_DIR}/B73_v4.genome
     "

# 2. Trim reads (per sample)
for R1 in ${FASTQ_DIR}/*_R1.fastq.gz; do
    sample=$(basename $R1 _R1.fastq.gz)
    R2="${FASTQ_DIR}/${sample}_R2.fastq.gz"
    
    bsub -J "trim_${sample}" \
         -n 4 \
         -R "rusage[mem=8GB]" \
         -o ${OUTPUT_DIR}/logs/trim_${sample}.out \
         -e ${OUTPUT_DIR}/logs/trim_${sample}.err \
         "
         trimmomatic PE -threads ${THREADS} \
             ${R1} ${R2} \
             ${OUTPUT_DIR}/trimmed/${sample}_R1_trimmed.fastq.gz \
             ${OUTPUT_DIR}/trimmed/${sample}_R1_unpaired.fastq.gz \
             ${OUTPUT_DIR}/trimmed/${sample}_R2_trimmed.fastq.gz \
             ${OUTPUT_DIR}/trimmed/${sample}_R2_unpaired.fastq.gz \
             ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 \
             SLIDINGWINDOW:4:15 MINLEN:36
         "
done

# 3. Map reads per chromosome (per sample)
for R1 in ${FASTQ_DIR}/*_R1.fastq.gz; do
    sample=$(basename $R1 _R1.fastq.gz)
    
    # Submit job array for mapping to each chromosome
    bsub -J "map_${sample}[1-10]" \
         -w "done(trim_${sample}) && done(prep_ref)" \
         -n 4 \
         -R "rusage[mem=16GB]" \
         -o ${OUTPUT_DIR}/logs/map_${sample}_%I.out \
         -e ${OUTPUT_DIR}/logs/map_${sample}_%I.err \
         "
         # Get chromosome from job array index
         idx=\$((LSB_JOBINDEX - 1))
         chrom=${CHROMOSOMES[\$idx]}
         
         # Map reads to chromosome-specific reference
         bwa mem -t ${THREADS} ${OUTPUT_DIR}/mapping/\${chrom}.fasta \
             ${OUTPUT_DIR}/trimmed/${sample}_R1_trimmed.fastq.gz \
             ${OUTPUT_DIR}/trimmed/${sample}_R2_trimmed.fastq.gz | \
             samtools view -b -o ${OUTPUT_DIR}/mapping/${sample}.\${chrom}.raw.bam -
         
         # Sort and mark duplicates
         samtools sort -@ ${THREADS} -o ${OUTPUT_DIR}/mapping/${sample}.\${chrom}.sorted.bam \
             ${OUTPUT_DIR}/mapping/${sample}.\${chrom}.raw.bam
         
         samtools markdup -@ ${THREADS} ${OUTPUT_DIR}/mapping/${sample}.\${chrom}.sorted.bam \
             ${OUTPUT_DIR}/mapping/${sample}.\${chrom}.bam
         
         samtools index ${OUTPUT_DIR}/mapping/${sample}.\${chrom}.bam
         
         # Clean up intermediate files
         rm ${OUTPUT_DIR}/mapping/${sample}.\${chrom}.raw.bam ${OUTPUT_DIR}/mapping/${sample}.\${chrom}.sorted.bam
         "
done

# 4. Call variants per chromosome using bcftools mpileup
for R1 in ${FASTQ_DIR}/*_R1.fastq.gz; do
    sample=$(basename $R1 _R1.fastq.gz)
    
    bsub -J "call_${sample}[1-10]" \
         -w "done(map_${sample}*)" \
         -n 4 \
         -R "rusage[mem=16GB]" \
         -o ${OUTPUT_DIR}/logs/call_${sample}_%I.out \
         -e ${OUTPUT_DIR}/logs/call_${sample}_%I.err \
         "
         # Get chromosome from job array index
         idx=\$((LSB_JOBINDEX - 1))
         chrom=${CHROMOSOMES[\$idx]}
         
         # Call variants using mpileup
         bcftools mpileup -Ou -f ${OUTPUT_DIR}/mapping/\${chrom}.fasta \
             ${OUTPUT_DIR}/mapping/${sample}.\${chrom}.bam | \
             bcftools call -mv -Oz -o ${OUTPUT_DIR}/variants/${sample}.\${chrom}.vcf.gz
         
         bcftools index ${OUTPUT_DIR}/variants/${sample}.\${chrom}.vcf.gz
         "
done

# 5. Filter variants using extended HapMap3 per chromosome
for R1 in ${FASTQ_DIR}/*_R1.fastq.gz; do
    sample=$(basename $R1 _R1.fastq.gz)
    
    bsub -J "filter_${sample}[1-10]" \
         -w "done(call_${sample}*)" \
         -n 4 \
         -R "rusage[mem=16GB]" \
         -o ${OUTPUT_DIR}/logs/filter_${sample}_%I.out \
         -e ${OUTPUT_DIR}/logs/filter_${sample}_%I.err \
         "
         # Get chromosome from job array index
         idx=\$((LSB_JOBINDEX - 1))
         chrom=${CHROMOSOMES[\$idx]}
         
         # Extract chromosome from extended HapMap
         bcftools view -r \${chrom} ${EXTENDED_HAPMAP} -Oz -o ${OUTPUT_DIR}/hapmap_filtered/hapmap3.\${chrom}.vcf.gz
         bcftools index ${OUTPUT_DIR}/hapmap_filtered/hapmap3.\${chrom}.vcf.gz
         
         # Intersect with HapMap positions
         bcftools view -R ${OUTPUT_DIR}/hapmap_filtered/hapmap3.\${chrom}.vcf.gz \
             ${OUTPUT_DIR}/variants/${sample}.\${chrom}.vcf.gz \
             -Oz -o ${OUTPUT_DIR}/hapmap_filtered/${sample}.\${chrom}.hapmap.vcf.gz
         
         bcftools index ${OUTPUT_DIR}/hapmap_filtered/${sample}.\${chrom}.hapmap.vcf.gz
         
         # Extract non-reference positions
         bcftools view -i 'GT=\"alt\"' ${OUTPUT_DIR}/hapmap_filtered/${sample}.\${chrom}.hapmap.vcf.gz \
             -Oz -o ${OUTPUT_DIR}/hapmap_filtered/${sample}.\${chrom}.nonref.vcf.gz
         
         bcftools index ${OUTPUT_DIR}/hapmap_filtered/${sample}.\${chrom}.nonref.vcf.gz
         
         # Convert to tabix format for bin calculation
         bcftools query -f '%CHROM\\t%POS0\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' \
             ${OUTPUT_DIR}/hapmap_filtered/${sample}.\${chrom}.hapmap.vcf.gz \
             > ${OUTPUT_DIR}/hapmap_filtered/${sample}.\${chrom}.bed
         "
done

# 6. Create bins and calculate bin frequencies per chromosome
for R1 in ${FASTQ_DIR}/*_R1.fastq.gz; do
    sample=$(basename $R1 _R1.fastq.gz)
    
    bsub -J "bins_${sample}[1-10]" \
         -w "done(filter_${sample}*)" \
         -n 4 \
         -R "rusage[mem=16GB]" \
         -o ${OUTPUT_DIR}/logs/bins_${sample}_%I.out \
         -e ${OUTPUT_DIR}/logs/bins_${sample}_%I.err \
         "
         # Get chromosome from job array index
         idx=\$((LSB_JOBINDEX - 1))
         chrom=${CHROMOSOMES[\$idx]}
         
         # Create bins for this chromosome
         chr_size=\$(grep \"^\${chrom}\" ${OUTPUT_DIR}/B73_v4.genome | cut -f2)
         bedtools makewindows -g <(echo -e \"\${chrom}\\t\${chr_size}\") -w ${BIN_SIZE} \
             > ${OUTPUT_DIR}/bins/\${chrom}_bins_${BIN_SIZE}bp.bed
         
         # Calculate reference and non-reference read counts in each bin
         bedtools map -a ${OUTPUT_DIR}/bins/\${chrom}_bins_${BIN_SIZE}bp.bed \
             -b ${OUTPUT_DIR}/hapmap_filtered/${sample}.\${chrom}.bed \
             -c 4,5 -o count,count > ${OUTPUT_DIR}/bins/${sample}.\${chrom}.bin_counts.txt
         
         # Calculate frequencies using Python
         python3 -c \"
import pandas as pd

# Read bin counts
counts = pd.read_csv('${OUTPUT_DIR}/bins/${sample}.\${chrom}.bin_counts.txt', sep='\\t', 
                    header=None, names=['chrom', 'start', 'end', 'ref_count', 'alt_count'])

# Calculate frequencies
counts['total'] = counts['ref_count'] + counts['alt_count']
counts['ref_freq'] = counts['ref_count'] / counts['total'].replace(0, float('nan'))
counts['alt_freq'] = counts['alt_count'] / counts['total'].replace(0, float('nan'))

# Save results
counts.to_csv('${OUTPUT_DIR}/bins/${sample}.\${chrom}.bin_frequencies.txt', sep='\\t', index=False)
\"
         "
done

# 7. Calculate Jaccard Index with HapMap panel per chromosome
for R1 in ${FASTQ_DIR}/*_R1.fastq.gz; do
    sample=$(basename $R1 _R1.fastq.gz)
    
    bsub -J "haplo_${sample}[1-10]" \
         -w "done(filter_${sample}*)" \
         -n 4 \
         -R "rusage[mem=16GB]" \
         -o ${OUTPUT_DIR}/logs/haplo_${sample}_%I.out \
         -e ${OUTPUT_DIR}/logs/haplo_${sample}_%I.err \
         "
         # Get chromosome from job array index
         idx=\$((LSB_JOBINDEX - 1))
         chrom=${CHROMOSOMES[\$idx]}
         
         # Extract positions of non-reference variants
         bcftools query -f '%CHROM\\t%POS\\n' \
             ${OUTPUT_DIR}/hapmap_filtered/${sample}.\${chrom}.nonref.vcf.gz \
             > ${OUTPUT_DIR}/haplotypes/${sample}.\${chrom}.positions.txt
         
         # Extract these positions from HapMap3 panel for this chromosome
         bcftools view -r \${chrom} ${HAPMAP_PANEL} \
             -Oz -o ${OUTPUT_DIR}/haplotypes/hapmap_panel.\${chrom}.vcf.gz
         bcftools index ${OUTPUT_DIR}/haplotypes/hapmap_panel.\${chrom}.vcf.gz
         
         bcftools view -R ${OUTPUT_DIR}/haplotypes/${sample}.\${chrom}.positions.txt \
             ${OUTPUT_DIR}/haplotypes/hapmap_panel.\${chrom}.vcf.gz \
             -Oz -o ${OUTPUT_DIR}/haplotypes/${sample}.\${chrom}.hapmap_panel_subset.vcf.gz
         
         # Filter by MAF
         bcftools view -q 0.05:minor \
             ${OUTPUT_DIR}/haplotypes/${sample}.\${chrom}.hapmap_panel_subset.vcf.gz \
             -Oz -o ${OUTPUT_DIR}/haplotypes/${sample}.\${chrom}.hapmap_panel_filtered.vcf.gz
         
         bcftools index ${OUTPUT_DIR}/haplotypes/${sample}.\${chrom}.hapmap_panel_filtered.vcf.gz
         
         # Run Jaccard calculation with Python
         python3 -c \"
import pandas as pd
import numpy as np
from sklearn.metrics import jaccard_score
import subprocess
import io

# Read sample variants - encode as 0/1
sample_cmd = 'bcftools query -f \\\"%CHROM\\\\t%POS\\\\t%REF\\\\t%ALT[\\\\t%GT]\\\\n\\\" ${OUTPUT_DIR}/hapmap_filtered/${sample}.\${chrom}.nonref.vcf.gz'
sample_vars = pd.read_csv(io.StringIO(subprocess.check_output(sample_cmd, shell=True).decode()), 
                         sep='\\t', header=None, names=['chrom', 'pos', 'ref', 'alt', 'gt'])
sample_vars['value'] = 1  # Non-reference positions set to 1

# Get list of accessions in HapMap3 panel
panel_cmd = 'bcftools query -l ${OUTPUT_DIR}/haplotypes/${sample}.\${chrom}.hapmap_panel_filtered.vcf.gz'
accessions = subprocess.check_output(panel_cmd, shell=True).decode().strip().split('\\n')

# Calculate Jaccard indices
results = []
for acc in accessions:
    # Extract genotypes for this accession
    acc_cmd = f'bcftools query -s {acc} -f \\\"%CHROM\\\\t%POS\\\\t%REF\\\\t%ALT[\\\\t%GT]\\\\n\\\" ${OUTPUT_DIR}/haplotypes/${sample}.\${chrom}.hapmap_panel_filtered.vcf.gz'
    acc_vars = pd.read_csv(io.StringIO(subprocess.check_output(acc_cmd, shell=True).decode()),
                          sep='\\t', header=None, names=['chrom', 'pos', 'ref', 'alt', 'gt'])
    
    # Convert genotypes to binary (1 for alt allele present, 0 for ref only)
    acc_vars['value'] = acc_vars['gt'].apply(lambda x: 1 if '1' in x else 0)
    
    # Merge and calculate Jaccard index
    merged = pd.merge(sample_vars, acc_vars, on=['chrom', 'pos'], suffixes=('_sample', '_acc'))
    if merged.shape[0] > 0:
        jaccard = jaccard_score(merged['value_sample'], merged['value_acc'])
    else:
        jaccard = 0
    
    results.append({'accession': acc, 'chromosome': '\${chrom}', 'jaccard_index': jaccard})

# Save results
pd.DataFrame(results).sort_values('jaccard_index', ascending=False).to_csv(
    '${OUTPUT_DIR}/haplotypes/${sample}.\${chrom}.jaccard_similarities.txt', sep='\\t', index=False)
\"
         "
done

# 8. Merge results across chromosomes
for R1 in ${FASTQ_DIR}/*_R1.fastq.gz; do
    sample=$(basename $R1 _R1.fastq.gz)
    
    bsub -J "merge_${sample}" \
         -w "done(bins_${sample}*) && done(haplo_${sample}*)" \
         -n 4 \
         -R "rusage[mem=16GB]" \
         -o ${OUTPUT_DIR}/logs/merge_${sample}.out \
         -e ${OUTPUT_DIR}/logs/merge_${sample}.err \
         "
         # Merge bin frequencies
         awk 'FNR>1 || NR==1' ${OUTPUT_DIR}/bins/${sample}.chr*.bin_frequencies.txt > ${OUTPUT_DIR}/merged/${sample}.all_bins.txt
         
         # Merge Jaccard similarities (keeping top matches per chromosome)
         for chrom in ${CHROMOSOMES[@]}; do
             head -n 11 ${OUTPUT_DIR}/haplotypes/${sample}.\${chrom}.jaccard_similarities.txt >> ${OUTPUT_DIR}/merged/${sample}.top_matches.txt
         done
         
         # Merge VCFs
         chr_vcfs=\$(ls ${OUTPUT_DIR}/hapmap_filtered/${sample}.chr*.hapmap.vcf.gz)
         bcftools concat \$chr_vcfs -Oz -o ${OUTPUT_DIR}/merged/${sample}.all_hapmap.vcf.gz
         bcftools index ${OUTPUT_DIR}/merged/${sample}.all_hapmap.vcf.gz
         
         # Create comprehensive report
         echo \"WideSeq Analysis Report for ${sample}\" > ${OUTPUT_DIR}/merged/${sample}_report.txt
         echo \"===================================\" >> ${OUTPUT_DIR}/merged/${sample}_report.txt
         echo \"\" >> ${OUTPUT_DIR}/merged/${sample}_report.txt
         
         echo \"Top Haplotype Matches by Chromosome:\" >> ${OUTPUT_DIR}/merged/${sample}_report.txt
         echo \"-----------------------------------\" >> ${OUTPUT_DIR}/merged/${sample}_report.txt
         for chrom in ${CHROMOSOMES[@]}; do
             echo \"\\n\${chrom}:\" >> ${OUTPUT_DIR}/merged/${sample}_report.txt
             head -n 6 ${OUTPUT_DIR}/haplotypes/${sample}.\${chrom}.jaccard_similarities.txt | tail -n 5 >> ${OUTPUT_DIR}/merged/${sample}_report.txt
         done
         
         echo \"\\nVariant Statistics:\" >> ${OUTPUT_DIR}/merged/${sample}_report.txt
         echo \"------------------\" >> ${OUTPUT_DIR}/merged/${sample}_report.txt
         for chrom in ${CHROMOSOMES[@]}; do
             count=\$(bcftools view -H ${OUTPUT_DIR}/hapmap_filtered/${sample}.\${chrom}.hapmap.vcf.gz | wc -l)
             echo \"\${chrom}: \$count SNPs\" >> ${OUTPUT_DIR}/merged/${sample}_report.txt
         done
         total=\$(bcftools view -H ${OUTPUT_DIR}/merged/${sample}.all_hapmap.vcf.gz | wc -l)
         echo \"Total: \$total SNPs\" >> ${OUTPUT_DIR}/merged/${sample}_report.txt
         "
done

# 9. Create visualization of bin frequencies (using R)
for R1 in ${FASTQ_DIR}/*_R1.fastq.gz; do
    sample=$(basename $R1 _R1.fastq.gz)
    
    bsub -J "viz_${sample}" \
         -w "done(merge_${sample})" \
         -n 4 \
         -R "rusage[mem=16GB]" \
         -o ${OUTPUT_DIR}/logs/viz_${sample}.out \
         -e ${OUTPUT_DIR}/logs/viz_${sample}.err \
         "
         # Create R script for visualization
         cat > ${OUTPUT_DIR}/merged/plot_${sample}.R << 'EOF'
# Load libraries
library(ggplot2)
library(data.table)
library(gridExtra)

# Read merged bin data
sample_name <- \"${sample}\"
bins <- fread(\"${OUTPUT_DIR}/merged/${sample}.all_bins.txt\")

# Plot reference allele frequency across all chromosomes
p <- ggplot(bins, aes(# Create R script for visualization
         cat > ${OUTPUT_DIR}/merged/plot_${sample}.R << 'EOF'
# Load libraries
library(ggplot2)
library(data.table)
library(gridExtra)

# Read merged bin data
sample_name <- "${sample}"
bins <- fread("${OUTPUT_DIR}/merged/${sample}.all_bins.txt")

# Plot reference allele frequency across all chromosomes
p <- ggplot(bins, aes(x=start/1e6, y=ref_freq)) +
  geom_line() +
  geom_smooth(span=0.2, se=FALSE, color="blue") +
  facet_wrap(~chrom, scales="free_x", ncol=2) +
  labs(x="Position (Mb)", y="Reference Allele Frequency",
       title=paste("Reference Allele Frequency for", sample_name)) +
  theme_minimal() +
  theme(strip.background=element_rect(fill="lightblue", color="black"),
        strip.text=element_text(face="bold"))

# Save plot
ggsave(paste0("${OUTPUT_DIR}/merged/", sample_name, "_allele_freq_plot.pdf"), 
       p, width=12, height=10)

# Create heatmap-style visualization of all chromosomes
# Reshape data for heatmap
bins$bin_id <- 1:nrow(bins)
bins$alt_freq[is.na(bins$alt_freq)] <- 0

# Cap at 100 bins per chromosome for visualization
heat_data <- list()
for (chr in unique(bins$chrom)) {
  chr_bins <- bins[bins$chrom == chr,]
  n_bins <- nrow(chr_bins)
  
  if (n_bins > 100) {
    # Aggregate bins if more than 100
    bin_size <- ceiling(n_bins/100)
    agg_bins <- data.table()
    for (i in 1:ceiling(n_bins/bin_size)) {
      start_idx <- (i-1)*bin_size + 1
      end_idx <- min(i*bin_size, n_bins)
      if (start_idx <= n_bins) {
        sub_bins <- chr_bins[start_idx:end_idx,]
        agg_bins <- rbind(agg_bins, data.table(
          chrom = chr,
          bin_group = i,
          start = min(sub_bins$start),
          end = max(sub_bins$end),
          alt_freq = mean(sub_bins$alt_freq, na.rm=TRUE)
        ))
      }
    }
    heat_data[[chr]] <- agg_bins
  } else {
    # Keep as is if fewer than 100 bins
    chr_bins$bin_group <- 1:n_bins
    heat_data[[chr]] <- chr_bins[,.(chrom, bin_group, start, end, alt_freq)]
  }
}

# Combine all chromosomes
heat_data <- rbindlist(heat_data)

# Create heatmap
h <- ggplot(heat_data, aes(x=bin_group, y=chrom, fill=alt_freq)) +
  geom_tile() +
  scale_fill_gradientn(
    colors=c("navy", "white", "firebrick"), 
    values=c(0, 0.5, 1),
    name="Alt Allele\nFrequency") +
  labs(x="Bin Position", y="Chromosome",
       title=paste("Genome-wide Ancestry Heatmap for", sample_name)) +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Save heatmap
ggsave(paste0("${OUTPUT_DIR}/merged/", sample_name, "_ancestry_heatmap.pdf"), 
       h, width=10, height=8)

# Create top matches visualization
# Read Jaccard data
jaccard_data <- fread("${OUTPUT_DIR}/merged/${sample}.top_matches.txt")

# Keep top 15 unique accessions
top_matches <- jaccard_data[order(-jaccard_index)]
top_matches <- top_matches[!duplicated(accession)][1:15]

# Plot top matches
j <- ggplot(top_matches, aes(x=reorder(accession, jaccard_index), y=jaccard_index)) +
  geom_bar(stat="identity", fill="steelblue") +
  geom_text(aes(label=round(jaccard_index, 3)), hjust=-0.1, size=3) +
  labs(x="Accession", y="Jaccard Similarity Index",
       title=paste("Top Haplotype Matches for", sample_name)) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.y=element_text(size=8))

# Save plot
ggsave(paste0("${OUTPUT_DIR}/merged/", sample_name, "_top_matches.pdf"), 
       j, width=8, height=6)
EOF

         # Run R script
         Rscript ${OUTPUT_DIR}/merged/plot_${sample}.R
         "
done

# 10. Create final summary report across all samples
bsub -J "final_summary" \
     -w "done(viz_*)" \
     -n 4 \
     -R "rusage[mem=16GB]" \
     -o ${OUTPUT_DIR}/logs/final_summary.out \
     -e ${OUTPUT_DIR}/logs/final_summary.err \
     "
     # Generate summary of all samples
     echo 'WideSeq Analysis Summary' > ${OUTPUT_DIR}/merged/all_samples_summary.txt
     echo '=======================' >> ${OUTPUT_DIR}/merged/all_samples_summary.txt
     echo '' >> ${OUTPUT_DIR}/merged/all_samples_summary.txt
     
     echo 'Sample Statistics:' >> ${OUTPUT_DIR}/merged/all_samples_summary.txt
     echo '-----------------' >> ${OUTPUT_DIR}/merged/all_samples_summary.txt
     
     # Process each sample
     for R1 in ${FASTQ_DIR}/*_R1.fastq.gz; do
         sample=\$(basename \$R1 _R1.fastq.gz)
         
         # Get total SNPs
         total_snps=\$(bcftools view -H ${OUTPUT_DIR}/merged/\${sample}.all_hapmap.vcf.gz 2>/dev/null | wc -l)
         
         # Get top match across all chromosomes
         top_match=\$(sort -t \$'\\t' -k3,3nr ${OUTPUT_DIR}/merged/\${sample}.top_matches.txt | head -n2 | tail -n1)
         
         echo \"\${sample}:\" >> ${OUTPUT_DIR}/merged/all_samples_summary.txt
         echo \"  - Total SNPs: \${total_snps}\" >> ${OUTPUT_DIR}/merged/all_samples_summary.txt
         echo \"  - Top Match: \${top_match}\" >> ${OUTPUT_DIR}/merged/all_samples_summary.txt
         echo \"\" >> ${OUTPUT_DIR}/merged/all_samples_summary.txt
     done
     
     echo 'All analysis complete. Summary reports and visualizations available in ${OUTPUT_DIR}/merged/'
     "