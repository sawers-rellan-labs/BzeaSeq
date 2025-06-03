#!/bin/bash
# calculate_variant_stats.sh
#
# Calculates statistics for the filtered teosinte variants.

INPUT_DIR="./wideseq_ref"
OUTPUT_DIR="./stats"

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Calculate total variants across all chromosomes
total_variants=0

# Process each chromosome
for chr in {1..10}; do
    echo "Calculating statistics for chromosome ${chr}..."
    
    # Input file
    INPUT_VCF="${INPUT_DIR}/wideseq_chr${chr}.vcf.gz"
    
    # Calculate statistics
    bcftools stats ${INPUT_VCF} > ${OUTPUT_DIR}/wideseq_chr${chr}.stats
    
    # Count variants
    variants=$(grep "number of SNPs:" ${OUTPUT_DIR}/wideseq_chr${chr}.stats | awk '{print $6}')
    total_variants=$((total_variants + variants))
    
    echo "Chromosome ${chr}: ${variants} variants"
done

echo "Total variants across all chromosomes: ${total_variants}"

# Calculate average variant density per 100kb
genome_size=2131846805  # Zm-B73-REFERENCE-NAM-5.0 nuclear chromosome length in bp
density_per_100kb=$(echo "scale=2; ${total_variants} * 100000 / ${genome_size}" | bc)

echo "Average variant density: ${density_per_100kb} variants per 100kb"
