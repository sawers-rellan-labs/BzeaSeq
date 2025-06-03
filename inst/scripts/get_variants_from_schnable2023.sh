#!/bin/bash
# get_variants_from_schnable2023.sh
#
# Extracts teosinte and Tripsacum samples from Schnable2023 VCF files
# and applies MAF filtering.

INPUT_DIR="schnable2023"
OUTPUT_DIR="wideseq_ref"
SAMPLE_LIST="wideseq_ref_id.list"

# This will bias against huehue 
MAF_THRESHOLD="0.05"

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Process each chromosome
for chr in {1..10}; do
    echo "Processing chromosome ${chr}..."
    
    # Input and output file paths
    INPUT_VCF="${INPUT_DIR}/schnable2023_chr${chr}.vcf.gz"
    SAMPLE_FILTER_VCF="${OUTPUT_DIR}/wideseq_taxa_chr${chr}.vcf.gz"
    OUTPUT_VCF="${OUTPUT_DIR}/wideseq_chr${chr}.vcf.gz"
    
    # Extract samples and 
    bcftools view -S ${SAMPLE_LIST} -m2 -M2 -v snps ${INPUT_VCF} -Oz -o ${SAMPLE_FILTER_VCF}
    
    # apply MAF filtering   
    bcftools view --min-af ${MAF_THRESHOLD}:minor -v snps ${SAMPLE_FILTER_VCF} -Oz -o ${OUTPUT_VCF}
    
    # Index the output file
    bcftools index ${OUTPUT_VCF}
    
    # gatk index
    gatk IndexFeatureFile -I ${OUTPUT_VCF}
    
    echo "Completed chromosome ${chr}"
done

echo "All chromosomes processed."
