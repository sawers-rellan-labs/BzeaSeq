#!/bin/bash

# Set file paths (adjust these to your file locations)
ALLELIC_COUNT_FILE="/rsstu/users/r/rrellan/BZea/bzeaseq/allelic_counts/PN15_SID1353_chr2.allelicCounts.tsv"
WIDESEQ_REF_FILE="/rsstu/users/r/rrellan/BZea/bzeaseq/wideseq_ref/wideseq_chr2.vcf.gz"
OUTPUT_DIR="./SNP_COUNT_DIFF"

# Create output directory
mkdir -p ${OUTPUT_DIR}

# 1. Extract positions from the allelic counts file (skip header line)
echo "Extracting positions from allelic counts file..."
grep -v -P "^@|^C" ${ALLELIC_COUNT_FILE}| cut -f2 > ${OUTPUT_DIR}/allelic_positions.txt
ALLELIC_COUNT=$(wc -l ${OUTPUT_DIR}/allelic_positions.txt | awk '{print $1}')
echo "Found ${ALLELIC_COUNT} positions in allelic counts file"

# 2. Extract positions from the wideseq reference VCF
echo "Extracting positions from wideseq reference VCF..."
bcftools query -f '%POS\n' ${WIDESEQ_REF_FILE} > ${OUTPUT_DIR}/wideseq_positions.txt
WIDESEQ_COUNT=$(wc -l ${OUTPUT_DIR}/wideseq_positions.txt | awk '{print $1}')
echo "Found ${WIDESEQ_COUNT} positions in wideseq reference VCF"

# 3. Find positions in allelic counts that aren't in wideseq
echo "Finding unique positions in allelic counts file..."
sort -n ${OUTPUT_DIR}/allelic_positions.txt | uniq > ${OUTPUT_DIR}/allelic_positions_sorted.txt
sort -n ${OUTPUT_DIR}/wideseq_positions.txt | uniq > ${OUTPUT_DIR}/wideseq_positions_sorted.txt
comm -23 ${OUTPUT_DIR}/allelic_positions_sorted.txt ${OUTPUT_DIR}/wideseq_positions_sorted.txt > ${OUTPUT_DIR}/positions_only_in_allelic.txt
UNIQUE_COUNT=$(wc -l ${OUTPUT_DIR}/positions_only_in_allelic.txt | awk '{print $1}')
echo "Found ${UNIQUE_COUNT} positions that are in allelic counts but not in wideseq reference"

# 4. Optional: Extract a sample of these differences for review
head -n 20 ${OUTPUT_DIR}/positions_only_in_allelic.txt > ${OUTPUT_DIR}/sample_differences.txt
echo "Sample of positions only in allelic counts file saved to ${OUTPUT_DIR}/sample_differences.txt"

# 5. Report summary
echo ""
echo "=== Summary ==="
echo "Allelic counts positions: ${ALLELIC_COUNT}"
echo "Wideseq reference positions: ${WIDESEQ_COUNT}"
echo "Difference: $((ALLELIC_COUNT - WIDESEQ_COUNT))"
echo "Positions only in allelic counts: ${UNIQUE_COUNT}"
echo "=================="

# 6. Optionally, check if these positions exist in the original filtered VCF
# This helps determine if the issue is in the SelectVariants step
if [ -n "$FILTERED_VCF" ]; then
    echo ""
    echo "Checking if unique positions exist in the original filtered VCF..."
    bcftools query -f '%POS\n' ${FILTERED_VCF} > ${OUTPUT_DIR}/filtered_positions.txt
    sort -n ${OUTPUT_DIR}/filtered_positions.txt | uniq > ${OUTPUT_DIR}/filtered_positions_sorted.txt
    
    # Check overlap with the unique positions
    comm -12 ${OUTPUT_DIR}/positions_only_in_allelic.txt ${OUTPUT_DIR}/filtered_positions_sorted.txt > ${OUTPUT_DIR}/positions_in_filtered_not_wideseq.txt
    OVERLAP_COUNT=$(wc -l ${OUTPUT_DIR}/positions_in_filtered_not_wideseq.txt | awk '{print $1}')
    echo "Of the ${UNIQUE_COUNT} unique positions, ${OVERLAP_COUNT} exist in the original filtered VCF"
fi
