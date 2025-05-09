#!/bin/bash
# todo(frodig4) This script does not do the ancestry call by bins correctly
#               Rerun on test data and correct!!!!!!
# run_wideseq.sh
#
# Description: This script collects allelic counts at teosinte-specific variant positions
# for all samples across all chromosomes using GATK's CollectAllelicCounts tool.
# It follows the workflow described in the WideSeq documentation.
# After collecting allelic counts, it performs ancestry calling on the merged results
# and cleans up temporary files including BAM files and indices.
#
# Usage: ./run_wideseq.sh [config_file]
# If no config file is specified, the script will use default values.

# activate conda environment
source ~/.bashrc
conda activate /share/maize/frodrig4/conda/env/bzeaseq

# Default configuration (can be overridden by config file)
REF_GENOME="../ref/Zm-B73-REFERENCE-NAM-5.0.fa"
BAM_BASE_DIR="/rsstu/users/r/rrellan/DOE_CAREER/BZea/mapped_bwa"
WIDESEQ_REF_DIR="./wideseq_ref"
OUTPUT_DIR="./allelic_counts"
TEMP_BAM_DIR="./bam"
SCRIPTS_DIR="./jobs"
ANCESTRY_DIR="./ancestry"
ANCESTRY_SCRIPT="./get_ancestry_calls.R"
BIN_SIZE=1000000  # Default bin size for ancestry calling (1Mb)
CHROMOSOMES=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10")
LSF_QUEUE="sara"
LSF_MEMORY="16GB"
LSF_CORES="4"

# Load configuration file if provided
if [ ! -z "$1" ]; then
    if [ -f "$1" ]; then
        echo "Loading configuration from $1"
        source "$1"
    else
        echo "Error: Configuration file $1 not found."
        exit 1
    fi
fi

# Create output and temporary directories if they don't exist
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${TEMP_BAM_DIR}"
mkdir -p "${SCRIPTS_DIR}"
mkdir -p "${ANCESTRY_DIR}"

# Function to display progress with timestamp
progress() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# Function to check if a tool is available
check_tool() {
    command -v $1 >/dev/null 2>&1 || { echo "Error: $1 is required but not installed. Aborting."; exit 1; }
}

# Check required tools
check_tool "gatk"
check_tool "picard"
check_tool "samtools"
check_tool "Rscript"

# Check if ancestry calling script exists
if [ ! -f "${ANCESTRY_SCRIPT}" ]; then
    echo "Error: Ancestry calling script ${ANCESTRY_SCRIPT} not found. Aborting."
    exit 1
fi

# Step 1: Create sequence dictionary if it doesn't exist
if [ ! -f "${REF_GENOME%.*}.dict" ]; then
    progress "Creating sequence dictionary for reference genome"
    picard CreateSequenceDictionary \
        -R="${REF_GENOME}" \
        -O="${REF_GENOME%.*}.dict"

    if [ $? -ne 0 ]; then
        echo "Error: Failed to create sequence dictionary. Aborting."
        exit 1
    else
        progress "Sequence dictionary created successfully: ${REF_GENOME%.*}.dict"
    fi
else
    progress "Sequence dictionary already exists: ${REF_GENOME%.*}.dict"
fi

# Find all sample plates (filtered_S* directories)
SAMPLE_PLATES=($(find "${BAM_BASE_DIR}" -type d -name "filtered_S*" | sort))

if [ ${#SAMPLE_PLATES[@]} -eq 0 ]; then
    echo "Error: No sample plate directories found in ${BAM_BASE_DIR}"
    exit 1
fi

progress "Found ${#SAMPLE_PLATES[@]} sample plates for processing"

# Index VCF files for all chromosomes if needed
for CHR in "${CHROMOSOMES[@]}"; do
    CHR_NUM=${CHR#chr}  # Remove 'chr' prefix to get just the number

    # Check if VCF file exists
    if [ ! -f "${WIDESEQ_REF_DIR}/wideseq_chr${CHR_NUM}.vcf.gz" ]; then
        echo "Error: VCF file ${WIDESEQ_REF_DIR}/wideseq_chr${CHR_NUM}.vcf.gz not found."
        continue
    fi

    # Check if VCF index exists
    if [ ! -f "${WIDESEQ_REF_DIR}/wideseq_chr${CHR_NUM}.vcf.gz.tbi" ]; then
        progress "Indexing VCF file for ${CHR}"
        gatk IndexFeatureFile -I "${WIDESEQ_REF_DIR}/wideseq_chr${CHR_NUM}.vcf.gz"

        if [ $? -ne 0 ]; then
            echo "Error: Failed to index VCF file for ${CHR}. Aborting."
            exit 1
        fi
    else
        progress "VCF index already exists for ${CHR}"
    fi
done

# Create a master submission script
SUBMIT_SCRIPT="${SCRIPTS_DIR}/submit_all_jobs.sh"
echo '#!/bin/bash' > "${SUBMIT_SCRIPT}"
echo '' >> "${SUBMIT_SCRIPT}"
echo '# Master script to submit all WideSeq jobs to LSF' >> "${SUBMIT_SCRIPT}"
echo '# Generated on '"$(date)" >> "${SUBMIT_SCRIPT}"
echo '' >> "${SUBMIT_SCRIPT}"

# Process each sample plate and find BAM files
TOTAL_SAMPLES=0
for PLATE_DIR in "${SAMPLE_PLATES[@]}"; do
    # Extract plate number from directory name
    PLATE_NUM=$(echo "${PLATE_DIR}" | grep -o "S[0-9]*")

    progress "Processing sample plate ${PLATE_NUM} in ${PLATE_DIR}"

    # Find all *sorted_alignment.bam files in this plate directory
    BAM_FILES=($(find "${PLATE_DIR}" -name "*sorted_alignment.bam" | grep -v "sorted_rg.bam"))

    if [ ${#BAM_FILES[@]} -eq 0 ]; then
        echo "Warning: No *sorted*alignment.bam files found in ${PLATE_DIR}"
        continue
    fi

    progress "Found ${#BAM_FILES[@]} *sorted*alignment.bam files in plate ${PLATE_NUM}"
    TOTAL_SAMPLES=$((TOTAL_SAMPLES + ${#BAM_FILES[@]}))

    # Process each BAM file in this plate - create individual job scripts
    for BAM_FILE in "${BAM_FILES[@]}"; do
        # Extract sample name from BAM filename
        SAMPLE=$(basename "${BAM_FILE}" .bam)
        SAMPLE=${SAMPLE/_sorted_alignment/}  # Remove _sorted_alignment suffix if present

        progress "Creating job script for sample: ${SAMPLE} (Plate ${PLATE_NUM})"

        # Create job script for this sample
        JOB_SCRIPT="${SCRIPTS_DIR}/wideseq_${SAMPLE}.sh"

        cat > "${JOB_SCRIPT}" << EOF
#!/bin/bash

#BSUB -J wideseq_${SAMPLE}
#BSUB -q ${LSF_QUEUE}
#BSUB -n ${LSF_CORES}
#BSUB -W 180
#BSUB -R "rusage[mem=${LSF_MEMORY}]"
#BSUB -o ${SCRIPTS_DIR}/wideseq_${SAMPLE}_%J.out
#BSUB -e ${SCRIPTS_DIR}/wideseq_${SAMPLE}_%J.err

# activate conda environment
source ~/.bashrc
conda activate /share/maize/frodrig4/conda/env/bzeaseq

# WideSeq Pipeline for ${SAMPLE} (Plate ${PLATE_NUM})
# Generated on $(date)
# This script performs both allelic count collection and ancestry calling

# Configuration
REF_GENOME="${REF_GENOME}"
BAM_FILE="${BAM_FILE}"
WIDESEQ_REF_DIR="${WIDESEQ_REF_DIR}"
OUTPUT_DIR="${OUTPUT_DIR}"
TEMP_BAM_DIR="${TEMP_BAM_DIR}"
ANCESTRY_DIR="${ANCESTRY_DIR}"
ANCESTRY_SCRIPT="${ANCESTRY_SCRIPT}"
BIN_SIZE="${BIN_SIZE}"
SAMPLE="${SAMPLE}"
PLATE_NUM="${PLATE_NUM}"
CHROMOSOMES=($(printf '"%s" ' "${CHROMOSOMES[@]}"))

# Function to display progress with timestamp
progress() {
    echo "\$(date '+%Y-%m-%d %H:%M:%S') - \$1"
}

progress "Starting WideSeq pipeline for \${SAMPLE} (Plate \${PLATE_NUM})"

# Step 1: Add read groups if not already present
if ! samtools view -H "\${BAM_FILE}" | grep -q "@RG"; then
    progress "Adding read groups to \${SAMPLE} BAM file"
    picard AddOrReplaceReadGroups \\
        -I "\${BAM_FILE}" \\
        -O "\${TEMP_BAM_DIR}/\${SAMPLE}_sorted_rg.bam" \\
        --RGID "\${PLATE_NUM}" \\
        --RGPL illumina \\
        --RGLB "\${SAMPLE}" \\
        --RGPU NONE \\
        --RGSM "\${SAMPLE}"

    if [ \$? -ne 0 ]; then
        echo "Error: Failed to add read groups for \${SAMPLE}. Aborting."
        exit 1
    else
        progress "Read groups added successfully for \${SAMPLE}"
    fi

    # Step 2: Index the read group BAM
    progress "Indexing read group BAM for \${SAMPLE}"
    samtools index "\${TEMP_BAM_DIR}/\${SAMPLE}_sorted_rg.bam"

    if [ \$? -ne 0 ]; then
        echo "Error: Failed to index BAM file for \${SAMPLE}. Aborting."
        exit 1
    fi

    # Use the new BAM file for downstream processing
    PROCESSED_BAM="\${TEMP_BAM_DIR}/\${SAMPLE}_sorted_rg.bam"
else
    progress "Read groups already present in \${SAMPLE} BAM file"

    # Check if BAM is indexed
    if [ ! -f "\${BAM_FILE}.bai" ]; then
        progress "Indexing BAM file for \${SAMPLE}"
        samtools index "\${BAM_FILE}"

        if [ \$? -ne 0 ]; then
            echo "Error: Failed to index BAM file for \${SAMPLE}. Aborting."
            exit 1
        fi
    fi

    # Use the original BAM file for downstream processing
    PROCESSED_BAM="\${BAM_FILE}"
fi

# Create merged output file
MERGED_OUTPUT="\${OUTPUT_DIR}/\${SAMPLE}.allelicCounts.tsv"

# Process all chromosomes for this sample
for CHR in "\${CHROMOSOMES[@]}"; do
    CHR_NUM=\${CHR#chr}  # Remove 'chr' prefix to get just the number

    progress "Processing chromosome \${CHR} for sample \${SAMPLE}"

    # Step 3: Collect allelic counts
    progress "Collecting allelic counts for \${SAMPLE} on \${CHR}"
    gatk CollectAllelicCounts \\
        -I "\${PROCESSED_BAM}" \\
        -R "\${REF_GENOME}" \\
        -L "\${WIDESEQ_REF_DIR}/wideseq_chr\${CHR_NUM}.vcf.gz" \\
        -O "\${OUTPUT_DIR}/\${SAMPLE}_chr\${CHR_NUM}.allelicCounts.tsv"

    if [ \$? -ne 0 ]; then
        echo "Error: Failed to collect allelic counts for \${SAMPLE} on \${CHR}. Aborting."
        exit 1
    else
        progress "Completed allelic count collection for \${SAMPLE} on \${CHR}"
    fi
done

# Step 4: Merge all chromosome results
progress "Merging allelic counts from all chromosomes for \${SAMPLE}"

# Get the header from the first chromosome file
if [ -f "\${OUTPUT_DIR}/\${SAMPLE}_chr1.allelicCounts.tsv" ]; then
    # Extract header (lines starting with @)
    grep "^@" "\${OUTPUT_DIR}/\${SAMPLE}_chr1.allelicCounts.tsv" > "\${MERGED_OUTPUT}"

    # Add the column header line (first non-@ line)
    grep -v "^@" "\${OUTPUT_DIR}/\${SAMPLE}_chr1.allelicCounts.tsv" | head -n 1 >> "\${MERGED_OUTPUT}"
else
    echo "Error: No output file found for chromosome 1. Cannot create merged file."
    exit 1
fi

# Append data from all chromosome files (skipping headers)
for CHR in "\${CHROMOSOMES[@]}"; do
    CHR_NUM=\${CHR#chr}
    if [ -f "\${OUTPUT_DIR}/\${SAMPLE}_chr\${CHR_NUM}.allelicCounts.tsv" ]; then
        grep -v "^@" "\${OUTPUT_DIR}/\${SAMPLE}_chr\${CHR_NUM}.allelicCounts.tsv" | tail -n +2 >> "\${MERGED_OUTPUT}"
    fi
done

# Check if merge was successful
if [ ! -f "\${MERGED_OUTPUT}" ]; then
    echo "Error: Failed to create merged output file for \${SAMPLE}."
    exit 1
fi

# Calculate some statistics for allelic counts
VARIANT_POSITIONS=\$(grep -v "^@" "\${MERGED_OUTPUT}" | wc -l)
VARIANT_POSITIONS=\$((VARIANT_POSITIONS - 1))  # Subtract 1 for header

echo "=== WideSeq Allelic Count Collection Summary for \${SAMPLE} ==="
echo "Sample name: \${SAMPLE}"
echo "Plate: \${PLATE_NUM}"
echo "Total positions collected: \${VARIANT_POSITIONS}"
echo "Results available in: \${MERGED_OUTPUT}"
echo "=========================================================="

# Cleanup temporary BAM files and indexes after allelic counts are collected
progress "Cleaning up temporary BAM files and indexes"

# If we created a temporary BAM with read groups, delete it
if [ "\${PROCESSED_BAM}" != "\${BAM_FILE}" ]; then
    # Check if the temporary BAM was in the temp directory
    if [[ "\${PROCESSED_BAM}" == "\${TEMP_BAM_DIR}"* ]]; then
        progress "Removing temporary BAM file: \${PROCESSED_BAM}"
        rm -f "\${PROCESSED_BAM}"
        rm -f "\${PROCESSED_BAM}.bai"
    fi
fi

# Step 5: Run ancestry calling on the merged file
progress "Running ancestry calling on merged allelic counts for \${SAMPLE}"
SAMPLE_PREFIX="\${ANCESTRY_DIR}/\${SAMPLE}"
Rscript "\${ANCESTRY_SCRIPT}" "\${MERGED_OUTPUT}" "\${SAMPLE_PREFIX}" "\${BIN_SIZE}"

if [ \$? -ne 0 ]; then
    echo "Error: Ancestry calling failed for \${SAMPLE}. Check logs for details."
    exit 1
else
    progress "Ancestry calling completed successfully for \${SAMPLE}"
fi

# Clean up if ancestry calling was successful
if [ -f "\${SAMPLE_PREFIX}_bin_genotypes.tsv" ]; then
    progress "Cleaning up intermediate files"
    
    # Clean up per-chromosome allelic count files
    for CHR in "\${CHROMOSOMES[@]}"; do
        CHR_NUM=\${CHR#chr}
        # Remove per-chromosome allelic count files
        rm -f "\${OUTPUT_DIR}/\${SAMPLE}_chr\${CHR_NUM}.allelicCounts.tsv"
        # Remove intermediate VCF files
        rm -f "\${OUTPUT_DIR}/\${SAMPLE}_wideseq_selected_chr\${CHR_NUM}.vcf"
        rm -f "\${OUTPUT_DIR}/\${SAMPLE}_wideseq_selected_chr\${CHR_NUM}.vcf.idx"
    done
    
    # Calculate and display ancestry statistics
    if [ -f "\${SAMPLE_PREFIX}_summary.tsv" ]; then
        REF_PCT=\$(awk -F'\t' '{print \$6}' "\${SAMPLE_PREFIX}_summary.tsv" | tail -n 1)
        HET_PCT=\$(awk -F'\t' '{print \$7}' "\${SAMPLE_PREFIX}_summary.tsv" | tail -n 1)
        ALT_PCT=\$(awk -F'\t' '{print \$8}' "\${SAMPLE_PREFIX}_summary.tsv" | tail -n 1)
        NON_REF_PCT=\$(awk -F'\t' '{print \$9}' "\${SAMPLE_PREFIX}_summary.tsv" | tail -n 1)
        
        echo "=== WideSeq Ancestry Calling Summary for \${SAMPLE} ==="
        echo "Sample name: \${SAMPLE}"
        echo "Plate: \${PLATE_NUM}"
        echo "Bin size: \${BIN_SIZE} bp"
        echo "REF bins: \${REF_PCT}%"
        echo "HET bins: \${HET_PCT}%"
        echo "ALT bins: \${ALT_PCT}%"
        echo "NON-REF bins: \${NON_REF_PCT}%"
        echo "Results available in: \${SAMPLE_PREFIX}_*"
        echo "=========================================================="
    fi
else
    echo "Warning: Ancestry calling output not found. Intermediate files not cleaned up."
fi

progress "WideSeq pipeline completed for \${SAMPLE} (Plate \${PLATE_NUM})"
EOF

        # Make the job script executable
        chmod +x "${JOB_SCRIPT}"

        # Add this job to the master submission script
        echo "bsub < ${JOB_SCRIPT}" >> "${SUBMIT_SCRIPT}"

        progress "Created job script for ${SAMPLE} (Plate ${PLATE_NUM}): ${JOB_SCRIPT}"
    done
done

# Make the submission script executable
chmod +x "${SUBMIT_SCRIPT}"

# Print summary
echo ""
echo "=== WideSeq Pipeline Summary ==="
echo "Total sample plates found: ${#SAMPLE_PLATES[@]}"
echo "Total samples: ${TOTAL_SAMPLES}"
echo "Job scripts created in: ${SCRIPTS_DIR}"
echo "To submit all jobs, run: ${SUBMIT_SCRIPT}"
echo "===================================="
