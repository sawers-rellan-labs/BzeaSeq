#!/bin/bash
# ./generate_wideseq_jobs.sh [config_file]
# WideSeq Job Generator for LSF
# This script generates individual job scripts for each sample and creates
# a submission script to send all jobs to the LSF scheduler.
# It accounts for different sample plate numbers in the BAM file structure.
#
# Usage: ./generate_wideseq_jobs.sh [config_file]
# If no config file is specified, the script will use default values.

# Default configuration (can be overridden by config file)
REF_GENOME="../ref/Zm-B73-REFERENCE-NAM-5.0.fa"
BAM_BASE_DIR="/rsstu/users/r/rrellan/DOE_CAREER/BZea/mapped_bwa"
WIDESEQ_REF_DIR="./wideseq_ref"
FILTERED_VCF_DIR="/rsstu/users/r/rrellan/DOE_CAREER/BZea/joint_genotype/all_samps/9_final_samples/more_filtered"
OUTPUT_DIR="./bzea"
TEMP_BAM_DIR="./bam"
SCRIPTS_DIR="./jobs"
CHROMOSOMES=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10")
LSF_QUEUE="normal"
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

# Create output and script directories if they don't exist
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${TEMP_BAM_DIR}"
mkdir -p "${SCRIPTS_DIR}"

# Function to display progress
progress() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# Create sequence dictionary if it doesn't exist
if [ ! -f "${REF_GENOME%.*}.dict" ]; then
    progress "Creating sequence dictionary for reference genome"
    picard CreateSequenceDictionary \
        --REFERENCE="${REF_GENOME}" \
        --OUTPUT="${REF_GENOME%.*}.dict"
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

# Index VCF files for all chromosomes
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
    else
        progress "VCF index already exists for ${CHR}"
    fi
done

# Create a master submission script
SUBMIT_SCRIPT="${SCRIPTS_DIR}/submit_all_jobs.sh"
echo '#!/bin/bash' > "${SUBMIT_SCRIPT}"
echo '' >> "${SUBMIT_SCRIPT}"
echo '# Master script to submit all WideSeq read count jobs to LSF' >> "${SUBMIT_SCRIPT}"
echo '# Generated on '"$(date)" >> "${SUBMIT_SCRIPT}"
echo '' >> "${SUBMIT_SCRIPT}"

# Process each sample plate and find BAM files
TOTAL_SAMPLES=0
for PLATE_DIR in "${SAMPLE_PLATES[@]}"; do
    # Extract plate number from directory name
    PLATE_NUM=$(echo "${PLATE_DIR}" | grep -o "S[0-9]*")
    
    progress "Processing sample plate ${PLATE_NUM} in ${PLATE_DIR}"
    
    # Find all BAM files in this plate directory
    BAM_FILES=($(find "${PLATE_DIR}" -name "*.bam" | grep -v "sorted_rg.bam"))
    
    if [ ${#BAM_FILES[@]} -eq 0 ]; then
        echo "Warning: No BAM files found in ${PLATE_DIR}"
        continue
    fi
    
    progress "Found ${#BAM_FILES[@]} BAM files in plate ${PLATE_NUM}"
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
#BSUB -n ${LSF_CORES}
#BSUB -q ${LSF_QUEUE}
#BSUB -R "rusage[mem=${LSF_MEMORY}]"
#BSUB -o ${SCRIPTS_DIR}/wideseq_${SAMPLE}_%J.out
#BSUB -e ${SCRIPTS_DIR}/wideseq_${SAMPLE}_%J.err

# WideSeq Read Count Collection for ${SAMPLE} (Plate ${PLATE_NUM})
# Generated on $(date)

# Configuration
REF_GENOME="${REF_GENOME}"
BAM_FILE="${BAM_FILE}"
WIDESEQ_REF_DIR="${WIDESEQ_REF_DIR}"
FILTERED_VCF_DIR="${FILTERED_VCF_DIR}"
OUTPUT_DIR="${OUTPUT_DIR}"
TEMP_BAM_DIR="${TEMP_BAM_DIR}"
SAMPLE="${SAMPLE}"
PLATE_NUM="${PLATE_NUM}"
CHROMOSOMES=($(printf '"%s" ' "${CHROMOSOMES[@]}"))

# Function to display progress with timestamp
progress() {
    echo "\$(date '+%Y-%m-%d %H:%M:%S') - \$1"
}

progress "Starting WideSeq read count collection for \${SAMPLE} (Plate \${PLATE_NUM})"

# Step 1: Add read groups if not already present
if ! samtools view -H "\${BAM_FILE}" | grep -q "@RG"; then
    progress "Adding read groups to \${SAMPLE}"
    picard AddOrReplaceReadGroups \\
        -I "\${BAM_FILE}" \\
        -O "\${TEMP_BAM_DIR}/\${SAMPLE}_sorted_rg.bam" \\
        --RGID "\${PLATE_NUM}" \\
        --RGPL illumina \\
        --RGLB "\${SAMPLE}" \\
        --RGPU NONE \\
        --RGSM "\${SAMPLE}"
        
    # Index the read group BAM
    progress "Indexing read group BAM for \${SAMPLE}"
    samtools index "\${TEMP_BAM_DIR}/\${SAMPLE}_sorted_rg.bam"
    
    # Use the new BAM file for downstream processing
    PROCESSED_BAM="\${TEMP_BAM_DIR}/\${SAMPLE}_sorted_rg.bam"
else
    progress "Read groups already present in \${SAMPLE}"
    
    # Check if BAM is indexed
    if [ ! -f "\${BAM_FILE}.bai" ]; then
        progress "Indexing BAM file for \${SAMPLE}"
        samtools index "\${BAM_FILE}"
    fi
    
    # Use the original BAM file for downstream processing
    PROCESSED_BAM="\${BAM_FILE}"
fi

# Find the corresponding VCF file for this sample
SAMPLE_VCF="\${FILTERED_VCF_DIR}/\${SAMPLE}_more_filtered.vcf.gz"

if [ ! -f "\${SAMPLE_VCF}" ]; then
    echo "Error: VCF file \${SAMPLE_VCF} not found for sample \${SAMPLE}. Exiting."
    exit 1
fi

# Check if VCF is indexed
if [ ! -f "\${SAMPLE_VCF}.tbi" ]; then
    progress "Indexing sample VCF file for \${SAMPLE}"
    gatk IndexFeatureFile -I "\${SAMPLE_VCF}"
fi

# Process all chromosomes for this sample
for CHR in "\${CHROMOSOMES[@]}"; do
    CHR_NUM=\${CHR#chr}  # Remove 'chr' prefix to get just the number
    
    progress "Processing chromosome \${CHR} for sample \${SAMPLE}"
    
    # Step 2: Select variants using teosinte reference as intervals
    progress "Selecting variants for \${SAMPLE} on \${CHR}"
    gatk SelectVariants \\
        -R "\${REF_GENOME}" \\
        -V "\${SAMPLE_VCF}" \\
        -L "\${WIDESEQ_REF_DIR}/wideseq_chr\${CHR_NUM}.vcf.gz" \\
        -O "\${OUTPUT_DIR}/\${SAMPLE}_wideseq_selected_chr\${CHR_NUM}.vcf"
    
    # Step 3: Collect allelic counts
    progress "Collecting allelic counts for \${SAMPLE} on \${CHR}"
    gatk CollectAllelicCounts \\
        -I "\${PROCESSED_BAM}" \\
        -R "\${REF_GENOME}" \\
        -L "\${OUTPUT_DIR}/\${SAMPLE}_wideseq_selected_chr\${CHR_NUM}.vcf" \\
        -O "\${OUTPUT_DIR}/\${SAMPLE}_chr\${CHR_NUM}.allelicCounts.tsv"
    
    progress "Completed allelic count collection for \${SAMPLE} on \${CHR}"
done

# Merge results from all chromosomes for this sample
progress "Merging allelic counts for all chromosomes of \${SAMPLE}"

# Get the header from the first chromosome file
head -n 1 "\${OUTPUT_DIR}/\${SAMPLE}_chr1.allelicCounts.tsv" > "\${OUTPUT_DIR}/\${SAMPLE}.allelicCounts.tsv"

# Append data from all chromosome files (skipping headers)
for CHR in "\${CHROMOSOMES[@]}"; do
    CHR_NUM=\${CHR#chr}
    if [ -f "\${OUTPUT_DIR}/\${SAMPLE}_chr\${CHR_NUM}.allelicCounts.tsv" ]; then
        tail -n +2 "\${OUTPUT_DIR}/\${SAMPLE}_chr\${CHR_NUM}.allelicCounts.tsv" >> "\${OUTPUT_DIR}/\${SAMPLE}.allelicCounts.tsv"
    fi
done

# Clean up temporary files
progress "Cleaning up temporary files"
for CHR in "\${CHROMOSOMES[@]}"; do
    CHR_NUM=\${CHR#chr}
    rm -f "\${OUTPUT_DIR}/\${SAMPLE}_wideseq_selected_chr\${CHR_NUM}.vcf"
done

progress "WideSeq read count collection completed for \${SAMPLE} (Plate \${PLATE_NUM})"

# Calculate some statistics
VARIANT_POSITIONS=\$(wc -l "\${OUTPUT_DIR}/\${SAMPLE}.allelicCounts.tsv" | awk '{print \$1-1}')  # Subtract 1 for header

echo "=== WideSeq Read Count Collection Summary for \${SAMPLE} ==="
echo "Sample name: \${SAMPLE}"
echo "Plate: \${PLATE_NUM}"
echo "Total positions collected: \${VARIANT_POSITIONS}"
echo "Results available in: \${OUTPUT_DIR}/\${SAMPLE}.allelicCounts.tsv"
echo "=========================================================="
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
echo "=== WideSeq Job Generation Summary ==="
echo "Total sample plates found: ${#SAMPLE_PLATES[@]}"
echo "Total samples: ${TOTAL_SAMPLES}"
echo "Job scripts created in: ${SCRIPTS_DIR}"
echo "To submit all jobs, run: ${SUBMIT_SCRIPT}"
echo "===================================="
