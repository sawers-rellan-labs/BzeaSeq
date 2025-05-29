#!/bin/bash
#
# get_all_sample_ancestry.sh
#
# Description: This script processes all merged allelic count files and runs
# BinAncestry.R on each one to generate bin-level ancestry calls.
#
# Usage: ./process_all_samples.sh [bin_size]
# If no bin_size is specified, 1000000 (1Mb) will be used.

# Configuration
ALLELIC_COUNTS_DIR="./allelic_counts"
ANCESTRY_DIR="./ancestry"
ANCESTRY_SCRIPT="./get_ancestry_calls.R"
BIN_SIZE=${1:-1000000}  # Default bin size is 1Mb unless specified as first argument
LOG_FILE="get_all_sample_ancestry.log"

# Function to log messages with timestamp
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "${LOG_FILE}"
}

# Check if BinAncestry.R script exists
if [ ! -f "${ANCESTRY_SCRIPT}" ]; then
    echo "Error: BinAncestry.R script not found at ${ANCESTRY_SCRIPT}. Aborting."
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "${ANCESTRY_DIR}"

# Initialize log file
echo "Starting processing at $(date)" > "${LOG_FILE}"
log_message "Using bin size of ${BIN_SIZE} bp"
log_message "Finding merged allelic count files..."

# Find all merged allelic count files (excluding chromosome-specific files)
ALLELIC_FILES=$(ls -1 ${ALLELIC_COUNTS_DIR}/*.allelicCounts.tsv | grep -v "chr")

# Count files
FILE_COUNT=$(echo "${ALLELIC_FILES}" | wc -l)
log_message "Found ${FILE_COUNT} merged allelic count files to process"

# Initialize counters
PROCESSED=0
SUCCEEDED=0
FAILED=0

# Process each file
for FILE in ${ALLELIC_FILES}; do
    # Extract sample name from file name
    SAMPLE=$(basename "${FILE}" .allelicCounts.tsv)

    # Update progress
    PROCESSED=$((PROCESSED + 1))
    log_message "[${PROCESSED}/${FILE_COUNT}] Processing sample ${SAMPLE}"

    # Run BinAncestry.R
    log_message "Running: Rscript ${ANCESTRY_SCRIPT} ${FILE} ${ANCESTRY_DIR}/${SAMPLE} ${BIN_SIZE}"

    Rscript "${ANCESTRY_SCRIPT}" "${FILE}" "${ANCESTRY_DIR}/${SAMPLE}" "${BIN_SIZE}"

    # Check if the script succeeded
    if [ $? -eq 0 ]; then
        SUCCEEDED=$((SUCCEEDED + 1))
        log_message "Successfully processed ${SAMPLE}"
    else
        FAILED=$((FAILED + 1))
        log_message "ERROR: Failed to process ${SAMPLE}"
    fi

    log_message "Progress: ${PROCESSED}/${FILE_COUNT} (${SUCCEEDED} succeeded, ${FAILED} failed)"
    echo "--------------------------------------------------------------------------------" >> "${LOG_FILE}"
done

# Final summary
log_message "Processing completed"
log_message "Total files: ${FILE_COUNT}"
log_message "Succeeded: ${SUCCEEDED}"
log_message "Failed: ${FAILED}"

# Print instructions for viewing results
echo ""
echo "Processing completed. See ${LOG_FILE} for details."
echo "Results are available in ${ANCESTRY_DIR}/"
echo ""

if [ ${FAILED} -gt 0 ]; then
    echo "WARNING: ${FAILED} files failed to process. Check the log for details."
    exit 1
else
    echo "All files processed successfully."
    exit 0
fi
