#!/bin/bash
#
# get_ancestry_parallel.sh
#
# Description: This script processes all merged allelic count files and runs
# get_ancestry_calls.R on each one in parallel using GNU parallel.
#
# Usage: ./get_ancestry_parallel.sh [bin_size] [max_jobs]
# - bin_size: Optional, default is 1000000 (1Mb)
# - max_jobs: Optional, maximum number of parallel jobs, default is number of cores
# chmod +x parallel_process_samples.sh
# ./get_ancestry_parallel.sh 1000000 8  # Process with bin size 1Mb, max 8 parallel jobs


# Configuration
ALLELIC_COUNTS_DIR="./allelic_counts"
ANCESTRY_DIR="./ancestry"
ANCESTRY_SCRIPT="./get_ancestry_calls.R"
BIN_SIZE=${1:-1000000}  # Default bin size is 1Mb unless specified as first argument
MAX_JOBS=${2:-0}  # Default to number of available cores (0 means use all available)
LOG_DIR="./logs"
MAIN_LOG="${LOG_DIR}/get_ancestry_parallel.log"

# Function to log messages with timestamp
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "${MAIN_LOG}"
}

# Check if GNU parallel is installed
if ! command -v parallel &>/dev/null; then
    echo "Error: GNU parallel is not installed. Please install it first."
    echo "On most systems, you can install it with: apt-get install parallel"
    exit 1
fi

# Check if BinAncestry.R script exists
if [ ! -f "${ANCESTRY_SCRIPT}" ]; then
    echo "Error: get_ancestry_calls.R script not found at ${ANCESTRY_SCRIPT}. Aborting."
    exit 1
fi

# Create output directories
mkdir -p "${ANCESTRY_DIR}"
mkdir -p "${LOG_DIR}"

# Initialize log file
echo "Starting parallel processing at $(date)" > "${MAIN_LOG}"
log_message "Using bin size of ${BIN_SIZE} bp"
log_message "Finding merged allelic count files..."

# Find all merged allelic count files (excluding chromosome-specific files)
ALLELIC_FILES=$(ls -1 ${ALLELIC_COUNTS_DIR}/*.allelicCounts.tsv | grep -v "chr")

# Create a temporary file with the list of samples to process
TEMP_FILE=$(mktemp)
echo "${ALLELIC_FILES}" > "${TEMP_FILE}"

# Count files
FILE_COUNT=$(echo "${ALLELIC_FILES}" | wc -l)
log_message "Found ${FILE_COUNT} merged allelic count files to process"

# Function to process a single sample
process_sample() {
    local FILE=$1
    local SAMPLE=$(basename "${FILE}" .allelicCounts.tsv)
    local LOG_FILE="${LOG_DIR}/${SAMPLE}.log"

    echo "$(date '+%Y-%m-%d %H:%M:%S') - Processing sample ${SAMPLE}" > "${LOG_FILE}"

    # Run BinAncestry.R
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Running: Rscript ${ANCESTRY_SCRIPT} ${FILE} ${ANCESTRY_DIR}/${SAMPLE} ${BIN_SIZE}" >> "${LOG_FILE}"

    Rscript "${ANCESTRY_SCRIPT}" "${FILE}" "${ANCESTRY_DIR}/${SAMPLE}" "${BIN_SIZE}" >> "${LOG_FILE}" 2>&1

    # Check if the script succeeded
    if [ $? -eq 0 ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') - Successfully processed ${SAMPLE}" >> "${LOG_FILE}"
        echo "SUCCESS: ${SAMPLE}"
    else
        echo "$(date '+%Y-%m-%d %H:%M:%S') - ERROR: Failed to process ${SAMPLE}" >> "${LOG_FILE}"
        echo "FAILED: ${SAMPLE}"
    fi
}

# Export the function and variables so parallel can use them
export -f process_sample
export ANCESTRY_SCRIPT
export ANCESTRY_DIR
export BIN_SIZE
export LOG_DIR

# Determine number of CPU cores if MAX_JOBS not specified
if [ ${MAX_JOBS} -eq 0 ]; then
    if command -v nproc &>/dev/null; then
        MAX_JOBS=$(nproc)
    elif [ -f /proc/cpuinfo ]; then
        MAX_JOBS=$(grep -c ^processor /proc/cpuinfo)
    else
        # Default to 4 if we can't determine
        MAX_JOBS=4
    fi
fi

log_message "Starting parallel processing with up to ${MAX_JOBS} simultaneous jobs"

# Run parallel processing and capture results
RESULTS=$(cat "${TEMP_FILE}" | parallel --jobs ${MAX_JOBS} --eta process_sample)

# Count successes and failures
SUCCEEDED=$(echo "${RESULTS}" | grep "SUCCESS" | wc -l)
FAILED=$(echo "${RESULTS}" | grep "FAILED" | wc -l)

# Final summary
log_message "Parallel processing completed"
log_message "Total files: ${FILE_COUNT}"
log_message "Succeeded: ${SUCCEEDED}"
log_message "Failed: ${FAILED}"

# Print list of failed samples if any
if [ ${FAILED} -gt 0 ]; then
    log_message "Failed samples:"
    echo "${RESULTS}" | grep "FAILED" | cut -d' ' -f2 | tee -a "${MAIN_LOG}"
fi

# Clean up
rm "${TEMP_FILE}"

# Print instructions for viewing results
echo ""
echo "Processing completed. See ${MAIN_LOG} for details."
echo "Individual sample logs are available in ${LOG_DIR}/"
echo "Results are available in ${ANCESTRY_DIR}/"
echo ""

if [ ${FAILED} -gt 0 ]; then
    echo "WARNING: ${FAILED} files failed to process. Check the logs for details."
    exit 1
else
    echo "All files processed successfully."
    exit 0
fi
