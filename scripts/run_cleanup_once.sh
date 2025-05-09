#!/bin/bash
#
# run_cleanup_once.sh
#
# Description: One-time script to clean up chromosome-specific allelic count files
# for samples where the merged file is greater than 590 MB.
#
# Usage: ./run_cleanup_once.sh
#        (Run manually or as a cron job)

# Configuration
ALLELIC_COUNTS_DIR="/rsstu/users/r/rrellan/BZea/bzeaseq/allelic_counts"
LOG_FILE="./cleanup_allelic_counts.log"
SIZE_THRESHOLD=590   # Size threshold in MB

# Function to log messages with timestamp
log_message() {
  echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "${LOG_FILE}"
}

# Function to clean up chromosome files for a sample
cleanup_chromosome_files() {
  local sample=$1
  local base_name=$(basename "${sample}" .allelicCounts.tsv)
  local deleted_count=0
  local deleted_size=0
  
  # Find all chromosome-specific files for this sample
  for chr_file in "${ALLELIC_COUNTS_DIR}/${base_name}_chr"*".allelicCounts.tsv"; do
  if [ -f "${chr_file}" ]; then
  # Get file size before deletion
  local file_size=$(du -k "${chr_file}" | cut -f1)
  
  # Delete the file
  rm -f "${chr_file}"
  
  if [ $? -eq 0 ]; then
  deleted_count=$((deleted_count + 1))
  deleted_size=$((deleted_size + file_size))
  fi
  fi
  done
  
  # Convert KB to MB for reporting
  deleted_size_mb=$(echo "scale=2; ${deleted_size}/1024" | bc)
  
  log_message "Cleaned up ${deleted_count} chromosome files for ${base_name}, freed ${deleted_size_mb} MB"
}

# Create log file directory if it doesn't exist
mkdir -p "$(dirname "${LOG_FILE}")"

# Initial log entry
log_message "Starting one-time allelic counts cleanup"
log_message "Monitoring directory: ${ALLELIC_COUNTS_DIR}"
log_message "Size threshold: ${SIZE_THRESHOLD} MB"

# Find all merged allelic count files (not chromosome-specific)
merged_files=$(find "${ALLELIC_COUNTS_DIR}" -name "*.allelicCounts.tsv" -not -name "*_chr*.allelicCounts.tsv")

total_files=0
processed_files=0

for file in ${merged_files}; do
total_files=$((total_files + 1))

# Get file size in MB
file_size_mb=$(du -m "${file}" | cut -f1)

# Check if file size exceeds threshold
if [ ${file_size_mb} -gt ${SIZE_THRESHOLD} ]; then
log_message "Processing: $(basename "${file}") (${file_size_mb} MB)"

# Clean up chromosome files for this sample
cleanup_chromosome_files "${file}"

processed_files=$((processed_files + 1))
fi
done

log_message "Cleanup completed: ${processed_files} of ${total_files} merged files processed"
