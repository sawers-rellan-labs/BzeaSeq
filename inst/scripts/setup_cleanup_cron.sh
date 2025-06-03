#!/bin/bash
#
# setup_cleanup_cron.sh
#
# Description: Sets up a cron job to run the cleanup script every 20 minutes
#
# Usage: ./setup_cleanup_cron.sh

# Configuration
CLEANUP_SCRIPT="./run_cleanup_once.sh"
SCRIPT_DIR="$(pwd)"

# Make sure the cleanup script is executable
chmod +x "${CLEANUP_SCRIPT}"

# Create a temporary file for the new crontab
TEMP_CRONTAB=$(mktemp)

# Export current crontab
crontab -l > "${TEMP_CRONTAB}" 2>/dev/null

# Check if the cleanup job is already in crontab
if ! grep -q "${CLEANUP_SCRIPT}" "${TEMP_CRONTAB}"; then
    # Add the cleanup job to run every 30 minutes
    echo "*/30 * * * * cd ${SCRIPT_DIR} && ${SCRIPT_DIR}/${CLEANUP_SCRIPT} >> ${SCRIPT_DIR}/cron_cleanup.log 2>&1" >> "${TEMP_CRONTAB}"

    # Install the new crontab
    crontab "${TEMP_CRONTAB}"

    echo "Cron job set up successfully. The cleanup script will run every 30 minutes."
else
    echo "Cron job already exists for the cleanup script."
fi

# Remove the temporary file
rm "${TEMP_CRONTAB}"

# Display current crontab to verify
echo "Current crontab:"
crontab -l
