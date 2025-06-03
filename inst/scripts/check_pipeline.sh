#!/bin/bash
#check_pipeline.sh

# Job status monitoring script
OUTPUT_DIR="/path/to/output"
LOG_DIR="${OUTPUT_DIR}/logs"

# Check status of all jobs
function check_status() {
  echo "=== Current Job Status ==="
  bjobs | awk '{print $1, $2, $3}' | column -t
  
  echo -e "\n=== Job Distribution ==="
  bjobs | awk '{print $3}' | sort | uniq -c | sort -nr
  
  echo -e "\n=== Recently Completed Jobs ==="
  bjobs -d -w | tail -n 10
}

# Kill all jobs in a specific state
function kill_jobs() {
  state=$1
  if [ -z "$state" ]; then
  echo "Please specify a job state (e.g., PEND, RUN, ZOMBI)"
  return 1
  fi
  
  echo "Killing all jobs in $state state..."
  jobids=$(bjobs | grep $state | awk '{print $1}')
  if [ -n "$jobids" ]; then
  bkill $jobids
  echo "Killed jobs: $jobids"
  else
    echo "No jobs in $state state found."
  fi
}

# Restart failed jobs
function restart_failed() {
  echo "Looking for failed jobs..."
  find $LOG_DIR -name "*.err" -type f -size +0c | while read errfile; do
  jobname=$(basename $errfile .err | cut -d '_' -f 1)
  echo "Potential failure in job: $jobname"
  echo "Error file: $errfile"
  echo "First few lines of error:"
  head -n 5 $errfile
  echo "-----"
  done
  
  echo -e "\nTo resubmit a specific job array, use: bsub < job_script.sh"
}

# Show resource usage
function show_resources() {
  echo "=== Resource Usage ==="
  lsload | sort -k 7 -nr | head -n 10
  
  echo -e "\n=== Jobs with High Memory Usage ==="
  bjobs -l | grep -A 5 "MEMORY USAGE" | grep -v "^$" | head -n 20
}

# Main menu
case "$1" in
status)
check_status
;;
kill)
kill_jobs $2
;;
failed)
restart_failed
;;
resources)
show_resources
;;
*)
echo "Usage: $0 {status|kill [state]|failed|resources}"
;;
esac