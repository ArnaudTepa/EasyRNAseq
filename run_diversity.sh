#!/bin/bash

# Check arguments
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 /path/to/Variance-sliding.pl input.pileup output_dir pool_size window_sizes step_sizes MAX_JOBS"
    echo "Example: $0 /path/to/Variance-sliding.pl input.pileup /path/to/output 40 5000,10000 1000,2000 4"
    exit 1
fi

Variance_tool=$1
input_file=$2
output_dir=$3
pool_size=$4
window_sizes_raw=$5
step_sizes_raw=$6
MAX_JOBS=$7

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Extract base filename
filename=$(basename "$input_file" .pileup)

# Function to limit parallel jobs
function wait_for_jobs {
    while [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; do
        sleep 1
    done
}

# Convert comma-separated strings to arrays
IFS=',' read -r -a WINDOW_SIZES <<< "$window_sizes_raw"
IFS=',' read -r -a STEP_SIZES <<< "$step_sizes_raw"

# Loop over combinations
for window_size in "${WINDOW_SIZES[@]}"; do
  for step_size in "${STEP_SIZES[@]}"; do
    wait_for_jobs

    tajima_out="${output_dir}/${filename}_tajima_${window_size}_step${step_size}.txt"
    pi_out="${output_dir}/${filename}_pi_${window_size}_step${step_size}.txt"

    echo "?? Running Tajima's D: window=$window_size, step=$step_size"
    perl "$Variance_tool" \
      --input "$input_file" \
      --output "$tajima_out" \
      --measure d \
      --min-count 2 \
      --min-coverage 5 \
      --max-coverage 1000 \
      --window-size "$window_size" \
      --step-size "$step_size" \
      --pool-size "$pool_size" \
      --fastq-type illumina &

    echo "?? Running Pi: window=$window_size, step=$step_size"
    perl "$Variance_tool" \
      --input "$input_file" \
      --output "$pi_out" \
      --measure pi \
      --min-count 2 \
      --min-coverage 5 \
      --max-coverage 1000 \
      --window-size "$window_size" \
      --step-size "$step_size" \
      --pool-size "$pool_size" \
      --fastq-type illumina &
  done
done

wait
echo "? All diversity calculations completed."

