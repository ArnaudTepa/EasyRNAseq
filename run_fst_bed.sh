#!/bin/bash

# Check arguments
if [ "$#" -lt 8 ] || [ "$#" -gt 9 ]; then
    echo "Usage: $0 perl_script input.sync output_dir pool_size window_sizes step_sizes MAX_JOBS sample_file [bed_file]"
    echo "Example: $0 fst-sliding.pl input.sync outdir 40 100000,250000 1000 4 sample.txt regions.bed"
    exit 1
fi

perl_file=$1
input_file=$2
output_dir=$3
pool_size=$4
window_sizes_raw=$5
step_sizes_raw=$6
MAX_JOBS=$7
sample_file=$8
bed_file=${9:-}

# Validate sample file
if [ ! -f "$sample_file" ]; then
    echo "? Error: Sample file '$sample_file' not found."
    exit 1
fi

# Validate BED file if provided
if [ -n "$bed_file" ] && [ ! -f "$bed_file" ]; then
    echo "? Error: BED file '$bed_file' not found."
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Extract base filename
filename=$(basename "$input_file" .sync)

# Create a temporary file for cleaned input
TMP_DATA="${filename}_cleaned_chrom.sync"

# Remove prefix before underscore from chromosome names
awk 'BEGIN {FS=OFS="	"}
     !/^#/ {
         split($1, parts, "_");
         $1 = parts[2];
         print
     }' "$input_file" > "$TMP_DATA"

# If BED file provided, filter TMP_DATA by regions
if [ -n "$bed_file" ]; then
    echo "? Filtering sync file by BED regions..."
    FILTERED_DATA="${filename}_filtered.sync"
    awk 'NR==FNR {start[$1][++count[$1]]=$2; end[$1][count[$1]]=$3; next}
    {
      chrom=$1; pos=$2;
      if (chrom in start) {
        for (i=1; i<=count[chrom]; i++) {
          if (pos >= start[chrom][i] && pos <= end[chrom][i]) {
            print $0;
            break;
          }
        }
      }
    }' "$bed_file" "$TMP_DATA" > "$FILTERED_DATA"
    TMP_DATA="$FILTERED_DATA"
fi

# Reheader function
function reheader_fst_file {
    local input_file="$1"
    local output_file="$2"

    # Read sample names
    mapfile -t samples < "$sample_file"

    # Generate pairwise labels
    pairwise_labels=()
    for ((i=0;i<${#samples[@]};i++)); do
        for ((j=i+1;j<${#samples[@]};j++)); do
            pairwise_labels+=("${samples[i]}:${samples[j]}")
        done
    done

    # Write header and cleaned data
    {
        echo -e "CHROM	POS	N_VAR	FST	COV	$(IFS=$'	'; echo "${pairwise_labels[*]}")"
        awk '{
            printf "%s	%s	%s	%s	%s", $1, $2, $3, $4, $5;
            for (i=6;i<=NF;i++) {
                split($i, val, "=");
                printf "	%s", val[2];
            }
            printf "
";
        }' "$input_file"
    } > "$output_file"
}

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

    output_file="${output_dir}/${filename}_fst${window_size}_step${step_size}.txt"
    cleaned_file="${output_dir}/${filename}_fst${window_size}_step${step_size}_cleaned.txt"
    echo "? Running: window=$window_size, step=$step_size"

    (
      perl "$perl_file"         --input "$TMP_DATA"         --output "$output_file"         --min-count 2         --min-coverage 5         --max-coverage 1000         --window-size "$window_size"         --step-size "$step_size"         --pool-size "$pool_size"

      # Reheader and clean
      if [ -s "$output_file" ]; then
        reheader_fst_file "$output_file" "$cleaned_file"
        rm -f "$output_file"
        echo "? Cleaned file saved: $cleaned_file"
      else
        echo "?? Warning: $output_file is empty or missing. Skipping reheader."
      fi
    ) &
  done
done

wait
rm -f "$TMP_DATA"
echo "? All FST calculations and reheadering completed."
