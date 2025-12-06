#!/bin/bash

# Set variables
INPUT_SYNC="$1"
SCRIPT_PATH="$2"
OUTPUT_DIR="$3"
SAMPLE_FILE="$4"

# Extract base name without extension
BASENAME=$(basename "$INPUT_SYNC" .sync)
FINAL_OUTPUT="${OUTPUT_DIR}/${BASENAME}_fisher_results.txt"
CLEANED_FILE="${OUTPUT_DIR}/${BASENAME}_fisher_results_cleaned.txt"
TMP_DATA="${OUTPUT_DIR}/${BASENAME}_cleaned_chrom.sync"

# Step 0: Check input file
if [ ! -f "$INPUT_SYNC" ]; then
  echo "? Error: Input file '$INPUT_SYNC' not found."
  exit 1
fi

# Step 1: Prepare output directory
mkdir -p "$OUTPUT_DIR"

# Step 2: Validate sample file
if [ ! -f "$SAMPLE_FILE" ]; then
  echo "? Error: Sample file '$SAMPLE_FILE' not found."
  exit 1
fi

# Step 3: Clean chromosome names and create cleaned input
awk 'BEGIN {FS=OFS="\t"}
     !/^#/ {
         split($1, parts, "_");
         $1 = parts[2];
         print
     }' "$INPUT_SYNC" > "$TMP_DATA"

# Step 4: Run fisher-test.pl
echo "?? Running fisher-test.pl on cleaned sync file..."
perl "$SCRIPT_PATH" \
  --input "$TMP_DATA" \
  --output "$FINAL_OUTPUT" \
  --min-count 2 \
  --min-coverage 5 \
  --max-coverage 1000  \
  --suppress-noninformative

# Step 5: Reheader function
reheader_fst_file() {
    local input_file="$1"
    local output_file="$2"

    # Read sample names
    mapfile -t samples < "$SAMPLE_FILE"

    # Generate pairwise labels
    pairwise_labels=()
    for ((i=0; i<${#samples[@]}; i++)); do
        for ((j=i+1; j<${#samples[@]}; j++)); do
            pairwise_labels+=("${samples[i]}:${samples[j]}")
        done
    done

    # Write header and cleaned data
    {
        echo -e "CHROM\tPOS\tN_VAR\tfrac_cov\tmin_cov\t$(IFS=$'\t'; echo "${pairwise_labels[*]}")"
        awk '{
            printf "%s\t%s\t%s\t%s\t%s", $1, $2, $3, $4, $5;
            for (i=6; i<=NF; i++) {
                split($i, val, "=");
                printf "\t%s", val[2];
            }
            printf "\n";
        }' "$input_file"
    } > "$output_file"
}

# Step 6: Reheader and clean output
if [ -s "$FINAL_OUTPUT" ]; then
  echo "?? Reheadering and cleaning output..."
  reheader_fst_file "$FINAL_OUTPUT" "$CLEANED_FILE"
  rm -f "$FINAL_OUTPUT"
  echo "? Cleaned file saved: $CLEANED_FILE"
else
  echo "?? Warning: '$FINAL_OUTPUT' is empty or missing. Skipping reheader."
fi
