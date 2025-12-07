#!/bin/bash

# Usage:
# ./diversity_grenedalf_sync.sh sync_file output_dir pool_size threads prefix grenedalf_path bed_file file_suffix sample_file
# Example:
# ./diversity_grenedalf_sync.sh /path/to/input.sync /path/to/output 40 40 MyPrefix /home/arnaud.tepa/pipeline/grenedalf/bin/grenedalf /path/to/intervals.bed .snp.diversity sample.txt

if [ "$#" -ne 9 ]; then
    echo "Usage: $0 sync_file output_dir pool_size threads prefix grenedalf_path bed_file file_suffix sample_file"
    exit 1
fi

SYNC_FILE=$1
OUT_DIR=$2
POOL_SIZE=$3
THREADS=$4
PREFIX=$5
GRENDALF_BIN=$6
BED_FILE=$7
FILE_SUFFIX=$8
SAMPLE_FILE=$9

# Validate sample file
if [ ! -f "$SAMPLE_FILE" ]; then
    echo "? Error: Sample file '$SAMPLE_FILE' not found."
    exit 1
fi

# Create output directory if not exists
mkdir -p "$OUT_DIR"

# Log file
LOG_FILE="${OUT_DIR}/${PREFIX}${FILE_SUFFIX}.log"

echo "? Starting Grenedalf Diversity calculation..."
echo "Input SYNC file: $SYNC_FILE"
echo "Output directory: $OUT_DIR"
echo "Pool size: $POOL_SIZE"
echo "Threads: $THREADS"
echo "Prefix: $PREFIX"
echo "Grenedalf path: $GRENDALF_BIN"
echo "BED file: $BED_FILE"
echo "File suffix: $FILE_SUFFIX"
echo "Sample file: $SAMPLE_FILE"

# Run Grenedalf Diversity
$GRENDALF_BIN diversity \
  --sync-path "$SYNC_FILE" \
  --multi-file-locus-set intersection \
  --filter-sample-min-count 2 \
  --filter-sample-min-read-depth 10 \
  --filter-total-min-read-depth 10 \
  --filter-total-snp-min-frequency 0.01 \
  --window-type regions \
  --window-region-bed "$BED_FILE" \
  --window-average-policy window-length \
  --pool-sizes "$POOL_SIZE" \
  --separator-char tab \
  --out-dir "$OUT_DIR" \
  --file-prefix "$PREFIX" \
  --file-suffix "$FILE_SUFFIX" \
  --allow-file-overwriting \
  --threads "$THREADS" \
  --log-file "$LOG_FILE"

# ? Auto-detect Grenedalf output file (csv, tsv, or txt)
for ext in csv tsv txt; do
    if [ -f "${OUT_DIR}/${PREFIX}diversity${FILE_SUFFIX}.${ext}" ]; then
        OUTPUT_FILE="${OUT_DIR}/${PREFIX}diversity${FILE_SUFFIX}.${ext}"
        break
    fi
done

if [ -z "$OUTPUT_FILE" ]; then
    echo "? Error: No Grenedalf output file found with prefix ${PREFIX}diversity${FILE_SUFFIX} in ${OUT_DIR}"
    exit 1
fi


# ? Clean and reheader function
function clean_diversity_file {
    local input_file="$1"
    local output_file="$2"

    # Read sample names into array
    mapfile -t samples < "$SAMPLE_FILE"

    # Extract basename of sync file (without extension)
    base=$(basename "$SYNC_FILE")
    base=${base%%.*}

    # Extract header
    header=$(head -1 "$input_file")

    # Keep only chrom, start, end and diversity columns
    keep_cols="chrom|start|end|theta_pi|theta_watterson|tajimas_d"
    new_header=$(echo "$header" | tr '\t' '\n' | grep -E "$keep_cols" | tr '\n' '\t')

    # Replace basename.N with sample names
    for i in "${!samples[@]}"; do
        n=$((i+1))
        new_header=$(echo "$new_header" | sed "s/${base}.${n}/${samples[$i]}/g")
    done

    # Write new header
    echo -e "${new_header%$'\t'}" > "$output_file"

    # Filter data rows
    awk -v keep="$keep_cols" '
        NR==1 {
            for (i=1; i<=NF; i++) {
                if ($i ~ keep) col[i]=1
            }
            next
        }
        {
            out=""
            for (i=1; i<=NF; i++) {
                if (col[i]) out = out ? out OFS $i : $i
            }
            print out
        }
    ' OFS="\t" "$input_file" >> "$output_file"
}

# ? Apply cleaning if output exists
if [ -s "$OUTPUT_FILE" ]; then
    FINAL_FILE="${OUT_DIR}/${PREFIX}diversity${FILE_SUFFIX}_cleaned.txt"
    clean_diversity_file "$OUTPUT_FILE" "$FINAL_FILE"
    echo "? Cleaned file saved: $FINAL_FILE"
else
    echo "?? Warning: $OUTPUT_FILE is empty. Skipping cleaning."
fi

echo "? Diversity calculation and cleaning completed!"
echo "Results saved in: $OUT_DIR"
echo "Log file: $LOG_FILE"