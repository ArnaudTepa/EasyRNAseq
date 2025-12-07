
#!/bin/bash

# Usage:
# ./fst_grenedalf.sh "bam1 bam2 ..." output_dir pool_size threads prefix grenedalf_path ref_genome window_queue_count file_suffix
# Example:
# ./fst_grenedalf.sh "/path/to/sample1.bam /path/to/sample2.bam" /path/to/output 40 40 MyPrefix /home/arnaud.tepa/pipeline/grenedalf/bin/grenedalf /home/arnaud.tepa/reference/AGAP61.fasta 100000 .100ksnp.fst

if [ "$#" -ne 9 ]; then
    echo "Usage: $0 \"bam_files\" output_dir pool_size threads prefix grenedalf_path ref_genome window_queue_count file_suffix"
    exit 1
fi

BAM_FILES=$1
OUT_DIR=$2
POOL_SIZE=$3
THREADS=$4
PREFIX=$5
GRENDALF_BIN=$6
REF_GENOME=$7
WINDOW_QUEUE_COUNT=$8
FILE_SUFFIX=$9

# Create output directory if not exists
mkdir -p "$OUT_DIR"

# Log file
LOG_FILE="${OUT_DIR}/${PREFIX}${FILE_SUFFIX}.log"

echo "? Starting Grenedalf FST calculation..."
echo "Input BAM files: $BAM_FILES"
echo "Output directory: $OUT_DIR"
echo "Pool size: $POOL_SIZE"
echo "Threads: $THREADS"
echo "Prefix: $PREFIX"
echo "Grenedalf path: $GRENDALF_BIN"
echo "Reference genome: $REF_GENOME"
echo "Window queue count: $WINDOW_QUEUE_COUNT"
echo "File suffix: $FILE_SUFFIX"

$GRENDALF_BIN fst \
  --sam-path $BAM_FILES \
  --sam-min-map-qual 10 \
  --sam-min-base-qual 10 \
  --multi-file-locus-set intersection \
  --reference-genome-fasta "$REF_GENOME" \
  --filter-total-min-read-depth 10 \
  --filter-total-only-biallelic-snps \
  --filter-total-snp-min-frequency 0.05 \
  --window-type queue \
  --window-queue-count "$WINDOW_QUEUE_COUNT" \
  --window-queue-stride 0 \
  --window-average-policy valid-snps \
  --method kofler \
  --pool-sizes "$POOL_SIZE" \
  --separator-char tab \
  --out-dir "$OUT_DIR" \
  --file-prefix "$PREFIX" \
  --file-suffix "$FILE_SUFFIX" \
  --allow-file-overwriting \
  --threads "$THREADS" \
  --log-file "$LOG_FILE"

echo "? FST calculation completed!"
echo "Results saved in: $OUT_DIR"
echo "Log file: $LOG_FILE"
