#!/bin/bash
set -e

# ---------------------------
# User-defined variables
# ---------------------------
mpileup2sync_PATH=$1     # Path to mpileup2sync jar
REF=$2                   # Reference FASTA
FAI_FILE=$3              # Reference .fai file
BAM_DIR=$4               # Directory containing BAM files
WORK_DIR=$5              # Output directory

REGIONS_FILE="$WORK_DIR/regions.txt"

# ---------------------------
# Create output directory
# ---------------------------
mkdir -p "$WORK_DIR"


# ---------------------------
# Generate regions from .fai (if not exists)
# ---------------------------
if [ ! -f "$REGIONS_FILE" ]; then
    echo "Generating regions from $FAI_FILE..."
    cut -f 1,2 "$FAI_FILE" | awk '
    {
        chrom=$1; size=$2; region_size=1000000000;
        for (start=1; start<=size; start+=region_size) {
            end=start+region_size-1;
            if (end>size) end=size;
            print chrom":"start"-"end;
        }
    }' > "$REGIONS_FILE"
fi

# ---------------------------
# Find all BAM files
# ---------------------------

ls "$BAM_DIR"/*.bam > bam_files.txt

# ---------------------------
# Run pileup per-region in parallel
# ---------------------------
echo "Starting parallel variant calling..."
cat "$REGIONS_FILE" | while read region; do
  (
    # Loop each bam
    for bam in "$BAM_DIR"/*.bam; do
        sample=$(basename "$bam" .bam)

    # windowed pileup file
    out_pileup="$WORK_DIR/${sample}_${region//[:\-]/_}.pileup"

    # Run pileup

        echo "Processing $sample for region $region..."
        samtools mpileup -Q 30 -f "$REF" "$bam" -r "$region" > "$out_pileup"
    done
  ) &
done

wait  # Ensure all parallel jobs are finished

echo "Merging pileup files by sample..."
for bam in "$BAM_DIR"/*.bam; do
    sample=$(basename "$bam" .bam)
    final_pileup="$WORK_DIR/${sample}.pileup"

    # Merge all region-specific pileups for this sample
    cat "$WORK_DIR/${sample}_*.pileup" > "$final_pileup"

    # Remove intermediate files
    rm "$WORK_DIR/${sample}_*.pileup"
done

echo "Final merged pileup files saved in: $WORK_DIR"
