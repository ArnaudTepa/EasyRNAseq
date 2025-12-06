#!/bin/bash
set -e

# ---------------------------
# Usage: ./run_varscan.sh reference.fasta bam_dir output.vcf varscan.jar
# Example: ./run_varscan.sh ref.fa /path/to/bams merged.vcf /path/to/VarScan.jar
# ---------------------------

REF=$1             # Reference FASTA
BAM_DIR=$2         # Directory containing BAM files
OUTPUT_VCF=$3      # Final VCF output
VARSCAN_PATH=$4    # Path to VarScan jar

# Validate inputs
if [ $# -ne 4 ]; then
    echo "Usage: $0 reference.fasta bam_dir output.vcf varscan.jar"
    exit 1
fi

if [ ! -f "$REF" ]; then
    echo "Error: Reference FASTA not found!"
    exit 1
fi

if [ ! -d "$BAM_DIR" ]; then
    echo "Error: BAM directory not found!"
    exit 1
fi

# Collect BAM files
BAMS=$(ls "$BAM_DIR"/*.bam | sort)
if [ -z "$BAMS" ]; then
    echo "Error: No BAM files found in $BAM_DIR"
    exit 1
fi

echo "Found BAM files:"
echo "$BAMS"

# ---------------------------
# Generate sample list automatically
# ---------------------------
SAMPLE_LIST="${OUTPUT_VCF%.vcf}_samples.txt"
echo "Generating sample list: $SAMPLE_LIST"
rm -f "$SAMPLE_LIST"
for bam in $BAMS; do
    basename "$bam" .bam >> "$SAMPLE_LIST"
done

echo "Sample names:"
cat "$SAMPLE_LIST"

# ---------------------------
# Run mpileup and pipe to VarScan
# ---------------------------
echo "Running samtools mpileup and VarScan..."
samtools mpileup -f "$REF" $BAMS | \
java -jar "$VARSCAN_PATH" mpileup2snp \
    --min-avg-qual 30 \
    --min-coverage 10 \
    --min-reads2 5 \
    --min-var-freq 0.25 \
    --p-value 0.01 \
    --min-freq-for-hom 0.75 \
    --output-vcf 1 \
    --vcf-sample-list "$SAMPLE_LIST" \
> "$OUTPUT_VCF"

echo "? VCF generated: $OUTPUT_VCF"

# ---------------------------
# Compress and index while keeping original VCF
# ---------------------------
echo "Compressing and indexing VCF..."
bgzip -c "$OUTPUT_VCF" > "${OUTPUT_VCF}.gz"
tabix -p vcf "${OUTPUT_VCF}.gz"

echo "? Original VCF kept: $OUTPUT_VCF"
echo "? Compressed and indexed VCF: ${OUTPUT_VCF}.gz"
echo "? Sample list file: $SAMPLE_LIST"