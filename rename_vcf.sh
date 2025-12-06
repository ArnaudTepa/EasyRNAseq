#!/bin/bash

# Usage: ./rename_vcf.sh <vcf_file> <sample_file> [output_vcf]
VCF_FILE=$1
SAMPLE_FILE=$2
OUTPUT_VCF=$3

# Check arguments
if [ -z "$VCF_FILE" ] || [ -z "$SAMPLE_FILE" ]; then
    echo "Usage: $0 <vcf_file> <sample_file> [output_vcf]"
    exit 1
fi

# Set default output name if not provided
if [ -z "$OUTPUT_VCF" ]; then
    OUTPUT_VCF="merged_renamed.vcf"
fi

# Check if bcftools is installed
if ! command -v bcftools &> /dev/null; then
    echo "Error: bcftools is not installed. Please install bcftools and retry."
    exit 1
fi

# Validate input files
if [ ! -f "$VCF_FILE" ]; then
    echo "Error: VCF file '$VCF_FILE' not found."
    exit 1
fi

if [ ! -f "$SAMPLE_FILE" ]; then
    echo "Error: Sample file '$SAMPLE_FILE' not found."
    exit 1
fi

# Extract old sample names from VCF
bcftools query -l "$VCF_FILE" > old_names.txt

# Count old and new names
OLD_COUNT=$(wc -l < old_names.txt)
NEW_COUNT=$(wc -l < "$SAMPLE_FILE")

if [ "$OLD_COUNT" -ne "$NEW_COUNT" ]; then
    echo "Error: Number of samples in VCF ($OLD_COUNT) does not match sample.txt ($NEW_COUNT)."
    exit 1
fi

# Create rename.map file
paste old_names.txt "$SAMPLE_FILE" > rename.map

# Perform renaming
bcftools reheader -s rename.map "$VCF_FILE" -o "$OUTPUT_VCF"

# Verify and print new header
if [ $? -eq 0 ]; then
    echo "Renaming completed successfully. New VCF: $OUTPUT_VCF"
    grep '#CHROM' "$OUTPUT_VCF"
else
    echo "Error: Renaming failed."
    exit 1
fi
