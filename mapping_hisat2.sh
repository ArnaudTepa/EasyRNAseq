#!/bin/bash  

# Specify the HISAT2 index directory and reference genome  
HISAT2_INDEX=$1  
OUTPUT_DIR=$3  
READS_DIR=$2  

# Number of threads to use  
threads=40  

# Log file to track processed samples  
LOG_FILE="$OUTPUT_DIR/alignment_log.txt"  

# Create the output directory and log file if they don't exist  
mkdir -p "$OUTPUT_DIR"  
touch "$LOG_FILE"  

# Iterate over all directories in the reads directory  
for dir in "$READS_DIR"/*; do  
    if [ -d "$dir" ]; then  
        echo "Processing directory: $dir"  

        # Find paired-end reads (assuming they end with _1.fq.gz and _2.fq.gz)  
        read1=$(find "$dir" -name "*_1.fq.gz")  
        read2=$(find "$dir" -name "*_2.fq.gz")  

        # Check if both read files exist  
        if [[ -f "$read1" && -f "$read2" ]]; then  
            # Generate output filenames  
            sample_name=$(basename "$dir")  
            output_prefix="$OUTPUT_DIR/${sample_name}"  

            # Check if the sorted BAM file already exists  
            if [[ ! -f "${output_prefix}.sorted.bam" ]]; then  
                # Check if this sample has already been processed    
                if ! grep -q "$sample_name" "$LOG_FILE"; then  
                    echo "Aligning $read1 and $read2"  

                    # Align reads with HISAT2  
                    hisat2 -x "$HISAT2_INDEX" -1 "$read1" -2 "$read2" -p "$threads" -t --summary-file "${output_prefix}.txt" --dta \
                    | samtools view -bS \
                    | samtools sort -o "${output_prefix}.sorted.bam"  
                    
                    # Index the sorted BAM file  
                    samtools index "${output_prefix}.sorted.bam"  
                    
                    echo "Alignment completed for $sample_name."  

                    # Log the completed sample name  
                    echo "$sample_name" >> "$LOG_FILE"  
                else  
                    echo "Alignment already done for $sample_name, skipping."  
                fi  
            else  
                echo "Sorted BAM file already exists for $sample_name, skipping."  
            fi  
        else  
            echo "Warning: Paired-end reads not found in $dir"  
        fi  
    fi  
done  

echo "All directories processed."
