#!/bin/bash  

# Set the base directory for input files  
BASE_DIR=$1  

# Set the output directory for FastQC results  
OUTPUT_DIR=$2  

# Create the output directory if it does not exist  
mkdir -p "$OUTPUT_DIR"  

# Loop through each subdirectory within the base directory  
for SUBDIR in "$BASE_DIR"/*; do  
    # Check if it is a directory  
    if [ -d "$SUBDIR" ]; then  
        echo "Processing directory: $SUBDIR"  

        # Find .fastq.gz files in the current subdirectory  
        for FILE in "$SUBDIR"/*.fq.gz; do  
            # Check if the file exists  
            if [ -e "$FILE" ]; then  
                echo "Running fastqc on $FILE"  
                fastqc "$FILE" -o "$OUTPUT_DIR"   # Save output in the specified output directory  
            else  
                echo "No .fq.gz files found in $SUBDIR"  
            fi  
        done  
    fi  
done  

echo "Processing complete."
