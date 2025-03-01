#!/bin/bash

# This script performs gene expression quantification using RSEM for paired-end RNA-Seq data.
# It processes multiple samples, aligning them to a reference genome using Bowtie2 and calculating transcript abundance.

# ---------------------
# Command-line Arguments
# ---------------------
# $1 - RAW_DATA_DIR: Directory containing raw paired-end FASTQ files.
# $2 - OUTPUT_DIR: Directory where RSEM output will be stored.
# $3 - REFERENCE_PREFIX: Prefix of the pre-built RSEM reference (must be prepared in advance using rsem-prepare-reference).
# $4 - NUM_THREADS: Number of CPU threads to use for parallel processing.

RAW_DATA_DIR="$1"      # Input directory containing raw FASTQ files
OUTPUT_DIR="$2"        # Output directory for RSEM results
REFERENCE_PREFIX="$3"  # Prefix for pre-built RSEM reference files (prepared earlier using rsem-prepare-reference)
NUM_THREADS="$4"       # Number of threads for parallel execution

# Ensure that the output directory exists
mkdir -p "$OUTPUT_DIR"

# ---------------------
# Ensure that the RSEM reference was built earlier
# ---------------------
if [[ ! -f "${REFERENCE_PREFIX}.grp" || ! -f "${REFERENCE_PREFIX}.seq" ]]; then
    echo "Error: RSEM reference files not found. Ensure you have run 'rsem-prepare-reference' before executing this script."
    exit 1  # Exit if reference files are missing
fi

# ---------------------
# Loop through all forward read files and process them
# ---------------------
for FORWARD_READ in "$RAW_DATA_DIR"/*_1.fastq; do
    # Extract the sample name from the forward read file name
    BASENAME=$(basename "$FORWARD_READ" _1.fastq)

    # Identify the corresponding reverse read file
    REVERSE_READ="$RAW_DATA_DIR/${BASENAME}_2.fastq"

    # Define output prefix for RSEM results
    OUTPUT_PREFIX="$OUTPUT_DIR/${BASENAME}_rsem"

    # Check if the corresponding reverse read file exists
    if [[ ! -f "$REVERSE_READ" ]]; then
        echo "Reverse read file for $BASENAME not found. Skipping."
        continue  # Skip this sample and move to the next one
    fi

    # ---------------------
    # Run RSEM for expression quantification
    # ---------------------
    rsem-calculate-expression --paired-end \  # Specifies paired-end sequencing data
                              --bowtie2 \      # Uses Bowtie2 as the alignment tool
                              --num-threads "$NUM_THREADS" \  # Uses multiple threads for faster processing
                              "$FORWARD_READ" "$REVERSE_READ" \  # Input FASTQ files
                              "$REFERENCE_PREFIX" \  # Pre-built reference genome prefix
                              "$OUTPUT_PREFIX"  # Output prefix for RSEM results

    # ---------------------
    # Error Handling
    # ---------------------
    if [[ $? -ne 0 ]]; then
        echo "Error in rsem-calculate-expression for sample $BASENAME. Exiting."
        exit 1  # Exit the script if RSEM fails
    fi

    echo "Quantification for sample $BASENAME completed."
done

# ---------------------
# Completion Message
# ---------------------
echo "All samples processed successfully."
