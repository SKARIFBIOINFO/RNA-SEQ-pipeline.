#!/bin/bash

# This script aligns paired-end FASTQ files to a pre-indexed E. coli genome using Bowtie2.
# It assumes that the genome has already been indexed using bowtie2-build.
# The script processes all input files, aligns them, converts them to BAM format,
# sorts and indexes the output for downstream analysis.

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input_directory> <output_directory> <bowtie2_index_prefix> <threads>"
    exit 1
fi

# Assign input arguments to variables
input_dir="$1"      # Directory containing paired-end FASTQ files
output_dir="$2"     # Directory where output files will be stored
index_dir="$3"      # Directory containing the Bowtie2 index files
no_of_threads="$4"  # Number of threads to use for processing

# Create the output directory if it does not exist
mkdir -p "$output_dir"

# Define the Bowtie2 index file prefix
index_prefix="$index_dir/E_coli_genome_index"

# Loop through all forward read files in the input directory
for forward_read in "$input_dir"/*_1_paired.fastq; do
    # Extract the base name without the _1_paired.fastq suffix
    BASENAME=$(basename "$forward_read" _1_paired.fastq)

    # Define the corresponding reverse read file
    reverse_read="$input_dir/${BASENAME}_2_paired.fastq"

    # Define output file paths
    sam_output="$output_dir/${BASENAME}.sam"               # SAM file (alignment output)
    bam_output="$output_dir/${BASENAME}.bam"               # BAM file (compressed SAM)
    sorted_bam_output="$output_dir/${BASENAME}_sorted.bam" # Sorted BAM file

    echo "Processing input sample: $BASENAME"

    # Run Bowtie2 alignment using the pre-built genome index
    bowtie2 -x "$index_prefix" -1 "$forward_read" -2 "$reverse_read" -S "$sam_output" -p "$no_of_threads"

    # Convert SAM to BAM to reduce file size and enable faster downstream analysis
    samtools view -@ "$no_of_threads" -bS "$sam_output" > "$bam_output"

    # Sort BAM file to organize alignments
    samtools sort -@ "$no_of_threads" -o "$sorted_bam_output" "$bam_output"

    # Index the sorted BAM file for efficient access in downstream analysis
    samtools index "$sorted_bam_output"

    echo "Alignment completed for: $BASENAME"
done

echo "Bowtie2 alignment pipeline finished. Results are in '$output_dir'."
