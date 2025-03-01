#!/bin/bash

# This script performs quality trimming on paired-end FASTQ files using Trimmomatic.
# It processes all files in the specified input directory that match the pattern *_1.fastq,
# assuming paired files are named with _1.fastq and _2.fastq suffixes.

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_directory> <output_directory> <threads>"
    exit 1
fi

# Assign command-line arguments to variables
input_dir="$1"   # Directory containing input FASTQ files
output_dir="$2"  # Directory to save trimmed output files
threads="$3"     # Number of threads to use for Trimmomatic

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Define paths for Trimmomatic and adapter files
# Update these paths based on your Anaconda environment installation
trimmomatic="/home/arif/anaconda3/envs/rna_seq/share/trimmomatic-0.39-2/trimmomatic"
adapters="/home/arif/anaconda3/envs/rna_seq/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa"

# Loop through all forward read files in the input directory
for file in "$input_dir"/*_1.fastq; do
    # Extract the base name without the _1.fastq suffix
    BASENAME=$(basename "$file" _1.fastq)

    # Define file paths for forward and reverse reads
    forward_read="$input_dir/${BASENAME}_1.fastq"
    reverse_read="$input_dir/${BASENAME}_2.fastq"

    # Define file paths for output paired and unpaired reads
    paired_forward="$output_dir/${BASENAME}_1_paired.fastq"
    unpaired_forward="$output_dir/${BASENAME}_1_unpaired.fastq"
    paired_reverse="$output_dir/${BASENAME}_2_paired.fastq"
    unpaired_reverse="$output_dir/${BASENAME}_2_unpaired.fastq"

    # Run Trimmomatic for paired-end read trimming
    "$trimmomatic" PE -threads "$threads" -phred33 \
        "$forward_read" "$reverse_read" \
        "$paired_forward" "$unpaired_forward" \
        "$paired_reverse" "$unpaired_reverse" \
        ILLUMINACLIP:"$adapters":2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

    # Output message indicating completion of trimming for the current pair
    echo "Trimmed: $forward_read and $reverse_read using $threads threads"
done

# Final message indicating completion of the trimming process
echo "Trimmomatic processing complete. Results saved in '$output_dir'."
