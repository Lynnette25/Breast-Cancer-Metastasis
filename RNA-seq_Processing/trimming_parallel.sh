#!/bin/bash

# This script trims RNA-seq reads in parallel using Trimmomatic and GNU Parallel.

# Set paths
RNA_HOME=/home/fheredia/RNS-seq/Metastasis_Project
RNA_REFS_DIR=/home/fheredia/RNS-seq/RNA_refs
adapters=$RNA_REFS_DIR/illumina_multiplex.fa
input_dir=/opt/data_hf/metastasis/RNA_seq
output_dir=$RNA_HOME/RNA_seq/trimmed
mkdir -p "$output_dir"

# Define trimming function
trim_sample() {
    base=$1
    r1="$input_dir/${base}_1.fq.gz"
    r2="$input_dir/${base}_2.fq.gz"
    out_p1="$output_dir/${base}_1_paired.fq.gz"
    out_u1="$output_dir/${base}_1_unpaired.fq.gz"
    out_p2="$output_dir/${base}_2_paired.fq.gz"
    out_u2="$output_dir/${base}_2_unpaired.fq.gz"

    echo "Trimming $base..."
    trimmomatic PE -threads 4 -phred33 \
        "$r1" "$r2" \
        "$out_p1" "$out_u1" \
        "$out_p2" "$out_u2" \
        ILLUMINACLIP:$adapters:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
}

export -f trim_sample
export input_dir output_dir adapters

# Run in parallel: 4 jobs at a time
parallel -j 4 trim_sample :::: sample_names.txt

echo "All samples trimmed in parallel."