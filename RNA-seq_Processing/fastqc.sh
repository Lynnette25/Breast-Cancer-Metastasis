#!/bin/bash

# This  script runs FastQC on all .fq.gz files in the specified directory
# and outputs the results to the current working directory.

# Get the current directory
current_dir=$(pwd)

# Loop through each .fq.gz file in the specified directory
for filename in /opt/data_hf/metastasis/RNA_seq/*.fq.gz
do
  fastqc -o "$current_dir" "$filename"
done
echo "FastQC processing complete."