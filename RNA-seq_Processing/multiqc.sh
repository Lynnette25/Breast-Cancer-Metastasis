#!/bin/bash

# This script runs MultiQC to aggregate FastQC results in the current directory

# Get the current directory
current_dir=$(pwd)

# Run MultiQC to aggregate the FastQC results
multiqc "$current_dir" -o "$current_dir"

echo "FastQC and MultiQC analysis completed."
