################################################################################
# Frances Heredia
# Franco Lab
# Description: This script corrects for sequence bias in ATAC-seq data using 
# TOBIAS ATACorrect.
################################################################################

# Load necessary libraries
library(parallel)

## Setup

# Define the samples to be processed
samples <- c("Breast", "Liver", "Lung")

# Function to run TOBIAS ATACorrect for a given sample
run_tobias <- function(sample) {
    bam_file <- file.path("Footprints/merged/bam/", sample, paste0(sample, "_merged.bam"))
    peaks_file <- file.path("Footprints/merged/peaks/", sample, paste0(sample, "_merged.narrowPeak"))
    
    command <- sprintf("~/miniconda3/bin/TOBIAS ATACorrect --bam %s --genome ~/pepatac/genome_folder/alias/hg38/fasta/default/hg38.fa --peaks %s --outdir Footprints --prefix %s",
                       bam_file, peaks_file, sample)
                       
    # Execute TOBIAS command
    system(command)
}

## Execution

# Run the correction function in parallel using available cores
mclapply(samples, run_tobias, mc.cores = detectCores() - 4) # Use one less than total cores

## Notes
# - Ensure the TOBIAS package is properly installed in the conda environment specified.
# - Check that the file paths for BAM files and peak files are correct.
# - Adjust the number of cores based on your system's capacity for optimal performance.
