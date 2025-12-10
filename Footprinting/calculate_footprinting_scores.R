################################################################################
# Frances Heredia
# Franco Lab
# Description: This script calculates footprinting scores using TOBIAS ScoreBigWig.
################################################################################

## Setup

# Define the samples to be processed
samples <- c("Breast", "Liver", "Lung")

# Loop through each sample and calculate footprinting scores
for (sample in samples) {
    signal_file <- paste0("Footprints/", sample, "_corrected.bw")
    regions_file <- file.path("Footprints/merged/peaks/", sample, paste0(sample, "_merged.narrowPeak"))
    output_file <- paste0("Footprints/", sample, "_footprints.bw")
    
    # Run TOBIAS ScoreBigwig command
    command <- sprintf("/home/fheredia/miniconda3/bin/TOBIAS ScoreBigwig --signal %s --regions %s --cores 16 --output %s",
                       signal_file, regions_file, output_file)
    
    system(command)
}

## Notes
# Ensure that:
# - The corrected BigWig files are located in the specified directory.
# - The peak region files are also correctly placed in the specified path.
# - Adjust the number of cores if needed based on your system's performance.
