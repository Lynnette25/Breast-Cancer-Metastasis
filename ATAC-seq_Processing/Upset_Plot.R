################################################################################
# Frances Heredia
# Franco Lab
# Description: This script performs the following tasks  
# 	  1) Adds annotations using specified gene bed and GTF files.
#     2) Retrieves bed files based on specified peaks.
#     3) Combines results from retrieved bed files.
#     4) Creates an UpSet plot of the combined results and saves it as a PNG file.
################################################################################



# Load necessary library
library(cisDynet)

# Add annotations using specified gene bed and GTF files.
addAnnotation(gene_bed = "beds/50/hg38_gene_standard.bed",
              gtf = "beds/50/Modified_Homo_sapiens.GRCh38.108.gtf",
              genome_size = "beds/50/hg38.chrom.size")

# Define peak files to be analyzed
peaks <- c("Early_Breast_consensus_04.bed", "Breast_consensus_04.bed", 
           "Liver_consensus_04.bed", "Lung_consensus_04.bed")

# Retrieve bed files based on specified peaks 
bed <- get_beds(peaks)

# Combine results from retrieved bed files
inputs <- get_combine_result(bed)

# Create an UpSet plot of the combined results
output_file <- "Back_to_Peaks/cisDynet_consensus_UpSet_plus_Early.png"
png(output_file, width = 10, height = 4, units = "in", res = 300)
plotUpset(inputs, bed)
dev.off()

# Notify user of completion
cat("Analysis complete. Output saved to", output_file, "\n")
