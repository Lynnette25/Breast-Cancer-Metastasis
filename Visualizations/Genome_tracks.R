################################################################################
# Frances Heredia
# Franco Lab
# Description: This script performs the following tasks  
# 1. Reads a CSV file containing genomic peak data.
# 2. For each gene, it calculates the minimum start and maximum end positions of the
#    peaks, expands these positions by 500 base pairs on both sides, and ensures that
#    the start position is not less than 1.
# 3. It then saves the calculated ranges to a new CSV file.
# 4. Finally, it reads the new CSV file and generates genome track plots for each
#    gene using the pyGenomeTracks tool, saving each plot as a PNG file.
################################################################################


library(dplyr)
library(tidyr)

# Load your data
df <- read.csv("Correlated_peak.csv", sep = "\t",stringsAsFactors = FALSE)


gene_ranges <- df_parsed %>%
  group_by(Gene, chr) %>%
  summarise(
    peak_start = min(start) - 500,
    peak_end = max(end) + 500,
    .groups = "drop"
  ) %>%
  # Ensure start is not less than 1
  mutate(peak_start = ifelse(peak_start < 1, 1, peak_start),
         range = paste0(chr, ":", peak_start, "-", peak_end))


# Print result
print(gene_ranges)

# Save to file
write.csv(gene_ranges, "Correlated_peak_ranges_expanded500bp.csv", row.names = FALSE)

# Load necessary libraries
library(dplyr)
library(readr)

# Read the CSV
df <- read_csv("Back_to_Peaks/Correlated_peak_ranges_expanded500bp.csv")


for (i in 1:nrow(df)) {
  gene_name <- df$Gene[i]
  # Use the 'range' column from the CSV for the main region
  region_str <- df$range[i]
  # Remove quotes if present
  region_str <- gsub('"', '', region_str)
  cat("Gene:", gene_name, "Region:", region_str, "\n")
  
  # Call pyGenomeTracks for each gene with the region
  cmd <- sprintf("pyGenomeTracks --region %s --tracks 3_Organ_tracks.ini -o Tracks_%s.png ",
                 region_str, gene_name)
  cat("Running:", cmd, "\n")
  system(cmd)
}
