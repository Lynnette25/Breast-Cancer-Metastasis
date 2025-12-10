################################################################################
# Frances Heredia
# Franco Lab
# Description: This file contains code snippets for data wrangling tasks related to downstream peak-to-gene analysis.
# Tasks include:
# 1. Creating a BED file of peaks from a counts CSV file.
# 2. Generating a range CSV file for genome tracks visualization.
# 3. Creating a short BED file from the Top 10 Peaks per Gene file for genome tracks.   
################################################################################



# Creating the bed file of peaks to use in downstream analysis

```python
import pandas as pd

df = pd.read_csv('Back_to_Peaks/counts_file_peaks.csv')

# Adjust 'start' to 0-based for BED
df['start_bed'] = df['start'] - 1

# Prepare the BED columns
bed_df = df[['seqnames', 'start_bed', 'end']]

# Save to BED file
bed_df.to_csv('Back_to_Peaks/Peaks.bed', sep='\t', header=False, index=False)
```


############################################################################################
# Creating the range csv to be used in genome tracks
```R
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

```


############################################################################################

# Creating the short bed file from the Top 10 Peaks per Gene file, for genome tracks
```R
library(dplyr)
library(tidyr)

# Load the Top 10 Peaks per Gene file
df <- read.csv("Top10Peaks_Per_Gene.csv", stringsAsFactors = FALSE)

# Parse the Peak column into chr, start, end
bed_df <- df %>%
  separate(Peak, into = c("chr", "coords"), sep = ":", remove = FALSE) %>%
  separate(coords, into = c("start", "end"), sep = "-", convert = TRUE) %>%
  select(chr, start, end, Gene, correlations)

# Write as a BED file (tab-separated, no headers, no row names)
write.table(
  bed_df,
  file = "Top10Peaks_Per_Gene.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

```

