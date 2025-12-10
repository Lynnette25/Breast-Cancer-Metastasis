################################################################################
# Frances Heredia
# Franco Lab
# Description: This script performs the following tasks  
# 1. Reads in ATAC-seq counts data
# 2. Normalizes the data to log2(CPM + 1)
# 3. Saves the normalized data to a TSV file	
################################################################################

library(data.table)
library(dplyr)

df <- fread("Back_to_Peaks/counts_file_peaks.csv")

# Identify sample columns (assuming they are all columns after the annotation columns)
# Adjust the column indices based on your data structure
annotation_cols <- c("seqnames" , "start" , "end" , "width" , "strand" , "X32")
sample_cols <- setdiff(names(df), annotation_cols)

# Calculate total counts per sample
total_counts <- df[, lapply(.SD, sum), .SDcols = sample_cols]

# Calculate CPM for each sample
# CPM = (counts / total counts) * 1,000,000
for (sample in sample_cols) {
  total <- total_counts[[sample]]
  df[, paste0(sample, "_CPM") := get(sample) / total * 1e6]
}
  
# Select 'seqnames  start    end' and all columns with '_CPM' in their names
filtered_df <- df %>% select(seqnames,  start,    end, matches("_CPM"))

colnames(filtered_df) <- c("seqnames" , "start" , "end", gsub("_CPM$", "", colnames(filtered_df)[-c(1, 2, 3)]))

# View the resulting data frame
head(filtered_df)

# Save coordinates in a single col

filtered_df$PeakID <- paste0(filtered_df$seqnames, ":", filtered_df$start, "-", filtered_df$end)
filtered_df$seqnames <- NULL
filtered_df$start <- NULL
filtered_df$end <- NULL

filtered_df <- as.data.frame(filtered_df)


# Set 'Peak_ID' as row names
numeric_cols <- sapply(filtered_df, is.numeric)
numeric_data <- filtered_df[, numeric_cols]

# Calculate log2 CPM
cpm_log2 <- log2(numeric_data + 1)
cpm_log2$Peak_ID <- filtered_df$PeakID

# Move PeakID to the front
cpm_log2_df <- cpm_log2 %>%
  select(Peak_ID, everything())

write.table(cpm_log2_df, file="Back_to_Peaks/ATAC_lod2CPM.tsv", row.names=FALSE, quote=FALSE, sep=" ")