################################################################################
# Frances Heredia
# Franco Lab
# Description: This script combines all bound TFs and their scores into a
# single matrix with peaks as rows and TFs as columns. This matrix can be used
# as input for TFSEE.
# It workes for a single tissue, in this case Breast. It has to be modified for other tissues. 
################################################################################


library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(stringr)

# Path to all your footprint bed files
footprint_dir <- "Footprints/Breast_Liver"

# Find all bound.bed files
files <- list.files(footprint_dir, pattern = "_Breast_bound\\.bed$", recursive = TRUE, full.names = TRUE)

# Function to read each file and extract relevant info
read_tf_file <- function(file) {
  df <- read_tsv(file, col_names = FALSE, show_col_types = FALSE)
  
  # Columns:
  # X7, X8, X9 -> peak coordinates
  # X4 -> TF name
  # X15 -> score (last column)
  df %>%
    mutate(
      TF = X4,
      peak_id = paste0(X7, ":", X8, "-", X9),
      score = X15
    ) %>%
    select(peak_id, TF, score)
}

# Read all files and combine
all_data <- map_dfr(files, read_tf_file)

# Collapse duplicates: keep max score per peak_id + TF
all_data_clean <- all_data %>%
  group_by(peak_id, TF) %>%
  summarise(score = max(score, na.rm = TRUE), .groups = "drop")

# Now pivot safely
wide_matrix <- all_data_clean %>%
  pivot_wider(
    names_from = TF,
    values_from = score,
    values_fill = list(score = 0)
  )


# Arrange rows by peak_id
wide_matrix <- wide_matrix %>% arrange(peak_id)

# Write final matrix
write_csv(wide_matrix, "Merged_Breast_Footprints.csv")

# Repeat the above steps for Liver and Lung tissues, changing the pattern in list.files accordingly.


################################################################################
# Combine multiple tissues into a rds list of matrices for TFSEE

################################################################################

# Combine footprinting files into one peak matrix
```R
library(data.table)
library(dplyr)
library(tidyr)

# Load your merged footprint matrices (each already peak × TF)
breast <- fread("Merged_Breast_Footprints.csv")
liver  <- fread("Merged_Liver_Footprints.csv")
lung   <- fread("Merged_Lung_Footprints.csv")

# Align TF columns across tissues
common_tfs <- Reduce(intersect, list(colnames(breast), colnames(liver), colnames(lung)))
common_tfs <- setdiff(common_tfs, "peak_id")

breast <- breast[, c("peak_id", common_tfs), with = FALSE]
liver  <- liver[,  c("peak_id", common_tfs), with = FALSE]
lung   <- lung[,   c("peak_id", common_tfs), with = FALSE]

# Convert each to matrix with peak_id as rownames
fp_list <- list(
  Breast = as.data.frame(breast),
  Liver  = as.data.frame(liver),
  Lung   = as.data.frame(lung)
)

fp_list <- lapply(fp_list, function(x) {
  rownames(x) <- x$peak_id
  x$peak_id <- NULL
  as.matrix(x)
})

# Save
saveRDS(fp_list, file = "footprinting_matrices_list.rds")

