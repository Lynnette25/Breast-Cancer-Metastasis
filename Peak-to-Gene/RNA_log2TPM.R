################################################################################
# Frances Heredia
# Franco Lab
# Description: This script performs the following tasks  
# 1. Reads in RNA-seq counts data
# 2. Filters for important genes based on a provided CSV file
# 3. Normalizes the data to log2(TPM + 1)
# 4. Saves the normalized data to a TSV file
################################################################################


library(data.table)
library(dplyr)


df <- read.csv("~/RNA_seq/counts_PE.csv", header = TRUE, sep = ",")

rownames(df) <- df$Geneid  

################################################################################
# This section filters the RNA-seq data for clinicaly relevant genes. If you wish to
# use all genes, skip to the normalization section.  
################################################################################
# Read the CSV file
data <- read.csv("~/PE_counts/Up-genes-three-way-gene_intersections_p_val.csv", header = TRUE)

selected_columns <- data %>% select(1, 3, 5, 7 )

# Flatten the data into a single list, ignoring the row and column structure
flattened_list <- unlist(selected_columns, use.names = FALSE)

# All the importtan genes for Met/Liver/Lung
clean_vector <- flattened_list[!is.na(flattened_list)]

df <- df[rownames(df) %in% clean_vector, ]

################################################################################
# Normalization to log2(TPM + 1)
################################################################################

sample_cols <- grep("^R", names(df), value = TRUE)

df[sample_cols] <- lapply(df[sample_cols], as.numeric)

for (sample in sample_cols) {
  length_kb <- df$Length / 1000
  df[[paste0(sample, "_RPK")]] <- df[[sample]] / length_kb
  sum_rpk <- sum(df[[paste0(sample, "_RPK")]], na.rm = TRUE)
  df[[paste0(sample, "_TPM")]] <- (df[[paste0(sample, "_RPK")]] / sum_rpk) * 1e6
}

# Supprimer les colonnes de comptage brut et RPK
df <- df %>% select(Geneid, starts_with("R") & ends_with("TPM"))

# Redéfinir les rownames à partir de la colonne GeneID

df$Geneid <- NULL

# Appliquer log2(TPM + 1)
tpm_log2 <- log2(df + 1)


write.table(tpm_log2,
            file = "~/ATAC_seq/Third_run/results_pipeline/Peak_to_Gene/Filtered_RNA_MET_Liver_Lung_Cancer_all_samples_PE.tsv",
            row.names = TRUE, quote = FALSE, sep = " ")



df  <- read.csv("~/ATAC_seq/Third_run/results_pipeline/Peak_to_Gene/Filtered_RNA_MET_Liver_Lung_Cancer_all_samples_PE.tsv", header = TRUE, sep = " ")

colnames(df) <- gsub("_TPM$", "", colnames(df))
df <- df[, !colnames(df) %in% "R49BR_3"]

desired_order <- c(
  "R8BR_1", "R8BR_2", "R49BR_1", "R49BR_2", "R39BR_1",
  "R8LIV_1", "R8LIV_2", "R39LIV_1", "R39LIV_2", "R49LIV_1", "R49LIV_2",
  "R18LIV_1", "R18LIV_2", "R23LIV_1", "R23LIV_2", "R33LIV_1", "R33LIV_2", "R36LIV_1", "R36LIV_2",
  "R8LUN_1", "R8LUN_2", "R49LUN_1", "R49LUN_2", "R18LUN_1", "R18LUN_2",
  "R23LUN_1", "R23LUN_2", "R33LUN_1", "R33LUN_2", "R36LUN_1", "R36LUN_2"
)

df <- df[, desired_order]

write.table(df,
            file = "~/ATAC_seq/Third_run/results_pipeline/Peak_to_Gene/Filtered_RNA_MET_Liver_Lung_Cancer_all_samples_PE.tsv",
            row.names = TRUE, quote = FALSE, sep = " ")



################################################################################
# This section splits the RNA-seq data into separate files for each sample type
# (Liver, Lung, Breast) for downstream analysis.			
# If you do not need this, you can skip it.
################################################################################


                      
rna_matrix <- fread("Peak_to_Gene/Filtered_RNA_MET_Liver_Lung_Cancer_all_samples_PE.tsv")

rna_LIV <- rna_matrix[, c("V1", grep("LIV", names(rna_matrix), value = TRUE)), with = FALSE]
rna_LUN <- rna_matrix[, c("V1", grep("LUN", names(rna_matrix), value = TRUE)), with = FALSE]
rna_BR <- rna_matrix[, c("V1", grep("BR", names(rna_matrix), value = TRUE)), with = FALSE]

fwrite(rna_LIV, "Peak_to_Gene/Filtered_RNA_MET_Liver_Lung_Cancer_LIV_only.tsv", sep = "\t")
fwrite(rna_LUN, "Peak_to_Gene/Filtered_RNA_MET_Liver_Lung_Cancer_LUN_only.tsv", sep = "\t")
fwrite(rna_BR, "Peak_to_Gene/Filtered_RNA_MET_Liver_Lung_Cancer_BR_only.tsv", sep = "\t")

