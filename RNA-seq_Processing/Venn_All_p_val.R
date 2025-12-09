################################################################################
# Frances Heredia
# Franco Lab
# Description: This script performs the following tasks  
# 1. Reads differential expression results from CSV files for Liver, Breast, and Lung tissues.
# 2. Merges p-values from a separate CSV file into each tissue's data.
# 3. Filters genes based on specified thresholds for P.Value, tissue-specific p-values, and logFC.
# 4. Creates a Venn diagram to visualize the overlap of upregulated genes across the three tissues.
# 5. Extracts and saves the gene intersections to a CSV file. 
################################################################################


library(VennDiagram)
library(readr)
library(dplyr)
library(tidyr)


P_value <- read_csv("FULL_gene_data_PE_MET_GTEX_p_value.csv")

#P-value cols:  "Breast_p_value"     "Liver_p_value"      "Lung_p_value"      

Liver_data <- read_csv("PE_counts/differential_expression_Liver_results_Met_Gtex_clean.csv")
Liver_data <- as.data.frame(Liver_data)
rownames(Liver_data) <- Liver_data$...1
head(Liver_data)

Liver_data$ID <- as.character(Liver_data$ID)
P_value$Geneid <- as.character(P_value$Geneid)

# Merge the p-value into Liver_data
Liver_data_merged <- Liver_data %>%left_join(P_value %>% select(Geneid, Liver_p_value), by = c("ID" = "Geneid"))

Breast_data <- read_csv("PE_counts/differential_expression_Breast_results_Met_Gtex_clean.csv")  # Adjust file name as needed
Breast_data <- as.data.frame(Breast_data)
rownames(Breast_data) <- Breast_data$...1
head(Breast_data)

Breast_data$ID <- as.character(Breast_data$ID)
# Merge the p-value into Liver_data
Breast_data_merged <- Breast_data %>%left_join(P_value %>% select(Geneid, Breast_p_value), by = c("ID" = "Geneid"))

Lung_data <- read_csv("PE_counts/differential_expression_Lung_results_Met_Gtex_clean.csv")  # Adjust file name as needed
Lung_data <- as.data.frame(Lung_data)
rownames(Lung_data) <- Lung_data$...1
head(Lung_data)

Lung_data$ID <- as.character(Lung_data$ID)
# Merge the p-value into Liver_data
Lung_data_merged <- Lung_data %>%left_join(P_value %>% select(Geneid, Lung_p_value), by = c("ID" = "Geneid"))

threshold_pvalue <- 0.05  # Example threshold for P.Value
threshold_logFC <- 0.1

Liver_genes <- Liver_data_merged %>%  filter(P.Value < threshold_pvalue & Liver_p_value < threshold_pvalue & logFC > threshold_logFC)
Breast_genes <- Breast_data_merged %>%  filter(P.Value < threshold_pvalue & Breast_p_value < threshold_pvalue & logFC > threshold_logFC)
Lung_genes <- Lung_data_merged %>%  filter(P.Value < threshold_pvalue & Lung_p_value < threshold_pvalue & logFC > threshold_logFC)

Liver_genes <-  Liver_genes$ID
Breast_genes <- Breast_genes$ID
Lung_genes <- Lung_genes$ID


length(Liver_genes)
length(Breast_genes)
length(Lung_genes)

head(Liver_genes)
head(Breast_genes)
head(Lung_genes)

venn_file_path <- "PE_counts/Up-genes-three-way-venn_diagram_p_val.png"


venn.plot <- venn.diagram(
  x = list(
    Set1 = Liver_genes,
    Set2 = Breast_genes,
    Set3 = Lung_genes
  ),
  category.names = c("Liver", "Breast", "Lung"),
  filename = venn_file_path,
  lwd = 2,
  fill = c("red", "blue", "green"),
  cex = 1.5,  # Size of the text
  cat.cex = 1.5,  # Size of the category label
  cat.col = c("red", "blue", "green"),
  scaled = TRUE,  # Scale the circles based on the data
  cat.pos = c(285, 135, 90),  # Position of category labels
  margin = 0.1
)



# Extracting the intersections
intersect1_2 <- setdiff(intersect(Liver_genes, Breast_genes), Lung_genes)
intersect1_3 <- setdiff(intersect(Liver_genes, Lung_genes), Breast_genes)
intersect2_3 <- setdiff(intersect(Breast_genes, Lung_genes), Liver_genes)
intersect_all <-Reduce(intersect, list(Liver_genes, Breast_genes, Lung_genes))

unique_liver <- setdiff(Liver_genes, union(Breast_genes, Lung_genes))
unique_breast <- setdiff(Breast_genes, union(Liver_genes, Lung_genes))
unique_lung <- setdiff(Lung_genes, union(Liver_genes, Breast_genes))

cat("Number in intersect1_2:", length(intersect1_2), "\n")
cat("Number in intersect1_3:", length(intersect1_3), "\n")
cat("Number in intersect2_3:", length(intersect2_3), "\n")
cat("Number in all three:", length(intersect_all), "\n")
cat("Unique to Liver:", length(unique_liver), "\n")
cat("Unique to Breast:", length(unique_breast), "\n")
cat("Unique to Lung:", length(unique_lung), "\n")

# Updated list to include unique genes
intersect_data <- list(
  Unique_Liver = unique_liver,
  Unique_Breast = unique_breast,
  Unique_Lung = unique_lung,
  Common_Liver_Breast = intersect1_2,
  Common_Liver_Lung = intersect1_3,
  Common_Breast_Lung = intersect2_3,
  Common_All = intersect_all
)

# Convert to data frame for CSV
intersect_df <- do.call(data.frame, lapply(intersect_data, function(x) {
  as.data.frame(matrix(c(x, rep(NA, max(sapply(intersect_data, length)) - length(x))), nrow=max(sapply(intersect_data, length))))
}))

for (col in names(intersect_df)) {
  non_na_count <- sum(!is.na(intersect_df[[col]]))
  cat("Number of non-NA entries in", col, ":", non_na_count, "\n")
}

# Write to CSV
write_csv(intersect_df, "PE_counts/Up-genes-three-way-gene_intersections_p_val.csv")
