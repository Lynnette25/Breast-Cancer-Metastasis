################################################################################
# Frances Heredia
# Franco Lab
# Description: This script performs the following tasks  
#       This script performs differential expression analysis using DESeq2 with SVA correction.
#       It generates PCA plots, volcano plots, MA plots, and heatmaps of significant genes.
#	   It takes as input a counts file and a design file, and outputs results to a specified directory.
# 	   Usage: Rscript Differential_analisis_SV.R <counts_file> <design_file>
################################################################################

library(pheatmap)
library(DESeq2)
library(sva)
library(ggplot2)
library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(rafalib)
library(scales)
library(amap)
library(EnhancedVolcano)
library(ggpubr)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
file <- args[1]         # First argument is the file (counts data)
design_file <- args[2]  # Second argument is the design file

# Step 1: Read the count data and design matrix
raw_counts <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
#raw_counts <- raw_counts %>% mutate(gene_name = ifelse(gene_name == "" | is.na(gene_name), Geneid, gene_name))

print(colnames(raw_counts))

# Identify duplicate row names
duplicates <- raw_counts$gene_name[duplicated(raw_counts$gene_name) | duplicated(raw_counts$gene_name, fromLast = TRUE)]
print(duplicates)
# Modify row names to make them unique
for (dup in duplicates) {
  # Create a new name with an underscore and the count of occurrences
  count <- sum(raw_counts$gene_name == dup)
  new_names <- paste(dup, 1:count, sep = "_")
  # Replace the old names with new unique names
  raw_counts$gene_name[raw_counts$gene_name == dup] <- new_names
}

xp_design <- read.csv(design_file, header = TRUE, stringsAsFactors = FALSE)

print(head(xp_design$sample))
output_dir <- paste0("NEW_", tools::file_path_sans_ext(basename(design_file)))
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Prepare the experimental design
xp_design$organ <- make.names(xp_design$organ)
xp_design$patient <- make.names(xp_design$patient)
xp_design$batch <- make.names(xp_design$patient)
xp_design$met <- make.names(xp_design$met)

# Select only the columns in raw_counts that exist in xp_design$sample
common_samples <- intersect(xp_design$sample, colnames(raw_counts))

counts <- raw_counts %>%
    column_to_rownames("gene_name") %>%
    select(one_of(common_samples))

filter_threshold <- 10  # Minimum count threshold for expression
sample_cutoff <- 3      # Minimum number of samples where the count is above the threshold

counts <- counts[rowSums(counts >= filter_threshold) >= sample_cutoff, ]

print("counts after filtering")
print(nrow(counts))
print(paste("Number of common samples:", length(common_samples)))

# If no common samples are found, skip this file
if (length(common_samples) == 0) {
    print(paste("Skipping", file, ": No matching samples between design and counts file."))
    next
}

# Prepare the design matrix and DESeq2 object
group <- factor(xp_design$organ[xp_design$sample %in% common_samples])
xp_design_filtered <- xp_design[xp_design$sample %in% common_samples, ]
design <- model.matrix(~ group)
dex <- DESeqDataSetFromMatrix(countData = counts, colData = xp_design_filtered, design = ~ batch + met)



# Prepare your design matrix and DESeq2 object (modify as needed)

dat <- as.matrix(assay(dex))
mm <- model.matrix(~ batch + met, colData(dex))
mm0 <- model.matrix(~ 1, colData(dex))
n.sv <- num.sv(dat, mm, method="leek")
#n.sv <- 2
if (n.sv > 4) {
    n.sv <- 4}

o <- svaseq(dat, mm, mm0, n.sv=n.sv) 

if (n.sv == 1) {
    design <- ~ V1 + met
} else if (n.sv == 2) {
    design <- ~ V1 + V2 + met
} else if (n.sv == 3) {
    design <- ~ V1 + V2 + V3 + met
} else if (n.sv == 4){
    design <- ~ V1 + V2 + V3 + V4+ met  # Default case for more than 3 SVs
}

# Extract the surrogate variables
W <- o$sv  # Surrogate variables

#  Create the design matrix combining surrogate variables and your variable of interest (e.g., "met")
# Make sure to convert W to a data frame
W_df <- as.data.frame(W)


# Combine your design data with the surrogate variables
colData_combined <- cbind(xp_design, W_df)

# Creation of the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData_combined, design = design  )


# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract results for differential expression analysis
res <- results(dds)

# Print the DESeqResults object
print(res)

# Explicitly specify the comparison of interest: treated vs untreated
all_genes_results <- results(dds, contrast = c("met","Primary", "Metastasis"))

# Add a column for gene names
all_genes_results <- as.data.frame(all_genes_results)
all_genes_results$gene <- rownames(all_genes_results)

print(head(all_genes_results))

# Perform PCA on the regularized log-transformed data (without SVs)
rld_no_svs <- rlog(dds, blind = FALSE)
pcaData_no_svs <- plotPCA(rld_no_svs, intgroup = c("met", "patient"), returnData = TRUE)
percentVar_no_svs <- round(100 * attr(pcaData_no_svs, "percentVar"))

# Save the PCA plot without SVs to a PNG file
output_file <- file.path(output_dir, paste0("PCA_plot_no_correction_", tools::file_path_sans_ext(basename(design_file)), ".png"))
png(output_file, width = 800, height = 600)
ggplot(pcaData_no_svs, aes(PC1, PC2, color = met)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar_no_svs[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_no_svs[2], "% variance")) +
  ggtitle("PCA with no Correction") +
  coord_fixed()

dev.off()

pcaData_no_svs <- plotPCA(rld_no_svs, intgroup = c("organ", "patient"), returnData = TRUE)
percentVar_no_svs <- round(100 * attr(pcaData_no_svs, "percentVar"))

# Save the PCA plot without SVs to a PNG file
output_file <- file.path(output_dir, paste0("PCA_plot_no_correction_", tools::file_path_sans_ext(basename(design_file)), "_organ.png"))
png(output_file, width = 800, height = 600)
ggplot(pcaData_no_svs, aes(PC1, PC2, color = organ)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar_no_svs[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_no_svs[2], "% variance")) +
  ggtitle("PCA with no Correction") +
  coord_fixed()

dev.off()
# Create a basic volcano plot and save it to a file
output_file <- file.path(output_dir, paste0("Volcano_plot_", tools::file_path_sans_ext(basename(design_file)), ".png"))
png(output_file, width = 800, height = 600)
EnhancedVolcano(toptable = all_genes_results,
                x = "log2FoldChange",
                y = "padj",
                lab = all_genes_results$gene,
                xlim = c(-7, 7),
                ylim = c(0, 7.5),
                pCutoff = 0.1,  # Adjusted p-value cutoff
                FCcutoff = 0.5)
dev.off()

output_file <- file.path(output_dir, paste0("MA_plot_", tools::file_path_sans_ext(basename(design_file)), ".png"))
png(output_file, width = 800, height = 600)
resLFC <- lfcShrink(dds = dds, 
                  res = res,
                  type = "apeglm",
                  coef ="met_Primary_vs_Metastasis") # corresponds to "infected_Pseudomonas_syringae_DC3000_vs_mock" comparison
plotMA(resLFC, alpha = 0.01)
dev.off()


# Specify output file for the MA plot
output_file <- file.path(output_dir, paste0("GGMA_plot_", tools::file_path_sans_ext(basename(design_file)), ".png"))
# Create the MA plot using ggmaplot from ggpubr
png(output_file, width = 800, height = 600)

ggmaplot(
  data = all_genes_results,           # Input data
  fdr = 0.05,                         # False discovery rate
  fc = 1.5,                           # Fold change threshold
  genenames = all_genes_results$gene_name,  # Gene names for labeling
  palette = c("#B31B21", "#1465AC", "darkgray"),  # Custom color palette
  top = 15,                           # Number of top genes to label
  main = "MA-Plot: Primary vs Metastasis",  # Plot title
  xlab = "Log2 mean expression",      # X-axis label
  ylab = "Log2 fold change",          # Y-axis label
  size = 0.5, 
  legend = "top", 
  font.label = c("bold", 11),
  font.legend = "bold",
  font.main = "bold",
  ggtheme = ggplot2::theme_minimal()
)

dev.off()


# Define the cutoffs
padj_cutoff <- 0.05
log2FC_cutoff <- 1

# Filter the differential genes based on the cutoffs
significant_genes <- res %>% 
  as.data.frame() %>% 
  rownames_to_column("gene_name") %>% 
  filter(padj < padj_cutoff & abs(log2FoldChange) >= log2FC_cutoff) %>% 
  arrange(desc(log2FoldChange), desc(padj))

output_file <- file.path(output_dir, paste0("Adjusted_pvalue_distribution_", tools::file_path_sans_ext(basename(design_file)), ".png"))
png(output_file, width = 800, height = 600)
# Distribution of adjusted p-values
hist(res$padj, col="lightblue", main = "Adjusted p-value distribution")
dev.off()


# Distribution of non-adjusted p-values
output_file <- file.path(output_dir, paste0("Non-adjusted_pvalue_distribution_", tools::file_path_sans_ext(basename(design_file)), ".png"))
png(output_file, width = 800, height = 600)
hist(res$pvalue, col="grey", main = "Non-adjusted p-value distribution")
# Save the non-adjusted p-value distribution plot
dev.off()

# Save the table of significant genes to a CSV file
output_file <- file.path(output_dir, paste0("Significant_genes_", tools::file_path_sans_ext(basename(design_file)), ".csv"))
write.csv(significant_genes, output_file, row.names = FALSE)



# Step 2: Perform SVA to identify 3 surrogate variables
dex <- estimateSizeFactors(dex)
norm.cts <- counts(dex, normalized=TRUE)
mm <- model.matrix(~ group, colData(dex))
mm0 <- model.matrix(~ 1, colData(dex))
norm.cts <- norm.cts[rowSums(norm.cts) > 0,]

# Perform SVA to detect  SVs
fit <- svaseq(norm.cts, mod = mm, mod0 = mm0, n.sv = n.sv)

# Step 3: Define the SVA correction function
svaBatchCor <- function(dat, mmi, mm0, n.sv=NULL){
  dat <- as.matrix(dat)
  Y <- t(dat)
  if(is.null(n.sv))   n.sv <- num.sv(dat, mmi, method="leek")
  o <- svaseq(dat, mmi, mm0, n.sv=n.sv)
  W <- o$sv
  alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
  o$corrected <- t(Y - W %*% alpha)
  return(o)
}

# Step 4: Use the correction function
moda <- model.matrix(~ 1 + met, colData(dex))  # Assuming 'group' is your main variable
modb <- model.matrix(~ 1, colData(dex))
correcteddata <- svaBatchCor(assay(dex), moda, modb, n.sv = n.sv)
setest <- SummarizedExperiment(correcteddata$corrected, colData = colData(dex))

correctedse <- DESeqTransform(setest)

rownames(correctedse) <- rownames(counts)

# Step 5: Filter significant DEGs
res <- results(DESeq(dex))
res <- res[!is.na(res$padj), ]  # Remove NA padj entries
# Add a column for gene names
corrected_results <- as.data.frame(res)
corrected_results$gene <- rownames(corrected_results)

# Print summary to check results
print(head(corrected_results))



# Ensure correctedse has valid row names
if (any(is.na(rownames(correctedse)))) {
    warning("correctedse contains NA row names.")
    correctedse <- correctedse[!is.na(rownames(correctedse)), ]  # Remove NA row names
}

# Check for significant DEGs
if (nrow(res) > 0 && any(res$padj < 0.0005)) {
    correctedse <- correctedse[rownames(correctedse) %in% rownames(res[res$padj < 0.0005,]), ]
} else {
    warning("No significant DEGs found with the specified threshold or all entries are NA.")
}


pcaData <- plotPCA(correctedse, intgroup = c( c("organ", "batch")), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
# Apply pseudo-log transformation to the PCA data directly
pcaData <- pcaData %>%
  mutate(PC1 = pseudo_log_trans(base = 10)$transform(PC1),
         PC2 = pseudo_log_trans(base = 10)$transform(PC2))

# Now plot the transformed PCA data without any additional axis scaling

# Save the PCA plot without SVs to a PNG file
output_file <- file.path(output_dir, paste0("PCA_plot_met_vs_batch_", tools::file_path_sans_ext(basename(design_file)), ".png"))
png(output_file, width = 800, height = 600)
ggplot(pcaData, aes(PC1, PC2, color = organ)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle(paste0("PCA after SVA Correction with ", n.sv," SVs ")) +
  coord_fixed()

dev.off()


pcaData <- plotPCA(correctedse, intgroup = c( c("organ", "patient")), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
# Apply pseudo-log transformation to the PCA data directly
pcaData <- pcaData %>%
  mutate(PC1 = pseudo_log_trans(base = 10)$transform(PC1),
         PC2 = pseudo_log_trans(base = 10)$transform(PC2))

# Now plot the transformed PCA data without any additional axis scaling

# Save the PCA plot without SVs to a PNG file
output_file <- file.path(output_dir, paste0("PCA_plot_", tools::file_path_sans_ext(basename(design_file)), "_organ.png"))
png(output_file, width = 800, height = 600)
ggplot(pcaData, aes(PC1, PC2, color = organ, shape = patient)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle(paste0("PCA after SVA Correction with ", n.sv," SVs ")) +
  coord_fixed()

dev.off()

#rownames(correctedse) <- rownames(counts)
print("correctedse")
print(nrow(correctedse))
print(colnames(correctedse))
logcts <- assay(correctedse)

# Select top variable genes and prepare for heatmap

rv <- rowVars(logcts)
last_row <- nrow(correctedse)
o <- order(rv, decreasing = TRUE)[1:last_row]
qlfvals <- assay(correctedse)

dists <- Dist(t(qlfvals[o, ]), method = 'pearson')
cols <- colorRampPalette(c("blue",  "white",  "red"))(255)
distmat <- as.matrix(dists)
df <- data.frame(condition = colnames(qlfvals[o, ]), row.names = colnames(distmat))

# Convert to data frame for heatmap
qlfvals_df <- as.data.frame(qlfvals[o, ])
#qlfvals_df <- qlfvals_df[, !colnames(qlfvals_df) %in% c('R18LIV_1', 'R18LIV_2')]
#qlfvals_df <- apply(qlfvals_df, 2, function(x) pseudo_log_trans(base = 10)$transform(x))

# Create a dataframe that labels each sample as Primary or Metastasis
sample_annotation <- data.frame(
  Condition = colData_combined$met  # Assuming 'met' contains 'Primary' and 'Metastasis'
)

# Set row names of the annotation to be the sample names, which must match the column names of your data
rownames(sample_annotation) <- colnames(qlfvals_df)

# Define colors for the annotation (e.g., Primary in blue, Metastasis in red)
annotation_colors <- list(
  Condition = c(Primary = "blue", Metastasis = "red")
)


# Save heatmap
output_file <- file.path(output_dir, paste0("Heatmap_", tools::file_path_sans_ext(basename(design_file)), ".png"))
png(output_file, width = 800, height = 600)
pheatmap(na.omit(qlfvals_df), scale = 'row',
            color = cols,
            cluster_rows = TRUE,
            show_rownames = FALSE, 
            show_colnames = TRUE,
            annotation_col = sample_annotation,   # Add annotation for Primary/Metastasis
            annotation_colors = annotation_colors)
dev.off()





# Create a  volcano plot and save it to a file
output_file <- file.path(output_dir, paste0("Corrected_Volcano_plot_", tools::file_path_sans_ext(basename(design_file)), ".png"))
png(output_file, width = 800, height = 600)
EnhancedVolcano(toptable = corrected_results,
                x = "log2FoldChange",
                y = "padj",
                lab = corrected_results$gene,
                xlim = c(-7, 7),
                ylim = c(0, 7.5),
                pCutoff = 0.05,  # Adjusted p-value cutoff
                FCcutoff = 1.5)
dev.off()




# Specify output file for the MA plot
output_file <- file.path(output_dir, paste0("Corrected_GGMA_plot_", tools::file_path_sans_ext(basename(design_file)), ".png"))
# Create the MA plot using ggmaplot from ggpubr
png(output_file, width = 800, height = 600)

ggmaplot(
  data = corrected_results,           # Input data
  fdr = 0.05,                         # False discovery rate
  fc = 1.5,                           # Fold change threshold
  genenames = corrected_results$gene_name,  # Gene names for labeling
  palette = c("#B31B21", "#1465AC", "darkgray"),  # Custom color palette
  top = 15,                           # Number of top genes to label
  main = "MA-Plot: Primary vs Metastasis",  # Plot title
  xlab = "Log2 mean expression",      # X-axis label
  ylab = "Log2 fold change",          # Y-axis label
  size = 0.5, 
  legend = "top", 
  font.label = c("bold", 11),
  font.legend = "bold",
  font.main = "bold",
  ggtheme = ggplot2::theme_minimal()
)

dev.off()
