################################################################################
# Frances Heredia
# Franco Lab
# Description: This script performs the following tasks  
# 	  1) Differential expression analysis per organ (Breast, Liver, Lung)
# 	  2) PCA plots per organ
# 	  3) Volcano plots per organ
# 	  4) Heatmaps of differentially expressed genes per organ
################################################################################

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)  # Load the limma package for batch effect correction
library(pheatmap)  # Load the pheatmap package

# Step 1: Read and prepare data
full_data <- read_csv("FULL_gene_data_PE_MET_GTEX_one_to_one.csv")
full_data <- as.data.frame(full_data)

# Define the organs to analyze
organs_to_analyze <- c("Breast", "Liver", "Lung")

# Loop through each organ
for (organ_name in organs_to_analyze) {
    

    cat("Processing organ:", organ_name, "\n")

    organ_columns <- colnames(full_data)[grepl(organ_name, colnames(full_data))]

    tumor_columns <- organ_columns[grepl("Tumor", organ_columns)]
    
    normal_columns <- organ_columns[grepl("Normal", organ_columns)]
    print(tumor_columns)
    print(normal_columns)
        # Select only normal and tumor columns for PCA
    pca_data <- full_data %>%
    dplyr::select(gene_name, Geneid, all_of(c(normal_columns, tumor_columns))) %>% drop_na()

    count_data <- as.matrix(pca_data[, -(1:2)])  # remove gene name columns
    filtered_data <- pca_data[, c("gene_name", "Geneid")]

    # Normalize and log-transform
    normalized_counts <- t(apply(count_data, 1, function(x) x / sum(x) * median(rowSums(count_data, na.rm = TRUE))))
    log_counts <- log2(normalized_counts + 1)

    # Remove low variance rows
    row_variances <- apply(log_counts, 1, var)
    low_variance_rows <- row_variances < 1e-2 | is.na(row_variances)

    log_counts <- log_counts[!low_variance_rows, ]
    filtered_data <- filtered_data[!low_variance_rows, ]

    # Use both gene_name and Geneid for clarity (or just gene_name if you prefer)
    rownames(log_counts) <- paste0(filtered_data$gene_name, " (", filtered_data$Geneid, ")")


    sample_names <- colnames(log_counts)
    group <- ifelse(grepl("Normal", sample_names, ignore.case = TRUE), "Normal", 
                    ifelse(grepl("Tumor", sample_names, ignore.case = TRUE), "Tumor", NA))

    


    # Create design matrix
    design <- model.matrix(~ 0 + group ) 
    colnames(design) <- gsub("group", "", colnames(design))  
    
    if (nrow(design) != ncol(log_counts)) {
        stop("Error: Design matrix dimensions do not match log count data.")
    }

    # Fit the model and apply batch effect correction
    fit <- lmFit(log_counts, design)


    # Perform PCA on the corrected counts
    pca_result <- prcomp(t(log_counts), center = TRUE, scale. = TRUE)
    explained_variance <- summary(pca_result)$importance[2, ]
    variance_pc1 <- round(explained_variance[1] * 100, 2)
    variance_pc2 <- round(explained_variance[2] * 100, 2)

    short_sample_names <- sub("^(R\\d+).*", "\\1", sample_names)

    # Create a data frame for PCA results
    pca_df <- data.frame(Sample = short_sample_names, 
                         PC1 = pca_result$x[, 1], 
                         PC2 = pca_result$x[, 2], 
                         Group = group)

    

    # Create a second PCA plot (corrected)
    png(paste0("PE_counts/pca_plot_", organ_name, "_short-label.png"), width = 10, height = 8)
    ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
        geom_point() +
        geom_text(vjust = 1.5, hjust = 0.5, size = 3) +  
        labs(title = paste("PCA of Normal and Tumor Samples -", organ_name),
             x = paste("PC1 (", variance_pc1, "%)", sep = ""),
             y = paste("PC2 (", variance_pc2, "%)", sep = "")) +
        theme_minimal() +
        scale_color_manual(values = c("Normal" = "blue", "Tumor" = "red")) + 
        theme(legend.title = element_blank()) 
    dev.off()

    cat("Finished processing:PCA")

    # Adding Differential expression and volcano plots
    design <- model.matrix(~ 0 + group )  # Include organ in the model
    colnames(design) <- gsub("group", "", colnames(design))  # Clean up column names

    fit <- lmFit(log_counts, design)
    contrast_matrix <- makeContrasts(Tumor = Tumor - Normal, levels = design)
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)

    # Get the results for differential expression
    results <- topTable(fit2, adjust = "fdr", sort.by = "P", number = Inf)
    
    # Create Volcano Plot
    results$Significance <- -log10(results$P.Value)
    results$Threshold <- 0.05  # Set threshold for significance

    png(paste0("PE_counts/volcano_plot_", organ_name, ".png"), width = 8, height = 6)
    ggplot(results, aes(x = logFC, y = Significance)) +
        geom_point(aes(color = (P.Value < 0.05)), alpha = 0.5) +
        scale_color_manual(values = c("grey", "red")) +
        geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
        geom_vline(xintercept = c(-1, 1), linetype = "dotted") +
        labs(title = paste("Volcano Plot of Differential Expression -", organ_name),
             x = "Log Fold Change",
             y = "-Log10(p-value)") +
        theme_minimal() 
    dev.off()

    # Create a  volcano plot and save it to a file
    output_file <- file.path(output_dir, paste0("Enhanced_Volcano_plot_", organ_name, ".png"))
    png(output_file, width = 800, height = 600)
    EnhancedVolcano(toptable = results,
                    x = "log2FoldChange",
                    y = "padj",
                    lab = results$gene_name,
                    xlim = c(-7, 7),
                    ylim = c(0, 7.5),
                    pCutoff = 0.05,  # Adjusted p-value cutoff
                    FCcutoff = 1.5)
    dev.off()


    cat("Finished processing:Volcano")

	results <- results[ , -1]

# Extract only the Ensembl ID from the ID column
	resultsf$ID <- sub(".*\\((ENSG[0-9]+)\\)", "\\1", results$ID)

    write.csv(results, paste0("PE_counts/differential_expression_", organ_name, "_results_Met_Gtex_corrected.csv"), row.names = TRUE)

    results$log2FoldChange <- results$logFC  # Rename or ensure the column for log2FoldChange
    results$padj <- results$adj.P.Val  

    # Create a new data frame with necessary columns
    corrected_results <- data.frame(
        log2FoldChange = results$log2FoldChange, 
        padj = results$padj, 
        gene = rownames(results)  # Use row names as gene identifiers
    )

    # Set NA for adjusted p-values when thresholding, optional
    corrected_results[is.na(corrected_results$padj), "padj"] <- 1

    library(EnhancedVolcano)

    output_file <-  paste0("PE_counts/Corrected_Volcano_plot_", organ_name, ".png")
    png(output_file, width = 800, height = 600)
    EnhancedVolcano(toptable = corrected_results,
                x = "log2FoldChange",
                y = "padj",
                lab = corrected_results$gene,
                xlim = c(-4, 7),
                ylim = c(0, 40),
                pCutoff = 0.05,  # Adjusted p-value cutoff
                FCcutoff = 1.5)
    dev.off()
    cat("Finished processing:Volcano 2")
    threshold <- 0.05

    # Create a heatmap of significant results
    significant_results <- results[results$P.Value < threshold, ]
    top_genes <- significant_results %>% 
        arrange(desc(abs(logFC)))

    # Get gene names
    selected_genes <- rownames(top_genes)

    # Subset original log counts data for selected genes
    heatmap_data <- log_counts[selected_genes, ]

    # Normalize for heatmap
    heatmap_data <- t(scale(t(heatmap_data)))

    custom_colors <- colorRampPalette(c("blue", "white", "red"))(50)

    png(paste0("PE_counts/heatmap_", organ_name, "_V2.png"), width = 1400, height = 800)
    pheatmap(heatmap_data,
             show_rownames = FALSE,
             show_colnames = TRUE,
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             clustering_method = "complete",
             main = paste("Heatmap of Top Differentially Expressed Genes -", organ_name),
             color = custom_colors)
    dev.off()
    
    cat("Finished processing organ:", organ_name, "\n")
}
