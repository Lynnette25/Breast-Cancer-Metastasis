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
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)

# Output dir
output_dir <- "PE_counts"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Step 1: Read data
full_data <- read_csv("FULL_gene_data_PE_MET_GTEX_one_to_one.csv", show_col_types = FALSE)
full_data <- as.data.frame(full_data)

# Clean up gene IDs (remove version numbers if present)
full_data$Geneid <- sub("\\..*$", "", full_data$Geneid)

# Ensure gene_name exists
if (!"gene_name" %in% colnames(full_data)) {
    full_data$gene_name <- full_data$Geneid
}

# Define organs to analyze
organs_to_analyze <- c("Breast", "Liver", "Lung")

for (organ_name in organs_to_analyze) {
    cat("Processing organ:", organ_name, "\n")

    # Identify columns for this organ
    organ_columns <- colnames(full_data)[grepl(organ_name, colnames(full_data))]
    tumor_columns <- organ_columns[grepl("Tumor", organ_columns, ignore.case = TRUE)]
    normal_columns <- organ_columns[grepl("Normal", organ_columns, ignore.case = TRUE)]

    if (length(tumor_columns) == 0 || length(normal_columns) == 0) {
        warning(paste0("No tumor or normal columns found for ", organ_name, " — skipping."))
        next
    }

    # Select gene id/name + the organ count columns
    sub_df <- full_data %>%
        dplyr::select(gene_name, Geneid, all_of(c(normal_columns, tumor_columns))) %>%
        drop_na()

    # Build count matrix
    count_data <- as.matrix(sub_df[, -(1:2)])  # remove gene_name & Geneid
    # Ensure integer counts (DESeq2 requires integers)
    storage.mode(count_data) <- "double"
    count_data <- round(count_data)            # round just in case
    storage.mode(count_data) <- "integer"

    # Set rownames as Geneid (unique) — keep gene_name separately for labels
    rownames(count_data) <- sub_df$Geneid

    # Build colData (sample metadata)
    sample_names <- colnames(count_data)
    condition <- ifelse(grepl("Normal", sample_names, ignore.case = TRUE), "Normal",
                        ifelse(grepl("Tumor", sample_names, ignore.case = TRUE), "Tumor", NA))
    # Source: assume Normal == GTEx, Tumor == Metastasis
    source <- ifelse(condition == "Normal", "GTEx",
                     ifelse(condition == "Tumor", "Metastasis", NA))

    coldata <- data.frame(
        sample = sample_names,
        condition = factor(condition, levels = c("Normal", "Tumor")),
        source = factor(source)
    )
    rownames(coldata) <- coldata$sample

    # Filter genes with very low counts across all samples (pre-filter)
    keep <- rowSums(count_data >= 10) >= 2   # keep genes with >=10 counts in >=2 samples
    count_data_filt <- count_data[keep, ]
    if (nrow(count_data_filt) < 50) {
        warning("Very few genes pass the filter for ", organ_name)
    }

    # Create DESeq2 dataset with batch/source correction in design
    dds <- DESeqDataSetFromMatrix(countData = count_data_filt,
                                  colData = coldata,
                                  design = ~ source + condition)

    # Relevel condition to ensure results use Tumor vs Normal as desired (Tumor as second level)
    dds$condition <- relevel(dds$condition, ref = "Normal")

    # Run DESeq
    dds <- DESeq(dds, quiet = TRUE)

    # Extract results for Tumor vs Normal (adjusted for source)
    res <- results(dds, contrast = c("condition", "Tumor", "Normal"), alpha = 0.05)
    res_df <- as.data.frame(res)
    res_df$Geneid <- rownames(res_df)

    # Add gene_name (if available)
    gene_map <- sub_df %>% dplyr::select(Geneid, gene_name) %>% distinct()
    res_df <- left_join(res_df, gene_map, by = "Geneid")

    # Save full results
    write.csv(res_df, file.path(output_dir, paste0("differential_expression_", organ_name, "_DESeq2_results.csv")), row.names = FALSE)

    # Variance stabilizing transform for PCA and heatmaps
    vst_dds <- vst(dds, blind = FALSE)
    vst_mat <- assay(vst_dds)

    # PCA
    pca_res <- prcomp(t(vst_mat), center = TRUE, scale. = TRUE)
    explained_var <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
    pc1_pct <- round(explained_var[1] * 100, 2)
    pc2_pct <- round(explained_var[2] * 100, 2)

    pca_df <- data.frame(Sample = rownames(coldata),
                         PC1 = pca_res$x[, 1],
                         PC2 = pca_res$x[, 2],
                         Condition = coldata$condition,
                         Source = coldata$source,
                         Short = sub("^(R\\d+).*", "\\1", rownames(coldata)))

    png(file.path(output_dir, paste0("pca_plot_", organ_name, "_short-label.png")), width = 1000, height = 800)
    print(
        ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition, shape = Source, label = Short)) +
            geom_point(size = 3) +
            geom_text(vjust = -1, hjust = 0.5, size = 3) +
            labs(title = paste("PCA of VST counts -", organ_name),
                 x = paste0("PC1 (", pc1_pct, "%)"),
                 y = paste0("PC2 (", pc2_pct, "%)")) +
            theme_minimal() +
            scale_color_manual(values = c("Normal" = "blue", "Tumor" = "red")) +
            theme(legend.title = element_blank())
    )
    dev.off()

    cat("PCA done for", organ_name, "\n")

    # Volcano plot with EnhancedVolcano
    # Ensure we have log2FoldChange and padj columns
    res_df$log2FoldChange <- res_df$log2FoldChange
    res_df$padj <- res_df$padj
    # Some packages expect column names exactly; also remove NA genes
    volcano_tbl <- res_df %>% dplyr::filter(!is.na(log2FoldChange))

    png(file.path(output_dir, paste0("Enhanced_Volcano_plot_", organ_name, ".png")), width = 900, height = 700)
    EnhancedVolcano(toptable = volcano_tbl,
                    lab = volcano_tbl$gene_name,
                    x = "log2FoldChange",
                    y = "padj",
                    pCutoff = 0.05,
                    FCcutoff = 1.5,
                    title = paste("Tumor vs Normal -", organ_name),
                    subtitle = "DESeq2 (adjusted for source)")
    dev.off()

    # Basic ggplot2 volcano (for raw p-values too)
    volcano_gg <- volcano_tbl %>%
        mutate(Significance = -log10(pvalue),
               Sig = ifelse(padj < 0.05, "FDR<0.05", "NS"))

    png(file.path(output_dir, paste0("volcano_plot_", organ_name, ".png")), width = 800, height = 600)
    print(
        ggplot(volcano_gg, aes(x = log2FoldChange, y = Significance, color = Sig)) +
            geom_point(alpha = 0.6) +
            geom_vline(xintercept = c(-1, 1), linetype = "dotted") +
            geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
            labs(title = paste("Volcano Plot -", organ_name),
                 x = "log2(Fold Change)",
                 y = "-log10(p-value)") +
            theme_minimal() +
            scale_color_manual(values = c("NS" = "grey", "FDR<0.05" = "red"))
    )
    dev.off()

    cat("Volcano done for", organ_name, "\n")

    # Heatmap of top significant genes (padj < 0.05). If none, choose top 50 by padj.
    sig_genes <- volcano_tbl %>% filter(!is.na(padj) & padj < 0.05) %>% arrange(padj)
    if (nrow(sig_genes) == 0) {
        # fallback: top 50 genes by smallest pvalue
        warning("No genes passed padj < 0.05 for ", organ_name, "- using top 50 by p-value for heatmap.")
        top_genes <- head(volcano_tbl %>% arrange(pvalue), 50)$Geneid
    } else {
        top_genes <- head(sig_genes$Geneid, 200)  # select up to 200 top genes for the heatmap
    }

    # Subset vst matrix and scale for heatmap
    ht_mat <- vst_mat[top_genes, , drop = FALSE]
    # If single gene or small number, pheatmap may require a matrix; ensure dimensions
    if (nrow(ht_mat) > 1) {
        ht_scaled <- t(scale(t(ht_mat)))
    } else {
        ht_scaled <- ht_mat
    }

    png(file.path(output_dir, paste0("heatmap_", organ_name, "_DESeq2.png")), width = 1400, height = 900)
    pheatmap(ht_scaled,
             show_rownames = FALSE,
             show_colnames = TRUE,
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             clustering_method = "complete",
             main = paste("Heatmap of DE genes -", organ_name),
             annotation_col = data.frame(Condition = coldata$condition, Source = coldata$source))
    dev.off()

    cat("Heatmap done for", organ_name, "\n")
    cat("Finished processing organ:", organ_name, "\n\n")
}
