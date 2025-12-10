################################################################################
# Frances Heredia
# Franco Lab
# Description: This script performs the following tasks  
# 1. Reads in ATAC-seq counts data
# 2. Normalizes the data using variance stabilizing transformation (VST)
# 3. Conducts PCA analysis to visualize sample clustering
# 4. Applies batch correction using limma's removeBatchEffect function
# 5. Visualizes PCA results before and after batch correction	
################################################################################

# Read the counts table
counts <- read.csv("Back_to_Peaks/All_samples_counts_peaks.csv", header = TRUE)

# First 3 columns are coordinates, rest are counts
peak_coords <- counts[,1:3]
count_matrix <- counts[,-c(1:3)]
rownames(count_matrix) <- paste(peak_coords$CHR, peak_coords$START, peak_coords$END, sep = "_")

head(count_matrix[,1:5])

coldata <- read.csv("ATACdata.csv", row.names = 1)
coldata <- coldata[colnames(count_matrix), , drop=FALSE]  # make sure order matches

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                              colData   = coldata,
                              design    = ~ Tissue)

dds <- estimateSizeFactors(dds)
vst_mat <- assay(vst(dds, blind=TRUE))

library(ggplot2)
library(matrixStats)

# ensure vst_mat is a numeric matrix
vst_mat <- as.matrix(vst_mat)

# compute row variances
rv <- rowVars(vst_mat, useNames = FALSE)

# select top variable peaks
top5000 <- order(rv, decreasing = TRUE)[1:5000]
mat <- vst_mat[top5000, ]


pca <- prcomp(t(mat), scale.=TRUE)

percentVar <- pca$sdev^2 / sum(pca$sdev^2)

pca_df <- data.frame(PC1 = pca$x[,1],
                     PC2 = pca$x[,2],
                     PC3 = pca$x[,3],
                     Tissue = coldata$Tissue,
                     Sample = rownames(pca$x))

# 2D PCA
png("Back_to_Peaks/Early_Breast_PCA_BreastLiverLung.png", width=600, height=400, res=150)
ggplot(pca_df, aes(x=PC1, y=PC2, color=Tissue, label=Sample)) +
  geom_point(size=3) +
  theme_minimal() +
  xlab(sprintf("PC1 (%.1f%%)", percentVar[1]*100)) +
  ylab(sprintf("PC2 (%.1f%%)", percentVar[2]*100))
dev.off()

dds <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                              colData   = coldata,
                              design    = ~ Status)

pca_df <- data.frame(PC1 = pca$x[,1],
                     PC2 = pca$x[,2],
                     PC3 = pca$x[,3],
                     Tissue = coldata$Status,
                     Sample = rownames(pca$x))

png("Back_to_Peaks/Early_Breast_PCA_EarlyPriMets.png", width=600, height=400, res=150)
ggplot(pca_df, aes(x=PC1, y=PC2, color=Tissue, label=Sample)) +
  geom_point(size=3) +
  theme_minimal() +
  xlab(sprintf("PC1 (%.1f%%)", percentVar[1]*100)) +
  ylab(sprintf("PC2 (%.1f%%)", percentVar[2]*100))
dev.off()	

# select top variable peaks
top1000 <- order(rv, decreasing = TRUE)[1:1000]
mat <- vst_mat[top1000, ]

# select top variable peaks
top500 <- order(rv, decreasing = TRUE)[1:500]
mat <- vst_mat[top500, ]


coldata <- read.csv("ATACdata.csv", row.names = 1)
coldata <- coldata[colnames(count_matrix), , drop=FALSE]  # make sure order matches

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                              colData   = coldata,
                              design    = ~ Status)

dds <- estimateSizeFactors(dds)

library(matrixStats)
library(limma)
library(ggplot2)

# vst transform
vst <- vst(dds)
vst_mat <- assay(vst)

# select top variable peaks
vst_mat <- as.matrix(vst_mat)
rv <- rowVars(vst_mat, useNames = FALSE)
top1000 <- order(rv, decreasing = TRUE)[1:1000]
mat <- vst_mat[top1000, ]

# batch correction (if you have a 'Batch' column)
mat_corrected <- removeBatchEffect(
  mat,
  batch = coldata$Batch1,
  batch2 = coldata$Batch2
)

# PCA
pca <- prcomp(t(mat_corrected))

df <- data.frame(pca$x, 
                 Tissue = vst$Status, 
                 Batch1 = vst$Batch2,
                 Batch2 = vst$Batch2,
                 Sample = colnames(vst_mat))
library(ggrepel)

png("Back_to_Peaks/Early_Breast_PCA_EarlyPriMets_Batchcorrection.png", width=600, height=400, res=150)
ggplot(df, aes(PC1, PC2, color = Tissue, label = Sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, show.legend = FALSE) +
  theme_bw() +
  labs(title = "ATAC-seq PCA (Batch-corrected)",
       x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)"))



library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                              colData   = coldata,
                              design    = ~ Status)

dds <- estimateSizeFactors(dds)

library(matrixStats)
library(limma)
library(ggplot2)

# vst transform
vst <- vst(dds)
vst_mat <- assay(vst)

# select top variable peaks
vst_mat <- as.matrix(vst_mat)
rv <- rowVars(vst_mat, useNames = FALSE)
top1000 <- order(rv, decreasing = TRUE)[1:1000]
mat <- vst_mat[top1000, ]

# batch correction (if you have a 'Batch' column)
mat_corrected <- removeBatchEffect(
  mat,
  batch = coldata$Batch1,
  batch2 = coldata$Batch2,
  batch3 = coldata$Batch3
)

# PCA
pca <- prcomp(t(mat_corrected))

df <- data.frame(pca$x, 
                 Tissue = vst$Status, 
                 Batch1 = vst$Batch1,
                 Batch2 = vst$Batch2,
                 Batch3 = vst$Batch3,
                 Sample = colnames(vst_mat))
library(ggrepel)

png("Back_to_Peaks/Early_Breast_PCA_EarlyPriMets_Batchcorrection3.png", width=600, height=400, res=150)
ggplot(df, aes(PC1, PC2, color = Tissue, label = Sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, show.legend = FALSE) +
  theme_bw() +
  labs(title = "ATAC-seq PCA (Batch-corrected)",
       x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)"))
dev.off()


```