################################################################################
# Frances Heredia
# Franco Lab
# Description: This script implements the TFSEE pipeline to compute TFSEE scores
# from footprinting data and bulk RNA-seq data across multiple tissues.		
# It generates heatmaps and deltaZ plots to visualize TF activity differences.
################################################################################



library(data.table) 
library(stringr)

rna_counts <- read.csv("~/RNA_seq/counts_PE.csv", header = TRUE, sep = ",")
rna_counts <- rna_counts[, !colnames(rna_counts) %in% c("Chr", "Start", "End", "Strand", "Length")]
# Fix for : Error: unexpected numeric constant in "RNA_final$49"  when running RNA_final$49BR_3 <- NULL
rna_counts$"R49BR_3" <- NULL
colnames(rna_counts) <- sub("^R", "", colnames(rna_counts))

nonid_cols <- setdiff(colnames(rna_counts), c("gene_name", "Geneid"))
RNA_samples <- as.matrix(sapply(rna_counts[, nonid_cols, drop=FALSE], as.numeric))
# set rownames to gene symbols (you could use Geneid if you prefer)
rownames(RNA_samples) <- rna_counts$gene_name
dup_indices <- which(duplicated(rownames(RNA_samples)) | duplicated(rownames(RNA_samples), fromLast = TRUE))
RNA_samples <- RNA_samples[-dup_indices, ]

sample_info_atac <- fread("sample_info") 
sample_info_atac <- as.data.frame(sample_info_atac)
rownames(sample_info_atac) <- sample_info_atac$Samples
sample_info_rna <- sample_info_atac 
group_col <- "Organ"

# Remove everything after the first underscore ("_")

organs <- c("Breast", "Liver" , "Lung")
for ( organ in organs) {
	colnames(fp_list[[organ]]) <- sub("_.*", "", colnames(fp_list[[organ]]))
}


human_tfs <- readLines("human_tfs.txt")

# Find overlapping TF between RNA and Footprints
common_tfs <- intersect(colnames(fp_list[["Liver"]]), rownames(RNA_samples))
message("TFs in common (motif matrix vs RNA genes): ", length(common_tfs),
        " / ", ncol(fp_list[["Liver"]]), " (T) ; ", nrow(RNA_samples), " (RNA)")
# TFs in common (motif matrix vs RNA genes): 643 / 2346 (T) ; 60981 (RNA)       

RNA_final <- RNA_samples[common_tfs, , drop = FALSE]


library(edgeR)


## ---------- Helpers ----------
aggregate_by_group <- function(mat, sample_info, group_col) {
  if (is.null(group_col)) {
    return(mat)  # no aggregation; keep each sample as its own condition
  }
  stopifnot(group_col %in% colnames(sample_info))
  groups <- split(seq_len(ncol(mat)), sample_info[[group_col]])
  out <- sapply(groups, function(idx) rowSums(as.matrix(mat[, idx, drop = FALSE])))
  if (is.null(dim(out))) out <- matrix(out, nrow = nrow(mat), dimnames = list(rownames(mat), names(groups)))
  out
}

## ---------- 2) Build R from BULK RNA ----------
## Aggregate replicates (if any), normalize to log2 CPM, filter to TFs, and z-score per TF across conditions.
build_R_from_bulk_rna <- function(rna_counts, sample_info_rna, group_col = NULL, tf_genes = NULL) {
  stopifnot(all(colnames(rna_counts) %in% rownames(sample_info_rna)))
  sample_info_rna <- sample_info_rna[colnames(rna_counts), , drop = FALSE]
  pb <- aggregate_by_group(rna_counts, sample_info_rna, group_col) # [genes x conditions]

  dge <- DGEList(counts = pb)
  dge <- calcNormFactors(dge, method = "TMM")
  logCPM <- cpm(dge, log = TRUE, prior.count = 1)  # [genes x conditions]
  # z-score per TF (row) across conditions (columns)
  R_raw <- logCPM[rownames(logCPM) %in% human_tfs, , drop = FALSE]
  R_scaled <- zscore_rows(R_raw)  # [TF x conditions]
  R_scaled
}


R_scaled <- build_R_from_bulk_rna(RNA_final, sample_info_rna, group_col = "Organ", tf_genes = human_tfs)

# -------------- 3) Compute Z from Footprints and R --------------
compute_Z_from_footprints <- function(fp_list, R_scaled) {
  # Ensure organ names in fp_list match R_scaled column names
  orgs <- intersect(names(fp_list), colnames(R_scaled))
  fp_list <- fp_list[orgs]
  R_scaled <- R_scaled[, orgs, drop = FALSE]
  
  # Identify common TFs across footprint and RNA
  common_tfs <- Reduce(intersect, lapply(fp_list, colnames))
  common_tfs <- intersect(common_tfs, rownames(R_scaled))
  
  # Filter RNA and footprints
  R_scaled <- R_scaled[common_tfs, , drop = FALSE]
  fp_list  <- lapply(fp_list, function(x) x[, common_tfs, drop = FALSE])
  
  message("TFs in common across organs: ", length(common_tfs))
  
  # Compute per-organ aggregate footprint signal
  Z_list <- lapply(orgs, function(org) {
    fp_mat <- fp_list[[org]]   # [peaks x TF]
    expr_vec <- R_scaled[, org]  # [TF]
    # Elementwise multiply each TF’s footprint vector by its expression
    weighted_fp <- sweep(fp_mat, 2, expr_vec, "*")
    # Aggregate across peaks → one Z-score per TF per organ
    z_scores <- colMeans(weighted_fp, na.rm = TRUE)
    z_scores
  })
  
  Z <- do.call(rbind, Z_list)
  rownames(Z) <- orgs
  return(Z)
}

Z <- compute_Z_from_footprints(fp_list, R_scaled)

# Scale -1 to 1 for comparability
Z_min <- min(Z)
Z_max <- max(Z)
Z_scaled <- 2 * (Z - Z_min) / (Z_max - Z_min) - 1

saveRDS(Z_scaled, "Z_scaled_footprinting_organs.rds")

library(pheatmap)

# Compute variance for Breast vs Liver
Z_BL <- Z_scaled[c("Breast", "Liver"), , drop = FALSE]
var_BL <- apply(Z_BL, 2, var)

# Top 75 most variable TFs between Breast and Liver
top_BL <- names(sort(var_BL, decreasing = TRUE))[1:75]

# Compute variance for Breast vs Lung
Z_BLu <- Z_scaled[c("Breast", "Lung"), , drop = FALSE]
var_BLu <- apply(Z_BLu, 2, var)

# Top 75 most variable TFs between Breast and Lung
top_BLu <- names(sort(var_BLu, decreasing = TRUE))[1:75]

# Combine and remove duplicates
top_tfs <- unique(c(top_BL, top_BLu))

# Subset Z_scaled to these TFs
Z_top <- Z_scaled[, top_tfs]

# Check
dim(Z_top)
head(Z_top[, 1:10])

# Plot the heatmap

png("Footprinting_heatmap_Z_scaled_top_var_pairs.png", width=4000, height=1500, res=300)
pheatmap(Z_top,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         scale = "none",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Top variable TFs between Breast–Liver and Breast–Lung",
		 show_rownames=TRUE,
		 show_colnames=TRUE,
		 fontsize=6)
dev.off()

# Compute deltaZ = Liver - Breast
deltaZ <- Z_top["Liver", ] - Z_top["Breast", ]

# Sort deltaZ decreasingly
deltaZ_sorted <- sort(deltaZ, decreasing = FALSE)
genes <- colnames(Z_top)
# Get the gene names
genes_sorted <- names(deltaZ_sorted)

# Create x-axis: 1 to 50
x <- 1:length(deltaZ_sorted)

png("deltaZ__Liver_vs_Breast_no_label.png", width=600, height=600)
plot(x, deltaZ_sorted,
     type = "p",      # dots only
     pch = 21,        # circle with fill color
     col = "#1961c5ff",          # outline color
     bg = "#c3e8f8",             # fill color inside the dots
     cex = 2,                     # bigger size
     xlab = "Transcription Factors in Liver ranked",
     ylab = expression(Delta * " scaled TFSEE score"),
     main = "DeltaZ: Liver vs Breast",
     xaxt = "n")

# Add a loess trendline
trend <- loess(deltaZ_sorted ~ x, span = 0.3)
x_pred <- seq(min(x), max(x), length.out = 200)
y_pred <- predict(trend, newdata = data.frame(x = x_pred))
lines(x_pred, y_pred, col = "darkorange", lwd = 2)

# Add dotted lines at ± 0.25
abline(h = c(-0.25, 0.25), col = "gray", lty = 2)

# Add x-axis numbers
axis(1, at = x)

dev.off()

png("deltaZ_Liver_vs_Breast_shaded.png", width=600, height=600)

# Set up the base plot (no points yet)
plot(x, deltaZ_sorted,
     type = "n",   # empty plot first (to add shading)
     xlab = "Transcription Factors in Liver ranked",
     ylab = expression(Delta * " scaled TFSEE score"),
     main = "DeltaZ: Liver vs Breast",
     xaxt = "n",
     ylim = range(deltaZ_sorted))

# Add shaded regions above and below ±0.25
usr <- par("usr")  # get plot limits: c(xmin, xmax, ymin, ymax)
rect(usr[1], 0.25, usr[2], usr[4], col = "gray95", border = NA)  # top region
rect(usr[1], usr[3], usr[2], -0.25, col = "gray95", border = NA) # bottom region

# Now add points
points(x, deltaZ_sorted,
       pch = 21,
       col = "#1961c5ff",
       bg = "#c3e8f8",
       cex = 2)

# Add loess trendline
trend <- loess(deltaZ_sorted ~ x, span = 0.3)
x_pred <- seq(min(x), max(x), length.out = 200)
y_pred <- predict(trend, newdata = data.frame(x = x_pred))
lines(x_pred, y_pred, col = "darkorange", lwd = 2)

# Add dotted lines at ± 0.25
abline(h = c(-0.25, 0.25), col = "gray50", lty = 2)

# Add x-axis
axis(1, at = x)

dev.off()

# Compute deltaZ = Lung - Breast
deltaZ <- Z_top["Lung", ] - Z_top["Breast", ]

# Sort deltaZ decreasingly
deltaZ_sorted <- sort(deltaZ, decreasing = FALSE)
genes <- colnames(Z_top)
# Get the gene names
genes_sorted <- names(deltaZ_sorted)

# Create x-axis: 1 to 50
x <- 1:length(deltaZ_sorted)

png("deltaZ_Lung_vs_Breast_shaded.png", width=600, height=600)
# Set up the base plot (no points yet)
plot(x, deltaZ_sorted,
     type = "n",   # empty plot first (to add shading)
     xlab = "Transcription Factors in Lung ranked",
     ylab = expression(Delta * " scaled TFSEE score"),
     main = "DeltaZ: Lung vs Breast",
     xaxt = "n",
     ylim = range(deltaZ_sorted))

# Add shaded regions above and below ±0.25
usr <- par("usr")  # get plot limits: c(xmin, xmax, ymin, ymax)
rect(usr[1], 0.25, usr[2], usr[4], col = "gray95", border = NA)  # top region
rect(usr[1], usr[3], usr[2], -0.25, col = "gray95", border = NA) # bottom region

# Now add points
points(x, deltaZ_sorted,
       pch = 21,
       col = "#006633",
       bg = "#cbf2d1",
       cex = 2)

# Add loess trendline
trend <- loess(deltaZ_sorted ~ x, span = 0.3)
x_pred <- seq(min(x), max(x), length.out = 200)
y_pred <- predict(trend, newdata = data.frame(x = x_pred))
lines(x_pred, y_pred, col = "darkorange", lwd = 2)

# Add dotted lines at ± 0.25
abline(h = c(-0.25, 0.25), col = "gray50", lty = 2)

# Add x-axis
axis(1, at = x)

dev.off()
