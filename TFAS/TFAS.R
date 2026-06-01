##############################################################################################################
# Frances Heredia
# Franco Lab
# Description: This script performs the following tasks
# 1. Reads TF footprints data
# 2. Intersects with ReMap TF binding data to identify which TFs have footprints in your ATAC peaks
# 3. Computes a per-sample TF activity score by weighting footprint presence by RNA expression of the TF	
##############################################################################################################


library(data.table)
library(GenomicRanges)
library(edgeR)

# ---------- USER PATHS / SETTINGS ----------
fp_list_rds <- "footprinting_matrices_list.rds"   # list: names = organs, each matrix peaks x TF (numeric)
remap_tf_beds_dir <- "beds_by_TF"                  # per-TF merged ReMap BEDs (TF_merged.bed)
sample_info_file <- "sample_info"                  # sample info table (must have column 'Samples' and group_col)
rna_counts_file <- "/Metastasis_Project/RNA_seq/counts_PE.csv"
group_col <- "Organ"                               # grouping column name in sample_info
human_tfs_file <- "tf_list.txt"                    # list of TF gene symbols (one per line)
out_Z_scaled_rds <- "Z_scaled_remap_footprints.rds"

# ---------- 0) Load inputs ----------
fp_list <- readRDS(fp_list_rds)   # list of matrices: peaks x TF (rownames "chr:start-end")
stopifnot(length(fp_list) > 0)

sample_info <- fread(sample_info_file)
sample_info <- as.data.frame(sample_info)
rownames(sample_info) <- sample_info$Samples

rna_counts <- fread(rna_counts_file)
# remove unwanted cols if present
rna_counts <- as.data.frame(rna_counts)
rna_counts <- rna_counts[, !colnames(rna_counts) %in% c("Chr","Start","End","Strand","Length")]
# fix a problematic column if present (adjust if not needed)
if("R49BR_3" %in% colnames(rna_counts)) rna_counts$R49BR_3 <- NULL
colnames(rna_counts) <- sub("^R", "", colnames(rna_counts))

human_tfs <- if(file.exists(human_tfs_file)) readLines(human_tfs_file) else NULL

# ---------- 1) Build unified peak set from footprint matrices ----------
all_peaks <- unique(unlist(lapply(fp_list, rownames)))
df <- do.call(rbind, strsplit(all_peaks, "[:-]"))
peaks_df <- data.frame(chr = df[,1], start = as.integer(df[,2]), end = as.integer(df[,3]), stringsAsFactors = FALSE)

gr_all <- GRanges(seqnames = peaks_df$chr, ranges = IRanges(start = peaks_df$start, end = peaks_df$end))
seqlevelsStyle(gr_all) <- "UCSC"
gr_all <- keepSeqlevels(gr_all, paste0("chr", c(1:22,"X","Y")), pruning.mode="coarse")
gr_peaks <- reduce(gr_all, ignore.strand = TRUE)   # merged footprint peak set

# ---------- 2) Build ReMap presence matrix aligned to gr_peaks ----------
# For each TF BED in remap_tf_beds_dir, count overlaps to gr_peaks (0/1)
tf_files <- list.files(remap_tf_beds_dir, pattern = "_merged\\.bed$", full.names = TRUE)
TF_names <- sub("_merged.bed$", "", basename(tf_files))

T_on_peaks <- matrix(0L, nrow = length(gr_peaks), ncol = length(TF_names),
                     dimnames = list(paste0(seqnames(gr_peaks), ":", start(gr_peaks), "-", end(gr_peaks)), TF_names))

wanted <- paste0("chr", c(1:22, "X", "Y"))

for(i in seq_along(tf_files)) {
  remap_gr <- rtracklayer::import(tf_files[i])
  seqlevelsStyle(remap_gr) <- "UCSC"
  common_sl <- intersect(seqlevels(remap_gr), wanted)
  if(length(common_sl) == 0) next
  remap_gr <- keepSeqlevels(remap_gr, common_sl, pruning.mode = "coarse")
  # if any intervals are outside ranges (NA starts/ends) skip
  ov <- countOverlaps(gr_peaks, remap_gr, ignore.strand = TRUE)
  T_on_peaks[, i] <- as.integer(ov > 0)
}


# ---------- 3) Prepare footprint matrices: ensure same peaks order as gr_peaks ----------

# For each organ, map its footprint rows onto gr_peaks and compute aggregated footprint per gr_peak.
# If an organ's original peak row overlaps multiple gr_peaks, assign by overlap and average values into gr_peak.
map_aggregate_fp <- function(fp_mat, gr_orig, gr_target) {
  # fp_mat: peaks x TF (rownames are "chr:start-end")
  # gr_orig: GRanges for fp_mat rows
  ov <- findOverlaps(gr_orig, gr_target, ignore.strand = TRUE)
  q <- queryHits(ov); s <- subjectHits(ov)
  # for each target index, aggregate (mean) values across overlapping original peaks
  TFs <- colnames(fp_mat)
  out <- matrix(0, nrow = length(gr_target), ncol = length(TFs), dimnames = list(names(gr_target), TFs))
  if(length(ov)==0) return(out)
  # use data.table for speed
  dt_idx <- data.table(q = q, s = s)
  agg_list <- dt_idx[, .(rows = list(q)), by = s]
  for(i in seq_len(nrow(agg_list))) {
    sidx <- agg_list$s[i]
    ridxs <- agg_list$rows[[i]]
    out[sidx, ] <- colMeans(fp_mat[ridxs, , drop = FALSE], na.rm = TRUE)
  }
  out
}

# parse original fp_list rownames to GRanges once
gr_fp_list <- lapply(fp_list, function(m) {
  df <- do.call(rbind, strsplit(rownames(m), "[:-]"))
  GRanges(seqnames = df[,1], ranges = IRanges(start = as.integer(df[,2]), end = as.integer(df[,3])), names = rownames(m))
})
# build aggregated footprint matrices aligned to gr_peaks
fp_on_peaks_list <- list()
for(org in names(fp_list)) {
  fp_mat <- fp_list[[org]]
  gr_orig <- gr_fp_list[[org]]
  seqlevelsStyle(gr_orig) <- "UCSC"
  gr_orig <- keepSeqlevels(gr_orig, paste0("chr", c(1:22,"X","Y")), pruning.mode="coarse")
  fp_on_peaks_list[[org]] <- map_aggregate_fp(fp_mat, gr_orig, gr_peaks)
}

# ---------- 4) Build R_scaled from RNA (TFs x organs) ----------
# Prepare RNA matrix: rows = genes, cols = samples
rna_mat <- as.matrix(sapply(rna_counts[, setdiff(colnames(rna_counts), c("gene_name","Geneid")), drop = FALSE], as.numeric))
rownames(rna_mat) <- rna_counts$gene_name
# Remove duplicated gene names
dup <- which(duplicated(rownames(rna_mat)) | duplicated(rownames(rna_mat), fromLast = TRUE))
if(length(dup) > 0) rna_mat <- rna_mat[-dup, , drop = FALSE]

# aggregate by organ (sum across samples per organ)
library(edgeR)

aggregate_by_group <- function(mat, sample_info, group_col) {
  if (is.null(group_col)) return(mat)
  stopifnot(group_col %in% colnames(sample_info))
  sample_info <- sample_info[colnames(mat), , drop = FALSE]
  groups <- split(seq_len(ncol(mat)), sample_info[[group_col]])
  out <- sapply(groups, function(idx) rowSums(as.matrix(mat[, idx, drop = FALSE])))
  if (is.null(dim(out))) out <- matrix(out, nrow = nrow(mat), dimnames = list(rownames(mat), names(groups)))
  out
}

zscore_rows <- function(m) {
  z <- t(scale(t(m)))
  z[is.na(z)] <- 0
  z
}

# then your code
stopifnot(all(colnames(rna_mat) %in% rownames(sample_info)))
si <- sample_info[colnames(rna_mat), , drop = FALSE]
pb <- aggregate_by_group(rna_mat, si, group_col)
dge <- DGEList(counts = pb); dge <- calcNormFactors(dge, method = "TMM")
logCPM <- cpm(dge, log = TRUE, prior.count = 1)


# Restrict to TF genes and z-scale rows
if(!is.null(human_tfs)) tf_genes <- human_tfs else tf_genes <- rownames(logCPM)
R_raw <- logCPM[rownames(logCPM) %in% tf_genes, , drop = FALSE]
zscore_rows <- function(m) { z <- t(scale(t(m))); z[is.na(z)] <- 0; z }
R_scaled <- zscore_rows(R_raw)   # rows = TF, cols = organs (group_col)

# ---------- 5) Compute Z: for each organ, weight footprint signal by TF expression and aggregate across peaks ----------
orgs <- intersect(names(fp_on_peaks_list), colnames(R_scaled))
fp_on_peaks_list <- fp_on_peaks_list[orgs]
R_scaled <- R_scaled[, orgs, drop = FALSE]

# Identify TFs common to footprints-on-peaks and R_scaled
common_tfs <- Reduce(intersect, lapply(fp_on_peaks_list, colnames))
common_tfs <- intersect(common_tfs, rownames(R_scaled))
if(length(common_tfs) == 0) stop("No TFs in common between footprints and RNA")

# Filter footprints and R_scaled to common TFs
fp_on_peaks_list <- lapply(fp_on_peaks_list, function(m) m[, common_tfs, drop = FALSE])
R_scaled <- R_scaled[common_tfs, , drop = FALSE]

# Compute Z: organ × TF
Z <- matrix(NA_real_, nrow = length(orgs), ncol = length(common_tfs),
            dimnames = list(orgs, common_tfs))
for(org in orgs) {
  fp_mat <- fp_on_peaks_list[[org]]           # peaks x TF
  expr_vec <- R_scaled[, org]                 # TF
  weighted_fp <- sweep(fp_mat, 2, expr_vec, FUN = "*")
  Z[org, ] <- colMeans(weighted_fp, na.rm = TRUE)
}

# ---------- 6) Scale Z to -1..1 and save ----------
Z_min <- min(Z, na.rm = TRUE); Z_max <- max(Z, na.rm = TRUE)
Z_scaled <- 2 * (Z - Z_min) / (Z_max - Z_min) - 1
saveRDS(Z_scaled, out_Z_scaled_rds)

# print basic info
cat("Saved Z_scaled to", out_Z_scaled_rds, "\n")
cat("Dimensions (orgs x TFs):", dim(Z_scaled), "\n")

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

# 3️⃣ Combine and remove duplicates
top_tfs <- unique(c(top_BL, top_BLu))

# Subset Z_scaled to these TFs
Z_top <- Z_scaled[, top_tfs]

# Plot the heatmap
png("Footprinting_Remap_heatmap_Z_scaled_top_var_pairs.png", width=4000, height=1500, res=300)
pheatmap(Z_scaled,
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
 deltaZ <- Z_scaled["Liver", ] - Z_scaled["Breast", ]
 
 # Sort deltaZ decreasingly
 deltaZ_sorted <- sort(deltaZ, decreasing = FALSE)
 genes <- colnames(Z_scaled)
 # Get the gene names
 genes_sorted <- names(deltaZ_sorted)
 
 # Create x-axis: 1 to 50
 x <- 1:length(deltaZ_sorted)

  png("deltaZ_Remap_Liver_vs_Breast_shaded.png", width=600, height=600)
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
 deltaZ <- Z_scaled["Lung", ] - Z_scaled["Breast", ]
 
 # Sort deltaZ decreasingly
 deltaZ_sorted <- sort(deltaZ, decreasing = FALSE)
 genes <- colnames(Z_scaled)
 # Get the gene names
 genes_sorted <- names(deltaZ_sorted)
 
 # Create x-axis: 1 to 50
 x <- 1:length(deltaZ_sorted)
 
 png("deltaZ_Remap_Lung_vs_Breast_shaded.png", width=600, height=600)
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



