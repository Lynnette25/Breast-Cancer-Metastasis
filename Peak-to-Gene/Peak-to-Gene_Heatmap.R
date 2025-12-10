################################################################################
# Frances Heredia
# Franco Lab
# Description: This script performs the following tasks  
# 1. Loads necessary libraries and functions for peak-to-gene analysis.
# 2. Defines a function to plot heatmaps of peak-to-gene linkages based on ATAC-seq and RNA-seq data.
# 3. Generates and saves a heatmap visualization of significant peak-to-gene linkages.		
################################################################################

library(cisDynet)
library(GenomicRanges)
library(IRanges)
library(dplyr)
library(progress)
library(BSDA)
library(tibble)
library(paletteer)

addAnnotation(gene_bed = "beds/50/hg38_gene_standard.bed",
              gtf = "beds/50/Modified_Homo_sapiens.GRCh38.108.gtf",
              genome_size = "beds/50/hg38.chrom.size")

plotP2GHeatmap <- function(p2g_res, cor_cutoff, atac_matrix, rna_matrix, cluster_N =4, raster = F,
                           palATAC = NA, palRNA= NA, save_path=NA, fig_width=12,fig_height=12){
  checkGeAnno()
  blueYellow <- c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")
  solarExtra <- c("5"='#3361A5', "7"='#248AF3', "1"='#14B3FF', "8"='#88CEEF', "9"='#C1D5DC', "4"='#EAD397', "3"='#FDB31A',"2"= '#E42A2A', "6"='#A31D1D')
  p2g <- readRDS(p2g_res)
  p2g <- p2g[p2g$correlations >= cor_cutoff, ]
  logfile(sprintf("Depending on the cutoff. The Peak2Gene links number: %s",dim(p2g)[1]))
  #logfile("This step may take a while, please wait...")
  rna <- read.table(rna_matrix, header = T, row.names = 1)
  p2g$idx_atac <- sprintf("ATAC_%s",1:nrow(p2g))
  p2g$idx_rna <- sprintf("RNA_%s",1:nrow(p2g))
  p2g_atac <- p2g[,c("Peak","idx_atac")]
  rownames(p2g_atac) <- p2g_atac$idx_atac
  atac <- read.table(atac_matrix, row.names = 1, header = T)
  atac_mat <- merge(p2g_atac, atac, by.x="Peak",by.y=0,all=F)
  rownames(atac_mat) <- atac_mat$idx_atac
  atac_mat <- atac_mat[,-c(1,2)]
  atac_scale <- rowZscores(as.matrix(atac_mat),limit=T)

  df1 <- rowZscores(as.matrix(atac_scale), limit = TRUE) %>% as.data.frame()
    logfile("Calculating ATAC matrix.")
    set.seed(2023)

    if (cluster_N > 1) {
    # Perform hierarchical clustering only if cluster_N > 1
    row_dend <- hclust(dist(df1))
    # Ensure k does not exceed the number of actual clusters
    n_clusters_possible <- length(unique(cutree(row_dend, k=cluster_N)))
    if (cluster_N > n_clusters_possible) {
        warning(sprintf("Requested cluster_N=%d exceeds the number of possible clusters=%d. Adjusting cluster_N to %d.", cluster_N, n_clusters_possible, n_clusters_possible))
        cluster_N <- n_clusters_possible
    } 

    group <- data.frame(C=cutree(row_dend, k=cluster_N))
    group$Cluster <- paste0("Cluster", group$C)
    group <- group[order(group$Cluster),]
    group$Cluster <- factor(group$Cluster, levels = sort(unique(group$Cluster)))
    mat <- df1[rownames(group),]
    group <- subset(group, select=-(C))
    } else {
    # For cluster_N = 1, assign all to a single cluster named "All"
    group <- data.frame(C=rep(1, nrow(df1)))
    group$Cluster <- "All"
    mat <- df1
    }
  annotation <- mat[match(rownames(group),rownames(as.data.frame(mat))),]
  mycol <- c(paletteer::paletteer_d("ggthemes::Tableau_10"), paletteer::paletteer_d("ggthemes::Tableau_20"))
  gcols <- setNames(as.character(mycol[1:cluster_N]),unique(group$Cluster))
  gcol <- list(Cluster=gcols)

  N <- cluster_N - 1
  if(!is.na(palATAC)){
    p3 <- ComplexHeatmap::pheatmap(as.matrix(annotation), show_rownames=F,cluster_row=F,cluster_col=F,border_color=NA,use_raster=raster,
                                   annotation_colors =gcol, color = colorRampPalette(palATAC)(256),
                                   annotation_row = group, gaps_row = cumsum(as.numeric(table(group$Cluster)))[1:N],
                                   main="ATAC-seq",name="ATAC Z-score")
  } else{
    p3 <- ComplexHeatmap::pheatmap(as.matrix(annotation), show_rownames=F,cluster_row=F,cluster_col=F,border_color=NA,use_raster=raster,
                                   annotation_colors =gcol, color = colorRampPalette(blueYellow)(256),
                                   annotation_row = group, gaps_row = cumsum(as.numeric(table(group$Cluster)))[1:N],
                                   main="ATAC-seq",name="ATAC Z-score")
  }
  # Set the RNA matrix gene order
  rownames(p2g) <- p2g$idx_atac
  logfile("Calculating RNA matrix.")
  p2gnew <- p2g[rownames(group),]
  p2g_rna <- p2gnew[,c("Gene","idx_rna")]
  rna_order <- merge(p2g_rna, rna, by.x="Gene",by.y=0,all.x =T)
  rownames(rna_order) <- rna_order$idx_rna
  rna_order <- rna_order[p2gnew$idx_rna,-c(1,2)]
  colnames(rna_order) <- gsub("^R", "", colnames(rna_order))
  colnames(annotation) <- gsub("^A", "", colnames(annotation))
  print(colnames(annotation))
  print(colnames(rna_order))
  rna_order <- rna_order[,colnames(annotation)]
  rna_scale <- rowZscores(as.matrix(rna_order),limit=T)
  if(!is.na(palRNA)){
    p4 <- ComplexHeatmap::pheatmap(as.matrix(rna_scale), show_rownames=F,cluster_row=F,cluster_col=F,border_color=NA,use_raster=raster,
                                   annotation_colors =gcol, color = colorRampPalette(palRNA)(256),
                                   annotation_row = group, gaps_row = cumsum(as.numeric(table(group$Cluster)))[1:N],
                                   main="RNA-seq", name="RNA Z-score")
  }else{
    p4 <- ComplexHeatmap::pheatmap(as.matrix(rna_scale), show_rownames=F,cluster_row=F,cluster_col=F,border_color=NA,use_raster=raster,
                                   annotation_colors =gcol, color = colorRampPalette(solarExtra)(256),
                                   annotation_row = group, gaps_row = cumsum(as.numeric(table(group$Cluster)))[1:N],
                                   main="RNA-seq", name="RNA Z-score")
  }
  pcor <- ComplexHeatmap::Heatmap(p2gnew$correlations, cluster_rows=F, name = "Correlations", col = paletteer::paletteer_d("RColorBrewer::Oranges")[1:7])
  ppvalue <- ComplexHeatmap::Heatmap(-log10(p2gnew$FDR), cluster_rows=F, name = "-log10(FDR)", col = paletteer::paletteer_d("RColorBrewer::Greens")[1:7])
  gene_type_col = setNames(as.character(paletteer_d("ggthemes::Classic_10")[1:length(unique(p2gnew$Type))]),unique(p2gnew$Type))
  ptype <- ComplexHeatmap::Heatmap(p2gnew$Type, cluster_rows=F, name = "Type", col = gene_type_col)
  pdis <- ComplexHeatmap::Heatmap(log10(abs(p2gnew$Summit2TSS)+1), cluster_rows=F, name = "log10(Summit2TSS + 1)", col = paletteer::paletteer_d("grDevices::blues9")[1:7])
  if(!is.na(save_path)){
    write.table(rna_scale, sprintf("%s/Peak2Gene_Links_RNA_matrix.tsv",save_path),sep='\t',quote=F)
    write.table(annotation, sprintf("%s/Peak2Gene_Links_ATAC_matrix.tsv",save_path),sep='\t',quote=F)
    write.table(p2g, sprintf("%s/Peak2Gene_Links_matrix.tsv",save_path),sep='\t',quote=F)
    pdf(sprintf("%s/ATAC_RNA_All_Links_Heatmap.pdf",save_path),width=fig_width,height=fig_height)
    ComplexHeatmap::draw(p3 + p4 + pdis + ptype + pcor + ppvalue, column_title = sprintf("%s Peak-to-Gene linkages",dim(p2g)[1]))
    dev.off()
  }
  return(ComplexHeatmap::draw(p3 + p4 + pdis + ptype + pcor + ppvalue, column_title = sprintf("%s Peak-to-Gene linkages",dim(p2g)[1])))
}

png("Peak_to_Gene/Peak2Gene_Heatmap.png", width = 800, height = 600)

# Plot heatmap of significant peak-to-gene linkages
# Fill in the file names as needed
plotP2GHeatmap(
  p2g_res = " .rds",
  cluster_N =4,
  cor_cutoff = 0.7,
  atac_matrix = " .tsv",
  rna_matrix = " .tsv"
)

# Close the PNG device
dev.off()
cat("Heatmap saved to Peak2Gene_Heatmap.png\n")