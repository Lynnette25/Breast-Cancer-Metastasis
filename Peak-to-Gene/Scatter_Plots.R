################################################################################
# Frances Heredia
# Franco Lab
# Description: This script performs the following tasks 
# plot all scatterplots for genes that have positive corr in liver and lung.

################################################################################

```R
# Load libraries
library(ggplot2)
library(readr)
library(ggpubr)

# Create output directory for plots if it doesn't exist
if (!dir.exists("Back_to_Peaks/plots_by_organ")) {
  dir.create("Back_to_Peaks/plots_by_organ")
}

# Gene data
my_data <- data.frame(
  Description = c(
    "SPINT1-AS1", "NR2C2AP", "OVOL1", "TATDN1",
    "RCOR2", "WDR5-DT", "HOXC11", "RBM12B-AS1", 
    "MORF4L2-AS1", "SLC2A1-AS1", "SYCP2", "DMBX1",
    "HLTF",  "RECQL4", "CCNF", "PPP1R14D",
    "ATP2C2-AS1", "KCNMB2-AS1", "TPD52",  "OXR1-AS1",
    "GOSR2-DT", "PADI3", "PRAME"
  ),
  Geneid = c(
    "ENSG00000261183", "ENSG00000184162", "ENSG00000172818", "ENSG00000147687", 
    "ENSG00000167771", "ENSG00000273249", "ENSG00000123388", "ENSG00000279331", 
    "ENSG00000231154", "ENSG00000227533", "ENSG00000196074", "ENSG00000197587", 
    "ENSG00000071794", "ENSG00000160957", "ENSG00000162063", "ENSG00000166143", 
    "ENSG00000261286", "ENSG00000237978", "ENSG00000076554", "ENSG00000253582", 
    "ENSG00000261886", "ENSG00000142619", "ENSG00000185686"
  ),
  stringsAsFactors = FALSE
)


# Sample color map
color_dict <- list(
  "A8BR_1" = "#ffbaca", "A8BR_2" = "#ffbaca", "A49BR_1" = "#ffbaca", "A49BR_2" = "#ffbaca", "A39BR_1" = "#ffbaca",

  "A8LIV_1" = "#c3e8f8", "A8LIV_2" = "#c3e8f8", "A39LIV_1" = "#c3e8f8", "A39LIV_2" = "#c3e8f8", 
  "A49LIV_1" = "#c3e8f8", "A49LIV_2" = "#c3e8f8", "A18LIV_1" = "#c3e8f8", "A18LIV_2" = "#c3e8f8", 
  "A23LIV_1" = "#c3e8f8", "A23LIV_2" = "#c3e8f8", "A33LIV_1" = "#c3e8f8", "A33LIV_2" = "#c3e8f8", 
  "A36LIV_1" = "#c3e8f8", "A36LIV_2" = "#c3e8f8",

  "A8LUN_1" = "#cbf2d1", "A8LUN_2" = "#cbf2d1", "A49LUN_1" = "#cbf2d1", "A49LUN_2" = "#cbf2d1", 
  "A18LUN_1" = "#cbf2d1", "A18LUN_2" = "#cbf2d1", "A23LUN_1" = "#cbf2d1", "A23LUN_2" = "#cbf2d1", 
  "A33LUN_1" = "#cbf2d1", "A33LUN_2" = "#cbf2d1", "A36LUN_1" = "#cbf2d1", "A36LUN_2" = "#cbf2d1"
)



# Read matrices once (cached)            
atac_paired_norm <- read.table("Back_to_Peaks/MET_regions_counts_peaks_log2CPM_plus_1.tsv", row.names=1, header=TRUE)
rna_paired <- read.table("Peak_to_Gene/Filtered_RNA_MET_Liver_Lung_Cancer_all_samples_PE.tsv", row.names=1, header=TRUE)

# Define plotting function
# Define a plotting function that separates LIV, BR, and LUN samples
plotPeakGeneBySampleType <- function(peak, gene, gene_name, atac_data, rna_data, color_dict = NULL, legend = TRUE) {
  if (!(peak %in% rownames(atac_data))) {
    message(sprintf("⚠️ Peak '%s' not found in ATAC matrix", peak))
    return(NULL)
  }
  if (!(gene %in% rownames(rna_data))) {
    message(sprintf("⚠️ Gene '%s' not found in RNA matrix", gene))
    return(NULL)
  }

  samples <- colnames(rna_data)
  
  # Determine the sample type
  sample_type <- sapply(samples, function(sample) {
    if (grepl("BR", sample)) "BR"
    else if (grepl("LIV", sample)) "LIV"
    else if (grepl("LUN", sample)) "LUN"
    else "Other"
  })
  
  df <- data.frame(
    ATAC = as.numeric(as.matrix(atac_data)[peak, ]),
    RNA = as.numeric(rna_data[gene, ]),
    SampleType = factor(sample_type, levels = c("BR", "LIV", "LUN")),
    sample = samples
  )
  
  ggplot(df, aes(x = ATAC, y = RNA, color = SampleType, shape = SampleType)) +
    geom_point(size = 3) +
    geom_smooth(method = lm, se = FALSE, aes(group = SampleType, color = SampleType), linetype = "solid") +
    scale_color_manual(values = c("BR" = "#ffbaca", "LIV" = "#c3e8f8", "LUN" = "#cbf2d1")) +
    scale_shape_manual(values = c("BR" = 16, "LIV" = 18, "LUN" = 17)) +
    ggpubr::theme_pubr(legend = ifelse(legend, "right", "none")) +
    theme(panel.border = element_rect(colour = "black", fill = NA)) +
    xlab("ATAC-seq log2(CPM)") +
    ylab("RNA-seq log2(TPM+1)") +
    ggtitle(sprintf("%s:%s\n%s", gene_name, gene, peak))
}

# Loop through genes and their peaks with the new plotting function
peak_gene <- read_rds("Back_to_Peaks/Peaks2Gene_Met_Liver_Lung_Cancer_Links_MET_only_peaks_filt_corr_0.7_ovl.rds")

for (i in 1:nrow(my_data)) {
  row <- my_data[i, ]
  gene_name <- row$Description
  cat(paste0("Plotting ", gene_name, "\n"))
  gene_id <- row$Geneid

  peaks_for_gene <- peak_gene$Peak[peak_gene$Gene == gene_id]

  for (peak in peaks_for_gene) {
    print(peak)
    result_plot <- plotPeakGeneBySampleType(
      peak = peak,
      gene = gene_id,
      gene_name = gene_name,
      atac_data = atac_paired_norm,
      rna_data = rna_paired,
      color_dict = color_dict,
      legend = TRUE
    )

    if (!is.null(result_plot)) {
      filename <- sprintf("Back_to_Peaks/plots_by_organ/%s_%s.png", gene_name, gsub("[:]", "_", peak))
      ggsave(filename, plot = result_plot, width = 5, height = 5)
    }
  }
}
