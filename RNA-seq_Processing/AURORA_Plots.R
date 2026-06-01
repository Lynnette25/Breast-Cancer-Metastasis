# AURORA VST Plots
################################################################################
# Frances Heredia
# Franco Lab
# Description: This script performs the following tasks  
# 1. Reads in AURORA RNA-seq counts data
# 2. Normalizes the data using VST
# 3. Creates boxplots comparing expression in primary breast tumors vs liver and lung metastases	
################################################################################


library(readr)
library(DESeq2)
library(dplyr)
library(stringr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(tidyr)

AURORA <- read_csv("AURORA/AURORA_Pri_vs_Mets_p-vals.csv")

AURORA <- AURORA %>%
  select(Name, contains("Tumor"))

count_matrix <- AURORA %>%
  column_to_rownames("Name") %>%
  as.matrix()
count_matrix <- round(count_matrix)

sample_info <- data.frame(
  sample = colnames(count_matrix)
)

sample_info <- sample_info %>%
  mutate(
    patient = str_extract(sample, "R[0-9]+"),
    tissue = case_when(
      grepl("Breast", sample) ~ "Breast",
      grepl("Liver", sample) ~ "Liver",
      grepl("Lung", sample)  ~ "Lung"
    ),
    mets = case_when(
      tissue == "Breast" ~ "Primary",
      tissue %in% c("Liver","Lung") ~ "Metastasis"
    )
  )


rownames(sample_info) <- sample_info$sample

sample_info$mets <- factor(sample_info$mets, levels = c("Primary","Metastasis"))

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ patient + mets
)
dds <- DESeq(dds)
res_Primary_Metastasis <- results(dds, contrast = c("mets","Metastasis","Primary"))


PM <- as.data.frame(res_Primary_Metastasis) %>%
  rownames_to_column("Gene")

PM_sig <- PM %>%
  filter(log2FoldChange > 1, padj < 0.05)

# -------------------------
# Compute VST values
# -------------------------

vsd <- vst(dds, blind = FALSE)

vst_mat <- assay(vsd)

vst_df <- as.data.frame(vst_mat) %>%
  rownames_to_column("Gene")

# -------------------------
# Convert to long format
# -------------------------

vst_long <- vst_df %>%
  pivot_longer(cols = -Gene, names_to = "sample", values_to = "expr")

# Add sample metadata
vst_long <- vst_long %>%
  left_join(sample_info, by = "sample")

# Define tissue groups for plotting
vst_long <- vst_long %>%
  mutate(
    tissue_plot = case_when(
      tissue == "Breast" ~ "A_Pri_Breast",
      tissue == "Liver"  ~ "Liver_Mets",
      tissue == "Lung"   ~ "Lung_Mets"
    )
  )

# -------------------------
# Genes to plot
# -------------------------
gene_list <- c("ISYNA1","RIPK4","RAB3D")

# -------------------------
# Plot loop
# -------------------------

tissue_colors <- c(
  "A_Pri_Breast" = "#ffbaca",
  "Liver_Mets"   = "#c3e8f8",
  "Lung_Mets"    = "#cbf2d1"
)



for (gene in gene_list) {

  gene_df <- vst_long %>%
    filter(Gene == gene)

  comparisons <- list(
    c("A_Pri_Breast","Liver_Mets"),
    c("A_Pri_Breast","Lung_Mets")
  )

  p <- ggplot(gene_df, aes(x = tissue_plot, y = expr, color = tissue_plot, fill = tissue_plot)) +
    geom_boxplot(alpha = 0.4, outlier.shape = NA, color = "gray40") +
    geom_jitter(width = 0.2, size = 2, shape = 21, aes(fill = tissue_colors), color = "gray40") +
    stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.format",method.args =list(alternative = "less")) +
    stat_summary(fun = median, geom = "crossbar", width = 0.6, colour = "gray40", size = 0.6, fatten = 0) +
    scale_color_manual(values = tissue_colors) +
    scale_fill_manual(values = tissue_colors) +
    theme_classic() +
    labs(
      title = paste0(gene, " expression in AURORA"),
      x = "Tissue",
      y = "VST expression"
    )

  ggsave( paste0("AURORA_VST_",gene,".png"),
    p,
    width = 6,
    height = 4,
    dpi = 300
  )
}
