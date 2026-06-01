##############################################################################################################
# Frances Heredia
# Franco Lab
# Description: This script performs the following tasks
# 1. Reads RNA-seq data and normal breast, liver metastasis, and lung metastasis gene lists
# 2. Performs differential expression analysis to identify upregulated genes in each tissue
# 3. Performs hypergeometric tests to assess the significance of overlap between upregulated genes and tissue-specific gene lists
##############################################################################################################


```R
library(tidyverse)
library(DESeq2)

counts <- read.csv("/home/fheredia/RNS-seq/Metastasis_Project/RNA_seq/FULL_gene_data.csv")
Liver_genes <- readLines("Liver_genes.txt")
Lung_genes <- readLines("Lung_genes.txt")
Breast_genes <- readLines("Breast_genes.txt")

length(Breast_genes)
length(Liver_genes)
length(Lung_genes)

sum(duplicated(counts$Geneid)) 
dupes <- unique(counts$Geneid[duplicated(counts$Geneid)])

counts_unique <- counts[!duplicated(counts$Geneid), ]
rownames(counts_unique) <- counts_unique$Geneid


annotation_cols <- c("gene_name","Geneid","Chr","Start",
                     "End","Strand","Length",
                     "Description.x","Description.y","Description")

counts_filtered <- counts_unique %>%
  select(-any_of(annotation_cols))



sample_info <- tibble(
  sample = colnames(counts_filtered)
) %>%
  mutate(
    tissue = case_when(
      str_detect(sample, "Breast") ~ "Breast",
      str_detect(sample, "Liver")  ~ "Liver",
      str_detect(sample, "Lung")   ~ "Lung"
    ),
    condition = case_when(
      str_detect(sample, "Normal") ~ "Normal",
      str_detect(sample, "Tumor")  ~ "Tumor"
    )
  ) %>%
  filter(!is.na(tissue)) %>%   
  column_to_rownames("sample")


counts_filtered[is.na(counts_filtered)] <- 0


run_deseq_by_tissue <- function(tissue_name) {

  samples_keep <- rownames(sample_info)[sample_info$tissue == tissue_name]

  dds <- DESeqDataSetFromMatrix(
    countData = counts_filtered[, samples_keep],
    colData   = sample_info[samples_keep, ],
    design    = ~ condition
  )

  dds <- dds[rowSums(counts(dds)) > 10, ]
  dds <- DESeq(dds)

  res <- results(dds, contrast = c("condition", "Tumor", "Normal"))
  res <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj))

  return(list(
    res = res,
    universe = rownames(dds)
  ))
}


breast <- run_deseq_by_tissue("Breast")
liver  <- run_deseq_by_tissue("Liver")
lung   <- run_deseq_by_tissue("Lung")

get_upregulated <- function(res) {
  res %>%
    filter(padj < 0.05, log2FoldChange > 0) %>%
    pull(gene)
}

up_breast <- get_upregulated(breast$res)
up_liver  <- get_upregulated(liver$res)
up_lung   <- get_upregulated(lung$res)


hypergeom_test <- function(setA, setB, universe) {

  N <- length(universe)
  K <- length(setA)
  M <- length(setB)
  k <- length(intersect(setA, setB))

  pval <- phyper(
    q = k - 1,
    m = K,
    n = N - K,
    k = M,
    lower.tail = FALSE
  )

  tibble(
    overlap = k,
    expected = (K * M) / N,
    fold_enrichment = k / ((K * M) / N),
    p_value = pval
  )
}

# Shared tumor programs across tissues
hypergeom_test(Breast_genes, up_breast, breast$universe)
hypergeom_test(Liver_genes, up_liver,  liver$universe)
hypergeom_test(Lung_genes,  up_lung,  lung$universe)

