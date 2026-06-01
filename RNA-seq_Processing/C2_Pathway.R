################################################################################
# Frances Heredia
# Franco Lab
# Description: This script performs the following tasks
#	1. Reads in DESeq2 results 
#	2. Maps ENSEMBL IDs to ENTREZ IDs
#	3. Creates a ranked gene list for GSEA
#	4. Runs GSEA using clusterProfiler with MSigDB C2 gene sets
################################################################################

library(readr)
library(dplyr)
library(tibble)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)

outdir <- "GSEA_results"
dir.create(outdir, showWarnings = FALSE)

# ---- Read DESeq2 results produced in the other env ----
res_df <- read_csv(file.path(outdir, "DESeq2_results_full.csv"), show_col_types = FALSE)


# remove version suffix from ENSEMBL if present
res_df$ENSEMBL <- sub("\\..*$", "", res_df$ENSEMBL)

# ---- Map ENSEMBL -> ENTREZ ----
ensembl_ids <- unique(res_df$ENSEMBL)
mapped <- mapIds(org.Hs.eg.db,
                 keys = ensembl_ids,
                 column = "ENTREZID",
                 keytype = "ENSEMBL",
                 multiVals = "first")
# Add mapping to results
res_df$ENTREZ <- as.character(mapped[res_df$ENSEMBL])

# ---- Remove NA Entrez and collapse duplicates keeping row with max |stat| ----
if(!"stat" %in% colnames(res_df)) stop("Column 'stat' not found in DESeq2 results. Use the Wald test statistic column name 'stat'.")
# ---- Ensure unique ENTREZ -> pick row with max |stat| per ENTREZ ----
res_map <- res_df %>% filter(!is.na(ENTREZ))

res_map_collapsed <- res_map %>%
  group_by(ENTREZ) %>%
  summarize(stat = stat[which.max(abs(stat))], .groups = "drop")

# ---- Create ranked geneList (named numeric vector) ----
geneList <- res_map_collapsed %>%
  filter(!is.na(stat)) %>%
  arrange(desc(stat)) %>%
  pull(stat)
names(geneList) <- res_map_collapsed %>% arrange(desc(stat)) %>% pull(ENTREZ)
# Ensure names are unique (safety)
if(any(duplicated(names(geneList)))){
  stop("Duplicate ENTREZ IDs remain in geneList after collapsing.")
}
geneList <- sort(geneList, decreasing = TRUE)


# Save ranked list
write.table(data.frame(ENTREZ=names(geneList), stat=geneList),
            file = file.path(outdir,"ranked_gene_list_by_stat.rnk"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ---- Load MSigDB C2 (curated gene sets) ----
m_df <- msigdbr(species = "Homo sapiens", category = "C2")
term2gene_c2 <- m_df %>% dplyr::select(gs_name, entrez_gene) %>%
  dplyr::filter(!is.na(entrez_gene)) %>%
  dplyr::rename(term = gs_name, gene = entrez_gene)

# ---- Run GSEA (clusterProfiler::GSEA with TERM2GENE) ----
gsea_c2 <- GSEA(geneList,
                TERM2GENE = term2gene_c2,
                pvalueCutoff = 0.05,
                minGSSize = 15,
                maxGSSize = 500,
                verbose = TRUE)

# Save results and plots
if(!is.null(gsea_c2) && nrow(gsea_c2@result) > 0){
  write.csv(as.data.frame(gsea_c2@result), file = file.path(outdir, "GSEA_C2_results.csv"), row.names = FALSE)

  pdf(file.path(outdir, "GSEA_top30_dotplot.pdf"), width = 8, height = 10)
  print(dotplot(gsea_c2, showCategory = 30) + ggtitle("GSEA C2 (top 30)"))+
    theme(axis.text.y = element_text(size = 6)) +     # smaller pathway labels
    theme(legend.title = element_text(size = 8),
          legend.text  = element_text(size = 7))
  dev.off()

  top_terms <- head(gsea_c2@result$ID, 5)
  pdf(file.path(outdir, "GSEA_top_terms_plots.pdf"), width = 8, height = 6)
  for(t in top_terms){
    print(gseaplot2(gsea_c2, geneSetID = t))
  }
  dev.off()
} else {
  message("No significant C2 pathways found (or GSEA returned no results).")
}


message("Done. Results in: ", normalizePath(outdir))
