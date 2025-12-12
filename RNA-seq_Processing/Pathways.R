################################################################################
# Frances Heredia
# Franco Lab
# Description: This script performs the following tasks  
#     1) Pathway enrichment analysis using Hallmark pathways
#        for upregulated genes in different components	
################################################################################

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(msigdbr)

# -------------------------------------------------------------R
# Load Hallmark gene sets (once, outside the loop)
# -------------------------------------------------------------

m_df <- msigdbr(species = "Homo sapiens", category = "H")

head(m_df)

hallmark_df <- m_df %>%
  dplyr::select(gs_name, entrez_gene) %>%
  dplyr::filter(!is.na(entrez_gene)) %>%
  dplyr::rename(TERM = gs_name, GENE = entrez_gene)

# -------------------------------------------------------------
# Components for looping
# -------------------------------------------------------------
# Load gene list
UP_genes <- read_csv("PE_counts/Up-genes-three-way-gene_intersections_p_val.csv")
colnames(UP_genes)[1] <- "V1.0"
components = c(1, 5, 3, 7)

for (component in components) {

    # Map component → composition label
    if (component == 1) {
        composition <- "liver"
		#color <- "#c3e8f8"
    } else if (component == 3) {
        composition <- "lung"
		#color <- "#cbf2d1"
    } else if (component == 5) {
        composition <- "metastasis"
		#color <- "#6a44a8"
    } else if (component == 7) {
        composition <- "cancer"
		#color <-"#6a44a8"
    }

    cat("Processing ", composition, "\n")

    
    column_name <- colnames(UP_genes)[component]

    upregulated_genes <- UP_genes %>%
        pivot_longer(cols =all_of(column_name),
                     names_to = "column",
                     values_to = "EnsemblID")

    # Extract and clean Ensembl IDs
    upregulated_genes_ids <- as.character(upregulated_genes$EnsemblID)
    upregulated_genes_ids <- upregulated_genes_ids[
        !is.na(upregulated_genes_ids) & upregulated_genes_ids != ""
    ]

    # Convert ENSEMBL → ENTREZ
    upregulated_genes_entrez_ids <- mapIds(
        org.Hs.eg.db,
        keys = upregulated_genes_ids,
        column = "ENTREZID",
        keytype = "ENSEMBL",
        multiVals = "list"
    )

    upregulated_genes_entrez_ids <- unlist(upregulated_genes_entrez_ids)
	upregulated_genes_entrez_ids <- na.omit(unique(upregulated_genes_entrez_ids))


    # ---------------------------------------------------------
    # Perform enrichment using Hallmark pathways
    # ---------------------------------------------------------
    enrichment_results <- enricher(
        gene = upregulated_genes_entrez_ids,
        TERM2GENE = hallmark_df)

    # ---------------------------------------------------------
    # If there are enriched pathways, generate plots
    # ---------------------------------------------------------
    if (!is.null(enrichment_results) && nrow(enrichment_results@result) > 0) {
        enrichment_results@result$Description <- gsub("^HALLMARK_", "", enrichment_results@result$Description)

        # -----------------------------------------------------
        # Dotplot
        # -----------------------------------------------------
        png(paste0("Genes_in_", composition, "_hallmark_pathways_p_val_color.png"),
            width = 2000,   height = 2600,  res = 300 )
        print(
            dotplot(enrichment_results, showCategory = 15) +
				geom_point(aes(size = Count)) +
                ggtitle(paste("Upregulated Hallmark Pathways in", composition))
        )
        dev.off()

    } else {
        cat("No enriched Hallmark pathways found for", composition, "\n")
    }
    

}
