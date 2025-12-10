################################################################################
# Frances Heredia
# Franco Lab
# Description: This script creates a peak annotation file for the peak-to-gene analysis
################################################################################

library(cisDynet)

# Add gene annotations using specified gene bed, GTF, and genome size files
addAnnotation(gene_bed = "beds/50/hg38_gene_standard.bed",
              gtf = "beds/50/Modified_Homo_sapiens.GRCh38.108.gtf",
              genome_size = "beds/50/hg38.chrom.size")

# Annotate merged peaks with specified cutoff and TSS flank parameters
anno <- annoMergedPeaks(quant_data = "Peak_to_Gene/ATAC_log2CPM.tsv", 
                        cutoff = 50000,
                        tss_flank = 1000, 
                        save_path = "Peak_to_Gene", 
                        save_name = "Merged_Peak_Anno")

# Write the annotated peaks to a TSV file
write.table(anno, "Peak_to_Gene/Merged_Peak_Anno.tsv", sep = "\t", quote = FALSE, row.names = FALSE) 
