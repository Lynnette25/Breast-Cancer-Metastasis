################################################################################
# Frances Heredia
# Franco Lab
# Description: This script performs the following tasks  
# 1. Load consensus peak BED files for multiple conditions into GRanges objects.
# 2. Create annotation data from a TxDb object.
# 3. Annotate each peak set using ChIPseeker.
# 4. Generate and save a bar plot visualizing the annotation results.	
################################################################################

library(GenomicRanges)
library(ChIPseeker)
library(org.Hs.eg.db)
library(rtracklayer)
library(ChIPpeakAnno)

library(GenomicRanges)
library(ChIPpeakAnno)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene) # adjust genome build
library(org.Hs.eg.db)

# 1. Load peaks into GRanges
peak_files <- c("nonMet_Breast_consensus_04.bed",
                "Breast_consensus_04.bed",
                "Liver_consensus_04.bed",
                "Lung_consensus_04.bed")

library(GenomicRanges)

load_clean_bed <- function(file) {
  df <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
  # Keep first 4 columns (chr, start, end, score)
  df <- df[, 1:4]
  colnames(df) <- c("chr", "start", "end", "score")
  
  GRanges(seqnames = df$chr,
          ranges = IRanges(start = df$start, end = df$end),
          score = df$score)
}

peakGRs <- lapply(peak_files, load_clean_bed)
names(peakGRs) <- names(peak_files)


# 2. Create annotation data
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
annoData <- toGRanges(txdb, feature="gene")

# 3. Annotate each peak set
peakAnnoList <- lapply(peakGRs, function(x) {
  annotatePeak(x,
               TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
               tssRegion = c(-3000, 3000),
               annoDb = "org.Hs.eg.db")
})

png("Genomic_regions_Consensus_Organs.png", width=800, height=600)
plotAnnoBar(peakAnnoList)
dev.off()

