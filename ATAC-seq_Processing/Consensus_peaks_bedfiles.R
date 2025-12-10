################################################################################
# Frances Heredia
# Franco Lab
# Description: This script performs the following tasks  
# 1. Loads ATAC-seq peak data for multiple samples using DiffBind
# 2. Creates consensus peak sets for different conditions (Breast, Lung, Metastasis)
# 3. Outputs consensus peak sets as BED files	
################################################################################



library(DiffBind)

breast_dba <- dba.peakset(NULL, 
                     peaks="A8BR_1/peak_calling_hg38/A8BR_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A8BR_1", 
                     condition="Breast", 
                     replicate=1, 
                     bamReads="A8BR_1/aligned_hg38/A8BR_1_sort_dedup.bam")

# Add the second replicate for A8BR
breast_dba <- dba.peakset(breast_dba, 
                     peaks="A8BR_2/peak_calling_hg38/A8BR_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A8BR_2", 
                     condition="Breast", 
                     replicate=2, 
                     bamReads="A8BR_2/aligned_hg38/A8BR_2_sort_dedup.bam")

# Add the first replicate for A49BR
breast_dba <- dba.peakset(breast_dba, 
                     peaks="A49BR_1/peak_calling_hg38/A49BR_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A49BR_1", 
                     condition="Breast", 
                     replicate=1, 
                     bamReads="A49BR_1/aligned_hg38/A49BR_1_sort_dedup.bam")

# Add the second replicate for A49BR
breast_dba <- dba.peakset(breast_dba, 
                     peaks="A49BR_2/peak_calling_hg38/A49BR_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A49BR_2", 
                     condition="Breast", 
                     replicate=2, 
                     bamReads="A49BR_2/aligned_hg38/A49BR_2_sort_dedup.bam")

# Add the first replicate for A39BR
breast_dba <- dba.peakset(breast_dba, 
                     peaks="A39BR_1/peak_calling_hg38/A39BR_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A39BR_1", 
                     condition="Breast", 
                     replicate=1, 
                     bamReads="A39BR_1/aligned_hg38/A39BR_1_sort_dedup.bam")


nonMet_breast_dba <- dba.peakset(nonMet_breast_dba, 
                     peaks="A_1D91L/peak_calling_hg38/A_1D91L_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A_1D91L", 
                     condition="Breast", 
                     replicate=1, 
                     bamReads="A_1D91L/aligned_hg38/A_1D91L_sort_dedup.bam")


nonMet_breast_dba <- dba.peakset(nonMet_breast_dba, 
                     peaks="A_345CDL/peak_calling_hg38/A_345CDL_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A_345CDL", 
                     condition="Breast", 
                     replicate=1, 
                     bamReads="A_345CDL/aligned_hg38/A_345CDL_sort_dedup.bam")


nonMet_breast_dba <- dba.peakset(nonMet_breast_dba, 
                     peaks="A35B65L1/peak_calling_hg38/A35B65L1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A35B65L1", 
                     condition="Breast", 
                     replicate=1, 
                     bamReads="A35B65L1/aligned_hg38/A35B65L1_sort_dedup.bam")

nonMet_breast_dba <- dba.peakset(nonMet_breast_dba, 
                     peaks="A35B65L2/peak_calling_hg38/A35B65L2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A35B65L2", 
                     condition="Breast", 
                     replicate=2, 
                     bamReads="A35B65L2/aligned_hg38/A35B65L2_sort_dedup.bam")


nonMet_breast_dba <- dba.peakset(nonMet_breast_dba, 
                     peaks="A_3E0C8L/peak_calling_hg38/A_3E0C8L_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A_3E0C8L", 
                     condition="Breast", 
                     replicate=1, 
                     bamReads="A_3E0C8L/aligned_hg38/A_3E0C8L_sort_dedup.bam")


liver_dba <- dba.peakset(NULL, 
                     peaks="A8LIV_1/peak_calling_hg38/A8LIV_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A8LIV_1", 
                     condition="Metastasis", 
                     replicate=1, 
                     bamReads="A8LIV_1/aligned_hg38/A8LIV_1_sort_dedup.bam")

liver_dba <- dba.peakset(liver_dba, 
                     peaks="A8LIV_2/peak_calling_hg38/A8LIV_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A8LIV_2", 
                     condition="Metastasis", 
                     replicate=2, 
                     bamReads="A8LIV_2/aligned_hg38/A8LIV_2_sort_dedup.bam")
liver_dba <- dba.peakset(liver_dba, 
                     peaks="A39LIV_1/peak_calling_hg38/A39LIV_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A39LIV_1", 
                     condition="Metastasis", 
                     replicate=1, 
                     bamReads="A39LIV_1/aligned_hg38/A39LIV_1_sort_dedup.bam")

liver_dba <- dba.peakset(liver_dba, 
                     peaks="A39LIV_2/peak_calling_hg38/A39LIV_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A39LIV_2", 
                     condition="Metastasis", 
                     replicate=2, 
                     bamReads="A39LIV_2/aligned_hg38/A39LIV_2_sort_dedup.bam")

liver_dba <- dba.peakset(liver_dba, 
                     peaks="A49LIV_1/peak_calling_hg38/A49LIV_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A49LIV_1", 
                     condition="Metastasis", 
                     replicate=1, 
                     bamReads="A49LIV_1/aligned_hg38/A49LIV_1_sort_dedup.bam")

liver_dba <- dba.peakset(liver_dba, 
                     peaks="A49LIV_2/peak_calling_hg38/A49LIV_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A49LIV_2", 
                     condition="Metastasis", 
                     replicate=2, 
                     bamReads="A49LIV_2/aligned_hg38/A49LIV_2_sort_dedup.bam")

liver_dba <- dba.peakset(liver_dba, 
                     peaks="A18LIV_1/peak_calling_hg38/A18LIV_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A18LIV_1", 
                     condition="Metastasis", 
                     replicate=1, 
                     bamReads="A18LIV_1/aligned_hg38/A18LIV_1_sort_dedup.bam")

liver_dba <- dba.peakset(liver_dba, 
                     peaks="A18LIV_2/peak_calling_hg38/A18LIV_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A18LIV_2", 
                     condition="Metastasis", 
                     replicate=2, 
                     bamReads="A18LIV_2/aligned_hg38/A18LIV_2_sort_dedup.bam")

liver_dba <- dba.peakset(liver_dba, 
                     peaks="A23LIV_1/peak_calling_hg38/A23LIV_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A23LIV_1", 
                     condition="Metastasis", 
                     replicate=1, 
                     bamReads="A23LIV_1/aligned_hg38/A23LIV_1_sort_dedup.bam")

liver_dba <- dba.peakset(liver_dba, 
                     peaks="A23LIV_2/peak_calling_hg38/A23LIV_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A23LIV_2", 
                     condition="Metastasis", 
                     replicate=2, 
                     bamReads="A23LIV_2/aligned_hg38/A23LIV_2_sort_dedup.bam")

liver_dba <- dba.peakset(liver_dba, 
                     peaks="A33LIV_1/peak_calling_hg38/A33LIV_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A33LIV_1", 
                     condition="Metastasis", 
                     replicate=1, 
                     bamReads="A33LIV_1/aligned_hg38/A33LIV_1_sort_dedup.bam")

liver_dba <- dba.peakset(liver_dba, 
                     peaks="A33LIV_2/peak_calling_hg38/A33LIV_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A33LIV_2", 
                     condition="Metastasis", 
                     replicate=2, 
                     bamReads="A33LIV_2/aligned_hg38/A33LIV_2_sort_dedup.bam")

liver_dba <- dba.peakset(liver_dba, 
                     peaks="A36LIV_1/peak_calling_hg38/A36LIV_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A36LIV_1", 
                     condition="Metastasis", 
                     replicate=1, 
                     bamReads="A36LIV_1/aligned_hg38/A36LIV_1_sort_dedup.bam")

liver_dba <- dba.peakset(liver_dba, 
                     peaks="A36LIV_2/peak_calling_hg38/A36LIV_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A36LIV_2", 
                     condition="Metastasis", 
                     replicate=2, 
                     bamReads="A36LIV_2/aligned_hg38/A36LIV_2_sort_dedup.bam")


# Add peak sets for the "Lung" condition with appropriate replicates

lung_dba <- dba.peakset(NULL, 
                     peaks="A8LUN_1/peak_calling_hg38/A8LUN_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A8LUN_1", 
                     condition="Lung", 
                     replicate=1, 
                     bamReads="A8LUN_1/aligned_hg38/A8LUN_1_sort_dedup.bam")

lung_dba <- dba.peakset(lung_dba, 
                     peaks="A8LUN_2/peak_calling_hg38/A8LUN_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A8LUN_2", 
                     condition="Lung", 
                     replicate=2, 
                     bamReads="A8LUN_2/aligned_hg38/A8LUN_2_sort_dedup.bam")

lung_dba <- dba.peakset(lung_dba, 
                     peaks="A49LUN_1/peak_calling_hg38/A49LUN_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A49LUN_1", 
                     condition="Lung", 
                     replicate=1, 
                     bamReads="A49LUN_1/aligned_hg38/A49LUN_1_sort_dedup.bam")

lung_dba <- dba.peakset(lung_dba, 
                     peaks="A49LUN_2/peak_calling_hg38/A49LUN_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A49LUN_2", 
                     condition="Lung", 
                     replicate=2, 
                     bamReads="A49LUN_2/aligned_hg38/A49LUN_2_sort_dedup.bam")

lung_dba <- dba.peakset(lung_dba, 
                     peaks="A18LUN_1/peak_calling_hg38/A18LUN_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A18LUN_1", 
                     condition="Lung", 
                     replicate=1, 
                     bamReads="A18LUN_1/aligned_hg38/A18LUN_1_sort_dedup.bam")

lung_dba <- dba.peakset(lung_dba, 
                     peaks="A18LUN_2/peak_calling_hg38/A18LUN_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A18LUN_2", 
                     condition="Lung", 
                     replicate=2, 
                     bamReads="A18LUN_2/aligned_hg38/A18LUN_2_sort_dedup.bam")

lung_dba <- dba.peakset(lung_dba, 
                     peaks="A23LUN_1/peak_calling_hg38/A23LUN_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A23LUN_1", 
                     condition="Lung", 
                     replicate=1, 
                     bamReads="A23LUN_1/aligned_hg38/A23LUN_1_sort_dedup.bam")

lung_dba <- dba.peakset(lung_dba, 
                     peaks="A23LUN_2/peak_calling_hg38/A23LUN_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A23LUN_2", 
                     condition="Lung", 
                     replicate=2, 
                     bamReads="A23LUN_2/aligned_hg38/A23LUN_2_sort_dedup.bam")

lung_dba <- dba.peakset(lung_dba, 
                     peaks="A33LUN_1/peak_calling_hg38/A33LUN_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A33LUN_1", 
                     condition="Lung", 
                     replicate=1, 
                     bamReads="A33LUN_1/aligned_hg38/A33LUN_1_sort_dedup.bam")

lung_dba <- dba.peakset(lung_dba, 
                     peaks="A33LUN_2/peak_calling_hg38/A33LUN_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A33LUN_2", 
                     condition="Lung", 
                     replicate=2, 
                     bamReads="A33LUN_2/aligned_hg38/A33LUN_2_sort_dedup.bam")

lung_dba <- dba.peakset(lung_dba, 
                     peaks="A36LUN_1/peak_calling_hg38/A36LUN_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A36LUN_1", 
                     condition="Lung", 
                     replicate=1, 
                     bamReads="A36LUN_1/aligned_hg38/A36LUN_1_sort_dedup.bam")

lung_dba <- dba.peakset( lung_dba, 
                     peaks="A36LUN_2/peak_calling_hg38/A36LUN_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A36LUN_2", 
                     condition="Lung", 
                     replicate=2, 
                     bamReads="A36LUN_2/aligned_hg38/A36LUN_2_sort_dedup.bam")

# Function to create consensus bed for a condition
create_condition_consensus <- function(dba_subset, condition_name) {
  # Subset

  # Count reads if not done
  dba_subset <- dba.count(dba_subset, bParallel=TRUE)
  # Create consensus peaks
  consensus <- dba.peakset(dba_subset, consensus=TRUE, minOverlap=2)
  # Retrieve and write to BED
  consensus_gr <- dba.peakset(consensus, bRetrieve=TRUE)
  write.table(as.data.frame(consensus_gr),
              file=paste0(condition_name, "_consensus.bed"),
              quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}

# For "Breast"

create_condition_consensus(breast_dba, "Breast")

create_condition_consensus(nonMet_breast_dba, "nonMet_Breast")

# For "Lung"

create_condition_consensus(lung_dba, "Lung")

# For "Metastasis" (Liver)

create_condition_consensus(liver_dba, "Liver")

