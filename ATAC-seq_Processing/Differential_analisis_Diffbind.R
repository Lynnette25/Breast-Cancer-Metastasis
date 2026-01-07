################################################################################
# Frances Heredia
# Franco Lab
# Description: This script performs the following tasks  
# 	  1) Differential accessibility analysis between Breast Tumor vs Metastasis (Liver and Lung)
# 	  2) PCA plots
# 	  3) Volcano plots
# 	  4) Heatmaps of differentially accessible regions	
################################################################################

library(DiffBind)

### Load the peak files from MACS2 to DiffBind
# Initialize the MET_df object with the first peakset
MET_df <- dba.peakset(NULL, 
                     peaks="A8BR_1/peak_calling_hg38/A8BR_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A8BR_1", 
                     condition="Breast", 
                     replicate=1, 
                     bamReads="A8BR_1/aligned_hg38/A8BR_1_sort_dedup.bam")

# Add the second replicate for A8BR
MET_df <- dba.peakset(MET_df, 
                     peaks="A8BR_2/peak_calling_hg38/A8BR_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A8BR_2", 
                     condition="Breast", 
                     replicate=2, 
                     bamReads="A8BR_2/aligned_hg38/A8BR_2_sort_dedup.bam")

# Add the first replicate for A49BR
MET_df <- dba.peakset(MET_df, 
                     peaks="A49BR_1/peak_calling_hg38/A49BR_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A49BR_1", 
                     condition="Breast", 
                     replicate=1, 
                     bamReads="A49BR_1/aligned_hg38/A49BR_1_sort_dedup.bam")

# Add the second replicate for A49BR
MET_df <- dba.peakset(MET_df, 
                     peaks="A49BR_2/peak_calling_hg38/A49BR_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A49BR_2", 
                     condition="Breast", 
                     replicate=2, 
                     bamReads="A49BR_2/aligned_hg38/A49BR_2_sort_dedup.bam")

# Add the first replicate for A39BR
MET_df <- dba.peakset(MET_df, 
                     peaks="A39BR_1/peak_calling_hg38/A39BR_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A39BR_1", 
                     condition="Breast", 
                     replicate=1, 
                     bamReads="A39BR_1/aligned_hg38/A39BR_1_sort_dedup.bam")



# Add peak sets for the "Liver" condition with appropriate replicates

MET_df <- dba.peakset(MET_df, 
                     peaks="A8LIV_1/peak_calling_hg38/A8LIV_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A8LIV_1", 
                     condition="Metastasis", 
                     replicate=1, 
                     bamReads="A8LIV_1/aligned_hg38/A8LIV_1_sort_dedup.bam")

MET_df <- dba.peakset(MET_df, 
                     peaks="A8LIV_2/peak_calling_hg38/A8LIV_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A8LIV_2", 
                     condition="Metastasis", 
                     replicate=2, 
                     bamReads="A8LIV_2/aligned_hg38/A8LIV_2_sort_dedup.bam")

MET_df <- dba.peakset(MET_df, 
                     peaks="A39LIV_1/peak_calling_hg38/A39LIV_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A39LIV_1", 
                     condition="Metastasis", 
                     replicate=1, 
                     bamReads="A39LIV_1/aligned_hg38/A39LIV_1_sort_dedup.bam")

MET_df <- dba.peakset(MET_df, 
                     peaks="A39LIV_2/peak_calling_hg38/A39LIV_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A39LIV_2", 
                     condition="Metastasis", 
                     replicate=2, 
                     bamReads="A39LIV_2/aligned_hg38/A39LIV_2_sort_dedup.bam")

MET_df <- dba.peakset(MET_df, 
                     peaks="A49LIV_1/peak_calling_hg38/A49LIV_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A49LIV_1", 
                     condition="Metastasis", 
                     replicate=1, 
                     bamReads="A49LIV_1/aligned_hg38/A49LIV_1_sort_dedup.bam")

MET_df <- dba.peakset(MET_df, 
                     peaks="A49LIV_2/peak_calling_hg38/A49LIV_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A49LIV_2", 
                     condition="Metastasis", 
                     replicate=2, 
                     bamReads="A49LIV_2/aligned_hg38/A49LIV_2_sort_dedup.bam")

MET_df <- dba.peakset(MET_df, 
                     peaks="A18LIV_1/peak_calling_hg38/A18LIV_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A18LIV_1", 
                     condition="Metastasis", 
                     replicate=1, 
                     bamReads="A18LIV_1/aligned_hg38/A18LIV_1_sort_dedup.bam")

MET_df <- dba.peakset(MET_df, 
                     peaks="A18LIV_2/peak_calling_hg38/A18LIV_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A18LIV_2", 
                     condition="Metastasis", 
                     replicate=2, 
                     bamReads="A18LIV_2/aligned_hg38/A18LIV_2_sort_dedup.bam")

MET_df <- dba.peakset(MET_df, 
                     peaks="A23LIV_1/peak_calling_hg38/A23LIV_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A23LIV_1", 
                     condition="Metastasis", 
                     replicate=1, 
                     bamReads="A23LIV_1/aligned_hg38/A23LIV_1_sort_dedup.bam")

MET_df <- dba.peakset(MET_df, 
                     peaks="A23LIV_2/peak_calling_hg38/A23LIV_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A23LIV_2", 
                     condition="Metastasis", 
                     replicate=2, 
                     bamReads="A23LIV_2/aligned_hg38/A23LIV_2_sort_dedup.bam")

MET_df <- dba.peakset(MET_df, 
                     peaks="A33LIV_1/peak_calling_hg38/A33LIV_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A33LIV_1", 
                     condition="Metastasis", 
                     replicate=1, 
                     bamReads="A33LIV_1/aligned_hg38/A33LIV_1_sort_dedup.bam")

MET_df <- dba.peakset(MET_df, 
                     peaks="A33LIV_2/peak_calling_hg38/A33LIV_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A33LIV_2", 
                     condition="Metastasis", 
                     replicate=2, 
                     bamReads="A33LIV_2/aligned_hg38/A33LIV_2_sort_dedup.bam")

MET_df <- dba.peakset(MET_df, 
                     peaks="A36LIV_1/peak_calling_hg38/A36LIV_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A36LIV_1", 
                     condition="Metastasis", 
                     replicate=1, 
                     bamReads="A36LIV_1/aligned_hg38/A36LIV_1_sort_dedup.bam")

MET_df <- dba.peakset(MET_df, 
                     peaks="A36LIV_2/peak_calling_hg38/A36LIV_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A36LIV_2", 
                     condition="Metastasis", 
                     replicate=2, 
                     bamReads="A36LIV_2/aligned_hg38/A36LIV_2_sort_dedup.bam")




# Add peak sets for the "Lung" condition with appropriate replicates

MET_df <- dba.peakset(MET_df, 
                     peaks="A8LUN_1/peak_calling_hg38/A8LUN_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A8LUN_1", 
                     condition="Metastasis", 
                     replicate=1, 
                     bamReads="A8LUN_1/aligned_hg38/A8LUN_1_sort_dedup.bam")

MET_df <- dba.peakset(MET_df, 
                     peaks="A8LUN_2/peak_calling_hg38/A8LUN_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A8LUN_2", 
                     condition="Metastasis", 
                     replicate=2, 
                     bamReads="A8LUN_2/aligned_hg38/A8LUN_2_sort_dedup.bam")

MET_df <- dba.peakset(MET_df, 
                     peaks="A49LUN_1/peak_calling_hg38/A49LUN_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A49LUN_1", 
                     condition="Metastasis", 
                     replicate=1, 
                     bamReads="A49LUN_1/aligned_hg38/A49LUN_1_sort_dedup.bam")

MET_df <- dba.peakset(MET_df, 
                     peaks="A49LUN_2/peak_calling_hg38/A49LUN_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A49LUN_2", 
                     condition="Metastasis", 
                     replicate=2, 
                     bamReads="A49LUN_2/aligned_hg38/A49LUN_2_sort_dedup.bam")

MET_df <- dba.peakset(MET_df, 
                     peaks="A18LUN_1/peak_calling_hg38/A18LUN_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A18LUN_1", 
                     condition="Metastasis", 
                     replicate=1, 
                     bamReads="A18LUN_1/aligned_hg38/A18LUN_1_sort_dedup.bam")

MET_df <- dba.peakset(MET_df, 
                     peaks="A18LUN_2/peak_calling_hg38/A18LUN_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A18LUN_2", 
                     condition="Metastasis", 
                     replicate=2, 
                     bamReads="A18LUN_2/aligned_hg38/A18LUN_2_sort_dedup.bam")

MET_df <- dba.peakset(MET_df, 
                     peaks="A23LUN_1/peak_calling_hg38/A23LUN_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A23LUN_1", 
                     condition="Metastasis", 
                     replicate=1, 
                     bamReads="A23LUN_1/aligned_hg38/A23LUN_1_sort_dedup.bam")

MET_df <- dba.peakset(MET_df, 
                     peaks="A23LUN_2/peak_calling_hg38/A23LUN_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A23LUN_2", 
                     condition="Metastasis", 
                     replicate=2, 
                     bamReads="A23LUN_2/aligned_hg38/A23LUN_2_sort_dedup.bam")

MET_df <- dba.peakset(MET_df, 
                     peaks="A33LUN_1/peak_calling_hg38/A33LUN_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A33LUN_1", 
                     condition="Metastasis", 
                     replicate=1, 
                     bamReads="A33LUN_1/aligned_hg38/A33LUN_1_sort_dedup.bam")

MET_df <- dba.peakset(MET_df, 
                     peaks="A33LUN_2/peak_calling_hg38/A33LUN_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A33LUN_2", 
                     condition="Metastasis", 
                     replicate=2, 
                     bamReads="A33LUN_2/aligned_hg38/A33LUN_2_sort_dedup.bam")

MET_df <- dba.peakset(MET_df, 
                     peaks="A36LUN_1/peak_calling_hg38/A36LUN_1_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A36LUN_1", 
                     condition="Metastasis", 
                     replicate=1, 
                     bamReads="A36LUN_1/aligned_hg38/A36LUN_1_sort_dedup.bam")

MET_df <- dba.peakset(MET_df, 
                     peaks="A36LUN_2/peak_calling_hg38/A36LUN_2_peaks.xls", 
                     peak.caller="macs", 
                     sampID="A36LUN_2", 
                     condition="Metastasis", 
                     replicate=2, 
                     bamReads="A36LUN_2/aligned_hg38/A36LUN_2_sort_dedup.bam")




#### Affinity analysis - differential accessibility analysis
cat("counts")
MET_df_counts<-dba.count(MET_df, bParallel = TRUE, score=DBA_SCORE_READS, summits=TRUE)
MET_df_counts

cat("normalize")
MET_df_normalize <- dba.normalize(MET_df_counts, method=DBA_DESEQ2, normalize=DBA_NORM_NATIVE, library=DBA_LIBSIZE_BACKGROUND, background=TRUE)
#counts<-dba.count(MET_df)
cat("contrast")
counts<-dba.contrast(MET_df_normalize, categories=DBA_CONDITION,minMembers = 2)

cat("analyzed")

MET_df_analyzed<-dba.analyze(counts)
cat("plots")

png('BR_vs_MET/All_volcano_BR_vs_MET.png', units="in", width=7, heigh=7, res=800)
dba.plotVolcano(MET_df_analyzed)
dev.off()

png('BR_vs_MET/All_Heatmap_BR_vs_MET.png', units="in", width=7, heigh=7, res=800)
dba.plotHeatmap(MET_df_analyzed)
dev.off()

png('BR_vs_MET/All_MA_BR_vs_MET.png', units="in", width=7, heigh=7, res=800)
dba.plotMA(MET_df_analyzed)
dev.off()

png('BR_vs_MET/All_PCA_BR_vs_MET.png', units="in", width=7, heigh=7, res=800)
dba.plotPCA(MET_df_analyzed,contrast=1,label=DBA_CONDITION)
dev.off()

profiles <- dba.plotProfile(MET_df_analyzed)
png('BR_vs_MET/All_Profile_BR_vs_MET.png', units="in", width=7, heigh=7, res=800)
dba.plotProfile(profiles)
dev.off()


MET_df_DE_peaks <- dba.report(MET_df_analyzed)
MET_df_DE_peaks_DF<-as.data.frame(MET_df_DE_peaks)
# Filtering
write.table(MET_df_DE_peaks_DF[abs(MET_df_DE_peaks_DF$Fold)>0.5 & MET_df_DE_peaks_DF$FDR<0.05,], "BR_vs_MET/All_diff_regions_BR_vs_MET.txt", sep = "\t", quote = F)
