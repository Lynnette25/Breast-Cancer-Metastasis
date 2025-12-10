################################################################################
# Frances Heredia
# Franco Lab
# Description: This file contains scripts to:
#   1) An R script to move BAM and peak files to directories by condition
#   2) A bash script to merge BAM and peak files by condition
################################################################################

```bash
# Create directories for merged files
mkdir -p Footprints/merged/bam/Breast \
         Footprints/merged/bam/Liver \
		 Footprints/merged/bam/Lung \
         Footprints/merged/peaks/Breast \
         Footprints/merged/peaks/Liver \
		 Footprints/merged/peaks/Lung 
```
 
```R
samples <- c( "A8LIV_2", "A8LIV_1", 
             "A49LIV_1", "A49LIV_2")

# Copying Bams to folder
for (sample in samples) {
  source_path <- file.path(sample, "aligned_hg38", paste0(sample, "_sort_dedup.bam"))
  dest_path <- file.path("Footprints", "merged", "bam", "Liver", paste0(sample, "_sort_dedup.bam"))
  system(paste("cp", shQuote(source_path), shQuote(dest_path)))}

# Copiying Peaks to folder
for (sample in samples) {
  source_path <- file.path(sample, "peak_calling_hg38", paste0(sample, "_peaks.narrowPeak"))
  dest_path <- file.path("Footprints", "merged", "peaks", "Liver", paste0(sample, "_peaks.narrowPeak"))
  system(paste("cp", shQuote(source_path), shQuote(dest_path)))}


samples <- c(  "A8BR_2", "A8BR_1", 
             "A49BR_1", "A49BR_2")

# Copying Bams to folder
for (sample in samples) {
  source_path <- file.path(sample, "aligned_hg38", paste0(sample, "_sort_dedup.bam"))
  dest_path <- file.path("Footprints", "merged", "bam", "Breast", paste0(sample, "_sort_dedup.bam"))
  system(paste("cp", shQuote(source_path), shQuote(dest_path)))}

# Copiying Peaks to folder
for (sample in samples) {
  source_path <- file.path(sample, "peak_calling_hg38", paste0(sample, "_peaks.narrowPeak"))
  dest_path <- file.path("Footprints", "merged", "peaks", "Breast", paste0(sample, "_peaks.narrowPeak"))
  system(paste("cp", shQuote(source_path), shQuote(dest_path)))}

samples <- c( "A49LUN_1", "A49LUN_2",
			  "A8LUN_1", "A8LUN_2")		

# Copying Bams to folder
for (sample in samples) {
  source_path <- file.path(sample, "aligned_hg38", paste0(sample, "_sort_dedup.bam"))
  dest_path <- file.path("Footprints", "merged", "bam", "Lung", paste0(sample, "_sort_dedup.bam"))
  system(paste("cp", shQuote(source_path), shQuote(dest_path)))}

# Copiying Peaks to folder
for (sample in samples) {
  source_path <- file.path(sample, "peak_calling_hg38", paste0(sample, "_peaks.narrowPeak"))
  dest_path <- file.path("Footprints", "merged", "peaks", "Lung", paste0(sample, "_peaks.narrowPeak"))
  system(paste("cp", shQuote(source_path), shQuote(dest_path)))}	 
```

```bash
# Merge the BAM files 
samtools merge -@ 16 Footprints/merged/bam/Breast/Breast_merged.bam  Footprints/merged/bam/Breast/*.bam
samtools merge -@ 16 Footprints/merged/bam/Liver/Liver_merged.bam  Footprints/merged/bam/Liver/*.bam
samtools merge -@ 16 Footprints/merged/bam/Lung/Lung_merged.bam  Footprints/merged/bam/Lung/*.bam

rm Footprints/merged/bam/Breast/A*
rm Footprints/merged/bam/Liver/A*
rm Footprints/merged/bam/Lung/A*

cat Footprints/merged/peaks/Breast/*.narrowPeak | sort -k1,1 -k2,2n | bedtools merge > Footprints/merged/peaks/Breast/Breast_merged.narrowPeak
cat Footprints/merged/peaks/Liver/*.narrowPeak | sort -k1,1 -k2,2n | bedtools merge > Footprints/merged/peaks/Liver/Liver_merged.narrowPeak
cat Footprints/merged/peaks/Lung/*.narrowPeak | sort -k1,1 -k2,2n | bedtools merge >  Footprints/merged/peaks/Lung/Lung_merged.narrowPeak

rm Footprints/merged/peaks/Breast/A*
rm Footprints/merged/peaks/Liver/A*
rm Footprints/merged/peaks/Lung/A*
```
