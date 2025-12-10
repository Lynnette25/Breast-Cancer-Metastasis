#!/bin/bash

################################################################################
# Frances Heredia
# Franco Lab
# Description: This script creates a combined set of peaks for comparative analysis
# between Breast, Liver, and Lung conditions, annotates the peaks, and filters the results.
################################################################################

# Create directories for comparative analysis
mkdir -p Footprints/Breast_Liver Footprints/Breast_Lung

# Combine peaks from Breast and Liver into a unified set
cat Footprints/merged/peaks/Breast/Breast_merged.narrowPeak \
    Footprints/merged/peaks/Liver/Liver_merged.narrowPeak | \
    sort -k1,1 -k2,2n | bedtools merge > Footprints/Breast_Liver/Breast_Liver_combined.narrowPeak

# Combine peaks from Breast and Lung into a unified set
cat Footprints/merged/peaks/Breast/Breast_merged.narrowPeak \
    Footprints/merged/peaks/Lung/Lung_merged.narrowPeak | \
    sort -k1,1 -k2,2n | bedtools merge > Footprints/Breast_Lung/Breast_Lung_combined.narrowPeak

# Annotate combined peaks with nearby genes for Breast_Liver
uropa \
    --bed Footprints/Breast_Liver/Breast_Liver_combined.narrowPeak \
    --gtf gencode.v45.annotation.gtf \
    --show_attributes gene_id gene_name \
    --feature_anchor start \
    --distance 20000 10000 \
    --feature gene \
    -o Footprints/Breast_Liver

# Annotate combined peaks with nearby genes for Breast_Lung
uropa \
    --bed Footprints/Breast_Lung/Breast_Lung_combined.narrowPeak \
    --gtf gencode.v45.annotation.gtf \
    --show_attributes gene_id gene_name \
    --feature_anchor start \
    --distance 20000 10000 \
    --feature gene \
    -o Footprints/Breast_Lung

# Extract header for annotated peaks in Breast_Liver
cut -f 1-6,16-17 Footprints/Breast_Liver/Breast_Liver_combined_finalhits.txt | \
    head -n 1 > Footprints/Breast_Liver/Breast_Liver_combined_finalhits_header.txt

# Match the header with annotated peaks for Breast_Liver
cut -f 1-6,16-17 Footprints/Breast_Liver/Breast_Liver_combined_finalhits.txt > \
    Footprints/Breast_Liver/Breast_Liver_combined_finalhits_2.txt

# Filter for standard chromosomes only (removes contigs and alternative assemblies) for Breast_Liver
grep -E '^chr[0-9XY]+\s' Footprints/Breast_Liver/Breast_Liver_combined_finalhits_2.txt > \
    Footprints/Breast_Liver/Breast_Liver_combined_finalhits_3.bed

# Extract header for annotated peaks in Breast_Lung
cut -f 1-6,16-17 Footprints/Breast_Lung/Breast_Lung_combined_finalhits.txt | \
    head -n 1 > Footprints/Breast_Lung/Breast_Lung_combined_finalhits_header.txt

# Match the header with annotated peaks for Breast_Lung
cut -f 1-6,16-17 Footprints/Breast_Lung/Breast_Lung_combined_finalhits.txt > \
    Footprints/Breast_Lung/Breast_Lung_combined_finalhits_2.txt

# Filter for standard chromosomes only (removes contigs and alternative assemblies) for Breast_Lung
grep -E '^chr[0-9XY]+\s' Footprints/Breast_Lung/Breast_Lung_combined_finalhits_2.txt > \
    Footprints/Breast_Lung/Breast_Lung_combined_finalhits_3.bed
