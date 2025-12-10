#!/bin/bash

################################################################################
# Frances Heredia
# Franco Lab
# Description: This script identifies differential TF binding sites between GroupA 
# and GroupB using TOBIAS BINDetect.

# IMPORTANT NOTE: This script runs using the modified version of TOBIAS BinDetect
# that is provided in this repository. This version only runs with a single core.
# Multicore functionality will be added in future versions.

################################################################################

# Define the path to TOBIAS
TOBIAS_CMD="~/miniconda3/bin/TOBIAS"

# Define paths to motifs and genome files
MOTIFS_FILE="JASPAR2024_CORE_non-redundant_pfms.jaspar"
GENOME_FILE="~/pepatac/genome_folder/alias/hg38/fasta/default/hg38.fa"


# Identify differential TF binding sites between Breast and Liver 
$TOBIAS_CMD BINDetect \
    --motifs $MOTIFS_FILE \
    --signals Footprints/Breast_footprints.bw Footprints/Liver_footprints.bw \
    --genome $GENOME_FILE \
    --peaks Footprints/Breast_Liver/Breast_Liver_combined_finalhits_3.bed \
    --peak_header Footprints/Breast_Liver/Breast_Liver_combined_finalhits_header.txt \
    --outdir Footprints/Breast_Liver \
    --cond_names Breast Liver \
    --cores 1

# Identify differential TF binding sites between Breast and Lung 
$TOBIAS_CMD BINDetect \
    --motifs $MOTIFS_FILE \
    --signals Footprints/Breast_footprints.bw Footprints/Lung_footprints.bw \
    --genome $GENOME_FILE \
    --peaks Footprints/Breast_Lung/Breast_Lung_combined_finalhits_3.bed \
    --peak_header Footprints/Breast_Lung/Breast_Lung_combined_finalhits_header.txt \
    --outdir Footprints/Breast_Lung \
    --cond_names Breast Lung \
    --cores 1
