#!/bin/bash

# Description: This script uses featureCounts to count reads mapped to genes from paired-end RNA-seq data.

# Set working directory and reference paths
RNA_HOME=/home/fheredia/RNS-seq/Metastasis_Project/RNA_seq
ALIGN_DIR=$RNA_HOME/PE_mapped
GTF=/home/fheredia/RNS-seq/Homo_sapiens.GRCh38.108.gtf
COUNTS_OUT=$RNA_HOME/counts_PE.txt

# Run featureCounts in paired-end mode
featureCounts \
    -p --countReadPairs -B -C \ 
    -O -t exon -g gene_id \ 
    -s 0 \ 
    -T 8 \
    -a $GTF \
    -o $COUNTS_OUT \
    $ALIGN_DIR/*.bam

# Check success
if [ $? -ne 0 ]; then
    echo "Error: featureCounts failed."
    exit 1
else
    echo " featureCounts completed successfully. Output: $COUNTS_OUT"
fi