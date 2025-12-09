#!/bin/bash

# This script performs paired-end alignment of RNA-seq reads using STAR aligner.

# Set paths
RNA_HOME=/home/fheredia/RNS-seq/Metastasis_Project/RNA_seq
GENOME_DIR=$RNA_HOME/genomeIndex
TRIMMED_DIR=$RNA_HOME/trimmed
ALIGN_DIR=$RNA_HOME/PE_mapped
mkdir -p $ALIGN_DIR

GTF=/home/fheredia/RNS-seq/Homo_sapiens.GRCh38.108.gtf

# Loop through paired trimmed files
for r1 in $TRIMMED_DIR/*_1_paired.fq.gz; do
    base=$(basename "$r1" _1_paired.fq.gz)
    r2="$TRIMMED_DIR/${base}_2_paired.fq.gz"
    
    echo "Aligning $base..."

    STAR --genomeDir $GENOME_DIR \
         --readFilesIn "$r1" "$r2" \
         --readFilesCommand zcat \
         --runThreadN 8 \
         --outFileNamePrefix $ALIGN_DIR/${base}_ \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped None \
         --outFilterMismatchNmax 3 \
         --outFilterMultimapNmax 1 \
         --outSAMattributes All \
         --sjdbGTFfile $GTF
done

echo "Paired-end STAR alignment complete."