#!/bin/bash

# This script generates a genome index for STAR alignment

# Set your working directory and adapter file paths
RNA_HOME=/home/fheredia/RNS-seq/Metastasis_Project/RNA_seq
RNA_REF_GTF=/home/fheredia/RNS-seq/Homo_sapiens.GRCh38.108.gtf
adapters=$RNA_REFS_DIR/illumina_multiplex.fa
output_dir=$RNA_HOME/trimmed
genome_dir=$RNA_HOME/genomeIndex
mkdir -p $genome_dir

# Create the output directory for STAR alignment if it doesn't exist
alignment_output_dir=$RNA_HOME/STAR_alignment
mkdir -p $alignment_output_dir

# Generate genome index (if not already generated)
STAR --runMode genomeGenerate --genomeDir $genome_dir --genomeFastaFiles $RNA_REFS_DIR/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --sjdbGTFfile $RNA_REF_GTF --runThreadN 2
