# Breast Cancer Matastasis 2025

Repository to store code used for the analyses presented in "Dissecting The Enhancer Logic in Breast Cancer Metastasis Through Transcriptional, Chromatin Accessibility Profiling and Footprint-Inferred Transcription Factor Activity". 

Please cite our paper published in Molecular Cancer Research: https://

### Data availability:
Processed data available at GEO (https://www.ncbi.nlm.nih.gov/geo/). 

Raw data available with controlled access at dbGaP ([https://www.ncbi.nlm.nih.gov/gap/](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs003363.v1.p1)).


## RNA-seq Processing

1. RNA-seq_Processing/fastqc.sh
2. RNA-seq_Processing/multiqc.sh
3. RNA-seq_Processing/trimming_parallel.sh
4. RNA-seq_Processing/genome_index.sh
5. RNA-seq_Processing/alignment_PE.sh
6. RNA-seq_Processing/counts_PE.sh
7. RNA-seq_Processing/Wranglin_counts_file_for_analisis.py
8. RNA-seq_Processing/Differential_analisis_SV.R
9. RNA-seq_Processing/Merging_MET_GTEX_datasets.R
10. RNA-seq_Processing/Differential_analisis_per_organ.R
11. RNA-seq_Processing/Venn_All_p_val.R
12. RNA-seq_Processing/Hazards_ratio_filtering.R
13. RNA-seq_Processing/Box_Plots.py
14. RNA-seq_Processing/Pathways.R

## ATAC-seq Processing

1. ATAC-seq_Processing/pepatac.yaml
2. ATAC-seq_Processing/project_pipeline_interface.yaml
3. ATAC-seq_Processing/sample_pipeline_interface.yaml
4. ATAC-seq_Processing/Metastasis.yaml
5. ATAC-seq_Processing/Run_ATAC_seq_pipeline.sh
6. ATAC-seq_Processing/Consensus_bigwig.md
7. ATAC-seq Processing/Differential_analisis_Diffbind.R
8. ATAC-seq_Processing/Differential_analisis_Batch_Correction.R
9. ATAC-seq Processing/Consensus_peaks_bedfiles.R
10. ATAC-seq Processing/Genomic_Annotation.R

## Peak-to-Gene Association Analyses

1. Peak-to-Gene/ATAC_log2CPM.R
2. Peak-to-Gene/RNA_log2TPM.R
3. Peak-to-Gene/Merge_Peak_Annotation.R
4. Peak-to-Gene/Peak_to_Gene_Association.R
5. Peak-to-Gene/Peak-to-Gene_Heatmap.R
6. Peak-to-Gene/Scatter_Plots.R
7. Peak-to-Gene/Data_Wrangling_for_peak2gene.R

## Footprinting

1. Footprinting/merge_bam_peak_by_condition.md
2. Footprinting/correct_sequence_bias_atac_seq.R
3. Footprinting/calculate_footprinting_scores.R
4. Footprinting/prepare_unified_peak_set.R
5. Footprinting/differential_tf_binding.R
6. Footprinting/modified_bin_detect.R

## TFAS
1. TFAS/create_dataset_for_tfsee.R
2. TFAS/TFSEE.R
3. TFAS/human_tfs.txt

## Visualizations

1. Visualizations/Genome_tracks.ini
2. Visualizations/Genome_tracks.R
3. Visualizations/Forest_plots.R

