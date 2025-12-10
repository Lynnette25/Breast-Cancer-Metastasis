
# Breast_Cancer_Metastasis

This script presents an example of how to generate BigWig files from BAM files and create an average BigWig file for multiple samples using R and bash commands.

# In this example, we focus on processing BAM files from non metestatic breast cancer samples.

## Requirements

- R with required libraries
- bedGraphToBigWig tool
- wiggletools
- awk

Ensure all dependencies are installed before running the scripts.

## Generating BigWig Files from BAM Files

The following R code contributes to generating BigWig files from sorted bedGraph files:



```R
sample_names <- list(
   "A_1D91L", "A_345CDL", "A35B65L1",
  "A35B65L2", "A_3E0C8L")


# Generate file paths for BAM files
bam_files <- sapply(sample_names, function(x) {
  sorted_file <- paste0("BedGraph/", x, "_sorted.bedgraph")
  bigwig_file <- paste0("BigWigs/", x, ".bw")
  command3 <- paste( "bedGraphToBigWig",  sorted_file , "hg38.chrom.sizes", bigwig_file )
  cat(command3)
  system(command3)
})
```

# Creating the consensus BigWig
```
nonMets_breast_bw_list.txt
BigWigs/A_1D91L.bw
BigWigs/A_345CDL.bw
BigWigs/A35B65L1.bw
BigWigs/A35B65L2.bw
BigWigs/A_3E0C8L.bw
```

```bash
wiggletools mean $(cat early_breast_bw_list.txt) > early_breast_avg.bedGraph
```

```bash
cd BigWigs

awk 'BEGIN { OFS="\t" } /^fixedStep/ { for (i = 1; i <= NF; i++) { if ($i ~ /^chrom=/) chrom=substr($i,7); if ($i ~ /^start=/) start=substr($i,7); if ($i ~ /^step=/) step=substr($i,6); } i=0; next } /^[0-9.]+$/ { s = start + i * step; e = s + step; print chrom, s, e, $1; i++; next }' early_breast_avg.bedGraph > early_breast_avg_cleaned.bedGraph
```

```bash
awk 'NR==FNR { chrom_size[$1] = $2; next }($3 <= chrom_size[$1]) { print }'  ../hg38.chrom.sizes  early_breast_avg_cleaned.bedGraph > early_breast_avg_cleaned_filtered.bedGraph
```

```bash
bedGraphToBigWig early_breast_avg_cleaned_filtered.bedGraph ../hg38.chrom.sizes early_breast_avg.bw
```
