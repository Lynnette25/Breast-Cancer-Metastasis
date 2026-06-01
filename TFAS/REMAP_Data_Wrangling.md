##############################################################################################################
# Frances Heredia
# Franco Lab
# Description: This script performs the following tasks
# 1. Downloads Remap 2022 ChIP-seq peak data for all cell lines
# 2. Extracts peaks for breast cancer cell lines (MCF-7, T47D, ZR-75-1, BT-474)
# 3. Merges peaks for each TF across the four cell lines
# 4. Creates a binary matrix of TF binding (1 if peak overlaps footprint, 0 otherwise) for all TFs and all footprints
##############################################################################################################


# Remap Breast Cancer Cell lines
MCF-7 — luminal A; ER+, PR+, HER2-
T47D — luminal A; ER+, PR+, HER2-
ZR-75-1 — luminal B/luminal; ER+, PR+, variable HER2
BT-474 — luminal B; ER+, PR+, HER2+

# Looking at the Remap data from wget https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/remap2022_all_macs2_hg38_v1_0.bed.gz
https
~/RNS-seq/Metastasis_Project/ATAC_seq/Third_run/results_pipeline$ zcat remap2022_all_macs2_hg38_v1_0.bed.gz | head -n 5
chr1    9829    10459   GSE137250.SMARCA4.HeLa  34.47676        .       10084   10085   224,84,252
chr1    9880    10389   ENCSR387QUV.RELB.GM12878        56.61047        .       10091   10092   252,168,56
chr1    9882    10307   ENCSR189YYK.ZBTB40.GM12878      38.13214        .       10080   10081   73,149,46
chr1    9883    10270   ENCSR552XSN.MLLT1.GM12878       15.89233        .       10088   10089   140,28,84
chr1    9884    10289   ENCSR004PLU.ZBTB10.HEK293       21.35251        .       10106   10107   70,67,181


# Modify to keep only useful data
```bash
# Extract lines for the four cell lines (creates remap_breastlines.bed with same BED columns)
zcat remap2022_all_macs2_hg38_v1_0.bed.gz | awk 'BEGIN{IGNORECASE=1} $4 ~ /MCF-7|T47D|ZR-75-1|BT-474/ {print}' > remap_breastlines.bed
```
> head remap_breastlines.bed
chr1    9892    10307   ENCSR318LVG.ZBTB40.MCF-7        47.62095        .       10102   10103   73,149,46
chr1    9905    10243   ENCSR042TWZ.SNIP1.MCF-7 12.56300        .       10087   10088   208,49,131
chr1    9907    10132   ENCSR801SWX.TARDBP.MCF-7        11.44046        .       10063   10064   168,84,112
chr1    9918    10132   ENCSR337NQP.ESRRA.MCF-7 9.48870 .       10057   10058   28,224,56
chr1    9952    10220   ENCSR591EBL.NBN.MCF-7   36.91385        .       10090   10091   140,84,196


# Cellline counts
grep -E -o 'MCF-7|T47D|ZR-75-1|BT-474' remap_breastlines.bed | sort | uniq -c
 204053 BT-474
17224393 MCF-7
  75121 T47D
 189843 ZR-75-1

```bash
# Create folder and list TF names present
mkdir -p beds_by_TF
awk -F'\t' '{ split($4,a,"."); tf=a[2]; print tf }' remap_breastlines.bed | sort -u > tf_list.txt
wc -l tf_list.txt

# Make one merged BED per TF (union of peaks). Adjust chromosome fields if needed (here input bed has chr,start,end in cols 1-3).
while read tf; do
  awk -v TF="$tf" -F'\t' 'BEGIN{OFS="\t"} { split($4,a,"."); if(a[2]==TF) print $1,$2,$3,TF }' remap_breastlines.bed 
 | sort -k1,1 -k2,2n 
 | bedtools merge -i - -c 4 -o collapse \

"beds_by_TF/${tf}_merged.bed"
done < tf_list.txt


# Footprint Matrix
```R
library(GenomicRanges)
library(data.table)

# load if needed
fp_list <- readRDS("footprinting_matrices_list.rds")

# 1) collect all peak ids from all organs
all_peaks <- unique(unlist(lapply(fp_list, rownames)))

# 2) parse "chr:start-end" into a data.frame
df <- do.call(rbind, strsplit(all_peaks, "[:-]"))
df <- data.frame(chr = df[,1], start = as.integer(df[,2]), end = as.integer(df[,3]), stringsAsFactors = FALSE)

# 3) build GRanges and reduce (merge overlapping/adjacent peaks)
gr <- GRanges(seqnames = df$chr, ranges = IRanges(start = df$start, end = df$end))
# optionally restrict to canonical chromosomes:
seqlevelsStyle(gr) <- "UCSC"
gr <- keepSeqlevels(gr, paste0("chr", c(1:22, "X", "Y")), pruning.mode = "coarse")

gr_reduced <- reduce(gr, ignore.strand = TRUE)

# 4) write peaks.bed (BED format: chr, start, end)
peaks_df <- as.data.frame(gr_reduced)[, c("seqnames", "start", "end")]


write.table(peaks_df, file = "Footprints.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```
> head(peaks_df)
  seqnames  start     end
1     chr1 181311  181606
2     chr1 778419  779036
3     chr1 822567  823247
4     chr1 827309  827726
5     chr1 959102  959444
6     chr1 999800 1000862

# ReMap T matrix (peaks × TF)
```bash
mkdir -p tf_Footprints_overlap_counts
for f in beds_by_TF/*_merged.bed; do
  tf=$(basename "$f" _merged.bed)
  bedtools intersect -a Footprints.bed -b "$f" -c > tf_Footprints_overlap_counts/"${tf}.counts.bed"
done

# assemble matrix
#Start with peak IDs
awk '{print $1":"$2"-"$3}' Footprints.bed > ReMap_T_binary_matrix.tsv

#Append one TF column at a time
for f in tf_Footprints_overlap_counts/*.counts.bed; do
  cut -f4 "$f" > col.tmp
  paste ReMap_T_binary_matrix.tsv col.tmp > tmp.tsv
  mv tmp.tsv ReMap_T_binary_matrix.tsv
done
rm -f col.tmp

#Prepend header (TF names)
( printf "peak_id"; for f in tf_Footprints_overlap_counts/*.counts.bed; do tf=$(basename "$f" .counts.bed); printf "\t%s" "$tf"; done; printf "\n" ) | cat - ReMap_T_binary_matrix.tsv > ReMap_T_with_header.tsv

#Binarize (>1 -> 1)
awk 'NR==1{print;next} {for(i=2;i<=NF;i++) if($i>1) $i=1; print}' ReMap_T_with_header.tsv > ReMap_T_binary_matrix_bin.tsv

```
head ReMap_T_binary_matrix_bin.tsv
peak_id AFF4    AHR     AHRR    AR      ARID1A  ARID1B  ARID2   ARID3A  ARNT    ATF7    BMI1    BRD3    BRD4    CARM1   CDK9    CEBPB   CEBPG   CHD1    CLOCK   CREB1   CREBBP  CTBP1   CTBP2   CTCF    CUX1     DDX20   DPF1    DPF2    E2F1    E2F4    E2F8    E4F1    EGR1    ELF1    ELK1    EP300   ERG     ESR1    ESR1_D538G      ESR1_P118       ESR1_pS118      ESR1_Y537C      ESR1_Y537N      ESR1_Y537S       ESR2    ESRRA   ESRRG   FOS     FOSL2   FOXA1   FOXK2   FOXM1   GABPA   GATA3   GATAD2B GRHL1   GRHL2   GTF2F1  HCFC1   HDAC1   HDAC2   HDGF    HES1    HOXB7   HSF1    JUNB    JUN     JUND    KDM5B    KLF10   KLF4    KLF9    KMT2A   KMT2B   MAFK    MAX     MAZ     MBD2    MBD3    MED12   MED1    MEN1    MLLT1   MNT     MSX2    MTA1    MTA2    MTA3    MYC     NBN     NCOA1   NCOA2   NCOA3   NELFA    NEUROD1 NFIB    NFXL1   NIPBL   NONO    NR2F2   NR3C1   NR3C1_mut       NR5A2   NRF1    NRIP1   OVOL1   PAX8    PGR     PKNOX1  PML     PRKDC   RAD21   RAD51   RB1     RCOR1   RELA    REST    RFX1     RFX5    RUNX1   SIN3A   SIX2    SIX4    SMARCA4 SMARCA5 SMARCB1 SMARCC1 SMARCD3 SMARCE1 SMC1A   SNIP1   SP1     SPDEF   SREBF1  SRF     STAG1   STAT3   SUZ12   TAF1    TARDBP  TCF12   TCF7L2  TEAD4    TFAP2A  TFAP2C  TLE3    TP53    TRIM22  TRPS1   UTX     WDHD1   YAP1    YBX1    YY1AP1  ZBTB11  ZBTB1   ZBTB33  ZBTB40  ZBTB7B  ZFX     ZHX2    ZKSCAN1 ZMYM3   ZNF143  ZNF207  ZNF217  ZNF24   ZNF444   ZNF507  ZNF512B ZNF574  ZNF579  ZNF592  ZNF687  ZNF75A  ZNF8    ZSCAN2  ZXDC
chr1:9970-10370 0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       00       0       0       0       0       0       0       1       0       0       0       0       0       0       0       0       0       0       0       0       0       1       0       0       0       00       0       1       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       3 0       0       0       0       0       00       0       0       0       0       0       0       0       0       0       0       0       0       1       0       0       0       0       0       0       0       0       0       0       0       00       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       01       0       1       0       0       0       0       0       0       0       0       1       0       0       0       0       0       0       0       0       0       0       0       0       0       00       0       0       1       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0
chr1:181215-181615      0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0   
