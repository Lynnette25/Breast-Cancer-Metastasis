################################################################################
# Frances Heredia
# Franco Lab
# Description: This script performs the following tasks  
# 1. Annotates merged peaks with gene information	
# 2. Computes peak-to-gene associations based on correlation analysis	
################################################################################

library(cisDynet)
library(GenomicRanges)
library(IRanges)
library(dplyr)
library(progress)
library(BSDA)
library(tibble)

addAnnotation(gene_bed = "beds/50/hg38_gene_standard.bed",
              gtf = "beds/50/Modified_Homo_sapiens.GRCh38.108.gtf",
              genome_size = "beds/50/hg38.chrom.size")


getPeak2Gene <- function(atac_matrix, rna_matrix, peak_annotation,
                          max_distance=20000, N_permutation=10000, save_path=NA){
  checkGeAnno()
  print("atac_paired_normd")
  atac_paired_norm <- read.table(atac_matrix, head=T, row.names = 1)
  print("rna_paired")
  rna_paired <- read.table(rna_matrix,head=T, row.names = 1)
  logfile("Remove the gene with all expression value is 0.")
  rna_paired <- rna_paired[rowSums(rna_paired) > 0,]
  print("tss")
  tss <- CATAnno$tss %>% as.data.frame()
  tss <- tss %>% distinct(name, .keep_all = TRUE)
  rownames(tss) <- tss$name
  print("tss_df")
  tss_df <- tss[,c(1,2,3,4,6)]
  colnames(tss_df) <- c("chr","start","end","gene","strand")
  gr.tss <- GenomicRanges::makeGRangesFromDataFrame(tss_df,keep.extra.columns=TRUE)
  print("peak_bed")
  peak_bed <- read.table(peak_annotation,head=T)
  peak_bed$Peak <- sprintf("%s:%s-%s",peak_bed$Chromosome, peak_bed$Start, peak_bed$End)
  rownames(peak_bed) <- NULL
  print("peak_gr")
  peak_gr <- GenomicRanges::makeGRangesFromDataFrame(peak_bed, keep.extra.columns=TRUE)

  ## Make the pseudo data for permutation
  logfile("Make the pseudo data for permutation...")
  test_random <- atac_paired_norm[sample(nrow(atac_paired_norm), N_permutation), ]


  ## Calculate the correlation coefficient and p-value for a given gene.
  calculate_pvalue <- function(row_data, df, gene, test_random=test_random){
    corr_random <- apply(test_random, 1, function(row){
      corr <- cor(as.numeric(row_data), row)
    })
    correlations <- apply(df, 1, function(row){
      corr <- cor(as.numeric(row_data), row)
      if (is.numeric(corr)){
      p <- BSDA::z.test(corr_random, mu = corr, sigma.x = 15)$p.value
      } else{
        corr = 0
        p <- BSDA::z.test(corr_random, mu = corr, sigma.x = 15)$p.value
      }
      return(list(corr = corr, p = p))
    })
    result <- data.frame(Peak = rownames(df), Gene = gene)
    result$correlations <- sapply(correlations, function(x) x$corr)
    result$p.value <- sapply(correlations, function(x) x$p)
    return(result)
  }

  ## get the gene2peak links function
  getGenePeaks <- function(gene, max_distance){
    gr.bm_gs <- unique(gr.tss[gr.tss$gene %in% gene])
    GenomicRanges::start(gr.bm_gs) <- GenomicRanges::start(gr.bm_gs) - max_distance
    GenomicRanges::end(gr.bm_gs) <- GenomicRanges::end(gr.bm_gs) + max_distance
    ol <- GenomicRanges::findOverlaps(peak_gr, gr.bm_gs)
    peak_in_window <- peak_gr[unique(IRanges::from(ol))]
    if(length(peak_in_window)!=0){
      df1 <- data.frame(chr=GenomicRanges::seqnames(peak_in_window), start=GenomicRanges::start(peak_in_window), end=GenomicRanges::end(peak_in_window))
      df1$peak <- sprintf("%s:%s-%s", df1$chr, df1$start, df1$end)
      ATAC <- atac_paired_norm[df1$peak,]
      RNA <- rna_paired[gene,]
      res <- calculate_pvalue(RNA, ATAC, gene, test_random=test_random)
      return(res)
    }else{
      return(NA)
    }
  }
  ##  calculate the all gene peak2gene links
  gene_list <- rownames(rna_paired)
  pb <- progress::progress_bar$new(total = length(gene_list))
  result_list <- lapply(c(1:length(gene_list)), function(x) {
    result <- getGenePeaks(gene_list[x],max_distance=max_distance)
    pb$tick()
    return(result)
  })
  na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }
  res_list <- na.omit.list(result_list)
  all_res <- dplyr::bind_rows(res_list)
  final <- merge(all_res,peak_bed,by="Peak",all.x=T)
  final_res <- merge(final,tss_df, by.x="Gene", by.y="gene", all.x=T)
  final_res$Summit2TSS <- final_res$summit - final_res$start
  final_res$orientation <- ifelse((final_res$strand=="+" & final_res$Summit2TSS>=0)|(final_res$strand=="-" & final_res$Summit2TSS<0), "Downstream", "Upstream")
  final_res <- final_res[,c("Peak","Gene","correlations","p.value","Type","summit","start","Summit2TSS","strand","orientation")]
  colnames(final_res)[6:7] <- c("PeakSummit","TSS")
  final_res <- tibble::add_column(final_res, FDR = p.adjust(final_res$p.value, method = "fdr"), .after = 4)
  if(!is.na(save_path)){
    saveRDS(final_res, sprintf("%s/Peak2Gene_All_Links.rds",save_path))
  }
  return(final_res)
}

p2g_res <- getPeak2Gene(
  atac_matrix = "Peak_to_Gene/ATAC_lod2CPM.tsv",
  rna_matrix = "Peak_to_Gene/Filtered_RNA_MET_Liver_Lung_Cancer_all_samples_PE",
  peak_annotation = "Peak_to_Gene/Merged_Peak_Anno.tsv",
  max_distance = 50000,
  N_permutation = 10000,
  save_path = "Peak_to_Gene/"
)