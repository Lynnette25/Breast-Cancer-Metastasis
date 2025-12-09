################################################################################
# Frances Heredia
# Franco Lab
# Description: This script merges Metastasis RNA-seq data with GTEX normal tissue data
#	   to create a comprehensive dataset for further analysis.
################################################################################

library(dplyr)
library(readr)
library(tidyr)

counts_data <- read_csv("/home/fheredia/RNS-seq/Metastasis_Project/RNA_seq/counts_PE.csv")

col_names_dict <- c("R18LIV_1"= "R18_1-Liver_Tumor"
,"R18LIV_2"= "R18_2-Liver_Tumor"
,"R18LUN_1"= "R18_1-Lung_Tumor"
,"R18LUN_2"= "R18_2-Lung_Tumor"
,"R23LIV_1"= "R23_1-Liver_Tumor"
,"R23LIV_2"= "R23_2-Liver_Tumor"
,"R23LUN_1"= "R23_1-Lung_Tumor"
,"R23LUN_2"= "R23_2-Lung_Tumor"
,"R33LIV_1"= "R33_1-Liver_Tumor"
,"R33LIV_2"= "R33_2-Liver_Tumor"
,"R33LUN_1"= "R33_1-Lung_Tumor"
,"R33LUN_2"= "R33_2-Lung_Tumor"

,"R36LIV_1"= "R36_1-Liver_Tumor"
,"R36LIV_2"= "R36_2-Liver_Tumor"
,"R36LUN_1"= "R36_1-Lung_Tumor"
,"R36LUN_2"= "R36_2-Lung_Tumor"
,"R39BR_1"= "R39_1-Breast_Tumor"
,"R39LIV_1"= "R39_1-Liver_Tumor"
,"R39LIV_2"= "R39_2-Liver_Tumor"
,"R49BR_1"= "R49_1-Breast_Tumor"
,"R49BR_2"= "R49_2-Breast_Tumor"
,"R49BR_3"= "R49_3-Breast_Tumor"
,"R49LIV_1"= "R49_1-Liver_Tumor"

,"R49LIV_2"= "R49_4-Liver_Tumor"
,"R49LUN_1"= "R49_1-Lung_Tumor"
,"R49LUN_2"= "R49_2-Lung_Tumor"
,"R8BR_1"= "R8_1-Breast_Tumor"
,"R8BR_2"= "R8_2-Breast_Tumor"
,"R8LIV_1"= "R8_1-Liver_Tumor"
,"R8LIV_2"= "R8_2-Liver_Tumor"
,"R8LUN_1"= "R8_1-Lung_Tumor"
,"R8LUN_2"= "R8_2-Lung_Tumor" )

# Assuming short_gene_data is your data frame and col_names_dict is defined
counts_data <- counts_data %>% rename(!!!setNames(names(col_names_dict), col_names_dict))


normal_breast <- read_tsv("AURORA/gene_reads_v10_breast_mammary_tissue.tsv", locale = locale(encoding = "UTF-8"), skip = 2)

cols <-c(         
"Name", "Description",  
"GTEX-11NV4-2026-SM-5N9DG"
,"GTEX-11TUW-1826-SM-5BC5D"
,"GTEX-11ZTS-2926-SM-GHWOW"
,"GTEX-13FLV-2526-SM-G8W75"
,"GTEX-11DXZ-1926-SM-5GZZL"
,"GTEX-13PVQ-1026-SM-5KM3M"
,"GTEX-14AS3-1626-SM-5S2OY"
,"GTEX-11EMC-2026-SM-5A5JV"
,"GTEX-11EQ9-1826-SM-5Q5AJ"
,"GTEX-11GSP-0926-SM-9WYSG"
,"GTEX-17F96-2426-SM-7IGLN"
,"GTEX-17HGU-1326-SM-79OKB"
,"GTEX-1269C-2426-SM-5FQSN"
,"GTEX-145MO-0826-SM-5NQBL")

normal_breast <- normal_breast %>% dplyr::select(all_of(cols))

col_names_dict <-c(           
"GTEX-11NV4-2026-SM-5N9DG"= "R50-Normal_Breast"
,"GTEX-11TUW-1826-SM-5BC5D" = "R51-Normal_Breast"
,"GTEX-11ZTS-2926-SM-GHWOW" = "R52-Normal_Breast"
,"GTEX-13FLV-2526-SM-G8W75"= "R53-Normal_Breast"
,"GTEX-11DXZ-1926-SM-5GZZL" = "R54-Normal_Breast"
,"GTEX-13PVQ-1026-SM-5KM3M"= "R55-Normal_Breast"
,"GTEX-14AS3-1626-SM-5S2OY" = "R56-Normal_Breast"
,"GTEX-11EMC-2026-SM-5A5JV" = "R57-Normal_Breast"
,"GTEX-11EQ9-1826-SM-5Q5AJ"= "R58-Normal_Breast"
,"GTEX-11GSP-0926-SM-9WYSG" = "R59-Normal_Breast"
,"GTEX-17F96-2426-SM-7IGLN"= "R60-Normal_Breast"
,"GTEX-17HGU-1326-SM-79OKB"= "R61-Normal_Breast"
,"GTEX-1269C-2426-SM-5FQSN"= "R62-Normal_Breast"
,"GTEX-145MO-0826-SM-5NQBL"= "R63-Normal_Breast")

normal_breast <- normal_breast %>% rename(!!!setNames(names(col_names_dict), col_names_dict))

normal_liver <- read_tsv("AURORA/gene_reads_v10_liver.tsv", locale = locale(encoding = "UTF-8"), skip = 2)

cols <-c(  
"Name", "Description",     
"GTEX-11NV4-1326-SM-5HL6V",
"GTEX-11TUW-1726-SM-5BC5C",
"GTEX-11ZTS-1426-SM-5EQMM",
"GTEX-13FLV-0326-SM-5N9DJ",
"GTEX-11DXZ-0126-SM-5EGGY",
"GTEX-13PVQ-1526-SM-5IFEQ",
"GTEX-14AS3-0126-SM-5Q5F4",
"GTEX-11EMC-0326-SM-HAV2K",
"GTEX-11EQ9-0526-SM-5A5JZ",
"GTEX-11GSP-0626-SM-5986T",
"GTEX-17F96-1226-SM-79OK2",
"GTEX-17HGU-1826-SM-7IGQM",
"GTEX-1269C-0626-SM-5FQSS",
"GTEX-145MO-2326-SM-5NQ9K")

normal_liver <- normal_liver %>% dplyr::select(all_of(cols))




col_names_dict_liver <-c(  
"GTEX-11NV4-1326-SM-5HL6V" = "R50-Normal_Liver",
"GTEX-11TUW-1726-SM-5BC5C" = "R51-Normal_Liver",
"GTEX-11ZTS-1426-SM-5EQMM" = "R52-Normal_Liver",
"GTEX-13FLV-0326-SM-5N9DJ" = "R53-Normal_Liver",
"GTEX-11DXZ-0126-SM-5EGGY" = "R54-Normal_Liver",
"GTEX-13PVQ-1526-SM-5IFEQ" = "R55-Normal_Liver",
"GTEX-14AS3-0126-SM-5Q5F4" = "R56-Normal_Liver",
"GTEX-11EMC-0326-SM-HAV2K" = "R57-Normal_Liver",
"GTEX-11EQ9-0526-SM-5A5JZ" = "R58-Normal_Liver",
"GTEX-11GSP-0626-SM-5986T" = "R59-Normal_Liver",
"GTEX-17F96-1226-SM-79OK2" = "R60-Normal_Liver",
"GTEX-17HGU-1826-SM-7IGQM" = "R61-Normal_Liver",
"GTEX-1269C-0626-SM-5FQSS" = "R62-Normal_Liver",
"GTEX-145MO-2326-SM-5NQ9K" = "R63-Normal_Liver")

normal_liver <- normal_liver %>% rename(!!!setNames(names(col_names_dict_liver), col_names_dict_liver))
colnames(normal_liver)

normal_lung<- read_tsv("AURORA/gene_reads_v10_lung.tsv", locale = locale(encoding = "UTF-8"), skip = 2)

cols <-c(  
"Name", "Description",     
"GTEX-11NV4-1126-SM-5HL6J",
"GTEX-11TUW-0526-SM-5LU9A",
"GTEX-11ZTS-1226-SM-5EQMQ",
"GTEX-13FLV-0426-SM-5KLZA",
"GTEX-11DXZ-0726-SM-5N9C4",
"GTEX-13PVQ-0926-SM-5IJFD",
"GTEX-14AS3-0926-SM-5TDD6",
"GTEX-11EMC-0126-SM-5EGKV",
"GTEX-11EQ9-0226-SM-5A5JX",
"GTEX-11GSP-0726-SM-5986L",
"GTEX-17F96-0626-SM-793CC",
"GTEX-17HGU-0926-SM-79OKO",
"GTEX-1269C-0926-SM-5FQSR",
"GTEX-145MO-1326-SM-5Q5EF")


normal_lung <- normal_lung %>% dplyr::select(all_of(cols))


col_names_dict_lung <-c(  
"GTEX-11NV4-1126-SM-5HL6J" = "R50-Normal_Lung",
"GTEX-11TUW-0526-SM-5LU9A" = "R51-Normal_Lung",
"GTEX-11ZTS-1226-SM-5EQMQ" = "R52-Normal_Lung",
"GTEX-13FLV-0426-SM-5KLZA" = "R53-Normal_Lung",
"GTEX-11DXZ-0726-SM-5N9C4" = "R54-Normal_Lung",
"GTEX-13PVQ-0926-SM-5IJFD" = "R55-Normal_Lung",
"GTEX-14AS3-0926-SM-5TDD6" = "R56-Normal_Lung",
"GTEX-11EMC-0126-SM-5EGKV" = "R57-Normal_Lung",
"GTEX-11EQ9-0226-SM-5A5JX" = "R58-Normal_Lung",
"GTEX-11GSP-0726-SM-5986L" = "R59-Normal_Lung",
"GTEX-17F96-0626-SM-793CC" = "R60-Normal_Lung",
"GTEX-17HGU-0926-SM-79OKO" = "R61-Normal_Lung",
"GTEX-1269C-0926-SM-5FQSR" = "R62-Normal_Lung",
"GTEX-145MO-1326-SM-5Q5EF" = "R63-Normal_Lung")

normal_lung <- normal_lung %>% rename(!!!setNames(names(col_names_dict_lung), col_names_dict_lung))
colnames(normal_lung)

merged_data <- normal_breast %>% full_join(normal_liver, by = "Name") %>% full_join(normal_lung, by = "Name")

merged_data$Name <- gsub("\\..*", "", merged_data$Name)

full_data <- counts_data %>% left_join(merged_data, by = c("Geneid" = "Name"))

print(full_data, width = Inf)

write.csv(full_data , "FULL_gene_data_PE_MET_GTEX_one_to_one.csv", row.names = FALSE)


