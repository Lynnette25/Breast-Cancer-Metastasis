######################################################################################################
# Frances Heredia
# Franco Lab
# Description: This script performs the following tasks
# 1. Reads in ESR1 binding data from primary breast tumors, liver metastases, and lung metastases
# 2. Classifies ESR1 binding sites based on their presence/absence across tissues (rewiring classes)
# 3. Creates a bar plot of the number of ESR1 binding sites in each rewiring class
# 4. Creates a heatmap of ESR1 binding scores across tissues
######################################################################################################

library(dplyr)
library(readxl)
non_mets_Breast_df <- read_xlsx("Footprints/Breast_non-mets-Breast/ESR1_MA0112.4/ESR1_MA0112.4_overview.xlsx")
liver_df <- read_xlsx("Footprints/Breast_Liver/ESR1_MA0112.4/ESR1_MA0112.4_overview.xlsx")
lung_df   <- read_xlsx("Footprints/Breast_Lung/ESR1_MA0112.4/ESR1_MA0112.4_overview.xlsx")

lung <- lung_df %>%
  rename(
	Breast_score_lung = Breast_score,
	Lung_score = Lung_score,
	Breast_bound_lung = Breast_bound,
	Lung_bound = Lung_bound,
	Breast_Lung_log2fc = Breast_Lung_log2fc
  )


liver <- liver_df %>%
  rename(
	Breast_score_liver = Breast_score,
	Liver_score = Liver_score,
	Breast_bound_liver = Breast_bound,
	Liver_bound = Liver_bound,
	Breast_Liver_log2fc = Breast_Liver_log2fc
  )



non_mets_Breast <- non_mets_Breast_df %>%
  rename(
	Breast_score_primary = Breast_score,
	Primary_score = "non-mets-Breast_score",
	Breast_bound_primary = Breast_bound,
	Primary_bound = "non-mets-Breast_bound",
	Breast_Primary_log2fc = "Breast_non-mets-Breast_log2fc"
  )


merged_tfbs <- lung %>%
  full_join(liver,
			by = c("TFBS_chr", "TFBS_start", "TFBS_end")) %>%
  full_join(non_mets_Breast,
			by = c("TFBS_chr", "TFBS_start", "TFBS_end"))


merged_tfbs <- merged_tfbs %>%
  mutate(
	gene_id   = coalesce(gene_id.x, gene_id.y, gene_id),
	gene_name = coalesce(gene_name.x, gene_name.y, gene_name),
	peak_id   = coalesce(peak_id.x, peak_id.y, peak_id),
	peak_chr  = coalesce(peak_chr.x, peak_chr.y, peak_chr),
	peak_start = coalesce(peak_start.x, peak_start.y, peak_start),
	peak_end   = coalesce(peak_end.x, peak_end.y, peak_end)
  ) %>%
  select(-matches("\\.x$"), -matches("\\.y$"))



# no duplicated TFBS rows?
merged_tfbs %>%
  count(TFBS_chr, TFBS_start, TFBS_end) %>%
  filter(n > 1)

# binding columns present?
grep("score|bound|log2fc", colnames(merged_tfbs), value = TRUE)

# NA patterns (expected but informative)
colSums(is.na(merged_tfbs))

tfbs_rewiring <- merged_tfbs %>%
  mutate(
	# replace NA with 0 for binding calls (not bound)
	Primary_bound = ifelse(is.na(Primary_bound), 0, Primary_bound),
	Lung_bound    = ifelse(is.na(Lung_bound), 0, Lung_bound),
	Liver_bound   = ifelse(is.na(Liver_bound), 0, Liver_bound),

	ER_rewiring_class = case_when(
	  # Shared everywhere
	  Primary_bound == 1 & Lung_bound == 1 & Liver_bound == 1 ~ "Shared across all",

	  # Primary only / lost in mets
	  Primary_bound == 1 & Lung_bound == 0 & Liver_bound == 0 ~ "Lost in metastasis",
	  Primary_bound == 1 & Lung_bound == 1 & Liver_bound == 0 ~ "Maintained in lung, lost in liver",
	  Primary_bound == 1 & Lung_bound == 0 & Liver_bound == 1 ~ "Maintained in liver, lost in lung",

	  # Gains
	  Primary_bound == 0 & Lung_bound == 1 & Liver_bound == 0 ~ "Gained in lung",
	  Primary_bound == 0 & Lung_bound == 0 & Liver_bound == 1 ~ "Gained in liver",
	  Primary_bound == 0 & Lung_bound == 1 & Liver_bound == 1 ~ "Gained in metastasis",

	  #Neither
	  Primary_bound == 0 & Lung_bound == 0 & Liver_bound == 0 ~ "Neither",

	  # Everything else
	  TRUE ~ "Other / mixed"
	)
  )

table(tfbs_rewiring$ER_rewiring_class)


library(ggplot2)

plot_df <- tfbs_rewiring %>%
  filter(ER_rewiring_class != "Neither") %>%
  count(ER_rewiring_class) %>%
  mutate(
	ER_rewiring_class = factor(
	  ER_rewiring_class,
	  levels = c(
		"Lost in metastasis",
		"Maintained in lung, lost in liver",
		"Maintained in liver, lost in lung",
		"Gained in lung",
		"Gained in liver",
		"Gained in metastasis",
		"Shared across all"
	  )
	)
  )

palette_aacr <- c(
  "Lost in metastasis"                 = "#004C6D", # deep teal/blue
  "Maintained in lung, lost in liver"  = "#2E8BC0", # medium blue
  "Maintained in liver, lost in lung"  = "#7FB3D5", # light blue
  "Gained in lung"                     = "#F29E4C", # warm orange
  "Gained in liver"                    = "#D9534F", # coral/red
  "Gained in metastasis"               = "#A61E22", # deep red
  "Shared across all"                  = "#6C757D"  # neutral grey
)

png(
  filename = "Footprints/ERE/ESR1_rewiring_classes_TFBS.png",
  width = 2400,
  height = 1600,
  res = 300
)

ggplot(plot_df, aes(x = ER_rewiring_class, y = n, fill = ER_rewiring_class)) +
  geom_col(width = 0.7, color = "black", size = 0.15) +
  scale_fill_manual(values = palette_aacr, guide = "none") +
  labs(
    x = NULL,
    y = "Number of ESR1 binding sites",
    title = "Estrogen receptor binding across primary and metastatic sites"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

dev.off()



tfbs_scores <- merged_tfbs %>%
  mutate(
	# Baseline ESR1 affinity at this TFBS (breast-defined)
	Breast_score_avg = rowMeans(
	  select(
		.,
		Breast_score_primary,
		Breast_score_lung,
		Breast_score_liver
	  ),
	  na.rm = TRUE
	)
  ) %>%
  select(
	TFBS_chr, TFBS_start, TFBS_end,
	Breast_score_avg,
	Primary_score, Lung_score, Liver_score
  )



library(tidyr)

tfbs_scores_long <- tfbs_scores %>%
  pivot_longer(
	cols = c(Primary_score,Breast_score_avg,  Liver_score, Lung_score),
	names_to = "Tissue",
	values_to = "Score"
  ) %>%
  mutate(
	Tissue = factor(
	  Tissue,
	  levels = c("Primary_score","Breast_score_avg" , "Liver_score", "Lung_score"),
	  labels = c("Primary", "Breast", "Liver","Lung")
	)
  )

tfbs_scores_long <- tfbs_scores_long %>%
  left_join(
	tfbs_scores %>%
	  select(TFBS_chr, TFBS_start, TFBS_end, Breast_score_avg),
	by = c("TFBS_chr", "TFBS_start", "TFBS_end")
  )

peak_order <- tfbs_scores_long %>%
  distinct(TFBS_chr, TFBS_start, TFBS_end, Breast_score_avg) %>%
  arrange(desc(Breast_score_avg)) %>%
  mutate(TFBS_id = paste(TFBS_chr, TFBS_start, TFBS_end, sep = ":")) %>%
  pull(TFBS_id)

tfbs_scores_long$TFBS_id <- paste(
  tfbs_scores_long$TFBS_chr,
  tfbs_scores_long$TFBS_start,
  tfbs_scores_long$TFBS_end,
  sep = ":"
)

tfbs_scores_long$TFBS_id <- factor(
  tfbs_scores_long$TFBS_id,
  levels = peak_order
)



p <- ggplot(tfbs_scores_long, aes(x = Tissue, y = TFBS_id, fill = Score)) +
  geom_tile() +
  scale_fill_viridis_c(name = "ESR1 binding score", option = "viridis", na.value = "grey90", limits = c(0,2.5)) +
  labs(
	x = NULL,
	y = NULL,
	title = "ESR1 binding strength across tissues"
  ) +
  theme_classic() +
  theme(
	axis.text.y  = element_blank(),
	axis.ticks.y = element_blank(),
	plot.title   = element_text(hjust = 0.5)
  )

ggsave(
  "Footprints/ERE/ESR1_binding_score_heatmap.png",
  plot = p,
  width = 6,
  height = 8,
  dpi = 300
)