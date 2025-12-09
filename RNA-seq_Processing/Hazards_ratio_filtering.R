################################################################################
# Frances Heredia
# Franco Lab
# Description: This script performs the following tasks  
# 	   1) Filter genes based on Hazard Ratios from Cox Proportional Hazards model
#	   2) Generate Kaplan-Meier plots for filtered genes
# 	   3) Generates summary table saved to CSV
################################################################################

# Preparing the list of genes and Hazards ratios:
library(dplyr)
library(readr)
library(tidyr)
library(RTCGA.clinical)
library(RTCGA.rnaseq)
library(survival)
library(survminer)

# Prepare the clinical data
clinical <- BRCA.clinical %>% 
  mutate(Sample = toupper(`patient.bcr_patient_barcode`)) %>% 
  mutate(OS_event = as.numeric(ifelse(`patient.vital_status` == 'alive', 0, ifelse(`patient.vital_status` == 'dead', 1, NA)))) %>% 
  mutate(OS_days = as.numeric(ifelse(OS_event == 0, `patient.follow_ups.follow_up.days_to_last_followup`, ifelse(OS_event == 1, `patient.follow_ups.follow_up.days_to_death`, NA)))) %>% 
  mutate(OS_years = OS_days / 365.25) %>% 
  select(Sample, OS_event, OS_days, OS_years)

# Modify column names
colnames(BRCA.rnaseq) <- ifelse(grepl("^\\?|bcr_patient_barcode", colnames(BRCA.rnaseq)),
                              colnames(BRCA.rnaseq),
                              sub("\\|.*", "", colnames(BRCA.rnaseq)))
                              

normal_lung<- read_tsv("gene_reads_v10_lung.tsv", locale = locale(encoding = "UTF-8"), skip = 2)
normal_lung <- normal_lung %>% dplyr::select("Name","Description")
normal_lung$Name <- gsub("\\..*", "", normal_lung$Name)
normal_lung$Description <- gsub("-.*$", "", normal_lung$Description)

UP_genes <- read_csv("PE_counts/Up-genes-three-way-gene_intersections.csv")
long_up_genes <- UP_genes %>%
  pivot_longer(cols = 3, names_to = "column", values_to = "EnsemblID")  %>% drop_na(EnsemblID)

nrow(long_up_genes)

merged <- long_up_genes %>% left_join(normal_lung, by = c("EnsemblID" = "Name")) 



TFs <- merged$Description

TFs <- Filter(Negate(is.na), TFs)


# Check which TFs are present in the data frame
existing_TFs <- TFs[TFs %in% colnames(BRCA.rnaseq)]

# Print the number of existing TFs
cat("Number of existing TFs:", length(existing_TFs), "\n")

colnames(BRCA.rnaseq) <- make.unique(colnames(BRCA.rnaseq))

# Prepare the expression data
expr <- BRCA.rnaseq %>% 
  filter(grepl("-01A-", bcr_patient_barcode)) %>% 
  mutate(Sample = sub("-01A-.*", "", bcr_patient_barcode)) %>% 
  mutate(across(all_of(existing_TFs), ~ log2(.), .names = "log2_{.col}")) %>% 
  mutate(across(starts_with("log2_"), ~ ifelse(. == -Inf, NA, .))) %>% 
  select(Sample, starts_with("log2_"))

all_data <- inner_join(clinical, expr, by = "Sample")

str(all_data)

results_list <- list()

cat("Cox proportional hazards regression\n")

# Function to create Kaplan-Meier plots with annotation for a given TF
create_km_plot <- function(data, expression_col, event_col, time_col, tf_name) {
  # Fit the Cox model
  cox_fit <- coxph(as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ", expression_col)), data = data)


  
  # Extract the necessary values from the model
  coef_val <- coef(cox_fit)
  exp_coef <- exp(coef_val)
  conf_int <- confint(cox_fit)
  exp_conf_int <- exp(conf_int)
  p_val <- summary(cox_fit)$coefficients[, "Pr(>|z|)"]
  

  # Calculate the C-index
  c_index <- summary(cox_fit)$concordance[1]
  
  # Format the text to be added to the plot
  output_text <- sprintf("%s Cox HR = %.4f\n95%% CI (%.4f, %.4f)\nCox P val = %.4f\nC-index = %.4f",
                         tf_name, exp_coef, exp_conf_int[1], exp_conf_int[2], p_val, c_index)
  
  # Create factor levels based on the median
  data[[tf_name]] <- factor(ifelse(data[[expression_col]] > median(data[[expression_col]], na.rm = TRUE), 
                                     "High expr", "Low expr" ), 
                            levels = c("Low expr", "High expr"))
  
  # Create the Kaplan-Meier plot with annotation
  fit <- surv_fit(as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ", tf_name)), data = data)
  plot <- ggsurvplot(fit, data = data, risk.table = TRUE, pval = TRUE,
                     xlab = "Time (years)", ylab = "Overall survival", 
                     legend.labs = c("Low expr", "High expr"),palette = c("blue", "red"))
  
  # Add the annotation to the plot
  plot$plot <- plot$plot + 
    annotate("text", x = Inf, y = Inf, label = output_text, hjust = 1, vjust = 1, size = 3.5, 
             fontface = "bold", color = "blue")
  
  return(plot)
}

for (tf_name in existing_TFs) {
    tryCatch({
        col_name <- paste0("log2_", tf_name)
        
        # Check if column exists
        if (!(col_name %in% colnames(all_data))) {
        cat("Column", col_name, "not found in data.\n")
        next
        }
        
        # Filter data
        current_data <- all_data %>% filter(!is.na(get(col_name)))
        if (nrow(current_data) < 3) {
        cat("Skipping", tf_name, "- not enough data.\n")
        next
        }
        
        # Fit the Cox model
        cox_fit <- coxph(as.formula(paste0("Surv(OS_years, OS_event) ~ ", col_name)), data = current_data)
        
        # Extract the necessary values from the model
        coef_val <- coef(cox_fit)
        exp_coef <- exp(coef_val)
        conf_int <- confint(cox_fit)
        exp_conf_int <- exp(conf_int)
        p_val <- summary(cox_fit)$coefficients[, "Pr(>|z|)"]
        
        
        
        # Only proceed if the p-value < 0.05 and exp(coefficient) > 1
        if (p_val >= 0.05 || exp_coef <= 1) {
            cat("Skipping", tf_name, "- does not meet criteria (p_val < 0.05 and exp_coef > 1).\n")
            next
        }
        
        # Store results in a data frame
        results_list[[tf_name]] <- data.frame(
            TF = tf_name,
            Coef = coef_val,
            Exp_Coef = exp_coef,
            Conf_Lower = exp_conf_int[1],
            Conf_Upper = exp_conf_int[2],
            P_value = p_val
        )
        # If criteria are met, create KM plot
        km_plot <- create_km_plot(current_data, col_name, "OS_event", "OS_years", tf_name)
        
        # Save the plot
        png(paste0("PE_counts/KM/KM_Lung_UP_Filtered_Kaplan_Meyer_", tf_name, ".png"), units="in", width=7, height=7, res=800)
        print(km_plot)
        dev.off()
        
    }, error = function(e) {
        # Catch errors, print message, and continue
        cat("Error processing", tf_name, ":", e$message, "\n")
    })
}



# Combine all into a single data frame
results_df <- do.call(rbind, results_list)

# Print or save the results table
write.csv(results_df, "PE_counts/KM_Lung_UP_Filtered_cox_results_PE.csv", row.names = FALSE)