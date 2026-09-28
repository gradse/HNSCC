# Analysis stage from archive/HNSCC_Analysis.Rmd (source lines 2068-3194).
# Run through scripts/run_r_analysis.R.

# Source notebook line 2068
ori_methyl <- read_parquet("methyl_calculated_removed_indeter_hpv.parquet")
ori_methyl <- as.data.frame(ori_methyl)

# Source notebook line 2075
gene_symbol_counts <- ori_methyl %>%
  group_by(Gene_Symbol) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

print(gene_symbol_counts)

# Source notebook line 2087
# This is the data that have been cleaned for both hpv "indeterminate" and "nan" rows
methyl_merge_calculated_cleaned <- open_dataset("methyl_calculated_removed_indeter_hpv_and_missing_beta.parquet")
methyl_merge_calculated_cleaned <- as.data.frame(methyl_merge_calculated_cleaned)

# Source notebook line 2095
methyl_merge_calculated_cleaned %>%
  group_by(Gene_Symbol) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

# Source notebook line 2150
methyl_gene <- open_dataset("methyl_summed_genes.parquet")
methyl_gene <- as.data.frame(methyl_gene)

# Source notebook line 2193
methyl_gene_cleaned <- open_dataset("methyl_cleaned_genes.parquet")
methyl_gene_cleaned <- as.data.frame(methyl_gene_cleaned)

methyl_gene_cleaned <- methyl_gene_cleaned[, -ncol(methyl_gene_cleaned)]

# Source notebook line 2203
transposed_methyl_gene_cleaned <- t(methyl_gene_cleaned[, -1])

# Convert it to a data frame
transposed_methyl_gene_cleaned <- as.data.frame(transposed_methyl_gene_cleaned)

# Set the column names to the values in the first row (Gene Symbols)
colnames(transposed_methyl_gene_cleaned) <- methyl_gene_cleaned$Gene_Symbol

# Set the row names to the original column names (Patient IDs)
rownames(transposed_methyl_gene_cleaned) <- colnames(methyl_gene_cleaned)[-1]

# Process clinical data (sort rownames)
sorted_clinical_focused_filtered_for_age <- clinical_focused_filtered_for_age[order(rownames(clinical_focused_filtered_for_age)), , drop = FALSE]
row_name_clinical <- rownames(sorted_clinical_focused_filtered_for_age)


# Merge them

# Check if the row names are consistent
all_row_names_match <- all(rownames(transposed_methyl_gene_cleaned) == row_name_clinical)
if (all_row_names_match) {
  gene_with_clinical_data <- cbind(sorted_clinical_focused_filtered_for_age, transposed_methyl_gene_cleaned)
  print("The row names are consistent and the data has been merged successfully.")
} else {
  print("The row names do not match. Please check the data.")
}

# Source notebook line 2238
# Identify indices for HPV positive and negative patients
hpv_positive_indices <- which(gene_with_clinical_data$hpv_status == "positive")
hpv_negative_indices <- which(gene_with_clinical_data$hpv_status == "negative")

# Function to perform Mann-Whitney U test for a single gene
perform_test <- function(gene_index) {
  gene_values_positive <- gene_with_clinical_data[hpv_positive_indices, gene_index]
  gene_values_negative <- gene_with_clinical_data[hpv_negative_indices, gene_index]
  
  test_result <- wilcox.test(gene_values_positive, gene_values_negative)
  
  return(c(
    gene = colnames(gene_with_clinical_data)[gene_index],
    p_value = test_result$p.value,
    W = test_result$statistic
  ))
}

# Use parallel processing to speed up the computation
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("gene_with_clinical_data", "hpv_positive_indices", "hpv_negative_indices", "perform_test"))

# Perform the test for all gene columns (10th to last)
test_results <- parLapply(cl, 10:ncol(gene_with_clinical_data), perform_test)
stopCluster(cl)

# Combine results into a data frame
test_results_df <- as.data.frame(do.call(rbind, test_results))
test_results_df$p_value <- as.numeric(as.character(test_results_df$p_value))
test_results_df$W <- as.numeric(as.character(test_results_df$W))

# Remove the duplicated W.W column
test_results_df$W.W <- NULL

# Format p-values to a specified number of significant digits
test_results_df$p_value <- format(test_results_df$p_value, digits = 4, scientific = TRUE)

# Save the results to a CSV file for easy checking
write.csv(test_results_df, "gene_comparison_results.csv", row.names = FALSE)

# Print a summary of significant results (e.g., p < 0.05)
significant_results <- test_results_df %>% filter(as.numeric(p_value) < 0.05)
print(significant_results)

# Save the results to a CSV file for easy checking
write.csv(significant_results, "gene_comparison_significant_results.csv", row.names = FALSE)

# Source notebook line 2294
methyl_gene <- open_dataset("methyl_summed_genes.parquet")
methyl_gene <- as.data.frame(methyl_gene)

# Source notebook line 2300
clinical_all_CDEs <- fread("clinical_All_CDEs.txt", sep = "\t", header = TRUE)

# Transpose the clinical data
clinical_transposed <- t(clinical_all_CDEs)

clinical_transposed <- as.data.frame(clinical_transposed)

colnames(clinical_transposed) <- clinical_transposed[1,]

clinical_transposed <- clinical_transposed[-1,]

# Select necessary variables
clinical_focused_for_age <- clinical_transposed[, c(9, 12, 21:23, 46, 48, 54, 114)]

# Filter data
clinical_focused_filtered_for_age <- clinical_focused_for_age %>%
  filter(hpv_status != "indeterminate", tumor_tissue_site == "head and neck")

# Source notebook line 2322
# HPV+ group
hpv_positive <- clinical_focused_filtered_for_age %>% filter(hpv_status == "positive")

# HPV- group
hpv_negative <- clinical_focused_filtered_for_age %>% filter(hpv_status == "negative")


hpv_pos_patient <- rownames(hpv_positive)
hpv_neg_patient <- rownames(hpv_negative)

# Find the maximum length
max_length <- max(length(hpv_pos_patient), length(hpv_neg_patient))

# Pad the shorter list with NA values
hpv_pos_patient <- c(hpv_pos_patient, rep(NA, max_length - length(hpv_pos_patient)))
hpv_neg_patient <- c(hpv_neg_patient, rep(NA, max_length - length(hpv_neg_patient)))

# Create a data frame with the two lists
patient_lists <- data.frame(HPV_Positive = hpv_pos_patient, HPV_Negative = hpv_neg_patient)

# Write the data frame to a CSV file
write.csv(patient_lists, "patient_lists.csv", row.names = FALSE)

# Source notebook line 2398
methyl_gene_cleaned_bygroup <- open_dataset("methyl_cleaned_bygroup_genes.parquet")
methyl_gene_cleaned_bygroup <- as.data.frame(methyl_gene_cleaned_bygroup)

methyl_gene_cleaned_bygroup <- methyl_gene_cleaned_bygroup[, -ncol(methyl_gene_cleaned_bygroup)]

# Source notebook line 2411
methyl_gene_cleaned_bygroup[, 2:ncol(methyl_gene_cleaned_bygroup)] <- log2(1 + methyl_gene_cleaned_bygroup[, 2:ncol(methyl_gene_cleaned_bygroup)])

# Source notebook line 2422
transposed_methyl_gene_cleaned_bygroup <- t(methyl_gene_cleaned_bygroup[, -1])

# Convert it to a data frame
transposed_methyl_gene_cleaned_bygroup <- as.data.frame(transposed_methyl_gene_cleaned_bygroup)

# Set the column names to the values in the first row (Gene Symbols)
colnames(transposed_methyl_gene_cleaned_bygroup) <- methyl_gene_cleaned_bygroup$Gene_Symbol

# Set the row names to the original column names (Patient IDs)
rownames(transposed_methyl_gene_cleaned_bygroup) <- colnames(methyl_gene_cleaned_bygroup)[-1]

# Process clinical data (sort rownames)
sorted_clinical_focused_filtered_for_age <- clinical_focused_filtered_for_age[order(rownames(clinical_focused_filtered_for_age)), , drop = FALSE]
row_name_clinical <- rownames(sorted_clinical_focused_filtered_for_age)


# Merge them

# Check if the row names are consistent
all_row_names_match <- all(rownames(transposed_methyl_gene_cleaned_bygroup) == row_name_clinical)
if (all_row_names_match) {
  gene_with_clinical_data_bygroup <- cbind(sorted_clinical_focused_filtered_for_age, transposed_methyl_gene_cleaned_bygroup)
  print("The row names are consistent and the data has been merged successfully.")
} else {
  print("The row names do not match. Please check the data.")
}

# Source notebook line 2458
# Identify indices for HPV positive and negative patients
hpv_positive_indices <- which(gene_with_clinical_data_bygroup$hpv_status == "positive")
hpv_negative_indices <- which(gene_with_clinical_data_bygroup$hpv_status == "negative")

# Function to perform Mann-Whitney U test for a single gene
perform_test_bygroup <- function(gene_index) {
  gene_values_positive <- gene_with_clinical_data_bygroup[hpv_positive_indices, gene_index]
  gene_values_negative <- gene_with_clinical_data_bygroup[hpv_negative_indices, gene_index]
  
  # Only perform the test if both groups have enough non-zero values
  if (length(gene_values_positive[gene_values_positive != 0]) >= 2 & length(gene_values_negative[gene_values_negative != 0]) >= 2) {
    test_result <- wilcox.test(gene_values_positive, gene_values_negative)
    difference <- median(gene_values_positive) - median(gene_values_negative)
    ratio <- 2^difference - 1
    return(c(
      gene = colnames(gene_with_clinical_data_bygroup)[gene_index],
      p_value = test_result$p.value,
      difference = difference,
      ratio = ratio
    ))
  } else {
    return(c(
      gene = colnames(gene_with_clinical_data_bygroup)[gene_index],
      p_value = NA,
      difference = NA,
      ratio = NA
    ))
  }
}

# Use parallel processing to speed up the computation
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("gene_with_clinical_data_bygroup", "hpv_positive_indices", "hpv_negative_indices", "perform_test"))

# Perform the test for all gene columns (10th to last)
test_results_bygroup <- parLapply(cl, 10:ncol(gene_with_clinical_data_bygroup), perform_test_bygroup)
stopCluster(cl)

# Combine results into a data frame
test_results_df_bygroup <- as.data.frame(do.call(rbind, test_results_bygroup))
test_results_df_bygroup$p_value <- as.numeric(as.character(test_results_df_bygroup$p_value))
test_results_df_bygroup$difference <- as.numeric(as.character(test_results_df_bygroup$difference))
test_results_df_bygroup$ratio <- as.numeric(as.character(test_results_df_bygroup$ratio))

# Calculate FDR
test_results_df_bygroup <- test_results_df_bygroup %>% 
  mutate(FDR = p.adjust(p_value, method = "fdr"))

# Format p-values to a specified number of significant digits
test_results_df_bygroup$p_value <- format(test_results_df_bygroup$p_value, digits = 4, scientific = TRUE)

# Sort the data frame by the difference column
test_results_df_bygroup <- test_results_df_bygroup %>% arrange(desc(difference))

# Save the results to a CSV file for easy checking
write.csv(test_results_df_bygroup, "gene_comparison_results_bygroup.csv", row.names = FALSE)

test_results_df_bygroup$p_value <- as.numeric(as.character(test_results_df_bygroup$p_value))
test_results_df_bygroup$FDR <- as.numeric(as.character(test_results_df_bygroup$FDR))
test_results_df_bygroup$difference <- as.numeric(as.character(test_results_df_bygroup$difference))

# Print the results
print(test_results_df_bygroup)

# Source notebook line 2533
# Filter genes based on p-value < 0.05
filtered_data_pvalue <- test_results_df_bygroup %>%
  filter(as.numeric(p_value) < 0.05) %>%
  mutate(
    negLog10PValue = -log10(as.numeric(p_value)),
    significant = as.numeric(p_value) < 0.05
  )

# Create the volcano plot for p-value < 0.05
volcano_plot_pvalue <- ggplot(filtered_data_pvalue, aes(x = difference, y = negLog10PValue)) +
  geom_point(aes(color = significant)) +
  scale_color_manual(values = c("black", "red")) +
  labs(
    title = "Volcano Plot (p-value < 0.05)",
    x = "Log2 Fold Change",
    y = "-Log10 P-Value"
  ) +
  theme_minimal()

# Print the first volcano plot
print(volcano_plot_pvalue)



ggplot(filtered_data_pvalue, aes(x = difference, y = -log10(p_value))) +
  geom_point(aes(color = p_value < 0.05)) +
  scale_color_manual(values = c("black", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "blue") +
  labs(title = "Volcano Plot (p-value < 0.05)", x = "Log2 Fold Change", y = "-Log10 P-Value") +
  theme_minimal()

# Source notebook line 2575
# Filter genes based on FDR < 0.001
filtered_data_fdr <- test_results_df_bygroup %>%
  filter(as.numeric(FDR) < 0.001) %>%
  mutate(
    negLog10PValue = -log10(as.numeric(p_value)),
    significant = as.numeric(FDR) < 0.001
  )

# Create the volcano plot for FDR < 0.001
volcano_plot_fdr <- ggplot(filtered_data_fdr, aes(x = difference, y = negLog10PValue)) +
  geom_point(aes(color = significant)) +
  scale_color_manual(values = c("black", "red")) +
  labs(
    title = "Volcano Plot (FDR < 0.001)",
    x = "Log2 Fold Change",
    y = "-Log10 P-Value"
  ) +
  theme_minimal()

# Print the second volcano plot
print(volcano_plot_fdr)



# Create the volcano plot for the FDR filtered data
ggplot(filtered_data_fdr, aes(x = difference, y = -log10(p_value))) +
  geom_point(aes(color = FDR < 0.001)) +
  scale_color_manual(values = c("black", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "blue") +
  labs(title = "Volcano Plot (FDR < 0.001)", x = "Log2 Fold Change", y = "-Log10 P-Value") +
  theme_minimal()


# Source notebook line 2622
# Identify indices for HPV positive and negative patients
hpv_positive_indices <- which(gene_with_clinical_data_bygroup$hpv_status == "positive")
hpv_negative_indices <- which(gene_with_clinical_data_bygroup$hpv_status == "negative")

# Function to perform t-test for a single gene
perform_test_bygroup_ttest <- function(gene_index) {
  gene_values_positive <- gene_with_clinical_data_bygroup[hpv_positive_indices, gene_index]
  gene_values_negative <- gene_with_clinical_data_bygroup[hpv_negative_indices, gene_index]
  
  # Only perform the test if both groups have enough non-zero values
  if (length(gene_values_positive[gene_values_positive != 0]) >= 2 & length(gene_values_negative[gene_values_negative != 0]) >= 2) {
    test_result <- t.test(gene_values_positive, gene_values_negative)
    difference <- test_result$estimate[1] - test_result$estimate[2]
    ratio <- 2^difference
    FC <- ifelse(ratio > 1, ratio, -1/ratio)
    
    return(c(
      gene = colnames(gene_with_clinical_data_bygroup)[gene_index],
      `difference_log2_scale (pos - neg)` = difference,
      `ratio (pos / neg)` = ratio,
      FC = FC,
      p_value = test_result$p.value
    ))
  } else {
    return(c(
      gene = colnames(gene_with_clinical_data_bygroup)[gene_index],
      `difference_log2_scale (pos - neg)` = NA,
      `ratio (pos / neg)` = NA,
      FC = NA,
      p_value = NA
    ))
  }
}

# Use parallel processing to speed up the computation
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("gene_with_clinical_data_bygroup", "hpv_positive_indices", "hpv_negative_indices", "perform_test_bygroup_ttest"))

# Perform the test for all gene columns (10th to last)
test_results_bygroup_ttest <- parLapply(cl, 10:ncol(gene_with_clinical_data_bygroup), perform_test_bygroup_ttest)
stopCluster(cl)

# Combine results into a data frame
test_results_df_bygroup_ttest <- as.data.frame(do.call(rbind, test_results_bygroup_ttest))
test_results_df_bygroup_ttest$p_value <- as.numeric(as.character(test_results_df_bygroup_ttest$p_value))
test_results_df_bygroup_ttest$`difference_log2_scale (pos - neg)` <- as.numeric(as.character(test_results_df_bygroup_ttest$`difference_log2_scale (pos - neg)`))
test_results_df_bygroup_ttest$`ratio (pos / neg)` <- as.numeric(as.character(test_results_df_bygroup_ttest$`ratio (pos / neg)`))
test_results_df_bygroup_ttest$FC <- as.numeric(as.character(test_results_df_bygroup_ttest$FC))

# Calculate FDR
test_results_df_bygroup_ttest <- test_results_df_bygroup_ttest %>% 
  mutate(FDR = p.adjust(p_value, method = "fdr"))

# Sort the data frame by the difference column
test_results_df_bygroup_ttest <- test_results_df_bygroup_ttest %>% arrange(desc(`difference_log2_scale (pos - neg)`))

# Remove 2nd to 4th columns that are character types
test_results_df_bygroup_ttest <- test_results_df_bygroup_ttest[, -c(2, 3, 4)]

# Save the results to a CSV file for easy checking
write.csv(test_results_df_bygroup_ttest, "gene_comparison_results_bygroup_ttest.csv", row.names = FALSE)

# Print the results
print(test_results_df_bygroup_ttest)

# Source notebook line 2695
# Determine genes outside of the vertical lines
significant_genes_difference <- test_results_df_bygroup_ttest %>%
  filter(`difference_log2_scale (pos - neg)` > log2(1.3) | `difference_log2_scale (pos - neg)` < -log2(1.3))

# Select the top 10 largest and lowest genes based on difference_log2_scale (pos - neg)
top_genes_diff <- significant_genes_difference %>%
  arrange(desc(`difference_log2_scale (pos - neg)`)) %>%
  head(10)

lowest_genes_diff <- significant_genes_difference %>%
  arrange(`difference_log2_scale (pos - neg)`) %>%
  head(10)

# Combine top and lowest genes
label_genes_diff <- bind_rows(top_genes_diff, lowest_genes_diff)

# Create the volcano plot
volcano_plot_difference <- ggplot(test_results_df_bygroup_ttest, aes(x = `difference_log2_scale (pos - neg)`, y = -log10(p_value))) +
  geom_point(aes(color = p_value < 0.05), alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("black", "red")) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "blue") +
  geom_text_repel(data = label_genes_diff, aes(label = gene), size = 3) +
  labs(title = "Volcano Plot (T-test Results)", x = "Difference_log2_scale (pos - neg)", y = "-Log10 P-Value") +
  theme_minimal()

# Print the plot
print(volcano_plot_difference)

# Source notebook line 2730
# Identify significant genes based on ratio
significant_genes_ratio <- test_results_df_bygroup_ttest %>%
  filter(`ratio (pos / neg)` > 1.3 | `ratio (pos / neg)` < 1/1.3)

# Select the top 10 largest and lowest genes based on ratio
top_genes_ratio <- significant_genes_ratio %>%
  arrange(desc(`ratio (pos / neg)`)) %>%
  head(10)

lowest_genes_ratio <- significant_genes_ratio %>%
  arrange(`ratio (pos / neg)`) %>%
  head(10)

# Combine top and lowest genes
label_genes_ratio <- bind_rows(top_genes_ratio, lowest_genes_ratio)

# Create volcano plot based on ratio (title = "Volcano Plot (T-test Results)")
volcano_plot_ratio <- ggplot(test_results_df_bygroup_ttest, aes(x = `ratio (pos / neg)`, y = -log10(p_value))) +
  geom_point(aes(color = p_value < 0.05), alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("black", "red")) +
  geom_vline(xintercept = c(1/1.3, 1.3), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "blue") +
  geom_text_repel(data = label_genes_ratio, aes(label = gene), size = 3) +
  labs(x = "Ratio (pos / neg)", y = "-Log10 P-Value") +
  theme_minimal()

# Save the plot as a PNG file with 300 DPI
ggsave("Figure_2A_colcano.png", plot = volcano_plot_ratio, width = 10, height = 8, units = "in", dpi = 300)

# Print the plot
print(volcano_plot_ratio)

# Source notebook line 2769
test_results_df_bygroup_ttest_copy <- test_results_df_bygroup_ttest

# Update the significant_genes_ratio data frame with new labels
test_results_df_bygroup_ttest_copy$Methylation <- "NO"
test_results_df_bygroup_ttest_copy$Methylation[test_results_df_bygroup_ttest_copy$`ratio (pos / neg)` > 1.3] <- "UP"
test_results_df_bygroup_ttest_copy$Methylation[test_results_df_bygroup_ttest_copy$`ratio (pos / neg)` < 1/1.3] <- "DOWN"

# Set colors for the new labels
mycolors <- c("green", "red", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")

# Identify significant genes based on ratio
significant_genes_ratio <- test_results_df_bygroup_ttest_copy %>%
  filter(`ratio (pos / neg)` > 1.3 | `ratio (pos / neg)` < 1/1.3)

# Select the top 10 largest and lowest genes based on ratio
top_genes_ratio <- significant_genes_ratio %>%
  arrange(desc(`ratio (pos / neg)`)) %>%
  head(10)

lowest_genes_ratio <- significant_genes_ratio %>%
  arrange(`ratio (pos / neg)`) %>%
  head(10)

# Combine top and lowest genes
label_genes_ratio <- bind_rows(top_genes_ratio, lowest_genes_ratio)

# Update the delabel column for the genes to be labeled
test_results_df_bygroup_ttest_copy$delabel <- NA
test_results_df_bygroup_ttest_copy$delabel[test_results_df_bygroup_ttest_copy$Methylation != "NO"] <- 
  test_results_df_bygroup_ttest_copy$gene[test_results_df_bygroup_ttest_copy$Methylation != "NO"]

# Create the volcano plot with updated labeling and coloring
volcano_plot_ratio <- ggplot(test_results_df_bygroup_ttest_copy, aes(x = `ratio (pos / neg)`, y = -log10(p_value), col = Methylation)) +
  geom_point(alpha = 0.8, size = 1.5) +
  geom_vline(xintercept = c(1/1.3, 1.3), col = "blue", linetype = "dashed") +
  geom_hline(yintercept = 1.3, col = "blue", linetype = "dashed") +
  geom_text_repel(data = label_genes_ratio, aes(label = gene), size = 3, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = mycolors, breaks = c("DOWN", "UP", "NO")) +
  labs(x = "Ratio (pos / neg)", y = "-Log10 P-Value", color = "Methylation") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

# Save the plot as a PNG file with 300 DPI
ggsave("Figure_2A_volcano_new.png", plot = volcano_plot_ratio, width = 10, height = 8, units = "in", dpi = 300)

# Print the plot
print(volcano_plot_ratio)

# Source notebook line 2829
# Identify significant genes based on FC
significant_genes_fc <- test_results_df_bygroup_ttest %>%
  filter(FC < -1.3 | FC > 1.3)

# Select the top 10 largest and lowest genes based on FC
top_genes_fc <- significant_genes_fc %>%
  arrange(desc(FC)) %>%
  head(10)

lowest_genes_fc <- significant_genes_fc %>%
  arrange(FC) %>%
  head(10)

# Combine top and lowest genes
label_genes_fc <- bind_rows(top_genes_fc, lowest_genes_fc)

# Create volcano plot based on FC
volcano_plot_fc <- ggplot(test_results_df_bygroup_ttest, aes(x = FC, y = -log10(p_value))) +
  geom_point(aes(color = p_value < 0.05), alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("black", "red")) +
  geom_vline(xintercept = c(-1.3, 1.3), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "blue") +
  geom_text_repel(data = label_genes_fc, aes(label = gene), size = 3) +
  labs(title = "Volcano Plot (T-test Results)", x = "Fold Change (FC)", y = "-Log10 P-Value") +
  theme_minimal()

print(volcano_plot_fc)

# Source notebook line 2862
# Determine genes outside of the vertical lines based on FDR
significant_genes_fdr <- test_results_df_bygroup_ttest %>%
  filter(FDR < 0.001)

# Create the volcano plot based on FDR
volcano_plot_fdr <- ggplot(test_results_df_bygroup_ttest, aes(x = `difference_log2_scale (pos - neg)`, y = -log10(FDR))) +
  geom_point(aes(color = FDR < 0.001), alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("black", "red")) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.001), linetype = "dashed", color = "blue") +
  labs(title = "Volcano Plot (FDR < 0.001)", x = "Difference_log2_scale (pos - neg)", y = "-Log10 FDR") +
  theme_minimal()

# Print the plot
print(volcano_plot_fdr)

# Source notebook line 2886
methyl_gene_cleaned_bygroup_ori <- open_dataset("methyl_cleaned_bygroup_genes.parquet")
methyl_gene_cleaned_bygroup_ori <- as.data.frame(methyl_gene_cleaned_bygroup_ori)

methyl_gene_cleaned_bygroup_ori <- methyl_gene_cleaned_bygroup_ori[, -ncol(methyl_gene_cleaned_bygroup_ori)]

# Source notebook line 2895
# Assuming `test_results_df_bygroup_ttest` is the test result data frame and `beta_value_data` is the beta-value data frame
test_results_df_bygroup_ttest <- as.data.frame(test_results_df_bygroup_ttest)
methyl_gene_cleaned_bygroup_ori <- as.data.frame(methyl_gene_cleaned_bygroup_ori)

# Set row names as gene names
rownames(test_results_df_bygroup_ttest) <- test_results_df_bygroup_ttest$gene
rownames(methyl_gene_cleaned_bygroup_ori) <- methyl_gene_cleaned_bygroup_ori$Gene_Symbol

# Remove the gene columns after setting them as row names
test_results_df_bygroup_ttest <- test_results_df_bygroup_ttest[, -1]
methyl_gene_cleaned_bygroup_ori <- methyl_gene_cleaned_bygroup_ori[, -1]

# Sort both data frames by row names (gene names)
test_results_df_bygroup_ttest <- test_results_df_bygroup_ttest[order(rownames(test_results_df_bygroup_ttest)), ]
methyl_gene_cleaned_bygroup_ori <- methyl_gene_cleaned_bygroup_ori[order(rownames(methyl_gene_cleaned_bygroup_ori)), ]

# Check if all row names are the same and in the same order
rownames_match <- all(rownames(test_results_df_bygroup_ttest) == rownames(methyl_gene_cleaned_bygroup_ori))

# If row names match, merge the data frames
if (rownames_match) {
  # Merge the data frames by row names
  merged_beta_ttest_table <- cbind(test_results_df_bygroup_ttest, methyl_gene_cleaned_bygroup_ori)
  print("The data frames have been merged successfully.")
} else {
  print("The row names do not match.")
}

# Sort the merged data by the `difference` column for clear heatmap visualization
merged_beta_ttest_table <- merged_beta_ttest_table[order(merged_beta_ttest_table$`difference_log2_scale (pos - neg)`, decreasing = TRUE), ]

# Save the merged data to a CSV file
write.csv(merged_beta_ttest_table, "merged_beta_ttest_table.csv", row.names = TRUE)

# Source notebook line 2938
# Get the list of significant gene names
significant_gene_names <- significant_genes_ratio$gene

# Filter the merged data to keep only the significant genes
merged_beta_ttest_table_filtered <- merged_beta_ttest_table[rownames(merged_beta_ttest_table) %in% significant_gene_names, ]

# Ensure all p-values are less than 0.05
merged_beta_ttest_table_filtered <- merged_beta_ttest_table_filtered[merged_beta_ttest_table_filtered$p_value < 0.05, ]

write.csv(merged_beta_ttest_table_filtered, "merged_beta_ttest_table_filtered.csv", row.names = TRUE)

# Source notebook line 2956
table2 <- merged_beta_ttest_table_filtered[, c(1,4,5)]

table2$Gene_Symbol <- rownames(table2)
rownames(table2) <- NULL

# Rearrange columns to have Gene_Symbol as the first column
table2 <- table2[, c("Gene_Symbol", "p_value", "FC", "FDR")]

# Adjust the columns' formats
# table2$p_value <- formatC(table2$p_value, format = "e", digits = 2)
# table2$FDR <- formatC(table2$FDR, format = "e", digits = 2)
# table2$FC <- formatC(table2$FC, format = "f", digits = 2)

write_xlsx(table2, "93_genes_table.xlsx")

largest_FC <- table2[order(-table2$FC), ][1:10, ]
smallest_FC <- table2[order(table2$FC), ][1:10, ]

# Flip the smallest_FC to be in ascending order
smallest_FC <- smallest_FC[order(-smallest_FC$FC), ]

# Adjust the columns' formats
largest_FC$p_value <- formatC(largest_FC$p_value, format = "e", digits = 2)
largest_FC$FDR <- formatC(largest_FC$FDR, format = "e", digits = 2)
largest_FC$FC <- formatC(largest_FC$FC, format = "f", digits = 2)

smallest_FC$p_value <- formatC(smallest_FC$p_value, format = "e", digits = 2)
smallest_FC$FDR <- formatC(smallest_FC$FDR, format = "e", digits = 2)
smallest_FC$FC <- formatC(smallest_FC$FC, format = "f", digits = 2)

# Combine the two subsets and add '...' in the middle
middle_row <- data.frame(Gene_Symbol = '...', p_value = '...', FC = '...', FDR = '...')
table2_new <- rbind(largest_FC, middle_row, smallest_FC)

table_grob <- tableGrob(table2_new, rows = NULL)

# Save the table as a PNG file
png(filename = "table2.png", width = 10, height = 12, units = "in", res = 100)
grid.draw(table_grob)
dev.off()

# Source notebook line 3010
# HPV+ group
hpv_positive <- clinical_focused_filtered_for_age %>% filter(hpv_status == "positive")

# HPV- group
hpv_negative <- clinical_focused_filtered_for_age %>% filter(hpv_status == "negative")


hpv_pos_patient <- rownames(hpv_positive)
hpv_neg_patient <- rownames(hpv_negative)

# Extract the patient columns from the merged data
patient_columns <- colnames(merged_beta_ttest_table)[6:ncol(merged_beta_ttest_table)]

# Verify if all patient columns are accounted for
all_match <- all(patient_columns %in% c(hpv_pos_patient, hpv_neg_patient))

if (all_match) {
  # Rearrange the columns: positive patients first, followed by negative patients
  ordered_patient_columns <- c(hpv_pos_patient, hpv_neg_patient)
  merged_beta_ttest_table_grouped <- merged_beta_ttest_table %>%
    dplyr::select(p_value, `difference_log2_scale (pos - neg)`, `ratio (pos / neg)`, FC, FDR, all_of(ordered_patient_columns))
  print("Patients are grouped.")
} else {
  print("Patient ID not match.")
}

write.csv(merged_beta_ttest_table_grouped, "merged_beta_ttest_table.csv", row.names = TRUE)

# Source notebook line 3043
# Extract the patient columns from the merged data
patient_columns <- colnames(merged_beta_ttest_table_filtered)[6:ncol(merged_beta_ttest_table_filtered)]

# Verify if all patient columns are accounted for
all_match <- all(patient_columns %in% c(hpv_pos_patient, hpv_neg_patient))

if (all_match) {
  # Rearrange the columns: positive patients first, followed by negative patients
  ordered_patient_columns <- c(hpv_pos_patient, hpv_neg_patient)
  merged_beta_ttest_table_filtered_grouped <- merged_beta_ttest_table_filtered %>%
    dplyr::select(p_value, `difference_log2_scale (pos - neg)`, `ratio (pos / neg)`, FC, FDR, all_of(ordered_patient_columns))
  print("Patients are grouped.")
} else {
  print("Patient ID not match.")
}

write.csv(merged_beta_ttest_table_filtered_grouped, "merged_beta_ttest_table_filtered.csv", row.names = TRUE)

# Source notebook line 3069
# Create a row indicating HPV status for each patient
hpv_status_row <- c("NA", "NA", "NA", "NA", "NA", 
                    ifelse(colnames(merged_beta_ttest_table_grouped)[6:ncol(merged_beta_ttest_table_grouped)] %in% hpv_pos_patient, "positive", "negative"))

# Combine this row with the existing data frame
merged_beta_ttest_table_grouped_with_status <- rbind(hpv_status_row, merged_beta_ttest_table_grouped)

# Rename the columns to avoid conflict
colnames(merged_beta_ttest_table_grouped_with_status) <- c("p_value", "difference_log2_scale (pos - neg)", "ratio (pos / neg)", "FC", "FDR",
                                                                    colnames(merged_beta_ttest_table_grouped)[6:ncol(merged_beta_ttest_table_grouped)])

rownames(merged_beta_ttest_table_grouped_with_status)[1] <- "HPV_status"

write.csv(merged_beta_ttest_table_grouped_with_status, "merged_beta_ttest_table_grouped_with_status.csv", row.names = TRUE)

# Source notebook line 3088
# Create a row indicating HPV status for each patient
hpv_status_row <- c("NA", "NA", "NA", "NA", "NA", 
                    ifelse(colnames(merged_beta_ttest_table_filtered_grouped)[6:ncol(merged_beta_ttest_table_filtered_grouped)] %in% hpv_pos_patient, "positive", "negative"))

# Combine this row with the existing data frame
merged_beta_ttest_table_filtered_grouped_with_status <- rbind(hpv_status_row, merged_beta_ttest_table_filtered_grouped)

# Rename the columns to avoid conflict
colnames(merged_beta_ttest_table_filtered_grouped_with_status) <- c("p_value", "difference_log2_scale (pos - neg)", "ratio (pos / neg)", "FC", "FDR",
                                                                    colnames(merged_beta_ttest_table_filtered_grouped)[6:ncol(merged_beta_ttest_table_filtered_grouped)])

rownames(merged_beta_ttest_table_filtered_grouped_with_status)[1] <- "HPV_status"

write.csv(merged_beta_ttest_table_filtered_grouped_with_status, "merged_beta_ttest_table_filtered_grouped_with_status.csv", row.names = TRUE)

# Source notebook line 3111
# Convert the data to numeric format for heatmap plotting
heatmap_data <- as.matrix(merged_beta_ttest_table_filtered_grouped_with_status[-1, -(1:5)])  # Exclude the first row (HPV status) and the first 5 columns (metadata)
heatmap_data <- apply(heatmap_data, 2, as.numeric)

# Create a data frame for annotation
annotation_col <- data.frame(
  HPV_status = factor(ifelse(colnames(heatmap_data) %in% hpv_pos_patient, "positive", "negative"))
)
rownames(annotation_col) <- colnames(heatmap_data)

# Plot the heatmap
pheatmap(heatmap_data, 
         annotation_col = annotation_col, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE, 
         show_colnames = FALSE, 
         scale = "row", 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         annotation_colors = list(HPV_status = c("positive" = "red", "negative" = "blue")),
         main = "Gene Methylation by HPV Status",
         fontsize_row = 6)

# Source notebook line 3141
# Set the row labels for the heatmap
row_labels <- rownames(merged_beta_ttest_table_filtered_grouped_with_status)[-1]

# Save the heatmap as a high-quality PNG file
png("high_quality_heatmap.png", width = 3600, height = 2400, res = 300)

# Plot the heatmap
pheatmap(heatmap_data, 
         annotation_col = annotation_col, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE, 
         show_colnames = FALSE, 
         scale = "row", 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         annotation_colors = list(HPV_status = c("positive" = "red", "negative" = "blue")),
         #main = "Gene Methylation by HPV Status",
         fontsize_row = 8,  # Adjusted for smaller font size
         labels_row = row_labels  # Ensures that row labels are displayed
)

# Close the PNG device
dev.off()

