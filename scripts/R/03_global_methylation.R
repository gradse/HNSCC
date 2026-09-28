# Analysis stage from archive/HNSCC_Analysis.Rmd (source lines 503-952).
# Run through scripts/run_r_analysis.R.

# Source notebook line 503
methyl_merge_calculated <- open_dataset("cleaned_methylation_data_final.parquet")

methyl_merge_calculated <- as.data.frame(methyl_merge_calculated)

sorted_clinical_focused <- clinical_focused[order(rownames(clinical_focused)), ]

methyl_sum_row <- fread("first_row_extracted.csv")

methyl_sum_row <- t(methyl_sum_row)

colnames(methyl_sum_row) <- "Summation"

sorted_methyl_sum <- methyl_sum_row[order(rownames(methyl_sum_row)), ]

sorted_methyl_sum <- as.data.frame(sorted_methyl_sum)

colnames(sorted_methyl_sum) <- "Summation"

same_names <- all(rownames(sorted_clinical_focused) == rownames(sorted_methyl_sum))

if (same_names) {
  clinical_focused_combined <- cbind(sorted_methyl_sum, sorted_clinical_focused)
} else {
  print("The row names of the datasets are not the same or not in the same order.")
}

# Source notebook line 532
# Convert the Summation column to numeric
clinical_focused_combined$Summation <- as.numeric(clinical_focused_combined$Summation)

# Check for any conversion warnings
if (any(is.na(clinical_focused_combined$Summation))) {
  warning("There are NA values in the Summation column.")
}

# Extract patient IDs where hpv_status is "indeterminate"
indeterminate_ids <- clinical_focused_combined %>% 
  filter(hpv_status == "indeterminate") %>% 
  rownames()

# Remove indeterminate value (2 rows removed)
clinical_focused_combined <- clinical_focused_combined %>% 
  filter(hpv_status != "indeterminate")

# Create the box plot
ggplot(clinical_focused_combined, aes(x = hpv_status, y = Summation)) +
  geom_boxplot() +
  labs(title = "Box Plot of Beta_value Summation by HPV Status",
       x = "HPV Status",
       y = "Beta-value Summation") +
  theme_minimal()

#write_parquet(clinical_focused_combined, "clinical_focused_combined_summation.parquet")

# Save the plot
ggsave("boxplot_summation_by_hpv_status_before_removal.png")

# Source notebook line 565
# Perform t-test between HPV positive and HPV negative groups
t_test_result <- t.test(Summation ~ hpv_status, data = clinical_focused_combined)

print(t_test_result)

# Source notebook line 573
wilcox_test_before <- wilcox.test(Summation ~ hpv_status, data = clinical_focused_combined)

print(wilcox_test_before)

# Source notebook line 582
# Remove columns with IDs in indeterminate_ids
methyl_merge_calculated <- methyl_merge_calculated %>%
  select(-all_of(indeterminate_ids))

# Save the cleaned dataframe to a new parquet file
write_parquet(methyl_merge_calculated, "methyl_calculated_removed_indeter_hpv.parquet")

# Source notebook line 623
# This is the data that have been cleaned for both hpv "indeterminate" and "nan" rows
methyl_merge_calculated_cleaned <- open_dataset("methyl_calculated_removed_indeter_hpv_and_missing_beta.parquet")
methyl_merge_calculated_cleaned <- as.data.frame(methyl_merge_calculated_cleaned)

# Source notebook line 657
methyl_sub <- open_dataset("subdataset_for_15_perc.parquet")
methyl_sub <- as.data.frame(methyl_sub)

# Source notebook line 666
# import pandas as pd
# import numpy as np
# 
# # Define the path to the Parquet file
# input_file = "subdataset_for_15_perc.parquet"
# 
# # Read the Parquet file into a DataFrame
# print("Reading Parquet file...")
# df = pd.read_parquet(input_file)
# 
# # Convert all values to numeric, coerce errors to NaN
# print("Converting values to numeric...")
# df = df.apply(pd.to_numeric, errors='coerce')
# 
# # Flatten the DataFrame to a single array, ignoring NaN values
# all_values = df.values.flatten()
# all_values = all_values[~np.isnan(all_values)]
# 
# # Calculate the 15th percentile value
# percentile_15th = np.percentile(all_values, 15)
# 
# print(f"The 15th percentile value is: {percentile_15th}")

# Source notebook line 736
methyl_calculated_patient <- open_dataset("methy_calculated_patient.parquet")
methyl_calculated_patient <- as.data.frame(methyl_calculated_patient)

summary(methyl_calculated_patient$Patients_Beta_Value_Count >= 447)

# Source notebook line 748
# import pandas as pd
# import numpy as np
# 
# # Define the input and output file paths
# input_file = "methy_calculated_patient.parquet"
# output_file = "methyl_removed_genes.parquet"
# 
# # Read the Parquet file into a DataFrame
# print("Reading Parquet file...")
# df = pd.read_parquet(input_file)
# 
# # Replace string 'nan' with actual NaN values
# df.replace('nan', np.nan, inplace=True)
# 
# # Preserve the first row (Summation)
# summation_row = df.iloc[0].copy()
# 
# # Remove rows where Patients_Beta_Value_Count >= 447, excluding the first row
# print("Removing rows where Patients_Beta_Value_Count >= 447...")
# df = df.iloc[1:]
# df = df[df["Patients_Beta_Value_Count"] < 447]
# 
# # Reinsert the Summation row at the top
# df = pd.concat([summation_row.to_frame().T, df]).reset_index(drop=True)
# 
# # Recalculate the first row "Summation"
# print("Recalculating the Summation row...")
# 
# # Specify the columns for summation (4th, 5th, and 7th onward are patients)
# columns_to_summate = ['Low', 'High'] + list(df.columns[6:])
# 
# # Convert the specified columns to numeric for summation, ignoring NaN values
# summation_values = df.loc[1:, columns_to_summate].apply(pd.to_numeric, errors='coerce').sum(axis=0, skipna=True)
# 
# # Convert the summation values back to strings with full decimals
# summation_values_str = summation_values.apply(lambda x: f"{x:.16f}")
# 
# # Update the Summation row in the DataFrame
# df.loc[0, columns_to_summate] = summation_values_str
# 
# # Find the Low and High values for the Summation row
# low_value = summation_values.min()
# high_value = summation_values.max()
# 
# # Update the Low and High columns in the Summation row
# df.loc[0, 'Low'] = f"{low_value:.16f}"
# df.loc[0, 'High'] = f"{high_value:.16f}"
# 
# # Save the updated DataFrame to a new Parquet file
# print("Saving the updated DataFrame to a new Parquet file...")
# df.to_parquet(output_file, index=False)
# 
# print("All operations completed and new Parquet file saved.")

# methyl_remove_genes.py:24: FutureWarning: The behavior of DataFrame concatenation with empty or all-NA entries is deprecated. In a future version, this will no longer exclude empty or all-NA columns when determining the result dtypes. To retain the old behavior, exclude the relevant entries before the concat operation.
#   df = pd.concat([summation_row.to_frame().T, df]).reset_index(drop=True)

# Source notebook line 811
methyl_removed_genes <- open_dataset("methyl_removed_genes.parquet")
methyl_removed_genes <- as.data.frame(methyl_removed_genes)

# Ready for binding new summation values
clinical_focused_combined_new <- clinical_focused_combined %>% 
  select(-Summation)

clinical_focused_combined_new <- clinical_focused_combined_new[order(rownames(clinical_focused_combined_new)), ]

# Read the CSV file
methyl_sum_row_new <- fread("extracted_first_row_new.csv", stringsAsFactors = FALSE)

# Transpose the data
methyl_sum_row_new <- t(methyl_sum_row_new)

# Set column name
colnames(methyl_sum_row_new) <- "Summation"

# Sort rows by their row names
sorted_methyl_sum_new <- methyl_sum_row_new[order(rownames(methyl_sum_row_new)), ]

# Convert to data frame
sorted_methyl_sum_new <- as.data.frame(sorted_methyl_sum_new, stringsAsFactors = FALSE)

# Set column name
colnames(sorted_methyl_sum_new) <- "Summation"

# Ensure the Summation column retains all decimal places
sorted_methyl_sum_new$Summation <- format(as.numeric(sorted_methyl_sum_new$Summation), digits = 16, scientific = FALSE)

same_names_new <- all(rownames(clinical_focused_combined_new) == rownames(sorted_methyl_sum_new))

if (same_names_new) {
  clinical_focused_combined_new <- cbind(sorted_methyl_sum_new, clinical_focused_combined_new)
} else {
  print("The row names of the datasets are not the same or not in the same order.")
}

clinical_focused_combined_new$Summation <- as.numeric(clinical_focused_combined_new$Summation)

# Create the box plot
ggplot(clinical_focused_combined_new, aes(x = hpv_status, y = Summation)) +
  geom_boxplot() +
  labs(title = "Box Plot of Beta_value Summation by HPV Status",
       x = "HPV Status",
       y = "Beta-value Summation") +
  theme_minimal()

# # Create the box plot with jittered dots
# ggplot(clinical_focused_combined_new, aes(x = hpv_status, y = Summation)) +
#   geom_boxplot(outlier.shape = NA) +  # Hide the default outliers to avoid overlapping with jittered points
#   geom_jitter(width = 0.2, alpha = 0.5, color = "blue") +
#   labs(title = "Box Plot of Beta-value Summation by HPV Status",
#        x = "HPV Status",
#        y = "Beta-value Summation") +
#   theme_minimal()

# Save the plot
ggsave("boxplot_summation_by_hpv_status_after_removal.png")

# Save the cleaned dataframe to a new csv file
write_csv(clinical_focused_combined_new, "clinical_focused_combined_new.csv")

# Source notebook line 878
clinical_ttest <- fread("clinical_focused_combined_new.csv")
clinical_ttest <- as.data.frame(clinical_ttest)

# Perform t-test between HPV positive and HPV negative groups
t_test_result_new <- t.test(Summation ~ hpv_status, data = clinical_ttest)

print(t_test_result_new)

# Source notebook line 890
wilcox_test <- wilcox.test(Summation ~ hpv_status, data = clinical_ttest, exact = FALSE)

print(wilcox_test)

# Source notebook line 897
subdivision_new <- as.factor(clinical_ttest$anatomic_neoplasm_subdivision)

is_hpv_new <- as.factor(clinical_ttest$hpv_status)

table1_new <- table(subdivision_new, is_hpv_new)

table1_new_margin <- addmargins(table1_new)

kable(table1_new_margin, format = "html",
      caption = "Cross-Tabulation of Anatomic Neoplasm Subdivision and HPV Status") %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                    full_width = FALSE)

# Source notebook line 918
# Filter the data for anatomic_neoplasm_subdivision == "tonsil"
clinical_ttest_filtered <- clinical_ttest %>%
  filter(anatomic_neoplasm_subdivision == "tonsil")

# Check for any conversion warnings
if (any(is.na(clinical_ttest_filtered$Summation))) {
  warning("There are NA values in the Summation column.")
}

# Perform t-test
t_test_result_tonsil <- t.test(Summation ~ hpv_status, data = clinical_ttest_filtered)

# Print the test result
print(t_test_result_tonsil)

# Source notebook line 936
# Create the box plot with jittered dots
ggplot(clinical_ttest_filtered, aes(x = hpv_status, y = Summation)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(width = 0.2, alpha = 0.5, color = "blue") +  
  labs(title = "Box Plot of Beta-value Summation by HPV Status",
       x = "HPV Status",
       y = "Beta-value Summation") +
  theme_minimal()

