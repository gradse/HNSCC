# Analysis stage from archive/HNSCC_Analysis.Rmd (source lines 1035-1732).
# Run through scripts/run_r_analysis.R.

# Source notebook line 1035
# BiocManager::install("doParallel")
# BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

# Source notebook line 1048
data(age_coefficients)

horvath_cpgs <- setdiff(names(ageCoefs$Horvath), "(Intercept)")
hannum_cpgs <- names(ageCoefs$Hannum)
pheno_cpgs <- setdiff(names(ageCoefs$PhenoAge), "(Intercept)")

write.csv(horvath_cpgs, file = "horvath_cpgs.csv", row.names = FALSE)
write.csv(hannum_cpgs, file = "hannum_cpgs.csv", row.names = FALSE)
write.csv(pheno_cpgs, file = "pheno_cpgs.csv", row.names = FALSE)

# Source notebook line 1116
horvath_data <- read_parquet("filtered_methylation_data_horvath.parquet")
hannum_data <- read_parquet("filtered_methylation_data_hannum.parquet")
phenoage_data <- read_parquet("filtered_methylation_data_phenoage.parquet")

# Delete the last column
horvath_data <- horvath_data[, -ncol(horvath_data)]
hannum_data <- hannum_data[, -ncol(hannum_data)]
phenoage_data <- phenoage_data[, -ncol(phenoage_data)]


# Function to convert data frame to numeric matrix
convert_to_numeric_matrix <- function(data) {
  if (ncol(data) < 2) {
    stop("Data must have at least two columns")
  }
  
  cpg_ids <- data[[1]]  # Extract the first column as CpG identifiers
  data <- data[, -1]  # Exclude the first column
  
  # Check for unique CpG identifiers
  if (length(unique(cpg_ids)) != length(cpg_ids)) {
    stop("CpG identifiers are not unique")
  }
  
  # Ensure all values are numeric while preserving row names
  data_matrix <- apply(data, 2, function(x) as.numeric(as.character(x)))
  rownames(data_matrix) <- cpg_ids
  
  if (length(rownames(data_matrix)) != nrow(data_matrix)) {
    print(paste("Number of row names:", length(rownames(data_matrix))))
    print(paste("Number of rows:", nrow(data_matrix)))
    stop("Mismatch between number of row names and number of rows")
  }
  return(data_matrix)
}

# Convert data to numeric matrices
horvath_betas <- convert_to_numeric_matrix(horvath_data)
hannum_betas <- convert_to_numeric_matrix(hannum_data)
phenoage_betas <- convert_to_numeric_matrix(phenoage_data)

# Check if the matrices contain any non-numeric values
check_numeric <- function(matrix_data, name) {
  if (!is.matrix(matrix_data) || mode(matrix_data) != "numeric") {
    stop(paste("Matrix for", name, "is not numeric"))
  }
}

check_numeric(horvath_betas, "Horvath")
check_numeric(hannum_betas, "Hannum")
check_numeric(phenoage_betas, "PhenoAge")

# Calculate the epigenetic clocks
horvath_age <- agep(horvath_betas, method = "horvath")
hannum_age <- agep(hannum_betas, method = "hannum")
phenoage_age <- agep(phenoage_betas, method = "phenoage")

# Check if row names are the same
same_rownames <- all(rownames(horvath_age) == rownames(hannum_age)) && all(rownames(horvath_age) == rownames(phenoage_age))

if (same_rownames) {
  # Combine results into a data frame with rownames
  clock_ages <- data.frame(Horvath = horvath_age[,1], Hannum = hannum_age[,1], PhenoAge = phenoage_age[,1])
  rownames(clock_ages) <- rownames(horvath_age)
  print(clock_ages)
} else {
  print("Row names do not match between the age calculations.")
}

# Source notebook line 1191
ori_data <- read_parquet("methyl_calculated_removed_indeter_hpv.parquet")
ori_data <- as.data.frame(ori_data)

# Source notebook line 1274
horvath_data_new <- read_parquet("filtered_methyl_horvath_new.parquet")
hannum_data_new <- read_parquet("filtered_methyl_hannum_new.parquet")
phenoage_data_new <- read_parquet("filtered_methyl_phenoage_new.parquet")

# Delete the last column
horvath_data_new <- horvath_data_new[, -ncol(horvath_data_new)]
hannum_data_new <- hannum_data_new[, -ncol(hannum_data_new)]
phenoage_data_new <- phenoage_data_new[, -ncol(phenoage_data_new)]


# Function to convert data frame to numeric matrix
convert_to_numeric_matrix <- function(data) {
  if (ncol(data) < 2) {
    stop("Data must have at least two columns")
  }
  
  cpg_ids <- data[[1]]  # Extract the first column as CpG identifiers
  data <- data[, -1]  # Exclude the first column
  
  # Check for unique CpG identifiers
  if (length(unique(cpg_ids)) != length(cpg_ids)) {
    stop("CpG identifiers are not unique")
  }
  
  # Ensure all values are numeric while preserving row names
  data_matrix <- apply(data, 2, function(x) as.numeric(as.character(x)))
  rownames(data_matrix) <- cpg_ids
  
  if (length(rownames(data_matrix)) != nrow(data_matrix)) {
    print(paste("Number of row names:", length(rownames(data_matrix))))
    print(paste("Number of rows:", nrow(data_matrix)))
    stop("Mismatch between number of row names and number of rows")
  }
  
  
  return(data_matrix)
}

# Convert data to numeric matrices
horvath_betas_new <- convert_to_numeric_matrix(horvath_data_new)
hannum_betas_new <- convert_to_numeric_matrix(hannum_data_new)
phenoage_betas_new <- convert_to_numeric_matrix(phenoage_data_new)

# Check if the matrices contain any non-numeric values
check_numeric <- function(matrix_data, name) {
  if (!is.matrix(matrix_data) || mode(matrix_data) != "numeric") {
    stop(paste("Matrix for", name, "is not numeric"))
  }
}

check_numeric(horvath_betas_new, "Horvath")
check_numeric(hannum_betas_new, "Hannum")
check_numeric(phenoage_betas_new, "PhenoAge")

# Calculate Horvath Clock
# horvath_age_new <- agep(horvath_betas_new, method = "horvath", n_missing = TRUE)
# print(horvath_age_new)

# Calculate Hannum Clock
# hannum_age_new <- agep(hannum_betas_new, method = "hannum", n_missing = TRUE)
# print(hannum_age_new)

# Calculate PhenoAge
# phenoage_new <- agep(phenoage_betas_new, method = "phenoage", n_missing = TRUE)
# print(phenoage_new)

# Calculate the epigenetic clocks
horvath_age <- agep(horvath_betas_new, method = "horvath")
hannum_age <- agep(hannum_betas_new, method = "hannum")
phenoage_age <- agep(phenoage_betas_new, method = "phenoage")

# Check if row names are the same
same_rownames <- all(rownames(horvath_age) == rownames(hannum_age)) && all(rownames(horvath_age) == rownames(phenoage_age))

if (same_rownames) {
  # Combine results into a data frame with rownames
  clock_ages <- data.frame(Horvath = horvath_age[,1], Hannum = hannum_age[,1], PhenoAge = phenoage_age[,1])
  rownames(clock_ages) <- rownames(horvath_age)
  print(clock_ages)
} else {
  print("Row names do not match between the age calculations.")
}

# Source notebook line 1365
new_sum_row <- read.csv("new_summation_row.csv")

methyl_imputed_new <- read_parquet("methyl_imputed_sum.parquet")
methyl_imputed_new <- as.data.frame(methyl_imputed_new)

# Source notebook line 1374
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

# Source notebook line 1395
# Process the summation row
rownames(new_sum_row) <- new_sum_row[,1]
new_sum_row <- new_sum_row[ , -1]

# Transpose the data
new_sum_row_transpose <- t(new_sum_row)

# Convert to data.frame
new_sum_row_transpose <- as.data.frame(new_sum_row_transpose)

# Sort the row names
sorted_new_sum_row_transpose <- new_sum_row_transpose[order(rownames(new_sum_row_transpose)), , drop = FALSE]
rownames(sorted_new_sum_row_transpose) <- gsub("\\.", "-", rownames(sorted_new_sum_row_transpose))
row_name_sum <- rownames(sorted_new_sum_row_transpose)

# Process clocks data (sort rownames)
sorted_clock_ages <- clock_ages[order(rownames(clock_ages)), , drop = FALSE]
row_name_clock <- rownames(sorted_clock_ages)


# Process clinical data (sort rownames)
sorted_clinical_focused_filtered_for_age <- clinical_focused_filtered_for_age[order(rownames(clinical_focused_filtered_for_age)), , drop = FALSE]
row_name_clinical <- rownames(sorted_clinical_focused_filtered_for_age)



# Merge them into clinical data
all_colnames_match <- all(row_name_sum == row_name_clock)

if (all_colnames_match) {
  # Merge the data frames
  merged_df <- cbind(sorted_new_sum_row_transpose, sorted_clock_ages)
  
  # Further merge with the clinical data
  all_colnames_match_clinical <- all(rownames(merged_df) == row_name_clinical)
  
  if (all_colnames_match_clinical) {
    sorted_clinical_focused_filtered_for_age <- cbind(merged_df, clinical_focused_filtered_for_age)
    
    # Save the merged data frame
    write.csv(sorted_clinical_focused_filtered_for_age, "sorted_clinical_focused_filtered_for_age.csv", row.names = TRUE)
    print("Merged data frame saved successfully.")
  } else {
    print("The row names of the merged data frame and clinical data frame do not match.")
  }
} else {
  print("The column names of the summation row transposed and clock ages data frames do not match.")
}

# Source notebook line 1453
# Filter out age is missing
clinical_filtered_age <- sorted_clinical_focused_filtered_for_age %>% 
  filter(!is.na(age_at_initial_pathologic_diagnosis))

# Ensure age_at_initial_pathologic_diagnosis is numeric
clinical_filtered_age$age_at_initial_pathologic_diagnosis <- as.numeric(clinical_filtered_age$age_at_initial_pathologic_diagnosis)

# Calculate the differences and add them as new columns
clinical_filtered_age$Horvath_diff <- abs(clinical_filtered_age$Horvath - clinical_filtered_age$age_at_initial_pathologic_diagnosis)
clinical_filtered_age$Hannum_diff <- abs(clinical_filtered_age$Hannum - clinical_filtered_age$age_at_initial_pathologic_diagnosis)
clinical_filtered_age$PhenoAge_diff <- abs(clinical_filtered_age$PhenoAge - clinical_filtered_age$age_at_initial_pathologic_diagnosis)

# Reorder columns to place the new columns in the 5th, 6th, and 7th positions
column_order <- c(names(clinical_filtered_age)[1:4], 'Horvath_diff', 'Hannum_diff', 'PhenoAge_diff', names(clinical_filtered_age)[5:(ncol(clinical_filtered_age)-3)])
clinical_filtered_age <- clinical_filtered_age[, column_order]


# Source notebook line 1479
# Create the boxplot
ggplot(clinical_filtered_age, aes(x = hpv_status, y = Summation)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  # Add points for individual data
  labs(title = "Box Plot of Summation by HPV Status",
       x = "HPV Status",
       y = "Summation") +
  theme_minimal()

# Check for normality (Shapiro-Wilk test)
shapiro_test <- shapiro.test(clinical_filtered_age$Summation)

if (shapiro_test$p.value < 0.05) {
  # Data is not normally distributed, use Wilcoxon rank-sum test
  wilcox_test <- wilcox.test(Summation ~ hpv_status, data = clinical_filtered_age)
  test_result <- wilcox_test
  test_name <- "Wilcoxon rank-sum test"
} else {
  # Data is normally distributed, use t-test
  t_test <- t.test(Summation ~ hpv_status, data = clinical_filtered_age)
  test_result <- t_test
  test_name <- "t-test"
}

# Print the results of the test
cat("Result of the", test_name, ":\n")
print(test_result)

# Source notebook line 1512
# Create the boxplot
ggplot(clinical_filtered_age, aes(x = hpv_status, y = Horvath_diff)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  # Add points for individual data
  labs(title = "Box Plot of Horvath_diff by HPV Status",
       x = "HPV Status",
       y = "Horvath_diff") +
  theme_minimal()

# Check for normality (Shapiro-Wilk test)
shapiro_test <- shapiro.test(clinical_filtered_age$Horvath_diff)

if (shapiro_test$p.value < 0.05) {
  # Data is not normally distributed, use Wilcoxon rank-sum test
  wilcox_test <- wilcox.test(Horvath_diff ~ hpv_status, data = clinical_filtered_age)
  test_result <- wilcox_test
  test_name <- "Wilcoxon rank-sum test"
} else {
  # Data is normally distributed, use t-test
  t_test <- t.test(Horvath_diff ~ hpv_status, data = clinical_filtered_age)
  test_result <- t_test
  test_name <- "t-test"
}

# Print the results of the test
cat("Result of the", test_name, ":\n")
print(test_result)

# Source notebook line 1543
# Create the boxplot
ggplot(clinical_filtered_age, aes(x = hpv_status, y = Horvath)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  # Add points for individual data
  labs(title = "Box Plot of Horvath by HPV Status",
       x = "HPV Status",
       y = "Horvath_diff") +
  theme_minimal()

# Check for normality (Shapiro-Wilk test)
shapiro_test <- shapiro.test(clinical_filtered_age$Horvath_diff)

if (shapiro_test$p.value < 0.05) {
  # Data is not normally distributed, use Wilcoxon rank-sum test
  wilcox_test <- wilcox.test(Horvath_diff ~ hpv_status, data = clinical_filtered_age)
  test_result <- wilcox_test
  test_name <- "Wilcoxon rank-sum test"
} else {
  # Data is normally distributed, use t-test
  t_test <- t.test(Horvath_diff ~ hpv_status, data = clinical_filtered_age)
  test_result <- t_test
  test_name <- "t-test"
}

# Print the results of the test
cat("Result of the", test_name, ":\n")
print(test_result)

# Source notebook line 1575
# Create the boxplot
ggplot(clinical_filtered_age, aes(x = hpv_status, y = Hannum_diff)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  # Add points for individual data
  labs(title = "Box Plot of Hannum_diff by HPV Status",
       x = "HPV Status",
       y = "Hannum_diff") +
  theme_minimal()

# Check for normality (Shapiro-Wilk test)
shapiro_test <- shapiro.test(clinical_filtered_age$Hannum_diff)

if (shapiro_test$p.value < 0.05) {
  # Data is not normally distributed, use Wilcoxon rank-sum test
  wilcox_test <- wilcox.test(Hannum_diff ~ hpv_status, data = clinical_filtered_age)
  test_result <- wilcox_test
  test_name <- "Wilcoxon rank-sum test"
} else {
  # Data is normally distributed, use t-test
  t_test <- t.test(Hannum_diff ~ hpv_status, data = clinical_filtered_age)
  test_result <- t_test
  test_name <- "t-test"
}

# Print the results of the test
cat("Result of the", test_name, ":\n")
print(test_result)

# Source notebook line 1608
# Create the boxplot
ggplot(clinical_filtered_age, aes(x = hpv_status, y = PhenoAge_diff)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  # Add points for individual data
  labs(title = "Box Plot of PhenoAge_diff by HPV Status",
       x = "HPV Status",
       y = "PhenoAge_diff") +
  theme_minimal()

# Check for normality (Shapiro-Wilk test)
shapiro_test <- shapiro.test(clinical_filtered_age$PhenoAge_diff)

if (shapiro_test$p.value < 0.05) {
  # Data is not normally distributed, use Wilcoxon rank-sum test
  wilcox_test <- wilcox.test(PhenoAge_diff ~ hpv_status, data = clinical_filtered_age)
  test_result <- wilcox_test
  test_name <- "Wilcoxon rank-sum test"
} else {
  # Data is normally distributed, use t-test
  t_test <- t.test(PhenoAge_diff ~ hpv_status, data = clinical_filtered_age)
  test_result <- t_test
  test_name <- "t-test"
}

# Print the results of the test
cat("Result of the", test_name, ":\n")
print(test_result)

# Source notebook line 1643
# Create the boxplot
ggplot(clinical_filtered_age, aes(x = hpv_status, y = age_at_initial_pathologic_diagnosis)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  # Add points for individual data
  labs(title = "Box Plot of age by HPV Status",
       x = "HPV Status",
       y = "Age") +
  theme_minimal()

# Check for normality (Shapiro-Wilk test)
shapiro_test <- shapiro.test(clinical_filtered_age$age_at_initial_pathologic_diagnosis)

if (shapiro_test$p.value < 0.05) {
  # Data is not normally distributed, use Wilcoxon rank-sum test
  wilcox_test <- wilcox.test(age_at_initial_pathologic_diagnosis ~ hpv_status, data = clinical_filtered_age)
  test_result <- wilcox_test
  test_name <- "Wilcoxon rank-sum test"
} else {
  # Data is normally distributed, use t-test
  t_test <- t.test(age_at_initial_pathologic_diagnosis ~ hpv_status, data = clinical_filtered_age)
  test_result <- t_test
  test_name <- "t-test"
}

# Print the results of the test
cat("Result of the", test_name, ":\n")
print(test_result)

# Source notebook line 1674
# HPV+ group
hpv_positive <- clinical_filtered_age %>% filter(hpv_status == "positive")

t_test_hannum_hpv_pos <- t.test(hpv_positive$Hannum, hpv_positive$age_at_initial_pathologic_diagnosis, paired = TRUE)
t_test_phenoage_hpv_pos <- t.test(hpv_positive$PhenoAge, hpv_positive$age_at_initial_pathologic_diagnosis, paired = TRUE)
t_test_horvath_hpv_pos <- t.test(hpv_positive$Horvath, hpv_positive$age_at_initial_pathologic_diagnosis, paired = TRUE)

# HPV- group
hpv_negative <- clinical_filtered_age %>% filter(hpv_status == "negative")

t_test_hannum_hpv_neg <- t.test(hpv_negative$Hannum, hpv_negative$age_at_initial_pathologic_diagnosis, paired = TRUE)
t_test_phenoage_hpv_neg <- t.test(hpv_negative$PhenoAge, hpv_negative$age_at_initial_pathologic_diagnosis, paired = TRUE)
t_test_horvath_hpv_neg <- t.test(hpv_negative$Horvath, hpv_negative$age_at_initial_pathologic_diagnosis, paired = TRUE)

# Print results
print(t_test_hannum_hpv_pos)
print(t_test_phenoage_hpv_pos)
print(t_test_horvath_hpv_pos)

print(t_test_hannum_hpv_neg)
print(t_test_phenoage_hpv_neg)
print(t_test_horvath_hpv_neg)

# Source notebook line 1704
# HPV+ group
wilcox_test_hannum_hpv_pos <- wilcox.test(hpv_positive$Hannum, hpv_positive$age_at_initial_pathologic_diagnosis, paired = TRUE)
wilcox_test_phenoage_hpv_pos <- wilcox.test(hpv_positive$PhenoAge, hpv_positive$age_at_initial_pathologic_diagnosis, paired = TRUE)
wilcox_test_horvath_hpv_pos <- wilcox.test(hpv_positive$Horvath, hpv_positive$age_at_initial_pathologic_diagnosis, paired = TRUE)

# HPV- group
wilcox_test_hannum_hpv_neg <- wilcox.test(hpv_negative$Hannum, hpv_negative$age_at_initial_pathologic_diagnosis, paired = TRUE)
wilcox_test_phenoage_hpv_neg <- wilcox.test(hpv_negative$PhenoAge, hpv_negative$age_at_initial_pathologic_diagnosis, paired = TRUE)
wilcox_test_horvath_hpv_neg <- wilcox.test(hpv_negative$Horvath, hpv_negative$age_at_initial_pathologic_diagnosis, paired = TRUE)

# Print results
print(wilcox_test_hannum_hpv_pos)
print(wilcox_test_phenoage_hpv_pos)
print(wilcox_test_horvath_hpv_pos)

print(wilcox_test_hannum_hpv_neg)
print(wilcox_test_phenoage_hpv_neg)
print(wilcox_test_horvath_hpv_neg)

