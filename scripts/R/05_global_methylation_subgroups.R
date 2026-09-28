# Analysis stage from archive/HNSCC_Analysis.Rmd (source lines 1733-2067).
# Run through scripts/run_r_analysis.R.

# Source notebook line 1733
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

# Source notebook line 1754
# only removed nan rows (396066 left)
methyl_sum_row <- fread("first_row_extracted.csv")
methyl_sum_row <- t(methyl_sum_row)
colnames(methyl_sum_row) <- "Summation"
sorted_methyl_sum <- methyl_sum_row[order(rownames(methyl_sum_row)), ]
sorted_methyl_sum <- as.data.frame(sorted_methyl_sum)
colnames(sorted_methyl_sum) <- "Summation"


# further cleaned (359376 left)
methyl_sum_row_new <- fread("extracted_first_row_new.csv")
methyl_sum_row_new <- t(methyl_sum_row_new)
colnames(methyl_sum_row_new) <- "Summation_cleaned"
sorted_methyl_sum_new <- methyl_sum_row_new[order(rownames(methyl_sum_row_new)), ]
sorted_methyl_sum_new <- as.data.frame(sorted_methyl_sum_new)
colnames(sorted_methyl_sum_new) <- "Summation_cleaned"


# Process clinical data (sort rownames)
sorted_clinical_focused_filtered_for_age <- clinical_focused_filtered_for_age[order(rownames(clinical_focused_filtered_for_age)), , drop = FALSE]
row_name_clinical <- rownames(sorted_clinical_focused_filtered_for_age)


# Merge them

# Remove rows in sorted_methyl_sum that are not in row_name_clinical
filtered_methyl_sum <- sorted_methyl_sum[row_name_clinical, , drop = FALSE]

# Check if the row names are consistent
all_row_names_match <- all(rownames(filtered_methyl_sum) == row_name_clinical) &&
                       all(rownames(sorted_methyl_sum_new) == row_name_clinical)
if (all_row_names_match) {
  final_clinical_data <- cbind(filtered_methyl_sum, sorted_methyl_sum_new, sorted_clinical_focused_filtered_for_age)
  print("The row names are consistent and the data has been merged successfully.")
} else {
  print("The row names do not match. Please check the data.")
}

# Source notebook line 1799
# Create the boxplot (title = "Box Plot of Global Methylation by HPV Status")
figure_1A <- ggplot(final_clinical_data, aes(x = hpv_status, y = Summation)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  
  labs(x = "HPV Status",
       y = "Global Methylation") +
  theme_minimal()

# Save the plot as a PNG file with 300 DPI
ggsave("1A-boxplot_global_methylation_by_hpv_status.png", plot = figure_1A, width = 8, height = 6, dpi = 300)

# Check for normality (Shapiro-Wilk test)
shapiro_test <- shapiro.test(final_clinical_data$Summation)

if (shapiro_test$p.value < 0.05) {
  # Data is not normally distributed, use Wilcoxon rank-sum test
  wilcox_test <- wilcox.test(Summation ~ hpv_status, data = final_clinical_data)
  test_result <- wilcox_test
  test_name <- "Wilcoxon rank-sum test"
} else {
  # Data is normally distributed, use t-test
  t_test <- t.test(Summation ~ hpv_status, data = final_clinical_data)
  test_result <- t_test
  test_name <- "t-test"
}

# Print the results of the test
cat("Result of the", test_name, ":\n")
print(test_result)

# Source notebook line 1833
# Create the boxplot
ggplot(final_clinical_data, aes(x = hpv_status, y = Summation_cleaned)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  
  labs(title = "Box Plot of Summation_cleaned by HPV Status",
       x = "HPV Status",
       y = "Summation_cleaned") +
  theme_minimal()

# Check for normality (Shapiro-Wilk test)
shapiro_test <- shapiro.test(final_clinical_data$Summation_cleaned)

if (shapiro_test$p.value < 0.05) {
  # Data is not normally distributed, use Wilcoxon rank-sum test
  wilcox_test <- wilcox.test(Summation_cleaned ~ hpv_status, data = final_clinical_data)
  test_result <- wilcox_test
  test_name <- "Wilcoxon rank-sum test"
} else {
  # Data is normally distributed, use t-test
  t_test <- t.test(Summation_cleaned ~ hpv_status, data = final_clinical_data)
  test_result <- t_test
  test_name <- "t-test"
}

# Print the results of the test
cat("Result of the", test_name, ":\n")
print(test_result)

# Source notebook line 1868
# male group
male_group <- final_clinical_data %>% filter(gender == "male")

# female group
female_group <- final_clinical_data %>% filter(gender == "female")

# Source notebook line 1879
# Create the boxplot (title = "Box Plot of Global Methylation by HPV Status in Male Group")
figure_1B1 <- ggplot(male_group, aes(x = hpv_status, y = Summation)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  
  labs(x = "HPV Status",
       y = "Global Methylation") +
  theme_minimal()

# Save the plot as a PNG file with 300 DPI
ggsave("1B1-boxplot_global_methylation_by_hpv_status_male_group.png", plot = figure_1B1, width = 8, height = 6, dpi = 300)

# Source notebook line 1895
# Create the boxplot (title = "Box Plot of Global Methylation by HPV Status in Female Group")
figure_1B2 <- ggplot(female_group, aes(x = hpv_status, y = Summation)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  
  labs(x = "HPV Status",
       y = "Global Methylation") +
  theme_minimal()

# Save the plot as a PNG file with 300 DPI
ggsave("1B2-boxplot_global_methylation_by_hpv_status_female_group.png", plot = figure_1B2, width = 8, height = 6, dpi = 300)

# Source notebook line 1914
# HPV+ group
hpv_positive <- final_clinical_data %>% filter(hpv_status == "positive")

# HPV- group
hpv_negative <- final_clinical_data %>% filter(hpv_status == "negative")

# Source notebook line 1925
# Create the boxplot for the HPV+ group
ggplot(hpv_positive, aes(x = gender, y = Summation)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  
  labs(title = "Box Plot of Summation by Gender in HPV positive group",
       x = "Gender",
       y = "Summation") +
  theme_minimal()

# Check for normality (Shapiro-Wilk test) for each gender group separately
shapiro_test_male <- shapiro.test(hpv_positive %>% filter(gender == "male") %>% pull(Summation))
shapiro_test_female <- shapiro.test(hpv_positive %>% filter(gender == "female") %>% pull(Summation))

# Use Wilcoxon rank-sum test if either group is not normally distributed
if (shapiro_test_male$p.value < 0.05 | shapiro_test_female$p.value < 0.05) {
  # Data is not normally distributed, use Wilcoxon rank-sum test
  wilcox_test <- wilcox.test(Summation ~ gender, data = hpv_positive)
  test_result <- wilcox_test
  test_name <- "Wilcoxon rank-sum test"
} else {
  # Data is normally distributed, use t-test
  t_test <- t.test(Summation ~ gender, data = hpv_positive)
  test_result <- t_test
  test_name <- "t-test"
}

# Print the results of the test
cat("Result of the", test_name, ":\n")
print(test_result)

# Source notebook line 1961
# Create the boxplot for the HPV+ group
ggplot(hpv_negative, aes(x = gender, y = Summation)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  
  labs(title = "Box Plot of Summation by Gender in HPV negative group",
       x = "Gender",
       y = "Summation") +
  theme_minimal()

# Check for normality (Shapiro-Wilk test) for each gender group separately
shapiro_test_male <- shapiro.test(hpv_negative %>% filter(gender == "male") %>% pull(Summation))
shapiro_test_female <- shapiro.test(hpv_negative %>% filter(gender == "female") %>% pull(Summation))

# Use Wilcoxon rank-sum test if either group is not normally distributed
if (shapiro_test_male$p.value < 0.05 | shapiro_test_female$p.value < 0.05) {
  # Data is not normally distributed, use Wilcoxon rank-sum test
  wilcox_test <- wilcox.test(Summation ~ gender, data = hpv_negative)
  test_result <- wilcox_test
  test_name <- "Wilcoxon rank-sum test"
} else {
  # Data is normally distributed, use t-test
  t_test <- t.test(Summation ~ gender, data = hpv_negative)
  test_result <- t_test
  test_name <- "t-test"
}

# Print the results of the test
cat("Result of the", test_name, ":\n")
print(test_result)

# Source notebook line 1997
# Create the boxplot for the HPV+ group
ggplot(hpv_positive, aes(x = gender, y = Summation_cleaned)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  
  labs(title = "Box Plot of Summation_cleaned by Gender in HPV positive group",
       x = "Gender",
       y = "Summation_cleaned") +
  theme_minimal()

# Check for normality (Shapiro-Wilk test) for each gender group separately
shapiro_test_male <- shapiro.test(hpv_positive %>% filter(gender == "male") %>% pull(Summation_cleaned))
shapiro_test_female <- shapiro.test(hpv_positive %>% filter(gender == "female") %>% pull(Summation_cleaned))

# Use Wilcoxon rank-sum test if either group is not normally distributed
if (shapiro_test_male$p.value < 0.05 | shapiro_test_female$p.value < 0.05) {
  # Data is not normally distributed, use Wilcoxon rank-sum test
  wilcox_test <- wilcox.test(Summation_cleaned ~ gender, data = hpv_positive)
  test_result <- wilcox_test
  test_name <- "Wilcoxon rank-sum test"
} else {
  # Data is normally distributed, use t-test
  t_test <- t.test(Summation_cleaned ~ gender, data = hpv_positive)
  test_result <- t_test
  test_name <- "t-test"
}

# Print the results of the test
cat("Result of the", test_name, ":\n")
print(test_result)

# Source notebook line 2031
# Create the boxplot for the HPV+ group
ggplot(hpv_negative, aes(x = gender, y = Summation_cleaned)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  
  labs(title = "Box Plot of Summation_cleaned by Gender in HPV negative group",
       x = "Gender",
       y = "Summation_cleaned") +
  theme_minimal()

# Check for normality (Shapiro-Wilk test) for each gender group separately
shapiro_test_male <- shapiro.test(hpv_negative %>% filter(gender == "male") %>% pull(Summation_cleaned))
shapiro_test_female <- shapiro.test(hpv_negative %>% filter(gender == "female") %>% pull(Summation_cleaned))

# Use Wilcoxon rank-sum test if either group is not normally distributed
if (shapiro_test_male$p.value < 0.05 | shapiro_test_female$p.value < 0.05) {
  # Data is not normally distributed, use Wilcoxon rank-sum test
  wilcox_test <- wilcox.test(Summation_cleaned ~ gender, data = hpv_negative)
  test_result <- wilcox_test
  test_name <- "Wilcoxon rank-sum test"
} else {
  # Data is normally distributed, use t-test
  t_test <- t.test(Summation_cleaned ~ gender, data = hpv_negative)
  test_result <- t_test
  test_name <- "t-test"
}

# Print the results of the test
cat("Result of the", test_name, ":\n")
print(test_result)

