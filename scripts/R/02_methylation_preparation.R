# Analysis stage from archive/HNSCC_Analysis.Rmd (source lines 227-408).
# Run through scripts/run_r_analysis.R.

# Source notebook line 227
# Open the Parquet file
dataset <- open_dataset('methylation_level3_within_bioassay.parquet')

# Read the first few rows
first_few_rows <- as.data.frame(dataset %>% slice_head(n = 10))

# Identify the columns of interest
initial_columns <- colnames(first_few_rows)[c(1, 3, 4)]

# Identify Beta_value columns
beta_value_columns <- names(first_few_rows)[which(first_few_rows[1, ] == "Beta_value")]

# Combine the columns to keep
columns_to_keep <- c(initial_columns, beta_value_columns)

# Filter the dataset to keep only the relevant columns
filtered_dataset <- dataset %>%
  dplyr::select(all_of(columns_to_keep)) 

# Convert to data frame
cleaned_data <- as.data.frame(filtered_dataset)

# Rename the second to third columns to have blank headers
colnames(cleaned_data)[2:3] <- c("Gene_Symbol", "Chromosome")

# Save the cleaned dataframe to a new Parquet file
write_parquet(cleaned_data, "cleaned_methylation_data.parquet")

# Source notebook line 267
extract_id <- function(header) {
  str_sub(header, 9, 12)
}

# Extract column headers starting from the 4th column
headers <- colnames(cleaned_data)[4:ncol(cleaned_data)]

# Extract the specific substring from each header
extracted_ids <- sapply(headers, extract_id)

# Check for duplicates
duplicate_ids <- extracted_ids[duplicated(extracted_ids)]

# Count the number of unique and duplicate IDs
unique_ids <- unique(extracted_ids)
num_unique_ids <- length(unique_ids)
num_duplicate_ids <- length(duplicate_ids)

# Print the results
cat("Number of unique IDs:", num_unique_ids, "\n")
cat("Number of duplicate IDs:", num_duplicate_ids, "\n")

# Source notebook line 295
# Function to extract the 14th and 15th digits from the column headers
extract_digits_14_15 <- function(header) {
  str_sub(header, 14, 15)
}

# Extract column headers starting from the 4th column
headers <- colnames(cleaned_data)[4:ncol(cleaned_data)]

# Extract the 14th and 15th digits from each header
digits_14_15 <- sapply(headers, extract_digits_14_15)

# Filter out columns where the 14th and 15th digits are not "01"
valid_columns <- headers[digits_14_15 == "01"]
removed_columns <- headers[digits_14_15 != "01"]

cleaned_data <- cleaned_data %>% select(all_of(c(colnames(cleaned_data)[1:3], valid_columns)))

# Print the removed columns
cat("Removed columns (sample type is not '01'):\n")
print(removed_columns)

# Extract the specific substring from each header
headers <- colnames(cleaned_data)[4:ncol(cleaned_data)]
extracted_ids <- sapply(headers, extract_id)

# Find duplicate IDs
duplicate_ids <- extracted_ids[duplicated(extracted_ids)]

# Print the results
cat("Number of unique IDs:", length(unique(extracted_ids)), "\n")
cat("Number of duplicate IDs:", length(duplicate_ids), "\n")

# If there are duplicates, print them
if (length(duplicate_ids) > 0) {
  cat("Duplicate IDs:\n")
  print(duplicate_ids)
}

# Source notebook line 339
remaining_headers <- colnames(cleaned_data)[4:ncol(cleaned_data)]
remaining_digits_14_15 <- sapply(remaining_headers, extract_digits_14_15)

if (all(remaining_digits_14_15 == "01")) {
  cat("All remaining columns have '01' as the sample type.\n")
} else {
  cat("Some remaining columns do not have '01' as the sample type:\n")
  print(remaining_headers[remaining_digits_14_15 != "01"])
}

# Save the cleaned dataframe to a new parquet file
write_parquet(cleaned_data, "cleaned_methylation_data.parquet")

# Source notebook line 382
# transpose the clinical data
clinical_focused_merge <- t(clinical_focused)

clinical_focused_merge <- as.data.frame(clinical_focused_merge)

cleaned_data <- open_dataset("cleaned_methylation_data_no_first_row.parquet")

cleaned_data <- as.data.frame(cleaned_data)

# Rename headers for each column starting from the 4th
rename_header <- function(header) {
  tolower(str_sub(header, 1, 12))
}

new_headers <- sapply(colnames(cleaned_data)[4:ncol(cleaned_data)], rename_header)

colnames(cleaned_data)[4:ncol(cleaned_data)] <- new_headers

# Save the new dataframe to cleaned_data_merge
cleaned_data_merge <- cleaned_data

# Save the cleaned dataframe to a new parquet file
write_parquet(cleaned_data_merge, "cleaned_methylation_data_merge.parquet")

