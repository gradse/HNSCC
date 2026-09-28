import pandas as pd
import numpy as np

# Define the input and output file paths
input_file = "methyl_calculated_removed_indeter_hpv_and_missing_beta.parquet"
output_file = "methy_calculated_patient.parquet"

# Define the 10th percentile value
percentile_15th_value = 0.0347440141330495

# Read the Parquet file into a DataFrame
print("Reading Parquet file...")
df = pd.read_parquet(input_file)

# Initialize the new column with NaN values
df.insert(5, "Patients_Beta_Value_Count", np.nan)

# Find the starting column index where the header is "tcga-4p-aa8j"
start_col_index = df.columns.get_loc("tcga-4p-aa8j")

# Calculate the count of beta_values <= percentile_15th_value for each row (starting from the second row)
print("Calculating the number of patients' beta_value <= 0.0347440141330495...")
for index, row in df.iterrows():
    if index > 0:  # Skip the first row
        beta_values = row[start_col_index:]
        # Convert beta_values to numeric for comparison, keep as strings in the DataFrame
        count = (pd.to_numeric(beta_values, errors='coerce') <= percentile_15th_value).sum()
        df.at[index, "Patients_Beta_Value_Count"] = count

# Save the updated DataFrame to a new Parquet file
print("Saving the updated DataFrame to a new Parquet file...")
df.to_parquet(output_file, index=False)
