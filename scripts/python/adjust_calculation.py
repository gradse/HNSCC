import pandas as pd
import pyarrow.parquet as pq
from tqdm import tqdm

# Define the path to the Parquet file
input_file = "cleaned_methylation_data_final.parquet"

# Read the Parquet file into a DataFrame
print("Reading Parquet file...")
methyl_merge_calculated = pq.read_table(input_file).to_pandas()

# Rename the first entry in the 'Hybridization REF' column to "Summation"
print("Renaming the first entry in 'Hybridization REF' column to 'Summation'...")
methyl_merge_calculated.at[0, 'Hybridization REF'] = 'Summation'

# Extract the first row, starting from the 6th column
first_row_data = methyl_merge_calculated.iloc[0, 5:]

# Convert the first row data to numeric for calculation
print("Converting first row data to numeric for calculation...")
numeric_first_row_data = pd.to_numeric(first_row_data, errors='coerce')

# Find the minimum and maximum values in the first row
print("Finding the minimum and maximum values in the first row...")
min_value = numeric_first_row_data.min()
max_value = numeric_first_row_data.max()

# Store the minimum and maximum values in the 'Low' and 'High' columns
print("Storing the minimum and maximum values in 'Low' and 'High' columns...")
methyl_merge_calculated.at[0, 'Low'] = str(min_value)
methyl_merge_calculated.at[0, 'High'] = str(max_value)

# Save the modified DataFrame back to Parquet file
output_file = "cleaned_methylation_data_final.parquet"
print("Saving the modified DataFrame back to Parquet file...")
methyl_merge_calculated.to_parquet(output_file, index=False)

print("All operations completed and new Parquet file saved.")
