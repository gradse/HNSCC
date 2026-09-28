import pandas as pd
from tqdm import tqdm

# Define the path to your Parquet file
input_file = "methyl_calculated_removed_indeter_hpv.parquet"
output_file = "methyl_calculated_removed_indeter_hpv_and_missing_beta.parquet"

# Read the Parquet file into a DataFrame
print("Reading Parquet file...")
df = pd.read_parquet(input_file)

# Ensure Low column is treated as a string to detect empty strings correctly
df['Low'] = df['Low'].astype(str)

# Remove rows where the Low column is empty or NaN
print("Removing rows with empty Low values...")
df_cleaned = df[df['Low'].notna() & (df['Low'] != '')]

# Save the cleaned DataFrame to a new Parquet file
print("Saving the cleaned DataFrame to a new Parquet file...")
df_cleaned.to_parquet(output_file, index=False)

print("All operations completed and new Parquet file saved.")
