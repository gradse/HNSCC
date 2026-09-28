import pandas as pd
import pyarrow.parquet as pq

# Define the path to your Parquet file
parquet_file = "cleaned_methylation_data.parquet"
output_file = "cleaned_methylation_data_no_first_row.parquet"

# Read the Parquet file into a DataFrame
df = pd.read_parquet(parquet_file)

# Remove the first row
df_no_first_row = df.iloc[1:].reset_index(drop=True)

# Save the modified DataFrame to a new Parquet file
df_no_first_row.to_parquet(output_file)

print("First row removed and new Parquet file saved.")
