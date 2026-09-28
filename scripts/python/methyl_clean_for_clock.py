import pandas as pd

# Define the input and output file paths
input_file = "methyl_removed_genes.parquet"
output_file = "methyl_clock.parquet"

# Read the Parquet file into a DataFrame
print("Reading Parquet file...")
df = pd.read_parquet(input_file, engine='pyarrow')

# Remove the first row
df = df.iloc[1:]

# Remove columns 2 to 6 
df = df.drop(df.columns[1:6], axis=1)

# Save the updated DataFrame to a new Parquet file
print("Saving the updated DataFrame to a new Parquet file...")
df.to_parquet(output_file, index=False, engine='pyarrow')

print("All operations completed and new Parquet file saved.")
