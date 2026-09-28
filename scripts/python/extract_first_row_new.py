import pandas as pd

# Define the input and output file paths
input_file = "methyl_removed_genes.parquet"
output_file = "extracted_first_row_new.csv"

# Read the Parquet file into a DataFrame
print("Reading Parquet file...")
df = pd.read_parquet(input_file)

# Extract the first row, starting from the 7th column
first_row = df.iloc[0, 6:]

# Convert the first row to a DataFrame to retain the header
first_row_df = first_row.to_frame().T

# Save the extracted row to a CSV file, keeping all decimals
first_row_df.to_csv(output_file, index=False)

print("First row extracted and saved to CSV file.")
