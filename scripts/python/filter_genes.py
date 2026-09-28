import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

# Load the Parquet file
input_file = "methyl_calculated_removed_indeter_hpv_and_missing_beta.parquet"
df = pd.read_parquet(input_file)

# Remove the first row (not headers)
df = df.iloc[1:, :]

# Remove the 1st, 3rd, 4th, and 5th columns
df = df.drop(columns=[df.columns[0], df.columns[2], df.columns[3], df.columns[4]])

# Rename columns for clarity
df.columns = ['Gene_Symbol'] + df.columns[1:].tolist()

# Remove rows where 'Gene_Symbol' is 'nan' or empty
df = df[df['Gene_Symbol'].notna() & df['Gene_Symbol'].str.strip().astype(bool) & (df['Gene_Symbol'] != 'nan')]

# Convert beta-values to numeric, keeping all decimals
for col in df.columns[1:]:
    df[col] = pd.to_numeric(df[col], errors='coerce')

# Sum beta-values for each gene across multiple rows
df_grouped = df.groupby('Gene_Symbol').sum()

# Reset the index to have 'Gene_Symbol' as a column again
df_grouped.reset_index(inplace=True)

# Save the resulting dataset to a new Parquet file
output_file = "methyl_summed_genes.parquet"
table = pa.Table.from_pandas(df_grouped)
pq.write_table(table, output_file)

print(f"Processed data saved to {output_file}")