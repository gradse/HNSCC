import pandas as pd
import numpy as np
from sklearn.impute import SimpleImputer
import pyarrow as pa
import pyarrow.parquet as pq

# Define the input and output file paths
input_file = "methyl_clock.parquet"
output_matrix_file = "methyl_imputed_matrix.parquet"
output_row_names_file = "methyl_row_names.parquet"

# Read the Parquet file into a DataFrame
print("Reading Parquet file...")
df = pd.read_parquet(input_file)

# Extract CpG identifiers
cpg_ids = df['Hybridization REF'].values

# Convert relevant columns to numeric, coercing errors to NaN
beta_values = df.iloc[:, 1:].apply(pd.to_numeric, errors='coerce')

# Impute missing values using mean imputation
print("Imputing missing values using mean imputation...")
imputer = SimpleImputer(strategy='mean')
beta_values_imputed = imputer.fit_transform(beta_values)

# Create DataFrame from imputed matrix
imputed_df = pd.DataFrame(beta_values_imputed, columns=df.columns[1:])

# Reattach the first column (CpG identifiers)
imputed_df.insert(0, 'Hybridization REF', cpg_ids)

# Save the imputed DataFrame to a new Parquet file
print("Saving the imputed DataFrame to a new Parquet file...")
imputed_df.to_parquet(output_matrix_file, index=False)

# Save the row names separately to Parquet
print("Saving the row names to a new Parquet file...")
row_names_df = pd.DataFrame({'Hybridization REF': cpg_ids})
row_names_df.to_parquet(output_row_names_file, index=False)

print("All operations completed and new Parquet files saved.")
