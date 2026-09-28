import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import numpy as np

# Load the Parquet file
input_file = "methyl_summed_genes.parquet"
df = pd.read_parquet(input_file)

# Calculate the 10th percentile among all beta-values (excluding the first column)
beta_values = df.iloc[:, 1:].values.flatten()
percentile_10th_value = np.percentile(beta_values, 10)
print(f"10th percentile value: {percentile_10th_value}")

# Filter out genes where >= 473 patients' beta-values are less than the 10th percentile value
def filter_genes(row, threshold, min_patients):
    return (row[1:] < threshold).sum() < min_patients

filtered_df = df[df.apply(filter_genes, axis=1, args=(percentile_10th_value, 473))]

# Save the resulting dataset to a new Parquet file
output_file = "methyl_cleaned_genes.parquet"
table = pa.Table.from_pandas(filtered_df)
pq.write_table(table, output_file)

print(f"Filtered data saved to {output_file}")
