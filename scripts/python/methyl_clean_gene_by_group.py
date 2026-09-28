import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import numpy as np

# Load the Parquet file
input_file = "methyl_summed_genes.parquet"
df = pd.read_parquet(input_file)

# Load the patient lists
patient_lists = pd.read_csv("patient_lists.csv")
hpv_pos_patient = patient_lists['HPV_Positive'].dropna().tolist()
hpv_neg_patient = patient_lists['HPV_Negative'].dropna().tolist()

# Calculate the 10th percentile among all beta-values (excluding the first column)
beta_values = df.iloc[:, 1:].values.flatten()
percentile_10th_value = np.percentile(beta_values, 10)
print(f"10th percentile value: {percentile_10th_value}")

# Filter out genes where >= 84 positive patients' beta-values are less than the 10th percentile value
# AND >= 363 negative patients' beta-values are less than the 10th percentile value
def filter_genes_by_group(row, pos_patients, neg_patients, threshold, pos_min, neg_min):
    pos_values = row[pos_patients].astype(float)
    neg_values = row[neg_patients].astype(float)
    pos_condition = (pos_values < threshold).sum() >= pos_min
    neg_condition = (neg_values < threshold).sum() >= neg_min
    return not (pos_condition and neg_condition)

filtered_df = df[df.apply(lambda row: filter_genes_by_group(row, hpv_pos_patient, hpv_neg_patient, percentile_10th_value, 84, 363), axis=1)]

# Save the resulting dataset to a new Parquet file
output_file = "methyl_cleaned_bygroup_genes.parquet"
table = pa.Table.from_pandas(filtered_df)
pq.write_table(table, output_file)

print(f"Filtered data saved to {output_file}")
