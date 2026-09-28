import pandas as pd

# Define the path to your Parquet file
input_file = "methyl_calculated_removed_indeter_hpv_and_missing_beta.parquet"
output_file = "subdataset_for_10_perc.parquet"

# Read the Parquet file into a DataFrame
print("Reading Parquet file...")
df = pd.read_parquet(input_file)

# Extract the sub-dataset from the 2nd row and starting from the 6th column
sub_df = df.iloc[1:, 5:]

# Save the sub-dataset to a new Parquet file
print("Saving the sub-dataset to a new Parquet file...")
sub_df.to_parquet(output_file, index=False)

print("All operations completed and new Parquet file saved.")
