import pandas as pd
import numpy as np

# Define the path to the Parquet file
input_file = "subdataset_for_15_perc.parquet"

# Read the Parquet file into a DataFrame
print("Reading Parquet file...")
df = pd.read_parquet(input_file)

# Convert all values to numeric, coerce errors to NaN
print("Converting values to numeric...")
df = df.apply(pd.to_numeric, errors='coerce')

# Flatten the DataFrame to a single array, ignoring NaN values
all_values = df.values.flatten()
all_values = all_values[~np.isnan(all_values)]

# Calculate the 15th percentile value
percentile_15th = np.percentile(all_values, 15)

print(f"The 15th percentile value is: {percentile_15th}")
