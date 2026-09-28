import pandas as pd
import numpy as np

# Define the input and output file paths
input_file = "methy_calculated_patient.parquet"
output_file = "methyl_removed_genes.parquet"

# Read the Parquet file into a DataFrame
print("Reading Parquet file...")
df = pd.read_parquet(input_file)

# Replace string 'nan' with actual NaN values
df.replace('nan', np.nan, inplace=True)

# Preserve the first row (Summation)
summation_row = df.iloc[0].copy()

# Remove rows where Patients_Beta_Value_Count >= 447, excluding the first row
print("Removing rows where Patients_Beta_Value_Count >= 447...")
df = df.iloc[1:]
df = df[df["Patients_Beta_Value_Count"] < 447]

# Reinsert the Summation row at the top
df = pd.concat([summation_row.to_frame().T, df]).reset_index(drop=True)

# Recalculate the first row "Summation"
print("Recalculating the Summation row...")

# Specify the columns for summation (4th, 5th, and 7th onward are patients)
columns_to_summate = ['Low', 'High'] + list(df.columns[6:])

# Convert the specified columns to numeric for summation, ignoring NaN values
summation_values = df.loc[1:, columns_to_summate].apply(pd.to_numeric, errors='coerce').sum(axis=0, skipna=True)

# Convert the summation values back to strings with full decimals
summation_values_str = summation_values.apply(lambda x: f"{x:.16f}")

# Update the Summation row in the DataFrame
df.loc[0, columns_to_summate] = summation_values_str

# Find the Low and High values for the Summation row
low_value = summation_values.min()
high_value = summation_values.max()

# Update the Low and High columns in the Summation row
df.loc[0, 'Low'] = f"{low_value:.16f}"
df.loc[0, 'High'] = f"{high_value:.16f}"

# Save the updated DataFrame to a new Parquet file
print("Saving the updated DataFrame to a new Parquet file...")
df.to_parquet(output_file, index=False)

print("All operations completed and new Parquet file saved.")
