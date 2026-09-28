import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from tqdm import tqdm

# Define the path to the Parquet file
input_file = "cleaned_methylation_data_merge.parquet"
output_file = "cleaned_methylation_data_final.parquet"

# Read the Parquet file into a DataFrame
df = pd.read_parquet(input_file)

# Replace "nan" with empty strings in the patient columns
patient_columns = df.columns[3:]
df[patient_columns] = df[patient_columns].apply(lambda col: col.map(lambda x: '' if pd.isna(x) else x))

# Initialize Low and High columns with NaN values
low_high_df = pd.DataFrame({'Low': [''] * len(df), 'High': [''] * len(df)})
df = pd.concat([df.iloc[:, :3], low_high_df, df.iloc[:, 3:]], axis=1)

# Calculate Low and High values for each row, using tqdm for progress tracking
for index, row in tqdm(df.iterrows(), total=df.shape[0], desc="Processing rows"):
    values = pd.to_numeric(row[5:], errors='coerce').dropna()
    if not values.empty:
        df.at[index, 'Low'] = values.min()
        df.at[index, 'High'] = values.max()

# Ensure 'Low' and 'High' columns are strings
df['Low'] = df['Low'].apply(lambda x: '' if pd.isna(x) else str(x))
df['High'] = df['High'].apply(lambda x: '' if pd.isna(x) else str(x))

# Add the 'Sum' row as the new first row
sum_row = [''] * 5 + df.iloc[:, 5:].apply(pd.to_numeric, errors='coerce').sum().tolist()
sum_df = pd.DataFrame([sum_row], columns=df.columns)
sum_df.index = ['Sum']
df = pd.concat([sum_df, df]).reset_index(drop=True)

# Ensure all patient columns are numeric before saving to Parquet
df[patient_columns] = df[patient_columns].apply(pd.to_numeric, errors='coerce')

# Convert all columns to strings to avoid type conversion issues
df = df.astype(str)

# Save the modified DataFrame to a new Parquet file
df.to_parquet(output_file)

print("All operations completed and new Parquet file saved.")
