import pandas as pd
import pyarrow.parquet as pq
from pyarrow import Table

# Load the Parquet file
input_file = "methyl_calculated_removed_indeter_hpv.parquet"
df = pd.read_parquet(input_file)

# Remove the first row and columns 2 to 5
df = df.iloc[1:]  # Remove the first row
df = df.drop(df.columns[1:5], axis=1)  # Remove columns 2 to 5

# Convert beta-values to numeric, keeping all decimals, and impute "nan" with the minimum value of the respective column
for col in df.columns[1:]:  # Skip the first column as it contains CpG names
    df[col] = pd.to_numeric(df[col], errors='coerce')
    median_value = df[col].min(skipna=True)
    df[col] = df[col].fillna(median_value)

# Add a new row called "Summation" as the first row
summation_row = df.iloc[:, 1:].sum()
summation_df = pd.DataFrame([summation_row], columns=df.columns[1:])
summation_df.insert(0, 'Hybridization REF', 'Summation')
df = pd.concat([summation_df, df]).reset_index(drop=True)

# Ensure no scientific notation
pd.set_option('display.float_format', lambda x: f'{x:.16f}')

# Save the entire methylation data into a Parquet file
output_file = "methyl_imputed_sum.parquet"
df_arrow = Table.from_pandas(df)
pq.write_table(df_arrow, output_file)
print(f"Entire methylation data saved to {output_file}")

# Extract the summation row into a CSV file without the "Hybridization REF" column
summation_row_with_header = summation_row.to_frame().T
summation_row_with_header.index = ['Summation']

# Save the CSV file with full decimal precision
summation_row_with_header.to_csv("new_summation_row.csv", index=True, header=True, float_format='%.16f')
print("Summation row saved to new_summation_row.csv")

# Function to filter and save dataset
def filter_and_save(df, cpg_list, cpg_name, output_file):
    missing_cpgs = [cpg for cpg in cpg_list if cpg not in df['Hybridization REF'].values]
    filtered_df = df[df['Hybridization REF'].isin(cpg_list)]
    
    if missing_cpgs:
        print(f"Missing CpGs in {cpg_name}: {missing_cpgs}")
    
    filtered_data_arrow = Table.from_pandas(filtered_df)
    pq.write_table(filtered_data_arrow, output_file)
    print(f"Filtered dataset for {cpg_name} saved to {output_file}")

# Load the CpG lists from CSV files
horvath_cpgs_df = pd.read_csv("horvath_cpgs.csv")
hannum_cpgs_df = pd.read_csv("hannum_cpgs.csv")
pheno_cpgs_df = pd.read_csv("pheno_cpgs.csv")

# Extract the CpG lists
horvath_cpgs = horvath_cpgs_df.iloc[:, 0].tolist()
hannum_cpgs = hannum_cpgs_df.iloc[:, 0].tolist()
pheno_cpgs = pheno_cpgs_df.iloc[:, 0].tolist()

# Filter and save the datasets for each clock
filter_and_save(df, horvath_cpgs, "Horvath", "filtered_methyl_horvath_new.parquet")
filter_and_save(df, hannum_cpgs, "Hannum", "filtered_methyl_hannum_new.parquet")
filter_and_save(df, pheno_cpgs, "PhenoAge", "filtered_methyl_phenoage_new.parquet")
