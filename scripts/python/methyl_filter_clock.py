import pandas as pd
import pyarrow.parquet as pq
from pyarrow import Table
from tqdm import tqdm

# Load the methylation dataset
input_file = "methyl_calculated_removed_indeter_hpv_and_missing_beta.parquet"
df = pd.read_parquet(input_file)

# Step 1: Remove the first row
df = df.iloc[1:].reset_index(drop=True)

# Step 2: Remove columns 2-5
df.drop(df.columns[1:5], axis=1, inplace=True)

# Function to check and print CSV columns
def check_csv_columns(file_path):
    df_csv = pd.read_csv(file_path)
    print(f"Columns in {file_path}: {df_csv.columns.tolist()}")
    return df_csv

# Check the columns in the CSV files
horvath_cpgs_df = check_csv_columns("horvath_cpgs.csv")
hannum_cpgs_df = check_csv_columns("hannum_cpgs.csv")
pheno_cpgs_df = check_csv_columns("pheno_cpgs.csv")

# Extract the CpG lists
horvath_cpgs = horvath_cpgs_df.iloc[:, 0].tolist()
hannum_cpgs = hannum_cpgs_df.iloc[:, 0].tolist()
pheno_cpgs = pheno_cpgs_df.iloc[:, 0].tolist()

# Function to filter and save dataset
def filter_and_save(df, cpg_list, cpg_name, output_file):
    missing_cpgs = [cpg for cpg in cpg_list if cpg not in df['Hybridization REF'].values]
    filtered_df = df[df['Hybridization REF'].isin(cpg_list)]
    
    if missing_cpgs:
        print(f"Missing CpGs in {cpg_name}: {missing_cpgs}")
    
    filtered_data_arrow = Table.from_pandas(filtered_df)
    pq.write_table(filtered_data_arrow, output_file)
    print(f"Filtered dataset for {cpg_name} saved to {output_file}")

# Filter and save the datasets
filter_and_save(df, horvath_cpgs, "Horvath", "filtered_methylation_data_horvath.parquet")
filter_and_save(df, hannum_cpgs, "Hannum", "filtered_methylation_data_hannum.parquet")
filter_and_save(df, pheno_cpgs, "PhenoAge", "filtered_methylation_data_phenoage.parquet")
