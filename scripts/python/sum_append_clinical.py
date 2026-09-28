import pandas as pd

# Define the path to the Parquet file
input_file = "cleaned_methylation_data_final.parquet"
output_file = "first_row_extracted.csv"

# Read the Parquet file into a DataFrame
print("Reading Parquet file...")
methyl_merge_calculated = pd.read_parquet(input_file)

# Extract the first row, starting from the 6th column
extracted_data = methyl_merge_calculated.iloc[0, 5:]

# Convert the extracted data to a DataFrame to include headers in the CSV
extracted_data_df = pd.DataFrame(extracted_data).transpose()

# Save the extracted data to a CSV file
print("Saving the extracted data to a CSV file...")
extracted_data_df.to_csv(output_file, index=False)

print("All operations completed and CSV file saved.")
