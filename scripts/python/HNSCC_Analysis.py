import pandas as pd
import os
import logging
import pyarrow  # Ensure pyarrow is imported
import fastparquet  # Ensure fastparquet is imported

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Define the file names (use absolute paths)
txt_file_name = 'methylation_level3_within_bioassay.txt'
parquet_file_name = 'methylation_level3_within_bioassay.parquet'

# Function to process and save chunks individually
def process_and_save_chunk(chunk, parquet_file, append=False):
    # Convert all columns to string type to avoid type conversion issues
    chunk = chunk.astype(str)
    
    # Write to Parquet file
    if append:
        chunk.to_parquet(parquet_file, engine='fastparquet', append=True)
    else:
        chunk.to_parquet(parquet_file, engine='fastparquet')

# Function to load the data in chunks and save as Parquet
def load_and_save_as_parquet(txt_file, parquet_file, chunk_size=100000):
    logging.info("Starting to load the dataset in chunks")

    try:
        first_chunk = True
        for chunk in pd.read_csv(txt_file, delimiter='\t', chunksize=chunk_size, low_memory=False):
            logging.info(f"Processed chunk with {len(chunk)} rows")
            process_and_save_chunk(chunk, parquet_file, append=not first_chunk)
            first_chunk = False
        
        logging.info("Dataset loaded and saved successfully")
    except Exception as e:
        logging.error(f"Error loading dataset: {e}")
        raise

# Check if the Parquet file exists
if os.path.exists(parquet_file_name):
    logging.info("Loading dataset from Parquet file")
    df = pd.read_parquet(parquet_file_name, engine='fastparquet')
else:
    logging.info("Parquet file not found. Loading dataset from text file and saving to Parquet.")
    load_and_save_as_parquet(txt_file_name, parquet_file_name)
    df = pd.read_parquet(parquet_file_name, engine='fastparquet')

# Get a list of all column names
logging.info("Column names:")
print(df.columns)

# Get a concise summary of the DataFrame
logging.info("DataFrame info:")
print(df.info())

# Get descriptive statistics for numerical columns
logging.info("Descriptive statistics:")
#print(df.describe())
