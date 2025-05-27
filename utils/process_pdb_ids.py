#!/usr/bin/env python3

import pandas as pd
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def process_pdb_ids():
    # Read the CSV file
    logging.info("Reading all_data_for_kfold.csv...")
    df = pd.read_csv('data/all_data_for_kfold.csv')
    
    # Print column names
    logging.info("Column names in the CSV file:")
    logging.info(df.columns.tolist())
    
    # Extract PDB IDs from the pdb_id column
    logging.info("Processing PDB IDs...")
    pdb_ids = df['pdb_id'].apply(lambda x: x.split('_')[0]).unique()
    
    # Sort the PDB IDs
    pdb_ids = sorted(pdb_ids)
    
    # Write to file
    output_file = 'official_training_pdbs.txt'
    logging.info(f"Writing {len(pdb_ids)} unique PDB IDs to {output_file}...")
    with open(output_file, 'w') as f:
        for pdb_id in pdb_ids:
            f.write(f"{pdb_id}\n")
    
    logging.info(f"Successfully wrote {len(pdb_ids)} unique PDB IDs to {output_file}")

if __name__ == "__main__":
    process_pdb_ids() 