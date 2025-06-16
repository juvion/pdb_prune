#!/usr/bin/env python3

import pandas as pd
import os
import shutil
import argparse
from pathlib import Path

def copy_files(csv_file: str, pdb_folder: str, fasta_folder: str, output_folder: str):
    """
    Copy PDB and FASTA files based on PDB_id and chain_id from CSV file.
    
    Args:
        csv_file (str): Path to CSV file containing PDB_id and chain_id_pdb
        pdb_folder (str): Path to folder containing PDB files
        fasta_folder (str): Path to folder containing FASTA files
        output_folder (str): Path to output folder for copied files
    """
    # Create output folders if they don't exist
    pdb_output = os.path.join(output_folder, 'pdb')
    fasta_output = os.path.join(output_folder, 'fasta')
    os.makedirs(pdb_output, exist_ok=True)
    os.makedirs(fasta_output, exist_ok=True)
    
    # Read CSV file
    df = pd.read_csv(csv_file)
    
    # Initialize counters and lists for missing files
    total_files = len(df)
    copied_pdb = 0
    copied_fasta = 0
    missing_pdb = []
    missing_fasta = []
    
    print(f"\nProcessing {total_files} entries from {csv_file}")
    print(f"Looking for PDB files in: {pdb_folder}")
    print(f"Looking for FASTA files in: {fasta_folder}")
    
    # Process each row
    for index, row in df.iterrows():
        pdb_id = row['PDB_id'].lower()  # Convert to lowercase for filename
        chain_id = row['chain_id_pdb']
        
        # Construct filenames
        pdb_filename = f"{pdb_id}_{chain_id}.pdb"
        fasta_filename = f"{pdb_id}_{chain_id}.fasta"
        
        # Source paths
        pdb_source = os.path.join(pdb_folder, pdb_filename)
        fasta_source = os.path.join(fasta_folder, fasta_filename)
        
        # Destination paths
        pdb_dest = os.path.join(pdb_output, pdb_filename)
        fasta_dest = os.path.join(fasta_output, fasta_filename)
        
        # Copy PDB file if it exists
        if os.path.exists(pdb_source):
            shutil.copy2(pdb_source, pdb_dest)
            copied_pdb += 1
        else:
            missing_pdb.append(pdb_filename)
        
        # Copy FASTA file if it exists
        if os.path.exists(fasta_source):
            shutil.copy2(fasta_source, fasta_dest)
            copied_fasta += 1
        else:
            missing_fasta.append(fasta_filename)
        
        # Print progress every 10 files
        if (index + 1) % 10 == 0:
            print(f"Processed {index + 1}/{total_files} entries...")
    
    # Print summary
    print("\nCopy Summary:")
    print(f"Total entries processed: {total_files}")
    print(f"PDB files copied: {copied_pdb}")
    print(f"FASTA files copied: {copied_fasta}")
    
    if missing_pdb:
        print(f"\nMissing PDB files ({len(missing_pdb)}):")
        for file in missing_pdb[:10]:  # Show first 10 missing files
            print(f"  - {file}")
        if len(missing_pdb) > 10:
            print(f"  ... and {len(missing_pdb) - 10} more")
    
    if missing_fasta:
        print(f"\nMissing FASTA files ({len(missing_fasta)}):")
        for file in missing_fasta[:10]:  # Show first 10 missing files
            print(f"  - {file}")
        if len(missing_fasta) > 10:
            print(f"  ... and {len(missing_fasta) - 10} more")
    
    print(f"\nFiles have been copied to:")
    print(f"  PDB files: {pdb_output}")
    print(f"  FASTA files: {fasta_output}")

def main():
    parser = argparse.ArgumentParser(description='Copy PDB and FASTA files based on CSV list')
    parser.add_argument('--csv', type=str, required=True,
                      help='Path to CSV file containing PDB_id and chain_id_pdb')
    parser.add_argument('--pdb_folder', type=str, required=True,
                      help='Path to folder containing PDB files')
    parser.add_argument('--fasta_folder', type=str, required=True,
                      help='Path to folder containing FASTA files')
    parser.add_argument('--output_folder', type=str, required=True,
                      help='Path to output folder for copied files')
    
    args = parser.parse_args()
    
    copy_files(
        csv_file=args.csv,
        pdb_folder=args.pdb_folder,
        fasta_folder=args.fasta_folder,
        output_folder=args.output_folder
    )

if __name__ == "__main__":
    main() 