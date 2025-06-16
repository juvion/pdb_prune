#!/usr/bin/env python3

import pandas as pd
import argparse

def match_csv_files(file1: str, file2: str, output_file: str = "merged_results.csv"):
    """
    Merge two CSV files with specific matching logic.
    
    Args:
        file1 (str): Path to first CSV file (test_pdb_meta.csv)
        file2 (str): Path to second CSV file (test_data_meta.csv)
        output_file (str): Path to output CSV file
    """
    # Read CSV files
    df_pdb = pd.read_csv(file1)  # test_pdb_meta.csv
    df_data = pd.read_csv(file2)  # test_data_meta.csv
    
    # Add uppercased keys for case-insensitive matching
    df_pdb['PDB_id_uc'] = df_pdb['PDB_id'].str.upper()
    df_data['PDB_id_uc'] = df_data['PDB_id'].str.upper()
    
    # Merge on PDB_id (uppercase version)
    merged_df = pd.merge(
        df_pdb,
        df_data,
        left_on=['PDB_id_uc'],
        right_on=['PDB_id_uc'],
        how='outer',
        suffixes=('_pdb', '_data')
    )
    
    # Initialize match_type column
    merged_df['match_type'] = 'no_match'
    
    # Apply matching logic
    # Perfect match: both chain_id and length match
    perfect_match = (
        (merged_df['chain_id_pdb'] == merged_df['chain_id_data']) & 
        (merged_df['length_pdb'] == merged_df['length_data'])
    )
    
    # Chain mismatch: chain_id matches but length doesn't
    chain_mismatch = (
        (merged_df['chain_id_pdb'] == merged_df['chain_id_data']) & 
        (merged_df['length_pdb'] != merged_df['length_data'])
    )
    
    # Length mismatch: chain_id doesn't match but length does
    length_mismatch = (
        (merged_df['chain_id_pdb'] != merged_df['chain_id_data']) & 
        (merged_df['length_pdb'] == merged_df['length_data'])
    )
    
    # PDB mismatch: both chain_id and length don't match
    pdb_mismatch = (
        (merged_df['chain_id_pdb'] != merged_df['chain_id_data']) & 
        (merged_df['length_pdb'] != merged_df['length_data'])
    )
    
    # Apply match types
    merged_df.loc[perfect_match, 'match_type'] = 'perfect'
    merged_df.loc[chain_mismatch, 'match_type'] = 'chain'
    merged_df.loc[length_mismatch, 'match_type'] = 'length'
    merged_df.loc[pdb_mismatch, 'match_type'] = 'pdb'
    
    # Create final result DataFrame with all columns
    result_df = pd.DataFrame({
        'PDB_id': merged_df['PDB_id_pdb'],
        'chain_id_pdb': merged_df['chain_id_pdb'],
        'length_pdb': merged_df['length_pdb'],
        'chain_id_data': merged_df['chain_id_data'],
        'length_data': merged_df['length_data'],
        'match_type': merged_df['match_type']
    })
    
    # Save results
    result_df.to_csv(output_file, index=False)
    
    # Print summary
    print("\nMatching Summary:")
    print(f"Total entries: {len(result_df)}")
    print("\nMatch types:")
    print(result_df['match_type'].value_counts())
    print(f"\nResults saved to {output_file}")
    
    # Print sample of each match type for verification
    print("\nSample of each match type:")
    for match_type in ['perfect', 'chain', 'length', 'pdb']:
        print(f"\n{match_type.upper()} matches:")
        print(result_df[result_df['match_type'] == match_type].head())

def main():
    parser = argparse.ArgumentParser(description='Merge and compare two CSV files with specific matching logic')
    parser.add_argument('--file1', type=str, required=True,
                      help='Path to first CSV file (test_pdb_meta.csv)')
    parser.add_argument('--file2', type=str, required=True,
                      help='Path to second CSV file (test_data_meta.csv)')
    parser.add_argument('--output', type=str, default='merged_results.csv',
                      help='Output CSV file path (default: merged_results.csv)')
    
    args = parser.parse_args()
    
    match_csv_files(file1=args.file1, file2=args.file2, output_file=args.output)

if __name__ == "__main__":
    main() 