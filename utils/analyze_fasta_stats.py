#!/usr/bin/env python3

import os
from pathlib import Path
from Bio import SeqIO
import pandas as pd
import argparse

def analyze_fasta_files(fasta_dir: str) -> pd.DataFrame:
    """Analyze FASTA files and return a DataFrame with statistics."""
    fasta_dir = Path(fasta_dir)
    stats = []
    
    for fasta_file in fasta_dir.glob("*.fasta"):
        try:
            # Read the FASTA file
            record = next(SeqIO.parse(fasta_file, "fasta"))
            
            # Extract PDB code and chain ID from filename
            pdb_code = fasta_file.stem.split('_')[0]  # Remove 'pdb' prefix
            chain_id = fasta_file.stem.split('_')[1]
            
            # Get sequence and length
            sequence = str(record.seq)
            length = len(sequence)
            
            stats.append({
                'sequence_id': f"{pdb_code}_{chain_id}",
                'sequence': sequence,
                'length': length
            })
            
        except Exception as e:
            print(f"Error processing {fasta_file}: {e}")
            continue
    
    # Create DataFrame and sort by length
    df = pd.DataFrame(stats)
    df = df.sort_values('length', ascending=False)
    return df

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Analyze FASTA files and generate statistics.')
    parser.add_argument('--fasta-dir', type=str, required=True,
                      help='Directory containing FASTA files')
    args = parser.parse_args()
    
    # Analyze FASTA files
    df = analyze_fasta_files(args.fasta_dir)
    
    # Print summary statistics
    print("\nFASTA File Statistics Summary:")
    print("-" * 80)
    print(f"Total number of sequences: {len(df)}")
    print(f"Average sequence length: {df['length'].mean():.2f}")
    print(f"Minimum sequence length: {df['length'].min()}")
    print(f"Maximum sequence length: {df['length'].max()}")
    print("\nTop 10 longest sequences:")
    print("-" * 80)
    print(df[['sequence_id', 'length', 'sequence']].head(10).to_string(index=False))
    
    # Save to CSV
    output_file = "fasta_statistics.csv"
    df.to_csv(output_file, index=False)
    print(f"\nDetailed statistics saved to {output_file}")

if __name__ == "__main__":
    main() 