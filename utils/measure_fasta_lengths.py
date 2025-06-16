#!/usr/bin/env python3

import os
from Bio import SeqIO
import pandas as pd
from pathlib import Path
import argparse

def measure_fasta_lengths(fasta_dir: str, output_csv: str = "fasta_lengths.csv"):
    """
    Measure the length of each FASTA file and export to CSV.
    
    Args:
        fasta_dir (str): Directory containing FASTA files
        output_csv (str): Output CSV file path
    """
    # Get all FASTA files
    fasta_files = list(Path(fasta_dir).glob("*.fasta"))
    if not fasta_files:
        print(f"No FASTA files found in {fasta_dir}")
        return
    
    # Collect file lengths
    file_lengths = []
    
    for fasta_file in fasta_files:
        try:
            # Read the first sequence from each file
            with open(fasta_file, 'r') as f:
                first_record = next(SeqIO.parse(f, "fasta"))
                file_lengths.append({
                    'file_name': fasta_file.name,
                    'length': len(first_record.seq)
                })
        except Exception as e:
            print(f"Error processing {fasta_file}: {str(e)}")
    
    if not file_lengths:
        print("No valid FASTA files found")
        return
    
    # Create DataFrame and save to CSV
    df = pd.DataFrame(file_lengths)
    df.to_csv(output_csv, index=False)
    print(f"\nFile lengths saved to {output_csv}")
    
    # Print summary
    print(f"\nProcessed {len(file_lengths)} FASTA files")
    print(f"Min length: {df['length'].min()}")
    print(f"Max length: {df['length'].max()}")
    print(f"Mean length: {df['length'].mean():.1f}")

def main():
    parser = argparse.ArgumentParser(description='Measure FASTA file lengths and export to CSV')
    parser.add_argument('--fasta-dir', type=str, required=True,
                      help='Directory containing input FASTA files')
    parser.add_argument('--output', type=str, default='fasta_lengths.csv',
                      help='Output CSV file path (default: fasta_lengths.csv)')
    
    args = parser.parse_args()
    
    measure_fasta_lengths(fasta_dir=args.fasta_dir, output_csv=args.output)

if __name__ == "__main__":
    main() 