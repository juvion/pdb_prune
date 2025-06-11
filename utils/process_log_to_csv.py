#!/usr/bin/env python3

import csv
import re

def process_log_file(input_file: str, output_file: str):
    """Process log file and extract PDB codes, chain IDs, and lengths."""
    results = []
    
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip lines starting with "("
            if line.startswith('('):
                continue
                
            # Extract PDB code (first 4 characters)
            pdb_code = line[:4]
            
            # Extract chain ID (between second "_" and ".")
            chain_id = line.split('_')[2].split('.')[0]
            
            # Extract length (number after "=")
            length_match = re.search(r'len=(\d+)', line)
            length = int(length_match.group(1)) if length_match else None
            
            results.append({
                'PDB_id': pdb_code,
                'chain_id': chain_id,
                'length': length
            })
    
    # Write to CSV
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['PDB_id', 'chain_id', 'length'])
        writer.writeheader()
        writer.writerows(results)

def main():
    input_file = "data/0607_output_log.txt"
    output_file = "data/pdb_chain_lengths.csv"
    
    process_log_file(input_file, output_file)
    print(f"Processed data saved to {output_file}")

if __name__ == "__main__":
    main() 