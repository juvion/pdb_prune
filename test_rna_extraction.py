#!/usr/bin/env python3

import sys
import os
from pathlib import Path

# Add utils directory to path so we can import the processor
sys.path.append(str(Path(__file__).parent / "utils"))

from pdb_rna_processor import PDBRNAProcessor

def main():
    # Initialize the processor
    processor = PDBRNAProcessor()
    
    # Define the PDB files to process
    pdb_files = [
        "/Users/xiaojuzhang/Github/pdb_prune/raw_pdbs/2xyz.pdb",
        "/Users/xiaojuzhang/Github/pdb_prune/raw_pdbs/3lob.pdb"
    ]
    
    # Extract individual RNA chains
    output_dir = "extracted_rna_chains"
    extracted_files = processor.extract_individual_rna_chains(pdb_files, output_dir)
    
    print(f"\nExtraction Summary:")
    print(f"Total extracted files: {len(extracted_files)}")
    print(f"Output directory: {output_dir}")
    
    if extracted_files:
        print("\nExtracted files:")
        for file_path in extracted_files:
            print(f"  - {file_path}")
    else:
        print("No RNA chains were extracted.")

if __name__ == "__main__":
    main()