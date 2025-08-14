#!/usr/bin/env python3
"""
Example script to run RNA clustering pipeline on the manuscript dataset.

This script demonstrates how to use the RNA clustering pipeline with the
specific data directories mentioned in the requirements.
"""

import os
import sys
from pathlib import Path

# Add utils to path
sys.path.append('utils')
from rna_clustering_pipeline import RNAClusteringPipeline

def main():
    # Data directories from the requirements
    pdb_dir = "data/experiments_data/exp2.2_manuscript/dataset_20250731/extracted_pdbs"
    fasta_dir = "data/experiments_data/exp2.2_manuscript/dataset_20250731/extracted_sequences"
    output_dir = "data/experiments_data/exp2.2_manuscript/dataset_20250731/clustering_results"
    
    # Verify directories exist
    for dir_path in [pdb_dir, fasta_dir]:
        if not os.path.exists(dir_path):
            print(f"Error: Directory does not exist: {dir_path}")
            print("Please check the paths and ensure data is available.")
            return
    
    print("Starting RNA clustering pipeline...")
    print(f"PDB directory: {pdb_dir}")
    print(f"FASTA directory: {fasta_dir}")
    print(f"Output directory: {output_dir}")
    print()
    
    # Initialize and run pipeline
    pipeline = RNAClusteringPipeline(
        pdb_dir=pdb_dir,
        fasta_dir=fasta_dir,
        output_dir=output_dir,
        max_seq_len=500  # Adjust as needed
    )
    
    try:
        pipeline.run_pipeline()
        print("\nPipeline completed successfully!")
        print(f"Results saved to: {output_dir}")
        
        # List output files
        output_path = Path(output_dir)
        if output_path.exists():
            print("\nGenerated files:")
            for file in sorted(output_path.glob("*.txt")):
                print(f"  {file.name}")
                
    except Exception as e:
        print(f"\nError running pipeline: {e}")
        print("\nTroubleshooting tips:")
        print("1. Ensure PSI-CD-HIT and US-align are installed and in PATH")
        print("2. Check that input directories contain valid PDB and FASTA files")
        print("3. Verify you have write permissions to the output directory")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())