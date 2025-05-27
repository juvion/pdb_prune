#!/usr/bin/env python3

import os
import sys
import logging
from pathlib import Path
from Bio import PDB
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
import re

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

class ShortRNAExtractor:
    def __init__(self, input_dir: str, output_dir: str, max_length: int = 25, max_pdbs: int = None):
        """Initialize the RNA extractor.
        
        Args:
            input_dir (str): Directory containing input PDB files
            output_dir (str): Directory to save extracted RNA PDB files
            max_length (int): Maximum sequence length to extract
            max_pdbs (int, optional): Maximum number of PDB files to generate. If None, no limit.
        """
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.max_length = max_length
        self.max_pdbs = max_pdbs
        self.parser = PDBParser(QUIET=True)
        self.io = PDBIO()
        
        # Create output directory if it doesn't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Validate input directory
        if not self.input_dir.exists():
            raise FileNotFoundError(f"Input directory not found: {input_dir}")

// ... existing code ...

    def process(self):
        """Process all PDB files in the input directory."""
        pdb_files = list(self.input_dir.glob("*.pdb")) + \
                   list(self.input_dir.glob("*.ent"))
        
        if not pdb_files:
            raise FileNotFoundError(f"No PDB files found in {self.input_dir}")
        
        logging.info(f"Found {len(pdb_files)} PDB files")
        
        total_extracted = 0
        for pdb_file in pdb_files:
            if self.max_pdbs is not None and total_extracted >= self.max_pdbs:
                logging.info(f"Reached maximum number of PDBs ({self.max_pdbs})")
                break
                
            logging.info(f"Processing {pdb_file.name}")
            extracted = self.process_pdb(pdb_file)
            total_extracted += extracted
        
        logging.info(f"\nExtraction Summary:")
        logging.info(f"Total PDB files processed: {len(pdb_files)}")
        logging.info(f"Total RNA chains extracted: {total_extracted}")
        if self.max_pdbs is not None:
            logging.info(f"Maximum PDBs limit: {self.max_pdbs}")
        logging.info(f"Output directory: {self.output_dir}")

def main():
    """Main function to run the RNA extractor."""
    if len(sys.argv) < 3:
        print("Usage: extract_short_rna.py <input_directory> <output_directory> [max_length] [max_pdbs]")
        sys.exit(1)
        
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    max_length = int(sys.argv[3]) if len(sys.argv) > 3 else 25
    max_pdbs = int(sys.argv[4]) if len(sys.argv) > 4 else None
    
    try:
        extractor = ShortRNAExtractor(input_dir, output_dir, max_length, max_pdbs)
        extractor.process()
    except Exception as e:
        logging.error(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()