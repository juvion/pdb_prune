#!/usr/bin/env python3

import os
import numpy as np
from Bio import PDB
from Bio.PDB import Structure, Model, Chain, Residue
from Bio.PDB.Atom import Atom
from pathlib import Path
import logging
import argparse
import sys

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

class SequenceExtractor:
    def __init__(self, input_dir: str = "raw_pdbs", output_dir: str = "sequences"):
        """Initialize the sequence extractor.
        
        Args:
            input_dir (str): Directory containing input PDB files
            output_dir (str): Directory to save extracted sequences
        """
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.pdb_parser = PDB.PDBParser(QUIET=True)
        
        # RNA residue mapping
        self.rna_residues = {
            'A': 'A',  # Adenosine
            'U': 'U',  # Uridine
            'G': 'G',  # Guanosine
            'C': 'C',  # Cytidine
            'I': 'I',  # Inosine
            'T': 'U',  # Thymidine (treated as U)
            'N': 'N'   # Unknown nucleotide
        }
    
    def extract_sequence(self, structure: Structure) -> str:
        """Extract RNA sequence from a structure.
        
        Args:
            structure (Structure): Input structure
            
        Returns:
            str: RNA sequence
        """
        sequence = []
        
        # Get the first chain
        chain = next(structure.get_chains())
        
        # Extract residues
        for res in chain:
            # Get residue name (3-letter code)
            resname = res.get_resname()
            
            # Convert to 1-letter code
            if resname in self.rna_residues:
                sequence.append(self.rna_residues[resname])
            else:
                # Skip non-RNA residues
                continue
        
        return ''.join(sequence)
    
    def save_sequence(self, sequence: str, output_path: Path) -> None:
        """Save a sequence to a FASTA file.
        
        Args:
            sequence (str): RNA sequence to save
            output_path (Path): Path to save the FASTA file
        """
        with open(output_path, 'w') as f:
            f.write(f">{output_path.stem}\n")
            f.write(sequence + "\n")
    
    def process_single_file(self, pdb_path: str) -> bool:
        """Process a single PDB file to extract sequence.
        
        Args:
            pdb_path (str): Path to PDB file
            
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            # Parse PDB file
            structure = self.pdb_parser.get_structure('RNA', pdb_path)
            
            # Extract sequence
            sequence = self.extract_sequence(structure)
            
            if not sequence:
                logging.warning(f"No RNA sequence found in {pdb_path}")
                return False
            
            # Set output path
            output_path = self.output_dir / f"{Path(pdb_path).stem}.fasta"
            
            # Save sequence
            self.save_sequence(sequence, output_path)
            
            logging.info(f"Successfully extracted sequence from {pdb_path} to {output_path}")
            return True
            
        except Exception as e:
            logging.error(f"Error processing {pdb_path}: {e}")
            return False
    
    def process_all_files(self) -> None:
        """Process all PDB files in the input directory."""
        # Get all PDB files
        pdb_files = list(self.input_dir.glob("*.pdb"))
        
        if not pdb_files:
            logging.error(f"No PDB files found in {self.input_dir}")
            return
        
        logging.info(f"Found {len(pdb_files)} PDB files to process")
        
        # Process each file
        success_count = 0
        for pdb_path in pdb_files:
            if self.process_single_file(str(pdb_path)):
                success_count += 1
        
        logging.info(f"\nProcessing complete. Successfully extracted sequences from {success_count} files.")
    
    def get_sequence_lengths(self) -> dict:
        """Get sequence lengths for all processed files.
        
        Returns:
            dict: Dictionary mapping PDB IDs to sequence lengths
        """
        lengths = {}
        
        # Get all FASTA files
        fasta_files = list(self.output_dir.glob("*.fasta"))
        
        for fasta_path in fasta_files:
            try:
                with open(fasta_path, 'r') as f:
                    # Skip header line
                    next(f)
                    # Get sequence length
                    sequence = next(f).strip()
                    lengths[fasta_path.stem] = len(sequence)
            except Exception as e:
                logging.error(f"Error reading {fasta_path}: {e}")
        
        return lengths

def main():
    """Main function to run the sequence extractor."""
    parser = argparse.ArgumentParser(description='Extract RNA sequences from PDB files')
    parser.add_argument('--input-dir', type=str, default='raw_pdbs',
                      help='Directory containing input PDB files')
    parser.add_argument('--output-dir', type=str, default='sequences',
                      help='Directory to save extracted sequences')
    
    args = parser.parse_args()
    
    try:
        extractor = SequenceExtractor(
            input_dir=args.input_dir,
            output_dir=args.output_dir
        )
        extractor.process_all_files()
        
        # Print sequence length statistics
        lengths = extractor.get_sequence_lengths()
        if lengths:
            print("\nSequence Length Statistics:")
            print(f"Total sequences: {len(lengths)}")
            print(f"Min length: {min(lengths.values())}")
            print(f"Max length: {max(lengths.values())}")
            print(f"Average length: {sum(lengths.values()) / len(lengths):.2f}")
            
    except Exception as e:
        logging.error(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 