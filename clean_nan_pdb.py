#!/usr/bin/env python3

import os
import numpy as np
import logging
from pathlib import Path
import argparse
import shutil
from typing import Dict, List, Tuple, Set

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

class PDBCleaner:
    def __init__(self, coords_dir: str, seqs_dir: str, 
                 new_coords_dir: str = "new_coords", 
                 new_seqs_dir: str = "new_seqs"):
        """Initialize the PDB cleaner.
        
        Args:
            coords_dir (str): Directory containing NumPy coordinate arrays
            seqs_dir (str): Directory containing FASTA sequence files
            new_coords_dir (str): Directory for cleaned coordinate files
            new_seqs_dir (str): Directory for cleaned sequence files
        """
        self.coords_dir = Path(coords_dir)
        self.seqs_dir = Path(seqs_dir)
        self.new_coords_dir = Path(new_coords_dir)
        self.new_seqs_dir = Path(new_seqs_dir)
        
        # Create output directories
        self.new_coords_dir.mkdir(parents=True, exist_ok=True)
        self.new_seqs_dir.mkdir(parents=True, exist_ok=True)
        
        # Statistics tracking
        self.stats = {
            'total_pdbs': 0,
            'valid_pdbs': 0,
            'removed_pdbs': 0,
            'errors': []
        }
    
    def check_coordinates(self, coords: np.ndarray) -> bool:
        """Check if coordinates are valid (only first P atom can have NaN).
        
        Args:
            coords (np.ndarray): Coordinate array
            
        Returns:
            bool: True if valid, False otherwise
        """
        # Check if any atom other than first P has NaN
        for i, residue_coords in enumerate(coords):
            for j, atom_coords in enumerate(residue_coords):
                # Skip the first P atom (first atom of first residue)
                if i == 0 and j == 0:
                    continue
                if np.isnan(atom_coords).any():
                    return False
        
        return True
    
    def process_file_pair(self, npy_file: Path, fasta_file: Path) -> None:
        """Process a single pair of .npy and .fasta files.
        
        Args:
            npy_file (Path): Path to NumPy array file
            fasta_file (Path): Path to FASTA sequence file
        """
        try:
            # Load coordinates
            coords = np.load(npy_file)
            
            # Check coordinates
            if self.check_coordinates(coords):
                # Copy files to new directories
                shutil.copy2(npy_file, self.new_coords_dir / npy_file.name)
                shutil.copy2(fasta_file, self.new_seqs_dir / fasta_file.name)
                self.stats['valid_pdbs'] += 1
                logging.info(f"Copied {npy_file.stem} (valid)")
            else:
                self.stats['removed_pdbs'] += 1
                logging.info(f"Skipped {npy_file.stem} (has NaN in non-first residues)")
            
        except Exception as e:
            error_msg = f"Error processing {npy_file.stem}: {str(e)}"
            logging.error(error_msg)
            self.stats['errors'].append(error_msg)
    
    def get_matching_files(self) -> Dict[str, Tuple[Path, Path]]:
        """Find matching .npy and .fasta files.
        
        Returns:
            Dict[str, Tuple[Path, Path]]: Dictionary mapping base names to (npy_file, fasta_file) tuples
        """
        npy_files = {f.stem: f for f in self.coords_dir.glob("*.npy")}
        fasta_files = {f.stem: f for f in self.seqs_dir.glob("*.fasta")}
        
        # Find common base names
        common_names = set(npy_files.keys()) & set(fasta_files.keys())
        
        return {name: (npy_files[name], fasta_files[name]) 
                for name in common_names}
    
    def process(self):
        """Process all matching file pairs."""
        # Find matching files
        file_pairs = self.get_matching_files()
        self.stats['total_pdbs'] = len(file_pairs)
        
        if not file_pairs:
            logging.error("No matching .npy and .fasta files found")
            return
        
        logging.info(f"Found {len(file_pairs)} matching file pairs")
        
        # Process each pair
        for base_name, (npy_file, fasta_file) in file_pairs.items():
            logging.info(f"Processing {base_name}")
            self.process_file_pair(npy_file, fasta_file)
        
        # Print summary
        self.print_summary()
    
    def print_summary(self):
        """Print cleaning summary."""
        logging.info("\nCleaning Summary:")
        logging.info(f"Total PDBs processed: {self.stats['total_pdbs']}")
        logging.info(f"Valid PDBs copied: {self.stats['valid_pdbs']}")
        logging.info(f"PDBs removed: {self.stats['removed_pdbs']}")
        if self.stats['errors']:
            logging.info("\nErrors encountered:")
            for error in self.stats['errors']:
                logging.info(f"- {error}")

def main():
    """Main function to run the PDB cleaner."""
    parser = argparse.ArgumentParser(description='Clean PDB files by removing those with NaN coordinates')
    parser.add_argument('--coords-dir', type=str, required=True,
                      help='Directory containing NumPy coordinate arrays (.npy files)')
    parser.add_argument('--seqs-dir', type=str, required=True,
                      help='Directory containing FASTA sequence files (.fasta files)')
    parser.add_argument('--new-coords-dir', type=str, default='new_coords',
                      help='Directory for cleaned coordinate files (default: new_coords)')
    parser.add_argument('--new-seqs-dir', type=str, default='new_seqs',
                      help='Directory for cleaned sequence files (default: new_seqs)')
    
    args = parser.parse_args()
    
    try:
        cleaner = PDBCleaner(args.coords_dir, args.seqs_dir, 
                           args.new_coords_dir, args.new_seqs_dir)
        cleaner.process()
    except Exception as e:
        logging.error(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 