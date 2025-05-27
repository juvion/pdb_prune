#!/usr/bin/env python3

import os
import sys
import logging
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from Bio import PDB
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Bio.PDB.PDBIO import PDBIO
from Bio import SeqIO
import argparse

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

class RNAReconstructor:
    def __init__(self, coords_dir: str, seqs_dir: str, output_dir: str = "reconstructed_pdbs"):
        """Initialize the RNA reconstructor.
        
        Args:
            coords_dir (str): Directory containing NumPy coordinate arrays
            seqs_dir (str): Directory containing FASTA sequence files
            output_dir (str): Directory for output PDB files
        """
        self.coords_dir = Path(coords_dir)
        self.seqs_dir = Path(seqs_dir)
        self.output_dir = Path(output_dir)
        
        # Create output directory if it doesn't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize Biopython objects
        self.io = PDBIO()
        
        # Statistics tracking
        self.stats = {
            'total_pairs': 0,
            'successful': 0,
            'errors': []
        }
        
        # Atom order in the NumPy array
        self.atom_order = ['P', "O5'", "C5'", "C4'", "C3'", "O3'"]
        
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
    
    def load_data(self, npy_file: Path, fasta_file: Path) -> Tuple[np.ndarray, str]:
        """Load coordinate array and sequence data.
        
        Args:
            npy_file (Path): Path to NumPy array file
            fasta_file (Path): Path to FASTA sequence file
            
        Returns:
            Tuple[np.ndarray, str]: (coordinates array, sequence)
        """
        # Load coordinates
        coords = np.load(npy_file)
        
        # Load sequence
        with open(fasta_file) as handle:
            record = next(SeqIO.parse(handle, "fasta"))
            sequence = str(record.seq)
        
        return coords, sequence
    
    def validate_data(self, coords: np.ndarray, sequence: str) -> bool:
        """Validate that coordinates and sequence match.
        
        Args:
            coords (np.ndarray): Coordinate array
            sequence (str): RNA sequence
            
        Returns:
            bool: True if valid, False otherwise
        """
        if len(sequence) != coords.shape[0]:
            return False
        if coords.shape[1] != 7:  # 7 atoms per residue
            return False
        if coords.shape[2] != 3:  # 3 coordinates per atom
            return False
        return True
    
    def create_structure(self, coords: np.ndarray, sequence: str, 
                        base_name: str) -> Optional[Structure]:
        """Create a Biopython Structure object from coordinates and sequence.
        
        Args:
            coords (np.ndarray): Coordinate array
            sequence (str): RNA sequence
            base_name (str): Base name for the structure
            
        Returns:
            Optional[Structure]: Created structure or None if failed
        """
        try:
            # Create structure hierarchy
            structure = Structure(base_name)
            model = Model(0)
            chain = Chain('A')
            
            # Add residues
            for i, (base, residue_coords) in enumerate(zip(sequence, coords)):
                # Create residue
                residue = Residue((' ', i + 1, ' '), base, ' ')
                
                # Add backbone atoms
                for j, atom_name in enumerate(self.atom_order):
                    atom_coords = residue_coords[j]
                    # Only skip the first P atom if it has NaN coordinates
                    if i == 0 and j == 0 and np.isnan(atom_coords).any():
                        logging.warning(f"Skipping first P atom in residue {i+1} due to NaN coordinates")
                        continue
                    atom = Atom(atom_name, atom_coords, 20.0, 1.0, ' ',
                              atom_name, i + 1, atom_name[0])
                    residue.add(atom)
                
                # Add base atom (N1 for pyrimidines, N9 for purines)
                base_atom_name = 'N1' if base in ['U', 'C'] else 'N9'
                base_coords = residue_coords[6]
                base_atom = Atom(base_atom_name, base_coords, 20.0, 1.0, ' ',
                               base_atom_name, i + 1, 'N')
                residue.add(base_atom)
                
                # Add residue to chain
                chain.add(residue)
            
            # Add chain to model and structure
            model.add(chain)
            structure.add(model)
            return structure
            
        except Exception as e:
            logging.error(f"Error creating structure for {base_name}: {str(e)}")
            return None
    
    def process_file_pair(self, base_name: str, 
                         npy_file: Path, fasta_file: Path) -> None:
        """Process a single pair of .npy and .fasta files.
        
        Args:
            base_name (str): Base name of the files
            npy_file (Path): Path to NumPy array file
            fasta_file (Path): Path to FASTA sequence file
        """
        try:
            # Load data
            coords, sequence = self.load_data(npy_file, fasta_file)
            
            # Validate data
            if not self.validate_data(coords, sequence):
                error_msg = f"Data validation failed for {base_name}: " \
                           f"sequence length {len(sequence)} != " \
                           f"coordinates length {coords.shape[0]}"
                logging.error(error_msg)
                self.stats['errors'].append(error_msg)
                return
            
            # Create structure
            structure = self.create_structure(coords, sequence, base_name)
            if structure is None:
                return
            
            # Save PDB file
            output_file = self.output_dir / f"{base_name}.pdb"
            self.io.set_structure(structure)
            self.io.save(str(output_file))
            
            self.stats['successful'] += 1
            logging.info(f"Successfully reconstructed {base_name}.pdb")
            
        except Exception as e:
            error_msg = f"Error processing {base_name}: {str(e)}"
            logging.error(error_msg)
            self.stats['errors'].append(error_msg)
    
    def process(self):
        """Process all matching file pairs."""
        # Find matching files
        file_pairs = self.get_matching_files()
        self.stats['total_pairs'] = len(file_pairs)
        
        if not file_pairs:
            logging.error("No matching .npy and .fasta files found")
            return
        
        logging.info(f"Found {len(file_pairs)} matching file pairs")
        
        # Process each pair
        for base_name, (npy_file, fasta_file) in file_pairs.items():
            logging.info(f"Processing {base_name}")
            self.process_file_pair(base_name, npy_file, fasta_file)
        
        # Print summary
        logging.info("\nReconstruction Summary:")
        logging.info(f"Total file pairs found: {self.stats['total_pairs']}")
        logging.info(f"Successfully reconstructed: {self.stats['successful']}")
        if self.stats['errors']:
            logging.info("\nErrors encountered:")
            for error in self.stats['errors']:
                logging.info(f"- {error}")

def main():
    """Main function to run the RNA reconstructor."""
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Reconstruct PDB files from NumPy coordinates and FASTA sequences')
    parser.add_argument('--coords-dir', type=str, required=True,
                      help='Directory containing NumPy coordinate arrays (.npy files)')
    parser.add_argument('--seqs-dir', type=str, required=True,
                      help='Directory containing FASTA sequence files (.fasta files)')
    parser.add_argument('--output-dir', type=str, default='reconstructed_pdbs',
                      help='Directory for output PDB files (default: reconstructed_pdbs)')
    
    args = parser.parse_args()
    
    try:
        reconstructor = RNAReconstructor(args.coords_dir, args.seqs_dir, args.output_dir)
        reconstructor.process()
    except Exception as e:
        logging.error(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 