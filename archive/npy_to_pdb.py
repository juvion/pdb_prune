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

class NumpyToPDBConverter:
    def __init__(self, npy_dir: str = "npy_files", output_dir: str = "reconstructed_pdbs"):
        """Initialize the NumPy to PDB converter.
        
        Args:
            npy_dir (str): Directory containing NumPy files
            output_dir (str): Directory to save reconstructed PDB files
        """
        self.npy_dir = Path(npy_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Define the atoms in the correct order
        self.atom_names = ['P', 'O5\'', 'C5\'', 'C4\'', 'C3\'', 'O3\'', 'N1/N9']
        
        # Define which base uses N1 vs N9
        self.base_atoms = {
            'A': 'N9',
            'U': 'N1',
            'G': 'N9',
            'C': 'N1'
        }
    
    def read_fasta(self, fasta_path: str) -> str:
        """Read sequence from FASTA file.
        
        Args:
            fasta_path (str): Path to FASTA file
            
        Returns:
            str: RNA sequence
        """
        with open(fasta_path, 'r') as f:
            lines = f.readlines()
            # Skip header line and join sequence lines
            sequence = ''.join(line.strip() for line in lines[1:])
        return sequence
    
    def read_coordinates(self, npy_path: str) -> np.ndarray:
        """Read coordinates from NPY file.
        
        Args:
            npy_path (str): Path to NPY file
            
        Returns:
            np.ndarray: Array of coordinates with shape (sequence_length, 7, 3)
        """
        return np.load(npy_path)
    
    def create_residue(self, resname: str, res_id: int, coords: np.ndarray) -> Residue:
        """Create a residue with atoms based on coordinates.
        
        Args:
            resname (str): Residue name (A, U, G, C)
            res_id (int): Residue ID
            coords (np.ndarray): Array of coordinates for all atoms (7, 3)
            
        Returns:
            Residue: Biopython Residue object
        """
        residue = Residue((' ', res_id, ' '), resname, ' ')
        
        # Create atoms with coordinates
        for i, atom_name in enumerate(self.atom_names):
            # Replace N1/N9 with the correct atom name for this base
            if atom_name == 'N1/N9':
                atom_name = self.base_atoms[resname]
            
            # Get coordinates for this atom
            atom_coords = coords[i]
            
            # Create atom with proper parameters
            atom = Atom(
                name=atom_name,
                coord=atom_coords,
                bfactor=1.0,
                occupancy=1.0,
                altloc=' ',
                fullname=atom_name,
                serial_number=res_id,
                element=atom_name[0]
            )
            residue.add(atom)
        
        return residue
    
    def convert_single_file(self, npy_path: str, fasta_path: str) -> bool:
        """Convert a single NumPy file to PDB format.
        
        Args:
            npy_path (str): Path to NumPy file
            fasta_path (str): Path to FASTA file
            
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            # Read sequence and coordinates
            sequence = self.read_fasta(fasta_path)
            coordinates = self.read_coordinates(npy_path)
            
            # Verify dimensions
            if len(sequence) != len(coordinates):
                logging.error(f"Sequence length ({len(sequence)}) does not match number of coordinates ({len(coordinates)})")
                return False
            
            if coordinates.shape[1:] != (7, 3):
                logging.error(f"Invalid coordinate shape: {coordinates.shape}. Expected (sequence_length, 7, 3)")
                return False
            
            # Create structure
            structure = Structure.Structure('RNA')
            model = Model.Model(0)
            chain = Chain.Chain('A')
            
            # Add residues
            for i, (resname, coords) in enumerate(zip(sequence, coordinates), start=1):
                residue = self.create_residue(resname, i, coords)
                chain.add(residue)
            
            model.add(chain)
            structure.add(model)
            
            # Set output path
            output_path = self.output_dir / f"{Path(npy_path).stem}.pdb"
            
            # Save PDB file
            io = PDB.PDBIO()
            io.set_structure(structure)
            io.save(str(output_path))
            
            logging.info(f"Successfully converted {npy_path} to {output_path}")
            return True
            
        except Exception as e:
            logging.error(f"Error converting {npy_path}: {e}")
            return False

    def convert_all_npy(self) -> None:
        """Convert all NumPy files in the input directory to PDB format."""
        # Get all NPY files
        npy_files = list(self.npy_dir.glob("*.npy"))
        
        if not npy_files:
            logging.error(f"No NPY files found in {self.npy_dir}")
            return
        
        logging.info(f"Found {len(npy_files)} NPY files to process")
        
        # Process each file
        success_count = 0
        for npy_path in npy_files:
            # Get corresponding FASTA file
            fasta_path = Path("sequences") / f"{npy_path.stem}.fasta"
            if not fasta_path.exists():
                logging.warning(f"No matching FASTA file found for {npy_path.name}")
                continue
            
            # Convert file
            if self.convert_single_file(str(npy_path), str(fasta_path)):
                success_count += 1
        
        logging.info(f"\nProcessing complete. Successfully converted {success_count} files.")

def main():
    """Main function to run the NumPy to PDB converter."""
    parser = argparse.ArgumentParser(description='Convert NumPy files to PDB format')
    parser.add_argument('--npy-dir', type=str, default='npy_files',
                      help='Directory containing NumPy files')
    parser.add_argument('--output-dir', type=str, default='reconstructed_pdbs',
                      help='Directory to save reconstructed PDB files')
    
    args = parser.parse_args()
    
    try:
        converter = NumpyToPDBConverter(
            npy_dir=args.npy_dir,
            output_dir=args.output_dir
        )
        converter.convert_all_npy()
    except Exception as e:
        logging.error(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 