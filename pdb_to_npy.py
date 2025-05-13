#!/usr/bin/env python3

import numpy as np
from Bio import PDB
from Bio.PDB import PDBParser
from pathlib import Path
from typing import List, Dict, Tuple
import os

class PDBToNumpyConverter:
    def __init__(self, processed_dir: str = "processed_pdbs", npy_dir: str = "npy_files"):
        """Initialize the PDB to NumPy converter.
        
        Args:
            processed_dir (str): Directory containing processed PDB files
            npy_dir (str): Directory to store NumPy files
        """
        self.processed_dir = Path(processed_dir)
        self.npy_dir = Path(npy_dir)
        self.npy_dir.mkdir(exist_ok=True)
        self.pdb_parser = PDBParser(QUIET=True)
        
        # Define the order of atoms we want to extract
        self.atom_order = ["P", "O5'", "C5'", "C4'", "C3'", "O3'", "N1", "N9"]
        
    def get_atom_position(self, residue: PDB.Residue, atom_name: str) -> np.ndarray:
        """Get the position of a specific atom in a residue.
        
        Args:
            residue (PDB.Residue): Biopython Residue object
            atom_name (str): Name of the atom to get position for
            
        Returns:
            np.ndarray: Array of [x, y, z] coordinates or [nan, nan, nan] if atom not found
        """
        try:
            atom = residue[atom_name]
            return np.array(atom.get_coord())
        except KeyError:
            return np.array([np.nan, np.nan, np.nan])
    
    def get_base_connecting_atom(self, residue: PDB.Residue) -> np.ndarray:
        """Get the position of the base connecting atom (N1 for U/C, N9 for A/G).
        
        Args:
            residue (PDB.Residue): Biopython Residue object
            
        Returns:
            np.ndarray: Array of [x, y, z] coordinates or [nan, nan, nan] if atom not found
        """
        resname = residue.get_resname()
        if resname in ['U', 'C']:
            return self.get_atom_position(residue, "N1")
        elif resname in ['A', 'G']:
            return self.get_atom_position(residue, "N9")
        return np.array([np.nan, np.nan, np.nan])
    
    def process_pdb(self, pdb_path: Path) -> np.ndarray:
        """Process a PDB file and convert it to a NumPy array of atom positions.
        
        Args:
            pdb_path (Path): Path to the processed PDB file
            
        Returns:
            np.ndarray: Array of shape (seq_length, 7, 3) containing atom positions
        """
        try:
            structure = self.pdb_parser.get_structure('RNA', str(pdb_path))
            chain = next(structure.get_chains())  # Get the first (and only) chain
            
            # Get sequence length
            seq_length = len(list(chain))
            
            # Initialize array with NaN values
            positions = np.full((seq_length, 7, 3), np.nan)
            
            # Fill in atom positions
            for i, residue in enumerate(chain):
                # Get backbone atoms
                for j, atom_name in enumerate(self.atom_order[:-2]):  # Exclude N1 and N9
                    positions[i, j] = self.get_atom_position(residue, atom_name)
                
                # Get base connecting atom (N1 or N9)
                positions[i, -1] = self.get_base_connecting_atom(residue)
            
            return positions
            
        except Exception as e:
            print(f"Error processing {pdb_path}: {e}")
            return None
    
    def convert_all_pdbs(self):
        """Convert all processed PDB files to NumPy arrays."""
        pdb_files = list(self.processed_dir.glob("*.pdb"))
        print(f"\nConverting {len(pdb_files)} PDB files to NumPy arrays...")
        
        for pdb_path in pdb_files:
            try:
                # Process PDB file
                positions = self.process_pdb(pdb_path)
                if positions is None:
                    continue
                
                # Save as NumPy array
                npy_path = self.npy_dir / f"{pdb_path.stem}.npy"
                np.save(str(npy_path), positions)
                print(f"Created NumPy file: {npy_path}")
                
            except Exception as e:
                print(f"Error converting {pdb_path}: {e}")
                continue

def main():
    converter = PDBToNumpyConverter()
    converter.convert_all_pdbs()

if __name__ == "__main__":
    main() 