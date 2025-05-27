#!/usr/bin/env python3

import os
import numpy as np
import logging
from pathlib import Path
import argparse
from collections import defaultdict
from typing import Dict, List, Tuple, Set
import csv

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

class PDBAnalyzer:
    def __init__(self, coords_dir: str):
        """Initialize the PDB analyzer.
        
        Args:
            coords_dir (str): Directory containing NumPy coordinate arrays
        """
        self.coords_dir = Path(coords_dir)
        self.stats = {
            'total_pdbs': 0,
            'first_p_missing': 0,
            'first_atom_missing': 0,  # Track if first atom has NaN
            'residues_with_nan': defaultdict(int),  # Count of PDBs with X residues having NaN
            'bin_size': 5,  # Size of residue count bins
            'high_nan_pdbs': set(),  # Set of PDB codes with >20 residues having NaN
            'bin_pdbs': defaultdict(set)  # Track PDB codes for each bin
        }
    
    def analyze_coordinates(self, coords: np.ndarray) -> Tuple[bool, bool, int]:
        """Analyze coordinates for NaN values.
        
        Args:
            coords (np.ndarray): Coordinate array
            
        Returns:
            Tuple[bool, bool, int]: (has_first_p_nan, has_first_atom_nan, count_of_residues_with_nan)
        """
        # Check if first P atom has NaN coordinates
        has_first_p_nan = np.isnan(coords[0, 0]).any()
        
        # Check if first atom (any type) has NaN coordinates
        has_first_atom_nan = np.isnan(coords[0]).any()
        
        # Count residues with any NaN coordinates
        residues_with_nan = 0
        for residue_coords in coords:
            if np.isnan(residue_coords).any():
                residues_with_nan += 1
        
        return has_first_p_nan, has_first_atom_nan, residues_with_nan
    
    def get_bin_index(self, count: int) -> str:
        """Get the bin label for a residue count.
        
        Args:
            count (int): Number of residues with NaN
            
        Returns:
            str: Bin label (e.g., "1-5", "6-10", etc.)
        """
        bin_size = self.stats['bin_size']
        lower = ((count - 1) // bin_size) * bin_size + 1
        upper = lower + bin_size - 1
        return f"{lower}-{upper}"
    
    def process_file(self, npy_file: Path) -> None:
        """Process a single .npy file.
        
        Args:
            npy_file (Path): Path to NumPy array file
        """
        try:
            # Load coordinates
            coords = np.load(npy_file)
            
            # Analyze coordinates
            has_first_p_nan, has_first_atom_nan, residues_with_nan = self.analyze_coordinates(coords)
            
            # Update statistics
            self.stats['total_pdbs'] += 1
            if has_first_p_nan:
                self.stats['first_p_missing'] += 1
            if has_first_atom_nan:
                self.stats['first_atom_missing'] += 1
            
            if residues_with_nan > 0:
                bin_label = self.get_bin_index(residues_with_nan)
                self.stats['residues_with_nan'][bin_label] += 1
                self.stats['bin_pdbs'][bin_label].add(npy_file.stem)
                
                # Track PDBs with more than 20 residues having NaN
                if residues_with_nan > 20:
                    self.stats['high_nan_pdbs'].add(npy_file.stem)
            
        except Exception as e:
            logging.error(f"Error processing {npy_file.name}: {str(e)}")
    
    def process(self):
        """Process all .npy files in the directory."""
        npy_files = list(self.coords_dir.glob("*.npy"))
        
        if not npy_files:
            logging.error(f"No .npy files found in {self.coords_dir}")
            return
        
        logging.info(f"Found {len(npy_files)} PDB files to analyze")
        
        # Process each file
        for npy_file in npy_files:
            logging.info(f"Processing {npy_file.name}")
            self.process_file(npy_file)
        
        # Print summary
        self.print_summary()
        
        # Export statistics to CSV
        self.export_statistics()
    
    def print_summary(self):
        """Print analysis summary."""
        logging.info("\nPDB Quality Analysis Summary:")
        logging.info(f"Total PDBs analyzed: {self.stats['total_pdbs']}")
        logging.info(f"PDBs with missing first P atom: {self.stats['first_p_missing']} "
                    f"({self.stats['first_p_missing']/self.stats['total_pdbs']*100:.1f}%)")
        logging.info(f"PDBs with missing first atom: {self.stats['first_atom_missing']} "
                    f"({self.stats['first_atom_missing']/self.stats['total_pdbs']*100:.1f}%)")
        
        logging.info("\nDistribution of PDBs by number of residues with NaN coordinates:")
        for bin_label in sorted(self.stats['residues_with_nan'].keys(), 
                              key=lambda x: int(x.split('-')[0])):
            count = self.stats['residues_with_nan'][bin_label]
            percentage = count / self.stats['total_pdbs'] * 100
            logging.info(f"{bin_label} residues: {count} PDBs ({percentage:.1f}%)")
        
        # Print PDBs with high NaN counts
        if self.stats['high_nan_pdbs']:
            logging.info("\nPDBs with more than 20 residues having NaN coordinates:")
            for pdb_code in sorted(self.stats['high_nan_pdbs']):
                logging.info(f"- {pdb_code}")
    
    def export_statistics(self):
        """Export statistics to CSV file."""
        csv_file = "pdb_nan_statistics.csv"
        logging.info(f"\nExporting statistics to {csv_file}")
        
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['bin_range', 'pdb_counts', 'pdbs'])
            
            for bin_label in sorted(self.stats['residues_with_nan'].keys(), 
                                  key=lambda x: int(x.split('-')[0])):
                count = self.stats['residues_with_nan'][bin_label]
                pdbs = ';'.join(sorted(self.stats['bin_pdbs'][bin_label]))
                writer.writerow([bin_label, count, pdbs])

def main():
    """Main function to run the PDB analyzer."""
    parser = argparse.ArgumentParser(description='Analyze PDB quality in coordinate files')
    parser.add_argument('--coords-dir', type=str, required=True,
                      help='Directory containing NumPy coordinate arrays (.npy files)')
    
    args = parser.parse_args()
    
    try:
        analyzer = PDBAnalyzer(args.coords_dir)
        analyzer.process()
    except Exception as e:
        logging.error(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 