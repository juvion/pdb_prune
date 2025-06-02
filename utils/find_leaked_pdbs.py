#!/usr/bin/env python3

import os
from pathlib import Path
import shutil
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def extract_pdb_codes(coords_dir: str) -> set:
    """Extract PDB codes from .npy files in the coords directory."""
    coords_dir = Path(coords_dir)
    pdb_codes = set()
    
    # Scan for .npy files
    for npy_file in coords_dir.glob("*.npy"):
        # Extract first 4 characters as PDB code
        pdb_code = npy_file.stem[:4].lower()
        pdb_codes.add(pdb_code)
    
    return pdb_codes

def find_and_move_leaked_pdbs(pdb_codes: set, source_dir: str, target_dir: str):
    """Find matching .ent files and move them to target directory."""
    source_dir = Path(source_dir)
    target_dir = Path(target_dir)
    
    # Create target directory if it doesn't exist
    target_dir.mkdir(parents=True, exist_ok=True)
    
    moved_files = []
    for pdb_code in pdb_codes:
        # Construct the expected filename
        ent_file = source_dir / f"pdb{pdb_code}.ent"
        
        if ent_file.exists():
            # Move file to target directory
            target_file = target_dir / ent_file.name
            shutil.move(str(ent_file), str(target_file))
            moved_files.append(ent_file.name)
            logging.info(f"Moved {ent_file.name} to {target_dir}")
    
    return moved_files

def main():
    # Define directories
    coords_dir = "competition/train/coords"
    source_dir = "downloaded_rna_pdbs"
    target_dir = "downloaded_rna_pdbs/info_leaked_in_train"
    
    # Extract PDB codes from coords directory
    logging.info(f"Scanning {coords_dir} for PDB codes...")
    pdb_codes = extract_pdb_codes(coords_dir)
    logging.info(f"Found {len(pdb_codes)} unique PDB codes")
    
    # Find and move leaked PDB files
    logging.info(f"Searching for matching .ent files in {source_dir}...")
    moved_files = find_and_move_leaked_pdbs(pdb_codes, source_dir, target_dir)
    
    # Print summary
    logging.info("\nSummary:")
    logging.info(f"Total PDB codes found: {len(pdb_codes)}")
    logging.info(f"Files moved to {target_dir}: {len(moved_files)}")
    if moved_files:
        logging.info("\nMoved files:")
        for file in sorted(moved_files):
            logging.info(f"- {file}")

if __name__ == "__main__":
    main() 