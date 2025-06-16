#!/usr/bin/env python3

import os
import shutil
from pathlib import Path
from Bio import PDB
import logging
from tqdm import tqdm
import argparse

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def extract_rna_sequence(pdb_file):
    """
    Extract RNA sequence from a PDB file.
    Returns the sequence length and a list of chain IDs.
    """
    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('RNA', pdb_file)
        
        # Get all RNA chains
        rna_chains = []
        for model in structure:
            for chain in model:
                # Check if chain contains RNA (has C4' atoms)
                if any(res.get_resname() in ['A', 'U', 'G', 'C'] for res in chain):
                    rna_chains.append(chain)
        
        if not rna_chains:
            return 0, []
        
        # Calculate total sequence length
        total_length = sum(len(chain) for chain in rna_chains)
        chain_ids = [chain.id for chain in rna_chains]
        
        return total_length, chain_ids
    
    except Exception as e:
        logger.error(f"Error processing {pdb_file}: {str(e)}")
        return 0, []

def filter_pdbs_by_length(input_dir, output_dir, min_length=30):
    """
    Filter PDB files based on RNA sequence length and copy them to output directory.
    """
    input_path = Path(input_dir).resolve()
    output_path = Path(output_dir).resolve()
    
    # Create output directory if it doesn't exist
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Get all PDB files
    pdb_files = list(input_path.glob('*.pdb'))
    logger.info(f"Found {len(pdb_files)} PDB files to process")
    
    # Process files with progress bar
    filtered_count = 0
    skipped_count = 0
    
    for pdb_file in tqdm(pdb_files, desc="Processing PDB files"):
        seq_length, chain_ids = extract_rna_sequence(pdb_file)
        
        if seq_length > min_length:
            # Copy file to output directory
            shutil.copy2(pdb_file, output_path / pdb_file.name)
            filtered_count += 1
            logger.debug(f"Copied {pdb_file.name} (length: {seq_length}, chains: {chain_ids})")
        else:
            skipped_count += 1
            logger.debug(f"Skipped {pdb_file.name} (length: {seq_length})")
    
    logger.info(f"Processing complete:")
    logger.info(f"- Total files processed: {len(pdb_files)}")
    logger.info(f"- Files copied: {filtered_count}")
    logger.info(f"- Files skipped: {skipped_count}")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Filter PDB files based on RNA sequence length')
    parser.add_argument('--input-dir', type=str, default='reconstructed_pdbs',
                      help='Input directory containing PDB files (default: reconstructed_pdbs)')
    parser.add_argument('--output-dir', type=str, default='reconstructed_pdbs_gt30',
                      help='Output directory for filtered PDB files (default: reconstructed_pdbs_gt30)')
    parser.add_argument('--min-length', type=int, default=30,
                      help='Minimum RNA sequence length (default: 30)')
    
    args = parser.parse_args()
    
    # Get the absolute path of the script's directory
    script_dir = Path(__file__).parent.resolve()
    
    # Convert relative paths to absolute paths
    input_dir = script_dir.parent / args.input_dir
    output_dir = script_dir.parent / args.output_dir
    
    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Minimum sequence length: {args.min_length}")
    
    logger.info("Starting PDB filtering process")
    filter_pdbs_by_length(input_dir, output_dir, args.min_length)
    logger.info("PDB filtering process completed")

if __name__ == "__main__":
    main() 