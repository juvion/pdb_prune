#!/usr/bin/env python3

import os
import sys
from pathlib import Path

# Add the parent directory to Python path
sys.path.append(str(Path(__file__).parent.parent))

from utils.pdb_to_fasta import SequenceExtractor
from utils.pdb_to_npy import PDBToNumpyConverter
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def process_directory(input_dir: str):
    """Process all PDB files in a directory using both extractors."""
    input_path = Path(input_dir)
    if not input_path.exists():
        logger.error(f"Directory {input_dir} does not exist")
        return
    
    # Create output directories with clear names
    base_output_dir = Path("processed_outputs")
    cutoff_name = input_path.name  # e.g., "extracted_loops_8.0"
    
    fasta_dir = base_output_dir / cutoff_name / "fasta"
    npy_dir = base_output_dir / cutoff_name / "npy"
    
    fasta_dir.mkdir(parents=True, exist_ok=True)
    npy_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Output directories created:")
    logger.info(f"FASTA files will be saved to: {fasta_dir}")
    logger.info(f"NumPy files will be saved to: {npy_dir}")
    
    # Initialize extractors
    seq_extractor = SequenceExtractor()
    npy_converter = PDBToNumpyConverter(
        processed_dir=str(input_path / "extracted_loop_pdbs"),
        npy_dir=str(npy_dir)
    )
    
    # Process all PDB files
    pdb_files = list(input_path.glob("extracted_loop_pdbs/*.pdb"))
    logger.info(f"Found {len(pdb_files)} PDB files in {input_dir}/extracted_loop_pdbs")
    
    # Extract sequences
    total_sequences = 0
    for pdb_file in pdb_files:
        if seq_extractor.extract_sequence(str(pdb_file), str(fasta_dir)):
            total_sequences += 1
    logger.info(f"Extracted {total_sequences} sequences from {len(pdb_files)} PDB files")
    
    # Convert to NumPy arrays
    npy_converter.convert_all_pdbs()
    
    logger.info(f"Processing complete for {cutoff_name}")
    logger.info(f"Results saved in: {base_output_dir / cutoff_name}")

def main():
    # Process each cutoff directory
    for cutoff in [8.0, 10.0, 12.0]:
        input_dir = f"reconstructed_pdbs_gt30/extracted_loops_{cutoff}"
        logger.info(f"\nProcessing directory: {input_dir}")
        process_directory(input_dir)
    
    logger.info("\nAll processing complete!")
    logger.info("Results are organized in the 'processed_outputs' directory:")
    logger.info("processed_outputs/")
    logger.info("├── extracted_loops_8.0/")
    logger.info("│   ├── fasta/")
    logger.info("│   └── npy/")
    logger.info("├── extracted_loops_10.0/")
    logger.info("│   ├── fasta/")
    logger.info("│   └── npy/")
    logger.info("└── extracted_loops_12.0/")
    logger.info("    ├── fasta/")
    logger.info("    └── npy/")

if __name__ == "__main__":
    main() 