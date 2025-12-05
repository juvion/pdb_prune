#!/usr/bin/env python3
"""
PDB Chain Extractor Script

This script extracts specific chains from PDB files based on a CSV specification file.
It reads a CSV file containing PDB IDs and chain specifications, then extracts
the specified chains from the corresponding PDB files and saves them as separate
files with their sequences in FASTA format.

Usage:
    python pdb_chain_extractor.py --csv_file validation_pdbid.csv --pdb_dir validation_downloads --output_dir validation_extracted

Author: METiS AiRNA Platform
Date: 2024
"""

import os
import sys
import argparse
import pandas as pd
from typing import Dict, List, Set, Tuple
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('pdb_chain_extraction.log')
    ]
)
logger = logging.getLogger(__name__)


class PDBChainExtractor:
    """Extracts specific chains from PDB files based on CSV specifications."""
    
    def __init__(self, csv_file: str, pdb_dir: str, output_dir: str):
        """
        Initialize the PDB chain extractor.
        
        Args:
            csv_file: Path to CSV file with PDB IDs and chain specifications
            pdb_dir: Directory containing PDB files
            output_dir: Output directory for extracted files
        """
        self.csv_file = csv_file
        self.pdb_dir = pdb_dir
        self.output_dir = output_dir
        self.pdb_chains = {}  # Dict to store chain specifications
        
        # Create output directories
        self.pdb_output_dir = os.path.join(output_dir, 'PDBs')
        self.seq_output_dir = os.path.join(output_dir, 'SEQs')
        
        os.makedirs(self.pdb_output_dir, exist_ok=True)
        os.makedirs(self.seq_output_dir, exist_ok=True)
        
        logger.info(f"Initialized PDB Chain Extractor")
        logger.info(f"CSV file: {csv_file}")
        logger.info(f"PDB directory: {pdb_dir}")
        logger.info(f"Output directory: {output_dir}")
    
    def parse_csv_file(self) -> Dict[str, List[str]]:
        """
        Parse the CSV file to extract PDB IDs and chain specifications.
        
        Returns:
            Dictionary mapping PDB IDs to lists of chain specifications
        """
        try:
            # Read CSV file
            df = pd.read_csv(self.csv_file, header=None, names=['pdb_id', 'model', 'chains'])
            
            pdb_chains = {}
            for _, row in df.iterrows():
                pdb_id = str(row['pdb_id']).strip().upper()
                chains_str = str(row['chains']).strip()
                
                # Parse chain specifications (handle ranges like "B-E", "A-B", single chains)
                chains = self._parse_chain_spec(chains_str)
                pdb_chains[pdb_id] = chains
                
                logger.debug(f"Parsed {pdb_id}: chains {chains}")
            
            logger.info(f"Successfully parsed {len(pdb_chains)} PDB entries from CSV")
            return pdb_chains
            
        except Exception as e:
            logger.error(f"Error parsing CSV file {self.csv_file}: {e}")
            raise
    
    def _parse_chain_spec(self, chains_str: str) -> List[str]:
        """
        Parse chain specification string into individual chain IDs.
        
        Args:
            chains_str: Chain specification (e.g., "B-E", "A-B", "W", "C-D")
            
        Returns:
            List of individual chain IDs
        """
        chains = []
        
        # Handle different formats
        if '-' in chains_str:
            # Handle ranges like "B-E", "A-B"
            parts = chains_str.split('-')
            if len(parts) == 2:
                start_chain = parts[0].strip()
                end_chain = parts[1].strip()
                
                # For single character chains (A-Z), generate range
                if len(start_chain) == 1 and len(end_chain) == 1:
                    start_ord = ord(start_chain)
                    end_ord = ord(end_chain)
                    if start_ord <= end_ord:
                        chains.extend([chr(i) for i in range(start_ord, end_ord + 1)])
                    else:
                        logger.warning(f"Invalid chain range: {chains_str}")
                        chains.extend([start_chain, end_chain])
                else:
                    # For multi-character chains, just add both
                    chains.extend([start_chain, end_chain])
            else:
                logger.warning(f"Complex chain specification: {chains_str}")
                chains.extend([part.strip() for part in parts])
        else:
            # Single chain or comma-separated
            if ',' in chains_str:
                chains.extend([chain.strip() for chain in chains_str.split(',')])
            else:
                chains.append(chains_str)
        
        return chains
    
    def extract_chains_from_pdb(self, pdb_file: str, target_chains: List[str]) -> Tuple[Dict[str, List[str]], Dict[str, str]]:
        """
        Extract specified chains from a PDB file.
        
        Args:
            pdb_file: Path to PDB file
            target_chains: List of chain IDs to extract
            
        Returns:
            Tuple of (pdb_records_by_chain, sequences_by_chain)
        """
        pdb_records = {}  # chain -> list of PDB records
        sequences = {}    # chain -> sequence string
        
        try:
            with open(pdb_file, 'r') as f:
                current_chain = None
                current_residue = None
                current_sequence = []
                
                for line in f:
                    line = line.rstrip()
                    
                    # Handle ATOM and HETATM records
                    if line.startswith(('ATOM', 'HETATM')):
                        if len(line) >= 22:
                            chain_id = line[21].strip()  # Chain ID is at position 22
                            residue_name = line[17:20].strip()  # Residue name
                            residue_num = line[22:26].strip()  # Residue number
                            
                            # Only process target chains
                            if chain_id in target_chains:
                                if chain_id not in pdb_records:
                                    pdb_records[chain_id] = []
                                
                                pdb_records[chain_id].append(line)
                                
                                # Build sequence (only for ATOM records, not HETATM)
                                if line.startswith('ATOM'):
                                    if current_chain != chain_id or current_residue != residue_num:
                                        # New residue
                                        if residue_name in ['A', 'U', 'G', 'C', 'T']:  # RNA nucleotides
                                            current_sequence.append(residue_name)
                                        elif residue_name in ['DA', 'DU', 'DG', 'DC', 'DT']:  # DNA nucleotides
                                            current_sequence.append(residue_name[1])  # Remove D prefix
                                        elif len(residue_name) == 1 and residue_name.isalpha():
                                            # Single letter amino acid or nucleotide
                                            current_sequence.append(residue_name)
                                        
                                        current_chain = chain_id
                                        current_residue = residue_num
                    
                    # Handle TER records
                    elif line.startswith('TER') and current_chain in target_chains:
                        if current_chain in pdb_records:
                            pdb_records[current_chain].append(line)
                    
                    # Handle END records
                    elif line.startswith('END'):
                        if current_chain in target_chains:
                            if current_chain not in sequences:
                                sequences[current_chain] = ''.join(current_sequence)
                        break
                
                # Finalize sequences
                for chain in current_sequence:
                    if current_chain not in sequences:
                        sequences[current_chain] = ''.join(current_sequence)
            
            logger.info(f"Extracted {len(pdb_records)} chains from {os.path.basename(pdb_file)}")
            return pdb_records, sequences
            
        except Exception as e:
            logger.error(f"Error extracting chains from {pdb_file}: {e}")
            return {}, {}
    
    def save_extracted_chains(self, pdb_id: str, pdb_records: Dict[str, List[str]], sequences: Dict[str, str]):
        """
        Save extracted chains as PDB and FASTA files.
        
        Args:
            pdb_id: PDB ID
            pdb_records: Dictionary of chain -> PDB records
            sequences: Dictionary of chain -> sequence
        """
        saved_files = []
        
        try:
            for chain_id, records in pdb_records.items():
                if not records:
                    continue
                
                # Create output filename
                output_pdb = os.path.join(self.pdb_output_dir, f"{pdb_id}_{chain_id}.pdb")
                output_fasta = os.path.join(self.seq_output_dir, f"{pdb_id}_{chain_id}.fasta")
                
                # Save PDB file
                with open(output_pdb, 'w') as f:
                    # Add header
                    f.write(f"HEADER    EXTRACTED CHAIN {chain_id} FROM {pdb_id}\n")
                    f.write(f"REMARK    Extracted chain {chain_id} from PDB {pdb_id}\n")
                    
                    # Write atom records
                    for record in records:
                        f.write(record + '\n')
                    
                    # Add END record
                    f.write("END\n")
                
                # Save FASTA file
                if chain_id in sequences and sequences[chain_id]:
                    with open(output_fasta, 'w') as f:
                        f.write(f">{pdb_id}_{chain_id}\n")
                        f.write(sequences[chain_id] + '\n')
                
                saved_files.extend([output_pdb, output_fasta])
                logger.info(f"Saved {pdb_id}_{chain_id} PDB and FASTA files")
            
            return saved_files
            
        except Exception as e:
            logger.error(f"Error saving extracted chains for {pdb_id}: {e}")
            return []
    
    def process_all_pdbs(self):
        """Process all PDB files according to the CSV specifications."""
        try:
            # Parse CSV file
            self.pdb_chains = self.parse_csv_file()
            
            total_processed = 0
            total_errors = 0
            
            for pdb_id, target_chains in self.pdb_chains.items():
                try:
                    # Find PDB file
                    pdb_file = os.path.join(self.pdb_dir, f"{pdb_id.lower()}.pdb")
                    
                    if not os.path.exists(pdb_file):
                        logger.warning(f"PDB file not found: {pdb_file}")
                        total_errors += 1
                        continue
                    
                    logger.info(f"Processing {pdb_id} with chains {target_chains}")
                    
                    # Extract chains
                    pdb_records, sequences = self.extract_chains_from_pdb(pdb_file, target_chains)
                    
                    if pdb_records:
                        # Save extracted chains
                        saved_files = self.save_extracted_chains(pdb_id, pdb_records, sequences)
                        total_processed += 1
                        logger.info(f"Successfully processed {pdb_id}: {len(saved_files)} files saved")
                    else:
                        logger.warning(f"No chains extracted from {pdb_id}")
                        total_errors += 1
                        
                except Exception as e:
                    logger.error(f"Error processing {pdb_id}: {e}")
                    total_errors += 1
            
            logger.info(f"Processing complete: {total_processed} successful, {total_errors} errors")
            
        except Exception as e:
            logger.error(f"Error in main processing: {e}")
            raise


def main():
    """Main function to run the PDB chain extractor."""
    parser = argparse.ArgumentParser(description='Extract specific chains from PDB files based on CSV specifications')
    parser.add_argument('--csv_file', required=True, help='CSV file with PDB IDs and chain specifications')
    parser.add_argument('--pdb_dir', required=True, help='Directory containing PDB files')
    parser.add_argument('--output_dir', required=True, help='Output directory for extracted files')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose logging')
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    try:
        # Create extractor
        extractor = PDBChainExtractor(
            csv_file=args.csv_file,
            pdb_dir=args.pdb_dir,
            output_dir=args.output_dir
        )
        
        # Process all PDBs
        extractor.process_all_pdbs()
        
        print(f"Chain extraction complete! Check {args.output_dir} for results.")
        print(f"Log file: pdb_chain_extraction.log")
        
    except Exception as e:
        logger.error(f"Fatal error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()