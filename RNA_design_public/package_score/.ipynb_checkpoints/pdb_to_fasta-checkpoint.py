#!/usr/bin/env python3

import os
from typing import List, Dict, Set, Tuple, Optional
from Bio import PDB
from Bio.PDB import PDBParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import logging
from pathlib import Path
import argparse

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

class SequenceExtractor:
    def __init__(self):
        """Initialize the sequence extractor."""
        self.pdb_parser = PDBParser(QUIET=True)
        
    def get_sequence(self, chain) -> Optional[str]:
        """Extract RNA sequence from a chain.
        
        Args:
            chain: Biopython Chain object
            
        Returns:
            Optional[str]: RNA sequence or None if no RNA residues found
        """
        sequence = []
        for residue in chain:
            resname = residue.get_resname()
            if resname in ['A', 'U', 'G', 'C']:
                sequence.append(resname)
        
        return ''.join(sequence) if sequence else None
    
    def extract_sequence(self, pdb_path: str, output_dir: str = "sequences") -> Optional[str]:
        """Extract sequence from a PDB file and save to FASTA.
        
        Args:
            pdb_path (str): Path to the PDB file
            output_dir (str): Directory to save FASTA files
            
        Returns:
            Optional[str]: Path to the saved FASTA file or None if extraction failed
        """
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Get PDB code from filename
        pdb_code = Path(pdb_path).stem
        
        try:
            # Parse structure
            structure = self.pdb_parser.get_structure('RNA', pdb_path)
            
            # Get first RNA chain
            rna_chain = None
            for model in structure:
                for chain in model:
                    if any(res.get_resname() in ['A', 'U', 'G', 'C'] for res in chain):
                        rna_chain = chain
                        break
                if rna_chain:
                    break
            
            if not rna_chain:
                logging.warning(f"No RNA chain found in {pdb_code}")
                return None
            
            # Extract sequence
            sequence = self.get_sequence(rna_chain)
            if not sequence:
                logging.warning(f"No RNA sequence found in {pdb_code}")
                return None
            
            # Create SeqRecord
            seq_record = SeqRecord(
                Seq(sequence),
                id=pdb_code,
                description=f"RNA sequence from {pdb_code}"
            )
            
            # Save to FASTA file
            output_file = os.path.join(output_dir, f"{pdb_code}.fasta")
            SeqIO.write(seq_record, output_file, "fasta")
            logging.info(f"Saved sequence from {pdb_code} to {output_file}")
            
            return output_file
            
        except Exception as e:
            logging.error(f"Error processing {pdb_code}: {str(e)}")
            return None

def main():
    """Main function to run the sequence extractor."""
    parser = argparse.ArgumentParser(description='Extract RNA sequences from PDB files')
    parser.add_argument('--pdb-dir', type=str, required=True,
                      help='Directory containing PDB files')
    parser.add_argument('--output-dir', type=str, default='sequences',
                      help='Directory to save FASTA files (default: sequences)')
    
    args = parser.parse_args()
    
    extractor = SequenceExtractor()
    
    # Process all PDB files in the directory
    pdb_files = [f for f in os.listdir(args.pdb_dir) if f.endswith('.pdb')]
    total_sequences = 0
    
    for pdb_file in pdb_files:
        pdb_path = os.path.join(args.pdb_dir, pdb_file)
        if extractor.extract_sequence(pdb_path, args.output_dir):
            total_sequences += 1
    
    logging.info(f"Extracted {total_sequences} sequences from {len(pdb_files)} PDB files")

if __name__ == "__main__":
    main() 