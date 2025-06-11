#!/usr/bin/env python3

import os
import sys
import logging
import random
import argparse
from pathlib import Path
from typing import List, Tuple, Dict, Set
from Bio import PDB
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import defaultdict

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

class RNAExtractor:
    def __init__(self, input_dir: str, min_length: int = 5, max_length: int = 20, 
                 coverage_rate: float = 0.1, generation_id: str = None):
        """Initialize the RNA segment extractor.
        
        Args:
            input_dir (str): Directory containing PDB files
            min_length (int): Minimum length of extracted RNA segments
            max_length (int): Maximum length of extracted RNA segments
            coverage_rate (float): Minimum sampling coverage rate (0-1)
            generation_id (str, optional): Unique identifier for this extraction run
        """
        self.input_dir = Path(input_dir)
        self.min_length = min_length
        self.max_length = max_length
        self.coverage_rate = coverage_rate
        
        # Create output directories
        if generation_id:
            self.output_dir = Path(f"extracted_rna_segments_{generation_id}")
        else:
            self.output_dir = Path("extracted_rna_segments")
            
        self.pdb_output_dir = self.output_dir / "extracted_rna_pdbs"
        self.fasta_output_dir = self.output_dir / "extracted_rna_seqs"
        
        # Create directories if they don't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.pdb_output_dir.mkdir(parents=True, exist_ok=True)
        self.fasta_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize Biopython objects
        self.parser = PDBParser(QUIET=True)
        self.io = PDBIO()
        
        # Statistics tracking
        self.stats = {
            'total_pdbs': 0,
            'processed_pdbs': 0,
            'total_chains': 0,
            'extracted_segments': 0,
            'errors': []
        }

    def calculate_num_samples(self, rna_size: int) -> int:
        """Calculate number of samples based on RNA size and coverage rate."""
        return int((self.max_length - self.min_length) * rna_size * self.coverage_rate)

    def identify_continuous_chains(self, structure) -> List[Tuple[str, List[PDB.Chain.Chain]]]:
        """Identify continuous RNA chains across chain breaks.
        
        Returns:
            List[Tuple[str, List[PDB.Chain.Chain]]]: List of (chain_ids, chain_objects) tuples
        """
        continuous_chains = []
        current_chain_group = []
        current_chain_ids = []
        
        # Sort chains by their first residue number
        all_chains = []
        for model in structure:
            for chain in model:
                residues = list(chain.get_residues())
                if residues:
                    first_res = residues[0].get_id()[1]  # Get residue number
                    all_chains.append((first_res, chain))
        
        all_chains.sort(key=lambda x: x[0])
        
        # Group chains that are continuous
        for _, chain in all_chains:
            if not current_chain_group:
                current_chain_group = [chain]
                current_chain_ids = [chain.id]
                continue
            
            # Check if this chain continues from the previous one
            prev_chain = current_chain_group[-1]
            prev_residues = list(prev_chain.get_residues())
            curr_residues = list(chain.get_residues())
            
            if prev_residues and curr_residues:
                prev_last = prev_residues[-1].get_id()[1]
                curr_first = curr_residues[0].get_id()[1]
                
                if curr_first == prev_last + 1:
                    current_chain_group.append(chain)
                    current_chain_ids.append(chain.id)
                else:
                    if current_chain_group:
                        continuous_chains.append((''.join(current_chain_ids), current_chain_group))
                    current_chain_group = [chain]
                    current_chain_ids = [chain.id]
        
        # Add the last group if it exists
        if current_chain_group:
            continuous_chains.append((''.join(current_chain_ids), current_chain_group))
        
        return continuous_chains

    def get_rna_sequence(self, residues: List[PDB.Residue.Residue]) -> str:
        """Extract RNA sequence from a list of residues."""
        sequence = ""
        for residue in residues:
            resname = residue.get_resname()
            if resname in ['A', 'U', 'G', 'C', 'RA', 'RU', 'RG', 'RC']:
                if resname.startswith('R'):
                    resname = resname[1]
                sequence += resname
        return sequence

    def find_continuous_segments(self, chain_group: List[PDB.Chain.Chain]) -> List[Tuple[int, int]]:
        """Find all continuous segments in a chain group.
        
        Returns:
            List[Tuple[int, int]]: List of (start_res, end_res) tuples
        """
        segments = []
        current_segment = None
        
        # Collect all residues from all chains in the group
        all_residues = []
        for chain in chain_group:
            all_residues.extend(list(chain.get_residues()))
        
        # Sort residues by their number
        all_residues.sort(key=lambda r: r.get_id()[1])
        
        # Find continuous segments
        for i, residue in enumerate(all_residues):
            res_num = residue.get_id()[1]
            
            if current_segment is None:
                current_segment = [res_num, res_num]
            elif res_num == current_segment[1] + 1:
                current_segment[1] = res_num
            else:
                segments.append(tuple(current_segment))
                current_segment = [res_num, res_num]
        
        if current_segment is not None:
            segments.append(tuple(current_segment))
        
        return segments

    def extract_random_segments(self, chain_group: List[PDB.Chain.Chain], 
                              segments: List[Tuple[int, int]], 
                              num_samples: int) -> List[Tuple[PDB.Structure.Structure, str, int, int]]:
        """Extract random segments from the chain group.
        
        Returns:
            List[Tuple[PDB.Structure.Structure, str, int, int]]: List of (structure, sequence, start_res, end_res)
        """
        # Filter segments by length
        valid_segments = [seg for seg in segments 
                         if self.min_length <= seg[1] - seg[0] + 1 <= self.max_length]
        if not valid_segments:
            return []
        
        # Generate all possible segments
        all_possible_segments = []
        for start_res, end_res in valid_segments:
            for length in range(self.min_length, min(self.max_length + 1, end_res - start_res + 2)):
                for start in range(start_res, end_res - length + 2):
                    all_possible_segments.append((start, start + length - 1))
        
        # Randomly sample without replacement
        if len(all_possible_segments) <= num_samples:
            selected_segments = all_possible_segments
        else:
            selected_segments = random.sample(all_possible_segments, num_samples)
        
        # Extract segments
        extracted_segments = []
        for start_res, end_res in selected_segments:
            # Create new structure
            new_structure = PDB.Structure.Structure("extracted")
            new_model = PDB.Model.Model(0)
            new_structure.add(new_model)
            
            # Extract residues for the segment
            segment_residues = []
            for chain in chain_group:
                for residue in chain:
                    res_num = residue.get_id()[1]
                    if start_res <= res_num <= end_res:
                        segment_residues.append(residue)
            
            # Create new chain and add residues
            new_chain = PDB.Chain.Chain("A")
            for residue in segment_residues:
                new_chain.add(residue.copy())
            new_model.add(new_chain)
            
            # Get sequence
            sequence = self.get_rna_sequence(segment_residues)
            
            extracted_segments.append((new_structure, sequence, start_res, end_res))
        
        return extracted_segments

    def process_pdb(self, pdb_file: Path) -> None:
        """Process a single PDB file."""
        try:
            structure = self.parser.get_structure(pdb_file.stem, pdb_file)
            pdb_id = pdb_file.stem[:4]  # Get 4-character PDB ID
            
            # Identify continuous chains
            continuous_chains = self.identify_continuous_chains(structure)
            self.stats['total_chains'] += len(continuous_chains)
            
            for chain_ids, chain_group in continuous_chains:
                # Find continuous segments
                segments = self.find_continuous_segments(chain_group)
                
                # Calculate total RNA size
                total_residues = sum(seg[1] - seg[0] + 1 for seg in segments)
                
                # Calculate number of samples
                num_samples = self.calculate_num_samples(total_residues)
                
                # Extract random segments
                extracted_segments = self.extract_random_segments(
                    chain_group, segments, num_samples)
                
                for i, (structure, sequence, start_res, end_res) in enumerate(extracted_segments):
                    # Check if segment ends at the last residue
                    is_terminal = any(end_res == seg[1] for seg in segments)
                    
                    # Generate output filenames
                    range_str = f"{start_res}-{end_res}{'T' if is_terminal else ''}"
                    base_name = f"{pdb_id}_{chain_ids}_{range_str}"
                    
                    pdb_filename = f"{base_name}.pdb"
                    fasta_filename = f"{base_name}.fasta"
                    
                    pdb_path = self.pdb_output_dir / pdb_filename
                    fasta_path = self.fasta_output_dir / fasta_filename
                    
                    # Save PDB file
                    self.io.set_structure(structure)
                    self.io.save(str(pdb_path))
                    
                    # Save FASTA file
                    fasta_header = f">{pdb_id}|{chain_ids}|{range_str}|Length={len(sequence)}"
                    record = SeqRecord(sequence, id=fasta_header, description="")
                    SeqIO.write(record, fasta_path, "fasta")
                    
                    self.stats['extracted_segments'] += 1
                    logging.info(f"Extracted segment from {pdb_id} chain {chain_ids}: "
                           f"{range_str} (length={len(sequence)})")
            
            self.stats['processed_pdbs'] += 1
            
        except Exception as e:
            error_msg = f"Error processing {pdb_file}: {str(e)}"
            logging.error(error_msg)
            self.stats['errors'].append(error_msg)

    def process(self):
        """Process all PDB files in the input directory."""
        pdb_files = list(self.input_dir.glob("*.pdb")) + \
                   list(self.input_dir.glob("*.ent"))
        
        self.stats['total_pdbs'] = len(pdb_files)
        logging.info(f"Found {self.stats['total_pdbs']} PDB files to process")
        
        for pdb_file in pdb_files:
            self.process_pdb(pdb_file)
        
        # Print statistics
        logging.info("\nProcessing complete!")
        logging.info(f"Processed {self.stats['processed_pdbs']} out of {self.stats['total_pdbs']} PDB files")
        logging.info(f"Found {self.stats['total_chains']} continuous chains")
        logging.info(f"Extracted {self.stats['extracted_segments']} segments")
        
        if self.stats['errors']:
            logging.warning(f"Encountered {len(self.stats['errors'])} errors:")
            for error in self.stats['errors']:
                logging.warning(error)

def main():
    parser = argparse.ArgumentParser(description='Extract RNA segments from PDB files')
    parser.add_argument('--input-dir', type=str, required=True,
                      help='Directory containing input PDB files')
    parser.add_argument('--min-length', type=int, default=5,
                      help='Minimum length of extracted segments')
    parser.add_argument('--max-length', type=int, default=20,
                      help='Maximum length of extracted segments')
    parser.add_argument('--coverage-rate', type=float, default=0.1,
                      help='Minimum sampling coverage rate (0-1)')
    parser.add_argument('--generation-id', type=str,
                      help='Optional identifier for this extraction run')
    
    args = parser.parse_args()
    
    # Validate coverage rate
    if not 0 < args.coverage_rate <= 1:
        parser.error("Coverage rate must be between 0 and 1")
    
    # Validate length parameters
    if args.min_length >= args.max_length:
        parser.error("Minimum length must be less than maximum length")
    
    try:
        extractor = RNAExtractor(
            input_dir=args.input_dir,
            min_length=args.min_length,
            max_length=args.max_length,
            coverage_rate=args.coverage_rate,
            generation_id=args.generation_id
        )
        extractor.process()
    except Exception as e:
        logging.error(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()