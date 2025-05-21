#!/usr/bin/env python3

import os
import sys
import logging
import random
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
    def __init__(self, input_dir: str, generation_id: str, max_length: int = 25, 
                 num_extractions_per_chain: int = 5):
        """Initialize the RNA extractor.
        
        Args:
            input_dir (str): Directory containing input PDB files
            generation_id (str): Unique identifier for this extraction run
            max_length (int): Maximum length of extracted RNA segments
            num_extractions_per_chain (int): Number of segments to extract per chain
        """
        self.input_dir = Path(input_dir)
        self.generation_id = generation_id
        self.max_length = max_length
        self.num_extractions_per_chain = num_extractions_per_chain
        
        # Create output directories
        self.output_dir = Path(f"extracted_rna_segments_{generation_id}")
        self.pdb_output_dir = self.output_dir / "extracted_rna_pdbs"
        self.fasta_output = self.output_dir / "all_extracted_rna_sequences.fasta"
        
        # Create directories if they don't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.pdb_output_dir.mkdir(parents=True, exist_ok=True)
        
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

    def extract_random_segment(self, chain_group: List[PDB.Chain.Chain], 
                             segments: List[Tuple[int, int]]) -> Tuple[PDB.Structure.Structure, str, int, int]:
        """Extract a random segment from the chain group.
        
        Returns:
            Tuple[PDB.Structure.Structure, str, int, int]: (structure, sequence, start_res, end_res)
        """
        # Filter segments by length
        valid_segments = [seg for seg in segments if seg[1] - seg[0] + 1 <= self.max_length]
        if not valid_segments:
            return None, "", 0, 0
        
        # Choose a random segment
        start_res, end_res = random.choice(valid_segments)
        
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
        
        return new_structure, sequence, start_res, end_res

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
                
                # Extract random segments
                for i in range(self.num_extractions_per_chain):
                    structure, sequence, start_res, end_res = self.extract_random_segment(
                        chain_group, segments)
                    
                    if structure and sequence:
                        # Generate output filenames
                        pdb_filename = f"{self.generation_id}_{pdb_id}_{chain_ids}_{i+1}.pdb"
                        pdb_path = self.pdb_output_dir / pdb_filename
                        
                        # Save PDB file
                        self.io.set_structure(structure)
                        self.io.save(str(pdb_path))
                        
                        # Create FASTA record
                        fasta_header = f">{pdb_id}|{chain_ids}|{start_res}-{end_res}|Length={len(sequence)}"
                        record = SeqRecord(sequence, id=fasta_header, description="")
                        
                        # Append to FASTA file
                        with open(self.fasta_output, 'a') as handle:
                            SeqIO.write(record, handle, "fasta")
                        
                        self.stats['extracted_segments'] += 1
                        logging.info(f"Extracted segment from {pdb_id} chain {chain_ids}: "
                                   f"{start_res}-{end_res} (length={len(sequence)})")
            
            self.stats['processed_pdbs'] += 1
            
        except Exception as e:
            error_msg = f"Error processing {pdb_file}: {str(e)}"
            logging.error(error_msg)
            self.stats['errors'].append(error_msg)

    def process(self):
        """Process all PDB files in the input directory."""
        pdb_files = list(self.input_dir.glob("*.pdb")) + \
                   list(self.input_dir.glob("*.ent"))
        
        if not pdb_files:
            raise FileNotFoundError(f"No PDB files found in {self.input_dir}")
        
        self.stats['total_pdbs'] = len(pdb_files)
        logging.info(f"Found {len(pdb_files)} PDB files")
        
        # Clear FASTA file if it exists
        if self.fasta_output.exists():
            self.fasta_output.unlink()
        
        for pdb_file in pdb_files:
            logging.info(f"Processing {pdb_file.name}")
            self.process_pdb(pdb_file)
        
        # Print summary
        logging.info("\nExtraction Summary:")
        logging.info(f"Total PDB files found: {self.stats['total_pdbs']}")
        logging.info(f"Successfully processed: {self.stats['processed_pdbs']}")
        logging.info(f"Total chains identified: {self.stats['total_chains']}")
        logging.info(f"Total segments extracted: {self.stats['extracted_segments']}")
        if self.stats['errors']:
            logging.info("\nErrors encountered:")
            for error in self.stats['errors']:
                logging.info(f"- {error}")

def main():
    """Main function to run the RNA extractor."""
    if len(sys.argv) < 3:
        print("Usage: extract_rna_segments.py <input_directory> <generation_id> [max_length] [num_extractions_per_chain]")
        sys.exit(1)
        
    input_dir = sys.argv[1]
    generation_id = sys.argv[2]
    max_length = int(sys.argv[3]) if len(sys.argv) > 3 else 25
    num_extractions = int(sys.argv[4]) if len(sys.argv) > 4 else 5
    
    try:
        extractor = RNAExtractor(input_dir, generation_id, max_length, num_extractions)
        extractor.process()
    except Exception as e:
        logging.error(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()