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
from typing import List
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

class LoopExtractor:
    def __init__(self, input_dir: str = "raw_pdbs", min_length: int = 3, max_length: int = 20, 
                 distance_cutoff: float = 10.0, atom_type: str = "C4'", output_dir: str = "extracted_loops"):
        """
        Initialize the LoopExtractor.
        
        Args:
            input_dir: Directory containing input PDB files
            min_length: Minimum length of loops to extract
            max_length: Maximum length of loops to extract
            distance_cutoff: Maximum distance between C4' atoms to consider residues as connected
            atom_type: Type of atom to use for distance calculations (default: C4')
            output_dir: Directory to save extracted loops
        """
        self.input_dir = Path(input_dir)
        self.min_length = min_length
        self.max_length = max_length
        self.distance_cutoff = distance_cutoff
        self.atom_type = atom_type
        self.output_dir = Path(output_dir)
        
        # Define which base uses N1 vs N9
        self.base_atoms = {
            'A': 'N9',
            'U': 'N1',
            'G': 'N9',
            'C': 'N1',
            'RA': 'N9',
            'RU': 'N1',
            'RG': 'N9',
            'RC': 'N1'
        }
        
        # Create output directories
        self.pdb_output_dir = self.output_dir / "extracted_loop_pdbs"
        self.seq_output_dir = self.output_dir / "extracted_loop_seqs"
        self.pdb_output_dir.mkdir(parents=True, exist_ok=True)
        self.seq_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize PDB parser
        self.parser = PDB.PDBParser(QUIET=True)
        
    def _is_rna_residue(self, residue):
        """Check if a residue is RNA."""
        resname = residue.get_resname()
        return resname in ['A', 'U', 'G', 'C', 'RA', 'RU', 'RG', 'RC']

    def _is_rna_chain(self, chain):
        """Check if a chain contains RNA residues."""
        return any(self._is_rna_residue(res) for res in chain)

    def _create_loop_structure(self, residues: List[PDB.Residue.Residue], chain_id: str) -> Structure:
        """Create a new structure containing the loop residues.
        
        Args:
            residues: List of residues to include in the loop
            chain_id: ID of the original chain
            
        Returns:
            Structure: New structure containing the loop
        """
        # Create new structure
        loop_structure = Structure.Structure('LOOP')
        model = Model.Model(0)
        chain = Chain.Chain(chain_id)
        
        # Add residues to the chain
        for res in residues:
            chain.add(res.copy())
            
        model.add(chain)
        loop_structure.add(model)
        return loop_structure

    def process(self):
        """Process all PDB files in the input directory."""
        pdb_files = list(self.input_dir.glob("*.pdb"))
        logging.info(f"Found {len(pdb_files)} PDB files to process")
        
        for pdb_file in pdb_files:
            try:
                logging.info(f"Processing {pdb_file.name}")
                structure = self.parser.get_structure(pdb_file.stem, pdb_file)
                
                # Process each model in the structure
                for model in structure:
                    # Process each chain in the model
                    for chain in model:
                        if not self._is_rna_chain(chain):
                            continue
                            
                        # Get all RNA residues in the chain
                        rna_residues = [r for r in chain if self._is_rna_residue(r)]
                        if len(rna_residues) < self.min_length:
                            continue
                            
                        # Find loops in the chain
                        loops = self._find_loops(rna_residues)
                        
                        # Extract and save each loop
                        for i, loop in enumerate(loops):
                            if self.min_length <= len(loop) <= self.max_length:
                                # Create loop structure using existing method
                                loop_structure = self.extract_loop(structure, loop[0].id[1], loop[-1].id[1])
                                
                                # Save PDB file
                                pdb_filename = f"{pdb_file.stem}_{chain.id}_{loop[0].id[1]}-{loop[-1].id[1]}"
                                if loop[-1].id[1] == rna_residues[-1].id[1]:
                                    pdb_filename += "T"
                                pdb_filename += ".pdb"
                                
                                io = PDB.PDBIO()
                                io.set_structure(loop_structure)
                                io.save(str(self.pdb_output_dir / pdb_filename))
                                
                                # Save FASTA file
                                seq = self.get_rna_sequence(loop)
                                seq_record = SeqRecord(
                                    seq=seq,
                                    id=f"{pdb_file.stem}|{chain.id}|{loop[0].id[1]}-{loop[-1].id[1]}|Length={len(loop)}",
                                    description=""
                                )
                                fasta_filename = pdb_filename.replace(".pdb", ".fasta")
                                SeqIO.write(seq_record, str(self.seq_output_dir / fasta_filename), "fasta")
                                
            except Exception as e:
                logging.error(f"Error processing {pdb_file.name}: {str(e)}")
                continue

    def _find_loops(self, residues):
        """Find loops in a list of residues based on distance cutoff."""
        if not residues:
            return []
            
        # Get coordinates and residue IDs
        coords = []
        res_ids = []
        for res in residues:
            resname = res.get_resname()
            if self.atom_type == "N1" or self.atom_type == "N9":
                atom_name = self.base_atoms.get(resname, self.atom_type)
            else:
                atom_name = self.atom_type
            if atom_name in res:
                coords.append(res[atom_name].get_coord())
                res_ids.append(res.id[1])
        
        if not coords:
            return []
            
        coords = np.array(coords)
        
        # Find regions where the distance between atoms of i and i+1 exceeds the cutoff
        loops = []
        start = None
        for i in range(len(coords) - 1):
            dist = np.linalg.norm(coords[i+1] - coords[i])
            if dist > self.distance_cutoff:
                if start is not None:
                    end = i
                    if self.min_length <= (end - start + 1) <= self.max_length:
                        loops.append(residues[start:end+1])
                    start = None
            else:
                if start is None:
                    start = i
                    
        # Handle last loop if it reaches the end
        if start is not None and (len(coords) - start) >= self.min_length:
            end = len(coords) - 1
            if self.min_length <= (end - start + 1) <= self.max_length:
                loops.append(residues[start:end+1])
                
        return loops

    def extract_loop(self, structure: Structure, start_res: int, end_res: int) -> Structure:
        # Create new structure
        loop_structure = Structure.Structure('LOOP')
        model = Model.Model(0)
        chain = Chain.Chain('A')
        orig_chain = next(structure.get_chains())
        for res in orig_chain:
            if start_res <= res.id[1] <= end_res:
                chain.add(res.copy())
        model.add(chain)
        loop_structure.add(model)
        return loop_structure

    def get_rna_sequence(self, residues: List[PDB.Residue.Residue]) -> str:
        """Extract RNA sequence from a list of residues.
        
        Args:
            residues: List of residues to extract sequence from
            
        Returns:
            str: RNA sequence
        """
        sequence = ""
        for residue in residues:
            resname = residue.get_resname()
            if resname in ['A', 'U', 'G', 'C', 'RA', 'RU', 'RG', 'RC']:
                if resname.startswith('R'):
                    resname = resname[1]
                sequence += resname
        return sequence

def main():
    parser = argparse.ArgumentParser(description='Extract RNA loops from PDB files using spatial distance cutoff')
    parser.add_argument('--input-dir', type=str, default='raw_pdbs', help='Directory containing input PDB files')
    parser.add_argument('--min-length', type=int, default=3, help='Minimum number of residues in a loop')
    parser.add_argument('--max-length', type=int, default=20, help='Maximum number of residues in a loop')
    parser.add_argument('--distance-cutoff', type=float, default=10.0, help='Distance cutoff (Å) for loop definition')
    parser.add_argument('--atom-type', type=str, default="C4'", help='Atom type to use for distance measurements (e.g., C4\', N1, N9)')
    parser.add_argument('--output-dir', type=str, default="extracted_loops", help='Directory to save extracted loops')
    
    args = parser.parse_args()
    try:
        extractor = LoopExtractor(
            input_dir=args.input_dir,
            min_length=args.min_length,
            max_length=args.max_length,
            distance_cutoff=args.distance_cutoff,
            atom_type=args.atom_type,
            output_dir=args.output_dir
        )
        extractor.process()
    except Exception as e:
        logging.error(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 