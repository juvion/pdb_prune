#!/usr/bin/env python3

import os
from typing import List, Dict, Set, Tuple, Optional, Any
from Bio import PDB
from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain
from Bio.PDB.Residue import Residue as BioResidue
from Bio.PDB.Atom import Atom as BioAtom
from Bio.PDB.NeighborSearch import NeighborSearch
import logging
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

class RNALoopExtractor:
    def __init__(self, distance_cutoff: float = 10.0, min_loop_size: int = 3, max_loop_size: int = 20):
        """Initialize the RNA loop extractor.
        
        Args:
            distance_cutoff (float): Distance cutoff in Angstroms for considering residues as neighbors
            min_loop_size (int): Minimum number of residues in a loop to be extracted
            max_loop_size (int): Maximum number of residues in a loop to be extracted
        """
        self.distance_cutoff = distance_cutoff
        self.min_loop_size = min_loop_size
        self.max_loop_size = max_loop_size
        self.pdb_parser = PDBParser(QUIET=True)
        
    def get_base_atom(self, residue: BioResidue) -> Optional[BioAtom]:
        """Get the N1 (for U/C) or N9 (for A/G) atom from a residue.
        
        Args:
            residue (BioResidue): Biopython Residue object
            
        Returns:
            Optional[BioAtom]: The N1 or N9 atom, or None if not found
        """
        resname = residue.get_resname()
        target_atom = 'N1' if resname in ['U', 'C'] else 'N9'
        
        # Get all atoms and find the target atom
        for atom in residue.get_atoms():
            if atom.get_name() == target_atom:
                return atom
        return None
    
    def get_sequential_residues(self, chain: Chain) -> List[Tuple[int, BioResidue]]:
        """Get residues in sequential order with their indices.
        
        Args:
            chain (Chain): Biopython Chain object
            
        Returns:
            List[Tuple[int, BioResidue]]: List of (index, residue) tuples
        """
        residues = []
        for i, residue in enumerate(chain):
            if residue.get_resname() in ['A', 'U', 'G', 'C']:
                residues.append((i, residue))
        return residues
    
    def find_neighbors(self, chain: Chain) -> Dict[int, Set[int]]:
        """Find 3D neighbors for each residue based on N1/N9 distances.
        
        Args:
            chain (Chain): Biopython Chain object
            
        Returns:
            Dict[int, Set[int]]: Dictionary mapping residue indices to sets of neighbor indices
        """
        # Get sequential residues
        seq_residues = self.get_sequential_residues(chain)
        if not seq_residues:
            return {}
        
        # Create atom list and index mapping
        atoms = []
        atom_to_index = {}
        for idx, residue in seq_residues:
            base_atom = self.get_base_atom(residue)
            if base_atom:
                atoms.append(base_atom)
                atom_to_index[base_atom] = idx
        
        # Use NeighborSearch to find neighbors
        ns = NeighborSearch(atoms)
        neighbors = {idx: set() for idx, _ in seq_residues}
        
        # Find neighbors for each atom
        for atom in atoms:
            idx = atom_to_index[atom]
            for neighbor in ns.search(atom.get_coord(), self.distance_cutoff):
                if neighbor != atom:
                    neighbors[idx].add(atom_to_index[neighbor])
        
        return neighbors
    
    def identify_loop_residues(self, neighbors: Dict[int, Set[int]]) -> Set[int]:
        """Identify residues that are part of loops.
        
        Args:
            neighbors (Dict[int, Set[int]]): Dictionary of residue neighbors
            
        Returns:
            Set[int]: Set of residue indices that are part of loops
        """
        loop_residues = set()
        
        for idx, neighbor_set in neighbors.items():
            # Check if any neighbor is more than 3 positions away
            is_loop = True
            for neighbor_idx in neighbor_set:
                if abs(idx - neighbor_idx) > 3:
                    is_loop = False
                    break
            if is_loop:
                loop_residues.add(idx)
        
        return loop_residues
    
    def extract_continuous_loops(self, loop_residues: Set[int]) -> List[Tuple[int, int]]:
        """Extract continuous segments of loop residues.
        
        Args:
            loop_residues (Set[int]): Set of loop residue indices
            
        Returns:
            List[Tuple[int, int]]: List of (start, end) indices for each continuous loop
        """
        if not loop_residues:
            return []
        
        # Sort residues
        sorted_residues = sorted(loop_residues)
        loops = []
        start = sorted_residues[0]
        prev = start
        
        for curr in sorted_residues[1:]:
            if curr != prev + 1:
                # End of current loop
                loops.append((start, prev))
                start = curr
            prev = curr
        
        # Add the last loop
        loops.append((start, prev))
        
        return loops
    
    def create_loop_structure(self, chain: Chain, start_idx: int, end_idx: int) -> Structure:
        """Create a new structure containing a loop fragment.
        
        Args:
            chain (Chain): Original chain
            start_idx (int): Start index of the loop
            end_idx (int): End index of the loop
            
        Returns:
            Structure: New structure containing the loop fragment
        """
        # Create new structure
        new_structure = Structure.Structure('RNA')
        new_model = Model.Model(0)
        new_chain = Chain.Chain(chain.id)
        
        # Add residues
        seq_residues = self.get_sequential_residues(chain)
        for idx, residue in seq_residues:
            if start_idx <= idx <= end_idx:
                new_residue = BioResidue(residue.id, residue.resname, residue.segid)
                for atom in residue:
                    new_residue.add(atom)
                new_chain.add(new_residue)
        
        new_model.add(new_chain)
        new_structure.add(new_model)
        return new_structure
    
    def extract_loops(self, pdb_path: str, output_dir: str = "extracted_loops") -> List[str]:
        """Extract loop fragments from a PDB file.
        
        Args:
            pdb_path (str): Path to the PDB file
            output_dir (str): Directory to save extracted loops
            
        Returns:
            List[str]: List of paths to extracted loop files
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
                return []
            
            # Find neighbors and identify loops
            neighbors = self.find_neighbors(rna_chain)
            loop_residues = self.identify_loop_residues(neighbors)
            continuous_loops = self.extract_continuous_loops(loop_residues)
            
            # Get sequential residues for terminal check
            seq_residues = self.get_sequential_residues(rna_chain)
            last_residue_idx = len(seq_residues) - 1
            
            # Filter loops based on size
            filtered_loops = []
            for start, end in continuous_loops:
                loop_size = end - start + 1
                if self.min_loop_size <= loop_size <= self.max_loop_size:
                    filtered_loops.append((start, end))
            
            # Extract and save loops
            output_files = []
            for start, end in filtered_loops:
                # Get original residue numbers
                start_res = seq_residues[start][1]
                end_res = seq_residues[end][1]
                start_num = start_res.get_id()[1]
                end_num = end_res.get_id()[1]
                
                # Check if this is the last residue in the chain
                is_terminal = end == last_residue_idx
                
                # Create loop structure
                loop_structure = self.create_loop_structure(rna_chain, start, end)
                
                # Save to file with residue range and 'T' suffix if terminal
                suffix = 'T' if is_terminal else ''
                output_file = os.path.join(output_dir, f"{pdb_code}_{start_num}-{end_num}{suffix}.pdb")
                io = PDBIO()
                io.set_structure(loop_structure)
                io.save(output_file)
                output_files.append(output_file)
                logging.info(f"Extracted loop from {pdb_code} (residues {start_num}-{end_num}{suffix})")
            
            return output_files
            
        except Exception as e:
            logging.error(f"Error processing {pdb_code}: {str(e)}")
            return []

def main():
    """Main function to run the RNA loop extractor."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Extract RNA loop fragments from PDB files')
    parser.add_argument('--pdb-dir', type=str, required=True,
                      help='Directory containing PDB files')
    parser.add_argument('--output-dir', type=str, default='extracted_loops',
                      help='Directory to save extracted loops (default: extracted_loops)')
    parser.add_argument('--distance-cutoff', type=float, default=10.0,
                      help='Distance cutoff in Angstroms (default: 10.0)')
    parser.add_argument('--min-loop-size', type=int, default=3,
                      help='Minimum number of residues in a loop (default: 3)')
    parser.add_argument('--max-loop-size', type=int, default=20,
                      help='Maximum number of residues in a loop (default: 20)')
    
    args = parser.parse_args()
    
    extractor = RNALoopExtractor(
        distance_cutoff=args.distance_cutoff,
        min_loop_size=args.min_loop_size,
        max_loop_size=args.max_loop_size
    )
    
    # Process all PDB files in the directory
    pdb_files = [f for f in os.listdir(args.pdb_dir) if f.endswith('.pdb')]
    total_loops = 0
    
    for pdb_file in pdb_files:
        pdb_path = os.path.join(args.pdb_dir, pdb_file)
        extracted_files = extractor.extract_loops(pdb_path, args.output_dir)
        total_loops += len(extracted_files)
    
    logging.info(f"Extracted {total_loops} loops from {len(pdb_files)} PDB files")

if __name__ == "__main__":
    main() 