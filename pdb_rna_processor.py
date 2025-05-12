#!/usr/bin/env python3

import os
import random
import json
from typing import List, Set, Dict, Optional
from Bio import PDB
from Bio.PDB import PDBParser, PPBuilder, PDBIO, Structure, Model, Chain, Residue
from Bio.PDB.Atom import Atom
import requests
from pathlib import Path

class PDBRNAProcessor:
    def __init__(self, original_dir: str = "original_pdbs", processed_dir: str = "processed_pdbs"):
        """Initialize the PDB RNA processor.
        
        Args:
            original_dir (str): Directory to store original PDB files
            processed_dir (str): Directory to store processed PDB files
        """
        self.original_dir = Path(original_dir)
        self.processed_dir = Path(processed_dir)
        self.original_dir.mkdir(exist_ok=True)
        self.processed_dir.mkdir(exist_ok=True)
        self.pdb_parser = PDBParser(QUIET=True)
        self.ppb = PPBuilder()
        
    def search_rna_structures(self, num_structures: int = 10) -> List[str]:
        """Search for RNA-containing structures using RCSB PDB API.
        
        Args:
            num_structures (int): Number of structures to retrieve
            
        Returns:
            List[str]: List of PDB IDs containing RNA
        """
        # Query to find structures containing RNA
        query = {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_polymer_entity.type",
                            "operator": "exact_match",
                            "value": "RNA"
                        }
                    }
                ]
            },
            "return_type": "entry",
            "request_options": {
                "results_content_type": ["experimental"],
                "pager": {
                    "start": 0,
                    "rows": num_structures
                }
            }
        }
        
        try:
            response = requests.post(
                "https://search.rcsb.org/rcsbsearch/v2/query",
                json=query,
                headers={"Content-Type": "application/json"}
            )
            response.raise_for_status()
            result = response.json()
            
            # Extract PDB IDs from the response
            pdb_ids = [entry["identifier"] for entry in result.get("result_set", [])]
            return pdb_ids
            
        except requests.exceptions.RequestException as e:
            print(f"Error searching for RNA structures: {e}")
            if hasattr(e.response, 'text'):
                print(f"Response: {e.response.text}")
            return []

    def download_pdb(self, pdb_id: str) -> Optional[str]:
        """Download a PDB file by its ID.
        
        Args:
            pdb_id (str): The PDB ID to download
            
        Returns:
            Optional[str]: Path to downloaded file if successful, None otherwise
        """
        pdb_id = pdb_id.lower()
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        output_path = self.original_dir / f"{pdb_id}.pdb"
        
        try:
            response = requests.get(url)
            response.raise_for_status()
            
            with open(output_path, 'w') as f:
                f.write(response.text)
            return str(output_path)
        except requests.exceptions.RequestException as e:
            print(f"Error downloading PDB {pdb_id}: {e}")
            return None

    def is_rna_chain(self, chain: Chain) -> bool:
        """Check if a chain contains RNA.
        
        Args:
            chain (Chain): Biopython Chain object
            
        Returns:
            bool: True if chain contains RNA, False otherwise
        """
        for residue in chain:
            if residue.get_resname() in ['A', 'U', 'G', 'C']:
                return True
        return False

    def get_rna_atoms(self, residue: Residue) -> List[Atom]:
        """Get the specified atoms from an RNA residue.
        
        Args:
            residue (Residue): Biopython Residue object
            
        Returns:
            List[Atom]: List of relevant atoms
        """
        backbone_atoms = ['P', 'O5\'', 'C5\'', 'C4\'', 'C3\'', 'O3\'', 'C2\'', 'O2\'', 'C1\'']
        base_atoms = ['N1'] if residue.get_resname() in ['U', 'C'] else ['N9']
        
        atoms_to_get = backbone_atoms + base_atoms
        return [atom for atom in residue if atom.get_name() in atoms_to_get]

    def get_chain_sequence(self, chain: Chain) -> str:
        """Get the sequence of a chain as a string of residue names (e.g., 'AUGC...')."""
        return ''.join(res.get_resname() for res in chain if res.get_resname() in ['A', 'U', 'G', 'C'])

    def process_pdb(self, pdb_path: str) -> List[str]:
        """Process a PDB file to extract unique RNA chains and their atoms.
        
        Args:
            pdb_path (str): Path to the PDB file
            
        Returns:
            List[str]: List of paths to processed PDB files
        """
        try:
            structure = self.pdb_parser.get_structure('RNA', pdb_path)
            processed_files = []
            seen_sequences = set()
            
            for model in structure:
                for chain in model:
                    if not self.is_rna_chain(chain):
                        continue
                    seq = self.get_chain_sequence(chain)
                    if not seq or seq in seen_sequences:
                        continue
                    seen_sequences.add(seq)
                    # Create a new structure for this chain
                    new_structure = Structure.Structure('RNA')
                    new_model = Model.Model(0)
                    new_chain = Chain.Chain(chain.id)
                    # Add only the specified atoms
                    for residue in chain:
                        new_residue = Residue.Residue(residue.id, residue.resname, residue.segid)
                        for atom in self.get_rna_atoms(residue):
                            new_residue.add(atom)
                        if len(new_residue) > 0:
                            new_chain.add(new_residue)
                    if len(new_chain) > 0:
                        new_model.add(new_chain)
                        new_structure.add(new_model)
                        output_path = self.processed_dir / f"{Path(pdb_path).stem}_{chain.id}.pdb"
                        io = PDBIO()
                        io.set_structure(new_structure)
                        io.save(str(output_path))
                        processed_files.append(str(output_path))
            return processed_files
        except Exception as e:
            print(f"Error processing PDB file {pdb_path}: {e}")
            return []

def main():
    processor = PDBRNAProcessor()
    
    # Search for RNA structures and get 10 PDB IDs
    print("Searching for RNA structures...")
    pdb_ids = processor.search_rna_structures(10)
    
    if not pdb_ids:
        print("No RNA structures found!")
        return
        
    print(f"Found {len(pdb_ids)} RNA structures: {', '.join(pdb_ids)}")
    
    # Process each PDB
    for pdb_id in pdb_ids:
        print(f"\nProcessing PDB ID: {pdb_id}")
        print("-" * 50)
        
        # Download the PDB
        pdb_path = processor.download_pdb(pdb_id)
        if not pdb_path:
            print(f"Failed to download PDB {pdb_id}")
            continue
            
        print(f"Successfully downloaded PDB to: {pdb_path}")
        
        # Process the PDB
        processed_files = processor.process_pdb(pdb_path)
        
        if not processed_files:
            print(f"No RNA chains found in {pdb_id}")
        else:
            print(f"Found {len(processed_files)} RNA chain(s):")
            for file_path in processed_files:
                print(f"- {file_path}")
                
            # Print some basic statistics
            for file_path in processed_files:
                try:
                    structure = processor.pdb_parser.get_structure('RNA', file_path)
                    chain = next(structure.get_chains())
                    num_residues = len(list(chain))
                    num_atoms = len(list(chain.get_atoms()))
                    print(f"  File: {file_path} | Residues: {num_residues} | Atoms: {num_atoms} | Chain: {chain.id}")
                except Exception as e:
                    print(f"Error analyzing {file_path}: {e}")

if __name__ == "__main__":
    main() 