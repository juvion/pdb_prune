#!/usr/bin/env python3
"""
Extract RNA chains from PDB files
Filters chains by length (20-280 nucleotides)
"""

import os
from Bio import PDB
from Bio.PDB import PDBIO
import warnings
import json
warnings.filterwarnings('ignore')

def extract_rna_chains(pdb_dir, output_dir):
    """
    Extract RNA chains from PDB files
    
    Args:
        pdb_dir: Directory containing PDB files
        output_dir: Directory to save extracted RNA chains
    
    Returns:
        List of dictionaries containing RNA structure information
    """
    os.makedirs(output_dir, exist_ok=True)
    parser = PDB.PDBParser(QUIET=True)
    io = PDBIO()
    
    rna_structures = []
    
    print(f"Scanning PDB files in: {pdb_dir}")
    pdb_files = [f for f in os.listdir(pdb_dir) if f.endswith('.pdb')]
    print(f"Found {len(pdb_files)} PDB files")
    
    for pdb_file in pdb_files:
        pdb_id = pdb_file.replace('.pdb', '')
        pdb_path = os.path.join(pdb_dir, pdb_file)
        
        try:
            structure = parser.get_structure(pdb_id, pdb_path)
            
            for model in structure:
                for chain in model:
                    # Check if chain contains RNA
                    is_rna = False
                    rna_residues = []
                    
                    for residue in chain:
                        resname = residue.get_resname().strip()
                        if resname in ['A', 'U', 'G', 'C']:
                            is_rna = True
                            rna_residues.append(resname)
                    
                    if is_rna:
                        seq_length = len(rna_residues)
                        
                        # Filter by length (20-280 nucleotides)
                        if 20 <= seq_length <= 280:
                            chain_id = f"{pdb_id}_{chain.id}"
                            output_file = os.path.join(output_dir, f"{chain_id}.pdb")
                            
                            io.set_structure(chain)
                            io.save(output_file)
                            
                            rna_structures.append({
                                'id': chain_id,
                                'pdb_id': pdb_id,
                                'chain': chain.id,
                                'length': seq_length,
                                'sequence': ''.join(rna_residues),
                                'file': output_file
                            })
                            print(f"✓ Extracted: {chain_id} (length: {seq_length})")
                        else:
                            print(f"✗ Skipped: {pdb_id}_{chain.id} (length: {seq_length}, out of range)")
        
        except Exception as e:
            print(f"✗ Error processing {pdb_file}: {str(e)}")
    
    # Save metadata
    metadata_file = os.path.join(output_dir, 'rna_chains_metadata.json')
    with open(metadata_file, 'w') as f:
        json.dump(rna_structures, f, indent=2)
    print(f"\nMetadata saved to: {metadata_file}")
    
    return rna_structures

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Extract RNA chains from PDB files')
    parser.add_argument('--pdb_dir', default='./pdb_files', help='Directory containing PDB files')
    parser.add_argument('--output_dir', default='./rna_chains', help='Output directory')
    
    args = parser.parse_args()
    
    structures = extract_rna_chains(args.pdb_dir, args.output_dir)
    print(f"\n{'='*60}")
    print(f"Total RNA chains extracted: {len(structures)}")
    print(f"Output directory: {args.output_dir}")
    print(f"{'='*60}")
