#!/usr/bin/env python3
"""
Extract RNA sequences from PDB files and save to FASTA format
"""

from Bio import PDB
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os

def extract_rna_sequences(pdb_dir, output_fasta):
    """
    Extract RNA sequences from PDB files and save to FASTA
    
    Args:
        pdb_dir: Directory containing RNA PDB files
        output_fasta: Output FASTA file path
    
    Returns:
        Number of sequences extracted
    """
    parser = PDB.PDBParser(QUIET=True)
    sequences = []
    
    print(f"Extracting sequences from: {pdb_dir}")
    pdb_files = sorted([f for f in os.listdir(pdb_dir) if f.endswith('.pdb')])
    print(f"Found {len(pdb_files)} PDB files")
    
    for pdb_file in pdb_files:
        chain_id = pdb_file.replace('.pdb', '')
        pdb_path = os.path.join(pdb_dir, pdb_file)
        
        try:
            structure = parser.get_structure(chain_id, pdb_path)
            
            # Extract sequence
            seq_list = []
            for model in structure:
                for chain in model:
                    for residue in chain:
                        resname = residue.get_resname().strip()
                        if resname in ['A', 'U', 'G', 'C']:
                            seq_list.append(resname)
            
            if seq_list:
                sequence = ''.join(seq_list)
                record = SeqRecord(
                    Seq(sequence),
                    id=chain_id,
                    description=f"RNA sequence from {chain_id}, length={len(sequence)}"
                )
                sequences.append(record)
                print(f"✓ {chain_id}: {len(sequence)} nt")
        
        except Exception as e:
            print(f"✗ Error processing {pdb_file}: {str(e)}")
    
    # Write to FASTA
    SeqIO.write(sequences, output_fasta, "fasta")
    print(f"\n{'='*60}")
    print(f"Extracted {len(sequences)} sequences to {output_fasta}")
    print(f"{'='*60}")
    
    return len(sequences)

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Extract RNA sequences to FASTA')
    parser.add_argument('--pdb_dir', default='./rna_chains', help='Directory containing RNA PDB files')
    parser.add_argument('--output', default='rna_sequences.fasta', help='Output FASTA file')
    
    args = parser.parse_args()
    
    extract_rna_sequences(args.pdb_dir, args.output)
