#!/usr/bin/env python3

import os
from typing import List, Set, Dict, Optional
from Bio import PDB
from Bio.PDB import PDBParser, PPBuilder, PDBIO, Structure, Model, Chain, Residue
from Bio.PDB.Atom import Atom
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

class PDBRNAProcessor:
    def __init__(self, original_dir: str = "original_pdbs", processed_dir: str = "processed_pdbs", fasta_dir: str = "rna_sequences"):
        """Initialize the PDB RNA processor.
        
        Args:
            original_dir (str): Directory containing original PDB files
            processed_dir (str): Directory to store processed PDB files
            fasta_dir (str): Directory to store RNA sequences in FASTA format
        """
        self.original_dir = Path(original_dir)
        self.processed_dir = Path(processed_dir)
        self.fasta_dir = Path(fasta_dir)
        self.processed_dir.mkdir(exist_ok=True)
        self.fasta_dir.mkdir(exist_ok=True)
        self.pdb_parser = PDBParser(QUIET=True)
        self.ppb = PPBuilder()
        
        # Statistics tracking
        self.stats = {
            "total_files": 0,
            "no_rna_chains": 0,
            "invalid_files": 0,
            "processed_files": 0,
            "chains_processed": 0
        }

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
            List[Atom]: List of relevant atoms (backbone + base connecting atom)
        """
        # Only keep specified backbone atoms
        backbone_atoms = ['P', 'O5\'', 'C5\'', 'C4\'', 'C3\'', 'O3\'']
        
        # Only keep N1 (for U/C) or N9 (for A/G) as base connecting atom
        base_atoms = ['N1'] if residue.get_resname() in ['U', 'C'] else ['N9']
        
        atoms_to_get = backbone_atoms + base_atoms
        return [atom for atom in residue if atom.get_name() in atoms_to_get]

    def get_chain_sequence(self, chain: Chain) -> str:
        """Get the sequence of a chain as a string of residue names (e.g., 'AUGC...')."""
        return ''.join(res.get_resname() for res in chain if res.get_resname() in ['A', 'U', 'G', 'C'])

    def process_pdb(self, pdb_path: str) -> List[str]:
        """Process a PDB file to extract RNA chains and their atoms.
        
        Args:
            pdb_path (str): Path to the PDB file
            
        Returns:
            List[str]: List of paths to processed PDB files
        """
        self.stats["total_files"] += 1
        processed_files = []
        
        try:
            structure = self.pdb_parser.get_structure('RNA', pdb_path)
            has_rna = False
            
            for model in structure:
                for chain in model:
                    if not self.is_rna_chain(chain):
                        continue
                    
                    has_rna = True
                    # Get sequence for logging
                    seq = self.get_chain_sequence(chain)
                    if not seq:
                        logging.warning(f"No sequence found in chain {chain.id} of {pdb_path}")
                        continue
                    
                    # Create a new structure for this chain
                    new_structure = Structure.Structure('RNA')
                    new_model = Model.Model(0)
                    new_chain = Chain.Chain(chain.id)
                    
                    # Add only the specified atoms
                    for residue in chain:
                        new_residue = Residue.Residue(residue.id, residue.resname, residue.segid)
                        for atom in self.get_rna_atoms(residue):
                            new_residue.add(atom)
                        if len(new_residue) > 0:  # Only add if it has atoms
                            new_chain.add(new_residue)
                    
                    if len(new_chain) > 0:  # Only add if it has residues
                        new_model.add(new_chain)
                        new_structure.add(new_model)
                        
                        # Save to file with chain ID in filename
                        output_path = self.processed_dir / f"{Path(pdb_path).stem}_{chain.id}.pdb"
                        io = PDBIO()
                        io.set_structure(new_structure)
                        io.save(str(output_path))
                        processed_files.append(str(output_path))
                        self.stats["chains_processed"] += 1
            
            if not has_rna:
                self.stats["no_rna_chains"] += 1
                logging.info(f"No RNA chains found in {pdb_path}")
            
        except Exception as e:
            self.stats["invalid_files"] += 1
            logging.error(f"Error processing PDB file {pdb_path}: {e}")
            return []
        
        self.stats["processed_files"] += len(processed_files)
        return processed_files

    def extract_sequences_to_fasta(self) -> List[str]:
        """Extract RNA sequences from processed PDB files and save them as FASTA files.
        
        Returns:
            List[str]: List of paths to generated FASTA files
        """
        fasta_files = []
        processed_files = list(self.processed_dir.glob("*.pdb"))
        
        logging.info(f"\nExtracting sequences from {len(processed_files)} processed PDB files...")
        
        for pdb_path in processed_files:
            try:
                structure = self.pdb_parser.get_structure('RNA', str(pdb_path))
                chain = next(structure.get_chains())  # Get the first (and only) chain
                
                # Get sequence
                sequence = self.get_chain_sequence(chain)
                if not sequence:
                    continue
                
                # Create FASTA file
                fasta_path = self.fasta_dir / f"{pdb_path.stem}.fasta"
                with open(fasta_path, 'w') as f:
                    # Write header with PDB ID and chain ID
                    f.write(f">{pdb_path.stem}\n")
                    # Write sequence in blocks of 60 characters
                    for i in range(0, len(sequence), 60):
                        f.write(f"{sequence[i:i+60]}\n")
                
                fasta_files.append(str(fasta_path))
                logging.info(f"Created FASTA file: {fasta_path}")
                
            except Exception as e:
                logging.error(f"Error extracting sequence from {pdb_path}: {e}")
                continue
        
        return fasta_files

def main():
    processor = PDBRNAProcessor()
    
    # Process all PDB files in the original_pdbs directory
    pdb_files = list(processor.original_dir.glob("*.pdb")) + list(processor.original_dir.glob("*.ent"))
    
    if not pdb_files:
        logging.error("No PDB files found in original_pdbs directory!")
        return
        
    logging.info(f"Found {len(pdb_files)} PDB files to process")
    
    # Process each PDB
    for pdb_path in pdb_files:
        logging.info(f"\nProcessing PDB: {pdb_path.name}")
        logging.info("-" * 50)
        
        # Process the PDB
        processed_files = processor.process_pdb(str(pdb_path))
        
        if not processed_files:
            logging.info(f"No RNA chains found in {pdb_path.name}")
        else:
            logging.info(f"Found {len(processed_files)} RNA chain(s):")
            for file_path in processed_files:
                logging.info(f"- {file_path}")
                
            # Print some basic statistics
            for file_path in processed_files:
                try:
                    structure = processor.pdb_parser.get_structure('RNA', file_path)
                    chain = next(structure.get_chains())
                    num_residues = len(list(chain))
                    num_atoms = len(list(chain.get_atoms()))
                    logging.info(f"  File: {file_path} | Residues: {num_residues} | Atoms: {num_atoms} | Chain: {chain.id}")
                except Exception as e:
                    logging.error(f"Error analyzing {file_path}: {e}")
    
    # Extract sequences to FASTA files
    fasta_files = processor.extract_sequences_to_fasta()
    logging.info(f"\nGenerated {len(fasta_files)} FASTA files in {processor.fasta_dir}")
    
    # Print final statistics
    logging.info("\nProcessing Statistics:")
    logging.info(f"Total files processed: {processor.stats['total_files']}")
    logging.info(f"Files with no RNA chains: {processor.stats['no_rna_chains']}")
    logging.info(f"Invalid files: {processor.stats['invalid_files']}")
    logging.info(f"Total RNA chains processed: {processor.stats['chains_processed']}")
    logging.info(f"Successfully processed files: {processor.stats['processed_files']}")

if __name__ == "__main__":
    main() 