#!/usr/bin/env python3

import os
from typing import List, Set, Dict, Optional
from Bio import PDB
from Bio.PDB import PDBParser, MMCIFParser, PPBuilder, PDBIO, Structure, Model, Chain, Residue
from Bio.PDB.Atom import Atom
from pathlib import Path
import logging
import argparse

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

class PDBRNAProcessor:
    def __init__(self, original_dir: str = "original_pdbs", processed_dir: str = "processed_pdbs", fasta_dir: str = "rna_sequences"):
        """Initialize the PDB RNA processor.
        
        Args:
            original_dir (str): Directory containing original PDB/CIF files
            processed_dir (str): Directory to store processed PDB files
            fasta_dir (str): Directory to store RNA sequences in FASTA format
        """
        self.original_dir = Path(original_dir)
        self.processed_dir = Path(processed_dir)
        self.fasta_dir = Path(fasta_dir)
        self.processed_dir.mkdir(exist_ok=True)
        self.fasta_dir.mkdir(exist_ok=True)
        self.pdb_parser = PDBParser(QUIET=True)
        self.mmcif_parser = MMCIFParser(QUIET=True)
        self.ppb = PPBuilder()
        
        # Chain ID mapping for long chain IDs
        self.chain_id_map = {}
        self.next_chain_id = 'A'
        
        # Statistics tracking
        self.stats = {
            "total_files": 0,
            "no_rna_chains": 0,
            "invalid_files": 0,
            "processed_files": 0,
            "chains_processed": 0,
            "pdb_files": 0,
            "cif_files": 0
        }

    def get_parser(self, file_path: str) -> PDBParser:
        """Get the appropriate parser based on file extension.
        
        Args:
            file_path (str): Path to the structure file
            
        Returns:
            PDBParser: Appropriate parser for the file type
        """
        if file_path.lower().endswith('.cif'):
            return self.mmcif_parser
        return self.pdb_parser

    def is_rna_chain(self, chain: Chain) -> bool:
        """Check if a chain contains RNA.
        
        Args:
            chain (Chain): Biopython Chain object
            
        Returns:
            bool: True if chain contains RNA, False otherwise
        """
        for residue in chain:
            # Handle both standard and modified RNA residues
            resname = residue.get_resname().strip()
            # Check for standard RNA residues and common modifications
            if resname in ['A', 'U', 'G', 'C', 'RA', 'RU', 'RG', 'RC']:
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
        
        # Handle both standard and modified RNA residues
        resname = residue.get_resname().strip()
        if resname.startswith('R'):
            resname = resname[1]  # Remove 'R' prefix if present
        
        # Only keep N1 (for U/C) or N9 (for A/G) as base connecting atom
        base_atoms = ['N1'] if resname in ['U', 'C'] else ['N9']
        
        atoms_to_get = backbone_atoms + base_atoms
        return [atom for atom in residue if atom.get_name() in atoms_to_get]

    def get_chain_sequence(self, chain: Chain) -> str:
        """Get the sequence of a chain as a string of residue names (e.g., 'AUGC...')."""
        sequence = []
        for res in chain:
            resname = res.get_resname().strip()
            if resname.startswith('R'):
                resname = resname[1]  # Remove 'R' prefix if present
            if resname in ['A', 'U', 'G', 'C']:
                sequence.append(resname)
        return ''.join(sequence)

    def get_single_char_chain_id(self, original_chain_id: str) -> str:
        """Map a chain ID to a single character for PDB format compatibility.
        
        Args:
            original_chain_id (str): Original chain ID from mmCIF
            
        Returns:
            str: Single character chain ID
        """
        if len(original_chain_id) == 1:
            return original_chain_id
            
        if original_chain_id not in self.chain_id_map:
            self.chain_id_map[original_chain_id] = self.next_chain_id
            # Move to next available character (A-Z, then a-z)
            if self.next_chain_id == 'Z':
                self.next_chain_id = 'a'
            elif self.next_chain_id == 'z':
                self.next_chain_id = 'A'
            else:
                self.next_chain_id = chr(ord(self.next_chain_id) + 1)
                
        return self.chain_id_map[original_chain_id]

    def process_structure(self, structure_path: str, atom_option: str = "backbone") -> List[str]:
        """Process a structure file (PDB or CIF) to extract RNA chains and their atoms.
        
        Args:
            structure_path (str): Path to the structure file
            
        Returns:
            List[str]: List of paths to processed PDB files
        """
        self.stats["total_files"] += 1
        if structure_path.lower().endswith('.cif'):
            self.stats["cif_files"] += 1
        else:
            self.stats["pdb_files"] += 1
            
        processed_files = []
        
        # Validate atom_option
        if atom_option not in {"backbone", "all"}:
            logging.warning(f"Invalid atom_option '{atom_option}'. Falling back to 'backbone'.")
            atom_option = "backbone"
        
        try:
            parser = self.get_parser(structure_path)
            structure = parser.get_structure('RNA', structure_path)
            has_rna = False
            
            for model in structure:
                for chain in model:
                    if not self.is_rna_chain(chain):
                        continue
                    
                    has_rna = True
                    # Get sequence for logging
                    seq = self.get_chain_sequence(chain)
                    if not seq:
                        logging.warning(f"No sequence found in chain {chain.id} of {structure_path}")
                        continue
                    
                    # Create a new structure for this chain
                    new_structure = Structure.Structure('RNA')
                    new_model = Model.Model(0)
                    
                    # Map chain ID to single character if needed
                    new_chain_id = self.get_single_char_chain_id(chain.id)
                    new_chain = Chain.Chain(new_chain_id)
                    
                    # Add only the specified atoms
                    for residue in chain:
                        # Handle residue name formatting
                        resname = residue.get_resname().strip()
                        if resname.startswith('R'):
                            resname = resname[1]  # Remove 'R' prefix if present
                        
                        new_residue = Residue.Residue(residue.id, resname, residue.segid)
                        # Choose atoms based on option
                        if atom_option == "backbone":
                            atoms_iter = self.get_rna_atoms(residue)
                        else:  # "all"
                            atoms_iter = list(residue.get_atoms())
                        for atom in atoms_iter:
                            new_residue.add(atom)
                        if len(new_residue) > 0:  # Only add if it has atoms
                            new_chain.add(new_residue)
                    
                    if len(new_chain) > 0:  # Only add if it has residues
                        new_model.add(new_chain)
                        new_structure.add(new_model)
                        
                        # Save to file with chain ID in filename
                        output_path = self.processed_dir / f"{Path(structure_path).stem}_{new_chain_id}.pdb"
                        io = PDBIO()
                        io.set_structure(new_structure)
                        io.save(str(output_path))
                        processed_files.append(str(output_path))
                        self.stats["chains_processed"] += 1
            
            if not has_rna:
                self.stats["no_rna_chains"] += 1
                logging.info(f"No RNA chains found in {structure_path}")
            
        except Exception as e:
            self.stats["invalid_files"] += 1
            logging.error(f"Error processing structure file {structure_path}: {e}")
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

    def extract_individual_rna_chains(self, pdb_files: List[str], output_dir: str = "extracted_rna_chains", atom_option: str = "backbone") -> List[str]:
        """
        Extract individual RNA chains from specific PDB/CIF files, with atom retention options.
        
        Options:
        - "backbone": 只保留骨架原子与 N1/N9
        - "all": 保留全部原子
        
        Args:
            pdb_files (List[str]): List of paths to PDB/CIF files to process
            output_dir (str): Directory to save extracted RNA chain PDB files
            atom_option (str): Atom selection option, either "backbone" or "all"
            
        Returns:
            List[str]: List of paths to extracted RNA chain PDB files
        """
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        extracted_files: List[str] = []
        
        if atom_option not in {"backbone", "all"}:
            logging.warning(f"Invalid atom_option '{atom_option}'. Falling back to 'backbone'.")
            atom_option = "backbone"
        
        for pdb_file in pdb_files:
            pdb_path = Path(pdb_file)
            if not pdb_path.exists():
                logging.error(f"PDB file not found: {pdb_file}")
                continue
                
            logging.info(f"Processing {pdb_path.name} for RNA chain extraction (atom_option={atom_option})...")
            
            try:
                parser = self.get_parser(str(pdb_path))
                structure = parser.get_structure('RNA', str(pdb_path))
                pdb_code = pdb_path.stem
                
                rna_chains_found = 0
                
                for model in structure:
                    for chain in model:
                        if not self.is_rna_chain(chain):
                            continue
                        
                        rna_chains_found += 1
                        
                        # Get sequence for logging
                        seq = self.get_chain_sequence(chain)
                        if not seq:
                            logging.warning(f"No valid RNA sequence found in chain {chain.id} of {pdb_path.name}")
                            continue
                        
                        # Create a new structure for this RNA chain
                        new_structure = Structure.Structure('RNA')
                        new_model = Model.Model(0)
                        
                        # Use original chain ID (preserve it as is)
                        new_chain = Chain.Chain(chain.id)
                        
                        # Add RNA residues with selected atoms
                        for residue in chain:
                            resname = residue.get_resname().strip()
                            
                            # Only process RNA residues
                            if resname in ['A', 'U', 'G', 'C', 'RA', 'RU', 'RG', 'RC']:
                                # Handle residue name formatting
                                if resname.startswith('R'):
                                    resname = resname[1]  # Remove 'R' prefix if present
                                
                                new_residue = Residue.Residue(residue.id, resname, residue.segid)
                                
                                # Choose atoms based on option
                                if atom_option == "backbone":
                                    atoms_iter = self.get_rna_atoms(residue)
                                else:  # "all"
                                    atoms_iter = list(residue.get_atoms())
                                
                                for atom in atoms_iter:
                                    new_residue.add(atom)
                                
                                if len(new_residue) > 0:  # Only add if it has atoms
                                    new_chain.add(new_residue)
                        
                        if len(new_chain) > 0:  # Only save if chain has residues
                            new_model.add(new_chain)
                            new_structure.add(new_model)
                            
                            # Save with naming convention: PDBCode_chainID.pdb
                            output_filename = f"{pdb_code}_{chain.id}.pdb"
                            output_file_path = output_path / output_filename
                            
                            io = PDBIO()
                            io.set_structure(new_structure)
                            io.save(str(output_file_path))
                            
                            extracted_files.append(str(output_file_path))
                            
                            # Log extraction details
                            num_residues = len(list(new_chain))
                            num_atoms = len(list(new_chain.get_atoms()))
                            logging.info(f"  Extracted chain {chain.id} ({atom_option} atoms): {output_filename} | Sequence: {seq} | Residues: {num_residues} | Atoms: {num_atoms}")
                
                if rna_chains_found == 0:
                    logging.info(f"  No RNA chains found in {pdb_path.name}")
                else:
                    logging.info(f"  Found and extracted {rna_chains_found} RNA chain(s) from {pdb_path.name}")
                    
            except Exception as e:
                logging.error(f"Error processing {pdb_file}: {e}")
                continue
        
        logging.info(f"\nExtraction complete. Total files extracted: {len(extracted_files)}")
        logging.info(f"Output directory: {output_path}")
        
        return extracted_files


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Process PDB/CIF files to extract RNA chains and sequences.')
    parser.add_argument('--original-dir', type=str, default="original_pdbs",
                      help='Directory containing original PDB/CIF files')
    parser.add_argument('--processed-dir', type=str, default="processed_pdbs",
                      help='Directory to store processed PDB files')
    parser.add_argument('--fasta-dir', type=str, default="rna_sequences",
                      help='Directory to store RNA sequences in FASTA format')
    parser.add_argument('--atom-option', type=str, default='backbone', choices=['backbone', 'all'],
                      help='Atom selection option for extraction: "backbone" keeps backbone + N1/N9; "all" keeps all atoms')
    args = parser.parse_args()
    
    processor = PDBRNAProcessor(
        original_dir=args.original_dir,
        processed_dir=args.processed_dir,
        fasta_dir=args.fasta_dir
    )
    
    # Process all structure files in the original directory
    structure_files = list(processor.original_dir.glob("*.pdb")) + \
                     list(processor.original_dir.glob("*.ent")) + \
                     list(processor.original_dir.glob("*.cif"))
    
    if not structure_files:
        logging.error(f"No structure files found in {args.original_dir} directory!")
        return
        
    logging.info(f"Found {len(structure_files)} structure files to process (atom_option={args.atom_option})")
    
    # Process each structure file
    for structure_path in structure_files:
        logging.info(f"\nProcessing structure: {structure_path.name} (atom_option={args.atom_option})")
        logging.info("-" * 50)
        
        # Process the structure
        processed_files = processor.process_structure(str(structure_path), atom_option=args.atom_option)
        
        if not processed_files:
            logging.info(f"No RNA chains found in {structure_path.name}")
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
    logging.info(f"PDB files: {processor.stats['pdb_files']}")
    logging.info(f"CIF files: {processor.stats['cif_files']}")
    logging.info(f"Files with no RNA chains: {processor.stats['no_rna_chains']}")
    logging.info(f"Invalid files: {processor.stats['invalid_files']}")
    logging.info(f"Total RNA chains processed: {processor.stats['chains_processed']}")
    logging.info(f"Successfully processed files: {processor.stats['processed_files']}")

if __name__ == "__main__":
    main()