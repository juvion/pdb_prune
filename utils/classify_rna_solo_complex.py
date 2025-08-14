#!/usr/bin/env python3

import os
import csv
from typing import List, Dict, Tuple
from Bio import PDB
from Bio.PDB import PDBParser, MMCIFParser, Chain
from pathlib import Path
import logging
import argparse

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

class RNAClassifier:
    def __init__(self, input_dir: str, output_dir: str):
        """Initialize the RNA classifier.
        
        Args:
            input_dir (str): Directory containing PDB files to analyze
            output_dir (str): Directory to save the classification results
        """
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.pdb_parser = PDBParser(QUIET=True)
        self.mmcif_parser = MMCIFParser(QUIET=True)
        
        # Statistics tracking
        self.stats = {
            "total_files": 0,
            "solo_rna": 0,
            "complex_rna": 0,
            "no_rna": 0,
            "invalid_files": 0
        }
        
        # Standard amino acid residues
        self.protein_residues = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
            'MSE', 'SEC', 'PYL'  # Modified amino acids
        }
        
        # Standard RNA residues
        self.rna_residues = {
            'A', 'U', 'G', 'C', 'RA', 'RU', 'RG', 'RC',
            # Common modified RNA residues
            'PSU', 'I', 'M2G', 'M7G', 'OMC', 'OMG', 'YYG', 'H2U', 'M1A',
            'M6A', 'T', 'M5C', 'M4C', 'M1G', 'M2A', 'QUO', 'YG', 'M1I',
            'M3C', 'UR3', 'M5U', 'S4U', 'M1Y', 'CBV', 'CCC', 'M7A'
        }
        
        # Standard DNA residues
        self.dna_residues = {
            'DA', 'DT', 'DG', 'DC', 'DI'
        }

    def get_parser(self, file_path: str):
        """Get the appropriate parser based on file extension."""
        if file_path.lower().endswith('.cif'):
            return self.mmcif_parser
        return self.pdb_parser

    def is_rna_chain(self, chain: Chain) -> bool:
        """Check if a chain contains RNA."""
        rna_residue_count = 0
        total_residues = 0
        
        for residue in chain:
            # Skip hetero residues (water, ions, etc.)
            if residue.id[0] != ' ':
                continue
                
            total_residues += 1
            resname = residue.get_resname().strip()
            
            if resname in self.rna_residues:
                rna_residue_count += 1
        
        # Consider it RNA if more than 50% of residues are RNA
        return total_residues > 0 and (rna_residue_count / total_residues) > 0.5

    def is_protein_chain(self, chain: Chain) -> bool:
        """Check if a chain contains protein."""
        protein_residue_count = 0
        total_residues = 0
        
        for residue in chain:
            # Skip hetero residues (water, ions, etc.)
            if residue.id[0] != ' ':
                continue
                
            total_residues += 1
            resname = residue.get_resname().strip()
            
            if resname in self.protein_residues:
                protein_residue_count += 1
        
        # Consider it protein if more than 50% of residues are protein
        return total_residues > 0 and (protein_residue_count / total_residues) > 0.5

    def is_dna_chain(self, chain: Chain) -> bool:
        """Check if a chain contains DNA."""
        dna_residue_count = 0
        total_residues = 0
        
        for residue in chain:
            # Skip hetero residues (water, ions, etc.)
            if residue.id[0] != ' ':
                continue
                
            total_residues += 1
            resname = residue.get_resname().strip()
            
            if resname in self.dna_residues:
                dna_residue_count += 1
        
        # Consider it DNA if more than 50% of residues are DNA
        return total_residues > 0 and (dna_residue_count / total_residues) > 0.5

    def get_chain_length(self, chain: Chain) -> int:
        """Get the number of standard residues in a chain (excluding hetero residues)."""
        count = 0
        for residue in chain:
            if residue.id[0] == ' ':  # Standard residue
                count += 1
        return count

    def analyze_structure(self, file_path: str) -> List[Dict]:
        """Analyze a PDB structure and classify RNA chains.
        
        Returns:
            List of dictionaries with classification results for each RNA chain
        """
        results = []
        pdb_code = Path(file_path).stem.upper()
        
        try:
            parser = self.get_parser(file_path)
            structure = parser.get_structure('structure', file_path)
            
            # Collect all chains and their types
            rna_chains = []
            protein_chains = []
            dna_chains = []
            other_chains = []
            
            for model in structure:
                for chain in model:
                    chain_length = self.get_chain_length(chain)
                    
                    # Skip very short chains (likely ligands or ions)
                    if chain_length < 3:
                        continue
                    
                    if self.is_rna_chain(chain):
                        rna_chains.append(chain)
                    elif self.is_protein_chain(chain):
                        protein_chains.append(chain)
                    elif self.is_dna_chain(chain):
                        dna_chains.append(chain)
                    else:
                        other_chains.append(chain)
            
            # Classify each RNA chain
            for rna_chain in rna_chains:
                chain_id = rna_chain.id
                pdb_code_chain_id = f"{pdb_code}_{chain_id}"
                
                # Determine if this is solo RNA
                # Solo RNA: exactly one RNA chain, no protein chains, no DNA chains
                is_solo = (
                    len(rna_chains) == 1 and 
                    len(protein_chains) == 0 and 
                    len(dna_chains) == 0
                )
                
                results.append({
                    'PDB_code': pdb_code,
                    'PDB_code_chainID': pdb_code_chain_id,
                    'is_solo': 'Y' if is_solo else 'N',
                    'chain_id': chain_id,
                    'rna_chains_count': len(rna_chains),
                    'protein_chains_count': len(protein_chains),
                    'dna_chains_count': len(dna_chains),
                    'other_chains_count': len(other_chains)
                })
                
                # Update statistics
                if is_solo:
                    self.stats["solo_rna"] += 1
                else:
                    self.stats["complex_rna"] += 1
            
            # If no RNA chains found
            if not rna_chains:
                self.stats["no_rna"] += 1
                logging.warning(f"No RNA chains found in {pdb_code}")
                
        except Exception as e:
            self.stats["invalid_files"] += 1
            logging.error(f"Error processing {file_path}: {e}")
            
        return results

    def process_all_files(self) -> str:
        """Process all PDB files in the input directory and generate classification CSV.
        
        Returns:
            Path to the output CSV file
        """
        # Find all PDB files
        pdb_files = list(self.input_dir.glob("*.pdb")) + \
                   list(self.input_dir.glob("*.ent")) + \
                   list(self.input_dir.glob("*.cif")) + \
                   list(self.input_dir.glob("*.PDB")) + \
                   list(self.input_dir.glob("*.ENT")) + \
                   list(self.input_dir.glob("*.CIF"))
        
        if not pdb_files:
            logging.error(f"No PDB files found in {self.input_dir}")
            return None
        
        logging.info(f"Found {len(pdb_files)} PDB files to process")
        
        all_results = []
        
        for i, pdb_file in enumerate(pdb_files, 1):
            self.stats["total_files"] += 1
            
            if i % 100 == 0:
                logging.info(f"Processed {i}/{len(pdb_files)} files...")
            
            results = self.analyze_structure(str(pdb_file))
            all_results.extend(results)
        
        # Write results to CSV
        output_file = self.output_dir / "rna_classification.csv"
        
        with open(output_file, 'w', newline='') as csvfile:
            fieldnames = ['PDB_code', 'PDB_code_chainID', 'is_solo']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            writer.writeheader()
            for result in all_results:
                writer.writerow({
                    'PDB_code': result['PDB_code'],
                    'PDB_code_chainID': result['PDB_code_chainID'],
                    'is_solo': result['is_solo']
                })
        
        # Write detailed results to a separate CSV for analysis
        detailed_output_file = self.output_dir / "rna_classification_detailed.csv"
        
        with open(detailed_output_file, 'w', newline='') as csvfile:
            fieldnames = [
                'PDB_code', 'PDB_code_chainID', 'is_solo', 'chain_id',
                'rna_chains_count', 'protein_chains_count', 'dna_chains_count', 'other_chains_count'
            ]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            writer.writeheader()
            for result in all_results:
                writer.writerow(result)
        
        # Print statistics
        logging.info("\nClassification Statistics:")
        logging.info(f"Total files processed: {self.stats['total_files']}")
        logging.info(f"Files with no RNA: {self.stats['no_rna']}")
        logging.info(f"Invalid files: {self.stats['invalid_files']}")
        logging.info(f"Solo RNA chains: {self.stats['solo_rna']}")
        logging.info(f"Complex RNA chains: {self.stats['complex_rna']}")
        logging.info(f"Total RNA chains analyzed: {len(all_results)}")
        
        if len(all_results) > 0:
            solo_percentage = (self.stats['solo_rna'] / len(all_results)) * 100
            logging.info(f"Solo RNA percentage: {solo_percentage:.1f}%")
        
        logging.info(f"\nResults saved to:")
        logging.info(f"Main output: {output_file}")
        logging.info(f"Detailed output: {detailed_output_file}")
        
        return str(output_file)

def download_sample_pdbs(pdb_ids_file: str, output_dir: str, max_samples: int = 10) -> List[str]:
    """Download a sample of PDB files for testing.
    
    Args:
        pdb_ids_file: Path to file containing PDB IDs
        output_dir: Directory to save downloaded PDB files
        max_samples: Maximum number of PDB files to download
        
    Returns:
        List of paths to downloaded PDB files
    """
    import requests
    from pathlib import Path
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Read PDB IDs from file
    pdb_ids = []
    try:
        with open(pdb_ids_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('pdb') and line.endswith('.ent'):
                    # Extract PDB ID from format like "pdb1al5.ent"
                    pdb_id = line[3:7]  # Extract 4-character PDB ID
                    pdb_ids.append(pdb_id)
                elif len(line) == 4 and line.isalnum():
                    # Direct 4-character PDB ID
                    pdb_ids.append(line)
    except FileNotFoundError:
        logging.error(f"PDB IDs file not found: {pdb_ids_file}")
        return []
    
    # Limit to max_samples
    pdb_ids = pdb_ids[:max_samples]
    
    downloaded_files = []
    for pdb_id in pdb_ids:
        pdb_id = pdb_id.lower()
        output_file = output_path / f"{pdb_id}.pdb"
        
        if output_file.exists():
            logging.info(f"File already exists: {output_file}")
            downloaded_files.append(str(output_file))
            continue
        
        try:
            # Try PDB format first
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            response = requests.get(url, timeout=30)
            
            if response.status_code == 200:
                with open(output_file, 'w') as f:
                    f.write(response.text)
                logging.info(f"Downloaded: {output_file}")
                downloaded_files.append(str(output_file))
            else:
                # Try CIF format
                url = f"https://files.rcsb.org/download/{pdb_id}.cif"
                response = requests.get(url, timeout=30)
                
                if response.status_code == 200:
                    output_file = output_path / f"{pdb_id}.cif"
                    with open(output_file, 'w') as f:
                        f.write(response.text)
                    logging.info(f"Downloaded: {output_file}")
                    downloaded_files.append(str(output_file))
                else:
                    logging.warning(f"Could not download {pdb_id}: HTTP {response.status_code}")
                    
        except Exception as e:
            logging.error(f"Error downloading {pdb_id}: {e}")
            continue
    
    return downloaded_files

def main():
    parser = argparse.ArgumentParser(
        description='Classify RNA structures as solo or complex based on PDB files.'
    )
    parser.add_argument(
        '--input-dir', 
        type=str, 
        default='/Users/xiaojuzhang/Dev/pdb_prune/data/raw_data/pdbs/raw_pdbs_full_download',
        help='Directory containing PDB files to analyze'
    )
    parser.add_argument(
        '--output-dir', 
        type=str, 
        default='/Users/xiaojuzhang/Dev/pdb_prune/data/experiments_data/exp3.1_is_solo',
        help='Directory to save classification results'
    )
    parser.add_argument(
        '--download-samples',
        action='store_true',
        help='Download sample PDB files for testing (uses PDB IDs from input directory)'
    )
    parser.add_argument(
        '--max-samples',
        type=int,
        default=10,
        help='Maximum number of sample PDB files to download (default: 10)'
    )
    
    args = parser.parse_args()
    
    # If download-samples is specified, download sample files first
    if args.download_samples:
        pdb_ids_file = Path(args.input_dir) / "downloaded_pdb_ids_20250731.txt"
        if not pdb_ids_file.exists():
            pdb_ids_file = Path(args.input_dir) / "found_pdb_ids.txt"
        
        if pdb_ids_file.exists():
            logging.info(f"Downloading sample PDB files from {pdb_ids_file}...")
            sample_dir = Path(args.input_dir) / "sample_pdbs"
            downloaded_files = download_sample_pdbs(str(pdb_ids_file), str(sample_dir), args.max_samples)
            
            if downloaded_files:
                logging.info(f"Downloaded {len(downloaded_files)} sample PDB files to {sample_dir}")
                # Update input directory to the sample directory
                args.input_dir = str(sample_dir)
            else:
                logging.error("No sample PDB files were downloaded")
                return
        else:
            logging.error(f"No PDB IDs file found in {args.input_dir}")
            return
    
    classifier = RNAClassifier(args.input_dir, args.output_dir)
    output_file = classifier.process_all_files()
    
    if output_file:
        logging.info(f"Classification complete. Results saved to {output_file}")
    else:
        logging.error("Classification failed.")

if __name__ == "__main__":
    main()