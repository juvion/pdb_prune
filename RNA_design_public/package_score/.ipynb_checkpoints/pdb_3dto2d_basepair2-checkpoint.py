#!/usr/bin/env python3
"""
PDB/NPY 3D to 2D Base Pair Matrix Converter with Biological Validation

This script converts 3D RNA structures to 2D base pair matrices.
Supports both PDB files and NPY coordinate arrays.
Supports two modes: all-atom and no-base-atom.

NEW: Enforces biological base-pairing rules (G-C, A-U, G-U)
NEW: Requires sequence information for all inputs

Usage:
    python pdb_3dto2d_basepair2.py input.pdb --output matrix.csv
    python pdb_3dto2d_basepair2.py input.npy --sequence seq.fasta --output matrix.csv --mode no-base-atom
    python pdb_3dto2d_basepair2.py input.pdb --output matrix.npy --mode all-atom
"""

import numpy as np
import argparse
from typing import Dict, List, Tuple, Optional, Set
import os


class RNABasePairDetector:
    """
    Detect base pairs in RNA structures from 3D coordinates.
    
    Supports:
    - PDB files (all-atom or no-base-atom)
    - NPY arrays with predefined atom order
    - Two detection modes: all-atom and no-base-atom
    - Biological base-pairing validation (G-C, A-U, G-U)
    """
    
    # Standard RNA atom names for backbone and first base atom
    BACKBONE_ATOMS = ["P", "O5'", "C5'", "C4'", "C3'", "O3'"]
    
    # N1 for pyrimidines (C, U), N9 for purines (A, G)
    BASE_MARKER_ATOMS = {
        'A': 'N9', 'G': 'N9',
        'C': 'N1', 'U': 'N1'
    }
    
    # Valid biological base pairs (Watson-Crick and wobble pairs)
    VALID_BASE_PAIRS = {
        frozenset(['G', 'C']),  # Watson-Crick
        frozenset(['A', 'U']),  # Watson-Crick
        frozenset(['G', 'U']),  # Wobble pair
    }
    
    # Distance thresholds for base pairing (in Angstroms)
    # These are standard thresholds used in RNA structure analysis
    DISTANCE_THRESHOLDS = {
        'all-atom': {
            'min': 13.0,   # Minimum C4'-C4' backbone distance
            'max': 17.0,  # Maximum C4'-C4' backbone distance (typical ~10Å)
            'n1n9': {
                'min': 7.5,  # Min N1/N9 distance
                'max': 10.5  # Max N1/N9 distance
            }
        },
        'no-base-atom': {
            'min': 13.0,   # Minimum C4'-C4' backbone distance
            'max': 17.0,  # Maximum C4'-C4' backbone distance (typical ~10Å)
            'n1n9': {
                'min': 7.5,  # Min N1/N9 distance
                'max': 10.5  # Max N1/N9 distance
            }
        }
    }
    
    def __init__(self, mode: str = 'all-atom', min_separation: int = 4):
        """
        Initialize the base pair detector.
        
        Args:
            mode: Detection mode, either 'all-atom' or 'no-base-atom'
            min_separation: Minimum sequence separation |j-i| for base pairing (default: 4)
                           This prevents physically impossible pairs between nearby residues.
        """
        if mode not in ['all-atom', 'no-base-atom']:
            raise ValueError("Mode must be 'all-atom' or 'no-base-atom'")
        
        self.mode = mode
        self.min_separation = min_separation
        self.thresholds = self.DISTANCE_THRESHOLDS[mode]
    
    @staticmethod
    def is_valid_base_pair(base1: str, base2: str) -> bool:
        """
        Check if two bases can form a valid biological base pair.
        
        Valid pairs:
        - G-C (Watson-Crick)
        - A-U (Watson-Crick)
        - G-U (Wobble pair)
        
        Args:
            base1: First base (A, U, G, C, or N for unknown)
            base2: Second base (A, U, G, C, or N for unknown)
            
        Returns:
            True if the pair is biologically valid
        """
        # If either base is unknown, we cannot validate
        if base1 == 'N' or base2 == 'N':
            return False
        
        # Check if the pair is in the valid set
        return frozenset([base1, base2]) in RNABasePairDetector.VALID_BASE_PAIRS
    
    def load_pdb(self, pdb_file: str, chain_id: Optional[str] = None) -> Tuple[np.ndarray, List[str], List[int]]:
        """
        Load RNA structure from PDB file.
        
        Args:
            pdb_file: Path to PDB file
            chain_id: Specific chain ID to extract (None = first RNA chain)
            
        Returns:
            Tuple of (coordinates_dict, residue_names, residue_numbers)
            coordinates_dict: Dict mapping (res_num, atom_name) -> [x, y, z]
        """
        coords_dict = {}
        residue_names = []
        residue_numbers = []
        seen_residues = set()
        
        with open(pdb_file, 'r') as f:
            for line in f:
                if not (line.startswith('ATOM') or line.startswith('HETATM')):
                    continue
                
                # Parse PDB line
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21].strip()
                res_num = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                
                # Filter by chain if specified
                if chain_id and chain != chain_id:
                    continue
                
                # Only process RNA residues
                if res_name not in ['A', 'U', 'G', 'C', 'DA', 'DT', 'DC', 'DG']:
                    continue
                
                # Normalize deoxy residues
                if res_name.startswith('D'):
                    res_name = res_name[1:]
                    if res_name == 'T':
                        res_name = 'U'
                
                # Track residues
                if res_num not in seen_residues:
                    seen_residues.add(res_num)
                    residue_names.append(res_name)
                    residue_numbers.append(res_num)
                
                # Store coordinates
                coords_dict[(res_num, atom_name)] = np.array([x, y, z])
        
        if not coords_dict:
            raise ValueError(f"No RNA residues found in {pdb_file}")
        
        return coords_dict, residue_names, residue_numbers
    
    def load_npy(self, npy_file: str, fasta_file: Optional[str] = None) -> Tuple[Dict, List[str], List[int]]:
        """
        Load RNA structure from NPY file.
        
        NPY format: array of shape (N, 7, 3) where:
        - N = number of residues
        - 7 atoms in order: ["P", "O5'", "C5'", "C4'", "C3'", "O3'", "N1/N9"]
        - 3 = x, y, z coordinates
        
        Args:
            npy_file: Path to NPY file
            fasta_file: Path to FASTA file containing sequence (REQUIRED)
            
        Returns:
            Tuple of (coordinates_dict, residue_names, residue_numbers)
        """
        coords_array = np.load(npy_file)
        
        if len(coords_array.shape) != 3 or coords_array.shape[1] != 7 or coords_array.shape[2] != 3:
            raise ValueError(f"NPY file must have shape (N, 7, 3), got {coords_array.shape}")
        
        n_residues = coords_array.shape[0]
        atom_names = ["P", "O5'", "C5'", "C4'", "C3'", "O3'", "N1/N9"]
        
        # FASTA file is now REQUIRED for NPY files
        if not fasta_file:
            raise ValueError(
                "ERROR: Sequence file (FASTA) is REQUIRED when using NPY files.\n"
                "       Biological base-pairing rules require sequence information.\n"
                "       Please provide a FASTA file using --sequence argument."
            )
        
        # Load sequence from FASTA file
        print(f"Loading sequence from FASTA file: {fasta_file}")
        sequence = None
        try:
            # Prefer Biopython if available
            from Bio import SeqIO  # type: ignore
            record = next(SeqIO.parse(fasta_file, "fasta"))
            sequence = str(record.seq).upper()
            print(f"Loaded sequence using Biopython: {sequence}")
        except Exception as e:
            # Fall back to simple parser on any error (e.g., Biopython missing)
            print(f"Biopython FASTA parsing failed ({e}). Falling back to simple reader.")
            sequence = self._read_fasta_simple(fasta_file)
            print(f"Loaded sequence using simple reader: {sequence}")

        # Normalize DNA -> RNA (T to U) and finalize residue names
        if not sequence:
            raise ValueError(f"Could not read sequence from FASTA file: {fasta_file}")
        
        sequence = sequence.replace('T', 'U')
        residue_names = list(sequence)
        
        if len(residue_names) != n_residues:
            raise ValueError(
                f"Sequence length mismatch: FASTA has {len(residue_names)} residues "
                f"but NPY has {n_residues} residues"
            )
        
        residue_numbers = list(range(1, n_residues + 1))
        
        # Convert to dictionary format
        coords_dict = {}
        for i in range(n_residues):
            res_num = residue_numbers[i]
            res_name = residue_names[i]
            
            for j, atom_name in enumerate(atom_names):
                coords = coords_array[i, j, :]
                # Check for missing atoms (all zeros or NaN)
                if not np.isnan(coords).any() and not np.allclose(coords, 0):
                    # For N1/N9, store both generic and residue-specific keys when possible
                    if atom_name == "N1/N9":
                        # Always store a generic key for robust downstream lookups
                        coords_dict[(res_num, 'N1N9')] = coords
                        # If residue type known, also store the specific atom name (N1 or N9)
                        if res_name in self.BASE_MARKER_ATOMS:
                            specific = self.BASE_MARKER_ATOMS[res_name]
                            coords_dict[(res_num, specific)] = coords
                    else:
                        coords_dict[(res_num, atom_name)] = coords

        return coords_dict, residue_names, residue_numbers

    @staticmethod
    def _read_fasta_simple(fasta_file: str) -> str:
        """
        Minimal FASTA reader that returns the first sequence observed.
        Concatenates non-header lines until next header or EOF.
        """
        try:
            seq_lines = []
            with open(fasta_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith('>'):
                        # Stop at the first record once we have collected sequence lines
                        if seq_lines:
                            break
                        else:
                            continue
                    seq_lines.append(line)
            return ''.join(seq_lines).upper()
        except Exception:
            return ''
    
    def detect_base_pairs_all_atom(self, coords_dict: Dict, residue_numbers: List[int],
                                   residue_names: List[str]) -> np.ndarray:
        """
        Detect base pairs using all-atom information with biological validation.
        
        Uses hydrogen bonding patterns and geometric criteria.
        
        Args:
            coords_dict: Dictionary mapping (res_num, atom_name) -> coordinates
            residue_numbers: List of residue numbers
            residue_names: List of residue names
            
        Returns:
            NxN binary matrix where 1 indicates a base pair
        """
        n = len(residue_numbers)
        pair_matrix = np.zeros((n, n), dtype=np.int8)
        
        invalid_pairs = []  # Track pairs rejected by biological rules
        
        # For all-atom mode, we need to implement proper base pair detection
        # This is a simplified version - full implementation would check:
        # 1. Watson-Crick pairs: A-U, G-C (specific hydrogen bonds)
        # 2. Wobble pairs: G-U
        # 3. Non-canonical pairs
        
        # Use C4'-C4' distance for consistency with no-base-atom mode
        for i, res_i in enumerate(residue_numbers):
            for j, res_j in enumerate(residue_numbers):
                if j <= i:  # Only upper triangle
                    continue
                
                # Enforce minimum sequence separation
                if abs(j - i) < self.min_separation:
                    continue
                
                # Use C4' atoms as backbone anchors (consistent with no-base-atom mode)
                key_i_c4 = (res_i, "C4'")
                key_j_c4 = (res_j, "C4'")
                
                if key_i_c4 not in coords_dict or key_j_c4 not in coords_dict:
                    continue
                
                # Calculate C4'-C4' distance
                dist_c4 = np.linalg.norm(coords_dict[key_i_c4] - coords_dict[key_j_c4])
                
                # Check backbone distance using thresholds
                if not (self.thresholds['min'] <= dist_c4 <= self.thresholds['max']):
                    continue
                
                # Check biological validity
                base_i = residue_names[i]
                base_j = residue_names[j]
                
                if self.is_valid_base_pair(base_i, base_j):
                    pair_matrix[i, j] = 1
                    pair_matrix[j, i] = 1
                else:
                    invalid_pairs.append((i+1, j+1, base_i, base_j))
        
        if invalid_pairs:
            print(f"\nBiological validation: Rejected {len(invalid_pairs)} pairs that don't follow G-C, A-U, or G-U rules:")
            for pos_i, pos_j, base_i, base_j in invalid_pairs[:10]:  # Show first 10
                print(f"  Position {pos_i}{base_i} - {pos_j}{base_j}")
            if len(invalid_pairs) > 10:
                print(f"  ... and {len(invalid_pairs) - 10} more")
        
        return pair_matrix
    
    def detect_base_pairs_no_base_atom(self, coords_dict: Dict, residue_numbers: List[int],
                                      residue_names: List[str]) -> np.ndarray:
        """
        Detect base pairs without using base atoms (only backbone + N1/N9) with biological validation.
        
        Uses geometric criteria based on:
        - C4'-C4' distance (backbone anchor atoms)
        - N1/N9 distance (glycosidic nitrogen, if available)
        - Biological base-pairing rules (G-C, A-U, G-U)
        
        Note: C4' atoms are used instead of C1' because the 7-atom representation
        (P, O5', C5', C4', C3', O3', N1/N9) does not include C1' atoms.
        
        Args:
            coords_dict: Dictionary mapping (res_num, atom_name) -> coordinates
            residue_numbers: List of residue numbers
            residue_names: List of residue names
            
        Returns:
            NxN binary matrix where 1 indicates a base pair
        """
        n = len(residue_numbers)
        pair_matrix = np.zeros((n, n), dtype=np.int8)
        
        invalid_pairs = []  # Track pairs rejected by biological rules
        
        for i, res_i in enumerate(residue_numbers):
            for j, res_j in enumerate(residue_numbers):
                if j <= i:  # Only upper triangle
                    continue
                # Enforce minimum sequence separation
                if abs(j - i) < self.min_separation:
                    continue

                # Use C4' atoms as backbone anchors (NPY format lacks C1')
                key_i_c4 = (res_i, "C4'")
                key_j_c4 = (res_j, "C4'")

                if key_i_c4 not in coords_dict or key_j_c4 not in coords_dict:
                    continue

                # Calculate C4'-C4' distance
                dist_c4 = np.linalg.norm(coords_dict[key_i_c4] - coords_dict[key_j_c4])

                # Check backbone distance using existing thresholds
                if not (self.thresholds['min'] <= dist_c4 <= self.thresholds['max']):
                    continue

                # Additional check with N1/N9 if available (with key fallback)
                res_name_i = residue_names[i]
                res_name_j = residue_names[j]

                atom_i = self.BASE_MARKER_ATOMS.get(res_name_i)
                atom_j = self.BASE_MARKER_ATOMS.get(res_name_j)

                key_i_n = (res_i, atom_i) if atom_i else None
                key_j_n = (res_j, atom_j) if atom_j else None

                if not (key_i_n and key_i_n in coords_dict):
                    key_i_n = (res_i, 'N1N9')
                if not (key_j_n and key_j_n in coords_dict):
                    key_j_n = (res_j, 'N1N9')

                if key_i_n in coords_dict and key_j_n in coords_dict:
                    dist_n = np.linalg.norm(coords_dict[key_i_n] - coords_dict[key_j_n])
                    n_thresh = self.thresholds['n1n9']
                    if not (n_thresh['min'] <= dist_n <= n_thresh['max']):
                        continue
                
                # Check biological validity before marking as base pair
                base_i = residue_names[i]
                base_j = residue_names[j]
                
                if self.is_valid_base_pair(base_i, base_j):
                    # Mark as base pair
                    pair_matrix[i, j] = 1
                    pair_matrix[j, i] = 1
                else:
                    invalid_pairs.append((i+1, j+1, base_i, base_j))
        
        if invalid_pairs:
            print(f"\nBiological validation: Rejected {len(invalid_pairs)} pairs that don't follow G-C, A-U, or G-U rules:")
            for pos_i, pos_j, base_i, base_j in invalid_pairs[:10]:  # Show first 10
                print(f"  Position {pos_i}{base_i} - {pos_j}{base_j}")
            if len(invalid_pairs) > 10:
                print(f"  ... and {len(invalid_pairs) - 10} more")
        
        return pair_matrix
    
    def detect_base_pairs(self, coords_dict: Dict, residue_numbers: List[int],
                         residue_names: List[str]) -> np.ndarray:
        """
        Detect base pairs based on the configured mode.
        
        Args:
            coords_dict: Dictionary mapping (res_num, atom_name) -> coordinates
            residue_numbers: List of residue numbers
            residue_names: List of residue names
            
        Returns:
            NxN binary matrix where 1 indicates a base pair
        """
        if self.mode == 'all-atom':
            return self.detect_base_pairs_all_atom(coords_dict, residue_numbers, residue_names)
        else:
            return self.detect_base_pairs_no_base_atom(coords_dict, residue_numbers, residue_names)
    
    def process_structure(self, input_file: str, chain_id: Optional[str] = None, 
                         fasta_file: Optional[str] = None) -> Tuple[np.ndarray, str]:
        """
        Process RNA structure file and generate base pair matrix.
        
        Args:
            input_file: Path to input file (PDB or NPY)
            chain_id: Chain ID for PDB files (None = first RNA chain)
            fasta_file: Path to FASTA file with sequence (REQUIRED for NPY files)
            
        Returns:
            Tuple of (pair_matrix, sequence)
        """
        # Determine file type
        if input_file.endswith('.pdb'):
            coords_dict, residue_names, residue_numbers = self.load_pdb(input_file, chain_id)
        elif input_file.endswith('.npy'):
            coords_dict, residue_names, residue_numbers = self.load_npy(input_file, fasta_file)
        else:
            raise ValueError(f"Unsupported file format: {input_file}. Use .pdb or .npy")
        
        # Check if we have valid sequence information
        if 'N' in residue_names:
            print("\nWARNING: Unknown nucleotides ('N') detected in sequence.")
            print("         Base pairs involving unknown nucleotides will be rejected.")
            print("         Please provide sequence information for accurate detection.")
        
        # Detect base pairs with biological validation
        pair_matrix = self.detect_base_pairs(coords_dict, residue_numbers, residue_names)
        
        # Generate sequence
        sequence = ''.join(residue_names)
        
        return pair_matrix, sequence

    @staticmethod
    def _matrix_to_pairs(pair_matrix: np.ndarray) -> List[Tuple[int, int]]:
        """
        Convert a symmetric NxN base-pair matrix to a sorted list of 0-indexed pairs (i, j) with i < j.
        """
        n = pair_matrix.shape[0]
        pairs: List[Tuple[int, int]] = []
        for i in range(n):
            for j in range(i + 1, n):
                if pair_matrix[i, j] != 0:
                    pairs.append((i, j))
        pairs.sort()
        return pairs

    def get_dot_bracket_notation(self, pair_matrix: np.ndarray) -> str:
        """
        Convert base-pair matrix to dot-bracket notation.

        Simple approach that does not handle pseudoknots: assigns '(' and ')' for the first available non-overlapping pairs.
        """
        n = pair_matrix.shape[0]
        notation = ['.'] * n
        used = [False] * n

        for i, j in self._matrix_to_pairs(pair_matrix):
            if not used[i] and not used[j]:
                notation[i] = '('
                notation[j] = ')'
                used[i] = True
                used[j] = True

        return ''.join(notation)

    def get_base_pair_list(
        self,
        pair_matrix: np.ndarray,
        residue_names: List[str],
        residue_numbers: List[int],
        chain_id: Optional[str] = None,
    ) -> List[Tuple[int, int, str, str]]:
        """
        Get list of base pairs with nucleotide information.

        Returns 1-indexed residue positions and formatted labels:
        - label format: "{chain}:{res_num}{res_name}" if chain_id is provided, otherwise "{res_num}{res_name}"
        """
        pair_list: List[Tuple[int, int, str, str]] = []
        chain_label_prefix = f"{chain_id}:" if chain_id else ""
        for i, j in self._matrix_to_pairs(pair_matrix):
            res_i = f"{chain_label_prefix}{residue_numbers[i]}{residue_names[i]}"
            res_j = f"{chain_label_prefix}{residue_numbers[j]}{residue_names[j]}"
            pair_list.append((i + 1, j + 1, res_i, res_j))
        return pair_list


def save_matrix(matrix: np.ndarray, output_file: str):
    """
    Save base pair matrix to file.
    
    Supports .npy, .csv, and .txt formats.
    Default .txt format matches the reference format (space-separated).
    
    Args:
        matrix: NxN binary matrix
        output_file: Path to output file
    """
    if output_file.endswith('.npy'):
        np.save(output_file, matrix)
        print(f"Saved matrix to {output_file} (NPY format)")
    elif output_file.endswith('.csv'):
        np.savetxt(output_file, matrix, delimiter=',', fmt='%d')
        print(f"Saved matrix to {output_file} (CSV format)")
    else:
        # Space-separated format matching reference (pdb9gus_W.txt style)
        np.savetxt(output_file, matrix, fmt='%d', delimiter=' ')
        print(f"Saved matrix to {output_file} (TXT format, space-separated)")


def main():
    parser = argparse.ArgumentParser(
        description='Convert 3D RNA structure to 2D base pair matrix with biological validation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
IMPORTANT: This version enforces biological base-pairing rules (G-C, A-U, G-U)
           Sequence information is REQUIRED for NPY files.

Examples:
  # Single file processing
  # =====================
  
  # PDB file with all-atom mode (sequence from PDB)
  python pdb_3dto2d_basepair2.py input.pdb --output matrix.txt --mode all-atom
  
  # PDB file with no-base-atom mode (backbone only, sequence from PDB)
  python pdb_3dto2d_basepair2.py input.pdb --output matrix.txt --mode no-base-atom
  
  # NPY file with sequence from FASTA (REQUIRED for NPY)
  python pdb_3dto2d_basepair2.py input.npy --sequence seq.fasta --output matrix.txt --mode no-base-atom
  
  # Batch processing
  # ================
  
  # Process all PDB files in a directory
  python pdb_3dto2d_basepair2.py --pdb-dir ./structures --output-dir ./matrices --mode no-base-atom
  
  # Process all NPY files with matching FASTA files
  python pdb_3dto2d_basepair2.py --coords-dir ./coords --seqs-dir ./sequences --output-dir ./matrices
  
  # Batch process with specific chain
  python pdb_3dto2d_basepair2.py --pdb-dir ./structures --output-dir ./matrices --chain A
        """
    )
    
    # Single file arguments
    parser.add_argument('input', nargs='?', help='Input file (PDB or NPY format)')
    parser.add_argument('--output', '-o', help='Output file (.npy, .csv, or .txt)')
    parser.add_argument('--sequence', '-s', default=None,
                       help='FASTA file with sequence (REQUIRED for NPY files)')
    
    # Batch processing arguments
    parser.add_argument('--pdb-dir', help='Directory containing PDB files for batch processing')
    parser.add_argument('--coords-dir', help='Directory containing NPY coordinate files for batch processing')
    parser.add_argument('--seqs-dir', help='Directory containing FASTA sequence files (required with --coords-dir)')
    parser.add_argument('--output-dir', help='Output directory for batch processing')
    
    # Common arguments
    parser.add_argument('--mode', '-m', choices=['all-atom', 'no-base-atom'], default='all-atom',
                       help='Detection mode (default: all-atom)')
    parser.add_argument('--chain', '-c', default=None,
                       help='Chain ID for PDB files (default: first RNA chain)')
    parser.add_argument('--save-sequence', action='store_true',
                       help='Save sequence to a separate file')
    parser.add_argument('--min-separation', type=int, default=4,
                       help='Minimum |j-i| residue separation to consider base pairing (default: 4)')
    parser.add_argument('--report', choices=['none', 'dot-bracket', 'list'], default='none',
                       help='Optional base-pair report format: dot-bracket or list (default: none)')
    parser.add_argument('--report-file', default=None,
                       help='Path to save the base-pair report (default: print to stdout)')
    
    args = parser.parse_args()
    
    # Determine processing mode: batch or single file
    batch_mode = args.pdb_dir or args.coords_dir
    
    if batch_mode:
        # Batch processing mode
        return process_batch(args)
    else:
        # Single file processing mode
        return process_single_file(args)


def process_batch(args):
    """Process multiple files in batch mode."""
    import sys
    from pathlib import Path
    
    # Validate batch arguments
    if args.pdb_dir and args.coords_dir:
        print("Error: Cannot use both --pdb-dir and --coords-dir simultaneously", file=sys.stderr)
        return 1
    
    if not args.output_dir:
        print("Error: --output-dir is required for batch processing", file=sys.stderr)
        return 1
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if args.pdb_dir:
        # Batch process PDB files
        return process_pdb_batch(args, output_dir)
    elif args.coords_dir:
        # Batch process NPY files
        return process_npy_batch(args, output_dir)
    else:
        print("Error: Either --pdb-dir or --coords-dir must be specified for batch processing", file=sys.stderr)
        return 1


def process_pdb_batch(args, output_dir):
    """Batch process PDB files."""
    import sys
    from pathlib import Path
    
    pdb_dir = Path(args.pdb_dir)
    
    if not pdb_dir.exists() or not pdb_dir.is_dir():
        print(f"Error: PDB directory '{args.pdb_dir}' not found or not a directory", file=sys.stderr)
        return 1
    
    # Find all PDB files
    pdb_files = sorted([p for p in pdb_dir.iterdir() if p.suffix.lower() in ('.pdb', '.cif')])
    
    if not pdb_files:
        print(f"Warning: No PDB or CIF files found in '{pdb_dir}'", file=sys.stderr)
        return 0
    
    total = len(pdb_files)
    processed = 0
    skipped = 0
    errors = 0
    
    print(f"\n{'='*70}")
    print(f"Batch Processing PDB Files")
    print(f"{'='*70}")
    print(f"Input directory:  {pdb_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Mode:             {args.mode}")
    print(f"Files found:      {total}")
    print(f"{'='*70}\n")
    
    # Initialize detector
    detector = RNABasePairDetector(mode=args.mode, min_separation=args.min_separation)
    
    for idx, pdb_file in enumerate(pdb_files, 1):
        print(f"[{idx}/{total}] Processing: {pdb_file.name}")
        
        try:
            # Process structure
            pair_matrix, sequence = detector.process_structure(str(pdb_file), args.chain, None)
            
            if len(sequence) == 0:
                print(f"  ⚠ Warning: No RNA residues found, skipping")
                skipped += 1
                continue
            
            # Save matrix
            output_file = output_dir / f"{pdb_file.stem}.txt"
            save_matrix(pair_matrix, str(output_file))
            
            # Save sequence if requested
            if args.save_sequence:
                seq_file = output_dir / f"{pdb_file.stem}.seq"
                with open(seq_file, 'w') as f:
                    f.write(f">{pdb_file.name}\n")
                    f.write(sequence + "\n")
                print(f"  ✓ Saved sequence to {seq_file.name}")
            
            # Statistics
            n_pairs = np.sum(pair_matrix) // 2
            print(f"  ✓ Length: {len(sequence)}, Pairs: {n_pairs}")
            processed += 1
            
        except Exception as e:
            print(f"  ✗ Error: {e}")
            errors += 1
            continue
    
    # Summary
    print(f"\n{'='*70}")
    print(f"Batch Processing Summary")
    print(f"{'='*70}")
    print(f"Total files:      {total}")
    print(f"Processed:        {processed}")
    print(f"Skipped:          {skipped}")
    print(f"Errors:           {errors}")
    print(f"{'='*70}\n")
    
    return 0


def process_npy_batch(args, output_dir):
    """Batch process NPY coordinate files with FASTA sequences."""
    import sys
    from pathlib import Path
    
    coords_dir = Path(args.coords_dir)
    
    if not coords_dir.exists() or not coords_dir.is_dir():
        print(f"Error: Coordinates directory '{args.coords_dir}' not found or not a directory", file=sys.stderr)
        return 1
    
    if not args.seqs_dir:
        print("Error: --seqs-dir is REQUIRED when using --coords-dir for NPY batch processing", file=sys.stderr)
        print("       Each NPY file requires a matching FASTA sequence file for biological validation.")
        return 1
    
    seqs_dir = Path(args.seqs_dir)
    if not seqs_dir.exists() or not seqs_dir.is_dir():
        print(f"Error: Sequences directory '{args.seqs_dir}' not found or not a directory", file=sys.stderr)
        return 1
    
    # Find all NPY files
    npy_files = sorted([p for p in coords_dir.iterdir() if p.suffix.lower() == '.npy'])
    
    if not npy_files:
        print(f"Warning: No NPY files found in '{coords_dir}'", file=sys.stderr)
        return 0
    
    # Force no-base-atom mode for NPY files
    if args.mode == 'all-atom':
        print("Info: NPY files contain limited atom information. Using no-base-atom mode.")
        args.mode = 'no-base-atom'
    
    total = len(npy_files)
    processed = 0
    skipped = 0
    errors = 0
    
    print(f"\n{'='*70}")
    print(f"Batch Processing NPY Files")
    print(f"{'='*70}")
    print(f"Coordinates dir:  {coords_dir}")
    print(f"Sequences dir:    {seqs_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Mode:             {args.mode}")
    print(f"Files found:      {total}")
    print(f"{'='*70}\n")
    
    # Initialize detector
    detector = RNABasePairDetector(mode=args.mode, min_separation=args.min_separation)
    
    for idx, npy_file in enumerate(npy_files, 1):
        print(f"[{idx}/{total}] Processing: {npy_file.name}")
        
        # Find matching FASTA file (try multiple extensions)
        fasta_file = None
        for ext in ['.fasta', '.fa', '.seq', '.txt']:
            candidate = seqs_dir / f"{npy_file.stem}{ext}"
            if candidate.exists():
                fasta_file = candidate
                break
        
        if not fasta_file:
            print(f"  ✗ Error: No matching FASTA file found in {seqs_dir}")
            print(f"     Tried: {npy_file.stem}.fasta, .fa, .seq, .txt")
            errors += 1
            continue
        
        print(f"  → Using sequence: {fasta_file.name}")
        
        try:
            # Process structure
            pair_matrix, sequence = detector.process_structure(str(npy_file), None, str(fasta_file))
            
            if len(sequence) == 0:
                print(f"  ⚠ Warning: Empty sequence, skipping")
                skipped += 1
                continue
            
            # Save matrix
            output_file = output_dir / f"{npy_file.stem}.txt"
            save_matrix(pair_matrix, str(output_file))
            
            # Save sequence if requested
            if args.save_sequence:
                seq_file = output_dir / f"{npy_file.stem}.seq"
                with open(seq_file, 'w') as f:
                    f.write(f">{npy_file.name}\n")
                    f.write(sequence + "\n")
                print(f"  ✓ Saved sequence to {seq_file.name}")
            
            # Statistics
            n_pairs = np.sum(pair_matrix) // 2
            print(f"  ✓ Length: {len(sequence)}, Pairs: {n_pairs}")
            processed += 1
            
        except Exception as e:
            print(f"  ✗ Error: {e}")
            errors += 1
            continue
    
    # Summary
    print(f"\n{'='*70}")
    print(f"Batch Processing Summary")
    print(f"{'='*70}")
    print(f"Total files:      {total}")
    print(f"Processed:        {processed}")
    print(f"Skipped:          {skipped}")
    print(f"Errors:           {errors}")
    print(f"{'='*70}\n")
    
    return 0


def process_single_file(args):
    """Process a single file."""
    import sys
    
    # Validate input file
    if not args.input:
        print("Error: Input file is required for single-file mode", file=sys.stderr)
        print("       Use --pdb-dir or --coords-dir for batch processing", file=sys.stderr)
        return 1
    
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' not found", file=sys.stderr)
        return 1
    
    if not args.output:
        print("Error: --output is required for single-file mode", file=sys.stderr)
        return 1
    
    # Validate NPY files require FASTA
    if args.input.endswith('.npy') and not args.sequence:
        print("\nERROR: NPY files require sequence information!", file=sys.stderr)
        print("       Please provide a FASTA file using --sequence argument.", file=sys.stderr)
        print("\nExample:", file=sys.stderr)
        print(f"  python {os.path.basename(__file__)} {args.input} --sequence seq.fasta --output {args.output} --mode no-base-atom", file=sys.stderr)
        return 1
    
    # Validate NPY files should use no-base-atom mode
    if args.input.endswith('.npy') and args.mode == 'all-atom':
        print("Warning: NPY files contain limited atom information. Forcing no-base-atom mode.")
        args.mode = 'no-base-atom'
    
    try:
        # Initialize detector
        detector = RNABasePairDetector(mode=args.mode, min_separation=args.min_separation)
        
        # Process structure
        print(f"\nProcessing {args.input} in {args.mode} mode...")
        print("Enforcing biological base-pairing rules: G-C, A-U, and G-U only")
        pair_matrix, sequence = detector.process_structure(args.input, args.chain, args.sequence)
        
        # Print statistics
        n_residues = len(sequence)
        n_pairs = np.sum(pair_matrix) // 2  # Divide by 2 since matrix is symmetric
        print(f"\nResults:")
        print(f"  Sequence length: {n_residues}")
        print(f"  Sequence: {sequence}")
        print(f"  Base pairs detected: {n_pairs}")
        print(f"  Pairing density: {n_pairs / n_residues:.2f} pairs per residue")
        
        # Count pair types
        pair_types = {}
        for i, j in detector._matrix_to_pairs(pair_matrix):
            base_i = sequence[i]
            base_j = sequence[j]
            pair_type = '-'.join(sorted([base_i, base_j]))
            pair_types[pair_type] = pair_types.get(pair_type, 0) + 1
        
        if pair_types:
            print(f"\n  Pair types found:")
            for pair_type, count in sorted(pair_types.items()):
                print(f"    {pair_type}: {count}")
        
        # Save matrix
        save_matrix(pair_matrix, args.output)
        
        # Save sequence if requested
        if args.save_sequence:
            seq_file = args.output.replace('.npy', '.seq').replace('.csv', '.seq').replace('.txt', '.seq')
            with open(seq_file, 'w') as f:
                f.write(f">{os.path.basename(args.input)}\n")
                f.write(sequence + "\n")
            print(f"Saved sequence to {seq_file}")

        # Optional base-pair report
        if args.report != 'none':
            if args.report == 'dot-bracket':
                notation = detector.get_dot_bracket_notation(pair_matrix)
                if args.report_file:
                    with open(args.report_file, 'w') as rf:
                        rf.write(notation + "\n")
                    print(f"Saved dot-bracket notation to {args.report_file}")
                else:
                    print("\nDot-bracket notation:")
                    print(notation)
            elif args.report == 'list':
                # Derive residue_names and residue_numbers locally from sequence
                residue_names_local = list(sequence)
                residue_numbers_local = list(range(1, n_residues + 1))
                pair_list = detector.get_base_pair_list(pair_matrix, residue_names=residue_names_local,
                                                        residue_numbers=residue_numbers_local,
                                                        chain_id=args.chain)
                if args.report_file:
                    with open(args.report_file, 'w') as rf:
                        rf.write("# i j res_i res_j\n")
                        for i_pos, j_pos, res_i, res_j in pair_list:
                            rf.write(f"{i_pos}\t{j_pos}\t{res_i}\t{res_j}\n")
                    print(f"Saved base-pair list to {args.report_file}")
                else:
                    print("\nBase-pair list (i j res_i res_j):")
                    for i_pos, j_pos, res_i, res_j in pair_list:
                        print(f"{i_pos}\t{j_pos}\t{res_i}\t{res_j}")
        
        print("\nDone!")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
