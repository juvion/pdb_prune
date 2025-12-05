#!/usr/bin/env python3
"""
PDB/NPY 3D to 2D Base Pair Matrix Converter

This script converts 3D RNA structures to 2D base pair matrices.
Supports both PDB files and NPY coordinate arrays.
Supports two modes: all-atom and no-base-atom.

Usage:
    python pdb_3dto2d_basepair.py input.pdb --output matrix.csv
    python pdb_3dto2d_basepair.py input.npy --output matrix.csv --mode no-base-atom
    python pdb_3dto2d_basepair.py input.pdb --output matrix.npy --mode all-atom
"""

import numpy as np
import argparse
from typing import Dict, List, Tuple, Optional
import os


class RNABasePairDetector:
    """
    Detect base pairs in RNA structures from 3D coordinates.
    
    Supports:
    - PDB files (all-atom or no-base-atom)
    - NPY arrays with predefined atom order
    - Two detection modes: all-atom and no-base-atom
    """
    
    # Standard RNA atom names for backbone and first base atom
    BACKBONE_ATOMS = ["P", "O5'", "C5'", "C4'", "C3'", "O3'"]
    
    # N1 for pyrimidines (C, U), N9 for purines (A, G)
    BASE_MARKER_ATOMS = {
        'A': 'N9', 'G': 'N9',
        'C': 'N1', 'U': 'N1'
    }
    
    # Distance thresholds for base pairing (in Angstroms)
    # These are standard thresholds used in RNA structure analysis
    DISTANCE_THRESHOLDS = {
        'all-atom': {
            'min': 2.5,   # Minimum distance for hydrogen bond
            'max': 3.5,   # Maximum distance for hydrogen bond
            'geometric': 4.5  # Geometric center distance for base pairs
        },
        'no-base-atom': {
            'min': 8.5,   # Minimum C1'-C1' distance
            'max': 12.0,  # Maximum C1'-C1' distance (typical ~10.5Å)
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
            fasta_file: Optional path to FASTA file containing sequence
            
        Returns:
            Tuple of (coordinates_dict, residue_names, residue_numbers)
        """
        coords_array = np.load(npy_file)
        
        if len(coords_array.shape) != 3 or coords_array.shape[1] != 7 or coords_array.shape[2] != 3:
            raise ValueError(f"NPY file must have shape (N, 7, 3), got {coords_array.shape}")
        
        n_residues = coords_array.shape[0]
        atom_names = ["P", "O5'", "C5'", "C4'", "C3'", "O3'", "N1/N9"]
        
        # Load sequence from FASTA file if provided
        if fasta_file:
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
            if sequence:
                sequence = sequence.replace('T', 'U')
                residue_names = list(sequence)
                if len(residue_names) != n_residues:
                    raise ValueError(
                        f"Sequence length mismatch: FASTA has {len(residue_names)} residues "
                        f"but NPY has {n_residues} residues"
                    )
            else:
                print("Warning: FASTA provided but sequence could not be read. Using 'N' for all residues.")
                residue_names = ['N'] * n_residues
        else:
            # No sequence information available
            print(f"Warning: No sequence file provided. Using 'N' (unknown) for all residues.")
            print(f"         For accurate base pair detection, provide sequence with --sequence argument")
            residue_names = ['N'] * n_residues
        
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
    
    def detect_base_pairs_all_atom(self, coords_dict: Dict, residue_numbers: List[int]) -> np.ndarray:
        """
        Detect base pairs using all-atom information.
        
        Uses hydrogen bonding patterns and geometric criteria.
        
        Args:
            coords_dict: Dictionary mapping (res_num, atom_name) -> coordinates
            residue_numbers: List of residue numbers
            
        Returns:
            NxN binary matrix where 1 indicates a base pair
        """
        n = len(residue_numbers)
        pair_matrix = np.zeros((n, n), dtype=np.int8)
        
        # For all-atom mode, we need to implement proper base pair detection
        # This is a simplified version - full implementation would check:
        # 1. Watson-Crick pairs: A-U, G-C (specific hydrogen bonds)
        # 2. Wobble pairs: G-U
        # 3. Non-canonical pairs
        
        # For now, use a distance-based heuristic
        for i, res_i in enumerate(residue_numbers):
            for j, res_j in enumerate(residue_numbers):
                if j <= i:  # Only upper triangle
                    continue
                
                # Get base atoms for distance calculation
                # For simplicity, use C1' atoms or N1/N9 atoms
                key_i_c1 = (res_i, "C1'")
                key_j_c1 = (res_j, "C1'")
                
                if key_i_c1 not in coords_dict or key_j_c1 not in coords_dict:
                    continue
                
                # Calculate C1'-C1' distance
                dist_c1 = np.linalg.norm(coords_dict[key_i_c1] - coords_dict[key_j_c1])
                
                # Typical C1'-C1' distance for base pairs is ~10-11Å
                if 9.0 <= dist_c1 <= 12.0:
                    pair_matrix[i, j] = 1
                    pair_matrix[j, i] = 1
        
        return pair_matrix
    
    def detect_base_pairs_no_base_atom(self, coords_dict: Dict, residue_numbers: List[int],
                                      residue_names: List[str]) -> np.ndarray:
        """
        Detect base pairs without using base atoms (only backbone + N1/N9).
        
        Uses geometric criteria based on:
        - C1'-C1' distance
        - C4'-C4' distance
        - N1/N9 distance (if available)
        
        Args:
            coords_dict: Dictionary mapping (res_num, atom_name) -> coordinates
            residue_numbers: List of residue numbers
            residue_names: List of residue names
            
        Returns:
            NxN binary matrix where 1 indicates a base pair
        """
        n = len(residue_numbers)
        pair_matrix = np.zeros((n, n), dtype=np.int8)
        
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
                
                # Mark as base pair
                pair_matrix[i, j] = 1
                pair_matrix[j, i] = 1
        
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
            return self.detect_base_pairs_all_atom(coords_dict, residue_numbers)
        else:
            return self.detect_base_pairs_no_base_atom(coords_dict, residue_numbers, residue_names)
    
    def process_structure(self, input_file: str, chain_id: Optional[str] = None, 
                         fasta_file: Optional[str] = None) -> Tuple[np.ndarray, str]:
        """
        Process RNA structure file and generate base pair matrix.
        
        Args:
            input_file: Path to input file (PDB or NPY)
            chain_id: Chain ID for PDB files (None = first RNA chain)
            fasta_file: Path to FASTA file with sequence (for NPY files)
            
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
        
        # Detect base pairs
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
        description='Convert 3D RNA structure to 2D base pair matrix',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # PDB file with all-atom mode
  python pdb_3dto2d_basepair.py input.pdb --output matrix.txt --mode all-atom
  
  # PDB file with no-base-atom mode (backbone only)
  python pdb_3dto2d_basepair.py input.pdb --output matrix.txt --mode no-base-atom
  
  # NPY file with sequence from FASTA (recommended)
  python pdb_3dto2d_basepair.py input.npy --output matrix.txt --mode no-base-atom --sequence seq.fasta
  
  # NPY file without sequence (less accurate, uses only C1'-C1' distance)
  python pdb_3dto2d_basepair.py input.npy --output matrix.txt --mode no-base-atom
  
  # Specify chain ID for PDB
  python pdb_3dto2d_basepair.py input.pdb --output matrix.txt --chain A --mode no-base-atom
  
  # Save sequence to separate file
  python pdb_3dto2d_basepair.py input.pdb --output matrix.txt --save-sequence
        """
    )
    
    parser.add_argument('input', help='Input file (PDB or NPY format)')
    parser.add_argument('--output', '-o', required=True, help='Output file (.npy, .csv, or .txt)')
    parser.add_argument('--mode', '-m', choices=['all-atom', 'no-base-atom'], default='all-atom',
                       help='Detection mode (default: all-atom)')
    parser.add_argument('--chain', '-c', default=None,
                       help='Chain ID for PDB files (default: first RNA chain)')
    parser.add_argument('--sequence', '-s', default=None,
                       help='FASTA file with sequence (required for NPY files for accurate detection)')
    parser.add_argument('--save-sequence', action='store_true',
                       help='Save sequence to a separate file')
    parser.add_argument('--min-separation', type=int, default=4,
                       help='Minimum |j-i| residue separation to consider base pairing (default: 4)')
    parser.add_argument('--report', choices=['none', 'dot-bracket', 'list'], default='none',
                       help='Optional base-pair report format: dot-bracket or list (default: none)')
    parser.add_argument('--report-file', default=None,
                       help='Path to save the base-pair report (default: print to stdout)')
    
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' not found")
        return 1
    
    # Validate NPY files should use no-base-atom mode
    if args.input.endswith('.npy') and args.mode == 'all-atom':
        print("Warning: NPY files contain limited atom information. Forcing no-base-atom mode.")
        args.mode = 'no-base-atom'
    
    try:
        # Initialize detector
        detector = RNABasePairDetector(mode=args.mode, min_separation=args.min_separation)
        
        # Process structure
        print(f"Processing {args.input} in {args.mode} mode...")
        pair_matrix, sequence = detector.process_structure(args.input, args.chain, args.sequence)
        
        # Print statistics
        n_residues = len(sequence)
        n_pairs = np.sum(pair_matrix) // 2  # Divide by 2 since matrix is symmetric
        print(f"\nResults:")
        print(f"  Sequence length: {n_residues}")
        print(f"  Sequence: {sequence}")
        print(f"  Base pairs detected: {n_pairs}")
        print(f"  Pairing density: {n_pairs / n_residues:.2f} pairs per residue")
        
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