#!/usr/bin/env python3
"""
RNA Base Pairing Matrix Tool v3.0
Extracts base pairing information from RNA PDB structures
Adds N1/N9 distance criterion alongside existing C1'-C1' approach
"""

import numpy as np
import argparse
import sys
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass
import warnings
from pathlib import Path

@dataclass
class Atom:
    """Represents an atom from a PDB file"""
    name: str
    x: float
    y: float
    z: float
    residue_name: str
    chain_id: str
    residue_num: int
    
    @property
    def coords(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z])

@dataclass
class Nucleotide:
    """Represents a nucleotide with its atoms"""
    residue_name: str
    chain_id: str
    residue_num: int
    atoms: Dict[str, Atom]
    
    def get_atom(self, atom_name: str) -> Optional[Atom]:
        """Get atom by name, handling common variations"""
        if atom_name in self.atoms:
            return self.atoms[atom_name]
        # Handle primed atoms
        if "'" in atom_name:
            alt_name = atom_name.replace("'", "*")
            if alt_name in self.atoms:
                return self.atoms[alt_name]
        return None
    
    def get_base_atoms(self) -> List[Atom]:
        """Get all base atoms (excluding sugar-phosphate backbone)"""
        base_atoms = []
        backbone_atoms = {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"}
        for atom_name, atom in self.atoms.items():
            if atom_name not in backbone_atoms and not atom_name.startswith("H"):
                base_atoms.append(atom)
        return base_atoms
    
    def get_glycosidic_atom(self) -> Optional[Atom]:
        """Get N1 for pyrimidines (C, U) or N9 for purines (A, G)"""
        if self.residue_name in ['C', 'U']:  # Pyrimidines
            return self.get_atom('N1')
        elif self.residue_name in ['A', 'G']:  # Purines
            return self.get_atom('N9')
        return None

class PDBParser:
    """Parse PDB/CIF files and extract nucleotide information"""
    
    RNA_RESIDUES = {"A", "U", "G", "C", "ADE", "URA", "GUA", "CYT", 
                    "rA", "rU", "rG", "rC", "DA", "DT", "DG", "DC"}
    
    @staticmethod
    def parse(filename: str) -> List[Nucleotide]:
        """Parse PDB or CIF file and return list of nucleotides"""
        # Check if it's a CIF file
        with open(filename, 'r') as f:
            first_line = f.readline()
            if first_line.startswith('#') or 'data_' in first_line:
                return PDBParser.parse_cif(filename)
            else:
                return PDBParser.parse_pdb(filename)
    
    @staticmethod
    def parse_cif(filename: str) -> List[Nucleotide]:
        """Parse CIF/mmCIF file format"""
        nucleotides_dict = {}
        
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        in_atom_site = False
        for line in lines:
            line = line.strip()
            
            if line.startswith('ATOM') or line.startswith('HETATM'):
                parts = line.split()
                if len(parts) >= 16:
                    try:
                        atom_name = parts[3]
                        residue_name = parts[5]
                        chain_id = parts[14] if parts[14] != '?' else 'A'
                        residue_num = int(parts[8]) if parts[8] != '?' else int(parts[7])
                        x = float(parts[10])
                        y = float(parts[11])
                        z = float(parts[12])
                        
                        # Check if it's an RNA residue
                        if residue_name not in PDBParser.RNA_RESIDUES:
                            continue
                        
                        # Normalize residue name
                        residue_name = PDBParser.normalize_residue_name(residue_name)
                        
                        atom = Atom(atom_name, x, y, z, residue_name, chain_id, residue_num)
                        
                        key = (chain_id, residue_num)
                        if key not in nucleotides_dict:
                            nucleotides_dict[key] = Nucleotide(
                                residue_name, chain_id, residue_num, {}
                            )
                        nucleotides_dict[key].atoms[atom_name] = atom
                    
                    except (ValueError, IndexError):
                        continue
        
        # Convert to list and sort
        nucleotides = list(nucleotides_dict.values())
        nucleotides.sort(key=lambda n: (n.chain_id, n.residue_num))
        
        return nucleotides
    
    @staticmethod
    def parse_pdb(filename: str) -> List[Nucleotide]:
        """Parse standard PDB file format"""
        nucleotides_dict = {}
        
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    try:
                        atom_name = line[12:16].strip()
                        residue_name = line[17:20].strip()
                        chain_id = line[21]
                        residue_num = int(line[22:26].strip())
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        
                        # Check if it's an RNA residue
                        if residue_name not in PDBParser.RNA_RESIDUES:
                            continue
                        
                        # Normalize residue name
                        residue_name = PDBParser.normalize_residue_name(residue_name)
                        
                        atom = Atom(atom_name, x, y, z, residue_name, chain_id, residue_num)
                        
                        key = (chain_id, residue_num)
                        if key not in nucleotides_dict:
                            nucleotides_dict[key] = Nucleotide(
                                residue_name, chain_id, residue_num, {}
                            )
                        nucleotides_dict[key].atoms[atom_name] = atom
                    
                    except (ValueError, IndexError):
                        continue
        
        # Convert to list and sort
        nucleotides = list(nucleotides_dict.values())
        nucleotides.sort(key=lambda n: (n.chain_id, n.residue_num))
        
        return nucleotides
    
    @staticmethod
    def normalize_residue_name(residue_name: str) -> str:
        """Normalize residue name to single letter"""
        if residue_name in ["ADE", "rA", "DA"]:
            return "A"
        elif residue_name in ["URA", "rU"]:
            return "U"
        elif residue_name in ["GUA", "rG", "DG"]:
            return "G"
        elif residue_name in ["CYT", "rC", "DC"]:
            return "C"
        elif residue_name == "DT":
            return "U"  # Treat thymine as uracil for RNA
        else:
            return residue_name

class BasePairDetector:
    """Detect base pairs from 3D structure using multiple criteria"""
    
    # Distance thresholds (in Angstroms)
    C1_C1_MIN = 9.0
    C1_C1_MAX = 12.5
    N1_N9_MIN = 7.5  # New: N1/N9 distance range
    N1_N9_MAX = 11.5  # New: N1/N9 distance range
    HBOND_MIN = 2.5
    HBOND_MAX = 3.5
    
    # Specific H-bond patterns for each base pair type
    HBOND_PATTERNS = {
        ('A', 'U'): [('N6', 'O4'), ('N1', 'N3')],
        ('U', 'A'): [('O4', 'N6'), ('N3', 'N1')],
        ('G', 'C'): [('O6', 'N4'), ('N1', 'N3'), ('N2', 'O2')],
        ('C', 'G'): [('N4', 'O6'), ('N3', 'N1'), ('O2', 'N2')],
        ('G', 'U'): [('O6', 'N3'), ('N1', 'O2')],  # Wobble pair
        ('U', 'G'): [('N3', 'O6'), ('O2', 'N1')],  # Wobble pair
    }
    
    @staticmethod
    def calculate_distance(atom1: Atom, atom2: Atom) -> float:
        """Calculate Euclidean distance between two atoms"""
        return np.linalg.norm(atom1.coords - atom2.coords)
    
    @staticmethod
    def check_coplanarity(nuc1: Nucleotide, nuc2: Nucleotide, threshold: float = 45.0) -> bool:
        """Check if two nucleotides are approximately coplanar"""
        base1_atoms = nuc1.get_base_atoms()
        base2_atoms = nuc2.get_base_atoms()
        
        if len(base1_atoms) < 3 or len(base2_atoms) < 3:
            return False
        
        # Fit planes to each base
        coords1 = np.array([a.coords for a in base1_atoms[:3]])
        coords2 = np.array([a.coords for a in base2_atoms[:3]])
        
        # Calculate normal vectors
        v1 = coords1[1] - coords1[0]
        v2 = coords1[2] - coords1[0]
        normal1 = np.cross(v1, v2)
        normal1 = normal1 / np.linalg.norm(normal1)
        
        v3 = coords2[1] - coords2[0]
        v4 = coords2[2] - coords2[0]
        normal2 = np.cross(v3, v4)
        normal2 = normal2 / np.linalg.norm(normal2)
        
        # Calculate angle between normals
        cos_angle = abs(np.dot(normal1, normal2))
        angle = np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))
        
        # Planes are coplanar if normals are parallel (angle ~0) or antiparallel (angle ~180)
        return angle < threshold or angle > (180 - threshold)
    
    @staticmethod
    def check_hbond_geometry(donor: Atom, acceptor: Atom, angle_threshold: float = 120.0) -> bool:
        """Check if H-bond geometry is reasonable (simplified check)"""
        # For a more complete check, we'd need the hydrogen position
        # Here we just check the distance
        distance = BasePairDetector.calculate_distance(donor, acceptor)
        return BasePairDetector.HBOND_MIN <= distance <= BasePairDetector.HBOND_MAX
    
    @staticmethod
    def check_n1_n9_distance(nuc1: Nucleotide, nuc2: Nucleotide) -> bool:
        """Check N1/N9 glycosidic atom distance criterion"""
        glyc1 = nuc1.get_glycosidic_atom()
        glyc2 = nuc2.get_glycosidic_atom()
        
        if not glyc1 or not glyc2:
            return False
        
        distance = BasePairDetector.calculate_distance(glyc1, glyc2)
        return BasePairDetector.N1_N9_MIN <= distance <= BasePairDetector.N1_N9_MAX
    
    @staticmethod
    def check_c1_c1_distance(nuc1: Nucleotide, nuc2: Nucleotide) -> bool:
        """Check C1'-C1' distance criterion"""
        c1_1 = nuc1.get_atom("C1'")
        c1_2 = nuc2.get_atom("C1'")
        
        if not c1_1 or not c1_2:
            return False
        
        distance = BasePairDetector.calculate_distance(c1_1, c1_2)
        return BasePairDetector.C1_C1_MIN <= distance <= BasePairDetector.C1_C1_MAX
    
    @staticmethod
    def is_base_paired(nuc1: Nucleotide, nuc2: Nucleotide, 
                       require_coplanarity: bool = True,
                       distance_method: str = 'both') -> bool:
        """Determine if two nucleotides are base paired
        
        Args:
            nuc1, nuc2: Nucleotides to check
            require_coplanarity: Whether to require coplanarity check
            distance_method: 'c1_c1', 'n1_n9', or 'both' (either criterion)
        """
        
        # Skip if too close in sequence (avoid backbone neighbors)
        if nuc1.chain_id == nuc2.chain_id and abs(nuc1.residue_num - nuc2.residue_num) <= 3:
            return False
        
        # Check base pair type compatibility
        pair_type = (nuc1.residue_name, nuc2.residue_name)
        if pair_type not in BasePairDetector.HBOND_PATTERNS:
            return False
        
        # Check distance criteria based on method
        distance_ok = False
        if distance_method == 'c1_c1':
            distance_ok = BasePairDetector.check_c1_c1_distance(nuc1, nuc2)
        elif distance_method == 'n1_n9':
            distance_ok = BasePairDetector.check_n1_n9_distance(nuc1, nuc2)
        elif distance_method == 'both':
            # Either criterion can satisfy the distance requirement
            distance_ok = (BasePairDetector.check_c1_c1_distance(nuc1, nuc2) or 
                          BasePairDetector.check_n1_n9_distance(nuc1, nuc2))
        else:
            raise ValueError(f"Unknown distance_method: {distance_method}")
        
        if not distance_ok:
            return False
        
        # Check coplanarity if required (with relaxed threshold)
        if require_coplanarity and not BasePairDetector.check_coplanarity(nuc1, nuc2, threshold=35.0):
            return False
        
        # Check specific H-bond patterns
        hbond_patterns = BasePairDetector.HBOND_PATTERNS[pair_type]
        hbonds_found = 0
        
        for donor_name, acceptor_name in hbond_patterns:
            donor = nuc1.get_atom(donor_name)
            acceptor = nuc2.get_atom(acceptor_name)
            
            if donor and acceptor:
                if BasePairDetector.check_hbond_geometry(donor, acceptor):
                    hbonds_found += 1
        
        # More lenient: require at least 1 H-bond for GU pairs, 2 for others
        min_hbonds = 1 if pair_type in [('G', 'U'), ('U', 'G')] else 2
        return hbonds_found >= min_hbonds

class BasePairingMatrix:
    """Represent and manipulate base pairing matrices"""
    
    def __init__(self, nucleotides: List[Nucleotide]):
        self.nucleotides = nucleotides
        self.n = len(nucleotides)
        self.matrix = np.zeros((self.n, self.n), dtype=int)
        self.base_pairs = set()
    
    def detect_base_pairs(self, require_coplanarity: bool = True, 
                         remove_isolated: bool = True,
                         distance_method: str = 'both'):
        """Detect all base pairs and build the matrix
        
        Args:
            require_coplanarity: Whether to require coplanarity check
            remove_isolated: Whether to remove isolated base pairs
            distance_method: 'c1_c1', 'n1_n9', or 'both' for distance criteria
        """
        print(f"Detecting base pairs among {self.n} nucleotides using {distance_method} distance method...")
        
        for i in range(self.n):
            for j in range(i+1, self.n):
                if BasePairDetector.is_base_paired(
                    self.nucleotides[i], 
                    self.nucleotides[j],
                    require_coplanarity,
                    distance_method
                ):
                    self.matrix[i, j] = 1
                    self.matrix[j, i] = 1
                    self.base_pairs.add((i, j))
        
        if remove_isolated:
            self.remove_isolated_pairs()
        
        print(f"Found {len(self.base_pairs)} base pairs")
    
    def remove_isolated_pairs(self):
        """Remove base pairs that have no neighboring pairs (more lenient)"""
        pairs_to_remove = set()
        
        for i, j in self.base_pairs:
            has_neighbor = False
            
            # Check for neighboring pairs with wider range
            for di in [-2, -1, 0, 1, 2]:
                for dj in [-2, -1, 0, 1, 2]:
                    if di == 0 and dj == 0:
                        continue
                    ni, nj = i + di, j + dj
                    if 0 <= ni < self.n and 0 <= nj < self.n:
                        if (ni, nj) in self.base_pairs or (nj, ni) in self.base_pairs:
                            has_neighbor = True
                            break
                if has_neighbor:
                    break
            
            if not has_neighbor:
                pairs_to_remove.add((i, j))
        
        # Remove isolated pairs
        for i, j in pairs_to_remove:
            self.matrix[i, j] = 0
            self.matrix[j, i] = 0
            self.base_pairs.discard((i, j))
        
        if pairs_to_remove:
            print(f"Removed {len(pairs_to_remove)} isolated base pairs")
    
    def get_dot_bracket_notation(self) -> str:
        """Convert matrix to dot-bracket notation"""
        notation = ['.'] * self.n
        used = [False] * self.n
        
        # Simple approach - doesn't handle pseudoknots
        for i, j in sorted(self.base_pairs):
            if not used[i] and not used[j]:
                notation[i] = '('
                notation[j] = ')'
                used[i] = used[j] = True
        
        return ''.join(notation)
    
    def get_base_pair_list(self) -> List[Tuple[int, int, str, str]]:
        """Get list of base pairs with nucleotide information"""
        pair_list = []
        for i, j in sorted(self.base_pairs):
            nuc1 = self.nucleotides[i]
            nuc2 = self.nucleotides[j]
            pair_list.append((
                i+1,  # 1-indexed
                j+1,  # 1-indexed
                f"{nuc1.chain_id}:{nuc1.residue_num}{nuc1.residue_name}",
                f"{nuc2.chain_id}:{nuc2.residue_num}{nuc2.residue_name}"
            ))
        return pair_list
    
    def save_matrix(self, filename: str):
        """Save matrix to file"""
        np.savetxt(filename, self.matrix, fmt='%d')
        print(f"Matrix saved to {filename}")
    
    def save_base_pairs(self, filename: str):
        """Save base pair list to file"""
        with open(filename, 'w') as f:
            f.write("# Base pairs (1-indexed)\n")
            f.write("# i\tj\tNuc1\tNuc2\n")
            for i, j, nuc1_str, nuc2_str in self.get_base_pair_list():
                f.write(f"{i}\t{j}\t{nuc1_str}\t{nuc2_str}\n")
        print(f"Base pair list saved to {filename}")
    
    def print_summary(self):
        """Print summary statistics"""
        print("\n=== Base Pairing Summary ===")
        print(f"Total nucleotides: {self.n}")
        print(f"Total base pairs: {len(self.base_pairs)}")
        
        # Count pair types
        pair_types = {}
        for i, j in self.base_pairs:
            pair_type = tuple(sorted([
                self.nucleotides[i].residue_name,
                self.nucleotides[j].residue_name
            ]))
            pair_types[pair_type] = pair_types.get(pair_type, 0) + 1
        
        print("\nBase pair types:")
        for pair_type, count in sorted(pair_types.items()):
            print(f"  {'-'.join(pair_type)}: {count}")
        
        print(f"\nDot-bracket notation:\n{self.get_dot_bracket_notation()}")

def main():
    parser = argparse.ArgumentParser(
        description='Extract RNA base pairing matrix from PDB file (v3.0 with N1/N9 distance option)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s input.pdb -o matrix.txt
  %(prog)s input.pdb --distance-method n1_n9
  %(prog)s input.pdb --distance-method both --no-coplanarity
  %(prog)s input.pdb -o matrix.txt -p pairs.txt -s
  %(prog)s --input-dir /path/to/pdbs --output-dir /path/to/matrices --distance-method both
        """
    )
    
    parser.add_argument('pdb_file', nargs='?', help='Input PDB/CIF file')
    parser.add_argument('-o', '--output-matrix', help='Output matrix file')
    parser.add_argument('-p', '--output-pairs', help='Output base pair list file')
    parser.add_argument('--distance-method', choices=['c1_c1', 'n1_n9', 'both'], 
                       default='both',
                       help='Distance criterion: c1_c1 (C1\'-C1\'), n1_n9 (N1/N9), or both (default: both)')
    parser.add_argument('--no-coplanarity', action='store_true',
                       help='Disable coplanarity check')
    parser.add_argument('--keep-isolated', action='store_true',
                       help='Keep isolated base pairs')
    parser.add_argument('-s', '--summary', action='store_true',
                       help='Print summary statistics')
    # New: directory processing arguments
    parser.add_argument('--input-dir', help='Input directory containing PDB/CIF files')
    parser.add_argument('--output-dir', help='Output directory to write matrix files (one per input)')
    
    args = parser.parse_args()
    
    # Directory mode
    if args.input_dir:
        in_dir = Path(args.input_dir)
        if not in_dir.exists() or not in_dir.is_dir():
            print(f"Error: Input directory '{args.input_dir}' not found or not a directory", file=sys.stderr)
            sys.exit(1)
        if not args.output_dir:
            print("Error: --output-dir is required when using --input-dir", file=sys.stderr)
            sys.exit(1)
        out_dir = Path(args.output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        
        pdb_files = sorted([p for p in in_dir.iterdir() if p.suffix.lower() in ('.pdb', '.cif')])
        if not pdb_files:
            print(f"Warning: No .pdb or .cif files found in '{in_dir}'", file=sys.stderr)
            sys.exit(0)
        
        total = len(pdb_files)
        processed = 0
        skipped = 0
        errors = 0
        print(f"Found {total} structure files in '{in_dir}'. Writing matrices to '{out_dir}'...")
        
        for pdb_path in pdb_files:
            try:
                print(f"\nParsing PDB/CIF file: {pdb_path}")
                nucleotides = PDBParser.parse(str(pdb_path))
                
                if not nucleotides:
                    print(f"Warning: No RNA nucleotides found in {pdb_path}, skipping.", file=sys.stderr)
                    skipped += 1
                    continue
                
                bp_matrix = BasePairingMatrix(nucleotides)
                bp_matrix.detect_base_pairs(
                    require_coplanarity=not args.no_coplanarity,
                    remove_isolated=not args.keep_isolated,
                    distance_method=args.distance_method
                )
                
                out_matrix = out_dir / f"{pdb_path.stem}.txt"
                bp_matrix.save_matrix(str(out_matrix))
                processed += 1
            except Exception as e:
                print(f"Error processing {pdb_path}: {e}", file=sys.stderr)
                errors += 1
                continue
        
        print("\n=== Batch Processing Summary ===")
        print(f"Total files: {total}")
        print(f"Processed: {processed}")
        print(f"Skipped (no RNA): {skipped}")
        print(f"Errors: {errors}")
        return
    
    # Single-file mode
    if not args.pdb_file:
        print("Error: PDB file is required in single-file mode. Provide a file path or use --input-dir.", file=sys.stderr)
        sys.exit(1)
    
    # Check if input file exists
    if not Path(args.pdb_file).exists():
        print(f"Error: PDB file '{args.pdb_file}' not found", file=sys.stderr)
        sys.exit(1)
    
    try:
        # Parse PDB file
        print(f"Parsing PDB file: {args.pdb_file}")
        nucleotides = PDBParser.parse(args.pdb_file)
        
        if not nucleotides:
            print("Error: No RNA nucleotides found in PDB file", file=sys.stderr)
            sys.exit(1)
        
        print(f"Found {len(nucleotides)} nucleotides")
        
        # Create base pairing matrix
        bp_matrix = BasePairingMatrix(nucleotides)
        
        # Detect base pairs
        bp_matrix.detect_base_pairs(
            require_coplanarity=not args.no_coplanarity,
            remove_isolated=not args.keep_isolated,
            distance_method=args.distance_method
        )
        
        # Save outputs
        if args.output_matrix:
            bp_matrix.save_matrix(args.output_matrix)
        
        if args.output_pairs:
            bp_matrix.save_base_pairs(args.output_pairs)
        
        # Print summary if requested
        if args.summary or (not args.output_matrix and not args.output_pairs):
            bp_matrix.print_summary()
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()