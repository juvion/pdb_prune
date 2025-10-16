#!/usr/bin/env python3
"""
RNA Base Pairing Matrix Tool
Extracts base pairing information from RNA PDB structures
Supports both PDB and CIF/mmCIF formats
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

class PDBParser:
    """Parse PDB and CIF files and extract nucleotide information"""
    
    RNA_RESIDUES = {"A", "U", "G", "C", "ADE", "URA", "GUA", "CYT", 
                    "rA", "rU", "rG", "rC", "DA", "DT", "DG", "DC"}
    
    @staticmethod
    def parse(filename: str) -> List[Nucleotide]:
        """Parse PDB or CIF file and return list of nucleotides"""
        # Determine file format
        with open(filename, 'r') as f:
            first_lines = f.read(100)
            if 'data_' in first_lines or '_atom_site' in first_lines:
                return PDBParser.parse_cif(filename)
            else:
                return PDBParser.parse_pdb(filename)
    
    @staticmethod
    def parse_cif(filename: str) -> List[Nucleotide]:
        """Parse CIF/mmCIF format file"""
        nucleotides_dict = {}
        
        with open(filename, 'r') as f:
            lines = f.readlines()
            
        # Find the atom_site loop
        in_atom_loop = False
        atom_data = []
        
        for i, line in enumerate(lines):
            if 'loop_' in line and i+1 < len(lines) and '_atom_site' in lines[i+1]:
                in_atom_loop = True
                # Find column headers
                headers = []
                j = i + 1
                while j < len(lines) and lines[j].startswith('_atom_site'):
                    headers.append(lines[j].strip())
                    j += 1
                # Start of data
                data_start = j
                continue
            
            if in_atom_loop:
                if line.strip() and not line.startswith('_') and not line.startswith('#') and 'loop_' not in line:
                    parts = line.split()
                    if len(parts) >= 15:  # Ensure we have enough fields
                        atom_data.append(parts)
                elif line.startswith('_') and '_atom_site' not in line:
                    in_atom_loop = False
        
        # Process atom data
        for parts in atom_data:
            try:
                # Standard positions for CIF format
                group = parts[0]  # ATOM or HETATM
                if group not in ['ATOM', 'HETATM']:
                    continue
                    
                atom_name = parts[3]
                residue_name = parts[5]
                chain_id = parts[6] if parts[6] != '.' else 'A'
                residue_num = int(parts[8]) if parts[8] != '?' else 1
                x = float(parts[10])
                y = float(parts[11])
                z = float(parts[12])
                
                # Check if it's an RNA residue
                if residue_name not in PDBParser.RNA_RESIDUES:
                    continue
                
                # Normalize residue name
                if residue_name in ["ADE", "rA", "DA"]:
                    residue_name = "A"
                elif residue_name in ["URA", "rU"]:
                    residue_name = "U"
                elif residue_name in ["GUA", "rG", "DG"]:
                    residue_name = "G"
                elif residue_name in ["CYT", "rC", "DC"]:
                    residue_name = "C"
                elif residue_name == "DT":
                    residue_name = "U"
                
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
        """Parse standard PDB format file"""
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
                        if residue_name in ["ADE", "rA", "DA"]:
                            residue_name = "A"
                        elif residue_name in ["URA", "rU"]:
                            residue_name = "U"
                        elif residue_name in ["GUA", "rG", "DG"]:
                            residue_name = "G"
                        elif residue_name in ["CYT", "rC", "DC"]:
                            residue_name = "C"
                        elif residue_name == "DT":
                            residue_name = "U"
                        
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

class BasePairDetector:
    """Detect base pairs from 3D structure"""
    
    # Adjusted distance thresholds (in Angstroms)
    C1_C1_MIN = 10.0  # Slightly relaxed
    C1_C1_MAX = 11.5  # Slightly relaxed
    HBOND_MIN = 2.5
    HBOND_MAX = 3.2   # Slightly relaxed
    
    # Specific H-bond patterns for each base pair type
    HBOND_PATTERNS = {
        ('A', 'U'): [('N6', 'O4', 3.0), ('N1', 'N3', 2.9)],
        ('U', 'A'): [('O4', 'N6', 3.0), ('N3', 'N1', 2.9)],
        ('G', 'C'): [('O6', 'N4', 2.9), ('N1', 'N3', 2.95), ('N2', 'O2', 2.9)],
        ('C', 'G'): [('N4', 'O6', 2.9), ('N3', 'N1', 2.95), ('O2', 'N2', 2.9)],
        ('G', 'U'): [('O6', 'N3', 2.85), ('N1', 'O2', 2.85)],  # Wobble pair
        ('U', 'G'): [('N3', 'O6', 2.85), ('O2', 'N1', 2.85)],  # Wobble pair
    }
    
    @staticmethod
    def calculate_distance(atom1: Atom, atom2: Atom) -> float:
        """Calculate Euclidean distance between two atoms"""
        return np.linalg.norm(atom1.coords - atom2.coords)
    
    @staticmethod
    def check_coplanarity(nuc1: Nucleotide, nuc2: Nucleotide, threshold: float = 35.0) -> bool:
        """Check if two nucleotides are approximately coplanar"""
        base1_atoms = nuc1.get_base_atoms()
        base2_atoms = nuc2.get_base_atoms()
        
        if len(base1_atoms) < 3 or len(base2_atoms) < 3:
            return False
        
        try:
            # Fit planes to each base
            coords1 = np.array([a.coords for a in base1_atoms[:3]])
            coords2 = np.array([a.coords for a in base2_atoms[:3]])
            
            # Calculate normal vectors
            v1 = coords1[1] - coords1[0]
            v2 = coords1[2] - coords1[0]
            normal1 = np.cross(v1, v2)
            norm1 = np.linalg.norm(normal1)
            if norm1 > 0:
                normal1 = normal1 / norm1
            else:
                return False
            
            v3 = coords2[1] - coords2[0]
            v4 = coords2[2] - coords2[0]
            normal2 = np.cross(v3, v4)
            norm2 = np.linalg.norm(normal2)
            if norm2 > 0:
                normal2 = normal2 / norm2
            else:
                return False
            
            # Calculate angle between normals
            cos_angle = abs(np.dot(normal1, normal2))
            angle = np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))
            
            # Planes are coplanar if normals are parallel (angle ~0) or antiparallel (angle ~180)
            return angle < threshold or angle > (180 - threshold)
        except:
            return False
    
    @staticmethod
    def check_hbond_geometry(donor: Atom, acceptor: Atom, max_dist: float = 3.2) -> bool:
        """Check if H-bond geometry is reasonable"""
        distance = BasePairDetector.calculate_distance(donor, acceptor)
        return BasePairDetector.HBOND_MIN <= distance <= max_dist
    
    @staticmethod
    def is_base_paired(nuc1: Nucleotide, nuc2: Nucleotide, 
                       require_coplanarity: bool = True,
                       strict: bool = False) -> bool:
        """Determine if two nucleotides are base paired"""
        
        # Skip if too close in sequence (avoid backbone neighbors)
        if nuc1.chain_id == nuc2.chain_id and abs(nuc1.residue_num - nuc2.residue_num) <= 3:
            return False
        
        # Check C1'-C1' distance
        c1_1 = nuc1.get_atom("C1'")
        c1_2 = nuc2.get_atom("C1'")
        
        if not c1_1 or not c1_2:
            return False
        
        c1_distance = BasePairDetector.calculate_distance(c1_1, c1_2)
        if not (BasePairDetector.C1_C1_MIN <= c1_distance <= BasePairDetector.C1_C1_MAX):
            return False
        
        # Check coplanarity if required (relaxed for wobble pairs)
        if require_coplanarity:
            pair_type = (nuc1.residue_name, nuc2.residue_name)
            threshold = 35.0 if pair_type in [('G', 'U'), ('U', 'G')] else 30.0
            if not BasePairDetector.check_coplanarity(nuc1, nuc2, threshold):
                return False
        
        # Check specific H-bond patterns
        pair_type = (nuc1.residue_name, nuc2.residue_name)
        if pair_type not in BasePairDetector.HBOND_PATTERNS:
            return False
        
        hbond_patterns = BasePairDetector.HBOND_PATTERNS[pair_type]
        hbonds_found = 0
        
        for pattern in hbond_patterns:
            donor_name = pattern[0]
            acceptor_name = pattern[1]
            max_dist = pattern[2] if len(pattern) > 2 else BasePairDetector.HBOND_MAX
            
            donor = nuc1.get_atom(donor_name)
            acceptor = nuc2.get_atom(acceptor_name)
            
            if donor and acceptor:
                if BasePairDetector.check_hbond_geometry(donor, acceptor, max_dist):
                    hbonds_found += 1
        
        # Require at least 2 H-bonds for WC pairs, 1 for wobble pairs
        min_hbonds = 2 if strict else (1 if pair_type in [('G', 'U'), ('U', 'G')] else 2)
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
                         strict: bool = False):
        """Detect all base pairs and build the matrix"""
        print(f"Detecting base pairs among {self.n} nucleotides...")
        
        for i in range(self.n):
            for j in range(i+1, self.n):
                if BasePairDetector.is_base_paired(
                    self.nucleotides[i], 
                    self.nucleotides[j],
                    require_coplanarity,
                    strict
                ):
                    self.matrix[i, j] = 1
                    self.matrix[j, i] = 1
                    self.base_pairs.add((i, j))
        
        if remove_isolated:
            self.remove_isolated_pairs()
        
        print(f"Found {len(self.base_pairs)} base pairs")
    
    def remove_isolated_pairs(self):
        """Remove base pairs that have no neighboring pairs"""
        pairs_to_remove = set()
        
        for i, j in self.base_pairs:
            has_neighbor = False
            
            # Check for neighboring pairs (more lenient for terminal pairs)
            for di in [-1, 0, 1]:
                for dj in [-1, 0, 1]:
                    if di == 0 and dj == 0:
                        continue
                    ni, nj = i + di, j + dj
                    if 0 <= ni < self.n and 0 <= nj < self.n:
                        if (ni, nj) in self.base_pairs or (nj, ni) in self.base_pairs:
                            has_neighbor = True
                            break
                if has_neighbor:
                    break
            
            # Be more lenient with terminal base pairs
            if not has_neighbor:
                # Check if it's a terminal pair
                is_terminal = (i <= 1 or i >= self.n-2) or (j <= 1 or j >= self.n-2)
                if not is_terminal:
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
            nuc1_name = self.nucleotides[i].residue_name
            nuc2_name = self.nucleotides[j].residue_name
            
            # Normalize pair type
            if (nuc1_name == 'A' and nuc2_name == 'U') or (nuc1_name == 'U' and nuc2_name == 'A'):
                pair_key = 'A-U'
            elif (nuc1_name == 'G' and nuc2_name == 'C') or (nuc1_name == 'C' and nuc2_name == 'G'):
                pair_key = 'C-G'
            elif (nuc1_name == 'G' and nuc2_name == 'U') or (nuc1_name == 'U' and nuc2_name == 'G'):
                pair_key = 'U-G'
            else:
                pair_key = '-'.join(sorted([nuc1_name, nuc2_name]))
            
            pair_types[pair_key] = pair_types.get(pair_key, 0) + 1
        
        print("\nBase pair types:")
        for pair_type in ['A-U', 'C-G', 'U-G']:  # Show in standard order
            if pair_type in pair_types:
                print(f"  {pair_type}: {pair_types[pair_type]}")
        
        print(f"\nDot-bracket notation:\n{self.get_dot_bracket_notation()}")

def main():
    parser = argparse.ArgumentParser(
        description='Extract RNA base pairing matrix from PDB/CIF file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s input.pdb -o matrix.txt
  %(prog)s input.cif --no-coplanarity --keep-isolated
  %(prog)s input.pdb -o matrix.txt -p pairs.txt -s
  %(prog)s input.pdb --relaxed  # Use relaxed parameters for difficult structures
        """
    )
    
    parser.add_argument('pdb_file', help='Input PDB or CIF file')
    parser.add_argument('-o', '--output-matrix', help='Output matrix file')
    parser.add_argument('-p', '--output-pairs', help='Output base pair list file')
    parser.add_argument('--no-coplanarity', action='store_true',
                       help='Disable coplanarity check')
    parser.add_argument('--keep-isolated', action='store_true',
                       help='Keep isolated base pairs')
    parser.add_argument('--strict', action='store_true',
                       help='Use strict base pairing criteria')
    parser.add_argument('--relaxed', action='store_true',
                       help='Use relaxed parameters (better for wobble pairs)')
    parser.add_argument('-s', '--summary', action='store_true',
                       help='Print summary statistics')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not Path(args.pdb_file).exists():
        print(f"Error: File '{args.pdb_file}' not found", file=sys.stderr)
        sys.exit(1)
    
    try:
        # Parse PDB/CIF file
        print(f"Parsing file: {args.pdb_file}")
        nucleotides = PDBParser.parse(args.pdb_file)
        
        if not nucleotides:
            print("Error: No RNA nucleotides found in file", file=sys.stderr)
            sys.exit(1)
        
        print(f"Found {len(nucleotides)} nucleotides")
        
        # Create base pairing matrix
        bp_matrix = BasePairingMatrix(nucleotides)
        
        # Adjust parameters based on flags
        if args.relaxed:
            args.no_coplanarity = True
            args.keep_isolated = True
        
        # Detect base pairs
        bp_matrix.detect_base_pairs(
            require_coplanarity=not args.no_coplanarity,
            remove_isolated=not args.keep_isolated,
            strict=args.strict
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
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()