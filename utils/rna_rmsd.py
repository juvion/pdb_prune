#!/usr/bin/env python3
"""
RNA RMSD Tool - Compute RMSD between two RNA PDB structures

This tool computes the Root Mean Square Deviation (RMSD) between two RNA
structures after optimal superposition using the Kabsch algorithm.
"""

import sys
import argparse
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select
from typing import List, Tuple, Optional, Set


# Standard RNA backbone atom names
BACKBONE_ATOMS = {'P', "O5'", "C5'", "C4'", "C3'", "O3'"}

# Nucleic acid residue names
RNA_RESIDUES = {'A', 'U', 'G', 'C', 'I'}
DNA_RESIDUES = {'DA', 'DG', 'DC', 'DT'}
NUCLEIC_RESIDUES = RNA_RESIDUES | DNA_RESIDUES


def log_message(message: str, quiet: bool = False):
    """Print diagnostic message to stderr unless quiet mode is enabled."""
    if not quiet:
        print(message, file=sys.stderr)


def get_atom_element(atom) -> str:
    """Get the element of an atom, handling different naming conventions."""
    if hasattr(atom, 'element') and atom.element:
        return atom.element.strip()
    # Fallback: parse from atom name
    name = atom.get_name().strip()
    if name[0].isalpha():
        return name[0]
    return ''


def select_best_altloc(residue):
    """
    For residues with alternate locations, keep only the highest occupancy
    conformer (or first if occupancy is equal/missing).
    """
    atoms_by_name = {}
    for atom in residue.get_atoms():
        name = atom.get_name()
        altloc = atom.get_altloc()
        
        if name not in atoms_by_name:
            atoms_by_name[name] = []
        atoms_by_name[name].append(atom)
    
    # For each atom name, select best altloc
    best_atoms = []
    for name, atom_list in atoms_by_name.items():
        if len(atom_list) == 1:
            best_atoms.append(atom_list[0])
        else:
            # Select by occupancy
            best_atom = max(atom_list, key=lambda a: a.get_occupancy())
            best_atoms.append(best_atom)
    
    return best_atoms


def filter_nucleic_residues(chain, quiet: bool = False):
    """Extract nucleic acid residues from a chain, handling altlocs."""
    nucleic_residues = []
    
    for residue in chain:
        # Skip hetero residues (water, ions, ligands)
        if residue.id[0] != ' ':
            continue
            
        resname = residue.get_resname().strip()
        if resname in NUCLEIC_RESIDUES:
            nucleic_residues.append(residue)
    
    return nucleic_residues


def select_chain(structure, chain_id: Optional[str], quiet: bool = False) -> Tuple[any, List]:
    """
    Select a chain and return its nucleic residues.
    If chain_id is None, auto-select the chain with most nucleic residues.
    """
    model = structure[0]  # Use first model only
    
    if chain_id:
        if chain_id not in model:
            raise ValueError(f"Chain {chain_id} not found in structure")
        chain = model[chain_id]
        residues = filter_nucleic_residues(chain, quiet)
        if not residues:
            raise ValueError(f"No nucleic residues found in chain {chain_id}")
        return chain, residues
    else:
        # Auto-select chain with most nucleic residues
        best_chain = None
        best_residues = []
        
        for chain in model:
            residues = filter_nucleic_residues(chain, quiet)
            if len(residues) > len(best_residues):
                best_chain = chain
                best_residues = residues
        
        if not best_residues:
            raise ValueError("No nucleic residues found in any chain")
        
        log_message(f"Auto-selected chain {best_chain.id} with {len(best_residues)} nucleic residues", quiet)
        return best_chain, best_residues


def apply_residue_range(residues: List, res_range: Optional[str], quiet: bool = False) -> List:
    """Apply residue range filter (1-based inclusive)."""
    if not res_range:
        return residues
    
    try:
        start, end = map(int, res_range.split(':'))
    except ValueError:
        raise ValueError(f"Invalid residue range format: {res_range}. Expected 'start:end'")
    
    # Convert to 0-based indexing
    start_idx = start - 1
    end_idx = end  # end is inclusive, so we don't subtract 1
    
    if start_idx < 0 or end_idx > len(residues):
        raise ValueError(f"Residue range {start}:{end} out of bounds (1-{len(residues)})")
    
    filtered = residues[start_idx:end_idx]
    log_message(f"Applied residue range {start}:{end}, selected {len(filtered)} residues", quiet)
    return filtered


def get_residue_atoms(residue, selection: str, strip_h: bool) -> List[Tuple[str, np.ndarray]]:
    """
    Get atoms from a residue based on selection criteria.
    Returns list of (atom_name, coordinates) tuples.
    """
    # Handle altlocs - select best conformer
    atoms = select_best_altloc(residue)
    
    result = []
    for atom in atoms:
        atom_name = atom.get_name().strip()
        
        # Filter by selection type
        if selection == 'backbone':
            if atom_name not in BACKBONE_ATOMS:
                continue
        # For 'all', include all atoms (already filtered to nucleic residues)
        
        # Strip hydrogens if requested
        if strip_h:
            element = get_atom_element(atom)
            if element == 'H':
                continue
        
        coord = atom.get_coord()
        result.append((atom_name, coord))
    
    return result


def collect_atom_pairs(residues1: List, residues2: List, selection: str, 
                       strip_h: bool, quiet: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    """
    Collect paired atom coordinates from two residue lists.
    Returns two arrays of shape (N, 3) where N is the total number of paired atoms.
    """
    if len(residues1) != len(residues2):
        raise ValueError(
            f"Chain length mismatch: structure 1 has {len(residues1)} residues, "
            f"structure 2 has {len(residues2)} residues"
        )
    
    X_list = []  # Coordinates from structure 1
    Y_list = []  # Coordinates from structure 2
    skipped_residues = 0
    total_residues = len(residues1)
    
    for i, (res1, res2) in enumerate(zip(residues1, residues2)):
        atoms1 = get_residue_atoms(res1, selection, strip_h)
        atoms2 = get_residue_atoms(res2, selection, strip_h)
        
        if not atoms1 or not atoms2:
            skipped_residues += 1
            continue
        
        # Create dictionaries for pairing by atom name
        atoms1_dict = {name: coord for name, coord in atoms1}
        atoms2_dict = {name: coord for name, coord in atoms2}
        
        # Find common atoms
        common_atoms = set(atoms1_dict.keys()) & set(atoms2_dict.keys())
        
        if not common_atoms:
            skipped_residues += 1
            continue
        
        # Sort atom names for deterministic ordering
        for atom_name in sorted(common_atoms):
            X_list.append(atoms1_dict[atom_name])
            Y_list.append(atoms2_dict[atom_name])
    
    if skipped_residues > 0:
        log_message(f"Skipped {skipped_residues}/{total_residues} residues with no common atoms", quiet)
    
    if not X_list:
        raise ValueError("No overlapping atoms found for selection")
    
    X = np.array(X_list)
    Y = np.array(Y_list)
    
    log_message(f"Paired {len(X)} atoms from {total_residues - skipped_residues} residues", quiet)
    
    return X, Y


def kabsch_rmsd(X: np.ndarray, Y: np.ndarray) -> Tuple[float, np.ndarray, np.ndarray]:
    """
    Compute RMSD between two coordinate sets using the Kabsch algorithm.
    
    Args:
        X: Reference coordinates (N, 3)
        Y: Coordinates to align (N, 3)
    
    Returns:
        rmsd: Root mean square deviation
        R: Rotation matrix (3, 3)
        t: Translation vector (3,)
    """
    # Center both coordinate sets
    centroid_X = np.mean(X, axis=0)
    centroid_Y = np.mean(Y, axis=0)
    
    X_centered = X - centroid_X
    Y_centered = Y - centroid_Y
    
    # Compute covariance matrix
    C = Y_centered.T @ X_centered
    
    # Singular Value Decomposition
    U, S, Vt = np.linalg.svd(C)
    
    # Compute rotation matrix
    R = Vt.T @ U.T
    
    # Handle reflection case
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    
    # Translation vector
    t = centroid_X - R @ centroid_Y
    
    # Apply transformation to Y
    Y_transformed = (Y_centered @ R.T) + centroid_X
    
    # Compute RMSD
    diff = X - Y_transformed
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    
    return rmsd, R, t


class SuperposedStructureSelect(Select):
    """Selector for writing superposed structure."""
    
    def __init__(self, chain_id: str):
        self.chain_id = chain_id
    
    def accept_chain(self, chain):
        return chain.id == self.chain_id


def write_superposed_pdb(structure, chain_id: str, residues: List, R: np.ndarray, 
                         t: np.ndarray, output_path: str, quiet: bool = False):
    """
    Write superposed structure to PDB file.
    Apply transformation R and t to all atoms in the selected residues.
    """
    # Get centroid of original coordinates (before transformation)
    all_coords = []
    for residue in residues:
        for atom in residue.get_atoms():
            all_coords.append(atom.get_coord())
    
    if not all_coords:
        raise ValueError("No atoms to write in superposed structure")
    
    centroid_original = np.mean(all_coords, axis=0)
    
    # Transform all atoms in the selected residues
    for residue in residues:
        for atom in residue.get_atoms():
            coord = atom.get_coord()
            # Apply transformation: (coord - centroid) @ R.T + centroid + t
            coord_centered = coord - centroid_original
            coord_transformed = coord_centered @ R.T + centroid_original + t
            atom.set_coord(coord_transformed)
    
    # Write to file
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_path, SuperposedStructureSelect(chain_id))
    
    log_message(f"Wrote superposed structure to {output_path}", quiet)


def main():
    parser = argparse.ArgumentParser(
        description='Compute RMSD between two RNA PDB structures',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('pdb1', help='First PDB file (reference)')
    parser.add_argument('pdb2', help='Second PDB file (to be aligned)')
    parser.add_argument('--selection', choices=['backbone', 'all'], default='backbone',
                        help='Atom selection: backbone (sugar-phosphate) or all atoms')
    parser.add_argument('--chain1', help='Chain ID to use from first structure')
    parser.add_argument('--chain2', help='Chain ID to use from second structure')
    parser.add_argument('--res-range', help='Residue range (1-based inclusive), e.g., 5:42')
    parser.add_argument('--output-superposed', help='Write superposed PDB2 onto PDB1')
    parser.add_argument('--strip-h', action='store_true', help='Ignore hydrogen atoms')
    parser.add_argument('--quiet', action='store_true', help='Suppress non-essential logs')
    
    args = parser.parse_args()
    
    try:
        # Parse PDB files
        parser_pdb = PDBParser(QUIET=True)
        
        log_message(f"Loading {args.pdb1}...", args.quiet)
        structure1 = parser_pdb.get_structure('struct1', args.pdb1)
        
        log_message(f"Loading {args.pdb2}...", args.quiet)
        structure2 = parser_pdb.get_structure('struct2', args.pdb2)
        
        # Select chains and filter to nucleic residues
        chain1, residues1 = select_chain(structure1, args.chain1, args.quiet)
        chain2, residues2 = select_chain(structure2, args.chain2, args.quiet)
        
        # Apply residue range if specified
        residues1 = apply_residue_range(residues1, args.res_range, args.quiet)
        residues2 = apply_residue_range(residues2, args.res_range, args.quiet)
        
        log_message(f"Selection: {args.selection}", args.quiet)
        if args.strip_h:
            log_message("Stripping hydrogen atoms", args.quiet)
        
        # Collect atom coordinate pairs
        X, Y = collect_atom_pairs(residues1, residues2, args.selection, 
                                   args.strip_h, args.quiet)
        
        # Compute RMSD using Kabsch algorithm
        rmsd, R, t = kabsch_rmsd(X, Y)
        
        # Print RMSD to stdout
        print(f"{rmsd:.6f}")
        
        # Write superposed structure if requested
        if args.output_superposed:
            write_superposed_pdb(structure2, chain2.id, residues2, R, t, 
                                args.output_superposed, args.quiet)
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
