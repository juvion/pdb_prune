#!/usr/bin/env python3
"""
Calculate pairwise TM-score matrix using US-align
"""

import os
import re
import shutil
import subprocess
import numpy as np
import pandas as pd
from tqdm import tqdm
import argparse
import time

def _parse_tm_scores(stdout: str) -> float:
    """
    Robustly parse TM-score values from US-align stdout.

    Returns the maximum TM-score found among the reported normalizations.
    If none found, returns 0.0.
    """
    tm_scores = []
    # Typical line: "TM-score= 0.5123 (normalized by length of Chain_1)"
    pattern = re.compile(r"TM-score=\s*([0-9]*\.?[0-9]+)")
    for line in stdout.splitlines():
        if 'TM-score' in line:
            m = pattern.search(line)
            if m:
                try:
                    tm_scores.append(float(m.group(1)))
                except ValueError:
                    continue
    return max(tm_scores) if tm_scores else 0.0


def _normalize_mol_arg(mol_type, verbose=False):
    """
    Normalize user-provided molecule type to what US-align expects.

    Accepts strings like 'auto', 'prot'/'protein', 'RNA', 'DNA', or legacy
    numeric values where 1=auto, 2=prot, 3=RNA.
    Returns one of: 'auto', 'prot', 'RNA'.
    """
    if mol_type is None:
        return 'RNA'

    # Handle integers or numeric strings
    try:
        val = int(mol_type)
        if val == 1:
            return 'auto'
        elif val == 2:
            return 'prot'
        elif val == 3:
            return 'RNA'
        else:
            if verbose:
                print(f"[DEBUG] Unknown numeric --mol '{mol_type}', defaulting to RNA")
            return 'RNA'
    except (ValueError, TypeError):
        pass

    # Handle string values and synonyms
    s = str(mol_type).strip().lower()
    if s in {'auto', 'a'}:
        return 'auto'
    if s in {'prot', 'protein', 'p'}:
        return 'prot'
    if s in {'rna', 'dna', 'nucleic', 'n', 'r', 'd'}:
        return 'RNA'

    if verbose:
        print(f"[DEBUG] Unrecognized --mol '{mol_type}', defaulting to RNA")
    return 'RNA'


def calculate_tm_score(pdb1, pdb2, timeout=60, mol_type='RNA', usalign_bin='USalign', verbose=False):
    """
    Calculate TM-score between two PDB structures using US-align
    
    Args:
        pdb1: Path to first PDB file
        pdb2: Path to second PDB file
        timeout: Maximum time to wait for US-align (seconds)
    
    Returns:
        TM-score (float)
    """
    try:
        mol_norm = _normalize_mol_arg(mol_type, verbose=verbose)
        cmd = [usalign_bin, '-mol', mol_norm, pdb1, pdb2]
        if verbose:
            print(f"[DEBUG] Running: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout
        )

        # Parse TM-score from output (robust to format changes)
        tm_score = _parse_tm_scores(result.stdout)
        if verbose and (result.returncode != 0 or tm_score == 0.0):
            print(f"[DEBUG] US-align returncode={result.returncode} tm={tm_score} for {os.path.basename(pdb1)} vs {os.path.basename(pdb2)}")
            if result.stderr:
                print(f"[DEBUG] stderr: {result.stderr.strip()[:200]}")
            if result.stdout:
                head = '\n'.join(result.stdout.splitlines()[:5])
                print(f"[DEBUG] stdout head:\n{head}")
        return tm_score
    
    except subprocess.TimeoutExpired:
        print(f"⚠ Timeout for {os.path.basename(pdb1)} vs {os.path.basename(pdb2)}")
        return 0.0
    except Exception as e:
        if verbose:
            print(f"✗ Error calculating TM-score: {str(e)}")
        return 0.0

def calculate_tm_matrix(pdb_dir, output_file, checkpoint_file=None, resume=False, limit=None, timeout=60, mol_type='RNA', usalign_bin='USalign', verbose=False):
    """
    Calculate pairwise TM-score matrix for all PDB files
    
    Args:
        pdb_dir: Directory containing PDB files
        output_file: Output file for TM-score matrix (CSV)
        checkpoint_file: File to save progress for resuming
        resume: Whether to resume from checkpoint
    
    Returns:
        tm_matrix: numpy array of TM-scores
        structure_ids: list of structure IDs
    """
    # Check US-align availability up-front
    usalign_path = shutil.which(usalign_bin)
    if usalign_path is None:
        raise RuntimeError(f"US-align executable '{usalign_bin}' not found in PATH. Install USalign or provide --usalign_bin.")

    # Get list of PDB files
    pdb_files = sorted([f for f in os.listdir(pdb_dir) if f.endswith('.pdb')])
    if limit is not None:
        pdb_files = pdb_files[:int(limit)]
    structure_ids = [f.replace('.pdb', '') for f in pdb_files]
    n = len(pdb_files)
    
    print(f"Found {n} PDB structures" + (f" (limited to first {limit})" if limit else ""))
    print(f"Total pairwise comparisons: {n * (n - 1) // 2}")
    
    # Initialize matrix
    tm_matrix = np.zeros((n, n))
    
    # Check for checkpoint
    start_i, start_j = 0, 0
    if resume and checkpoint_file and os.path.exists(checkpoint_file):
        print(f"Resuming from checkpoint: {checkpoint_file}")
        checkpoint_df = pd.read_csv(checkpoint_file, index_col=0)
        tm_matrix = checkpoint_df.values
        
        # Find where to resume
        for i in range(n):
            for j in range(i + 1, n):
                if tm_matrix[i, j] == 0 and i != j:
                    start_i, start_j = i, j
                    break
            if start_j > 0:
                break
    
    # Calculate pairwise TM-scores
    total_pairs = n * (n - 1) // 2
    if resume and checkpoint_file and os.path.exists(checkpoint_file):
        # Count completed pairs strictly above diagonal
        completed = int(np.sum(np.triu(tm_matrix, k=1) > 0))
    else:
        completed = 0

    print(f"Starting from ({start_i}, {start_j}), {completed} already completed")
    
    with tqdm(total=total_pairs, initial=completed) as pbar:
        for i in range(start_i, n):
            # Set diagonal
            tm_matrix[i, i] = 1.0

            # Only compute above the diagonal; skip self-pairs
            j_start = (start_j if (resume and i == start_i and start_j > i) else i + 1)

            for j in range(j_start, n):
                pdb1 = os.path.join(pdb_dir, pdb_files[i])
                pdb2 = os.path.join(pdb_dir, pdb_files[j])
                
                tm_score = calculate_tm_score(pdb1, pdb2, timeout=timeout, mol_type=mol_type, usalign_bin=usalign_bin, verbose=verbose)
                tm_matrix[i, j] = tm_score
                tm_matrix[j, i] = tm_score  # Symmetric
                
                pbar.update(1)
                
                # Save checkpoint every 100 comparisons
                if pbar.n % 100 == 0 and checkpoint_file:
                    df_checkpoint = pd.DataFrame(
                        tm_matrix, 
                        index=structure_ids, 
                        columns=structure_ids
                    )
                    df_checkpoint.to_csv(checkpoint_file)
    
    # Save final matrix
    df = pd.DataFrame(tm_matrix, index=structure_ids, columns=structure_ids)
    df.to_csv(output_file)
    print(f"\n✓ TM-score matrix saved to {output_file}")
    
    # Print statistics
    non_diag_scores = tm_matrix[np.triu_indices(n, k=1)]
    print(f"\nTM-score Statistics:")
    print(f"  Mean: {non_diag_scores.mean():.4f}")
    print(f"  Median: {np.median(non_diag_scores):.4f}")
    print(f"  Min: {non_diag_scores.min():.4f}")
    print(f"  Max: {non_diag_scores.max():.4f}")
    print(f"  Std: {non_diag_scores.std():.4f}")
    
    return tm_matrix, structure_ids

def validate_tm_matrix(tm_matrix_file):
    """
    Validate TM-score matrix
    
    Args:
        tm_matrix_file: Path to TM-score matrix CSV
    """
    df = pd.read_csv(tm_matrix_file, index_col=0)
    matrix = df.values
    n = len(matrix)
    
    print("Validating TM-score matrix...")
    
    # Check symmetry
    is_symmetric = np.allclose(matrix, matrix.T)
    print(f"  Symmetric: {is_symmetric}")
    
    # Check diagonal
    diagonal_ones = np.allclose(np.diag(matrix), 1.0)
    print(f"  Diagonal all 1.0: {diagonal_ones}")
    
    # Check range [0, 1]
    in_range = np.all((matrix >= 0) & (matrix <= 1))
    print(f"  All values in [0, 1]: {in_range}")
    
    # Check for missing values
    has_missing = np.any(np.isnan(matrix))
    print(f"  Has missing values: {has_missing}")
    
    if is_symmetric and diagonal_ones and in_range and not has_missing:
        print("\n✓ Matrix validation passed!")
    else:
        print("\n✗ Matrix validation failed!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Calculate TM-score matrix using US-align'
    )
    parser.add_argument(
        '--pdb_dir', 
        required=True,
        help='Directory containing PDB files'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output CSV file for TM-score matrix'
    )
    parser.add_argument(
        '--checkpoint',
        help='Checkpoint file to save progress'
    )
    parser.add_argument(
        '--resume',
        action='store_true',
        help='Resume from checkpoint'
    )
    parser.add_argument(
        '--validate',
        action='store_true',
        help='Validate existing TM-score matrix'
    )
    parser.add_argument(
        '--limit',
        type=int,
        help='Limit to first N PDBs for quick testing'
    )
    parser.add_argument(
        '--timeout',
        type=int,
        default=60,
        help='US-align timeout per comparison (seconds)'
    )
    parser.add_argument(
        '--mol',
        type=str,
        default='RNA',
        help="US-align molecule type: 'auto' | 'prot' | 'protein' | 'RNA' | 'DNA'. Legacy numeric values supported: 1=auto, 2=prot, 3=RNA"
    )
    parser.add_argument(
        '--usalign_bin',
        type=str,
        default='USalign',
        help='Path or name of US-align executable'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Print debug information when TM-score parsing fails or US-align errors'
    )
    
    args = parser.parse_args()
    
    print("="*60)
    print("TM-score Matrix Calculator")
    print("="*60)
    
    if args.validate:
        validate_tm_matrix(args.output)
    else:
        start_time = time.time()
        
        try:
            tm_matrix, ids = calculate_tm_matrix(
                args.pdb_dir,
                args.output,
                args.checkpoint,
                args.resume,
                limit=args.limit,
                timeout=args.timeout,
                mol_type=args.mol,
                usalign_bin=args.usalign_bin,
                verbose=args.verbose,
            )
        except RuntimeError as e:
            print(f"\n✗ {e}\n")
            print("Hint: Install USalign and ensure it's on PATH, or pass --usalign_bin /path/to/USalign")
            raise
        
        elapsed = time.time() - start_time
        print(f"\nTotal time: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")
    
    print("="*60)
