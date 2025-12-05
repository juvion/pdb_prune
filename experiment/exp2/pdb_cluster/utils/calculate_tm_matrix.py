#!/usr/bin/env python3
"""
Calculate pairwise TM-score matrix using US-align
"""

import os
import subprocess
import numpy as np
import pandas as pd
from tqdm import tqdm
import argparse
import time

def calculate_tm_score(pdb1, pdb2, timeout=60):
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
        result = subprocess.run(
            ['USalign', pdb1, pdb2],
            capture_output=True,
            text=True,
            timeout=timeout
        )
        
        # Parse TM-score from output
        # Look for line like: "TM-score= 0.xxxxx (normalized by length of Chain_1)"
        for line in result.stdout.split('\n'):
            if 'TM-score=' in line and 'Chain_1' in line:
                tm_score = float(line.split('TM-score=')[1].split()[0])
                return tm_score
        
        # If no TM-score found, return 0
        return 0.0
    
    except subprocess.TimeoutExpired:
        print(f"⚠ Timeout for {os.path.basename(pdb1)} vs {os.path.basename(pdb2)}")
        return 0.0
    except Exception as e:
        print(f"✗ Error calculating TM-score: {str(e)}")
        return 0.0

def calculate_tm_matrix(pdb_dir, output_file, checkpoint_file=None, resume=False):
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
    # Get list of PDB files
    pdb_files = sorted([f for f in os.listdir(pdb_dir) if f.endswith('.pdb')])
    structure_ids = [f.replace('.pdb', '') for f in pdb_files]
    n = len(pdb_files)
    
    print(f"Found {n} PDB structures")
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
    completed = start_i * (n - start_i) // 2 + (start_j - start_i - 1) if start_j > 0 else 0
    
    print(f"Starting from ({start_i}, {start_j}), {completed} already completed")
    
    with tqdm(total=total_pairs, initial=completed) as pbar:
        for i in range(start_i, n):
            # Set diagonal
            tm_matrix[i, i] = 1.0
            
            j_start = start_j if i == start_i else i + 1
            
            for j in range(j_start, n):
                pdb1 = os.path.join(pdb_dir, pdb_files[i])
                pdb2 = os.path.join(pdb_dir, pdb_files[j])
                
                tm_score = calculate_tm_score(pdb1, pdb2)
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
    
    args = parser.parse_args()
    
    print("="*60)
    print("TM-score Matrix Calculator")
    print("="*60)
    
    if args.validate:
        validate_tm_matrix(args.output)
    else:
        start_time = time.time()
        
        tm_matrix, ids = calculate_tm_matrix(
            args.pdb_dir,
            args.output,
            args.checkpoint,
            args.resume
        )
        
        elapsed = time.time() - start_time
        print(f"\nTotal time: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")
    
    print("="*60)
