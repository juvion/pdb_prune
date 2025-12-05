#!/usr/bin/env python3
"""
Alternative Sequence Clustering for Lower Thresholds (<0.8)

This script performs sequence clustering using pairwise alignment when
CD-HIT-EST cannot be used (threshold < 0.8). It's slower but more flexible.

Usage:
    python alternative_seq_clustering.py --input sequences.fasta --threshold 0.6 --output clusters.csv
"""

import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO, pairwise2
from Bio.SeqUtils import gc_fraction
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from collections import Counter
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

def calculate_sequence_identity(seq1, seq2):
    """
    Calculate sequence identity between two sequences using pairwise alignment
    
    Args:
        seq1, seq2: Bio.Seq objects or strings
    
    Returns:
        Identity score (0-1)
    """
    # Convert to strings
    s1, s2 = str(seq1), str(seq2)
    
    # Quick check for identical sequences
    if s1 == s2:
        return 1.0
    
    # Quick check for very different lengths (heuristic speedup)
    len_ratio = min(len(s1), len(s2)) / max(len(s1), len(s2))
    if len_ratio < 0.5:
        return 0.0
    
    # Perform alignment (using simple match/mismatch, no gaps penalty for speed)
    alignments = pairwise2.align.globalxx(s1, s2, one_alignment_only=True)
    
    if not alignments:
        return 0.0
    
    alignment = alignments[0]
    aligned_seq1, aligned_seq2 = alignment.seqA, alignment.seqB
    
    # Calculate identity
    matches = sum(a == b for a, b in zip(aligned_seq1, aligned_seq2) if a != '-' and b != '-')
    total = max(len(s1), len(s2))
    
    return matches / total if total > 0 else 0.0

def fast_sequence_identity_matrix(sequences, seq_ids):
    """
    Calculate pairwise sequence identity matrix with optimizations
    
    Args:
        sequences: List of Bio.Seq objects
        seq_ids: List of sequence IDs
    
    Returns:
        numpy array of identity scores
    """
    n = len(sequences)
    identity_matrix = np.eye(n)  # Diagonal is 1.0
    
    print(f"Calculating pairwise identities for {n} sequences...")
    total_pairs = n * (n - 1) // 2
    
    with tqdm(total=total_pairs, desc="Calculating identities") as pbar:
        for i in range(n):
            for j in range(i + 1, n):
                identity = calculate_sequence_identity(sequences[i], sequences[j])
                identity_matrix[i, j] = identity
                identity_matrix[j, i] = identity
                pbar.update(1)
    
    return identity_matrix

def cluster_by_threshold(identity_matrix, seq_ids, threshold):
    """
    Cluster sequences based on identity threshold using agglomerative clustering
    
    Args:
        identity_matrix: N×N matrix of sequence identities
        seq_ids: List of sequence IDs
        threshold: Identity threshold (0-1)
    
    Returns:
        Dictionary mapping sequence IDs to cluster IDs
    """
    # Convert similarity to distance
    distance_matrix = 1 - identity_matrix
    np.fill_diagonal(distance_matrix, 0)
    
    # Convert to condensed distance matrix
    condensed_dist = squareform(distance_matrix, checks=False)
    
    # Perform hierarchical clustering
    print(f"Performing hierarchical clustering at threshold {threshold}...")
    linkage_matrix = linkage(condensed_dist, method='average')
    
    # Cut tree at distance threshold
    distance_threshold = 1 - threshold
    cluster_labels = fcluster(linkage_matrix, distance_threshold, criterion='distance')
    
    # Create mapping
    clusters = {seq_ids[i]: int(cluster_labels[i]) for i in range(len(seq_ids))}
    
    return clusters

def greedy_clustering(identity_matrix, seq_ids, threshold):
    """
    Fast greedy clustering algorithm (similar to CD-HIT's approach)
    
    Args:
        identity_matrix: N×N matrix of sequence identities
        seq_ids: List of sequence IDs
        threshold: Identity threshold (0-1)
    
    Returns:
        Dictionary mapping sequence IDs to cluster IDs
    """
    n = len(seq_ids)
    
    # Sort sequences by length (longest first, like CD-HIT)
    seq_lengths = [len(str(seq)) for seq in sequences]
    sorted_indices = sorted(range(n), key=lambda i: seq_lengths[i], reverse=True)
    
    clusters = {}
    cluster_id = 0
    representatives = []
    
    print(f"Performing greedy clustering at threshold {threshold}...")
    
    for idx in tqdm(sorted_indices, desc="Clustering"):
        seq_id = seq_ids[idx]
        
        # Check if this sequence is similar to any representative
        assigned = False
        for rep_idx in representatives:
            if identity_matrix[idx, rep_idx] >= threshold:
                # Assign to this cluster
                clusters[seq_id] = clusters[seq_ids[rep_idx]]
                assigned = True
                break
        
        if not assigned:
            # Create new cluster
            clusters[seq_id] = cluster_id
            representatives.append(idx)
            cluster_id += 1
    
    return clusters

def main():
    parser = argparse.ArgumentParser(
        description='Alternative sequence clustering for thresholds < 0.8'
    )
    parser.add_argument(
        '--input',
        required=True,
        help='Input FASTA file'
    )
    parser.add_argument(
        '--threshold',
        type=float,
        required=True,
        help='Sequence identity threshold (0-1)'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output CSV file with cluster assignments'
    )
    parser.add_argument(
        '--method',
        choices=['greedy', 'hierarchical'],
        default='greedy',
        help='Clustering method (default: greedy, faster)'
    )
    parser.add_argument(
        '--dist',
        help='Output file for cluster size distribution'
    )
    
    args = parser.parse_args()
    
    print("="*70)
    print("Alternative Sequence Clustering")
    print("="*70)
    print(f"Input:      {args.input}")
    print(f"Threshold:  {args.threshold}")
    print(f"Method:     {args.method}")
    print(f"Output:     {args.output}")
    print("="*70)
    print()
    
    # Load sequences
    print("Loading sequences...")
    global sequences  # Make available to greedy_clustering
    sequences = []
    seq_ids = []
    
    for record in SeqIO.parse(args.input, "fasta"):
        sequences.append(record.seq)
        seq_ids.append(record.id)
    
    print(f"Loaded {len(sequences)} sequences")
    print(f"Length range: {min(len(s) for s in sequences)} - {max(len(s) for s in sequences)} nt")
    print()
    
    # Calculate identity matrix
    identity_matrix = fast_sequence_identity_matrix(sequences, seq_ids)
    
    # Perform clustering
    if args.method == 'greedy':
        clusters = greedy_clustering(identity_matrix, seq_ids, args.threshold)
    else:
        clusters = cluster_by_threshold(identity_matrix, seq_ids, args.threshold)
    
    # Analyze results
    cluster_sizes = Counter(clusters.values())
    n_clusters = len(cluster_sizes)
    
    print()
    print("="*70)
    print("Clustering Results:")
    print(f"  Total sequences:     {len(clusters)}")
    print(f"  Total clusters:      {n_clusters}")
    print(f"  Largest cluster:     {max(cluster_sizes.values())} sequences")
    print(f"  Smallest cluster:    {min(cluster_sizes.values())} sequence(s)")
    print(f"  Average cluster size: {np.mean(list(cluster_sizes.values())):.2f}")
    print(f"  Singleton clusters:  {sum(1 for s in cluster_sizes.values() if s == 1)}")
    print("="*70)
    
    # Save results
    df = pd.DataFrame([
        {'sequence_id': seq_id, 'cluster_id': cluster_id}
        for seq_id, cluster_id in clusters.items()
    ])
    df = df.sort_values(['cluster_id', 'sequence_id']).reset_index(drop=True)
    df.to_csv(args.output, index=False)
    print(f"\n✓ Cluster assignments saved to {args.output}")
    
    # Save distribution if requested
    if args.dist:
        size_dist = Counter(cluster_sizes.values())
        df_dist = pd.DataFrame([
            {'cluster_size': size, 'num_clusters': count}
            for size, count in sorted(size_dist.items())
        ])
        df_dist.to_csv(args.dist, index=False)
        print(f"✓ Cluster size distribution saved to {args.dist}")

if __name__ == "__main__":
    main()
