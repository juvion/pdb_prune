#!/usr/bin/env python3
"""
Perform agglomerative clustering on TM-score similarity matrix
"""

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
import argparse
from collections import Counter
import sys

def perform_structure_clustering(tm_matrix_file, thresholds=[0.6, 0.5, 0.4], 
                                 method='average', output_prefix='structure_clusters'):
    """
    Perform agglomerative clustering on TM-score matrix
    
    Args:
        tm_matrix_file: Path to TM-score matrix CSV
        thresholds: List of TM-score thresholds for clustering
        method: Linkage method ('single', 'complete', 'average', 'ward')
        output_prefix: Prefix for output files
    
    Returns:
        results: Dictionary of cluster assignments for each threshold
        linkage_matrix: Hierarchical linkage matrix
    """
    print(f"Loading TM-score matrix from: {tm_matrix_file}")
    
    # Load TM-score matrix
    df = pd.read_csv(tm_matrix_file, index_col=0)
    structure_ids = df.index.tolist()
    tm_matrix = df.values
    n = len(structure_ids)
    
    print(f"Matrix size: {n} x {n}")
    print(f"Linkage method: {method}")
    
    # Convert similarity to distance (1 - TM-score)
    distance_matrix = 1 - tm_matrix
    
    # Ensure diagonal is zero (numerical precision)
    np.fill_diagonal(distance_matrix, 0)
    
    # Convert to condensed distance matrix for scipy
    condensed_dist = squareform(distance_matrix, checks=False)
    
    # Perform hierarchical clustering
    print(f"\nPerforming agglomerative clustering...")
    linkage_matrix = linkage(condensed_dist, method=method)
    
    print(f"Linkage matrix shape: {linkage_matrix.shape}")
    
    # Cluster at different thresholds
    results = {}
    
    print(f"\nClustering at different thresholds:")
    print("-" * 60)
    
    for threshold in thresholds:
        # Convert TM-score threshold to distance threshold
        distance_threshold = 1 - threshold
        
        # Get cluster labels
        cluster_labels = fcluster(linkage_matrix, distance_threshold, criterion='distance')
        
        # Create mapping
        cluster_dict = {
            structure_ids[i]: int(cluster_labels[i])
            for i in range(len(structure_ids))
        }
        
        results[threshold] = cluster_dict
        
        # Analyze clusters
        cluster_sizes = Counter(cluster_labels)
        n_clusters = len(cluster_sizes)
        
        # Save results
        output_file = f"{output_prefix}_{threshold}.csv"
        df_out = pd.DataFrame([
            {'structure_id': sid, 'cluster_id': cid}
            for sid, cid in cluster_dict.items()
        ])
        df_out = df_out.sort_values(['cluster_id', 'structure_id']).reset_index(drop=True)
        df_out.to_csv(output_file, index=False)
        
        print(f"Threshold {threshold} (distance {distance_threshold:.2f}):")
        print(f"  Clusters: {n_clusters}")
        print(f"  Largest cluster: {max(cluster_sizes.values())} structures")
        print(f"  Smallest cluster: {min(cluster_sizes.values())} structure(s)")
        print(f"  Average cluster size: {np.mean(list(cluster_sizes.values())):.2f}")
        print(f"  Singleton clusters: {sum(1 for s in cluster_sizes.values() if s == 1)}")
        print(f"  Output: {output_file}")
        
        # Save cluster size distribution
        dist_file = f"{output_prefix}_{threshold}_distribution.csv"
        size_dist = Counter(cluster_sizes.values())
        df_dist = pd.DataFrame([
            {'cluster_size': size, 'num_clusters': count}
            for size, count in sorted(size_dist.items())
        ])
        df_dist.to_csv(dist_file, index=False)
        print(f"  Distribution: {dist_file}")
        print()
    
    return results, linkage_matrix, structure_ids

def plot_dendrogram(linkage_matrix, structure_ids, output_file='clustering_dendrogram.png',
                   figsize=(14, 8), max_d=None):
    """
    Plot hierarchical clustering dendrogram
    
    Args:
        linkage_matrix: Hierarchical linkage matrix
        structure_ids: List of structure IDs
        output_file: Output PNG file
        figsize: Figure size (width, height)
        max_d: Maximum distance to draw horizontal line
    """
    plt.figure(figsize=figsize)

    n = len(structure_ids)
    # Truncate labels if too many
    labels = structure_ids if n <= 50 else None

    # Increase recursion limit to handle large trees (SciPy dendrogram is recursive)
    try:
        if sys.getrecursionlimit() < 10000:
            sys.setrecursionlimit(10000)
    except Exception:
        pass

    # For very large N, use truncated dendrogram to avoid deep recursion
    truncate_kwargs = {}
    if n > 1000:
        truncate_kwargs = {
            'truncate_mode': 'level',  # show only the last p merged levels
            'p': 50,                   # adjust as needed
            'show_leaf_counts': True,
            'no_labels': True,
        }

    try:
        dendrogram(
            linkage_matrix,
            labels=labels,
            leaf_font_size=8 if labels else None,
            no_labels=(labels is None),
            **truncate_kwargs
        )
    except RecursionError:
        # Fallback: further truncate if recursion persists
        print("⚠ RecursionError while plotting dendrogram; falling back to aggressive truncation (p=20)")
        plt.clf()
        plt.figure(figsize=figsize)
        dendrogram(
            linkage_matrix,
            truncate_mode='level',
            p=20,
            show_leaf_counts=True,
            no_labels=True
        )
    
    plt.title('Hierarchical Clustering Dendrogram (Structure Similarity)', fontsize=14)
    plt.xlabel('Structure Index' if labels is None else 'Structure ID', fontsize=12)
    plt.ylabel('Distance (1 - TM-score)', fontsize=12)
    
    # Draw threshold line if specified
    if max_d:
        plt.axhline(y=max_d, c='red', linestyle='--', label=f'Threshold: {max_d:.2f}')
        plt.legend()
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Dendrogram saved to {output_file}")
    plt.close()

def plot_cluster_size_distribution(results, output_file='cluster_size_distribution.png'):
    """
    Plot cluster size distributions for different thresholds
    
    Args:
        results: Dictionary of cluster assignments for each threshold
        output_file: Output PNG file
    """
    fig, axes = plt.subplots(1, len(results), figsize=(5*len(results), 4))
    
    if len(results) == 1:
        axes = [axes]
    
    for idx, (threshold, cluster_dict) in enumerate(sorted(results.items())):
        cluster_sizes = Counter(cluster_dict.values())
        size_counts = Counter(cluster_sizes.values())
        
        sizes = sorted(size_counts.keys())
        counts = [size_counts[s] for s in sizes]
        
        axes[idx].bar(sizes, counts, alpha=0.7, edgecolor='black')
        axes[idx].set_xlabel('Cluster Size', fontsize=10)
        axes[idx].set_ylabel('Number of Clusters', fontsize=10)
        axes[idx].set_title(f'Threshold: {threshold}\n({len(cluster_sizes)} clusters)', fontsize=11)
        axes[idx].grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Cluster size distribution plot saved to {output_file}")
    plt.close()

def compare_thresholds(results, output_file='threshold_comparison.csv'):
    """
    Compare clustering results across thresholds
    
    Args:
        results: Dictionary of cluster assignments for each threshold
        output_file: Output CSV file
    """
    comparison = []
    
    for threshold in sorted(results.keys()):
        cluster_dict = results[threshold]
        cluster_sizes = Counter(cluster_dict.values())
        
        comparison.append({
            'threshold': threshold,
            'distance_threshold': 1 - threshold,
            'num_clusters': len(cluster_sizes),
            'num_structures': len(cluster_dict),
            'max_cluster_size': max(cluster_sizes.values()),
            'min_cluster_size': min(cluster_sizes.values()),
            'avg_cluster_size': np.mean(list(cluster_sizes.values())),
            'median_cluster_size': np.median(list(cluster_sizes.values())),
            'singleton_clusters': sum(1 for s in cluster_sizes.values() if s == 1)
        })
    
    df = pd.DataFrame(comparison)
    df.to_csv(output_file, index=False)
    print(f"✓ Threshold comparison saved to {output_file}")
    
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Agglomerative clustering on TM-score matrix'
    )
    parser.add_argument(
        '--input',
        required=True,
        help='Input TM-score matrix CSV file'
    )
    parser.add_argument(
        '--thresholds',
        nargs='+',
        type=float,
        default=[0.6, 0.5, 0.4],
        help='TM-score thresholds for clustering (default: 0.6 0.5 0.4)'
    )
    parser.add_argument(
        '--method',
        default='average',
        choices=['single', 'complete', 'average', 'ward'],
        help='Linkage method (default: average)'
    )
    parser.add_argument(
        '--output_prefix',
        default='structure_clusters',
        help='Prefix for output files (default: structure_clusters)'
    )
    parser.add_argument(
        '--plot',
        action='store_true',
        help='Generate plots (dendrogram and distributions)'
    )
    
    args = parser.parse_args()
    
    print("="*60)
    print("Agglomerative Clustering on Structure Similarity")
    print("="*60)
    
    # Perform clustering
    results, linkage_mat, struct_ids = perform_structure_clustering(
        args.input,
        args.thresholds,
        args.method,
        args.output_prefix
    )
    
    # Compare thresholds
    print("\nGenerating threshold comparison...")
    comparison_df = compare_thresholds(results, f"{args.output_prefix}_comparison.csv")
    print("\nThreshold Comparison:")
    print(comparison_df.to_string(index=False))
    
    # Generate plots if requested
    if args.plot:
        print("\nGenerating plots...")
        plot_dendrogram(linkage_mat, struct_ids, f"{args.output_prefix}_dendrogram.png")
        plot_cluster_size_distribution(results, f"{args.output_prefix}_size_distribution.png")
    
    print("\n" + "="*60)
    print("Clustering completed successfully!")
    print("="*60)
