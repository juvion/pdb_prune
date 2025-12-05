#!/usr/bin/env python3
"""
Parse CD-HIT clustering results from .clstr files
"""

import pandas as pd
import argparse

def parse_cdhit_clusters(clstr_file):
    """
    Parse CD-HIT cluster file (.clstr)
    
    Args:
        clstr_file: Path to .clstr file
    
    Returns:
        Dictionary mapping sequence IDs to cluster IDs
    """
    clusters = {}
    current_cluster = -1
    
    print(f"Parsing cluster file: {clstr_file}")
    
    with open(clstr_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>Cluster'):
                current_cluster = int(line.split()[1])
            elif line:
                # Extract sequence ID
                # Format: "0	500aa, >sequence_name... at 100.00%"
                try:
                    seq_id = line.split('>')[1].split('...')[0]
                    clusters[seq_id] = current_cluster
                except IndexError:
                    continue
    
    return clusters

def analyze_clusters(clusters):
    """
    Analyze cluster statistics
    
    Args:
        clusters: Dictionary mapping sequence IDs to cluster IDs
    
    Returns:
        Dictionary with statistics
    """
    from collections import Counter
    
    cluster_sizes = Counter(clusters.values())
    
    stats = {
        'total_sequences': len(clusters),
        'total_clusters': len(set(clusters.values())),
        'largest_cluster': max(cluster_sizes.values()),
        'smallest_cluster': min(cluster_sizes.values()),
        'average_cluster_size': sum(cluster_sizes.values()) / len(cluster_sizes),
        'singleton_clusters': sum(1 for size in cluster_sizes.values() if size == 1)
    }
    
    return stats

def save_cluster_mapping(clusters, output_file):
    """
    Save cluster mapping to CSV file
    
    Args:
        clusters: Dictionary mapping sequence IDs to cluster IDs
        output_file: Output CSV file path
    """
    df = pd.DataFrame([
        {'sequence_id': seq_id, 'cluster_id': cluster_id}
        for seq_id, cluster_id in clusters.items()
    ])
    
    # Sort by cluster_id then sequence_id
    df = df.sort_values(['cluster_id', 'sequence_id']).reset_index(drop=True)
    
    df.to_csv(output_file, index=False)
    print(f"\n✓ Saved cluster mapping to {output_file}")
    
    # Print statistics
    stats = analyze_clusters(clusters)
    print(f"\nCluster Statistics:")
    print(f"  Total sequences: {stats['total_sequences']}")
    print(f"  Total clusters: {stats['total_clusters']}")
    print(f"  Largest cluster: {stats['largest_cluster']} sequences")
    print(f"  Smallest cluster: {stats['smallest_cluster']} sequence(s)")
    print(f"  Average cluster size: {stats['average_cluster_size']:.2f}")
    print(f"  Singleton clusters: {stats['singleton_clusters']}")
    
    return df

def save_cluster_distribution(clusters, output_file):
    """
    Save cluster size distribution
    
    Args:
        clusters: Dictionary mapping sequence IDs to cluster IDs
        output_file: Output file path
    """
    from collections import Counter
    
    cluster_sizes = Counter(clusters.values())
    size_dist = Counter(cluster_sizes.values())
    
    df = pd.DataFrame([
        {'cluster_size': size, 'num_clusters': count}
        for size, count in sorted(size_dist.items())
    ])
    
    df.to_csv(output_file, index=False)
    print(f"✓ Saved cluster size distribution to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse CD-HIT cluster results')
    parser.add_argument('--input', required=True, help='Input .clstr file from CD-HIT')
    parser.add_argument('--output', required=True, help='Output CSV file')
    parser.add_argument('--dist', help='Output file for cluster size distribution')
    
    args = parser.parse_args()
    
    print("="*60)
    print("CD-HIT Cluster Parser")
    print("="*60)
    
    # Parse clusters
    clusters = parse_cdhit_clusters(args.input)
    
    # Save mapping
    save_cluster_mapping(clusters, args.output)
    
    # Save distribution if requested
    if args.dist:
        save_cluster_distribution(clusters, args.dist)
    
    print("="*60)
