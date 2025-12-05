#!/usr/bin/env python3
"""
FASTA Length Distribution Summary
Analyzes all FASTA files in the extracted_sequences directory
"""

import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

def parse_fasta(file_path):
    """Parse a FASTA file and return sequence lengths"""
    sequences = []
    current_seq = ""
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_seq:
                        sequences.append(len(current_seq))
                        current_seq = ""
                else:
                    current_seq += line
            
            # Don't forget the last sequence
            if current_seq:
                sequences.append(len(current_seq))
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return []
    
    return sequences

def analyze_directory(directory_path):
    """Analyze all FASTA files in a directory"""
    # Find all FASTA files
    fasta_patterns = ['*.fasta', '*.fa', '*.fas']
    fasta_files = []
    
    for pattern in fasta_patterns:
        fasta_files.extend(glob.glob(os.path.join(directory_path, pattern)))
    
    if not fasta_files:
        print(f"No FASTA files found in {directory_path}")
        return None
    
    print(f"Found {len(fasta_files)} FASTA files")
    
    # Collect all sequence lengths and file statistics
    all_lengths = []
    file_stats = []
    
    for fasta_file in sorted(fasta_files):
        lengths = parse_fasta(fasta_file)
        if lengths:
            all_lengths.extend(lengths)
            file_name = os.path.basename(fasta_file)
            file_stats.append({
                'file': file_name,
                'num_sequences': len(lengths),
                'min_length': min(lengths),
                'max_length': max(lengths),
                'mean_length': np.mean(lengths),
                'total_length': sum(lengths)
            })
            print(f"  {file_name}: {len(lengths)} sequences, lengths {min(lengths)}-{max(lengths)}")
    
    if not all_lengths:
        print("No sequences found in any FASTA files")
        return None
    
    return all_lengths, file_stats

def calculate_statistics(lengths):
    """Calculate comprehensive statistics"""
    lengths = np.array(lengths)
    
    stats = {
        'total_sequences': len(lengths),
        'min_length': np.min(lengths),
        'max_length': np.max(lengths),
        'mean_length': np.mean(lengths),
        'median_length': np.median(lengths),
        'std_length': np.std(lengths),
        'q25': np.percentile(lengths, 25),
        'q75': np.percentile(lengths, 75),
        'q10': np.percentile(lengths, 10),
        'q90': np.percentile(lengths, 90),
        'q95': np.percentile(lengths, 95),
        'q99': np.percentile(lengths, 99)
    }
    
    # Length categories
    categories = {
        'Very Short (1-10)': np.sum((lengths >= 1) & (lengths <= 10)),
        'Short (11-50)': np.sum((lengths >= 11) & (lengths <= 50)),
        'Medium (51-100)': np.sum((lengths >= 51) & (lengths <= 100)),
        'Long (101-500)': np.sum((lengths >= 101) & (lengths <= 500)),
        'Very Long (501-1000)': np.sum((lengths >= 501) & (lengths <= 1000)),
        'Extra Long (>1000)': np.sum(lengths > 1000)
    }
    
    return stats, categories

def create_visualizations(lengths, output_dir):
    """Create comprehensive visualizations"""
    plt.style.use('default')
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('FASTA Sequence Length Distribution Analysis', fontsize=16, fontweight='bold')
    
    lengths = np.array(lengths)
    
    # 1. Histogram
    axes[0, 0].hist(lengths, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
    axes[0, 0].set_xlabel('Sequence Length')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('Length Distribution (Linear Scale)')
    axes[0, 0].grid(True, alpha=0.3)
    
    # 2. Log-scale histogram
    axes[0, 1].hist(lengths, bins=50, alpha=0.7, color='lightcoral', edgecolor='black')
    axes[0, 1].set_xlabel('Sequence Length')
    axes[0, 1].set_ylabel('Frequency')
    axes[0, 1].set_title('Length Distribution (Log Scale)')
    axes[0, 1].set_yscale('log')
    axes[0, 1].grid(True, alpha=0.3)
    
    # 3. Box plot
    axes[1, 0].boxplot(lengths, vert=True, patch_artist=True,
                       boxprops=dict(facecolor='lightgreen', alpha=0.7))
    axes[1, 0].set_ylabel('Sequence Length')
    axes[1, 0].set_title('Length Distribution (Box Plot)')
    axes[1, 0].grid(True, alpha=0.3)
    
    # 4. Length categories
    categories = {
        'Very Short\n(1-10)': np.sum((lengths >= 1) & (lengths <= 10)),
        'Short\n(11-50)': np.sum((lengths >= 11) & (lengths <= 50)),
        'Medium\n(51-100)': np.sum((lengths >= 51) & (lengths <= 100)),
        'Long\n(101-500)': np.sum((lengths >= 101) & (lengths <= 500)),
        'Very Long\n(501-1000)': np.sum((lengths >= 501) & (lengths <= 1000)),
        'Extra Long\n(>1000)': np.sum(lengths > 1000)
    }
    
    cat_names = list(categories.keys())
    cat_counts = list(categories.values())
    colors = plt.cm.Set3(np.linspace(0, 1, len(cat_names)))
    
    bars = axes[1, 1].bar(cat_names, cat_counts, color=colors, alpha=0.8, edgecolor='black')
    axes[1, 1].set_ylabel('Number of Sequences')
    axes[1, 1].set_title('Sequences by Length Category')
    axes[1, 1].tick_params(axis='x', rotation=45)
    
    # Add value labels on bars
    for bar, count in zip(bars, cat_counts):
        if count > 0:
            axes[1, 1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(cat_counts)*0.01,
                           str(count), ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    
    # Save the plot
    output_file = os.path.join(output_dir, 'fasta_length_summary.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nVisualization saved to: {output_file}")
    
    return output_file

def main():
    # Target directory
    target_dir = "/Users/xiaojuzhang/Dev/pdb_prune/data/experiments_data/exp2.2_manuscript/dataset_20250731/extracted_sequences"
    
    print(f"Analyzing FASTA files in: {target_dir}")
    print("=" * 60)
    
    # Analyze the directory
    result = analyze_directory(target_dir)
    
    if result is None:
        print("No analysis could be performed.")
        return
    
    all_lengths, file_stats = result
    
    print(f"\n=== OVERALL SUMMARY ===")
    print(f"Total FASTA files analyzed: {len(file_stats)}")
    print(f"Total sequences: {len(all_lengths)}")
    
    # Calculate statistics
    stats, categories = calculate_statistics(all_lengths)
    
    print(f"\n=== LENGTH STATISTICS ===")
    print(f"Minimum length: {stats['min_length']}")
    print(f"Maximum length: {stats['max_length']}")
    print(f"Mean length: {stats['mean_length']:.2f}")
    print(f"Median length: {stats['median_length']:.2f}")
    print(f"Standard deviation: {stats['std_length']:.2f}")
    
    print(f"\n=== PERCENTILES ===")
    print(f"10th percentile: {stats['q10']:.1f}")
    print(f"25th percentile: {stats['q25']:.1f}")
    print(f"75th percentile: {stats['q75']:.1f}")
    print(f"90th percentile: {stats['q90']:.1f}")
    print(f"95th percentile: {stats['q95']:.1f}")
    print(f"99th percentile: {stats['q99']:.1f}")
    
    print(f"\n=== LENGTH CATEGORIES ===")
    total_seqs = len(all_lengths)
    for category, count in categories.items():
        percentage = (count / total_seqs) * 100
        print(f"{category}: {count} ({percentage:.1f}%)")
    
    # Create visualizations
    plot_file = create_visualizations(all_lengths, target_dir)
    
    # Save detailed statistics to CSV
    stats_df = pd.DataFrame([stats])
    stats_file = os.path.join(target_dir, 'overall_length_statistics.csv')
    stats_df.to_csv(stats_file, index=False)
    print(f"Statistics saved to: {stats_file}")
    
    # Save per-file statistics
    file_stats_df = pd.DataFrame(file_stats)
    file_stats_file = os.path.join(target_dir, 'per_file_statistics.csv')
    file_stats_df.to_csv(file_stats_file, index=False)
    print(f"Per-file statistics saved to: {file_stats_file}")
    
    # Save length categories
    categories_df = pd.DataFrame(list(categories.items()), columns=['Category', 'Count'])
    categories_df['Percentage'] = (categories_df['Count'] / total_seqs) * 100
    categories_file = os.path.join(target_dir, 'length_categories.csv')
    categories_df.to_csv(categories_file, index=False)
    print(f"Length categories saved to: {categories_file}")
    
    print(f"\n=== ANALYSIS COMPLETE ===")
    print(f"All results saved to: {target_dir}")

if __name__ == "__main__":
    main()