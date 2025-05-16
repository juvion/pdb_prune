#!/usr/bin/env python3

import os
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import pandas as pd

def analyze_sequence_lengths(fasta_dir: str = "rna_sequences", bin_size: int = 10):
    """
    Analyze the length distribution of RNA sequences from FASTA files.
    
    Args:
        fasta_dir (str): Directory containing FASTA files
        bin_size (int): Size of bins for length distribution (default: 10)
    """
    # Get all FASTA files
    fasta_files = list(Path(fasta_dir).glob("*.fasta"))
    if not fasta_files:
        print(f"No FASTA files found in {fasta_dir}")
        return
    
    # Collect sequence lengths
    lengths = []
    for fasta_file in fasta_files:
        for record in SeqIO.parse(fasta_file, "fasta"):
            lengths.append(len(record.seq))
    
    if not lengths:
        print("No sequences found in FASTA files")
        return
    
    # Calculate statistics
    lengths = np.array(lengths)
    stats = {
        "min": np.min(lengths),
        "max": np.max(lengths),
        "mean": np.mean(lengths),
        "median": np.median(lengths),
        "std": np.std(lengths)
    }
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Calculate histogram with specified bin size
    max_length = np.ceil(stats['max'] / bin_size) * bin_size
    bins = np.arange(0, max_length + bin_size, bin_size)
    hist, bin_edges = np.histogram(lengths, bins=bins)
    
    # Create distribution table
    distribution_df = pd.DataFrame({
        'Length_Range': [f"{int(bin_edges[i])}-{int(bin_edges[i+1])}" for i in range(len(bin_edges)-1)],
        'Count': hist,
        'Percentage': hist / len(lengths) * 100
    })
    
    # Save distribution table to CSV
    distribution_df.to_csv('sequence_length_distribution.csv', index=False)
    
    # Plot histogram
    ax1.hist(lengths, bins=bins, edgecolor='black')
    ax1.set_title(f'Distribution of RNA Sequence Lengths (Bin Size: {bin_size})')
    ax1.set_xlabel('Sequence Length')
    ax1.set_ylabel('Count')
    ax1.grid(True, alpha=0.3)
    
    # Add statistics as text
    stats_text = (
        f"Statistics:\n"
        f"Number of sequences: {len(lengths)}\n"
        f"Min length: {stats['min']:.1f}\n"
        f"Max length: {stats['max']:.1f}\n"
        f"Mean length: {stats['mean']:.1f}\n"
        f"Median length: {stats['median']:.1f}\n"
        f"Std dev: {stats['std']:.1f}"
    )
    ax1.text(0.95, 0.95, stats_text,
             transform=ax1.transAxes,
             verticalalignment='top',
             horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Box plot
    ax2.boxplot(lengths, vert=False)
    ax2.set_title('Box Plot of RNA Sequence Lengths')
    ax2.set_xlabel('Sequence Length')
    ax2.grid(True, alpha=0.3)
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig('sequence_length_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Print statistics
    print("\nSequence Length Statistics:")
    print(f"Number of sequences: {len(lengths)}")
    print(f"Min length: {stats['min']:.1f}")
    print(f"Max length: {stats['max']:.1f}")
    print(f"Mean length: {stats['mean']:.1f}")
    print(f"Median length: {stats['median']:.1f}")
    print(f"Standard deviation: {stats['std']:.1f}")
    print("\nPlot saved as 'sequence_length_distribution.png'")
    print("Distribution table saved as 'sequence_length_distribution.csv'")
    
    # Print distribution table
    print("\nLength Distribution Table:")
    print(distribution_df.to_string(index=False))

if __name__ == "__main__":
    analyze_sequence_lengths(bin_size=10) 