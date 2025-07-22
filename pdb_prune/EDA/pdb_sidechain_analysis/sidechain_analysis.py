"""
PDB Side Chain Analysis Script

Analyzes side chains in .ent PDB files, outputs statistics and visualizations.
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
from tqdm import tqdm

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

plt.style.use('seaborn')
sns.set_palette("husl")

# List of standard RNA side chain atom names (for nucleotides)
# For RNA, side chains are the base atoms (not the backbone: P, O5', C5', C4', C3', O3', etc.)
# We'll count all atoms not in the backbone as side chain atoms
BACKBONE_ATOMS = {'P', "O5'", "C5'", "C4'", "C3'", "O3'", "O4'", "C2'", "O2'", "C1'"}


def list_pdb_files(pdb_dir: str):
    """
    List all .ent PDB files in the given directory.
    """
    pdb_path = Path(pdb_dir)
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB directory not found: {pdb_dir}")
    pdb_files = list(pdb_path.glob("*.ent"))
    logger.info(f"Found {len(pdb_files)} PDB files in {pdb_dir}")
    return pdb_files

def count_side_chains_in_pdb(pdb_file: Path):
    """
    Count side chain atoms in a PDB file (excluding backbone atoms).
    Returns the number of side chain atoms, the total number of residues, and the set of unique residue types.
    """
    side_chain_count = 0
    residue_set = set()
    residue_types = set()
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    atom_name = line[12:16].strip()
                    res_id = line[22:26].strip()
                    res_name = line[17:20].strip()
                    if atom_name not in BACKBONE_ATOMS:
                        side_chain_count += 1
                    residue_set.add(res_id)
                    residue_types.add(res_name)
    except Exception as e:
        logger.error(f"Error reading {pdb_file.name}: {str(e)}")
    return side_chain_count, len(residue_set), residue_types

def analyze_side_chains(pdb_files):
    """
    Analyze side chain counts for all PDB files.
    """
    stats = []
    for pdb_file in tqdm(pdb_files, desc="Processing PDB files"):
        side_chain_count, n_residues, residue_types = count_side_chains_in_pdb(pdb_file)
        stats.append({
            'file_name': pdb_file.name,
            'side_chain_count': side_chain_count,
            'n_residues': n_residues,
            'n_unique_side_chains': len(residue_types)
        })
    stats_df = pd.DataFrame(stats)
    return stats_df

def plot_side_chain_histogram(stats_df: pd.DataFrame, output_dir: str):
    """
    Plot a histogram of side chain counts per PDB file.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(8, 5))
    sns.histplot(stats_df['side_chain_count'], bins=50, color='skyblue')
    plt.xlabel('Side Chain Atom Count')
    plt.ylabel('Number of PDB Files')
    plt.title('Histogram of Side Chain Atom Counts')
    plt.tight_layout()
    plt.savefig(output_path / 'side_chain_histogram.png', dpi=300)
    plt.close()

def plot_side_chain_vs_residues(stats_df: pd.DataFrame, output_dir: str):
    """
    Plot a scatter plot of side chain count vs. number of residues.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x='n_residues', y='side_chain_count', data=stats_df, alpha=0.7)
    plt.xlabel('Number of Residues')
    plt.ylabel('Side Chain Atom Count')
    plt.title('Side Chain Atom Count vs. Number of Residues')
    plt.tight_layout()
    plt.savefig(output_path / 'side_chain_vs_residues.png', dpi=300)
    plt.close()

def plot_unique_side_chain_histogram(stats_df: pd.DataFrame, output_dir: str):
    """
    Plot a histogram of the number of unique RNA side chains per PDB file.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(8, 5))
    sns.histplot(stats_df['n_unique_side_chains'], bins=20, color='purple')
    plt.xlabel('Number of Unique RNA Side Chains (Residue Types)')
    plt.ylabel('Number of PDB Files')
    plt.title('Distribution of Unique RNA Side Chains per PDB')
    plt.tight_layout()
    plt.savefig(output_path / 'unique_side_chain_histogram.png', dpi=300)
    plt.close()

def save_results(stats_df: pd.DataFrame, output_dir: str):
    """
    Save side chain statistics to CSV.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    stats_df.to_csv(output_path / 'side_chain_statistics.csv', index=False)
    logger.info(f"Results saved to {output_path}")

def print_side_chain_summary(stats_df: pd.DataFrame):
    """
    Print summary statistics for side chain analysis.
    """
    print("\n=== PDB Side Chain Analysis Summary ===")
    print(f"Files analyzed: {len(stats_df)}")
    print(f"Side chain atoms: min={stats_df['side_chain_count'].min()}, max={stats_df['side_chain_count'].max()}, mean={stats_df['side_chain_count'].mean():.2f}")
    print(f"Residues: min={stats_df['n_residues'].min()}, max={stats_df['n_residues'].max()}, mean={stats_df['n_residues'].mean():.2f}")
    print(f"Unique RNA side chains: min={stats_df['n_unique_side_chains'].min()}, max={stats_df['n_unique_side_chains'].max()}, mean={stats_df['n_unique_side_chains'].mean():.2f}")

def main():
    """
    Main function to run PDB side chain analysis.
    """
    pdb_dir = "competition/official_training_pdbs"
    output_dir = "EDA/pdb_sidechain_analysis"
    try:
        pdb_files = list_pdb_files(pdb_dir)
        stats_df = analyze_side_chains(pdb_files)
        print_side_chain_summary(stats_df)
        plot_side_chain_histogram(stats_df, output_dir)
        plot_side_chain_vs_residues(stats_df, output_dir)
        plot_unique_side_chain_histogram(stats_df, output_dir)
        save_results(stats_df, output_dir)
        print(f"\nAnalysis complete! Results saved to {output_dir}")
    except Exception as e:
        logger.error(f"Error in side chain analysis: {str(e)}")
        raise

if __name__ == "__main__":
    main() 