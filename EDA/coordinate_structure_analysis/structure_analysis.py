"""
Coordinate Structure Analysis Script

Analyzes the shape and structure of coordinate arrays in .npy files.
Outputs statistics on array dimensions and visualizations.
"""

import os
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

def list_npy_files(coords_dir: str):
    """
    List all .npy files in the given directory.
    """
    coords_path = Path(coords_dir)
    if not coords_path.exists():
        raise FileNotFoundError(f"Coords directory not found: {coords_dir}")
    npy_files = list(coords_path.glob("*.npy"))
    logger.info(f"Found {len(npy_files)} NPY files in {coords_dir}")
    return npy_files

def analyze_structure(npy_files):
    """
    Analyze the shape and structure of each .npy file.
    """
    structure_stats = []
    for npy_file in tqdm(npy_files, desc="Processing NPY files"):
        try:
            arr = np.load(npy_file)
            shape = arr.shape
            # Typical shapes: (residues, atoms, coords) or (residues, atoms)
            n_residues = shape[0] if len(shape) > 0 else 0
            n_atoms = shape[1] if len(shape) > 1 else 0
            n_coords = shape[2] if len(shape) > 2 else 0
            structure_stats.append({
                'file_name': npy_file.name,
                'n_residues': n_residues,
                'n_atoms': n_atoms,
                'n_coords': n_coords,
                'shape': shape
            })
        except Exception as e:
            logger.error(f"Error reading {npy_file.name}: {str(e)}")
            continue
    structure_stats_df = pd.DataFrame(structure_stats)
    return structure_stats_df

def plot_structure_histograms(structure_stats_df: pd.DataFrame, output_dir: str):
    """
    Plot histograms of residues, atoms, and coordinates per file.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(12, 4))
    plt.subplot(1, 3, 1)
    sns.histplot(structure_stats_df['n_residues'], bins=50, color='skyblue')
    plt.xlabel('Number of Residues')
    plt.title('Residues per File')
    plt.subplot(1, 3, 2)
    sns.histplot(structure_stats_df['n_atoms'], bins=10, color='orange')
    plt.xlabel('Number of Atoms')
    plt.title('Atoms per File')
    plt.subplot(1, 3, 3)
    sns.histplot(structure_stats_df['n_coords'], bins=5, color='green')
    plt.xlabel('Number of Coordinates')
    plt.title('Coordinates per File')
    plt.tight_layout()
    plt.savefig(output_path / 'structure_histograms.png', dpi=300)
    plt.close()

def save_results(structure_stats_df: pd.DataFrame, output_dir: str):
    """
    Save structure statistics to CSV.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    structure_stats_df.to_csv(output_path / 'structure_statistics.csv', index=False)
    logger.info(f"Results saved to {output_path}")

def print_structure_summary(structure_stats_df: pd.DataFrame):
    """
    Print summary statistics for structure analysis.
    """
    print("\n=== Coordinate Structure Analysis Summary ===")
    print(f"Files analyzed: {len(structure_stats_df)}")
    print(f"Residues: min={structure_stats_df['n_residues'].min()}, max={structure_stats_df['n_residues'].max()}, mean={structure_stats_df['n_residues'].mean():.2f}")
    print(f"Atoms: min={structure_stats_df['n_atoms'].min()}, max={structure_stats_df['n_atoms'].max()}, mean={structure_stats_df['n_atoms'].mean():.2f}")
    print(f"Coords: min={structure_stats_df['n_coords'].min()}, max={structure_stats_df['n_coords'].max()}, mean={structure_stats_df['n_coords'].mean():.2f}")

def main():
    """
    Main function to run coordinate structure analysis.
    """
    coords_dir = "competition/train/coords"
    output_dir = "pdb_prune/EDA/coordinate_structure_analysis"
    try:
        npy_files = list_npy_files(coords_dir)
        structure_stats_df = analyze_structure(npy_files)
        print_structure_summary(structure_stats_df)
        plot_structure_histograms(structure_stats_df, output_dir)
        save_results(structure_stats_df, output_dir)
        print(f"\nAnalysis complete! Results saved to {output_dir}")
    except Exception as e:
        logger.error(f"Error in structure analysis: {str(e)}")
        raise

if __name__ == "__main__":
    main() 