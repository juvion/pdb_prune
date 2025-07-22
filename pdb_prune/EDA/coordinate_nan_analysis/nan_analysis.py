"""
NaN Analysis in Coordinate Files

Analyzes NaN patterns in .npy coordinate files, excluding the first phosphorus atom (P).
Outputs per-file NaN statistics, histogram, and heatmap visualizations.
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

def analyze_nan_patterns(npy_files):
    """
    Analyze NaN patterns in each .npy file, excluding the first P atom.
    Also count the number of residues with at least one NaN per file.
    """
    nan_stats = []
    nan_positions = []
    residues_with_nan_counts = []
    for npy_file in tqdm(npy_files, desc="Processing NPY files"):
        try:
            arr = np.load(npy_file)
            arr_excl_first = arr[1:] if arr.shape[0] > 1 else arr
            nan_mask = np.isnan(arr_excl_first)
            nan_count = np.sum(nan_mask)
            # Count residues (rows) with at least one NaN
            if nan_mask.ndim >= 2:
                residues_with_nan = np.any(nan_mask, axis=tuple(range(1, nan_mask.ndim)))
            else:
                residues_with_nan = nan_mask
            n_residues_with_nan = np.sum(residues_with_nan)
            nan_stats.append({
                'file_name': npy_file.name,
                'nan_count': nan_count,
                'shape': arr.shape,
                'n_residues_with_nan': int(n_residues_with_nan)
            })
            nan_positions.append(np.sum(nan_mask, axis=0) if nan_mask.ndim > 1 else nan_mask.astype(int))
            residues_with_nan_counts.append(int(n_residues_with_nan))
        except Exception as e:
            logger.error(f"Error reading {npy_file.name}: {str(e)}")
            continue
    nan_stats_df = pd.DataFrame(nan_stats)
    nan_positions_arr = np.stack([x if x.shape == nan_positions[0].shape else np.zeros_like(nan_positions[0]) for x in nan_positions]) if nan_positions else np.array([])
    return nan_stats_df, nan_positions_arr, residues_with_nan_counts

def plot_nan_histogram(nan_stats_df: pd.DataFrame, output_dir: str):
    """
    Plot a histogram of NaN counts per file.
    """
    plt.figure(figsize=(8, 5))
    sns.histplot(nan_stats_df['nan_count'], bins=50, kde=False, color='skyblue')
    plt.xlabel('NaN Count per File')
    plt.ylabel('Number of Files')
    plt.title('Histogram of NaN Counts per Coordinate File')
    plt.tight_layout()
    plt.savefig(Path(output_dir) / 'nan_count_histogram.png', dpi=300)
    plt.close()

def plot_nan_heatmap(nan_positions_arr: np.ndarray, output_dir: str):
    """
    Plot a heatmap showing NaN positions across all files.
    """
    if nan_positions_arr.size == 0:
        logger.warning("No NaN position data to plot heatmap.")
        return
    # If 3D, sum over the last axis to get (files, atoms)
    if nan_positions_arr.ndim == 3:
        heatmap_data = nan_positions_arr.sum(axis=2)
    else:
        heatmap_data = nan_positions_arr
    plt.figure(figsize=(12, 6))
    sns.heatmap(heatmap_data, cmap='viridis', cbar_kws={'label': 'NaN Count'})
    plt.title('NaN Positions Heatmap (across all files, excluding first atom)')
    plt.xlabel('Atom Index')
    plt.ylabel('File Index')
    plt.tight_layout()
    plt.savefig(Path(output_dir) / 'nan_positions_heatmap.png', dpi=300)
    plt.close()

def save_results(nan_stats_df: pd.DataFrame, output_dir: str):
    """
    Save NaN statistics to CSV.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    nan_stats_df.to_csv(output_path / 'nan_statistics.csv', index=False)
    logger.info(f"Results saved to {output_path}")

def save_residues_with_nan_distribution(residues_with_nan_counts, output_dir: str):
    """
    Save a summary table: distribution of number of PDBs by number of residues with NaN.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    value_counts = pd.Series(residues_with_nan_counts).value_counts().sort_index()
    summary_df = pd.DataFrame({
        'n_residues_with_nan': value_counts.index,
        'n_pdbs': value_counts.values
    })
    summary_df.to_csv(output_path / 'residues_with_nan_distribution.csv', index=False)
    logger.info(f"Residues-with-NaN distribution saved to {output_path}")

def print_nan_summary(nan_stats_df: pd.DataFrame):
    """
    Print summary statistics for NaN analysis.
    """
    total_files = len(nan_stats_df)
    files_with_nan = (nan_stats_df['nan_count'] > 0).sum()
    percent_with_nan = (files_with_nan / total_files) * 100 if total_files > 0 else 0
    avg_nan_count = nan_stats_df['nan_count'].mean() if total_files > 0 else 0
    print(f"\n=== NaN Analysis Summary ===")
    print(f"Total files: {total_files}")
    print(f"Files with NaN: {files_with_nan} ({percent_with_nan:.2f}%)")
    print(f"Average NaN count per file: {avg_nan_count:.2f}")

def main():
    """
    Main function to run NaN analysis on coordinate files.
    """
    coords_dir = "competition/train/coords"
    output_dir = "EDA/coordinate_nan_analysis"
    try:
        npy_files = list_npy_files(coords_dir)
        nan_stats_df, nan_positions_arr, residues_with_nan_counts = analyze_nan_patterns(npy_files)
        print_nan_summary(nan_stats_df)
        plot_nan_histogram(nan_stats_df, output_dir)
        plot_nan_heatmap(nan_positions_arr, output_dir)
        save_results(nan_stats_df, output_dir)
        save_residues_with_nan_distribution(residues_with_nan_counts, output_dir)
        print(f"\nAnalysis complete! Results saved to {output_dir}")
    except Exception as e:
        logger.error(f"Error in NaN analysis: {str(e)}")
        raise

if __name__ == "__main__":
    main() 