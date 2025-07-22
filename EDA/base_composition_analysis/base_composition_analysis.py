"""
Base Composition Analysis Script

Analyzes the frequency of each RNA base (A, U, G, C) across all sequences in the FASTA files.
Outputs per-sequence and overall base composition, bar chart, and heatmap visualizations.
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

# Unified plot style settings
import seaborn as sns
sns.set_theme(style="whitegrid", palette="tab10")
import matplotlib as mpl
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.titlesize": 23,  # Title font size
    "axes.labelsize": 19,  # Label font size
    "xtick.labelsize": 17,
    "ytick.labelsize": 17,
    "legend.fontsize": 17,
    "figure.figsize": (8, 6),
    "axes.grid": True,
    "grid.alpha": 0.3,
    "axes.facecolor": "white",
    "savefig.dpi": 300,
    "lines.linewidth": 2,
})

sns.set_theme(style="whitegrid", palette="tab10")

BASES = ['A', 'U', 'G', 'C']

def read_fasta_files(fasta_dir: str) -> pd.DataFrame:
    """
    Read all FASTA files and extract sequence information.
    """
    fasta_path = Path(fasta_dir)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA directory not found: {fasta_dir}")
    logger.info(f"Reading FASTA files from {fasta_path}")
    fasta_files = list(fasta_path.glob("*.fasta"))
    logger.info(f"Found {len(fasta_files)} FASTA files")
    sequences_data = []
    for fasta_file in tqdm(fasta_files, desc="Processing FASTA files"):
        try:
            with open(fasta_file, 'r') as f:
                lines = f.readlines()
            sequence_id = None
            sequence = ""
            for line in lines:
                line = line.strip()
                if line.startswith('>'):
                    sequence_id = line[1:]
                elif line and not line.startswith('>'):
                    sequence += line
            if sequence_id and sequence:
                sequences_data.append({
                    'file_name': fasta_file.name,
                    'sequence_id': sequence_id,
                    'sequence': sequence
                })
            else:
                logger.warning(f"Invalid FASTA format in {fasta_file.name}")
        except Exception as e:
            logger.error(f"Error reading {fasta_file.name}: {str(e)}")
            continue
    sequences_df = pd.DataFrame(sequences_data)
    logger.info(f"Successfully processed {len(sequences_df)} sequences")
    return sequences_df

def calculate_base_composition(sequence: str) -> dict:
    """
    Calculate the percentage of each base in a sequence.
    """
    seq = sequence.upper()
    length = len(seq)
    if length == 0:
        return {base: 0.0 for base in BASES}
    base_counts = {base: seq.count(base) for base in BASES}
    base_percents = {base: (base_counts[base] / length) * 100 for base in BASES}
    return base_percents

def analyze_base_composition(sequences_df: pd.DataFrame) -> pd.DataFrame:
    """
    Add base composition columns to the DataFrame.
    """
    for base in BASES:
        sequences_df[f'percent_{base}'] = sequences_df['sequence'].apply(lambda seq: calculate_base_composition(seq)[base])
    return sequences_df

def overall_base_composition(sequences_df: pd.DataFrame) -> dict:
    """
    Calculate overall base composition for the dataset.
    """
    all_seq = ''.join(sequences_df['sequence'].tolist()).upper()
    total_length = len(all_seq)
    if total_length == 0:
        return {base: 0.0 for base in BASES}
    base_counts = {base: all_seq.count(base) for base in BASES}
    base_percents = {base: (base_counts[base] / total_length) * 100 for base in BASES}
    return base_percents

def plot_base_composition_bar(overall_comp: dict, output_dir: str):
    """
    Plot a bar chart of overall base composition using Tableau 10 blue for all bars.
    """
    import seaborn as sns
    sns.set_theme(style="whitegrid", palette="tab10")
    import matplotlib as mpl
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "axes.titlesize": 23,  # Title font size
        "axes.labelsize": 19,  # Label font size
        "xtick.labelsize": 17,
        "ytick.labelsize": 17,
        "legend.fontsize": 17,
        "figure.figsize": (8, 6),
        "axes.grid": True,
        "grid.alpha": 0.3,
        "axes.facecolor": "white",
        "savefig.dpi": 300,
        "lines.linewidth": 2,
    })

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(6, 5))
    palette = plt.get_cmap('tab10')
    color = palette(0)
    bars = plt.bar(list(overall_comp.keys()), list(overall_comp.values()), color=color, alpha=0.7)
    plt.ylabel('Percentage (%)')
    plt.xlabel('Base')
    plt.title('Overall Base Composition')
    plt.ylim(0, 100)
    for i, v in enumerate(overall_comp.values()):
        plt.text(i, v + 1, f"{v:.1f}", ha='center')
    plt.tight_layout()
    plt.savefig(output_path / 'overall_base_composition.png', dpi=300)
    plt.close()

def plot_base_composition_heatmap(sequences_df: pd.DataFrame, output_dir: str):
    """
    Plot a heatmap of per-sequence base composition.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    heatmap_data = sequences_df[[f'percent_{base}' for base in BASES]]
    plt.figure(figsize=(10, 8))
    sns.heatmap(heatmap_data, cmap='viridis', cbar_kws={'label': 'Percentage'})
    plt.title('Per-Sequence Base Composition Heatmap')
    plt.xlabel('Base')
    plt.ylabel('Sequence Index')
    plt.tight_layout()
    plt.savefig(output_path / 'base_composition_heatmap.png', dpi=300)
    plt.close()

def save_results(sequences_df: pd.DataFrame, overall_comp: dict, output_dir: str):
    """
    Save per-sequence and overall base composition to CSV.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    # Per-sequence
    sequences_df.to_csv(output_path / 'per_sequence_base_composition.csv', index=False)
    # Overall
    pd.DataFrame([overall_comp]).to_csv(output_path / 'overall_base_composition.csv', index=False)
    logger.info(f"Results saved to {output_path}")

def main():
    """
    Main function to run base composition analysis.
    """
    fasta_dir = "competition/train/seqs"
    output_dir = "EDA/base_composition_analysis"
    try:
        sequences_df = read_fasta_files(fasta_dir)
        sequences_df = analyze_base_composition(sequences_df)
        overall_comp = overall_base_composition(sequences_df)
        print("\n=== Overall Base Composition (%) ===")
        for base, percent in overall_comp.items():
            print(f"{base}: {percent:.2f}%")
        plot_base_composition_bar(overall_comp, output_dir)
        plot_base_composition_heatmap(sequences_df, output_dir)
        save_results(sequences_df, overall_comp, output_dir)
        print(f"\nAnalysis complete! Results saved to {output_dir}")
    except Exception as e:
        logger.error(f"Error in base composition analysis: {str(e)}")
        raise

if __name__ == "__main__":
    main() 