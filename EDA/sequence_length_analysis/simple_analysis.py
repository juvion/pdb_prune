"""
RNA Sequence Length Analysis Script

Simple script to analyze RNA sequence lengths from FASTA files.
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

# Set style for better plots
plt.style.use('seaborn')
sns.set_palette("husl")


def read_fasta_files(fasta_dir: str) -> pd.DataFrame:
    """Read all FASTA files and extract sequence information.
    
    Args:
        fasta_dir (str): Path to directory containing FASTA files
        
    Returns:
        pd.DataFrame: DataFrame containing sequence metadata and lengths
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
            
            # Extract sequence ID and sequence
            sequence_id = None
            sequence = ""
            
            for line in lines:
                line = line.strip()
                if line.startswith('>'):
                    sequence_id = line[1:]  # Remove '>' prefix
                elif line and not line.startswith('>'):
                    sequence += line
            
            if sequence_id and sequence:
                sequences_data.append({
                    'file_name': fasta_file.name,
                    'sequence_id': sequence_id,
                    'sequence': sequence,
                    'length': len(sequence),
                    'gc_content': calculate_gc_content(sequence)
                })
            else:
                logger.warning(f"Invalid FASTA format in {fasta_file.name}")
                
        except Exception as e:
            logger.error(f"Error reading {fasta_file.name}: {str(e)}")
            continue
    
    sequences_df = pd.DataFrame(sequences_data)
    logger.info(f"Successfully processed {len(sequences_df)} sequences")
    
    return sequences_df


def calculate_gc_content(sequence: str) -> float:
    """Calculate GC content of a sequence.
    
    Args:
        sequence (str): RNA sequence
        
    Returns:
        float: GC content as percentage
    """
    if not sequence:
        return 0.0
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return (gc_count / len(sequence)) * 100


def calculate_statistics(sequences_df: pd.DataFrame) -> dict:
    """Calculate comprehensive statistics on sequence lengths.
    
    Args:
        sequences_df (pd.DataFrame): DataFrame containing sequence data
        
    Returns:
        dict: Dictionary containing statistical summaries
    """
    if sequences_df.empty:
        raise ValueError("No sequence data available.")
    
    lengths = sequences_df['length'].values
    
    stats = {
        'total_sequences': len(sequences_df),
        'min_length': np.min(lengths),
        'max_length': np.max(lengths),
        'mean_length': np.mean(lengths),
        'median_length': np.median(lengths),
        'std_length': np.std(lengths),
        'q25_percentile': np.percentile(lengths, 25),
        'q75_percentile': np.percentile(lengths, 75),
        'iqr': np.percentile(lengths, 75) - np.percentile(lengths, 25),
        'skewness': calculate_skewness(lengths),
        'kurtosis': calculate_kurtosis(lengths)
    }
    
    return stats


def calculate_skewness(data: np.ndarray) -> float:
    """Calculate skewness of the data.
    """
    mean = np.mean(data)
    std = np.std(data)
    if std == 0:
        return 0
    return np.mean(((data - mean) / std) ** 3)


def calculate_kurtosis(data: np.ndarray) -> float:
    """Calculate kurtosis of the data.
    """
    mean = np.mean(data)
    std = np.std(data)
    if std == 0:
        return 0
    return np.mean(((data - mean) / std) ** 4) - 3


def create_visualizations(sequences_df: pd.DataFrame, output_dir: str):
    """Create comprehensive visualizations of sequence length distribution.
    
    Args:
        sequences_df (pd.DataFrame): DataFrame containing sequence data
        output_dir (str): Directory to save visualization files
    """
    if sequences_df.empty:
        raise ValueError("No sequence data available.")
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    lengths = sequences_df['length'].values
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('RNA Sequence Length Distribution Analysis', fontsize=16, fontweight='bold')
    
    # 1. Histogram with density curve
    axes[0, 0].hist(lengths, bins=50, alpha=0.7, density=True, color='skyblue', edgecolor='black')
    axes[0, 0].axvline(np.mean(lengths), color='red', linestyle='--', label=f'Mean: {np.mean(lengths):.1f}')
    axes[0, 0].axvline(np.median(lengths), color='green', linestyle='--', label=f'Median: {np.median(lengths):.1f}')
    axes[0, 0].set_xlabel('Sequence Length')
    axes[0, 0].set_ylabel('Density')
    axes[0, 0].set_title('Length Distribution Histogram')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # 2. Box plot
    axes[0, 1].boxplot(lengths, patch_artist=True, boxprops=dict(facecolor='lightblue'))
    axes[0, 1].set_xlabel('Sequence Length')
    axes[0, 1].set_title('Length Distribution Box Plot')
    axes[0, 1].grid(True, alpha=0.3)
    
    # 3. Cumulative distribution
    sorted_lengths = np.sort(lengths)
    cumulative_prob = np.arange(1, len(sorted_lengths) + 1) / len(sorted_lengths)
    axes[1, 0].plot(sorted_lengths, cumulative_prob, linewidth=2, color='purple')
    axes[1, 0].set_xlabel('Sequence Length')
    axes[1, 0].set_ylabel('Cumulative Probability')
    axes[1, 0].set_title('Cumulative Distribution Function')
    axes[1, 0].grid(True, alpha=0.3)
    
    # 4. Q-Q plot (normal distribution)
    from scipy import stats
    stats.probplot(lengths, dist='norm', plot=axes[1, 1])
    axes[1, 1].set_title('Q-Q Plot (Normal Distribution)')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path / 'sequence_length_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create additional detailed plots
    create_length_vs_gc_plot(sequences_df, output_path)
    create_length_distribution_by_category(sequences_df, output_path)
    
    logger.info(f"Visualizations saved to {output_path}")


def create_length_vs_gc_plot(sequences_df: pd.DataFrame, output_path: Path):
    """Create scatter plot of length vs GC content.
    """
    plt.figure(figsize=(10, 6))
    plt.scatter(sequences_df['length'], sequences_df['gc_content'], 
               alpha=0.6, s=20)
    plt.xlabel('Sequence Length')
    plt.ylabel('GC Content (%)')
    plt.title('Sequence Length vs GC Content')
    plt.grid(True, alpha=0.3)
    plt.savefig(output_path / 'length_vs_gc_content.png', dpi=300, bbox_inches='tight')
    plt.close()


def create_length_distribution_by_category(sequences_df: pd.DataFrame, output_path: Path):
    """
    Create length distribution categorized by length ranges.
    """
    # Define length categories
    sequences_df['length_category'] = pd.cut(
        sequences_df['length'], 
        bins=[0, 50, 100, 200, np.inf],
        labels=['0-50', '51-100', '101-200', '>200']
    )
    plt.figure(figsize=(10, 6))
    category_counts = sequences_df['length_category'].value_counts().sort_index()
    palette = plt.get_cmap('tab10')
    color = palette(0)  # Use the first color from tab10
    bars = plt.bar(
        range(len(category_counts)),
        category_counts.values,
        color=color,
        edgecolor='black',
        linewidth=1.5
    )
    plt.xlabel('Length Category', fontsize=13)
    plt.ylabel('Number of Sequences', fontsize=13)
    plt.title('Sequence Distribution by Length Category', fontsize=15, fontweight='bold')
    plt.xticks(range(len(category_counts)), category_counts.index, rotation=0, fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.3)
    # Add value labels on bars
    for bar, count in zip(bars, category_counts.values):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(category_counts.values)*0.01,
                 str(count), ha='center', va='bottom', fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_path / 'length_category_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()


def save_results(sequences_df: pd.DataFrame, stats: dict, output_dir: str):
    """Save analysis results to CSV files.
    
    Args:
        sequences_df (pd.DataFrame): DataFrame containing sequence data
        stats (dict): Dictionary containing statistics
        output_dir (str): Directory to save results
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Save sequence data
    sequences_df.to_csv(output_path / 'sequence_data.csv', index=False)
    
    # Save statistics
    stats_df = pd.DataFrame([stats])
    stats_df.to_csv(output_path / 'length_statistics.csv', index=False)
    
    # Save summary by length category
    if 'length_category' in sequences_df.columns:
        category_summary = sequences_df.groupby('length_category').agg({
            'length': ['count', 'mean', 'std'],
            'gc_content': ['mean', 'std']
        }).round(2)
        category_summary.to_csv(output_path / 'length_category_summary.csv')
    
    logger.info(f"Results saved to {output_path}")


def main():
    """Main function to run sequence length analysis.
    """
    # Configuration
    fasta_dir = "competition/train/seqs"
    output_dir = "pdb_prune/EDA/sequence_length_analysis"
    
    try:
        # Read FASTA files
        sequences_df = read_fasta_files(fasta_dir)
        
        # Calculate statistics
        stats = calculate_statistics(sequences_df)
        
        # Print summary statistics
        print("\n=== Sequence Length Analysis Summary ===")
        print(f"Total sequences: {stats['total_sequences']}")
        print(f"Length range: {stats['min_length']} - {stats['max_length']}")
        print(f"Mean length: {stats['mean_length']:.2f}")
        print(f"Median length: {stats['median_length']:.2f}")
        print(f"Standard deviation: {stats['std_length']:.2f}")
        print(f"Skewness: {stats['skewness']:.3f}")
        print(f"Kurtosis: {stats['kurtosis']:.3f}")
        
        # Create visualizations
        create_visualizations(sequences_df, output_dir)
        
        # Save results
        save_results(sequences_df, stats, output_dir)
        
        print(f"\nAnalysis complete! Results saved to {output_dir}")
        
    except Exception as e:
        logger.error(f"Error in sequence analysis: {str(e)}")
        raise


if __name__ == "__main__":
    main()
