import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Tuple

# Set seaborn style
sns.set(style="whitegrid", context="notebook")

# Output directory
OUTPUT_DIR = Path(__file__).parent
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Input files (relative to EDA root)
SEQ_CSV = OUTPUT_DIR.parent / "sequence_length_analysis/sequence_data.csv"
NAN_CSV = OUTPUT_DIR.parent / "coordinate_nan_analysis/nan_statistics.csv"
STRUCT_CSV = OUTPUT_DIR.parent / "coordinate_structure_analysis/structure_statistics.csv"


def extract_base_name(file_name: str) -> str:
    """
    Extract the base name (without extension) for joining.

    Args:
        file_name (str): File name (e.g., '7AOE_1_R.fasta' or '7AOE_1_R.npy')

    Returns:
        str: Base name (e.g., '7AOE_1_R')
    """
    return Path(file_name).stem


def load_and_merge_data() -> pd.DataFrame:
    """
    Load sequence, NaN, and structure statistics and merge on base name.

    Returns:
        pd.DataFrame: Merged dataframe with all relevant columns.
    """
    # Load sequence data
    seq_df = pd.read_csv(SEQ_CSV)
    seq_df["base_name"] = seq_df["file_name"].apply(extract_base_name)

    # Load NaN statistics
    nan_df = pd.read_csv(NAN_CSV)
    nan_df["base_name"] = nan_df["file_name"].apply(extract_base_name)

    # Load structure statistics
    struct_df = pd.read_csv(STRUCT_CSV)
    struct_df["base_name"] = struct_df["file_name"].apply(extract_base_name)

    # Merge all
    merged = seq_df.merge(nan_df, on="base_name", how="inner", suffixes=("_seq", "_nan"))
    merged = merged.merge(struct_df, on="base_name", how="inner", suffixes=("", "_struct"))

    # Clean up columns
    merged = merged.drop(columns=[
        "file_name_seq", "file_name_nan", "file_name", "shape_nan", "shape"
    ], errors="ignore")
    return merged


def plot_scatter_with_reg(x, y, data, xlabel, ylabel, title, out_path):
    """
    Scatter plot with regression line.
    """
    plt.figure(figsize=(8, 6))
    sns.regplot(x=x, y=y, data=data, scatter_kws={"s": 15, "alpha": 0.5}, line_kws={"color": "red"})
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def plot_correlation_heatmap(df, cols, out_path):
    """
    Plot a correlation heatmap for selected columns.
    """
    corr = df[cols].corr(method="pearson")
    plt.figure(figsize=(7, 6))
    sns.heatmap(corr, annot=True, fmt=".2f", cmap="coolwarm", square=True)
    plt.title("Pearson Correlation Heatmap")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def plot_box_by_category(df, category_col, value_col, ylabel, title, out_path):
    """
    Boxplot of value_col by category_col.
    """
    plt.figure(figsize=(10, 6))
    sns.boxplot(x=category_col, y=value_col, data=df, palette="tab10")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def main():
    """
    Main function to run cross-dataset correlation analysis.
    """
    merged = load_and_merge_data()
    merged.to_csv(OUTPUT_DIR / "merged_cross_correlation.csv", index=False)

    # Correlation: sequence length vs. n_residues (should be highly correlated)
    plot_scatter_with_reg(
        x="length", y="n_residues", data=merged,
        xlabel="Sequence Length (FASTA)", ylabel="n_residues (Coords)",
        title="Sequence Length vs. Coordinate Array Residue Count",
        out_path=OUTPUT_DIR / "length_vs_n_residues.png"
    )

    # Correlation: sequence length vs. nan_count
    plot_scatter_with_reg(
        x="length", y="nan_count", data=merged,
        xlabel="Sequence Length", ylabel="NaN Count (Coords)",
        title="Sequence Length vs. NaN Count in Coordinates",
        out_path=OUTPUT_DIR / "length_vs_nan_count.png"
    )

    # Correlation: GC content vs. nan_count
    plot_scatter_with_reg(
        x="gc_content", y="nan_count", data=merged,
        xlabel="GC Content (%)", ylabel="NaN Count (Coords)",
        title="GC Content vs. NaN Count in Coordinates",
        out_path=OUTPUT_DIR / "gc_content_vs_nan_count.png"
    )

    # Correlation: n_residues_with_nan vs. length
    plot_scatter_with_reg(
        x="length", y="n_residues_with_nan", data=merged,
        xlabel="Sequence Length", ylabel="# Residues with NaN",
        title="Sequence Length vs. # Residues with NaN",
        out_path=OUTPUT_DIR / "length_vs_n_residues_with_nan.png"
    )

    # Correlation heatmap
    corr_cols = ["length", "gc_content", "n_residues", "nan_count", "n_residues_with_nan"]
    plot_correlation_heatmap(
        merged, corr_cols, OUTPUT_DIR / "correlation_heatmap.png"
    )

    # Boxplot: length_category vs. nan_count
    plot_box_by_category(
        merged, "length_category", "nan_count",
        ylabel="NaN Count (Coords)",
        title="NaN Count by Sequence Length Category",
        out_path=OUTPUT_DIR / "nan_count_by_length_category.png"
    )

    # Boxplot: length_category vs. n_residues_with_nan
    plot_box_by_category(
        merged, "length_category", "n_residues_with_nan",
        ylabel="# Residues with NaN",
        title="# Residues with NaN by Sequence Length Category",
        out_path=OUTPUT_DIR / "n_residues_with_nan_by_length_category.png"
    )

    # Save summary statistics
    summary = merged[["length", "gc_content", "n_residues", "nan_count", "n_residues_with_nan"]].describe()
    summary.to_csv(OUTPUT_DIR / "cross_correlation_summary_statistics.csv")

    print("Cross-dataset correlation analysis complete. Outputs saved to:", OUTPUT_DIR)

if __name__ == "__main__":
    main() 