# RNA Clustering and Data Splitting Pipeline

A simple, effective pipeline for clustering RNA data based on sequence and structure similarity, followed by cluster-based splitting into train/validation/test sets.

## Overview

This pipeline implements the clustering methodology described in RNA structure prediction papers, providing:

- **Sequence clustering** using PSI-CD-HIT at thresholds 0.8, 0.9, 1.0
- **Structure clustering** using US-align TM-scores at thresholds 0.4, 0.45, 0.5, 0.6
- **Cluster-based data splitting** ensuring no data leakage between train/val/test sets
- **Length filtering** to manage computational resources

## Quick Start

### 1. Install Dependencies

```bash
# Install external tools
chmod +x utils/install_clustering_tools.sh
./utils/install_clustering_tools.sh

# Restart terminal or source bashrc
source ~/.bashrc

# Install Python dependencies
pip install networkx
```

### 2. Run Pipeline

```bash
# Using the example script (recommended)
python run_clustering_example.py

# Or run directly with custom parameters
python utils/rna_clustering_pipeline.py \
    --pdb_dir data/experiments_data/exp2.2_manuscript/dataset_20250731/extracted_pdbs \
    --fasta_dir data/experiments_data/exp2.2_manuscript/dataset_20250731/extracted_sequences \
    --output_dir data/experiments_data/exp2.2_manuscript/dataset_20250731/clustering_results \
    --max_seq_len 500
```

## Input Requirements

- **PDB files**: 3D RNA structures (`.pdb` format)
- **FASTA files**: Corresponding nucleotide sequences (`.fasta` format)
- **File naming**: One-to-one correspondence based on filename (e.g., `1ABC.pdb` ↔ `1ABC.fasta`)

## Output Files

The pipeline generates split files for each clustering method and threshold:

```
clustering_results/
├── seq_0.8_train.txt     # Sequence clustering at 80% similarity
├── seq_0.8_val.txt
├── seq_0.8_test.txt
├── seq_0.9_train.txt     # Sequence clustering at 90% similarity
├── seq_0.9_val.txt
├── seq_0.9_test.txt
├── seq_1.0_train.txt     # Sequence clustering at 100% similarity
├── seq_1.0_val.txt
├── seq_1.0_test.txt
├── struct_0.4_train.txt  # Structure clustering at TM-score 0.4
├── struct_0.4_val.txt
├── struct_0.4_test.txt
├── struct_0.45_train.txt # Structure clustering at TM-score 0.45
├── struct_0.45_val.txt
├── struct_0.45_test.txt
├── struct_0.5_train.txt  # Structure clustering at TM-score 0.5
├── struct_0.5_val.txt
├── struct_0.5_test.txt
├── struct_0.6_train.txt  # Structure clustering at TM-score 0.6
├── struct_0.6_val.txt
└── struct_0.6_test.txt
```

Each file contains PDB IDs (one per line) assigned to that split.

## Key Features

### Cluster-Based Splitting
- **No data leakage**: Entire clusters are assigned to single splits
- **Target ratios**: ~75% train, ~15% test, ~10% validation
- **Large cluster handling**: Clusters >30 samples go to training set
- **Length balancing**: Maintains similar average sequence lengths across splits

### Computational Efficiency
- **Length filtering**: Configurable maximum sequence length (default: 500nt)
- **Parallel processing**: Multi-threaded PSI-CD-HIT execution
- **Memory management**: Efficient TM-score matrix handling
- **Progress logging**: Detailed progress and statistics

### Error Handling
- **Robust file I/O**: Handles missing files gracefully
- **Tool validation**: Checks for required external tools
- **Comprehensive logging**: Detailed error messages and warnings

## External Tools

### PSI-CD-HIT
- **Purpose**: Fast sequence clustering at high identity
- **Installation**: Automated via `install_clustering_tools.sh`
- **Usage**: Clusters sequences based on similarity thresholds

### US-align
- **Purpose**: Universal structure alignment for TM-score calculation
- **Installation**: Automated via `install_clustering_tools.sh`
- **Usage**: Computes structural similarity between RNA 3D structures

## Troubleshooting

### Common Issues

1. **Tools not found**:
   ```bash
   # Verify installation
   which psicd-hit
   which US-align
   
   # Re-run installation if needed
   ./utils/install_clustering_tools.sh
   source ~/.bashrc
   ```

2. **Memory issues**:
   - Reduce `max_seq_len` parameter
   - Process smaller batches of structures
   - Increase system memory allocation

3. **Missing files**:
   - Check PDB/FASTA file correspondence
   - Verify file permissions
   - Review log messages for specific errors

### Performance Tips

- **Sequence length**: Lower `max_seq_len` for faster processing
- **Parallel processing**: Adjust PSI-CD-HIT thread count based on CPU cores
- **Disk space**: Ensure sufficient space for intermediate files
- **Memory**: Monitor memory usage for large datasets

## Customization

The pipeline can be easily customized by modifying:

- **Clustering thresholds**: Edit `seq_thresholds` and `struct_thresholds` in the pipeline class
- **Split ratios**: Adjust `test_ratio`, `val_ratio`, `train_ratio`
- **Large cluster threshold**: Modify the 30-sample threshold for training assignment
- **Tool parameters**: Customize PSI-CD-HIT and US-align command-line options

## Citation

If you use this pipeline in your research, please cite the relevant papers for the clustering methodologies and tools used.