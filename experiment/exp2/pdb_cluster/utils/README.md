# RNA PDB Data Clustering Pipeline

A comprehensive pipeline for clustering RNA structures from PDB database using both **sequence similarity** (CD-HIT) and **structure similarity** (TM-score with agglomerative clustering).

This implementation is based on Method 2 from the paper "RiboDiffusion: tertiary structure-based RNA inverse folding with generative diffusion models" (Huang et al., 2024).

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Detailed Usage](#detailed-usage)
- [Output Files](#output-files)
- [Troubleshooting](#troubleshooting)
- [References](#references)

## Overview

This pipeline performs dual clustering on RNA PDB structures:

1. **Sequence Similarity Clustering**: Uses PSI-CD-HIT to cluster based on nucleotide sequence identity
2. **Structure Similarity Clustering**: Uses US-align to compute TM-scores, followed by agglomerative clustering

### Features

- ✅ Automatic extraction of RNA chains from PDB files
- ✅ Length filtering (20-280 nucleotides)
- ✅ Multiple clustering thresholds (configurable)
- ✅ Comprehensive visualization (dendrograms, distributions)
- ✅ Resume capability for TM-score calculation
- ✅ Detailed summary reports

## Requirements

### System Requirements

- Linux/Unix or macOS
- 4GB+ RAM
- Python 3.7 or higher
- C++ compiler (gcc/g++)

### Software Dependencies

1. **CD-HIT** (for sequence clustering)
2. **US-align** (for structure alignment)
3. **Python packages**:
   - biopython
   - numpy
   - scipy
   - pandas
   - matplotlib
   - tqdm

## Installation

### Step 1: Install Python Packages

```bash
pip install biopython numpy scipy pandas matplotlib tqdm
```

Or using conda:

```bash
conda install -c conda-forge biopython numpy scipy pandas matplotlib tqdm
```

### Step 2: Install CD-HIT

**Option A: Using conda (recommended)**
```bash
conda install -c bioconda cd-hit
```

**Option B: From source**
```bash
git clone https://github.com/weizhongli/cdhit.git
cd cdhit
make
sudo make install  # or copy to ~/bin
```

Verify installation:
```bash
cd-hit-est -h
```

### Step 3: Install US-align

```bash
# Download and compile
git clone https://github.com/pylelab/USalign.git
cd USalign
g++ -static -O3 -ffast-math -lm -o USalign USalign.cpp

# Install
sudo cp USalign /usr/local/bin/
# Or for local installation:
mkdir -p ~/bin
cp USalign ~/bin/
export PATH=$PATH:~/bin
```

Verify installation:
```bash
USalign -h
```

### Step 4: Download Pipeline Scripts

Save all the provided Python scripts to a directory:

```bash
mkdir rna_clustering_pipeline
cd rna_clustering_pipeline

# Download or copy these files:
# - extract_rna_chains.py
# - extract_sequences.py
# - parse_cdhit_clusters.py
# - calculate_tm_matrix.py
# - agglomerative_clustering.py
# - run_clustering.sh

# Make shell script executable
chmod +x run_clustering.sh
```

## Quick Start

### Using the Automated Pipeline

1. **Prepare your data:**
   ```bash
   mkdir my_project
   cd my_project
   mkdir pdb_files
   
   # Place your PDB files in pdb_files/
   cp /path/to/your/*.pdb pdb_files/
   ```

2. **Run the complete pipeline:**
   ```bash
   # Basic usage
   ./run_clustering.sh pdb_files results
   
   # With custom thresholds
   ./run_clustering.sh pdb_files results 0.8,0.6,0.4 0.6,0.5,0.4
   ```

3. **Check results:**
   ```bash
   ls results/
   cat results/clustering_summary_report.txt
   ```

### Using Individual Scripts

If you prefer step-by-step execution:

```bash
# Step 1: Extract RNA chains
python3 extract_rna_chains.py --pdb_dir ./pdb_files --output_dir ./rna_chains

# Step 2: Extract sequences
python3 extract_sequences.py --pdb_dir ./rna_chains --output rna_sequences.fasta

# Step 3: Sequence clustering
cd-hit-est -i rna_sequences.fasta -o clusters_seq_0.8.fasta -c 0.8 -n 4 -M 2000 -T 4
python3 parse_cdhit_clusters.py --input clusters_seq_0.8.fasta.clstr --output sequence_clusters_0.8.csv

# Step 4: Calculate TM-score matrix
python3 calculate_tm_matrix.py --pdb_dir ./rna_chains --output tm_score_matrix.csv

# Step 5: Structure clustering
python3 agglomerative_clustering.py --input tm_score_matrix.csv --thresholds 0.6 0.5 0.4 --plot
```

## Detailed Usage

### Script 1: Extract RNA Chains

```bash
python3 extract_rna_chains.py --pdb_dir <input_dir> --output_dir <output_dir>
```

**Options:**
- `--pdb_dir`: Directory containing PDB files
- `--output_dir`: Output directory for extracted RNA chains

**Filters:**
- Sequence length: 20-280 nucleotides
- Only RNA residues (A, U, G, C)

### Script 2: Extract Sequences

```bash
python3 extract_sequences.py --pdb_dir <rna_chains_dir> --output <fasta_file>
```

### Script 3: Sequence Clustering (CD-HIT)

```bash
cd-hit-est -i <input.fasta> -o <output> -c <threshold> -n <word_size> -M <memory> -T <threads>
```

**Parameters:**
- `-c`: Sequence identity threshold (0.0-1.0)
- `-n`: Word size (4 for ≥0.6, 3 for ≥0.5, 2 for ≥0.4)
- `-M`: Memory limit in MB
- `-T`: Number of threads

**Example:**
```bash
# 80% identity threshold
cd-hit-est -i rna_sequences.fasta -o clusters_0.8.fasta -c 0.8 -n 4 -M 2000 -T 4

# 60% identity threshold
cd-hit-est -i rna_sequences.fasta -o clusters_0.6.fasta -c 0.6 -n 4 -M 2000 -T 4

# 40% identity threshold
cd-hit-est -i rna_sequences.fasta -o clusters_0.4.fasta -c 0.4 -n 2 -M 2000 -T 4
```

### Script 4: Parse CD-HIT Results

```bash
python3 parse_cdhit_clusters.py --input <clstr_file> --output <csv_file> [--dist <dist_file>]
```

### Script 5: Calculate TM-score Matrix

```bash
python3 calculate_tm_matrix.py --pdb_dir <dir> --output <matrix.csv> [--checkpoint <file>] [--resume]
```

**Options:**
- `--checkpoint`: Save progress to checkpoint file
- `--resume`: Resume from checkpoint
- `--validate`: Validate existing matrix

**Important:** This step can take a long time for large datasets (hours to days). Use checkpointing!

### Script 6: Agglomerative Clustering

```bash
python3 agglomerative_clustering.py --input <tm_matrix.csv> --thresholds <t1 t2 t3> [--plot]
```

**Options:**
- `--thresholds`: Space-separated TM-score thresholds (e.g., 0.6 0.5 0.4)
- `--method`: Linkage method (single/complete/average/ward)
- `--plot`: Generate dendrogram and distribution plots

## Output Files

After running the pipeline, you'll have the following directory structure:

```
results/
├── rna_chains/                              # Extracted RNA chain PDB files
│   ├── 1ABC_A.pdb
│   ├── 1ABC_B.pdb
│   ├── ...
│   └── rna_chains_metadata.json             # Metadata for all chains
│
├── rna_sequences.fasta                      # All RNA sequences in FASTA format
│
├── Sequence Clustering Results:
│   ├── clusters_seq_0.8.fasta               # CD-HIT clustered sequences (80%)
│   ├── clusters_seq_0.8.fasta.clstr         # CD-HIT cluster file
│   ├── sequence_clusters_0.8.csv            # Sequence→Cluster mapping
│   ├── sequence_clusters_0.8_distribution.csv
│   ├── (similar files for 0.6 and 0.4)
│
├── Structure Clustering Results:
│   ├── tm_score_matrix.csv                  # Pairwise TM-score matrix
│   ├── structure_clusters_0.6.csv           # Structure→Cluster mapping
│   ├── structure_clusters_0.6_distribution.csv
│   ├── (similar files for 0.5 and 0.4)
│   ├── structure_clusters_comparison.csv    # Compare thresholds
│   ├── structure_clusters_dendrogram.png    # Hierarchical clustering tree
│   └── structure_clusters_size_distribution.png
│
└── clustering_summary_report.txt            # Comprehensive summary report
```

### Key Output Files

1. **sequence_clusters_X.csv**
   - Mapping of sequence IDs to cluster IDs
   - Format: `sequence_id, cluster_id`

2. **structure_clusters_X.csv**
   - Mapping of structure IDs to cluster IDs
   - Format: `structure_id, cluster_id`

3. **tm_score_matrix.csv**
   - N×N matrix of pairwise TM-scores
   - Can be very large for many structures

4. **clustering_summary_report.txt**
   - Summary of all clustering results
   - Statistics for each threshold

## Performance Considerations

### For Large Datasets (>500 structures)

1. **TM-score calculation is the bottleneck:**
   - For N structures, requires N×(N-1)/2 comparisons
   - Example: 1000 structures = ~500,000 comparisons
   - Time estimate: 2-10 seconds per comparison = 28-139 hours

2. **Use checkpointing:**
   ```bash
   python3 calculate_tm_matrix.py \
       --pdb_dir ./rna_chains \
       --output tm_score_matrix.csv \
       --checkpoint checkpoint.csv \
       --resume
   ```

3. **Parallel processing (advanced):**
   - Modify `calculate_tm_matrix.py` to use multiprocessing
   - Split structures into batches

4. **Test on subset first:**
   ```bash
   # Test with first 50 structures
   mkdir test_subset
   ls rna_chains/*.pdb | head -50 | xargs -I {} cp {} test_subset/
   python3 calculate_tm_matrix.py --pdb_dir test_subset --output test_matrix.csv
   ```

## Troubleshooting

### Problem: "USalign: command not found"

**Solution:**
```bash
# Add to PATH
export PATH=$PATH:~/bin

# Or add permanently to ~/.bashrc
echo 'export PATH=$PATH:~/bin' >> ~/.bashrc
source ~/.bashrc
```

### Problem: CD-HIT runs out of memory

**Solution:**
```bash
# Increase memory limit
cd-hit-est -i input.fasta -o output.fasta -c 0.8 -n 4 -M 4000 -T 4
```

### Problem: TM-score calculation too slow

**Solution:**
```bash
# Use checkpoint and run in background
nohup python3 calculate_tm_matrix.py \
    --pdb_dir ./rna_chains \
    --output tm_score_matrix.csv \
    --checkpoint checkpoint.csv > tm_calc.log 2>&1 &

# Check progress
tail -f tm_calc.log
```

### Problem: "No module named 'Bio'"

**Solution:**
```bash
pip install biopython
```

### Problem: Structure clustering fails

**Solution:**
```bash
# Validate TM-score matrix first
python3 calculate_tm_matrix.py --validate --output tm_score_matrix.csv

# Check for NaN or invalid values
python3 -c "import pandas as pd; import numpy as np; df = pd.read_csv('tm_score_matrix.csv', index_col=0); print(f'NaN values: {df.isna().sum().sum()}'); print(f'Out of range: {((df < 0) | (df > 1)).sum().sum()}')"
```

## Example Workflow

Here's a complete example with a small test dataset:

```bash
# 1. Setup
mkdir rna_test_project
cd rna_test_project
mkdir pdb_files

# 2. Download some test PDB files (example)
# (Replace with your actual PDB files)
cd pdb_files
wget https://files.rcsb.org/download/1EHZ.pdb
wget https://files.rcsb.org/download/1F7U.pdb
# ... add more PDB files
cd ..

# 3. Run pipeline
../run_clustering.sh pdb_files results

# 4. View results
cat results/clustering_summary_report.txt
ls -lh results/
```

## Advanced Options

### Custom Threshold Sets

```bash
# Test many sequence thresholds
for thresh in 0.9 0.8 0.7 0.6 0.5 0.4 0.3; do
    cd-hit-est -i rna_sequences.fasta -o clusters_${thresh}.fasta \
        -c $thresh -n 4 -M 2000 -T 4
done

# Test many structure thresholds
python3 agglomerative_clustering.py \
    --input tm_score_matrix.csv \
    --thresholds 0.7 0.65 0.6 0.55 0.5 0.45 0.4 0.35 0.3 \
    --plot
```

### Different Linkage Methods

```bash
# Try different methods
for method in single complete average ward; do
    python3 agglomerative_clustering.py \
        --input tm_score_matrix.csv \
        --thresholds 0.6 0.5 0.4 \
        --method $method \
        --output_prefix clusters_${method}
done
```

## References

1. **CD-HIT**: Fu et al. (2012) "CD-HIT: accelerated for clustering the next-generation sequencing data"
   - GitHub: https://github.com/weizhongli/cdhit

2. **US-align**: Zhang et al. (2022) "US-align: universal structure alignments of proteins, nucleic acids, and macromolecular complexes"
   - GitHub: https://github.com/pylelab/USalign

3. **TM-score**: Zhang & Skolnick (2005) "TM-align: a protein structure alignment algorithm based on the TM-score"

4. **RiboDiffusion**: Huang et al. (2024) "RiboDiffusion: tertiary structure-based RNA inverse folding with generative diffusion models"

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@article{huang2024ribodiffusion,
  title={RiboDiffusion: tertiary structure-based RNA inverse folding with generative diffusion models},
  author={Huang, Han and others},
  journal={Nature Computational Science},
  year={2024}
}
```

## License

This pipeline combines multiple tools with different licenses:
- CD-HIT: GPLv2
- US-align: Custom (free for academic use)
- Python scripts: MIT License

## Contact & Support

For issues and questions:
1. Check the [Troubleshooting](#troubleshooting) section
2. Verify all dependencies are correctly installed
3. Test on a small subset first

---

**Last Updated:** November 2025
