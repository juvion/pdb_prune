# RNA Sequence Scoring Guide

This guide shows how to score a predicted RNA sequence against a reference PDB structure using two Python utilities.

## Overview

The workflow involves:
1. Extract base pair matrix from reference PDB structure
2. Score predicted sequence against reference using the base pair matrix

## Prerequisites

- Python 3.x
- Required Python packages: `numpy`, `biopython`, `editdistance`
- Input files:
  - Reference PDB file: `/Users/ju/Documents/Dev/pdb_prune/data/raw_data/pdbs/raw_pdbs_full_download/pdb1aqo.ent`
  - Backbone PDB file: `/Users/ju/Documents/Dev/pdb_prune/data/experiments_data/exp2.2_manuscript/extracted_pdbs/pdb1aqo_A.pdb`
  - Predicted sequence FASTA file: `pred_1aqo.fasta` (you need to create this)

## Step 1: Extract Base Pair Matrix

Extract the base pair matrix from the reference PDB structure:

```bash
python3 PATH_TO_SCRITPS/pdb_3dto2d_basepair.py \
  PATH_TO_RAW_PDB/pdb1aqo.ent \
  --distance-method c1_c1 \
  --no-coplanarity \
  -o pdb1aqo.mat
```

**Parameters explained:**
- `--distance-method c1_c1`: Use C1'-C1' distance criterion for base pair detection
- `--no-coplanarity`: Disable coplanarity check for more lenient base pair detection
- `-o pdb1aqo.mat`: Output matrix file name

**Expected output:**
```
Parsing PDB file: /Users/ju/Documents/Dev/pdb_prune/data/raw_data/pdbs/raw_pdbs_full_download/pdb1aqo.ent
Found X nucleotides
Detecting base pairs among X nucleotides using c1_c1 distance method...
Found Y base pairs
Matrix saved to pdb1aqo.mat
```

## Step 2: Extract Reference Sequence

Extract the reference sequence from the backbone PDB:

```bash
python3 utils/pdb_to_fasta.py -h
usage: pdb_to_fasta.py [-h] --pdb-dir PDB_DIR [--output-dir OUTPUT_DIR]

Extract RNA sequences from PDB files

options:
  -h, --help            show this help message and exit
  --pdb-dir PDB_DIR     Directory containing PDB files
  --output-dir OUTPUT_DIR
                        Directory to save FASTA files (default: sequences)
```

*Note: If `pdb_to_fasta.py` doesn't exist, you can manually create the reference sequence FASTA file based on the PDB structure.*

## Step 3: Create Predicted Sequence File

Create your predicted sequence FASTA file (`pred_1aqo.fasta`):

```
>predicted_sequence_1aqo
GGAGUGCUUCAACCAGUGCUUGGACGCUCC
```

*Replace the sequence above with your actual ML model prediction.*

## Step 4: Score Predicted Sequence

Score the predicted sequence against the reference:

```bash
python3 PATH_TO_SCRITPS/sequence_similarity_score.py \
  pred_1aqo.fasta \
  pdb1aqo_A.fasta \
  pdb1aqo.mat \
  --lambda_param 0.5 \
  --format summary
```

**Parameters explained:**
- `pred_1aqo.fasta`: Your predicted sequence FASTA file
- `pdb1aqo_A.fasta`: Reference sequence FASTA file
- `pdb1aqo.mat`: Base pair matrix from Step 1
- `--lambda_param 0.5`: Weight parameter (0.5 = equal weight for base pair and edit distance scores)
- `--format summary`: Output format (summary, detailed, json, csv)

**Expected output:**
```
Combined Score: 0.8500
Base Pair Score: 0.9200
Edit Distance Score: 0.7800
```

## Alternative Output Formats

### Detailed Output
```bash
python3 PATH_TO_SCRITPS/sequence_similarity_score.py \
  pred_1aqo.fasta pdb1aqo_A.fasta pdb1aqo.mat \
  --lambda_param 0.5 --format detailed
```

### JSON Output
```bash
python3 PATH_TO_SCRITPS/sequence_similarity_score.py \
  pred_1aqo.fasta pdb1aqo_A.fasta pdb1aqo.mat \
  --lambda_param 0.5 --format json --output results.json
```

## Troubleshooting

1. **File not found errors**: Ensure all file paths are correct and files exist
2. **Sequence length mismatch**: Predicted and reference sequences must have the same length
3. **Invalid RNA sequence**: Only A, U, C, G nucleotides are allowed
4. **Matrix dimension mismatch**: Matrix dimensions must match sequence length

## Score Interpretation

- **Combined Score**: Weighted combination of base pair and edit distance scores (0-1, higher is better)
- **Base Pair Score**: How well predicted sequence maintains reference base pairing (0-1, higher is better)
- **Edit Distance Score**: Sequence similarity (0-1, higher is better)
- **Lambda Parameter**: Controls weighting (0.7 = 70% base pair, 30% edit distance)

## Example Complete Workflow

```bash
# Step 1: Extract base pair matrix
python3 /Users/ju/Documents/Dev/pdb_prune/utils/pdb_3dto2d_basepair.py \
  /Users/ju/Documents/Dev/pdb_prune/data/raw_data/pdbs/raw_pdbs_full_download/pdb1aqo.ent \
  --distance-method c1_c1 --no-coplanarity -o pdb1aqo.mat

# Step 2: Score prediction (assuming you have pred_1aqo.fasta and pdb1aqo_A.fasta)
python3 /Users/ju/Documents/Dev/pdb_prune/utils/sequence_similarity_score.py \
  pred_1aqo.fasta pdb1aqo_A.fasta pdb1aqo.mat \
  --lambda_param 0.5 --format summary
```