# PDB to 2D Base Pairing Analysis

A simple Python tool to analyze 3D nucleic acid structures from PDB files and generate 2D base pairing matrices based on spatial proximity and canonical base pairing rules.

## Overview

This tool processes PDB files containing nucleic acid structures (RNA/DNA) and identifies potential base pairs based on:
- **Sequence distance**: Residues must be at least 5 positions apart (|j-i| > 4)
- **Base pairing rules**: Canonical Watson-Crick pairs (A-U, U-A, G-C, C-G) and wobble pairs (G-U, U-G)
- **Spatial distance**: C4' atoms of residues within 10Å of each other

The output is a 2D matrix showing which residues form base pairs.

## Features

- Parse PDB files to extract nucleic acid residues
- Calculate C4' atom distances between all residue pairs
- Apply canonical base pairing rules (G-C, A-U, G-U)
- Generate 2D base pairing matrix
- Simple, efficient implementation without over-engineering

## Requirements

```
numpy
biopython  # for PDB parsing
matplotlib  # for visualization (optional)
```

## Installation

```bash
pip install numpy biopython matplotlib
```

## Usage

### Basic Usage

```python
from pdb_basepair import analyze_base_pairs

# Analyze a PDB file
matrix, residue_info = analyze_base_pairs('structure.pdb')

# Display results
print("Base pairing matrix:")
print(matrix)
```

### Advanced Usage

```python
from pdb_basepair import PDBBasePairAnalyzer

# Initialize analyzer
analyzer = PDBBasePairAnalyzer(
    distance_cutoff=10.0,  # Angstroms
    chains=['A', 'B']      # Specific chains to analyze
)

# Load and analyze PDB
analyzer.load_pdb('structure.pdb')
matrix = analyzer.get_base_pair_matrix()

# Get detailed pairing information
pairs = analyzer.get_base_pairs()
for pair in pairs:
    print(f"Residue {pair['res1']} paired with {pair['res2']} (distance: {pair['distance']:.2f}Å)")

# Visualize matrix
analyzer.plot_matrix(save_path='basepair_matrix.png')
```

## Implementation Details

### Core Algorithm

1. **PDB Parsing**: Extract nucleic acid residues and their C4' atom coordinates
2. **Sequence Filtering**: Only consider residue pairs where |j-i| > 4 (at least 5 positions apart)
3. **Distance Calculation**: Compute pairwise distances between filtered C4' atoms
4. **Base Pair Filtering**: Apply distance cutoff (≤10Å) and canonical pairing rules
5. **Matrix Generation**: Create symmetric boolean matrix indicating base pairs

### Base Pairing Rules

Residues i and j are considered paired only if they satisfy ALL three criteria:

1. **Sequence Distance**: |j-i| > 4 (residues must be at least 5 positions apart)
2. **Base Compatibility**: Only the following pairs are allowed:
   - **A-U**: Standard Watson-Crick pair (2 hydrogen bonds)
   - **U-A**: Standard Watson-Crick pair (2 hydrogen bonds)
   - **G-C**: Standard Watson-Crick pair (3 hydrogen bonds)
   - **C-G**: Standard Watson-Crick pair (3 hydrogen bonds)
   - **G-U**: Wobble pair (2 hydrogen bonds)
   - **U-G**: Wobble pair (2 hydrogen bonds)
3. **Spatial Distance**: C4' atoms must be within 10Å of each other

All other combinations are excluded.

### Data Structures

**Input**: PDB file with nucleic acid structure
**Output**: 
- `matrix`: NxN boolean array where `matrix[i,j] = True` indicates base pair between residues i and j
- `residue_info`: List of dictionaries containing residue metadata (chain, number, base type)

## File Structure

```
pdb-basepair-analysis/
├── README.md
├── pdb_basepair.py          # Main implementation
├── examples/
│   ├── example_usage.py     # Usage examples
│   └── sample_output.png    # Sample matrix visualization
├── tests/
│   └── test_basepair.py     # Unit tests
└── requirements.txt
```

## Example Output

For a structure with 20 nucleotides, the output matrix might look like:

```
     A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T
A  [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 ]
B  [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0 ]
...
T  [ 1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ]
```

Where 1 indicates a valid base pair and 0 indicates no pairing.

## Design Philosophy

This implementation prioritizes:
- **Simplicity**: Clear, readable code without unnecessary complexity
- **Efficiency**: Fast distance calculations using vectorized operations
- **Flexibility**: Configurable parameters for different analysis needs
- **Reliability**: Robust PDB parsing and error handling

## Limitations

- Only considers C4' atom distances (not full geometric criteria)
- Limited to canonical base pairs (A-U, U-A, G-C, C-G, G-U, U-G)
- Requires minimum sequence separation of 5 residues (excludes local interactions)
- Does not account for tertiary structure constraints beyond distance
- Assumes standard nucleotide naming conventions in PDB files

## Future Enhancements

- Support for modified bases
- Geometric angle criteria for base pairing
- Export to common formats (CSV, JSON)
- Integration with structure visualization tools
- Batch processing capabilities

## Contributing

This is designed as a simple, focused tool. When contributing:
- Maintain the principle of simplicity
- Add comprehensive tests for new features
- Update documentation accordingly
- Follow PEP 8 style guidelines
