# PDB RNA Structure Processor

A Python tool for processing RNA structures from PDB files, extracting sequences, and converting them to various formats.

## Features

- Process PDB files containing RNA structures
- Extract RNA sequences and save them as FASTA files
- Convert PDB structures to NumPy arrays for machine learning applications
- Handle both `.pdb` and `.ent` file formats
- Support for multiple RNA chains in a single structure
- Automatic handling of missing atoms with NaN values

## Requirements

- Python 3.6+
- Biopython
- NumPy
- Requests

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/pdb_prune.git
cd pdb_prune
```

2. Install the required packages:
```bash
pip install -r requirements.txt
```

## Usage

### 1. Processing PDB Files

Place your PDB files (`.pdb` or `.ent`) in the `original_pdbs` directory. The script will:
- Extract RNA chains
- Save processed structures in `processed_pdbs`
- Generate FASTA sequences in `rna_sequences`
- Create NumPy arrays in `npy_files`

Run the processing script:
```bash
python pdb_rna_processor.py
```

### 2. Converting to NumPy Arrays

The script will automatically convert processed PDB files to NumPy arrays with the following structure:
- Shape: `(sequence_length, 7, 3)`
- 7 atoms per residue: P, O5', C5', C4', C3', O3', and N1/N9
- Missing atoms are filled with NaN values

Run the conversion script:
```bash
python pdb_to_npy.py
```

## Directory Structure

```
pdb_prune/
├── original_pdbs/     # Input PDB files
├── processed_pdbs/    # Processed RNA structures
├── rna_sequences/     # FASTA sequence files
├── npy_files/        # NumPy array files
├── pdb_rna_processor.py
├── pdb_to_npy.py
└── requirements.txt
```

## File Formats

### Processed PDB Files
- Located in `processed_pdbs/`
- Named as `pdb{id}_{chain}.pdb`
- Contain only RNA chains with specified atoms

### FASTA Files
- Located in `rna_sequences/`
- Named as `pdb{id}_{chain}.fasta`
- Contain RNA sequences in FASTA format

### NumPy Arrays
- Located in `npy_files/`
- Named as `pdb{id}_{chain}.npy`
- Shape: `(sequence_length, 7, 3)`
- Atom order: P, O5', C5', C4', C3', O3', N1/N9

## Contributing

Feel free to submit issues and enhancement requests!

## License

This project is licensed under the MIT License - see the LICENSE file for details. 