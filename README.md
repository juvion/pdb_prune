# PDB RNA Processor

This script processes PDB files containing RNA structures using Biopython. It can:
1. Download PDB files (either randomly or by specific ID)
2. Isolate RNA chains
3. Extract specific atoms from each RNA chain
4. Output the data in PDB format

## Requirements

- Python 3.6 or higher
- Required packages (install using `pip install -r requirements.txt`):
  - biopython
  - requests

## Usage

1. Install the required dependencies:
```bash
pip install -r requirements.txt
```

2. Run the script:
```bash
python pdb_rna_processor.py
```

The script will:
- Create a `processed_pdbs` directory to store the output files
- Process an example PDB (1ehz) and 3 random PDBs
- For each PDB file, it will:
  - Download the structure
  - Extract RNA chains
  - Save each RNA chain as a separate PDB file with only the specified atoms

## Output

The processed PDB files will be saved in the `processed_pdbs` directory with the naming format:
`{original_pdb_id}_{chain_id}.pdb`

Each output file contains only the RNA chain with the following atoms:
- Backbone atoms: P, O5', C5', C4', C3', O3', C2', O2', C1'
- Base atoms: N1 (for U and C) or N9 (for A and G)

## Customization

You can modify the script to:
- Change the number of random PDBs to process
- Use different PDB IDs
- Modify the output directory
- Change the set of atoms to extract 