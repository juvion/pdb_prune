# RNA Structure Processing Pipeline

A comprehensive toolkit for processing and analyzing RNA structures from the Protein Data Bank (PDB).

## Overview

This pipeline provides tools for:
- Downloading and processing RNA structures from PDB
- Converting between PDB and NumPy formats
- Extracting RNA segments and loops
- Generating FASTA sequences
- Analyzing RNA structure properties

## Installation

### Prerequisites
- Python 3.8+
- pip (Python package installer)

### Dependencies
```bash
pip install -r requirements.txt
```

Required packages:
- biopython
- numpy
- pandas
- matplotlib
- jupyter
- requests

## Directory Structure

```
pdb_prune/
├── README.md
├── requirements.txt
├── tutorial.ipynb
├── utils/
│   ├── __init__.py
│   ├── download_pdb.py
│   ├── pdb_to_npy.py
│   ├── npy_to_pdb.py
│   ├── extract_loops.py
│   ├── extract_rna_segments.py
│   └── extract_sequence.py
├── data/
│   ├── raw_pdbs/
│   ├── processed_pdbs/
│   └── npy_files/
└── examples/
    └── output/
```

## Usage Examples

### 1. Downloading PDB Files
```python
from utils.download_pdb import PDBDownloader

# Download specific PDB IDs
downloader = PDBDownloader()
downloader.download_pdbs(['1ABC', '2XYZ'])
```

### 2. Converting PDB to NumPy
```python
from utils.pdb_to_npy import PDBToNumpyConverter

# Convert PDB files to NumPy arrays
converter = PDBToNumpyConverter(
    processed_dir="processed_pdbs",
    npy_dir="npy_files"
)
converter.convert_all_pdbs()
```

### 3. Extracting RNA Loops
```python
from utils.extract_loops import LoopExtractor

# Extract loops with 10Å cutoff
extractor = LoopExtractor(
    input_dir="processed_pdbs",
    output_dir="extracted_loops",
    distance_cutoff=10.0
)
extractor.process()
```

### 4. Extracting RNA Segments
```python
from utils.extract_rna_segments import RNAExtractor

# Extract segments of length 5-20
extractor = RNAExtractor(
    input_dir="processed_pdbs",
    generation_id="run1",
    min_length=5,
    max_length=20,
    num_extractions_per_chain=5
)
extractor.process()
```

### 5. Generating FASTA Files
```python
from utils.extract_sequence import SequenceExtractor

# Extract sequences from PDB files
extractor = SequenceExtractor(
    pdb_dir="processed_pdbs",
    output_dir="sequences"
)
extractor.process()
```

## Input/Output Formats

### PDB Files
- Input: Standard PDB format files
- Output: Processed PDB files with RNA chains only

### NumPy Arrays
- Shape: (sequence_length, 7, 3)
- Contains coordinates for: P, O5', C5', C4', C3', O3', and base connecting atom (N1/N9)

### FASTA Files
- Standard FASTA format
- Header format: `>PDB_ID_chain description`

## Common Use Cases

1. **RNA Structure Analysis**
   - Download RNA structures
   - Convert to NumPy arrays for analysis
   - Extract specific regions of interest

2. **Loop Detection**
   - Identify and extract RNA loops
   - Analyze loop properties
   - Generate loop datasets

3. **Sequence Analysis**
   - Extract RNA sequences
   - Generate FASTA files
   - Analyze sequence patterns

## Best Practices

1. **File Organization**
   - Keep raw and processed files separate
   - Use consistent naming conventions
   - Maintain backup copies of important data

2. **Error Handling**
   - Check file existence before processing
   - Validate input formats
   - Handle missing atoms gracefully

3. **Performance**
   - Process files in batches
   - Use appropriate data structures
   - Monitor memory usage

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details. 