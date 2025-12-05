# RNA PDB Data Clustering - Technical Workflow
## Method 2: Dual Clustering Approach (Sequence + Structure Similarity)

This guide provides a complete technical workflow for clustering RNA PDB structures using both sequence similarity (PSI-CD-HIT) and structure similarity (US-align + agglomerative clustering).

---

## Table of Contents
1. [Prerequisites](#prerequisites)
2. [Software Installation](#software-installation)
3. [Data Preparation](#data-preparation)
4. [Sequence Similarity Clustering](#sequence-similarity-clustering)
5. [Structure Similarity Clustering](#structure-similarity-clustering)
6. [Complete Python Pipeline](#complete-python-pipeline)

---

## Prerequisites

### System Requirements
- Linux/Unix or macOS
- Python 3.7+
- 4GB+ RAM
- C++ compiler (gcc/g++)

### Python Packages
```bash
pip install biopython numpy scipy pandas matplotlib
```

---

## Software Installation

### 1. Install PSI-CD-HIT (for sequence clustering)

**Option A: Using conda (recommended)**
```bash
conda install -c bioconda cd-hit
```

**Option B: From source**
```bash
# Download CD-HIT
cd /tmp
git clone https://github.com/weizhongli/cdhit.git
cd cdhit

# Compile
make

# Install (may need sudo)
sudo make install

# Or install locally
make
mkdir -p ~/bin
cp cd-hit cd-hit-est psi-cd-hit ~/bin/
export PATH=$PATH:~/bin
```

**Verify installation:**
```bash
psi-cd-hit.pl -h
```

### 2. Install US-align (for structure alignment)

```bash
# Download US-align
cd /tmp
git clone https://github.com/pylelab/USalign.git
cd USalign

# Compile
g++ -static -O3 -ffast-math -lm -o USalign USalign.cpp

# Install
sudo cp USalign /usr/local/bin/
# Or for local installation
cp USalign ~/bin/
export PATH=$PATH:~/bin
```

**Verify installation:**
```bash
USalign -h
```

### 3. Install BioPython (if not already installed)

```bash
pip install biopython
```

---

## Data Preparation

### Step 1: Extract RNA chains from PDB files

Create `extract_rna_chains.py`:

```python
import os
from Bio import PDB
from Bio.PDB import PDBIO
import warnings
warnings.filterwarnings('ignore')

def extract_rna_chains(pdb_dir, output_dir):
    """
    Extract RNA chains from PDB files
    
    Args:
        pdb_dir: Directory containing PDB files
        output_dir: Directory to save extracted RNA chains
    """
    os.makedirs(output_dir, exist_ok=True)
    parser = PDB.PDBParser(QUIET=True)
    io = PDBIO()
    
    rna_structures = []
    
    for pdb_file in os.listdir(pdb_dir):
        if not pdb_file.endswith('.pdb'):
            continue
            
        pdb_id = pdb_file.replace('.pdb', '')
        pdb_path = os.path.join(pdb_dir, pdb_file)
        
        try:
            structure = parser.get_structure(pdb_id, pdb_path)
            
            for model in structure:
                for chain in model:
                    # Check if chain contains RNA
                    is_rna = False
                    for residue in chain:
                        if residue.get_resname().strip() in ['A', 'U', 'G', 'C']:
                            is_rna = True
                            break
                    
                    if is_rna:
                        # Extract sequence length
                        seq_length = sum(1 for res in chain if res.get_resname().strip() in ['A', 'U', 'G', 'C'])
                        
                        # Filter by length (20-280 nucleotides)
                        if 20 <= seq_length <= 280:
                            chain_id = f"{pdb_id}_{chain.id}"
                            output_file = os.path.join(output_dir, f"{chain_id}.pdb")
                            
                            io.set_structure(chain)
                            io.save(output_file)
                            
                            rna_structures.append({
                                'id': chain_id,
                                'pdb_id': pdb_id,
                                'chain': chain.id,
                                'length': seq_length,
                                'file': output_file
                            })
                            print(f"Extracted: {chain_id} (length: {seq_length})")
        
        except Exception as e:
            print(f"Error processing {pdb_file}: {str(e)}")
    
    return rna_structures

if __name__ == "__main__":
    pdb_directory = "./pdb_files"
    output_directory = "./rna_chains"
    
    structures = extract_rna_chains(pdb_directory, output_directory)
    print(f"\nTotal RNA chains extracted: {len(structures)}")
```

### Step 2: Extract sequences to FASTA format

Create `extract_sequences.py`:

```python
from Bio import PDB
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os

def extract_rna_sequences(pdb_dir, output_fasta):
    """
    Extract RNA sequences from PDB files and save to FASTA
    
    Args:
        pdb_dir: Directory containing RNA PDB files
        output_fasta: Output FASTA file path
    """
    parser = PDB.PDBParser(QUIET=True)
    sequences = []
    
    for pdb_file in sorted(os.listdir(pdb_dir)):
        if not pdb_file.endswith('.pdb'):
            continue
        
        chain_id = pdb_file.replace('.pdb', '')
        pdb_path = os.path.join(pdb_dir, pdb_file)
        
        try:
            structure = parser.get_structure(chain_id, pdb_path)
            
            # Extract sequence
            seq_list = []
            for model in structure:
                for chain in model:
                    for residue in chain:
                        resname = residue.get_resname().strip()
                        if resname in ['A', 'U', 'G', 'C']:
                            seq_list.append(resname)
            
            if seq_list:
                sequence = ''.join(seq_list)
                record = SeqRecord(
                    Seq(sequence),
                    id=chain_id,
                    description=f"RNA sequence from {chain_id}"
                )
                sequences.append(record)
        
        except Exception as e:
            print(f"Error processing {pdb_file}: {str(e)}")
    
    # Write to FASTA
    SeqIO.write(sequences, output_fasta, "fasta")
    print(f"Extracted {len(sequences)} sequences to {output_fasta}")

if __name__ == "__main__":
    pdb_directory = "./rna_chains"
    output_file = "rna_sequences.fasta"
    
    extract_rna_sequences(pdb_directory, output_file)
```

---

## Sequence Similarity Clustering

### Using PSI-CD-HIT

PSI-CD-HIT clusters sequences based on sequence identity. We'll test three thresholds: 0.8, 0.6, and 0.4.

**Command:**
```bash
# Clustering with 80% identity threshold
psi-cd-hit.pl -i rna_sequences.fasta -o clusters_seq_0.8 -c 0.8

# Clustering with 60% identity threshold
psi-cd-hit.pl -i rna_sequences.fasta -o clusters_seq_0.6 -c 0.6

# Clustering with 40% identity threshold
psi-cd-hit.pl -i rna_sequences.fasta -o clusters_seq_0.4 -c 0.4
```

**Alternative using CD-HIT-EST (for RNA):**
```bash
# 80% identity
cd-hit-est -i rna_sequences.fasta -o clusters_seq_0.8.fasta -c 0.8 -n 4 -M 2000 -T 4

# 60% identity
cd-hit-est -i rna_sequences.fasta -o clusters_seq_0.6.fasta -c 0.6 -n 4 -M 2000 -T 4

# 40% identity
cd-hit-est -i rna_sequences.fasta -o clusters_seq_0.4.fasta -c 0.4 -n 3 -M 2000 -T 4
```

**Parameters:**
- `-i`: Input FASTA file
- `-o`: Output file
- `-c`: Sequence identity threshold (0.0-1.0)
- `-n`: Word size (4 for ≥0.6, 3 for ≥0.5, 2 for ≥0.4)
- `-M`: Memory limit in MB
- `-T`: Number of threads

### Parse CD-HIT results

Create `parse_cdhit_clusters.py`:

```python
def parse_cdhit_clusters(clstr_file):
    """
    Parse CD-HIT cluster file (.clstr)
    
    Args:
        clstr_file: Path to .clstr file
    
    Returns:
        Dictionary mapping sequence IDs to cluster IDs
    """
    clusters = {}
    current_cluster = -1
    
    with open(clstr_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>Cluster'):
                current_cluster = int(line.split()[1])
            elif line:
                # Extract sequence ID
                seq_id = line.split('>')[1].split('...')[0]
                clusters[seq_id] = current_cluster
    
    return clusters

def save_cluster_mapping(clusters, output_file):
    """Save cluster mapping to file"""
    import pandas as pd
    
    df = pd.DataFrame([
        {'sequence_id': seq_id, 'cluster_id': cluster_id}
        for seq_id, cluster_id in clusters.items()
    ])
    
    df.to_csv(output_file, index=False)
    print(f"Saved cluster mapping to {output_file}")
    print(f"Total clusters: {df['cluster_id'].nunique()}")
    print(f"Total sequences: {len(df)}")

if __name__ == "__main__":
    # Parse different thresholds
    for threshold in ['0.8', '0.6', '0.4']:
        clstr_file = f"clusters_seq_{threshold}.fasta.clstr"
        output_file = f"sequence_clusters_{threshold}.csv"
        
        clusters = parse_cdhit_clusters(clstr_file)
        save_cluster_mapping(clusters, output_file)
```

---

## Structure Similarity Clustering

### Step 1: Calculate TM-score matrix using US-align

Create `calculate_tm_matrix.py`:

```python
import os
import subprocess
import numpy as np
import pandas as pd
from tqdm import tqdm

def calculate_tm_score(pdb1, pdb2):
    """
    Calculate TM-score between two PDB structures using US-align
    
    Args:
        pdb1: Path to first PDB file
        pdb2: Path to second PDB file
    
    Returns:
        TM-score (float)
    """
    try:
        result = subprocess.run(
            ['USalign', pdb1, pdb2],
            capture_output=True,
            text=True,
            timeout=60
        )
        
        # Parse TM-score from output
        for line in result.stdout.split('\n'):
            if 'TM-score=' in line and 'Chain_1' in line:
                tm_score = float(line.split('TM-score=')[1].split()[0])
                return tm_score
        
        return 0.0
    
    except Exception as e:
        print(f"Error calculating TM-score: {str(e)}")
        return 0.0

def calculate_tm_matrix(pdb_dir, output_file):
    """
    Calculate pairwise TM-score matrix for all PDB files
    
    Args:
        pdb_dir: Directory containing PDB files
        output_file: Output file for TM-score matrix
    """
    # Get list of PDB files
    pdb_files = sorted([f for f in os.listdir(pdb_dir) if f.endswith('.pdb')])
    structure_ids = [f.replace('.pdb', '') for f in pdb_files]
    n = len(pdb_files)
    
    print(f"Calculating TM-scores for {n} structures...")
    
    # Initialize matrix
    tm_matrix = np.zeros((n, n))
    
    # Calculate pairwise TM-scores
    total_pairs = n * (n - 1) // 2
    with tqdm(total=total_pairs) as pbar:
        for i in range(n):
            tm_matrix[i, i] = 1.0  # Self-alignment
            
            for j in range(i + 1, n):
                pdb1 = os.path.join(pdb_dir, pdb_files[i])
                pdb2 = os.path.join(pdb_dir, pdb_files[j])
                
                tm_score = calculate_tm_score(pdb1, pdb2)
                tm_matrix[i, j] = tm_score
                tm_matrix[j, i] = tm_score  # Symmetric
                
                pbar.update(1)
    
    # Save matrix
    df = pd.DataFrame(tm_matrix, index=structure_ids, columns=structure_ids)
    df.to_csv(output_file)
    print(f"TM-score matrix saved to {output_file}")
    
    return tm_matrix, structure_ids

if __name__ == "__main__":
    pdb_directory = "./rna_chains"
    output_matrix = "tm_score_matrix.csv"
    
    tm_matrix, ids = calculate_tm_matrix(pdb_directory, output_matrix)
```

### Step 2: Perform agglomerative clustering

Create `agglomerative_clustering.py`:

```python
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram

def perform_structure_clustering(tm_matrix_file, thresholds=[0.6, 0.5, 0.4]):
    """
    Perform agglomerative clustering on TM-score matrix
    
    Args:
        tm_matrix_file: Path to TM-score matrix CSV
        thresholds: List of TM-score thresholds for clustering
    
    Returns:
        Dictionary of cluster assignments for each threshold
    """
    # Load TM-score matrix
    df = pd.read_csv(tm_matrix_file, index_col=0)
    structure_ids = df.index.tolist()
    tm_matrix = df.values
    
    # Convert similarity to distance (1 - TM-score)
    distance_matrix = 1 - tm_matrix
    
    # Convert to condensed distance matrix for scipy
    condensed_dist = squareform(distance_matrix, checks=False)
    
    # Perform hierarchical clustering
    print("Performing agglomerative clustering...")
    linkage_matrix = linkage(condensed_dist, method='average')
    
    # Cluster at different thresholds
    results = {}
    
    for threshold in thresholds:
        # Convert TM-score threshold to distance threshold
        distance_threshold = 1 - threshold
        
        # Get cluster labels
        cluster_labels = fcluster(linkage_matrix, distance_threshold, criterion='distance')
        
        # Create mapping
        cluster_dict = {
            structure_ids[i]: cluster_labels[i]
            for i in range(len(structure_ids))
        }
        
        results[threshold] = cluster_dict
        
        # Save results
        output_file = f"structure_clusters_{threshold}.csv"
        df_out = pd.DataFrame([
            {'structure_id': sid, 'cluster_id': cid}
            for sid, cid in cluster_dict.items()
        ])
        df_out.to_csv(output_file, index=False)
        
        n_clusters = len(set(cluster_labels))
        print(f"Threshold {threshold}: {n_clusters} clusters")
    
    # Plot dendrogram (optional)
    plt.figure(figsize=(12, 6))
    dendrogram(linkage_matrix, no_labels=True)
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('Structure Index')
    plt.ylabel('Distance (1 - TM-score)')
    plt.savefig('clustering_dendrogram.png', dpi=300, bbox_inches='tight')
    print("Dendrogram saved to clustering_dendrogram.png")
    
    return results, linkage_matrix

if __name__ == "__main__":
    tm_matrix_file = "tm_score_matrix.csv"
    thresholds = [0.6, 0.5, 0.4]
    
    clusters, linkage_mat = perform_structure_clustering(tm_matrix_file, thresholds)
```

---

## Complete Python Pipeline

Create `complete_clustering_pipeline.py`:

```python
"""
Complete RNA PDB Clustering Pipeline
Performs both sequence and structure similarity clustering
"""

import os
import sys
import argparse
from extract_rna_chains import extract_rna_chains
from extract_sequences import extract_rna_sequences
from parse_cdhit_clusters import parse_cdhit_clusters, save_cluster_mapping
from calculate_tm_matrix import calculate_tm_matrix
from agglomerative_clustering import perform_structure_clustering
import subprocess
import pandas as pd

class RNAClusteringPipeline:
    def __init__(self, pdb_dir, output_dir, seq_thresholds=[0.8, 0.6, 0.4], 
                 struct_thresholds=[0.6, 0.5, 0.4]):
        self.pdb_dir = pdb_dir
        self.output_dir = output_dir
        self.seq_thresholds = seq_thresholds
        self.struct_thresholds = struct_thresholds
        
        os.makedirs(output_dir, exist_ok=True)
    
    def run_full_pipeline(self):
        """Execute complete clustering pipeline"""
        
        print("="*60)
        print("RNA PDB Clustering Pipeline")
        print("="*60)
        
        # Step 1: Extract RNA chains
        print("\n[Step 1] Extracting RNA chains from PDB files...")
        rna_chains_dir = os.path.join(self.output_dir, "rna_chains")
        structures = extract_rna_chains(self.pdb_dir, rna_chains_dir)
        print(f"Extracted {len(structures)} RNA chains")
        
        # Step 2: Extract sequences to FASTA
        print("\n[Step 2] Extracting sequences to FASTA...")
        fasta_file = os.path.join(self.output_dir, "rna_sequences.fasta")
        extract_rna_sequences(rna_chains_dir, fasta_file)
        
        # Step 3: Sequence clustering
        print("\n[Step 3] Performing sequence clustering...")
        self.run_sequence_clustering(fasta_file)
        
        # Step 4: Calculate TM-score matrix
        print("\n[Step 4] Calculating TM-score matrix...")
        tm_matrix_file = os.path.join(self.output_dir, "tm_score_matrix.csv")
        tm_matrix, structure_ids = calculate_tm_matrix(rna_chains_dir, tm_matrix_file)
        
        # Step 5: Structure clustering
        print("\n[Step 5] Performing structure clustering...")
        os.chdir(self.output_dir)
        clusters, linkage_mat = perform_structure_clustering(
            "tm_score_matrix.csv", 
            self.struct_thresholds
        )
        
        # Step 6: Generate summary report
        print("\n[Step 6] Generating summary report...")
        self.generate_summary_report()
        
        print("\n" + "="*60)
        print("Pipeline completed successfully!")
        print(f"Results saved to: {self.output_dir}")
        print("="*60)
    
    def run_sequence_clustering(self, fasta_file):
        """Run CD-HIT sequence clustering"""
        for threshold in self.seq_thresholds:
            output_prefix = os.path.join(self.output_dir, f"clusters_seq_{threshold}")
            
            # Determine word size based on threshold
            if threshold >= 0.6:
                word_size = 4
            elif threshold >= 0.5:
                word_size = 3
            else:
                word_size = 2
            
            cmd = [
                'cd-hit-est',
                '-i', fasta_file,
                '-o', f"{output_prefix}.fasta",
                '-c', str(threshold),
                '-n', str(word_size),
                '-M', '2000',
                '-T', '4'
            ]
            
            print(f"  Running CD-HIT with threshold {threshold}...")
            try:
                subprocess.run(cmd, check=True, capture_output=True)
                
                # Parse results
                clstr_file = f"{output_prefix}.fasta.clstr"
                clusters = parse_cdhit_clusters(clstr_file)
                save_cluster_mapping(
                    clusters, 
                    os.path.join(self.output_dir, f"sequence_clusters_{threshold}.csv")
                )
            except subprocess.CalledProcessError as e:
                print(f"  Error running CD-HIT: {e}")
    
    def generate_summary_report(self):
        """Generate summary report of clustering results"""
        report_file = os.path.join(self.output_dir, "clustering_summary.txt")
        
        with open(report_file, 'w') as f:
            f.write("RNA PDB Clustering Summary Report\n")
            f.write("="*60 + "\n\n")
            
            # Sequence clustering results
            f.write("Sequence Clustering Results:\n")
            f.write("-"*60 + "\n")
            for threshold in self.seq_thresholds:
                csv_file = os.path.join(self.output_dir, f"sequence_clusters_{threshold}.csv")
                if os.path.exists(csv_file):
                    df = pd.read_csv(csv_file)
                    n_clusters = df['cluster_id'].nunique()
                    n_sequences = len(df)
                    f.write(f"Threshold {threshold}: {n_clusters} clusters, {n_sequences} sequences\n")
            
            f.write("\n")
            
            # Structure clustering results
            f.write("Structure Clustering Results:\n")
            f.write("-"*60 + "\n")
            for threshold in self.struct_thresholds:
                csv_file = os.path.join(self.output_dir, f"structure_clusters_{threshold}.csv")
                if os.path.exists(csv_file):
                    df = pd.read_csv(csv_file)
                    n_clusters = df['cluster_id'].nunique()
                    n_structures = len(df)
                    f.write(f"Threshold {threshold}: {n_clusters} clusters, {n_structures} structures\n")
        
        print(f"Summary report saved to {report_file}")

def main():
    parser = argparse.ArgumentParser(
        description='RNA PDB Clustering Pipeline'
    )
    parser.add_argument(
        '--pdb_dir',
        required=True,
        help='Directory containing PDB files'
    )
    parser.add_argument(
        '--output_dir',
        required=True,
        help='Output directory for results'
    )
    parser.add_argument(
        '--seq_thresholds',
        nargs='+',
        type=float,
        default=[0.8, 0.6, 0.4],
        help='Sequence similarity thresholds'
    )
    parser.add_argument(
        '--struct_thresholds',
        nargs='+',
        type=float,
        default=[0.6, 0.5, 0.4],
        help='Structure similarity thresholds (TM-score)'
    )
    
    args = parser.parse_args()
    
    pipeline = RNAClusteringPipeline(
        args.pdb_dir,
        args.output_dir,
        args.seq_thresholds,
        args.struct_thresholds
    )
    
    pipeline.run_full_pipeline()

if __name__ == "__main__":
    main()
```

---

## Usage Instructions

### Quick Start

1. **Prepare your data:**
   ```bash
   # Create directory structure
   mkdir -p rna_clustering_project
   cd rna_clustering_project
   mkdir pdb_files
   
   # Place your PDB files in pdb_files/
   ```

2. **Run the complete pipeline:**
   ```bash
   python complete_clustering_pipeline.py \
       --pdb_dir ./pdb_files \
       --output_dir ./results \
       --seq_thresholds 0.8 0.6 0.4 \
       --struct_thresholds 0.6 0.5 0.4
   ```

### Step-by-Step Execution

If you prefer to run each step separately:

```bash
# Step 1: Extract RNA chains
python extract_rna_chains.py

# Step 2: Extract sequences
python extract_sequences.py

# Step 3: Sequence clustering
cd-hit-est -i rna_sequences.fasta -o clusters_seq_0.8.fasta -c 0.8 -n 4 -M 2000 -T 4
python parse_cdhit_clusters.py

# Step 4: Calculate TM-score matrix
python calculate_tm_matrix.py

# Step 5: Structure clustering
python agglomerative_clustering.py
```

---

## Output Files

After running the pipeline, you'll have:

```
results/
├── rna_chains/                    # Extracted RNA chain PDB files
│   ├── 1ABC_A.pdb
│   ├── 1ABC_B.pdb
│   └── ...
├── rna_sequences.fasta            # RNA sequences in FASTA format
├── clusters_seq_0.8.fasta         # CD-HIT output (80% threshold)
├── clusters_seq_0.8.fasta.clstr   # CD-HIT cluster file
├── sequence_clusters_0.8.csv      # Sequence cluster mapping
├── sequence_clusters_0.6.csv
├── sequence_clusters_0.4.csv
├── tm_score_matrix.csv            # Pairwise TM-score matrix
├── structure_clusters_0.6.csv     # Structure cluster mapping
├── structure_clusters_0.5.csv
├── structure_clusters_0.4.csv
├── clustering_dendrogram.png      # Visualization
└── clustering_summary.txt         # Summary report
```

---

## Troubleshooting

### Common Issues

1. **US-align not found:**
   ```bash
   export PATH=$PATH:~/bin
   # Or add to ~/.bashrc
   echo 'export PATH=$PATH:~/bin' >> ~/.bashrc
   source ~/.bashrc
   ```

2. **CD-HIT memory error:**
   - Increase `-M` parameter: `-M 4000` (for 4GB)

3. **TM-score calculation timeout:**
   - Increase timeout in `calculate_tm_score()` function
   - Run on fewer structures for testing

4. **Python package missing:**
   ```bash
   pip install biopython numpy scipy pandas matplotlib tqdm
   ```

---

## Performance Optimization

For large datasets (>1000 structures):

1. **Parallelize TM-score calculation:**
   Use multiprocessing to calculate TM-scores in parallel

2. **Use US-align batch mode:**
   ```bash
   USalign -dir pdb_directory/ -suffix .pdb -outfmt 2 > tm_scores.txt
   ```

3. **Subsample for initial testing:**
   Test on a small subset first (50-100 structures)

---

## References

- CD-HIT: https://github.com/weizhongli/cdhit
- US-align: https://zhanggroup.org/US-align/
- Scipy agglomerative clustering: https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html
