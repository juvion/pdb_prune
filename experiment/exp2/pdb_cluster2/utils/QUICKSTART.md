# Quick Start Guide - RNA PDB Clustering Pipeline

## 🚀 5-Minute Setup

### 1. Install Dependencies

```bash
chmod +x install_dependencies.sh
./install_dependencies.sh
source ~/.bashrc  # If ~/bin was added to PATH
```

### 2. Prepare Your Data

```bash
mkdir my_clustering_project
cd my_clustering_project
mkdir pdb_files

# Copy your PDB files to pdb_files/
cp /path/to/your/pdbs/*.pdb pdb_files/
```

### 3. Run the Pipeline

```bash
chmod +x run_clustering.sh
./run_clustering.sh pdb_files results
```

That's it! Results will be in the `results/` directory.

---

## 📊 What You'll Get

After running the pipeline, you'll have:

1. **Sequence Clusters** (using CD-HIT)
   - `sequence_clusters_0.8.csv` - 80% identity threshold
   - `sequence_clusters_0.6.csv` - 60% identity threshold  
   - `sequence_clusters_0.4.csv` - 40% identity threshold

2. **Structure Clusters** (using TM-score)
   - `structure_clusters_0.6.csv` - TM-score ≥ 0.6
   - `structure_clusters_0.5.csv` - TM-score ≥ 0.5
   - `structure_clusters_0.4.csv` - TM-score ≥ 0.4

3. **Visualizations**
   - `structure_clusters_dendrogram.png` - Hierarchical clustering tree
   - `structure_clusters_size_distribution.png` - Cluster size plots

4. **Summary Report**
   - `clustering_summary_report.txt` - Complete statistics

---

## 🔧 Customization

### Custom Thresholds

```bash
# Sequence: 0.9, 0.7, 0.5
# Structure: 0.7, 0.6, 0.5
./run_clustering.sh pdb_files results 0.9,0.7,0.5 0.7,0.6,0.5
```

### Step-by-Step Execution

```bash
# 1. Extract RNA chains
python3 extract_rna_chains.py --pdb_dir pdb_files --output_dir rna_chains

# 2. Extract sequences
python3 extract_sequences.py --pdb_dir rna_chains --output rna_sequences.fasta

# 3. Sequence clustering
cd-hit-est -i rna_sequences.fasta -o clusters_0.8.fasta -c 0.8 -n 4 -M 2000 -T 4
python3 parse_cdhit_clusters.py --input clusters_0.8.fasta.clstr --output seq_clusters_0.8.csv

# 4. Calculate TM-scores (slow!)
python3 calculate_tm_matrix.py --pdb_dir rna_chains --output tm_matrix.csv

# 5. Structure clustering
python3 agglomerative_clustering.py --input tm_matrix.csv --thresholds 0.6 0.5 0.4 --plot
```

---

## ⚠️ Important Notes

### Time Estimates

| Number of Structures | TM-score Calculation Time |
|---------------------|---------------------------|
| 50                  | ~1 hour                   |
| 100                 | ~5 hours                  |
| 500                 | ~3-7 days                 |
| 1000                | ~2-4 weeks                |

### Memory Requirements

- CD-HIT: ~2-4 GB RAM
- TM-score calculation: ~4-8 GB RAM
- Clustering: ~2-4 GB RAM

### For Large Datasets (>500 structures)

1. **Use checkpointing:**
   ```bash
   python3 calculate_tm_matrix.py \
       --pdb_dir rna_chains \
       --output tm_matrix.csv \
       --checkpoint checkpoint.csv
   ```

2. **Run in background:**
   ```bash
   nohup python3 calculate_tm_matrix.py \
       --pdb_dir rna_chains \
       --output tm_matrix.csv \
       --checkpoint checkpoint.csv > tm_calc.log 2>&1 &
   
   # Check progress
   tail -f tm_calc.log
   ```

3. **Resume if interrupted:**
   ```bash
   python3 calculate_tm_matrix.py \
       --pdb_dir rna_chains \
       --output tm_matrix.csv \
       --checkpoint checkpoint.csv \
       --resume
   ```

---

## 🐛 Common Issues

### "Command not found: cd-hit-est"
```bash
export PATH=$PATH:~/bin
source ~/.bashrc
```

### "Command not found: USalign"
```bash
export PATH=$PATH:~/bin
source ~/.bashrc
```

### "No module named 'Bio'"
```bash
pip install biopython
```

### TM-score calculation stuck
```bash
# Check if it's actually running
ps aux | grep USalign

# Check progress (if using tqdm)
tail -f tm_calc.log
```

---

## 📚 Files in This Package

| File | Purpose |
|------|---------|
| `README.md` | Comprehensive documentation |
| `rna_clustering_workflow.md` | Detailed technical workflow |
| `install_dependencies.sh` | Install all dependencies |
| `run_clustering.sh` | Run complete pipeline |
| `extract_rna_chains.py` | Extract RNA from PDB files |
| `extract_sequences.py` | Convert to FASTA format |
| `parse_cdhit_clusters.py` | Parse CD-HIT results |
| `calculate_tm_matrix.py` | Calculate TM-scores |
| `agglomerative_clustering.py` | Structure clustering |

---

## 💡 Tips

1. **Test on small subset first:**
   ```bash
   mkdir test_subset
   cp pdb_files/*.pdb test_subset/ | head -20
   ./run_clustering.sh test_subset test_results
   ```

2. **Validate TM-score matrix:**
   ```bash
   python3 calculate_tm_matrix.py --validate --output tm_matrix.csv
   ```

3. **Check cluster quality:**
   - Look at `structure_clusters_comparison.csv`
   - Examine dendrogram for clear clusters
   - Check if singleton clusters are expected

---

## 📖 More Information

- Full documentation: `README.md`
- Technical details: `rna_clustering_workflow.md`
- Example usage in README.md

---

## 🎯 Example: Test Run

```bash
# Quick test with small dataset
mkdir test && cd test

# Get a few PDB files (example)
wget https://files.rcsb.org/download/1EHZ.pdb -P pdb_files/
wget https://files.rcsb.org/download/1F7U.pdb -P pdb_files/
wget https://files.rcsb.org/download/1FFZ.pdb -P pdb_files/

# Run pipeline
../run_clustering.sh pdb_files results

# Check results
cat results/clustering_summary_report.txt
```

---

**For questions and issues, check README.md or rna_clustering_workflow.md**
