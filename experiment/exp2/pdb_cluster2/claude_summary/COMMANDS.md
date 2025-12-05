# 🚀 QUICK COMMAND REFERENCE

## Copy-Paste Commands for Your System

Based on your file paths from the logs:

### ✅ Option A: Run with 0.8 threshold only (RECOMMENDED)

```bash
cd /Users/xiaojuzhang/Dev/pdb_prune/experiment/exp2/pdb_cluster

./run_clustering_py2.sh \
    '/Users/xiaojuzhang/Dev/pdb_prune/experiment/exp0/competition/train/pdbs' \
    '/Users/xiaojuzhang/Dev/pdb_prune/experiment/exp2/pdb_cluster/cluster_project' \
    0.8 \
    0.6,0.5,0.4
```

**Time:** ~3-6 hours  
**Result:** Sequence clustering (0.8) + Structure clustering (0.6, 0.5, 0.4)

---

### 🔧 Option B: Use fixed script (if you install PSI-CD-HIT)

1. **First, copy the fixed script to your project:**

```bash
cd /Users/xiaojuzhang/Dev/pdb_prune/experiment/exp2/pdb_cluster

# Copy from wherever you downloaded it
cp /path/to/downloaded/run_clustering_fixed.sh ./

# Make executable
chmod +x run_clustering_fixed.sh
```

2. **Edit the PROJECT_HOME variable (line 61):**

```bash
# Open in editor
nano run_clustering_fixed.sh

# Change line 61 to:
PROJECT_HOME='/Users/xiaojuzhang/Dev/pdb_prune/experiment/exp2/pdb_cluster'

# Save and exit
```

3. **Run it:**

```bash
./run_clustering_fixed.sh \
    '/Users/xiaojuzhang/Dev/pdb_prune/experiment/exp0/competition/train/pdbs' \
    '/Users/xiaojuzhang/Dev/pdb_prune/experiment/exp2/pdb_cluster/cluster_project'
```

---

### 🐍 Option C: Python alternative for lower thresholds

1. **First complete what works (0.8 + structure clustering):**

```bash
cd /Users/xiaojuzhang/Dev/pdb_prune/experiment/exp2/pdb_cluster

./run_clustering_py2.sh \
    '/Users/xiaojuzhang/Dev/pdb_prune/experiment/exp0/competition/train/pdbs' \
    '/Users/xiaojuzhang/Dev/pdb_prune/experiment/exp2/pdb_cluster/cluster_project' \
    0.8 \
    0.6,0.5,0.4
```

2. **Then add lower thresholds with Python:**

```bash
cd /Users/xiaojuzhang/Dev/pdb_prune/experiment/exp2/pdb_cluster

# For threshold 0.6
python alternative_seq_clustering.py \
    --input cluster_project/rna_sequences.fasta \
    --threshold 0.6 \
    --output cluster_project/sequence_clusters_0.6.csv \
    --dist cluster_project/sequence_clusters_0.6_distribution.csv

# For threshold 0.4
python alternative_seq_clustering.py \
    --input cluster_project/rna_sequences.fasta \
    --threshold 0.4 \
    --output cluster_project/sequence_clusters_0.4.csv \
    --dist cluster_project/sequence_clusters_0.4_distribution.csv
```

**Note:** Each Python clustering takes ~1-2 hours for 1133 sequences.

---

## 📊 Check Your Results

After completion:

```bash
cd /Users/xiaojuzhang/Dev/pdb_prune/experiment/exp2/pdb_cluster/cluster_project

# View summary
cat clustering_summary_report.txt

# List all output files
ls -lh *.csv *.png

# Check sequence clustering
wc -l sequence_clusters_*.csv

# Check structure clustering
wc -l structure_clusters_*.csv
```

---

## 🔍 Monitor Progress (If Running in Background)

```bash
# Run in background
nohup ./run_clustering_py2.sh \
    '/Users/xiaojuzhang/Dev/pdb_prune/experiment/exp0/competition/train/pdbs' \
    '/Users/xiaojuzhang/Dev/pdb_prune/experiment/exp2/pdb_cluster/cluster_project' \
    0.8 \
    0.6,0.5,0.4 \
    > clustering.log 2>&1 &

# Get process ID
echo $!

# Monitor progress
tail -f clustering.log

# Check if still running
ps aux | grep clustering

# Check USalign (TM-score calculation)
ps aux | grep USalign
```

---

## 🛑 Resume if Interrupted

If TM-score calculation stops:

```bash
cd /Users/xiaojuzhang/Dev/pdb_prune/experiment/exp2/pdb_cluster

python utils/calculate_tm_matrix.py \
    --pdb_dir cluster_project/rna_chains \
    --output cluster_project/tm_score_matrix.csv \
    --checkpoint cluster_project/tm_score_checkpoint.csv \
    --resume
```

Then complete the pipeline:

```bash
python utils/agglomerative_clustering.py \
    --input cluster_project/tm_score_matrix.csv \
    --thresholds 0.6 0.5 0.4 \
    --method average \
    --output_prefix cluster_project/structure_clusters \
    --plot
```

---

## 📝 Expected Output Files

After successful completion, you should have:

```
cluster_project/
├── rna_chains/                           # 1133 PDB files
│   └── rna_chains_metadata.json
├── rna_sequences.fasta                   # 1133 sequences
├── sequence_clusters_0.8.csv             # 470 clusters
├── sequence_clusters_0.8_distribution.csv
├── tm_score_matrix.csv                   # 1133×1133 matrix (~10-20MB)
├── structure_clusters_0.6.csv
├── structure_clusters_0.5.csv
├── structure_clusters_0.4.csv
├── structure_clusters_0.6_distribution.csv
├── structure_clusters_0.5_distribution.csv
├── structure_clusters_0.4_distribution.csv
├── structure_clusters_dendrogram.png
├── structure_clusters_size_distribution.png
├── structure_clusters_comparison.csv
└── clustering_summary_report.txt
```

---

## ⚡ One-Liner (Just Copy and Run)

For the impatient:

```bash
cd /Users/xiaojuzhang/Dev/pdb_prune/experiment/exp2/pdb_cluster && ./run_clustering_py2.sh '/Users/xiaojuzhang/Dev/pdb_prune/experiment/exp0/competition/train/pdbs' '/Users/xiaojuzhang/Dev/pdb_prune/experiment/exp2/pdb_cluster/cluster_project' 0.8 0.6,0.5,0.4
```

---

## 📖 Which Files Should I Read?

1. **[START_HERE.md](computer:///mnt/user-data/outputs/START_HERE.md)** - Overview and context ⭐
2. **This file** - Quick commands to run
3. **[ACTION_PLAN.md](computer:///mnt/user-data/outputs/ACTION_PLAN.md)** - Detailed comparison of options
4. **[TROUBLESHOOTING_CDHIT.md](computer:///mnt/user-data/outputs/TROUBLESHOOTING_CDHIT.md)** - If you encounter issues

---

## 🎯 Bottom Line

**Just run Option A now:**

```bash
cd /Users/xiaojuzhang/Dev/pdb_prune/experiment/exp2/pdb_cluster
./run_clustering_py2.sh \
    '/Users/xiaojuzhang/Dev/pdb_prune/experiment/exp0/competition/train/pdbs' \
    '/Users/xiaojuzhang/Dev/pdb_prune/experiment/exp2/pdb_cluster/cluster_project' \
    0.8 \
    0.6,0.5,0.4
```

**Then get coffee. It'll take 3-6 hours. ☕**

When it's done, check the summary report and proceed with your research! 🚀
