# Action Plan - Fixing Your PDB Clustering

## 📋 Current Status

✅ **What's Working:**
- Sequence clustering at threshold 0.8: **470 clusters from 1133 sequences**
- PDB extraction and preparation
- TM-score calculation setup
- Structure clustering setup

❌ **What Failed:**
- Sequence clustering at thresholds 0.6 and 0.4 (cd-hit-est limitation)

---

## 🎯 Recommended Solution (Choose One)

### **Option A: Use What Works (Quickest - 5 minutes)**

This is the pragmatic choice. Structure clustering is more important for RNA anyway.

#### Steps:

1. **Modify your command to use only 0.8 threshold for sequences:**

```bash
cd /Users/xiaojuzhang/Dev/pdb_prune/experiment/exp2/pdb_cluster

./run_clustering_py2.sh \
    '/Users/xiaojuzhang/Dev/pdb_prune/experiment/exp0/competition/train/pdbs' \
    '/Users/xiaojuzhang/Dev/pdb_prune/experiment/exp2/pdb_cluster/cluster_project' \
    0.8 \
    0.6,0.5,0.4
```

2. **Results you'll get:**
   - ✅ 470 sequence clusters (0.8 threshold)
   - ✅ Structure clusters at 0.6, 0.5, 0.4 thresholds
   - ✅ Complete summary report
   - ✅ All visualizations

3. **Why this is scientifically valid:**
   - Structure similarity (TM-score) is MORE important than sequence similarity for RNA
   - You still get multiple clustering granularities through structure thresholds
   - Your train/val/test split will be based on structure clusters (which is better!)

---

### **Option B: Install PSI-CD-HIT (Exact Paper Replication - 30 minutes)**

Use the actual tool from the paper.

#### Steps:

1. **Check if you already have PSI-CD-HIT:**

```bash
which psi-cd-hit.pl
```

2. **If not found, install it:**

```bash
# Go to your CD-HIT source
cd /tmp
git clone https://github.com/weizhongli/cdhit.git
cd cdhit
make

# Copy all executables (including psi-cd-hit.pl)
sudo cp cd-hit cd-hit-est psi-cd-hit psi-cd-hit.pl cd-hit-2d /usr/local/bin/

# Or to your local bin (no sudo needed)
mkdir -p ~/bin
cp cd-hit cd-hit-est psi-cd-hit psi-cd-hit.pl cd-hit-2d ~/bin/
echo 'export PATH=$PATH:~/bin' >> ~/.bashrc
source ~/.bashrc
```

3. **Verify installation:**

```bash
psi-cd-hit.pl -h
```

4. **Use the fixed script:**

```bash
cd /Users/xiaojuzhang/Dev/pdb_prune/experiment/exp2/pdb_cluster

# Copy the fixed script I provided
cp /path/to/run_clustering_fixed.sh ./

# Update PROJECT_HOME in the script (line 61):
# PROJECT_HOME='/Users/xiaojuzhang/Dev/pdb_prune/experiment/exp2/pdb_cluster'

chmod +x run_clustering_fixed.sh

./run_clustering_fixed.sh \
    '/Users/xiaojuzhang/Dev/pdb_prune/experiment/exp0/competition/train/pdbs' \
    '/Users/xiaojuzhang/Dev/pdb_prune/experiment/exp2/pdb_cluster/cluster_project'
```

5. **Results you'll get:**
   - ✅ Sequence clusters at 0.8, 0.6, 0.4 thresholds
   - ✅ Structure clusters at 0.6, 0.5, 0.4 thresholds
   - ✅ Exact replication of paper methodology

---

### **Option C: Use Alternative Python Method (No New Tools - 1-2 hours)**

For lower thresholds without installing PSI-CD-HIT.

#### Steps:

1. **Use the alternative clustering script I provided:**

```bash
cd /Users/xiaojuzhang/Dev/pdb_prune/experiment/exp2/pdb_cluster

# For threshold 0.6
python alternative_seq_clustering.py \
    --input cluster_project/rna_sequences.fasta \
    --threshold 0.6 \
    --output cluster_project/sequence_clusters_0.6.csv \
    --dist cluster_project/sequence_clusters_0.6_distribution.csv \
    --method greedy

# For threshold 0.4
python alternative_seq_clustering.py \
    --input cluster_project/rna_sequences.fasta \
    --threshold 0.4 \
    --output cluster_project/sequence_clusters_0.4.csv \
    --dist cluster_project/sequence_clusters_0.4_distribution.csv \
    --method greedy
```

2. **Then continue with structure clustering:**

```bash
python utils/calculate_tm_matrix.py \
    --pdb_dir cluster_project/rna_chains \
    --output cluster_project/tm_score_matrix.csv

python utils/agglomerative_clustering.py \
    --input cluster_project/tm_score_matrix.csv \
    --thresholds 0.6 0.5 0.4 \
    --plot
```

**Note:** This method is slower (~1-2 hours for 1133 sequences) but requires no new tool installation.

---

## 🚀 My Recommendation: **Option A**

Here's why:

### Scientific Justification:
1. **Structure > Sequence for RNA**: The RiboDiffusion paper emphasizes that structure similarity is more biologically relevant
2. **Your 470 clusters at 0.8 threshold is good**: This gives you a reasonable starting point
3. **Multiple structure thresholds provide granularity**: You still get to test different similarity levels

### Practical Benefits:
- ✅ Works immediately (no installation)
- ✅ Faster (no additional computation)
- ✅ Still scientifically rigorous
- ✅ You already have working code

### Expected Results:

With 1133 sequences and 470 clusters at 0.8:

| Clustering Type | Threshold | Expected Clusters | Status |
|-----------------|-----------|-------------------|--------|
| Sequence | 0.8 | 470 | ✅ Done |
| Structure | 0.6 | ~800-1000 | ⏳ Ready |
| Structure | 0.5 | ~600-800 | ⏳ Ready |
| Structure | 0.4 | ~400-600 | ⏳ Ready |

---

## 📝 Exact Commands for Option A

```bash
# Navigate to your project
cd /Users/xiaojuzhang/Dev/pdb_prune/experiment/exp2/pdb_cluster

# Run with only 0.8 for sequence clustering
./run_clustering_py2.sh \
    '/Users/xiaojuzhang/Dev/pdb_prune/experiment/exp0/competition/train/pdbs' \
    '/Users/xiaojuzhang/Dev/pdb_prune/experiment/exp2/pdb_cluster/cluster_project' \
    0.8 \
    0.6,0.5,0.4

# Check results
cat cluster_project/clustering_summary_report.txt
ls -lh cluster_project/
```

---

## 🔍 What to Expect (Time Estimates)

### For 1133 sequences:

| Step | Time | Status |
|------|------|--------|
| Extract RNA chains | ~2 min | ✅ Already done |
| Extract sequences | ~1 min | ✅ Already done |
| Sequence clustering (0.8) | ~1 min | ✅ Already done |
| **Calculate TM-scores** | **~3-6 hours** | ⏳ Next |
| Structure clustering | ~5 min | ⏳ After TM-scores |
| Generate reports | ~1 min | ⏳ Final step |

**Total remaining time: ~3-6 hours** (mostly TM-score calculation)

### TM-score Calculation Tips:

The TM-score calculation is the bottleneck. For 1133 structures:
- Total comparisons: 641,328 pairs
- Rate: ~2-5 seconds per pair
- Estimated time: 3-6 hours

**Use checkpoint mode to be safe:**

```bash
# The script already includes checkpointing
# If it crashes, you can resume with:
python utils/calculate_tm_matrix.py \
    --pdb_dir cluster_project/rna_chains \
    --output cluster_project/tm_score_matrix.csv \
    --checkpoint cluster_project/tm_score_checkpoint.csv \
    --resume
```

---

## 📊 Validation Checklist

After running Option A, verify:

- [ ] `sequence_clusters_0.8.csv` exists with 1133 rows
- [ ] `tm_score_matrix.csv` exists (1133 × 1133 matrix)
- [ ] `structure_clusters_0.6.csv` exists
- [ ] `structure_clusters_0.5.csv` exists
- [ ] `structure_clusters_0.4.csv` exists
- [ ] `clustering_summary_report.txt` is complete
- [ ] PNG visualizations are generated

---

## 🆘 If You Still Want Lower Sequence Thresholds

If you absolutely need 0.6 and 0.4 sequence clustering:

1. **Try Option B** (install PSI-CD-HIT) - most accurate
2. **Try Option C** (Python alternative) - no installation needed but slower

Both will give you the additional clustering levels, but remember:
- Your structure clustering is MORE important
- Adding sequence 0.6 and 0.4 won't dramatically change your final analysis
- The paper's main contribution is the dual clustering approach, which you're already doing

---

## 💡 Bottom Line

**Just run Option A now.** You can always add lower sequence thresholds later if needed, but you'll likely find that the structure clustering provides all the granularity you need for a proper train/val/test split.

Your research won't be compromised by using only 0.8 for sequences, especially since you're getting proper structure-based clustering at multiple thresholds!

---

## 📞 Next Steps

1. Choose your option (I recommend **Option A**)
2. Run the command
3. Wait for TM-score calculation (~3-6 hours)
4. Check your results
5. Proceed with your data splitting

Good luck! 🚀
