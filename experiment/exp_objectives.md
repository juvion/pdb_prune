# Experiment Folders: Objectives and Quick Reference

This document summarizes the purpose and typical workflow of each experiment folder (exp0, exp1, exp2) in this repository. Use it as a quick reference to know what each experiment is about and how to run or extend it.

## exp0 — Data acquisition and preparation

Objective:
- Download, curate, and prepare RNA-related structural data (PDB) and corresponding sequence/array formats (FASTA, NPY).
- Perform basic quality checks and exploratory analysis to ensure downstream readiness.

Inputs:
- PDB ID lists and/or directories with raw PDB files.

Core tasks and typical tools:
- Download and resume downloads of PDBs (utils/pdb_downloader.py, utils/resume_download.py).
- Extract RNA chains and sequences to FASTA (utils/extract_rna_segments.py, utils/extract_sequence.py, utils/pdb_to_fasta.py).
- Convert to NPY arrays for coordinates and related features (utils/pdb_to_npy.py).
- Perform EDA and quality checks (EDA/*, utils/analyze_sequence_lengths.py, utils/analyze_pdb_quality.py).

Outputs:
- Curated PDB directory, rna_sequences.fasta, npy arrays, summary CSVs/plots.
- Notebook(s) for ad-hoc analysis (experiment/exp0/exp0.ipynb, EDA/*).

## exp1 — Competition training data reproduction (for CYX)

Objective:
- Reproduce the competition’s training dataset so CYX can run q9_1.py (and related q9_1_dev*.py variants).
- Align data format and preprocessing with the original competition code to ensure comparable results.

Inputs:
- Competition training data sources (IDs, sequences, metadata) and any accompanying labels.

Core tasks and typical tools:
- Parse and standardize training inputs (scripts in RNA_design_public/, exp1/train/).
- Execute/validate the q9_1.py workflow (see RNA_design_public/ for q9_1_dev* scripts).
- Log outcomes for reproducibility (experiment/exp1/log/).

Outputs:
- Prepared training datasets in expected formats (FASTA/CSV/NPY as needed).
- Logs and intermediate artifacts under experiment/exp1/.

## exp2 — PDB clustering (RIDiffusion-inspired)

Objective:
- Test RNA PDB clustering based on sequence and structure similarity following the approach described in the RIDiffusion paper.
- Produce clusters and summary reports that can be inspected and reused in later stages.

Entry points:
- Shell pipeline: experiment/exp2/pdb_cluster/utils/run_clustering_py2.sh (Python) and run_clustering.sh (Python 3).

Workflow (run_clustering_* scripts):
1) Extract RNA chains from PDBs (utils/extract_rna_chains.py) into rna_chains/.
2) Extract sequences to FASTA (utils/extract_sequences.py) into rna_sequences.fasta.
3) Sequence clustering with CD-HIT-EST (cd-hit-est):
   - Note: CD-HIT-EST requires identity threshold -c ≥ 0.8.
   - The script skips thresholds below 0.8 and logs a warning.
   - Outputs: clusters_seq_<thr>.fasta and sequence_clusters_<thr>.csv (+ distribution).
4) Structure similarity matrix using USalign (utils/calculate_tm_matrix.py) → tm_score_matrix.csv.
5) Agglomerative clustering on structure matrix (utils/agglomerative_clustering.py): outputs per-threshold CSVs, plots, and distributions.
6) Summary report aggregation (clustering_summary_report.txt) listing inputs, thresholds, and output files.

Outputs:
- rna_sequences.fasta, tm_score_matrix.csv, sequence_clusters_*.csv, structure_clusters_*.csv, plots, and a summary report under the chosen output directory.

Notes:
- The pipeline scripts have been updated to use an absolute PROJECT_HOME for utility paths and to skip invalid CD-HIT thresholds (< 0.8).
- Ensure required tools are available: python, cd-hit-est, USalign; and packages: biopython, numpy, pandas, scipy, matplotlib, tqdm.


## exp3 — PDB data preparation for RIDiffusion (Aborted)
PDB list (all_pdbid.txt) was provided by Xuke (From RIDifussion) 
The dataset now does not need to re-download from PDB, since the original data can be found in the paper's supplementary material.


---

If you add new experiments (e.g., exp3), please append a section here with objectives, inputs, core tasks, and outputs to keep this index up-to-date.


