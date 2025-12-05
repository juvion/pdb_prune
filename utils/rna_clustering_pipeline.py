#!/usr/bin/env python3
"""
RNA Clustering and Data Splitting Pipeline

A simple, effective script for clustering RNA data based on sequence and structure similarity,
followed by cluster-based splitting into train/validation/test sets.

Usage:
    python rna_clustering_pipeline.py --pdb_dir <path> --fasta_dir <path> --output_dir <path>
"""

import os
import sys
import glob
import subprocess
import argparse
import logging
from pathlib import Path
from collections import defaultdict
import networkx as nx

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class RNAClusteringPipeline:
    def __init__(self, pdb_dir, fasta_dir, output_dir, max_seq_len=1000):
        self.pdb_dir = Path(pdb_dir)
        self.fasta_dir = Path(fasta_dir)
        self.output_dir = Path(output_dir)
        self.max_seq_len = max_seq_len
        self.output_dir.mkdir(exist_ok=True)
        
        # Clustering thresholds
        self.seq_thresholds = [0.8, 0.9, 1.0]
        self.struct_thresholds = [0.4, 0.45, 0.5, 0.6]
        
        # Split ratios
        self.test_ratio = 0.15
        self.val_ratio = 0.10
        self.train_ratio = 0.75
        
    def collect_samples(self):
        """Collect and filter RNA samples by length."""
        logger.info("Collecting RNA samples...")
        
        pdb_files = list(self.pdb_dir.glob("*.pdb"))
        valid_samples = []
        
        for pdb_file in pdb_files:
            pdb_id = pdb_file.stem
            fasta_file = self.fasta_dir / f"{pdb_id}.fasta"
            
            if not fasta_file.exists():
                logger.warning(f"Missing FASTA file for {pdb_id}")
                continue
                
            # Check sequence length
            try:
                with open(fasta_file, 'r') as f:
                    lines = f.readlines()
                    sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
                    
                # Filter by both minimum and maximum length
                # US-align requires at least 3 residues to work properly
                if 3 <= len(sequence) <= self.max_seq_len:
                    valid_samples.append({
                        'pdb_id': pdb_id,
                        'pdb_file': pdb_file,
                        'fasta_file': fasta_file,
                        'length': len(sequence)
                    })
                elif len(sequence) < 3:
                    logger.debug(f"Skipping {pdb_id}: sequence too short ({len(sequence)} residues, minimum 3 required)")
                    
            except Exception as e:
                logger.error(f"Error reading {fasta_file}: {e}")
                
        logger.info(f"Found {len(valid_samples)} valid samples (max length: {self.max_seq_len})")
        return valid_samples
    
    def create_multifasta(self, samples):
        """Create multi-FASTA file for PSI-CD-HIT."""
        multifasta_file = self.output_dir / "all_sequences.fasta"
        
        with open(multifasta_file, 'w') as out_f:
            for sample in samples:
                with open(sample['fasta_file'], 'r') as in_f:
                    lines = in_f.readlines()
                    # Write header with PDB ID
                    out_f.write(f">{sample['pdb_id']}\n")
                    # Write sequence
                    for line in lines:
                        if not line.startswith('>'):
                            out_f.write(line)
                            
        return multifasta_file
    
    def run_sequence_clustering(self, multifasta_file):
        """Run PSI-CD-HIT for sequence clustering."""
        logger.info("Running sequence clustering with PSI-CD-HIT...")
        
        seq_clusters = {}
        
        for threshold in self.seq_thresholds:
            logger.info(f"Clustering at sequence similarity {threshold}")
            
            output_prefix = self.output_dir / f"seq_clusters_{threshold}"
            cluster_file = f"{output_prefix}.clstr"
            
            # Run CD-HIT (PSI-CD-HIT functionality)
            cmd = [
                "cd-hit-est",  # CD-HIT for nucleotide sequences
                "-i", str(multifasta_file),
                "-o", str(output_prefix),
                "-c", str(threshold),
                "-n", "8",  # word length for nucleotides
                "-M", "2000",  # memory limit in MB
                "-T", "4"  # number of threads
            ]
            
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                clusters = self.parse_psicd_hit_clusters(cluster_file)
                seq_clusters[threshold] = clusters
                logger.info(f"Found {len(clusters)} clusters at threshold {threshold}")
                
            except subprocess.CalledProcessError as e:
                logger.error(f"PSI-CD-HIT failed for threshold {threshold}: {e}")
                seq_clusters[threshold] = []
                
        return seq_clusters
    
    def parse_psicd_hit_clusters(self, cluster_file):
        """Parse PSI-CD-HIT cluster output."""
        clusters = []
        current_cluster = []
        
        try:
            with open(cluster_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>Cluster'):
                        if current_cluster:
                            clusters.append(current_cluster)
                        current_cluster = []
                    elif line and not line.startswith('>'):
                        # Extract PDB ID from cluster member line
                        parts = line.split('>')
                        if len(parts) > 1:
                            pdb_id = parts[1].split('...')[0]
                            current_cluster.append(pdb_id)
                            
                if current_cluster:
                    clusters.append(current_cluster)
                    
        except FileNotFoundError:
            logger.error(f"Cluster file not found: {cluster_file}")
            
        return clusters
    
    def calculate_tm_scores(self, samples):
        """Calculate TM-scores between all pairs using US-align."""
        logger.info("Calculating TM-scores with US-align...")
        
        tm_scores = {}
        pdb_ids = [s['pdb_id'] for s in samples]
        
        for i, sample1 in enumerate(samples):
            for j, sample2 in enumerate(samples[i+1:], i+1):
                pdb1, pdb2 = sample1['pdb_id'], sample2['pdb_id']
                
                try:
                    cmd = ["US-align", str(sample1['pdb_file']), str(sample2['pdb_file'])]
                    result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=30)
                    
                    # Parse TM-score from output
                    tm_score = self.parse_tm_score(result.stdout)
                    tm_scores[(pdb1, pdb2)] = tm_score
                    tm_scores[(pdb2, pdb1)] = tm_score  # Symmetric
                    
                except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
                    # Handle both regular errors and crashes (SIGSEGV, etc.)
                    if hasattr(e, 'returncode') and e.returncode < 0:
                        # Negative return codes indicate signals (e.g., SIGSEGV = -11)
                        logger.warning(f"US-align crashed for {pdb1}-{pdb2}: Signal {abs(e.returncode)} (likely due to problematic PDB structure)")
                    else:
                        logger.warning(f"US-align failed for {pdb1}-{pdb2}: {e}")
                    tm_scores[(pdb1, pdb2)] = 0.0
                    tm_scores[(pdb2, pdb1)] = 0.0
                except Exception as e:
                    # Catch any other unexpected errors
                    logger.warning(f"Unexpected error for {pdb1}-{pdb2}: {e}")
                    tm_scores[(pdb1, pdb2)] = 0.0
                    tm_scores[(pdb2, pdb1)] = 0.0
                    
        logger.info(f"Calculated {len(tm_scores)//2} pairwise TM-scores")
        return tm_scores
    
    def parse_tm_score(self, output):
        """Parse TM-score from US-align output."""
        for line in output.split('\n'):
            if 'TM-score=' in line:
                try:
                    tm_score = float(line.split('TM-score=')[1].split()[0])
                    return tm_score
                except (IndexError, ValueError):
                    continue
        return 0.0
    
    def run_structure_clustering(self, samples, tm_scores):
        """Run graph-based structure clustering."""
        logger.info("Running structure clustering...")
        
        struct_clusters = {}
        pdb_ids = [s['pdb_id'] for s in samples]
        
        for threshold in self.struct_thresholds:
            logger.info(f"Clustering at TM-score threshold {threshold}")
            
            # Create graph
            G = nx.Graph()
            G.add_nodes_from(pdb_ids)
            
            # Add edges based on TM-score threshold
            for (pdb1, pdb2), tm_score in tm_scores.items():
                if pdb1 < pdb2 and tm_score >= threshold:  # Avoid duplicate edges
                    G.add_edge(pdb1, pdb2)
                    
            # Find connected components (clusters)
            clusters = [list(component) for component in nx.connected_components(G)]
            struct_clusters[threshold] = clusters
            
            logger.info(f"Found {len(clusters)} clusters at threshold {threshold}")
            
        return struct_clusters
    
    def split_clusters(self, clusters, samples):
        """Split clusters into train/val/test sets."""
        # Calculate average sequence length
        avg_length = sum(s['length'] for s in samples) / len(samples)
        
        # Sort clusters by size (largest first)
        sorted_clusters = sorted(clusters, key=len, reverse=True)
        
        train_clusters, val_clusters, test_clusters = [], [], []
        train_size, val_size, test_size = 0, 0, 0
        total_samples = sum(len(cluster) for cluster in clusters)
        
        # Assign large clusters (>30 samples) to training
        for cluster in sorted_clusters[:]:
            if len(cluster) > 30:
                train_clusters.append(cluster)
                train_size += len(cluster)
                sorted_clusters.remove(cluster)
                
        # Distribute remaining clusters
        for cluster in sorted_clusters:
            cluster_size = len(cluster)
            
            # Calculate current ratios
            current_test_ratio = test_size / total_samples if total_samples > 0 else 0
            current_val_ratio = val_size / total_samples if total_samples > 0 else 0
            
            # Assign to the set that needs more samples
            if current_test_ratio < self.test_ratio:
                test_clusters.append(cluster)
                test_size += cluster_size
            elif current_val_ratio < self.val_ratio:
                val_clusters.append(cluster)
                val_size += cluster_size
            else:
                train_clusters.append(cluster)
                train_size += cluster_size
                
        return {
            'train': train_clusters,
            'val': val_clusters,
            'test': test_clusters
        }
    
    def save_splits(self, splits, method, threshold):
        """Save train/val/test splits to files."""
        for split_name, clusters in splits.items():
            filename = self.output_dir / f"{method}_{threshold}_{split_name}.txt"
            
            with open(filename, 'w') as f:
                for cluster in clusters:
                    for pdb_id in cluster:
                        f.write(f"{pdb_id}\n")
                        
            # Log statistics
            total_samples = sum(len(cluster) for cluster in clusters)
            logger.info(f"{method}_{threshold}_{split_name}: {len(clusters)} clusters, {total_samples} samples")
    
    def run_pipeline(self):
        """Run the complete clustering and splitting pipeline."""
        logger.info("Starting RNA clustering pipeline...")
        
        # Step 1: Collect and filter samples
        samples = self.collect_samples()
        if not samples:
            logger.error("No valid samples found!")
            return
            
        # Step 2: Sequence clustering
        multifasta_file = self.create_multifasta(samples)
        seq_clusters = self.run_sequence_clustering(multifasta_file)
        
        # Process sequence clustering results
        for threshold, clusters in seq_clusters.items():
            if clusters:
                splits = self.split_clusters(clusters, samples)
                self.save_splits(splits, "seq", threshold)
                
        # Step 3: Structure clustering
        tm_scores = self.calculate_tm_scores(samples)
        struct_clusters = self.run_structure_clustering(samples, tm_scores)
        
        # Process structure clustering results
        for threshold, clusters in struct_clusters.items():
            if clusters:
                splits = self.split_clusters(clusters, samples)
                self.save_splits(splits, "struct", threshold)
                
        logger.info("Pipeline completed successfully!")

def main():
    parser = argparse.ArgumentParser(description="RNA Clustering and Data Splitting Pipeline")
    parser.add_argument("--pdb_dir", required=True, help="Directory containing PDB files")
    parser.add_argument("--fasta_dir", required=True, help="Directory containing FASTA files")
    parser.add_argument("--output_dir", required=True, help="Output directory for results")
    parser.add_argument("--max_seq_len", type=int, default=500, help="Maximum sequence length")
    
    args = parser.parse_args()
    
    # Validate directories
    for dir_path in [args.pdb_dir, args.fasta_dir]:
        if not os.path.exists(dir_path):
            logger.error(f"Directory does not exist: {dir_path}")
            sys.exit(1)
            
    # Run pipeline
    pipeline = RNAClusteringPipeline(
        pdb_dir=args.pdb_dir,
        fasta_dir=args.fasta_dir,
        output_dir=args.output_dir,
        max_seq_len=args.max_seq_len
    )
    
    pipeline.run_pipeline()

if __name__ == "__main__":
    main()