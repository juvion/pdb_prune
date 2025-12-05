#!/usr/bin/env python3
"""
Parallel RNA Clustering and Data Splitting Pipeline

A highly optimized version that uses multiprocessing for parallel TM-score calculations,
making structure clustering feasible for large datasets.

Features:
- Parallel TM-score calculation using multiprocessing
- Intelligent work distribution and load balancing
- Progress tracking with shared state
- Memory-efficient batch processing
- Resumption capabilities with caching

Usage:
    python rna_clustering_pipeline_parallel.py --pdb_dir <path> --fasta_dir <path> --output_dir <path>
"""

import os
import sys
import glob
import subprocess
import argparse
import logging
import pickle
import random
import multiprocessing as mp
from pathlib import Path
from collections import defaultdict
from functools import partial
import networkx as nx
from tqdm import tqdm
import time
import itertools

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def calculate_tm_score_pair(pair_data):
    """
    Calculate TM-score for a single pair of PDB files.
    This function is designed to be used with multiprocessing.
    
    Args:
        pair_data: tuple containing (pdb1_path, pdb2_path, pdb1_id, pdb2_id)
    
    Returns:
        tuple: (pdb1_id, pdb2_id, tm_score)
    """
    pdb1_path, pdb2_path, pdb1_id, pdb2_id = pair_data
    
    try:
        cmd = ["US-align", str(pdb1_path), str(pdb2_path)]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=30)
        
        # Parse TM-score from output
        tm_score = parse_tm_score_from_output(result.stdout)
        return (pdb1_id, pdb2_id, tm_score)
        
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
        # Handle both regular errors and crashes (SIGSEGV, etc.)
        return (pdb1_id, pdb2_id, 0.0)
    except Exception as e:
        # Catch any other unexpected errors
        return (pdb1_id, pdb2_id, 0.0)

def parse_tm_score_from_output(output):
    """Parse TM-score from US-align output."""
    for line in output.split('\n'):
        if 'TM-score=' in line:
            try:
                tm_score = float(line.split('TM-score=')[1].split()[0])
                return tm_score
            except (IndexError, ValueError):
                continue
    return 0.0

class ParallelRNAClusteringPipeline:
    def __init__(self, pdb_dir, fasta_dir, output_dir, max_seq_len=1000, 
                 max_struct_samples=2000, skip_structure=False, n_processes=None):
        self.pdb_dir = Path(pdb_dir)
        self.fasta_dir = Path(fasta_dir)
        self.output_dir = Path(output_dir)
        self.max_seq_len = max_seq_len
        self.max_struct_samples = max_struct_samples
        self.skip_structure = skip_structure
        self.n_processes = n_processes or max(1, mp.cpu_count() - 1)  # Leave one core free
        self.output_dir.mkdir(exist_ok=True)
        
        # Clustering thresholds
        self.seq_thresholds = [0.8, 0.9, 1.0]
        self.struct_thresholds = [0.4, 0.45, 0.5, 0.6]
        
        # Split ratios
        self.test_ratio = 0.15
        self.val_ratio = 0.10
        self.train_ratio = 0.75
        
        # Cache files for resumption
        self.tm_scores_cache = self.output_dir / "tm_scores_cache_parallel.pkl"
        
        logger.info(f"Initialized parallel pipeline with {self.n_processes} processes")
        
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
                if 3 <= len(sequence) <= self.max_seq_len:
                    valid_samples.append({
                        'pdb_id': pdb_id,
                        'pdb_file': pdb_file,
                        'fasta_file': fasta_file,
                        'length': len(sequence)
                    })
                elif len(sequence) < 3:
                    logger.debug(f"Skipping {pdb_id}: sequence too short ({len(sequence)} residues)")
                    
            except Exception as e:
                logger.error(f"Error reading {fasta_file}: {e}")
                
        logger.info(f"Found {len(valid_samples)} valid samples (max length: {self.max_seq_len})")
        return valid_samples
    
    def create_multifasta(self, samples):
        """Create a multi-FASTA file for sequence clustering."""
        multifasta_file = self.output_dir / "all_rna_sequences.fasta"
        
        with open(multifasta_file, 'w') as outfile:
            for sample in samples:
                with open(sample['fasta_file'], 'r') as infile:
                    content = infile.read()
                    if not content.startswith('>'):
                        outfile.write(f">{sample['pdb_id']}\n")
                    outfile.write(content)
                    if not content.endswith('\n'):
                        outfile.write('\n')
                        
        return multifasta_file
    
    def run_sequence_clustering(self, multifasta_file):
        """Run CD-HIT for sequence clustering."""
        logger.info("Running sequence clustering with CD-HIT...")
        
        seq_clusters = {}
        
        for threshold in self.seq_thresholds:
            logger.info(f"Clustering at sequence similarity {threshold}")
            
            output_prefix = self.output_dir / f"seq_clusters_{threshold}"
            cluster_file = f"{output_prefix}.clstr"
            
            # Run CD-HIT with multiple threads
            cmd = [
                "cd-hit-est",
                "-i", str(multifasta_file),
                "-o", str(output_prefix),
                "-c", str(threshold),
                "-n", "8",
                "-M", "4000",  # Increased memory limit
                "-T", str(self.n_processes)  # Use all available cores
            ]
            
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                clusters = self.parse_cd_hit_clusters(cluster_file)
                seq_clusters[threshold] = clusters
                logger.info(f"Found {len(clusters)} clusters at threshold {threshold}")
                
            except subprocess.CalledProcessError as e:
                logger.error(f"CD-HIT failed for threshold {threshold}: {e}")
                seq_clusters[threshold] = []
                
        return seq_clusters
    
    def parse_cd_hit_clusters(self, cluster_file):
        """Parse CD-HIT cluster output."""
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
                        parts = line.split('>')
                        if len(parts) > 1:
                            pdb_id = parts[1].split('...')[0]
                            current_cluster.append(pdb_id)
                            
                if current_cluster:
                    clusters.append(current_cluster)
                    
        except FileNotFoundError:
            logger.error(f"Cluster file not found: {cluster_file}")
            
        return clusters
    
    def sample_for_structure_clustering(self, samples):
        """Sample a subset of structures for TM-score calculation."""
        if len(samples) <= self.max_struct_samples:
            return samples
            
        logger.info(f"Sampling {self.max_struct_samples} structures from {len(samples)} for structure clustering")
        
        # Stratified sampling by length to maintain diversity
        samples_by_length = defaultdict(list)
        for sample in samples:
            length_bin = sample['length'] // 50
            samples_by_length[length_bin].append(sample)
        
        sampled = []
        samples_per_bin = max(1, self.max_struct_samples // len(samples_by_length))
        
        for length_bin, bin_samples in samples_by_length.items():
            n_to_sample = min(samples_per_bin, len(bin_samples))
            sampled.extend(random.sample(bin_samples, n_to_sample))
        
        if len(sampled) < self.max_struct_samples:
            remaining = [s for s in samples if s not in sampled]
            additional_needed = self.max_struct_samples - len(sampled)
            if remaining:
                sampled.extend(random.sample(remaining, min(additional_needed, len(remaining))))
        
        logger.info(f"Selected {len(sampled)} samples for structure clustering")
        return sampled
    
    def calculate_tm_scores_parallel(self, samples):
        """Calculate TM-scores using parallel processing."""
        # Check if cached results exist
        if self.tm_scores_cache.exists():
            logger.info("Loading cached TM-scores...")
            with open(self.tm_scores_cache, 'rb') as f:
                return pickle.load(f)
        
        logger.info("Calculating TM-scores with parallel US-align...")
        
        # Sample structures if dataset is too large
        struct_samples = self.sample_for_structure_clustering(samples)
        
        # Prepare all pairs for parallel processing
        pairs_data = []
        for i, sample1 in enumerate(struct_samples):
            for j, sample2 in enumerate(struct_samples[i+1:], i+1):
                pairs_data.append((
                    sample1['pdb_file'], 
                    sample2['pdb_file'],
                    sample1['pdb_id'], 
                    sample2['pdb_id']
                ))
        
        total_pairs = len(pairs_data)
        logger.info(f"Computing {total_pairs} pairwise TM-scores using {self.n_processes} processes...")
        
        # Process pairs in parallel with progress tracking
        tm_scores = {}
        
        # Use multiprocessing Pool for parallel execution
        with mp.Pool(processes=self.n_processes) as pool:
            # Use imap for better memory efficiency and progress tracking
            results = list(tqdm(
                pool.imap(calculate_tm_score_pair, pairs_data, chunksize=max(1, total_pairs // (self.n_processes * 4))),
                total=total_pairs,
                desc="TM-score calculation"
            ))
        
        # Process results
        for pdb1_id, pdb2_id, tm_score in results:
            tm_scores[(pdb1_id, pdb2_id)] = tm_score
            tm_scores[(pdb2_id, pdb1_id)] = tm_score  # Symmetric
        
        # Cache results
        with open(self.tm_scores_cache, 'wb') as f:
            pickle.dump(tm_scores, f)
            
        logger.info(f"Calculated {len(tm_scores)//2} pairwise TM-scores")
        return tm_scores
    
    def run_structure_clustering(self, samples, tm_scores):
        """Run graph-based structure clustering."""
        logger.info("Running structure clustering...")
        
        struct_clusters = {}
        
        # Only include samples that were used in TM-score calculation
        tm_sample_ids = set()
        for (pdb1, pdb2) in tm_scores.keys():
            tm_sample_ids.add(pdb1)
            tm_sample_ids.add(pdb2)
        
        for threshold in self.struct_thresholds:
            logger.info(f"Clustering at TM-score threshold {threshold}")
            
            # Create graph
            G = nx.Graph()
            G.add_nodes_from(tm_sample_ids)
            
            # Add edges based on TM-score threshold
            for (pdb1, pdb2), tm_score in tm_scores.items():
                if pdb1 < pdb2 and tm_score >= threshold:
                    G.add_edge(pdb1, pdb2)
                    
            # Find connected components (clusters)
            clusters = [list(component) for component in nx.connected_components(G)]
            struct_clusters[threshold] = clusters
            
            logger.info(f"Found {len(clusters)} clusters at threshold {threshold}")
            
        return struct_clusters
    
    def split_clusters(self, clusters, samples):
        """Split clusters into train/val/test sets."""
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
        start_time = time.time()
        logger.info("Starting parallel RNA clustering pipeline...")
        
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
                
        # Step 3: Structure clustering (optional)
        if not self.skip_structure:
            logger.info(f"Starting parallel structure clustering with max {self.max_struct_samples} samples...")
            tm_scores = self.calculate_tm_scores_parallel(samples)
            struct_clusters = self.run_structure_clustering(samples, tm_scores)
            
            # Process structure clustering results
            for threshold, clusters in struct_clusters.items():
                if clusters:
                    splits = self.split_clusters(clusters, samples)
                    self.save_splits(splits, "struct", threshold)
        else:
            logger.info("Skipping structure clustering (--skip-structure flag set)")
        
        elapsed_time = time.time() - start_time
        logger.info(f"Pipeline completed successfully in {elapsed_time:.2f} seconds!")

def main():
    parser = argparse.ArgumentParser(description="Parallel RNA Clustering and Data Splitting Pipeline")
    parser.add_argument("--pdb_dir", required=True, help="Directory containing PDB files")
    parser.add_argument("--fasta_dir", required=True, help="Directory containing FASTA files")
    parser.add_argument("--output_dir", required=True, help="Output directory for results")
    parser.add_argument("--max_seq_len", type=int, default=500, help="Maximum sequence length")
    parser.add_argument("--max_struct_samples", type=int, default=2000, 
                       help="Maximum number of samples for structure clustering")
    parser.add_argument("--skip_structure", action="store_true", 
                       help="Skip structure clustering (sequence only)")
    parser.add_argument("--n_processes", type=int, default=None,
                       help="Number of processes for parallel execution (default: CPU count - 1)")
    
    args = parser.parse_args()
    
    # Validate directories
    for dir_path in [args.pdb_dir, args.fasta_dir]:
        if not os.path.exists(dir_path):
            logger.error(f"Directory does not exist: {dir_path}")
            sys.exit(1)
            
    # Run pipeline
    pipeline = ParallelRNAClusteringPipeline(
        pdb_dir=args.pdb_dir,
        fasta_dir=args.fasta_dir,
        output_dir=args.output_dir,
        max_seq_len=args.max_seq_len,
        max_struct_samples=args.max_struct_samples,
        skip_structure=args.skip_structure,
        n_processes=args.n_processes
    )
    
    pipeline.run_pipeline()

if __name__ == "__main__":
    main()