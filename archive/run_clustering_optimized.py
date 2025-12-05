#!/usr/bin/env python3
"""
Optimized RNA Clustering Example Script

This script runs the optimized RNA clustering pipeline with sensible defaults
for large datasets. It includes options for sampling and skipping expensive
structure clustering if needed.

Usage:
    python run_clustering_optimized.py [--structure-samples N] [--skip-structure]
"""

import os
import sys
import argparse
import logging
from pathlib import Path

# Add utils directory to path
sys.path.append(str(Path(__file__).parent / "utils"))

from rna_clustering_pipeline_optimized import OptimizedRNAClusteringPipeline

def main():
    parser = argparse.ArgumentParser(description="Run optimized RNA clustering pipeline")
    parser.add_argument("--structure-samples", type=int, default=1000,
                       help="Maximum samples for structure clustering (default: 1000)")
    parser.add_argument("--skip-structure", action="store_true",
                       help="Skip structure clustering entirely")
    parser.add_argument("--max-seq-len", type=int, default=500,
                       help="Maximum sequence length (default: 500)")
    
    args = parser.parse_args()
    
    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('clustering_optimized.log'),
            logging.StreamHandler()
        ]
    )
    
    logger = logging.getLogger(__name__)
    
    # Define paths
    base_dir = Path("/Users/xiaojuzhang/Dev/pdb_prune/data/experiments_data/exp2.2_manuscript/dataset_20250731")
    pdb_dir = base_dir / "extracted_pdbs"
    fasta_dir = base_dir / "extracted_sequences"
    output_dir = base_dir / "clustering_results_optimized"
    
    # Validate input directories
    if not pdb_dir.exists():
        logger.error(f"PDB directory not found: {pdb_dir}")
        sys.exit(1)
        
    if not fasta_dir.exists():
        logger.error(f"FASTA directory not found: {fasta_dir}")
        sys.exit(1)
    
    # Log configuration
    logger.info("=" * 60)
    logger.info("OPTIMIZED RNA CLUSTERING PIPELINE")
    logger.info("=" * 60)
    logger.info(f"PDB directory: {pdb_dir}")
    logger.info(f"FASTA directory: {fasta_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Max sequence length: {args.max_seq_len}")
    logger.info(f"Max structure samples: {args.structure_samples}")
    logger.info(f"Skip structure clustering: {args.skip_structure}")
    logger.info("=" * 60)
    
    # Create and run pipeline
    try:
        pipeline = OptimizedRNAClusteringPipeline(
            pdb_dir=pdb_dir,
            fasta_dir=fasta_dir,
            output_dir=output_dir,
            max_seq_len=args.max_seq_len,
            max_struct_samples=args.structure_samples,
            skip_structure=args.skip_structure
        )
        
        pipeline.run_pipeline()
        
        logger.info("\n" + "=" * 60)
        logger.info("PIPELINE COMPLETED SUCCESSFULLY!")
        logger.info("=" * 60)
        
        # Show output files
        if output_dir.exists():
            output_files = list(output_dir.glob("*.txt"))
            if output_files:
                logger.info("\nGenerated output files:")
                for file in sorted(output_files):
                    logger.info(f"  - {file.name}")
            else:
                logger.warning("No output files found!")
        
    except KeyboardInterrupt:
        logger.info("\nPipeline interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Pipeline failed with error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()