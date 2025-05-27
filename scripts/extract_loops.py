#!/usr/bin/env python3

import os
import sys
import logging
from pathlib import Path
from utils.extract_loops import LoopExtractor

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def main():
    # Example: Extract loops with different cutoffs
    for cutoff in [8.0, 10.0, 12.0]:
        extractor = LoopExtractor(
            input_dir="reconstructed_pdbs",
            min_length=3,
            max_length=20,
            distance_cutoff=cutoff,
            atom_type="C4'",
            generation_id=f"cutoff_{cutoff}"  # Use generation_id instead of output_dir
        )
        extractor.process()

if __name__ == "__main__":
    main() 