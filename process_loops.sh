#!/bin/bash

# Create output directories
# mkdir -p ../extracted_loops_8.0/extracted_loop_npys
# mkdir -p ../extracted_loops_10.0/extracted_loop_npys
# mkdir -p ../extracted_loops_12.0/extracted_loop_npys
mkdir -p ../extracted_rna_segments_run1/extracted_loop_npys
mkdir -p ../extracted_rna_segments_run2/extracted_loop_npys

# Process 8.0Å loops
# echo "Processing 8.0Å loops..."
# python utils/pdb_to_npy.py --input-dir ./extracted_loops_8.0/extracted_loop_pdbs --output-dir ./extracted_loops_8.0/extracted_loop_npys

# # Process 10.0Å loops
# echo "Processing 10.0Å loops..."
# python utils/pdb_to_npy.py --input-dir ./extracted_loops_10.0/extracted_loop_pdbs --output-dir ./extracted_loops_10.0/extracted_loop_npys

# # Process 12.0Å loops
# echo "Processing 12.0Å loops..."
# python utils/pdb_to_npy.py --input-dir ./extracted_loops_12.0/extracted_loop_pdbs --output-dir ./extracted_loops_12.0/extracted_loop_npys

# Process RNA segments run 1
echo "Processing RNA segments run 1..."
python utils/pdb_to_npy.py --input-dir ./extracted_rna_segments_run1/extracted_rna_pdbs --output-dir ./extracted_rna_segments_run1/extracted_loop_npys

# Process RNA segments run 2
echo "Processing RNA segments run 2..."
python utils/pdb_to_npy.py --input-dir ./extracted_rna_segments_run2/extracted_rna_pdbs --output-dir ./extracted_rna_segments_run2/extracted_loop_npys

echo "All processing complete!" 