RNA Loop and Segment Extraction Process
=====================================

This document describes the extraction process and contents of the following folders:
- extracted_loops_8.0/
- extracted_loops_10.0/
- extracted_loops_12.0/
- extracted_rna_segments_run1/
- extracted_rna_segments_run2/

1. Loop Extraction (extracted_loops_*.0/)
----------------------------------------
These folders contain RNA loops extracted using different distance cutoffs (8.0Å, 10.0Å, and 12.0Å).

Process:
- Input: PDB files from reconstructed_pdbs/
- Method: LoopExtractor class (utils/extract_loops.py)
- Parameters:
  * min_length: 3 residues
  * max_length: 20 residues
  * distance_cutoff: 8.0Å, 10.0Å, or 12.0Å
  * atom_type: "C4'" (sugar-phosphate backbone)
- Output Structure:
  * extracted_loop_pdbs/: Contains PDB files of extracted loops
  * extracted_loop_seqs/: Contains FASTA files of loop sequences

File Naming Convention:
- PDB files: [original_pdb]_[chain_id]_[start_res]-[end_res][T].pdb
  * T suffix indicates terminal loops
- FASTA files: [original_pdb]_[chain_id]_[start_res]-[end_res][T].fasta
  * Header format: >[pdb_id]|[chain_id]|[start_res]-[end_res]|Length=[length]

2. RNA Segment Extraction (extracted_rna_segments_run1/)
-----------------------------------------------------
This folder contains continuous RNA segments extracted from the PDB files.

Process:
- Input: PDB files from reconstructed_pdbs/
- Method: RNAExtractor class (utils/extract_rna_segments.py)
- Parameters:
  * min_length: 3 residues
  * max_length: 30 residues
  * coverage_rate: 0.02 (2% of possible segments are sampled)
- Output Structure:
  * extracted_segment_pdbs/: Contains PDB files of extracted segments
  * extracted_segment_seqs/: Contains FASTA files of segment sequences

File Naming Convention:
- PDB files: [pdb_id]_[chain_ids]_[start_res]-[end_res][T].pdb
  * T suffix indicates terminal segments
- FASTA files: [pdb_id]_[chain_ids]_[start_res]-[end_res][T].fasta
  * Header format: >[pdb_id]|[chain_ids]|[start_res]-[end_res]|Length=[length]

3. RNA Segment Extraction (extracted_rna_segments_run2/)
-----------------------------------------------------
This folder contains continuous RNA segments extracted from the PDB files.

Process:
- Input: PDB files from reconstructed_pdbs/
- Method: RNAExtractor class (utils/extract_rna_segments.py)
- Parameters:
  * min_length: 3 residues
  * max_length: 20 residues
  * coverage_rate: 0.05 (5% of possible segments are sampled)
- Output Structure:
  * extracted_segment_pdbs/: Contains PDB files of extracted segments
  * extracted_segment_seqs/: Contains FASTA files of segment sequences

File Naming Convention:
- PDB files: [pdb_id]_[chain_ids]_[start_res]-[end_res][T].pdb
  * T suffix indicates terminal segments
- FASTA files: [pdb_id]_[chain_ids]_[start_res]-[end_res][T].fasta
  * Header format: >[pdb_id]|[chain_ids]|[start_res]-[end_res]|Length=[length]

Notes:
1. The distance cutoff (8.0Å, 10.0Å, 12.0Å) affects how loops are identified:
   - Smaller cutoff (8.0Å): More conservative, shorter loops
   - Larger cutoff (12.0Å): More lenient, potentially longer loops
   - Medium cutoff (10.0Å): Balanced approach

2. Loop vs Segment Extraction:
   - Loops: Based on spatial distance between residues
   - Segments: Based on continuous chain regions with random sampling

3. Quality Control:
   - Only RNA residues (A, U, G, C) are included
   - Minimum length ensures meaningful structures
   - Maximum length prevents overly long fragments
   - Terminal markers (T) help identify end regions

4. File Organization:
   - Each folder contains both PDB and FASTA files
   - PDB files preserve 3D structure information
   - FASTA files provide sequence information
   - Consistent naming scheme across all extractions 