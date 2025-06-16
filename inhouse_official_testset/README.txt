Dataset Construction Process: Filtered Official Test In-house Dataset

This dataset was constructed through the following steps:

1. Initial PDB Download
   - Downloaded original PDB files from the official PDB list
   - Source: RCSB Protein Data Bank (PDB)
   - Files are stored in the 'raw_pdbs' directory

2. RNA Chain Extraction
   - Processed original PDB files to extract single RNA chains
   - Removed all non-RNA chains and multi-chain RNA structures
   - Preserved only single-stranded RNA chains
   - Extracted chains are stored in 'processed_pdbs' directory
   - Each file is named as 'PDB_id_chain_id.pdb'

3. Sequence Matching and Filtering
   - Matched PDB_id, length, and chain_id between datasets
   - Used case-insensitive matching for PDB_id and chain_id
   - Implemented matching logic:
     * Perfect match: both chain_id and length match
     * Chain mismatch: chain_id matches but length doesn't
     * Length mismatch: chain_id doesn't match but length does
     * PDB mismatch: both chain_id and length don't match
   - Matching results stored in 'merged_results_match_filtered.csv'

4. Manual Verification and Cleaning
   - Reviewed and resolved ambiguous cases:
     * One-to-many matches where PDB+length match but chain_id differs
     * Cases where PDB+chain match but length differs
   - Removed problematic entries to ensure data quality
   - Final filtered list stored in 'merged_results_match_filtered.csv'

5. Data Processing
   - Converted PDB files to numpy arrays (.npy format)
   - Preserved structural information and coordinates
   - Processed files stored in the dataset directory

Directory Structure:
- raw_pdbs/: Original PDB files
- processed_pdbs/: Extracted single RNA chains
- processed_rna_sequences/: FASTA files of RNA sequences
- destination_folder/: Final filtered dataset
  * pdb/: Processed PDB files
  * fasta/: Corresponding FASTA files

File Naming Convention:
- PDB files: PDB_id_chain_id.pdb
- FASTA files: PDB_id_chain_id.fasta
- All PDB_ids are in lowercase

Note: This dataset represents a carefully curated subset of RNA structures, ensuring high-quality single-chain RNA data for testing and validation purposes. 