


          
Here’s a clear prompt you can copy and use to describe exactly what you want the Python script to do:

Objective:
- Implement a Python script that extracts and merges specific chains from PDB files according to a CSV specification, then saves a single combined PDB file and a single combined FASTA file per CSV row.

Inputs:
- Source PDB directory: /Users/xiaojuzhang/Dev/pdb_prune/experiment/exp3/validation_downloads
- CSV file: /Users/xiaojuzhang/Dev/pdb_prune/experiment/exp3/validation_pdbid.csv
- Each CSV row has three columns: PDB_ID (e.g., 4KQ0), Model_Index (e.g., 1), Chain_Spec (e.g., B-E)
- The script should treat PDB_ID case-insensitively when locating files (e.g., open 4kq0.pdb for 4KQ0)

Chain specification rules:
- Chain_Spec can be a single chain (e.g., A), a comma-separated list (e.g., A,B,E), or an inclusive range with a hyphen (e.g., B-E)
- Inclusive range means: include every single-letter chain between the start and end, e.g., B-E means B, C, D, and E

Processing details:
- For each CSV row, open the corresponding PDB file from the source directory (e.g., 4KQ0 → /Users/xiaojuzhang/Dev/pdb_prune/experiment/exp3/validation_downloads/4kq0.pdb)
- If the file has MODEL records and a Model_Index is provided, only extract records from that model; otherwise, process the file as-is
- Extract ATOM and HETATM records whose chain ID (column 22 in PDB format) matches any chain from the Chain_Spec set
- Preserve TER records appropriately and maintain the original record order for extracted lines
- Preserve header/meta records (HEADER, TITLE, COMPND, REMARK) or at least add a minimal header/REMARK indicating which chains were extracted
- Case sensitivity: keep the original chain IDs from the PDB file; only PDB file names are case-insensitive

Sequence generation:
- Generate a FASTA file that includes the sequences for all extracted chains in the same order, as a multi-FASTA:
  - One FASTA entry per chain, with header “>PDBID_MODEL_CHAIN” (e.g., >4KQ0_1_B)
  - For protein chains: convert 3-letter residue names to standard 1-letter amino acid codes; skip unknowns or map to X
  - For nucleic acid chains: use A, U, G, C, T; map DNA residues (DA, DU, DG, DC, DT) to their single-letter equivalents; ignore modified bases or map them conservatively
  - Only count residues from ATOM records (skip HETATM for sequence unless clearly nucleotides)
- The single FASTA file per row should contain multiple entries (one per chain), not concatenated into one sequence

Outputs:
- Create output directories if they do not exist: PDBs and SEQs
- Write one combined PDB file per CSV row to PDBs with the filename: PDBID_MODEL_CHAINSPEC.pdb (e.g., 4KQ0_1_B-E.pdb)
- Write one combined FASTA file per CSV row to SEQs with the filename: PDBID_MODEL_CHAINSPEC.fasta (e.g., 4KQ0_1_B-E.fasta)
- Use uppercase PDB IDs in output filenames

Example:
- CSV row: 4KQ0,1,B-E
- Chain_Spec B-E means: extract chains B, C, D, and E
- Source file: /Users/xiaojuzhang/Dev/pdb_prune/experiment/exp3/validation_downloads/4kq0.pdb
- Output files:
  - /Users/xiaojuzhang/Dev/pdb_prune/experiment/exp3/PDBs/4KQ0_1_B-E.pdb containing ATOM/HETATM lines for chains B, C, D, and E
  - /Users/xiaojuzhang/Dev/pdb_prune/experiment/exp3/SEQs/4KQ0_1_B-E.fasta containing four FASTA entries: >4KQ0_1_B, >4KQ0_1_C, >4KQ0_1_D, >4KQ0_1_E

Error handling and logging:
- If a PDB file is missing, log a warning and continue
- If any specified chain does not exist in the PDB file, log it but still produce outputs for available chains
- Provide a summary at the end: total processed, total errors
- Log to both console and a file (e.g., pdb_chain_extraction.log)

Assumptions:
- The CSV has no header row and is in the format: PDB_ID,Model_Index,Chain_Spec
- The script runs in Python 3 and should avoid heavy dependencies; standard library is preferred, but Bio.PDB can be used if available for robustness

Acceptance criteria:
- For the example row 4KQ0,1,B-E, the script produces exactly two files named 4KQ0_1_B-E.pdb and 4KQ0_1_B-E.fasta, with the PDB containing extracted ATOM/HETATM records for B, C, D, and E, and the FASTA containing separate entries for each of those chains.
        