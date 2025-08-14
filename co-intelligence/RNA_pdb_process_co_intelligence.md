Please write a Python script that processes an RNA structure from a PDB file, identifies loop fragments based on a 3D contact analysis, and saves each identified loop fragment as a separate PDB file.

**Goal:** Identify continuous segments of the RNA chain that act as loops (lacking significant 3D contacts with distant residues) and save each loop as a new PDB file.

**Input:**
* A string representing the file path to an RNA structure in PDB format.

**Output:**
* Multiple PDB files, each containing the atomic coordinates for one identified loop fragment.
* The output files should be named following this convention: `[original_pdb_code]_[start_position]-[end_position].pdb`. The `start_position` and `end_position` should correspond to the 1-based residue numbering in the original PDB file.

**Algorithm to Implement:**

1.  **Load the PDB file:** Parse the input PDB file using a library like Biopython's `PDBParser`. Obtain the structure and identify the RNA chain(s). For simplicity, you can focus on the first model and the first RNA chain found, or make it adaptable if the user specifies a chain ID. Extract the original PDB code from the file name or header if possible.
2.  **Define a distance cutoff:** Set a constant variable for the neighbor distance cutoff (10 Angstroms).
3.  **Compute 3D Neighbors using N1/N9 distance:** For each residue in the RNA chain, determine which other residues are within the defined 10 Angstrom cutoff distance.
    * **Specific Distance Metric:** Calculate the distance between the **N1 atom** (for Adenine 'A' and Guanine 'G') or the **N9 atom** (for Cytosine 'C' and Uracil 'U') of the current residue and the corresponding N1/N9 atom of every other residue. If a residue is not A, G, C, or U, or the N1/N9 atom is missing, handle this appropriately (e.g., skip distance calculation for that residue pair or use an alternative representative atom like C4' if necessary, but prioritize N1/N9).
    * Store these neighbors for each residue. A dictionary where keys are residue objects (or their identifiers/indices) and values are lists of neighbor residue objects (or identifiers/indices) is suitable.
4.  **Identify "Loop-like" Residues:** Iterate through each residue in the RNA chain based on its sequential position. A residue at sequential position `i` (e.g., 0-indexed) is considered "loop-like" if and only if *none* of its 3D neighbors (identified in step 3) have a sequential position `j` such that the absolute difference `abs(i - j)` is greater than 3. The sequential position should reflect the order in the RNA chain.
5.  **Extract and Save Continuous Loop Fragments:** Identify consecutive runs of residues that were marked as "loop-like" in step 4. These continuous segments are the loop fragments.
    * Iterate through the chain's residues sequentially.
    * When you find a "loop-like" residue that is not the sequential successor of the previous "loop-like" residue, it marks the start of a new fragment.
    * Continue adding sequential "loop-like" residues to the current fragment.
    * When you encounter a non-"loop-like" residue or reach the end of the chain, the current fragment ends.
    * For each identified fragment containing one or more residues:
        * Create a new Biopython `Structure`, `Model`, and `Chain` containing only the residues belonging to this fragment. Include all atoms for these residues.
        * Get the 1-based residue sequence numbers from the original PDB for the start and end residues of the fragment.
        * Construct the output filename using the original PDB code (e.g., from the input filename like `1abc.pdb` -> `1abc`), the start position, and the end position: `[pdb_code]_[start]-[end].pdb`.
        * Save the newly created structure containing the fragment to a new PDB file with the generated filename using `Biopython.PDB.PDBIO`.

**Implementation Notes:**

* Use the Biopython library (`Bio.PDB`). You'll need `PDBParser`, `NeighborSearch`, and `PDBIO`.
* Carefully handle residue indexing. You'll likely need to map Biopython's potentially non-sequential residue numbering (`res.get_id()[1]`) to a 0-based or 1-based sequential index within the processed chain to correctly apply the `abs(i - j) <= 3` rule. The output filename requires the original PDB residue numbers.
* Ensure robust handling of finding the N1/N9 atoms based on residue name (A, G, C, U).
* Consider adding a check to make sure the input file is a PDB and contains an RNA chain.
* Specify where the output PDB files should be saved (e.g., in the current directory or a specified output directory). Saving in the current directory is a reasonable default.

Please provide the complete Python script implementing this algorithm.