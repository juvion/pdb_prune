Hello! I need a Python script that can **reconstruct RNA PDB files** from pairs of NumPy coordinate array files (.npy) and corresponding FASTA sequence files (.fasta). This script will essentially reverse a previous processing step, creating PDBs that contain only the specific atoms stored in the NumPy arrays.

Here are the detailed specifications for the script:

1.  **Input Data Structure:**
    * The script will take two command-line arguments:
        * **`coords_directory`**: Path to a directory containing NumPy array files (e.g., `/path/to/my_training_data/coords/`). These files will have a `.npy` extension.
        * **`seqs_directory`**: Path to a directory containing FASTA sequence files (e.g., `/path/to/my_training_data/seqs/`). These files will have a `.fasta` extension.
    * For each RNA molecule, there will be a corresponding `.npy` file and a `.fasta` file sharing the exact same base name (e.g., `rna_001.npy` and `rna_001.fasta`).

2.  **Input File Formats:**
    * **`.fasta` files:** Each `.fasta` file will contain a single RNA sequence. This sequence directly provides the base sequence (A, U, G, C) and thus the length (`i`) of the RNA molecule.
    * **`.npy` files:** Each `.npy` file contains a 3D NumPy array with the dimensions `(i, j, k)`:
        * `i`: Corresponds to the number of bases/residues in the RNA sequence. This length *must match* the length of the corresponding sequence in the `.fasta` file.
        * `j`: Is fixed to **7**. These 7 elements represent the coordinates for specific atoms in a defined order for each residue:
            1.  `P` (Phosphorus atom of the phosphate group)
            2.  `O5'` (Oxygen atom bonded to C5' in the sugar)
            3.  `C5'` (Carbon atom in the sugar)
            4.  `C4'` (Carbon atom in the sugar)
            5.  `C3'` (Carbon atom in the sugar)
            6.  `O3'` (Oxygen atom bonded to C3' in the sugar)
            7.  `N1` (for pyrimidines: U, C) **OR** `N9` (for purines: A, G) (the key base atom)
        * `k`: Is fixed to **3**. These represent the `(x, y, z)` spatial coordinates for each of the `j` atoms.

3.  **Output:**
    * The script should create a **new output directory** (e.g., `reconstructed_pdbs/`).
    * For each successfully reconstructed RNA molecule, it should generate a **new PDB file**.
    * The output PDB file should share the same base name as its input `.npy` and `.fasta` files (e.g., `rna_001.pdb`).

4.  **Core Reconstruction Logic:**
    * **File Pairing:** The script must iterate through the base names (e.g., `rna_001`) common to both directories, ensuring a `.npy` and a `.fasta` file exist for each.
    * **Data Loading:** For each pair, load the NumPy array and read the RNA sequence from the FASTA file.
    * **Length Validation:** Crucially, validate that the first dimension (`i`) of the NumPy array matches the length of the RNA sequence obtained from the FASTA file. If they don't match, log an error and skip the pair.
    * **PDB Structure Construction (using Biopython):**
        * Create a new `Bio.PDB.Structure.Structure` object.
        * Add a single `Model` (e.g., with ID 0).
        * Add a single `Chain` to this model. Let's make its ID configurable, with a default of **`'A'`**.
        * **For each residue** in the RNA sequence (from the FASTA file):
            * Create a `Bio.PDB.Residue.Residue` object.
            * Assign the correct residue name based on the base (e.g., 'A', 'U', 'G', 'C').
            * Assign a sequential residue ID (e.g., starting from 1, configurable).
            * **Add the 7 atoms** (`P`, `O5'`, `C5'`, `C4'`, `C3'`, `O3'`, and `N1/N9`) to this residue, taking their `(x, y, z)` coordinates directly from the corresponding slice of the NumPy array.
            * **Crucially, assign the correct atom name for the 7th atom:** `N1` if the base is 'U' or 'C', and `N9` if the base is 'A' or 'G'.
            * Assign appropriate atom IDs, element types (`'P'`, `'O'`, `'C'`, `'N'`), default occupancy (e.g., `1.0`), and default B-factor (e.g., `20.0`) for each atom.
    * **PDB File Writing:** Use `Bio.PDB.PDBIO` to write the newly constructed Biopython `Structure` object to a `.pdb` file in the output directory.

5.  **Important Limitation & Clarification:**
    * The reconstructed PDB files will **only contain the 7 specified backbone and base atoms** for each residue, as this is the only coordinate information available in the `.npy` files. They will not be full PDBs with all standard atoms (e.g., base carbons/nitrogens beyond N1/N9, sugar carbons/oxygens beyond the backbone atoms listed). This limitation should be implicitly understood in the output.

6.  **Error Handling and Informative Output:**
    * The script should gracefully handle cases where a `.npy` file doesn't have a matching `.fasta` file (and vice-versa).
    * It should log errors for length mismatches between `.npy` and `.fasta` files.
    * Informative messages should be printed to the console (e.g., "Processing rna_001...", "Reconstructing rna_001.pdb", "Error: Length mismatch for rna_xyz").
    * Provide a final summary of how many PDBs were successfully reconstructed.

7.  **Dependencies:**
    * The script will require the `biopython` library (specifically `Bio.PDB` for structure creation/writing and `Bio.SeqIO` for FASTA parsing).
    * It will require the `numpy` library for loading `.npy` files.
    * It will need Python's built-in `os` and `argparse` (for command-line arguments) modules.

Please provide a comprehensive, step-by-step implementation plan for this Python script. For each major phase (e.g., script setup, argument parsing, file pairing, data loading, structure construction, PDB writing, and error handling), include:

* A clear explanation of the logic and what the code aims to achieve.
* The corresponding Python code snippet.
* Detailed instructions on how to integrate these code components and how to execute the overall script from the command line.
* Specific guidance on implementing the Biopython `Structure`, `Model`, `Chain`, `Residue`, and `Atom` objects, and the `PDBIO` writer.

Ensure the code is clear, well-commented, and easy to understand.


python reconstruct_pdb.py competition/train/coords/ competition/train/seqs/ reconstructed_pdbs/