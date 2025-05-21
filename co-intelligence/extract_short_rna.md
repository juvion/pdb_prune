“Help me to implement this Python script to randomly extract short RNA pdbs from the existing RNA pdbs, which should follow these rules: 0. cover all the pdbs from the files, 1. Randomly extract different lengths of daughter RNAs from each PDB file, 2. residues should be within the same chain as a single strand, 3. residues' index should be continuous, 4.  The length of extracted RNA should be less than the given limit (default 25)， 5. Randomly choose the length. 5. output format should be one pdb and its Fasta corresponding file ”  Please help me to polish this prompt for cursor to be accurate and effective.


Hello! I need a Python script to **randomly extract short RNA segments (referred to as "daughter RNAs")** from a collection of existing RNA PDB files. The script must meticulously follow specific rules for extraction and, crucially, **output each extracted segment as a new PDB file AND contribute its sequence to a single, collective FASTA file**.

Here are the detailed requirements for the script:

1.  **Input Source:** The script should take two command-line arguments:
    * **`input_directory`**: The path to a directory that contains the original RNA PDB files (e.g., `/path/to/my_rna_pdbs/`).
    * **`generation_id`**: A unique string identifier for this specific run of the script (e.g., "run_20240520"). This ID will be incorporated into the output filenames and the FASTA header.

    The script must **process every single valid PDB file** (e.g., with `.pdb` or `.ent` extensions) found within the specified `input_directory`.

2.  **Extraction Logic & Rules for Daughter RNAs:**
    * **Processing All PDBs:** The script must iterate through all PDB files in the `input_directory`.
    * **Per-PDB Extraction Count:** From each processed PDB file, the script should attempt to extract a specified number of random daughter RNA segments per *eligible logical RNA chain*. Let's make this count a configurable parameter, defaulting to `num_extractions_per_chain = 5`. If fewer valid segments can be found, it should extract all possible ones.
    * **Random Length Selection:** For each daughter RNA segment being extracted, its exact length must be **randomly chosen**.
    * **Maximum Length Limit:** The length of any extracted daughter RNA segment must be **less than a specified maximum limit**. This limit should be a configurable parameter, with a default value of `max_length = 25` nucleotides.
    * **Residue Index Continuity:** The residues within any extracted segment must have **absolutely continuous PDB residue indexing**. This means that within the segment's span, all residue numbers (as annotated in the PDB file) must be present sequentially without any gaps (e.g., if a segment runs from residue `10` to `15`, residues `10, 11, 12, 13, 14, 15` must all exist in the original PDB model and chain).
    * **Single Logical Chain Constraint:** Each extracted daughter RNA segment must originate entirely from a **single, continuous logical RNA chain** within the PDB structure. This implies:
        * It must **not** span across different original chain IDs *unless* those chains are identified as part of a single biological RNA molecule due to **continuous residue numbering across `TER` records/chain breaks** (as per our previous discussions).
        * If multiple original chains (e.g., Q and R) are identified as one continuous logical chain, segments *can* span across their original boundaries, and the full range of original chain IDs should be reported.

3.  **Output Formats and Naming (Critical Dual Output):**
    * The script should create a **new main output directory**, ideally named based on the `generation_id` (e.g., `extracted_rna_segments_run_20240520`).
    * Inside this main output directory, it will contain:
        * **A subdirectory for individual PDB files:** (e.g., `extracted_rna_pdbs/` within the main output directory).
        * **A single, collective FASTA file:** (e.g., `all_extracted_rna_sequences.fasta` within the main output directory).

    **For each extracted daughter RNA segment, both of the following must happen:**

    * **A. Individual PDB File Output:**
        * The segment must be saved as a **separate PDB file** within the `extracted_rna_pdbs/` subdirectory.
        * The naming convention for each PDB file must be exact:
            `{generation_id}_{initial_pdb_code}_{chain_ID(s)}_{increment}.pdb`
            * **`{generation_id}`**: The unique string identifier provided as a command-line argument.
            * **`{initial_pdb_code}`**: The 4-character PDB ID (e.g., "1ABC") of the original PDB file from which the segment was extracted.
            * **`{chain_ID(s)}`**: Represents the original chain ID(s) from which the segment was extracted.
                * If the segment comes from a single original chain, use its ID (e.g., "A").
                * If it spans across multiple original chains that were identified as a continuous logical chain (e.g., Q and R), concatenate their chain IDs (e.g., "QR").
            * **`{increment}`**: A sequential counter, unique for each daughter RNA segment extracted from the *same original PDB file and the same logical chain*. For instance, if `1ABC`'s logical chain `A` yields three segments, they would be named `..._A_1.pdb`, `..._A_2.pdb`, `..._A_3.pdb`. This counter should reset for each new logical chain/PDB file.
        * **Example PDB Filenames:**
            * `run_20240520_1ABC_A_1.pdb`
            * `run_20240520_1A4D_QR_2.pdb`

    * **B. Single Collective FASTA File Output:**
        * The sequence of the extracted segment must be appended to the **single output FASTA file** (`all_extracted_rna_sequences.fasta`).
        * The header for each entry in this FASTA file must clearly identify its origin and details. The format should be:
            `>{PDB_ID}|{Chain_ID(s)}|{Start_Residue_Index}-{End_Residue_Index}|Length={Extracted_Length}`
            * **`{PDB_ID}`**: The 4-character PDB ID of the original file.
            * **`{Chain_ID(s)}`**: As defined for the PDB filename (single ID or concatenated IDs for continuous logical chains).
            * **`{Start_Residue_Index}-{End_Residue_Index}`**: The PDB residue index (not a 0-based Python index) of the first and last residue of the extracted segment.
            * **`{Extracted_Length}`**: The actual length of the extracted segment.
        * **Example FASTA Headers:**
            * `>1ABC|A|10-35|Length=26`
            * `>1A4D|QR|15-40|Length=26`

4.  **Error Handling and Informative Output:**
    * The script should gracefully handle PDB files that are invalid, do not contain any RNA chains, or whose RNA chains are too short to yield any segments meeting the length and continuity criteria.
    * Informative messages should be printed to the console throughout the process (e.g., "Processing PDB X...", "Identified Y logical RNA chains...", "Extracted Z segments from chain A, saved to ...").
    * Summarize how many PDBs were processed, how many daughter RNAs were extracted, and any encountered errors.

5.  **Dependencies:**
    * The script should utilize the `Biopython` library (specifically `Bio.PDB` for parsing PDB files, identifying RNA components, writing new PDB files, and `Bio.SeqIO` for FASTA output).
    * It will also need Python's built-in `random` module for random selection and `os` for file/directory operations.

Please provide a comprehensive, step-by-step implementation plan for this Python script. For each major phase (e.g., script setup, parsing arguments, PDB file iteration, identifying eligible continuous RNA chains, random segment extraction, constructing new PDB structures, writing new PDB files, appending to FASTA, and final cleanup), include:

* A clear explanation of the logic and what the code aims to achieve.
* The corresponding Python code snippet.
* Detailed instructions on how to integrate these code components and how to execute the overall script from the command line.
* Specific guidance on how to implement both the PDB writing using `Bio.PDB.PDBIO` and the FASTA appending using `Bio.SeqIO`, ensuring the exact filename/header generation as specified.

Ensure the code is clear, well-commented, and easy to understand.