Hello! I need help creating a Python script using the Biopython library to perform a specific bioinformatics task on PDB files containing RNA.

The script should accomplish the following steps:
0. **Retrieve PDB files in batch:** The script should be able to fetch all PDB file (by random sampling and download). It needs to handle cases where the PDB might contain proteins or other molecules in addition to RNA.
1.  OR **Retrieve PDB files with ID:** The script should be able to fetch a PDB file given its PDB ID. It needs to handle cases where the PDB might contain proteins or other molecules in addition to RNA.
2.  **Isolate RNA chains:** After parsing the PDB file, the script must identify and extract *only* the RNA chains present in the structure. It should be able to separate one single chain for each file, and remove the duplicate sequences.
3.  **Extract specific atoms from each RNA chain:** For each individual RNA chain identified, the script should iterate through its residues and extract a defined set of atoms. This set should include the backbone atoms (P, O5', C5', C4', C3', O3', C2', O2', C1') and the base atom N1 (for pyrimidines like Uridine 'U' and Cytidine 'C') or N9 (for purines like Adenosine 'A' and Guanine 'G').
4.  **Output the extracted atom data as the PDB format:** For each extracted atom, the script should output relevant information such as the chain ID, residue number, residue name, atom name, and its coordinates (x, y, z).

Please provide a step-by-step guide on how to write this script using Biopython's PDB module.

For each step, please:
* Explain the purpose of the step.
* Provide the Python code snippet required for that step.
* Explain how the code works.
* Include basic validation or error handling (e.g., how to check if a PDB ID is valid, or if no RNA is found).

Please keep the code simple, straightforward, and easy to understand. The focus is on functionality and readability, not over-engineering.

Thank you!