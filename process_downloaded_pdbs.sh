#!/bin/bash

# Create output directories
mkdir -p downloaded_coords
mkdir -p downloaded_seqs
mkdir -p temp_pdbs

# Convert .ent files to .pdb files
echo "Converting .ent files to .pdb format..."
for ent_file in downloaded_rna_pdbs/*.ent; do
    if [ -f "$ent_file" ]; then
        pdb_name=$(basename "$ent_file" .ent)
        cp "$ent_file" "temp_pdbs/${pdb_name}.pdb"
    fi
done

# Process PDB files to NPY format
echo "Converting PDB files to NPY format..."
python utils/pdb_to_npy.py --input-dir temp_pdbs --output-dir downloaded_coords

# Process PDB files to FASTA format
echo "Converting PDB files to FASTA format..."
python utils/pdb_to_fasta.py --pdb-dir temp_pdbs --output-dir downloaded_seqs

# Clean up temporary files
echo "Cleaning up temporary files..."
rm -rf temp_pdbs

echo "All processing complete!" 