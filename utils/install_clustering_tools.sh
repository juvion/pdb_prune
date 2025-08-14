#!/bin/bash

# Installation script for RNA clustering pipeline dependencies
# This script installs PSI-CD-HIT and US-align tools

set -e  # Exit on any error

echo "Installing RNA clustering pipeline dependencies..."

# Create tools directory
TOOLS_DIR="$HOME/rna_tools"
mkdir -p "$TOOLS_DIR"
cd "$TOOLS_DIR"

# Install CD-HIT (includes PSI-CD-HIT)
echo "Installing CD-HIT (includes PSI-CD-HIT)..."
if ! command -v psicd-hit &> /dev/null && ! command -v cd-hit &> /dev/null; then
    # Try multiple installation methods
    if command -v conda &> /dev/null; then
        echo "Installing CD-HIT via conda..."
        conda install -c bioconda cd-hit -y
    elif command -v brew &> /dev/null; then
        echo "Installing CD-HIT via homebrew..."
        brew install cd-hit
    else
        # Download and compile from source
        echo "Compiling CD-HIT from source..."
        if [ ! -d "cd-hit" ]; then
            git clone https://github.com/weizhongli/cdhit.git cd-hit
        fi
        cd cd-hit
        make
        
        # Add to PATH
        echo "export PATH=\$PATH:$TOOLS_DIR/cd-hit" >> ~/.bashrc
        echo "export PATH=\$PATH:$TOOLS_DIR/cd-hit" >> ~/.zshrc
        
        cd ..
    fi
    echo "CD-HIT installed successfully!"
else
    echo "CD-HIT already installed."
fi

# Install US-align
echo "Installing US-align..."
if ! command -v US-align &> /dev/null; then
    # Download US-align
    if [ ! -f "USalign" ]; then
        curl -O https://zhanggroup.org/US-align/bin/module/USalign.cpp
        # Use dynamic linking for macOS compatibility
        g++ -O3 -ffast-math -o USalign USalign.cpp
    fi
    
    # Create symlink with expected name
    ln -sf "$TOOLS_DIR/USalign" "$TOOLS_DIR/US-align"
    
    # Add to PATH
    echo "export PATH=\$PATH:$TOOLS_DIR" >> ~/.bashrc
    echo "export PATH=\$PATH:$TOOLS_DIR" >> ~/.zshrc
    
    echo "US-align installed successfully!"
else
    echo "US-align already installed."
fi

echo ""
echo "Installation complete!"
echo "Please run 'source ~/.bashrc' or restart your terminal to update PATH."
echo ""
echo "To verify installation, run:"
echo "  psicd-hit -h"
echo "  US-align -h"
echo ""
echo "To run the clustering pipeline:"
echo "  python utils/rna_clustering_pipeline.py --pdb_dir <pdb_dir> --fasta_dir <fasta_dir> --output_dir <output_dir>"