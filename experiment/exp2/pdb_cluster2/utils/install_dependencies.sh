#!/bin/bash
###############################################################################
# Installation Script for RNA PDB Clustering Pipeline
# 
# This script helps install all required dependencies
###############################################################################

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

print_step() {
    echo -e "${BLUE}[STEP]${NC} $1"
}

print_success() {
    echo -e "${GREEN}✓${NC} $1"
}

print_error() {
    echo -e "${RED}✗${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}⚠${NC} $1"
}

echo "================================================================================"
echo "          RNA PDB Clustering Pipeline - Dependency Installer"
echo "================================================================================"
echo ""

# Check if conda is available
if command -v conda &> /dev/null; then
    USE_CONDA=true
    print_success "Conda detected"
else
    USE_CONDA=false
    print_warning "Conda not found, will use pip for Python packages"
fi

# Install Python packages
print_step "Installing Python packages..."

if [ "$USE_CONDA" = true ]; then
    echo "Using conda..."
    conda install -y -c conda-forge biopython numpy scipy pandas matplotlib tqdm
else
    echo "Using pip..."
    pip install --user biopython numpy scipy pandas matplotlib tqdm
fi

if [ $? -eq 0 ]; then
    print_success "Python packages installed"
else
    print_error "Failed to install Python packages"
    exit 1
fi
echo ""

# Install CD-HIT
print_step "Installing CD-HIT..."

if [ "$USE_CONDA" = true ]; then
    echo "Using conda..."
    conda install -y -c bioconda cd-hit
    
    if [ $? -eq 0 ]; then
        print_success "CD-HIT installed via conda"
    else
        print_warning "Failed to install via conda, trying from source..."
        USE_CONDA=false
    fi
fi

if [ "$USE_CONDA" = false ]; then
    echo "Installing from source..."
    
    # Check if already installed
    if command -v cd-hit-est &> /dev/null; then
        print_warning "CD-HIT already installed"
    else
        mkdir -p /tmp/cdhit_install
        cd /tmp/cdhit_install
        
        git clone https://github.com/weizhongli/cdhit.git
        cd cdhit
        make
        
        # Try to install to user's bin
        mkdir -p ~/bin
        cp cd-hit cd-hit-est psi-cd-hit ~/bin/
        
        # Add to PATH
        if ! grep -q '~/bin' ~/.bashrc; then
            echo 'export PATH=$PATH:~/bin' >> ~/.bashrc
        fi
        export PATH=$PATH:~/bin
        
        cd ~
        rm -rf /tmp/cdhit_install
        
        print_success "CD-HIT installed to ~/bin"
        print_warning "Please run: source ~/.bashrc"
    fi
fi
echo ""

# Install US-align
print_step "Installing US-align..."

if command -v USalign &> /dev/null; then
    print_warning "US-align already installed"
else
    mkdir -p /tmp/usalign_install
    cd /tmp/usalign_install
    
    git clone https://github.com/pylelab/USalign.git
    cd USalign
    g++ -static -O3 -ffast-math -lm -o USalign USalign.cpp
    
    # Try to install to user's bin
    mkdir -p ~/bin
    cp USalign ~/bin/
    
    # Add to PATH
    if ! grep -q '~/bin' ~/.bashrc; then
        echo 'export PATH=$PATH:~/bin' >> ~/.bashrc
    fi
    export PATH=$PATH:~/bin
    
    cd ~
    rm -rf /tmp/usalign_install
    
    print_success "US-align installed to ~/bin"
    print_warning "Please run: source ~/.bashrc"
fi
echo ""

# Verify installations
print_step "Verifying installations..."

# Python packages
echo "Checking Python packages..."
python3 -c "import Bio; print(f'  BioPython {Bio.__version__}')" || print_error "BioPython not found"
python3 -c "import numpy; print(f'  NumPy {numpy.__version__}')" || print_error "NumPy not found"
python3 -c "import pandas; print(f'  Pandas {pandas.__version__}')" || print_error "Pandas not found"
python3 -c "import scipy; print(f'  SciPy {scipy.__version__}')" || print_error "SciPy not found"
python3 -c "import matplotlib; print(f'  Matplotlib {matplotlib.__version__}')" || print_error "Matplotlib not found"
python3 -c "import tqdm; print(f'  tqdm {tqdm.__version__}')" || print_error "tqdm not found"

echo ""
echo "Checking command-line tools..."

if command -v cd-hit-est &> /dev/null; then
    print_success "cd-hit-est found: $(which cd-hit-est)"
else
    print_error "cd-hit-est not found"
fi

if command -v USalign &> /dev/null; then
    print_success "USalign found: $(which USalign)"
else
    print_error "USalign not found"
fi

echo ""
echo "================================================================================"
echo "                         Installation Complete!"
echo "================================================================================"
echo ""
echo "Next steps:"
echo "  1. If you see warnings about ~/bin, run: source ~/.bashrc"
echo "  2. Verify tools: cd-hit-est -h && USalign -h"
echo "  3. Read README.md for usage instructions"
echo "  4. Run the pipeline: ./run_clustering.sh <pdb_dir> <output_dir>"
echo ""
echo "================================================================================"
