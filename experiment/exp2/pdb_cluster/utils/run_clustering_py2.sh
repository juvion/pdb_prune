#!/bin/bash
###############################################################################
# RNA PDB Clustering Pipeline - Complete Workflow Script
# 
# This script automates the entire process of clustering RNA structures
# based on both sequence similarity and structure similarity.
#
# Usage: ./run_clustering.sh <pdb_directory> <output_directory>
###############################################################################

set -e  # Exit on error

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_step() {
    echo -e "${BLUE}[STEP $1]${NC} $2"
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

# Check arguments
if [ $# -lt 2 ]; then
    echo "Usage: $0 <pdb_directory> <output_directory> [seq_thresholds] [struct_thresholds]"
    echo ""
    echo "Arguments:"
    echo "  pdb_directory     : Directory containing PDB files"
    echo "  output_directory  : Directory to save results"
    echo "  seq_thresholds    : Sequence similarity thresholds (default: 0.8,0.6,0.4)"
    echo "  struct_thresholds : Structure similarity thresholds (default: 0.6,0.5,0.4)"
    echo ""
    echo "Example:"
    echo "  $0 ./pdb_files ./results"
    echo "  $0 ./pdb_files ./results 0.8,0.6,0.4 0.6,0.5,0.4"
    exit 1
fi

PDB_DIR=$1
OUTPUT_DIR=$2
SEQ_THRESHOLDS=${3:-"0.8,0.6,0.4"}
STRUCT_THRESHOLDS=${4:-"0.6,0.5,0.4"}

# Absolute project home for clustering utilities
# Note: This replaces previous relative path references (../*.py)
PROJECT_HOME='/Users/xiaojuzhang/Dev/pdb_prune/experiment/exp2/pdb_cluster'

# Convert comma-separated thresholds to arrays
IFS=',' read -ra SEQ_THRESH_ARRAY <<< "$SEQ_THRESHOLDS"
IFS=',' read -ra STRUCT_THRESH_ARRAY <<< "$STRUCT_THRESHOLDS"

echo "================================================================================"
echo "                    RNA PDB CLUSTERING PIPELINE"
echo "================================================================================"
echo "PDB Directory:            $PDB_DIR"
echo "Output Directory:         $OUTPUT_DIR"
echo "Sequence Thresholds:      ${SEQ_THRESH_ARRAY[@]}"
echo "Structure Thresholds:     ${STRUCT_THRESH_ARRAY[@]}"
echo "================================================================================"
echo ""

# Check if input directory exists
if [ ! -d "$PDB_DIR" ]; then
    print_error "PDB directory does not exist: $PDB_DIR"
    exit 1
fi

# Check if required tools are available
print_step "0" "Checking required tools..."

check_tool() {
    if command -v $1 &> /dev/null; then
        print_success "$1 is installed"
    else
        print_error "$1 is not installed. Please install it first."
        exit 1
    fi
}

check_tool python
check_tool cd-hit-est
check_tool USalign

# Check Python packages
python -c "import Bio" 2>/dev/null || { print_error "BioPython not installed. Run: pip install biopython"; exit 1; }
python -c "import numpy" 2>/dev/null || { print_error "NumPy not installed. Run: pip install numpy"; exit 1; }
python -c "import pandas" 2>/dev/null || { print_error "Pandas not installed. Run: pip install pandas"; exit 1; }
python -c "import scipy" 2>/dev/null || { print_error "SciPy not installed. Run: pip install scipy"; exit 1; }
python -c "import matplotlib" 2>/dev/null || { print_error "Matplotlib not installed. Run: pip install matplotlib"; exit 1; }
python -c "import tqdm" 2>/dev/null || { print_error "tqdm not installed. Run: pip install tqdm"; exit 1; }

print_success "All required tools and packages are available"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# Step 1: Extract RNA chains
print_step "1" "Extracting RNA chains from PDB files..."
python "$PROJECT_HOME/utils/extract_rna_chains.py" --pdb_dir "$PDB_DIR" --output_dir ./rna_chains
if [ $? -eq 0 ]; then
    print_success "RNA chains extracted successfully"
else
    print_error "Failed to extract RNA chains"
    exit 1
fi
echo ""

# Step 2: Extract sequences to FASTA
print_step "2" "Extracting sequences to FASTA format..."
python "$PROJECT_HOME/utils/extract_sequences.py" --pdb_dir ./rna_chains --output rna_sequences.fasta
if [ $? -eq 0 ]; then
    print_success "Sequences extracted successfully"
else
    print_error "Failed to extract sequences"
    exit 1
fi
echo ""

# Step 3: Sequence clustering with CD-HIT
print_step "3" "Performing sequence clustering with CD-HIT..."
for threshold in "${SEQ_THRESH_ARRAY[@]}"; do
    echo "  Clustering at threshold: $threshold"
    # CD-HIT-EST requires identity cutoff (-c) to be >= 0.8
    if (( $(echo "$threshold < 0.8" | bc -l) )); then
        print_warning "CD-HIT-EST requires -c >= 0.8. Skipping threshold $threshold"
        continue
    fi
    
    # Determine word size based on threshold (CD-HIT-EST allows 4-10; using 4 for sensitivity)
    WORD_SIZE=4
    
    cd-hit-est \
        -i rna_sequences.fasta \
        -o "clusters_seq_${threshold}.fasta" \
        -c "$threshold" \
        -n "$WORD_SIZE" \
        -M 2000 \
        -T 4 \
        > "cdhit_${threshold}.log" 2>&1
    
    if [ $? -eq 0 ]; then
        # Parse clusters
        python "$PROJECT_HOME/utils/parse_cdhit_clusters.py" \
            --input "clusters_seq_${threshold}.fasta.clstr" \
            --output "sequence_clusters_${threshold}.csv" \
            --dist "sequence_clusters_${threshold}_distribution.csv"
        print_success "Sequence clustering at threshold $threshold completed"
    else
        print_error "Failed to cluster sequences at threshold $threshold"
    fi
done
echo ""

# Step 4: Calculate TM-score matrix
print_step "4" "Calculating TM-score matrix (this may take a while)..."
python "$PROJECT_HOME/utils/calculate_tm_matrix.py" \
    --pdb_dir ./rna_chains \
    --output tm_score_matrix.csv \
    --checkpoint tm_score_checkpoint.csv

if [ $? -eq 0 ]; then
    print_success "TM-score matrix calculated successfully"
else
    print_error "Failed to calculate TM-score matrix"
    exit 1
fi
echo ""

# Step 5: Structure clustering
print_step "5" "Performing structure clustering..."
python "$PROJECT_HOME/utils/agglomerative_clustering.py" \
    --input tm_score_matrix.csv \
    --thresholds ${STRUCT_THRESH_ARRAY[@]} \
    --method average \
    --output_prefix structure_clusters \
    --plot

if [ $? -eq 0 ]; then
    print_success "Structure clustering completed successfully"
else
    print_error "Failed to perform structure clustering"
    exit 1
fi
echo ""

# Step 6: Generate summary report
print_step "6" "Generating summary report..."

REPORT_FILE="clustering_summary_report.txt"

cat > "$REPORT_FILE" << EOF
===============================================================================
                    RNA PDB CLUSTERING SUMMARY REPORT
===============================================================================
Generated on: $(date)
PDB Directory: $PDB_DIR
Output Directory: $OUTPUT_DIR

===============================================================================
SEQUENCE CLUSTERING RESULTS (CD-HIT)
===============================================================================

EOF

for threshold in "${SEQ_THRESH_ARRAY[@]}"; do
    if [ -f "sequence_clusters_${threshold}.csv" ]; then
        NUM_SEQS=$(wc -l < "sequence_clusters_${threshold}.csv")
        NUM_SEQS=$((NUM_SEQS - 1))  # Subtract header
        NUM_CLUSTERS=$(tail -n +2 "sequence_clusters_${threshold}.csv" | cut -d',' -f2 | sort -u | wc -l)
        
        echo "Threshold: $threshold" >> "$REPORT_FILE"
        echo "  Total sequences: $NUM_SEQS" >> "$REPORT_FILE"
        echo "  Total clusters: $NUM_CLUSTERS" >> "$REPORT_FILE"
        echo "" >> "$REPORT_FILE"
    fi
done

cat >> "$REPORT_FILE" << EOF

===============================================================================
STRUCTURE CLUSTERING RESULTS (TM-score + Agglomerative)
===============================================================================

EOF

for threshold in "${STRUCT_THRESH_ARRAY[@]}"; do
    if [ -f "structure_clusters_${threshold}.csv" ]; then
        NUM_STRUCTS=$(wc -l < "structure_clusters_${threshold}.csv")
        NUM_STRUCTS=$((NUM_STRUCTS - 1))  # Subtract header
        NUM_CLUSTERS=$(tail -n +2 "structure_clusters_${threshold}.csv" | cut -d',' -f2 | sort -u | wc -l)
        
        echo "Threshold: $threshold" >> "$REPORT_FILE"
        echo "  Total structures: $NUM_STRUCTS" >> "$REPORT_FILE"
        echo "  Total clusters: $NUM_CLUSTERS" >> "$REPORT_FILE"
        echo "" >> "$REPORT_FILE"
    fi
done

cat >> "$REPORT_FILE" << EOF

===============================================================================
OUTPUT FILES
===============================================================================

Sequence Clustering:
EOF

for threshold in "${SEQ_THRESH_ARRAY[@]}"; do
    echo "  - sequence_clusters_${threshold}.csv" >> "$REPORT_FILE"
    echo "  - sequence_clusters_${threshold}_distribution.csv" >> "$REPORT_FILE"
done

cat >> "$REPORT_FILE" << EOF

Structure Clustering:
  - tm_score_matrix.csv
EOF

for threshold in "${STRUCT_THRESH_ARRAY[@]}"; do
    echo "  - structure_clusters_${threshold}.csv" >> "$REPORT_FILE"
    echo "  - structure_clusters_${threshold}_distribution.csv" >> "$REPORT_FILE"
done

cat >> "$REPORT_FILE" << EOF
  - structure_clusters_dendrogram.png
  - structure_clusters_size_distribution.png
  - structure_clusters_comparison.csv

Metadata:
  - rna_chains/rna_chains_metadata.json
  - rna_sequences.fasta

===============================================================================
EOF

print_success "Summary report generated: $REPORT_FILE"
echo ""

# Display summary
echo "================================================================================"
echo "                         PIPELINE COMPLETED!"
echo "================================================================================"
echo ""
echo "Summary:"
cat "$REPORT_FILE" | grep -A 100 "SEQUENCE CLUSTERING RESULTS" | head -n 50
echo ""
echo "Full report saved to: $OUTPUT_DIR/$REPORT_FILE"
echo "================================================================================"
