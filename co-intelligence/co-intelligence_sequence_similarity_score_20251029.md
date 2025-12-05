I'll walk you through the structure and functionality of the RNA sequence scoring code I created:

```python
import numpy as np
from typing import Tuple, Dict, List
import editdistance  # pip install editdistance
from Bio import SeqIO  # pip install biopython
import argparse

class RNASequenceScorer:
    """
    A class to score predicted RNA sequences against reference PDB structures.
    
    Combines base pairing score and sequence edit distance score for comprehensive evaluation.
    """
    
    def __init__(self):
        # Base pair scoring matrix as specified
        self.bp_score_matrix = {
            ('AU', 'AU'): 3, ('AU', 'UA'): 2,
            ('GC', 'GC'): 3, ('GC', 'CG'): 2,
            ('GU', 'GU'): 3, ('UG', 'GU'): 2,
            ('GC', 'GU'): 1, ('GC', 'AU'): 1,
            ('AU', 'GC'): 1, ('AU', 'GU'): 1,
            ('GU', 'GC'): 1, ('GU', 'AU'): 1
        }
        
        # Valid RNA base pairs
        self.valid_pairs = {'AU', 'UA', 'GC', 'CG', 'GU', 'UG'}
    
    def load_fasta_sequence(self, fasta_file: str) -> str:
        """
        Load predicted RNA sequence from FASTA file.
        
        Args:
            fasta_file: Path to FASTA file containing predicted sequence
            
        Returns:
            RNA sequence as string
        """
        try:
            record = next(SeqIO.parse(fasta_file, "fasta"))
            return str(record.seq).upper()
        except Exception as e:
            raise ValueError(f"Error reading FASTA file: {e}")
    
    def create_base_pair_matrix(self, sequence: str, reference_matrix: np.ndarray) -> np.ndarray:
        """
        Create base pair matrix for predicted sequence based on reference structure.
        
        Args:
            sequence: Predicted RNA sequence
            reference_matrix: NxN binary matrix indicating base pairs in reference structure
            
        Returns:
            NxN matrix with base pair information for predicted sequence
        """
        n = len(sequence)
        if reference_matrix.shape != (n, n):
            raise ValueError(f"Reference matrix shape {reference_matrix.shape} doesn't match sequence length {n}")
        
        bp_matrix = np.zeros((n, n), dtype='U2')  # Store base pair strings
        
        for i in range(n):
            for j in range(n):
                if reference_matrix[i, j] == 1:
                    pair = sequence[i] + sequence[j]
                    if pair in self.valid_pairs:
                        bp_matrix[i, j] = pair
                    else:
                        bp_matrix[i, j] = 'XX'  # Invalid pair
        
        return bp_matrix
    
    def calculate_base_pair_score(self, predicted_seq: str, reference_seq: str, 
                                 reference_matrix: np.ndarray) -> Tuple[float, Dict]:
        """
        Calculate base pairing score between predicted and reference sequences.
        
        Args:
            predicted_seq: Predicted RNA sequence
            reference_seq: Reference RNA sequence from PDB
            reference_matrix: NxN binary matrix indicating base pairs in reference structure
            
        Returns:
            Tuple of (normalized_score, detailed_stats)
        """
        if len(predicted_seq) != len(reference_seq):
            raise ValueError("Predicted and reference sequences must have the same length")
        
        n = len(predicted_seq)
        predicted_bp_matrix = self.create_base_pair_matrix(predicted_seq, reference_matrix)
        reference_bp_matrix = self.create_base_pair_matrix(reference_seq, reference_matrix)
        
        total_score = 0
        total_pairs = 0
        pair_counts = {}
        
        for i in range(n):
            for j in range(i + 1, n):  # Only consider upper triangle to avoid double counting
                if reference_matrix[i, j] == 1:
                    total_pairs += 1
                    pred_pair = predicted_bp_matrix[i, j]
                    ref_pair = reference_bp_matrix[i, j]
                    
                    # Get score for this pair comparison
                    score = self.bp_score_matrix.get((pred_pair, ref_pair), 0)
                    total_score += score
                    
                    # Track pair statistics
                    pair_key = f"{pred_pair}->{ref_pair}"
                    pair_counts[pair_key] = pair_counts.get(pair_key, 0) + 1
        
        # Normalize score (max possible score is 3 * total_pairs)
        max_possible_score = 3 * total_pairs
        normalized_score = total_score / max_possible_score if max_possible_score > 0 else 0
        
        stats = {
            'raw_score': total_score,
            'max_possible_score': max_possible_score,
            'normalized_score': normalized_score,
            'total_base_pairs': total_pairs,
            'pair_counts': pair_counts
        }
        
        return normalized_score, stats
    
    def calculate_edit_distance_score(self, predicted_seq: str, reference_seq: str) -> Tuple[float, Dict]:
        """
        Calculate normalized edit distance score between sequences.
        
        Args:
            predicted_seq: Predicted RNA sequence
            reference_seq: Reference RNA sequence
            
        Returns:
            Tuple of (normalized_score, detailed_stats)
        """
        if len(predicted_seq) != len(reference_seq):
            raise ValueError("Predicted and reference sequences must have the same length")
        
        edit_dist = editdistance.eval(predicted_seq, reference_seq)
        max_length = max(len(predicted_seq), len(reference_seq))
        
        # Convert edit distance to similarity score (1 - normalized_edit_distance)
        normalized_score = 1 - (edit_dist / max_length) if max_length > 0 else 1
        
        stats = {
            'edit_distance': edit_dist,
            'sequence_length': len(predicted_seq),
            'normalized_score': normalized_score,
            'accuracy': (max_length - edit_dist) / max_length if max_length > 0 else 1
        }
        
        return normalized_score, stats
    
    def calculate_combined_score(self, predicted_seq: str, reference_seq: str, 
                               reference_matrix: np.ndarray, bp_weight: float = 0.7, 
                               ed_weight: float = 0.3) -> Dict:
        """
        Calculate combined score using both base pairing and edit distance scores.
        
        Args:
            predicted_seq: Predicted RNA sequence
            reference_seq: Reference RNA sequence from PDB
            reference_matrix: NxN binary matrix indicating base pairs in reference structure
            bp_weight: Weight for base pairing score (default: 0.7)
            ed_weight: Weight for edit distance score (default: 0.3)
            
        Returns:
            Dictionary containing all scores and statistics
        """
        if abs(bp_weight + ed_weight - 1.0) > 1e-6:
            raise ValueError("Weights must sum to 1.0")
        
        # Calculate individual scores
        bp_score, bp_stats = self.calculate_base_pair_score(predicted_seq, reference_seq, reference_matrix)
        ed_score, ed_stats = self.calculate_edit_distance_score(predicted_seq, reference_seq)
        
        # Calculate combined score
        combined_score = bp_weight * bp_score + ed_weight * ed_score
        
        return {
            'combined_score': combined_score,
            'base_pair_score': bp_score,
            'edit_distance_score': ed_score,
            'weights': {'bp_weight': bp_weight, 'ed_weight': ed_weight},
            'base_pair_stats': bp_stats,
            'edit_distance_stats': ed_stats
        }
    
    def print_detailed_results(self, results: Dict, predicted_seq: str, reference_seq: str):
        """
        Print detailed scoring results in a formatted way.
        
        Args:
            results: Results dictionary from calculate_combined_score
            predicted_seq: Predicted RNA sequence
            reference_seq: Reference RNA sequence
        """
        print("="*60)
        print("RNA SEQUENCE SCORING RESULTS")
        print("="*60)
        
        print(f"Sequence Length: {len(predicted_seq)}")
        print(f"Reference Seq:  {reference_seq}")
        print(f"Predicted Seq:  {predicted_seq}")
        print()
        
        print(f"COMBINED SCORE: {results['combined_score']:.4f}")
        print(f"  Base Pair Score:    {results['base_pair_score']:.4f} (weight: {results['weights']['bp_weight']})")
        print(f"  Edit Distance Score: {results['edit_distance_score']:.4f} (weight: {results['weights']['ed_weight']})")
        print()
        
        # Base pair statistics
        bp_stats = results['base_pair_stats']
        print("BASE PAIR ANALYSIS:")
        print(f"  Total base pairs: {bp_stats['total_base_pairs']}")
        print(f"  Raw score: {bp_stats['raw_score']}/{bp_stats['max_possible_score']}")
        print("  Pair transitions:")
        for pair_transition, count in bp_stats['pair_counts'].items():
            print(f"    {pair_transition}: {count}")
        print()
        
        # Edit distance statistics
        ed_stats = results['edit_distance_stats']
        print("SEQUENCE ANALYSIS:")
        print(f"  Edit distance: {ed_stats['edit_distance']}")
        print(f"  Sequence accuracy: {ed_stats['accuracy']:.4f}")
        print("="*60)


def load_reference_matrix(matrix_file: str) -> np.ndarray:
    """
    Load reference base pair matrix from file.
    Expected format: CSV or NPY file with binary matrix
    
    Args:
        matrix_file: Path to matrix file
        
    Returns:
        NxN numpy array with binary base pair information
    """
    try:
        if matrix_file.endswith('.npy'):
            return np.load(matrix_file)
        elif matrix_file.endswith('.csv'):
            return np.loadtxt(matrix_file, delimiter=',', dtype=int)
        else:
            # Try to load as text file
            return np.loadtxt(matrix_file, dtype=int)
    except Exception as e:
        raise ValueError(f"Error loading reference matrix: {e}")


def main():
    """
    Main function to run RNA sequence scoring from command line.
    """
    parser = argparse.ArgumentParser(description='Score predicted RNA sequences against PDB reference')
    parser.add_argument('predicted_fasta', help='Path to predicted sequence FASTA file')
    parser.add_argument('reference_sequence', help='Reference RNA sequence (string or FASTA file)')
    parser.add_argument('reference_matrix', help='Path to reference base pair matrix file')
    parser.add_argument('--bp_weight', type=float, default=0.7, help='Weight for base pair score (default: 0.7)')
    parser.add_argument('--ed_weight', type=float, default=0.3, help='Weight for edit distance score (default: 0.3)')
    
    args = parser.parse_args()
    
    # Initialize scorer
    scorer = RNASequenceScorer()
    
    try:
        # Load predicted sequence
        predicted_seq = scorer.load_fasta_sequence(args.predicted_fasta)
        
        # Load reference sequence (could be string or file)
        if args.reference_sequence.endswith('.fasta') or args.reference_sequence.endswith('.fa'):
            reference_seq = scorer.load_fasta_sequence(args.reference_sequence)
        else:
            reference_seq = args.reference_sequence.upper()
        
        # Load reference matrix
        reference_matrix = load_reference_matrix(args.reference_matrix)
        
        # Calculate scores
        results = scorer.calculate_combined_score(
            predicted_seq, reference_seq, reference_matrix,
            args.bp_weight, args.ed_weight
        )
        
        # Print results
        scorer.print_detailed_results(results, predicted_seq, reference_seq)
        
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    # Example usage
    # Uncomment and modify the following lines for direct script usage
    
    # scorer = RNASequenceScorer()
    # 
    # # Example data
    # predicted_seq = "AUCGAUCGAUC"
    # reference_seq = "AUCGAUCGAUC"  
    # reference_matrix = np.array([
    #     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
    #     [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
    #     [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
    #     [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
    #     [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
    #     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    #     [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
    #     [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
    #     [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    #     [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    #     [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    # ])
    # 
    # results = scorer.calculate_combined_score(predicted_seq, reference_seq, reference_matrix)
    # scorer.print_detailed_results(results, predicted_seq, reference_seq)
    
    main()
```

## **Overall Architecture**

The code is organized around a main class `RNASequenceScorer` that encapsulates all scoring functionality, plus helper functions for file I/O and a command-line interface.

## **Core Class: RNASequenceScorer**

### **Initialization (`__init__`)**
```python
self.bp_score_matrix = {
    ('AU', 'AU'): 3, ('AU', 'UA'): 2,
    # ... more pairs
}
self.valid_pairs = {'AU', 'UA', 'GC', 'CG', 'GU', 'UG'}
```
- Sets up the scoring matrix exactly as you specified
- Defines valid RNA base pairs for validation

### **Key Methods:**

#### **1. File Loading (`load_fasta_sequence`)**
- Uses BioPython to parse FASTA files
- Converts sequences to uppercase for consistency
- Handles file reading errors gracefully

#### **2. Base Pair Matrix Creation (`create_base_pair_matrix`)**
```python
def create_base_pair_matrix(self, sequence: str, reference_matrix: np.ndarray)
```
- Takes a sequence and the reference structure matrix
- Creates an NxN matrix where each cell contains the base pair string (like "AU", "GC")
- Only fills positions where reference_matrix[i,j] = 1
- Marks invalid pairs as "XX"

#### **3. Base Pair Scoring (`calculate_base_pair_score`)**
This is the core of your `score_bp` requirement:

```python
for i in range(n):
    for j in range(i + 1, n):  # Upper triangle only
        if reference_matrix[i, j] == 1:  # Only score actual base pairs
            pred_pair = predicted_bp_matrix[i, j]
            ref_pair = reference_bp_matrix[i, j]
            score = self.bp_score_matrix.get((pred_pair, ref_pair), 0)
```

**How it works:**
- Iterates through the upper triangle of the matrix (avoids double-counting)
- Only scores positions where the reference indicates a base pair
- Looks up the score from your provided scoring table
- Normalizes the final score by dividing by maximum possible score (3 × number of base pairs)

#### **4. Edit Distance Scoring (`calculate_edit_distance_score`)**
This implements your `score_ed` requirement:
- Uses the `editdistance` library for accurate Levenshtein distance
- Converts distance to similarity: `1 - (edit_distance / sequence_length)`
- Returns normalized score between 0-1

#### **5. Combined Scoring (`calculate_combined_score`)**
- Combines both scores with customizable weights
- Default: 70% base pairing, 30% edit distance
- Validates that weights sum to 1.0
- Returns comprehensive statistics

## **Scoring Logic Details**

### **Base Pair Scoring Process:**
1. Create matrices for both predicted and reference sequences
2. For each position where reference_matrix[i,j] = 1:
   - Extract base pair from predicted sequence (e.g., "AU")
   - Extract base pair from reference sequence (e.g., "GC")  
   - Look up score in your table (e.g., AU→GC = 1 point)
3. Sum all scores and normalize

### **Normalization Strategy:**
- **Base pair score:** Divided by (3 × total_base_pairs) since 3 is maximum possible per pair
- **Edit distance score:** Converted to similarity score (1 - normalized_distance)
- **Final score:** Weighted combination, still 0-1 scale

## **Input/Output Design**

### **Inputs the code handles:**
- **FASTA files** for sequences (using BioPython)
- **Matrix files** in multiple formats (.npy, .csv, or plain text)
- **Command line arguments** for batch processing

### **Outputs provided:**
- **Numerical scores** (combined, base pair, edit distance)
- **Detailed statistics** (pair transition counts, accuracy metrics)
- **Formatted reports** for human readability

## **Error Handling & Validation**

The code includes several safety checks:
- Sequence length validation
- Matrix dimension validation  
- File format validation
- Weight validation (must sum to 1.0)
- Invalid base pair handling

## **Key Design Decisions**

### **Why upper triangle iteration?**
```python
for j in range(i + 1, n):  # Only upper triangle
```
Base pair matrices are symmetric, so this avoids double-counting the same pair.

### **Why normalize scores?**
Makes scores comparable across different:
- Sequence lengths
- Numbers of base pairs
- Different RNA structures

### **Why separate the scoring components?**
Allows you to:
- Analyze each component independently
- Adjust weights based on your priorities
- Debug scoring issues more easily

## **Extensibility Features**

The code is designed to be easily modified:
- **Scoring matrix:** Easy to update the `bp_score_matrix` dictionary
- **Weights:** Adjustable via parameters
- **File formats:** Additional formats can be added to loading functions
- **Additional metrics:** New scoring methods can be added as class methods

## **Performance Considerations**

- Uses NumPy for matrix operations (efficient)
- Only processes upper triangle of matrices
- Minimal memory footprint for sequence operations
- Suitable for typical RNA sequence lengths (hundreds to thousands of nucleotides)

The code balances functionality, readability, and performance while strictly implementing your specified requirements for both base pair scoring and edit distance calculation.