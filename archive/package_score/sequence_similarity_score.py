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
            # Perfect matches (3 points)
            ('AU', 'AU'): 3, ('UA', 'UA'): 3,
            ('GC', 'GC'): 3, ('CG', 'CG'): 3,
            ('GU', 'GU'): 3, ('UG', 'UG'): 3,
            # Reverse pairs (2 points)
            ('AU', 'UA'): 2, ('UA', 'AU'): 2,
            ('GC', 'CG'): 2, ('CG', 'GC'): 2,
            ('GU', 'UG'): 2, ('UG', 'GU'): 2,
            # Cross matches (1 point) - all other base pair to base pair combinations
            ('AU', 'GC'): 1, ('AU', 'CG'): 1, ('AU', 'GU'): 1, ('AU', 'UG'): 1,
            ('UA', 'GC'): 1, ('UA', 'CG'): 1, ('UA', 'GU'): 1, ('UA', 'UG'): 1,
            ('GC', 'AU'): 1, ('GC', 'UA'): 1, ('GC', 'GU'): 1, ('GC', 'UG'): 1,
            ('CG', 'AU'): 1, ('CG', 'UA'): 1, ('CG', 'GU'): 1, ('CG', 'UG'): 1,
            ('GU', 'AU'): 1, ('GU', 'UA'): 1, ('GU', 'GC'): 1, ('GU', 'CG'): 1,
            ('UG', 'AU'): 1, ('UG', 'UA'): 1, ('UG', 'GC'): 1, ('UG', 'CG'): 1
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
                               reference_matrix: np.ndarray, lambda_param: float = 0.7) -> Dict:
        """
        Calculate combined score using both base pairing and edit distance scores.
        
        Args:
            predicted_seq: Predicted RNA sequence
            reference_seq: Reference RNA sequence from PDB
            reference_matrix: NxN binary matrix indicating base pairs in reference structure
            lambda_param: Weight parameter (0-1) for base pairing score. 
                         bp_weight = lambda_param, ed_weight = 1 - lambda_param (default: 0.7)
            
        Returns:
            Dictionary containing all scores and statistics
        """
        if not (0 <= lambda_param <= 1):
            raise ValueError("lambda_param must be between 0 and 1")
        
        bp_weight = lambda_param
        ed_weight = 1 - lambda_param
        
        # Calculate individual scores
        bp_score, bp_stats = self.calculate_base_pair_score(predicted_seq, reference_seq, reference_matrix)
        ed_score, ed_stats = self.calculate_edit_distance_score(predicted_seq, reference_seq)
        
        # Calculate combined score
        combined_score = bp_weight * bp_score + ed_weight * ed_score
        
        return {
            'combined_score': combined_score,
            'base_pair_score': bp_score,
            'edit_distance_score': ed_score,
            'lambda_param': lambda_param,
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
        print(f"  Base Pair Score:    {results['base_pair_score']:.4f} (λ = {results['lambda_param']:.2f})")
        print(f"  Edit Distance Score: {results['edit_distance_score']:.4f} (1-λ = {1-results['lambda_param']:.2f})")
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
    parser = argparse.ArgumentParser(
        description='Score predicted RNA sequences against PDB reference',
        epilog="""Examples:
  %(prog)s pred.fasta ref.fasta matrix.mat
  %(prog)s pred.fasta AUCGAUCG matrix.mat --lambda_param 0.5
  %(prog)s pred.fasta ref.fasta matrix.mat --format json --output results.json

Supported matrix formats: .mat (MATLAB), .npy (NumPy), .txt/.csv (text)""",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('predicted_fasta', help='Path to predicted sequence FASTA file')
    parser.add_argument('reference_sequence', help='Reference RNA sequence (string or FASTA file)')
    parser.add_argument('reference_matrix', help='Path to reference base pair matrix file')
    
    # Optional arguments
    parser.add_argument('--lambda_param', '--weight', type=float, default=0.7, 
                       help='Lambda parameter (0-1) for score weighting. bp_weight=λ, ed_weight=1-λ (default: 0.7)')
    parser.add_argument('--output', '-o', help='Output file path (JSON format)')
    parser.add_argument('--quiet', '-q', action='store_true', help='Suppress verbose output')
    parser.add_argument('--format', choices=['detailed', 'summary', 'json', 'csv'], default='detailed',
                       help='Output format (default: detailed)')
    parser.add_argument('--validate-only', dest='validate_only', action='store_true',
                       help='Only validate inputs without computing scores')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
    args = parser.parse_args()
    
    # Validate lambda parameter
    if not 0 <= args.lambda_param <= 1:
        parser.error(f"Lambda parameter must be between 0 and 1, got: {args.lambda_param}")
    
    # Initialize scorer
    scorer = RNASequenceScorer()
    
    try:
        # Validate input files exist
        import os
        if not os.path.exists(args.predicted_fasta):
            parser.error(f"Predicted FASTA file not found: {args.predicted_fasta}")
        if not os.path.exists(args.reference_matrix):
            parser.error(f"Reference matrix file not found: {args.reference_matrix}")
        
        # Load predicted sequence
        predicted_seq = scorer.load_fasta_sequence(args.predicted_fasta)
        
        # Load reference sequence (could be string or file)
        if (args.reference_sequence.endswith('.fasta') or 
            args.reference_sequence.endswith('.fa') or 
            args.reference_sequence.endswith('.fas')):
            if not os.path.exists(args.reference_sequence):
                parser.error(f"Reference FASTA file not found: {args.reference_sequence}")
            reference_seq = scorer.load_fasta_sequence(args.reference_sequence)
        else:
            # Validate RNA sequence string
            reference_seq = args.reference_sequence.upper()
            valid_bases = set('AUCG')
            if not all(base in valid_bases for base in reference_seq):
                parser.error(f"Invalid RNA sequence. Only A, U, C, G allowed. Got: {args.reference_sequence}")
        
        # Load reference matrix
        reference_matrix = load_reference_matrix(args.reference_matrix)
        
        # Validate sequence lengths match matrix dimensions
        if len(predicted_seq) != reference_matrix.shape[0] or len(reference_seq) != reference_matrix.shape[0]:
            parser.error(f"Sequence length mismatch. Predicted: {len(predicted_seq)}, "
                        f"Reference: {len(reference_seq)}, Matrix: {reference_matrix.shape[0]}x{reference_matrix.shape[1]}")
        
        if args.validate_only:
            if not args.quiet:
                print("✓ All inputs validated successfully")
                print(f"  Predicted sequence length: {len(predicted_seq)}")
                print(f"  Reference sequence length: {len(reference_seq)}")
                print(f"  Matrix dimensions: {reference_matrix.shape}")
            return 0
        
        # Calculate scores
        results = scorer.calculate_combined_score(
            predicted_seq, reference_seq, reference_matrix,
            args.lambda_param
        )
        
        # Handle output based on format and options
        if args.format == 'json' or args.output:
            import json
            output_data = {
                'sequences': {
                    'predicted': predicted_seq,
                    'reference': reference_seq
                },
                'scores': results,
                'parameters': {
                    'lambda_param': args.lambda_param,
                    'matrix_file': args.reference_matrix
                }
            }
            
            if args.output:
                with open(args.output, 'w') as f:
                    json.dump(output_data, f, indent=2)
                if not args.quiet:
                    print(f"Results saved to {args.output}")
            else:
                print(json.dumps(output_data, indent=2))
                
        elif args.format == 'csv':
            print("metric,score")
            print(f"combined_score,{results['combined_score']:.4f}")
            print(f"base_pair_score,{results['base_pair_score']:.4f}")
            print(f"edit_distance_score,{results['edit_distance_score']:.4f}")
            
        elif args.format == 'summary' or args.quiet:
            print(f"Combined Score: {results['combined_score']:.4f}")
            print(f"Base Pair Score: {results['base_pair_score']:.4f}")
            print(f"Edit Distance Score: {results['edit_distance_score']:.4f}")
            
        else:  # detailed format
            scorer.print_detailed_results(results, predicted_seq, reference_seq)
        
    except KeyboardInterrupt:
        print("\nOperation cancelled by user")
        return 1
    except Exception as e:
        if args.quiet:
            print(f"Error: {e}")
        else:
            print(f"Error: {e}")
            print("\nUse --help for usage information")
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
    # results = scorer.calculate_combined_score(predicted_seq, reference_seq, reference_matrix, lambda_param=0.8)
    # scorer.print_detailed_results(results, predicted_seq, reference_seq)
    
    main()