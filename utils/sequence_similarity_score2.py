"""
RNA SEQUENCE SCORING MECHANISM - OPTIMIZED FOR MODEL TRAINING
==============================================================

PHILOSOPHY:
-----------
This scoring system evaluates predicted RNA sequences against reference structures
from PDB files. It combines two complementary metrics:

1. STRUCTURAL SCORE (Base Pairing Score):
   - Evaluates how well the prediction captures both PAIRED and UNPAIRED regions
   - Solves the "zero-pair problem": sequences without structure get meaningful scores
   - Rewards correct identification of non-structured regions (not just base pairs)

2. SEQUENCE SCORE (Edit Distance Score):
   - Measures sequence-level similarity regardless of structure
   - Provides complementary information about prediction quality

SCORING STRATEGY:
-----------------
Upper Triangle Scoring (avoids double-counting):
  - Only scores positions where i < j (N×(N-1)/2 positions total)
  - Excludes diagonal (self-pairing impossible)

For each position (i,j) in the upper triangle:

  IF REFERENCE HAS A PAIR at (i,j):
    - Correctly predicted pair with perfect match (AU→AU, GC→GC, GU→GU): 4 points
    - Correctly predicted pair with good match (AU→UA, GC→CG, UG→GU): 3 points  
    - Correctly predicted pair with weak match (GC→AU, GC→GU, etc.): 2 points
    - Invalid pair match or no pair predicted: 0 points
    
  IF REFERENCE HAS NO PAIR at (i,j):
    - When no predicted structure is provided (sequence-only scoring), all
      unpaired positions in the reference are treated as correct: 1 point
      each. This avoids penalizing positions that happen to be complementary
      in sequence but are not paired in the reference structure.

Normalization:
  - Max Score = 4 × (paired positions) + 1 × (unpaired positions)
  - Normalized Score = Raw Score / Max Score  (range: 0 to 1)

COMBINED SCORING:
-----------------
  Final Score = λ × Structure_Score + (1-λ) × Sequence_Score
  
  Where λ ∈ [0,1]:
    - λ = 1.0: Pure structure scoring (only base pairing matters)
    - λ = 0.7: Structure-emphasized (default, 70% structure, 30% sequence)
    - λ = 0.5: Balanced scoring
    - λ = 0.3: Sequence-emphasized  
    - λ = 0.0: Pure sequence scoring (only edit distance matters)

EDGE CASE HANDLING:
-------------------
  - If no base pairs exist (e.g., short sequences, poly-A/U/G/C), 
    automatically falls back to edit distance only (λ → 0)

PERFORMANCE OPTIMIZATIONS:
--------------------------
  - Numba JIT compilation for 10-100x speedup (optional, automatic if available)
  - Direct string/array input (no file I/O in scoring functions)
  - Pre-computed lookup tables for O(1) scoring
  - Vectorized operations where possible
"""

import numpy as np
from typing import Tuple, Dict, List
import editdistance  # pip install editdistance
from Bio import SeqIO  # pip install biopython
import argparse

# Try to import numba for JIT compilation (optional but recommended)
try:
    from numba import jit
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    def jit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator


class RNASequenceScorer:
    """
    A class to score predicted RNA sequences against reference PDB structures.
    
    Optimized for model training with:
    - Numba JIT compilation (optional, 10-100x speedup)
    - Direct string/array input (no file I/O in hot path)
    - Pre-computed lookup tables
    """
    
    def __init__(self, use_numba: bool = True, min_separation: int = 4, unpaired_weight: float = 0.0):
        """
        Initialize the RNA sequence scorer.
        
        Args:
            use_numba: Enable Numba JIT compilation if available (default: True)
                      First call will be slow (compilation), subsequent calls very fast
            min_separation: Minimum sequence separation |j-i| for base pairing (default: 4)
                           Positions with j-i < min_separation are excluded from scoring
                           to avoid bias from physically impossible pairs
        """
        self.use_numba = use_numba and NUMBA_AVAILABLE
        self.min_separation = min_separation
        # Weight applied to unpaired component within base-pair score. Default 0.0 (paired-only).
        self.unpaired_weight = max(0.0, min(1.0, float(unpaired_weight)))
        
        if self.use_numba:
            print("Numba JIT enabled - first scoring call will compile (slow), then fast")
        
        print(f"Minimum separation constraint: {self.min_separation} (positions with j-i < {self.min_separation} excluded from scoring)")
        if self.unpaired_weight > 0:
            print(f"Unpaired weight within base-pair score: {self.unpaired_weight:.2f} (paired weight {(1-self.unpaired_weight):.2f})")
        
        # Valid RNA base pairs
        self.valid_pairs = {'AU', 'UA', 'GC', 'CG', 'GU', 'UG'}
        
        # Base pair scoring matrix - complete with all valid pair combinations
        self.bp_score_matrix = {
            # Perfect matches (score 3 -> 4 points after +1)
            ('AU', 'AU'): 3, ('UA', 'UA'): 3,
            ('GC', 'GC'): 3, ('CG', 'CG'): 3,
            ('GU', 'GU'): 3, ('UG', 'UG'): 3,
            
            # Good matches - orientation flips (score 2 -> 3 points)
            ('AU', 'UA'): 2, ('UA', 'AU'): 2,
            ('GC', 'CG'): 2, ('CG', 'GC'): 2,
            ('GU', 'UG'): 2, ('UG', 'GU'): 2,
            
            # Weak matches - cross-type pairings (score 1 -> 2 points)
            ('AU', 'GC'): 1, ('AU', 'CG'): 1, ('AU', 'GU'): 1, ('AU', 'UG'): 1,
            ('UA', 'GC'): 1, ('UA', 'CG'): 1, ('UA', 'GU'): 1, ('UA', 'UG'): 1,
            ('GC', 'AU'): 1, ('GC', 'UA'): 1, ('GC', 'GU'): 1, ('GC', 'UG'): 1,
            ('CG', 'AU'): 1, ('CG', 'UA'): 1, ('CG', 'GU'): 1, ('CG', 'UG'): 1,
            ('GU', 'AU'): 1, ('GU', 'UA'): 1, ('GU', 'GC'): 1, ('GU', 'CG'): 1,
            ('UG', 'AU'): 1, ('UG', 'UA'): 1, ('UG', 'GC'): 1, ('UG', 'CG'): 1,
        }
        
        # Create fast lookup structures for Numba
        self.base_to_idx = {'A': 0, 'U': 1, 'G': 2, 'C': 3}
        self._init_score_array()
        # Precompute valid pair lookup (4x4) and cache for upper-triangle indices
        self.valid_pair_array = np.zeros((4, 4), dtype=np.int8)
        for pair in self.valid_pairs:
            i_idx = self.base_to_idx[pair[0]]
            j_idx = self.base_to_idx[pair[1]]
            self.valid_pair_array[i_idx, j_idx] = 1
        # Cache for np.triu_indices by (n, min_separation) to avoid recomputation in training loops
        self._triu_cache: Dict[Tuple[int, int], Tuple[np.ndarray, np.ndarray]] = {}
        # Lightweight cache for encoded sequences to reduce repeated work in training
        self._seq_encode_cache: Dict[str, np.ndarray] = {}
        self._seq_encode_cache_max = 128
    
    def _init_score_array(self):
        """Initialize 4D scoring array for fast lookup: score_array[pred_i, pred_j, ref_i, ref_j]"""
        self.score_array = np.zeros((4, 4, 4, 4), dtype=np.int8)
        
        for (pred_pair, ref_pair), score in self.bp_score_matrix.items():
            pred_i_idx = self.base_to_idx[pred_pair[0]]
            pred_j_idx = self.base_to_idx[pred_pair[1]]
            ref_i_idx = self.base_to_idx[ref_pair[0]]
            ref_j_idx = self.base_to_idx[ref_pair[1]]
            self.score_array[pred_i_idx, pred_j_idx, ref_i_idx, ref_j_idx] = score

    def _get_triu_indices(self, n: int) -> Tuple[np.ndarray, np.ndarray]:
        """Cached upper-triangle indices with min_separation offset."""
        key = (n, self.min_separation)
        cached = self._triu_cache.get(key)
        if cached is not None:
            return cached
        tri_i, tri_j = np.triu_indices(n, k=self.min_separation)
        # Ensure contiguous arrays for faster advanced indexing
        tri_i = np.ascontiguousarray(tri_i, dtype=np.int32)
        tri_j = np.ascontiguousarray(tri_j, dtype=np.int32)
        self._triu_cache[key] = (tri_i, tri_j)
        return tri_i, tri_j

    def _encode_sequence(self, seq: str) -> np.ndarray:
        """Encode RNA sequence to index array with caching and uppercase handling."""
        # Normalize to uppercase to avoid invalid indices for lowercase input
        seq_u = seq.upper()
        cached = self._seq_encode_cache.get(seq_u)
        if cached is not None:
            return cached
        arr = np.array([self.base_to_idx.get(b, -1) for b in seq_u], dtype=np.int8)
        arr = np.ascontiguousarray(arr)
        # Maintain small cache to avoid memory bloat
        if len(self._seq_encode_cache) >= self._seq_encode_cache_max:
            self._seq_encode_cache.clear()
        self._seq_encode_cache[seq_u] = arr
        return arr

    def _score_loop_vectorized(self, pred_indices: np.ndarray, ref_indices: np.ndarray,
                               reference_matrix: np.ndarray) -> Tuple[int, int, int, int, int]:
        """
        Fully vectorized scoring path (fast even without Numba).
        Mirrors the semantics of the JIT loop, including invalid-base skipping.
        """
        n = pred_indices.shape[0]
        tri_i, tri_j = self._get_triu_indices(n)

        # Valid bases mask (skip pairs with invalid indices), matches JIT behavior
        valid_pred = (pred_indices[tri_i] >= 0) & (pred_indices[tri_j] >= 0)
        valid_ref = (ref_indices[tri_i] >= 0) & (ref_indices[tri_j] >= 0)
        valid_mask = valid_pred & valid_ref
        if not np.any(valid_mask):
            return 0, 0, 0, 0, 0

        tri_i_v = tri_i[valid_mask]
        tri_j_v = tri_j[valid_mask]

        # Reference pair mask on valid positions
        ref_has_pair = (reference_matrix[tri_i_v, tri_j_v] == 1)
        paired_mask = ref_has_pair
        unpaired_mask = ~ref_has_pair

        paired_positions = int(np.count_nonzero(paired_mask))
        unpaired_positions = int(np.count_nonzero(unpaired_mask))

        total_score = 0
        paired_correct = 0
        unpaired_correct = 0

        if paired_positions > 0:
            pi = pred_indices[tri_i_v]
            pj = pred_indices[tri_j_v]
            ri = ref_indices[tri_i_v]
            rj = ref_indices[tri_j_v]
            # Select only paired positions
            pi_p = pi[paired_mask]
            pj_p = pj[paired_mask]
            ri_p = ri[paired_mask]
            rj_p = rj[paired_mask]
            base_scores = self.score_array[pi_p, pj_p, ri_p, rj_p]
            pos_mask = base_scores > 0
            paired_correct = int(np.count_nonzero(pos_mask))
            if paired_correct > 0:
                total_score += int(np.sum((base_scores[pos_mask] + 1).astype(np.int64)))

        if unpaired_positions > 0:
            # Sequence-only structural scoring:
            # Treat all unpaired positions in the reference as correct (1 point each),
            # since we do not have a predicted structure matrix to assert false positives.
            unpaired_correct = unpaired_positions
            total_score += unpaired_correct

        return total_score, paired_positions, unpaired_positions, paired_correct, unpaired_correct
    
    def load_fasta_sequence(self, fasta_file: str) -> str:
        """
        Load RNA sequence from FASTA file.
        
        Note: Avoid using this in training loops - use direct string input instead.
        
        Args:
            fasta_file: Path to FASTA file containing sequence
            
        Returns:
            RNA sequence as string
        """
        try:
            record = next(SeqIO.parse(fasta_file, "fasta"))
            return str(record.seq).upper()
        except Exception as e:
            raise ValueError(f"Error reading FASTA file: {e}")
    
    def calculate_base_pair_score(self, predicted_seq: str, reference_seq: str, 
                                 reference_matrix: np.ndarray) -> Tuple[float, Dict]:
        """
        Calculate base pairing score between predicted and reference sequences.
        
        OPTIMIZED for training: Uses JIT compilation if available, no file I/O.
        
        CRITICAL: Only scores positions where j-i >= min_separation to avoid bias
        from physically impossible pairs.
        
        Scoring strategy:
        - Excludes positions where j-i < min_separation (cannot physically pair)
        - Scores only upper triangle (i < j) to avoid double-counting
        - Unpaired → Unpaired: 1 point
        - Paired → Paired: 2-4 points based on match quality
        - Structure mismatch: 0 points
        - Normalizes by theoretical maximum score
        
        Args:
            predicted_seq: Predicted RNA sequence (string)
            reference_seq: Reference RNA sequence from PDB (string)
            reference_matrix: NxN binary matrix indicating base pairs (numpy array)
            
        Returns:
            Tuple of (normalized_score, detailed_stats)
        """
        if len(predicted_seq) != len(reference_seq):
            raise ValueError("Predicted and reference sequences must have the same length")
        
        n = len(predicted_seq)
        
        # Encode sequences (cached, uppercase)
        pred_indices = self._encode_sequence(predicted_seq)
        ref_indices = self._encode_sequence(reference_seq)
        
        # Ensure matrix is contiguous and integer-typed for Numba/vectorized path
        reference_matrix = np.ascontiguousarray(reference_matrix.astype(np.int8))

        # Use JIT-compiled version if available, otherwise vectorized fast path
        if self.use_numba:
            try:
                total_score, paired_positions, unpaired_positions, paired_correct, unpaired_correct = \
                    _score_loop_numba(
                        np.ascontiguousarray(pred_indices),
                        np.ascontiguousarray(ref_indices),
                        reference_matrix,
                        np.ascontiguousarray(self.score_array),
                        int(self.min_separation)
                    )
            except Exception as e:
                # Robust fallback: if Numba fails (e.g., typing/unboxing issues),
                # fall back to fast vectorized implementation.
                print(f"Warning: Numba JIT failed ({e}). Falling back to vectorized scoring.")
                total_score, paired_positions, unpaired_positions, paired_correct, unpaired_correct = \
                    self._score_loop_vectorized(pred_indices, ref_indices, reference_matrix)
        else:
            total_score, paired_positions, unpaired_positions, paired_correct, unpaired_correct = \
                self._score_loop_vectorized(pred_indices, ref_indices, reference_matrix)
        
        # Recompute paired-only score (structure-focused) to avoid unpaired dominance
        # Create masks and compute base scores for paired positions only
        tri_i, tri_j = self._get_triu_indices(n)
        valid_pred = (pred_indices[tri_i] >= 0) & (pred_indices[tri_j] >= 0)
        valid_ref = (ref_indices[tri_i] >= 0) & (ref_indices[tri_j] >= 0)
        valid_mask = valid_pred & valid_ref
        paired_total_score = 0
        if paired_positions > 0 and np.any(valid_mask):
            tri_i_v = tri_i[valid_mask]
            tri_j_v = tri_j[valid_mask]
            ref_has_pair = (reference_matrix[tri_i_v, tri_j_v] == 1)
            if np.any(ref_has_pair):
                pi = pred_indices[tri_i_v][ref_has_pair]
                pj = pred_indices[tri_j_v][ref_has_pair]
                ri = ref_indices[tri_i_v][ref_has_pair]
                rj = ref_indices[tri_j_v][ref_has_pair]
                base_scores = self.score_array[pi, pj, ri, rj]
                pos_mask = base_scores > 0
                if np.any(pos_mask):
                    paired_total_score = int(np.sum((base_scores[pos_mask] + 1).astype(np.int64)))

        # Paired-only normalization (structure-focused base pair score)
        max_possible_paired = 4 * paired_positions
        normalized_paired = paired_total_score / max_possible_paired if max_possible_paired > 0 else 1.0

        # Unpaired component (reported for analysis but not used in base_pair_score)
        unpaired_total_score = unpaired_positions  # sequence-only structural scoring
        max_possible_unpaired = unpaired_positions
        normalized_unpaired = (unpaired_total_score / max_possible_unpaired) if max_possible_unpaired > 0 else 1.0

        # Weighted base-pair score (default paired-only when unpaired_weight=0)
        normalized_score = (1.0 - self.unpaired_weight) * normalized_paired + self.unpaired_weight * normalized_unpaired

        # Additional metrics
        total_positions = paired_positions + unpaired_positions
        structure_accuracy = paired_correct / paired_positions if paired_positions > 0 else 1.0
        non_structure_accuracy = unpaired_correct / unpaired_positions if unpaired_positions > 0 else 1.0

        stats = {
            'normalized_score': normalized_score,  # base-pair score (paired-only by default)
            'normalized_paired_score': normalized_paired,
            'normalized_unpaired_score': normalized_unpaired,
            'raw_score_paired': paired_total_score,
            'max_possible_score_paired': max_possible_paired,
            'raw_score_unpaired': unpaired_total_score,
            'max_possible_score_unpaired': max_possible_unpaired,
            'total_positions_scored': total_positions,
            'paired_positions': paired_positions,
            'unpaired_positions': unpaired_positions,
            'paired_correct': paired_correct,
            'unpaired_correct': unpaired_correct,
            'structure_accuracy': structure_accuracy,
            'non_structure_accuracy': non_structure_accuracy,
            'pair_counts': {},  # Simplified for performance
            'total_base_pairs': paired_positions,
            'min_separation': self.min_separation,
            'unpaired_weight': self.unpaired_weight
        }
        
        return normalized_score, stats
    
    def _score_loop_python(self, predicted_seq: str, reference_seq: str, 
                          reference_matrix: np.ndarray) -> Tuple[int, int, int, int, int]:
        """Pure Python scoring loop (fallback when Numba not available)."""
        n = len(predicted_seq)
        total_score = 0
        paired_positions = 0
        unpaired_positions = 0
        paired_correct = 0
        unpaired_correct = 0
        
        for i in range(n):
            for j in range(i + 1, n):
                # CRITICAL: Skip positions that cannot physically pair
                if (j - i) < self.min_separation:
                    continue
                
                ref_has_pair = (reference_matrix[i, j] == 1)
                
                if ref_has_pair:
                    paired_positions += 1
                    pred_pair = predicted_seq[i] + predicted_seq[j]
                    ref_pair = reference_seq[i] + reference_seq[j]
                    
                    if pred_pair in self.valid_pairs:
                        base_score = self.bp_score_matrix.get((pred_pair, ref_pair), 0)
                        if base_score > 0:
                            total_score += base_score + 1
                            paired_correct += 1
                else:
                    # Sequence-only structural scoring:
                    # Without a predicted structure, do not penalize unpaired positions
                    # that happen to be complementary by sequence. Award 1 point always.
                    unpaired_positions += 1
                    total_score += 1
                    unpaired_correct += 1
        
        return total_score, paired_positions, unpaired_positions, paired_correct, unpaired_correct
    
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
        max_length = len(predicted_seq)
        
        # Convert edit distance to similarity score (1 - normalized_edit_distance)
        normalized_score = 1 - (edit_dist / max_length) if max_length > 0 else 1
        
        stats = {
            'edit_distance': edit_dist,
            'sequence_length': len(predicted_seq),
            'normalized_score': normalized_score,
            'accuracy': (max_length - edit_dist) / max_length if max_length > 0 else 1
        }
        
        return normalized_score, stats
    
    def calculate_recovery_rate_score(self, predicted_seq: str, reference_seq: str) -> Tuple[float, Dict]:
        """
        Calculate recovery rate (position-wise accuracy) between sequences.
        
        Recovery rate is the percentage of positions where the predicted base
        matches the reference base. This is a simpler, more intuitive metric
        than edit distance.
        
        Args:
            predicted_seq: Predicted RNA sequence
            reference_seq: Reference RNA sequence
            
        Returns:
            Tuple of (recovery_rate, detailed_stats)
        """
        if len(predicted_seq) != len(reference_seq):
            raise ValueError("Predicted and reference sequences must have the same length")
        
        # Vectorized character comparison (fast for long sequences)
        pred_arr = np.frombuffer(predicted_seq.upper().encode('ascii'), dtype='S1')
        ref_arr = np.frombuffer(reference_seq.upper().encode('ascii'), dtype='S1')
        correct_positions = int(np.sum(pred_arr == ref_arr))
        n = pred_arr.shape[0]
        recovery_rate = correct_positions / n if n > 0 else 1.0
        
        stats = {
            'correct_positions': correct_positions,
            'total_positions': n,
            'recovery_rate': recovery_rate,
            'accuracy': recovery_rate,
            'sequence_length': n
        }
        
        return recovery_rate, stats
    
    def calculate_combined_score(self, predicted_seq: str, reference_seq: str, 
                               reference_matrix: np.ndarray, lambda_param: float = 0.5,
                               seq_score_type: str = 'edit_distance') -> Dict:
        """
        Calculate combined score using both base pairing and sequence similarity scores.
        
        THIS IS THE MAIN FUNCTION FOR MODEL TRAINING.
        Pass sequences as strings and matrix as numpy array (no file I/O).
        
        Handles edge cases:
        - If no base pairs exist (e.g., short sequences, poly-A/U/G/C), 
          automatically falls back to sequence similarity only
        
        Args:
            predicted_seq: Predicted RNA sequence (string)
            reference_seq: Reference RNA sequence from PDB (string)
            reference_matrix: NxN binary matrix indicating base pairs (numpy array)
            lambda_param: Weight parameter (0-1) for base pairing score. 
                         bp_weight = lambda_param, ed_weight = 1 - lambda_param (default: 0.5)
            seq_score_type: Type of sequence similarity score to use:
                           'edit_distance' (default) or 'recovery_rate'
            
        Returns:
            Dictionary containing all scores and statistics
        """
        if not (0 <= lambda_param <= 1):
            raise ValueError("lambda_param must be between 0 and 1")
        
        if seq_score_type not in ['edit_distance', 'recovery_rate']:
            raise ValueError("seq_score_type must be 'edit_distance' or 'recovery_rate'")
        
        # Calculate base pair score
        bp_score, bp_stats = self.calculate_base_pair_score(predicted_seq, reference_seq, reference_matrix)
        
        # Calculate sequence similarity score based on type
        if seq_score_type == 'edit_distance':
            seq_score, seq_stats = self.calculate_edit_distance_score(predicted_seq, reference_seq)
            seq_stats_key = 'edit_distance_stats'
        else:  # recovery_rate
            seq_score, seq_stats = self.calculate_recovery_rate_score(predicted_seq, reference_seq)
            seq_stats_key = 'recovery_rate_stats'
        
        # Check for edge cases
        total_base_pairs = bp_stats['total_base_pairs']
        scoring_mode = 'combined'
        
        if total_base_pairs == 0:
            # No base pairs exist - use only sequence similarity score
            combined_score = seq_score
            scoring_mode = f'{seq_score_type}_only'
            effective_lambda = 0.0
            bp_weight = 0.0
            seq_weight = 1.0
        else:
            # Normal case - use weighted combination
            effective_lambda = lambda_param
            bp_weight = lambda_param
            seq_weight = 1 - lambda_param
            combined_score = bp_weight * bp_score + seq_weight * seq_score
        
        result = {
            'combined_score': combined_score,
            'base_pair_score': bp_score,
            'sequence_score': seq_score,
            'sequence_score_type': seq_score_type,
            'lambda_param': effective_lambda,
            'scoring_mode': scoring_mode,
            'weights': {'bp_weight': bp_weight, 'seq_weight': seq_weight},
            'base_pair_stats': bp_stats,
            seq_stats_key: seq_stats
        }
        
        # Add legacy key for backward compatibility
        if seq_score_type == 'edit_distance':
            result['edit_distance_score'] = seq_score
        
        return result
    
    def score_batch(self, predicted_seqs: List[str], reference_seqs: List[str], 
                   reference_matrices: List[np.ndarray], lambda_param: float = 0.5,
                   seq_score_type: str = 'recovery_rate') -> List[Dict]:
        """
        Score a batch of predictions efficiently.
        
        Useful for model training when you have multiple predictions to score.
        
        Args:
            predicted_seqs: List of predicted RNA sequences
            reference_seqs: List of reference RNA sequences
            reference_matrices: List of NxN numpy arrays
            lambda_param: Weight parameter
            seq_score_type: Type of sequence score ('edit_distance' or 'recovery_rate')
            
        Returns:
            List of result dictionaries
        """
        results = []
        for pred, ref, matrix in zip(predicted_seqs, reference_seqs, reference_matrices):
            result = self.calculate_combined_score(pred, ref, matrix, lambda_param, seq_score_type)
            results.append(result)
        return results
    
    def print_detailed_results(self, results: Dict, predicted_seq: str, reference_seq: str):
        """
        Print detailed scoring results in a formatted way.
        
        Args:
            results: Results dictionary from calculate_combined_score
            predicted_seq: Predicted RNA sequence
            reference_seq: Reference RNA sequence
        """
        print("="*70)
        print("RNA SEQUENCE SCORING RESULTS")
        print("="*70)
        
        print(f"\nSequence Length: {len(predicted_seq)}")
        print(f"Reference Seq:  {reference_seq}")
        print(f"Predicted Seq:  {predicted_seq}")
        
        # Show scoring mode
        seq_score_type = results.get('sequence_score_type', 'edit_distance')
        if 'only' in results['scoring_mode']:
            print(f"\n⚠ SCORING MODE: {results['scoring_mode'].replace('_', ' ').title()} (no base pairs in reference)")
        
        print(f"\n{'COMBINED SCORE:':<25} {results['combined_score']:.4f}")
        print(f"  {'Base Pair Score:':<23} {results['base_pair_score']:.4f} (λ = {results['lambda_param']:.2f})")
        print(f"  {'Sequence Score:':<23} {results['sequence_score']:.4f} (1-λ = {1-results['lambda_param']:.2f}) [{seq_score_type}]")
        
        bp_stats = results['base_pair_stats']
        print(f"\n{'BASE PAIR ANALYSIS:'}")
        print(f"  {'Total positions scored:':<30} {bp_stats['total_positions_scored']} (upper triangle only)")
        print(f"  {'Paired positions:':<30} {bp_stats['paired_positions']}")
        print(f"  {'Unpaired positions:':<30} {bp_stats['unpaired_positions']}")
        print(f"  {'Raw paired score:':<30} {bp_stats['raw_score_paired']}/{bp_stats['max_possible_score_paired']}")
        print(f"  {'Raw unpaired score:':<30} {bp_stats['raw_score_unpaired']}/{bp_stats['max_possible_score_unpaired']}")
        
        if bp_stats['paired_positions'] > 0:
            print(f"  {'Structure accuracy:':<30} {bp_stats['structure_accuracy']:.4f} ({bp_stats['paired_correct']}/{bp_stats['paired_positions']})")
        
        if bp_stats['unpaired_positions'] > 0:
            print(f"  {'Non-structure accuracy:':<30} {bp_stats['non_structure_accuracy']:.4f} ({bp_stats['unpaired_correct']}/{bp_stats['unpaired_positions']})")
        
        print(f"  {'Min separation used:':<30} {bp_stats['min_separation']} (positions with j-i < {bp_stats['min_separation']} excluded)")
        
        # Sequence analysis - handle both types
        print(f"\n{'SEQUENCE ANALYSIS:'} (using {seq_score_type})")
        if seq_score_type == 'edit_distance':
            seq_stats = results.get('edit_distance_stats', {})
            print(f"  {'Edit distance:':<30} {seq_stats.get('edit_distance', 'N/A')}")
            print(f"  {'Sequence accuracy:':<30} {seq_stats.get('accuracy', 0):.4f}")
        else:  # recovery_rate
            seq_stats = results.get('recovery_rate_stats', {})
            print(f"  {'Correct positions:':<30} {seq_stats.get('correct_positions', 'N/A')}/{seq_stats.get('total_positions', 'N/A')}")
            print(f"  {'Recovery rate:':<30} {seq_stats.get('recovery_rate', 0):.4f}")
        
        print("="*70)


# Numba JIT-compiled scoring loop (defined at module level for Numba compatibility)
@jit(nopython=True, cache=True)
def _score_loop_numba(pred_indices: np.ndarray, ref_indices: np.ndarray, 
                     reference_matrix: np.ndarray, score_array: np.ndarray,
                     min_separation: int) -> Tuple[int, int, int, int, int]:
    """
    Numba JIT-compiled scoring loop for maximum performance.
    
    This function is compiled to native code and runs 10-100x faster than Python.
    First call will be slow (compilation), subsequent calls will be very fast.
    
    CRITICAL: Excludes positions where j-i < min_separation to avoid bias.
    """
    n = len(pred_indices)
    total_score = 0
    paired_positions = 0
    unpaired_positions = 0
    paired_correct = 0
    unpaired_correct = 0
    
    # Valid pairs (encoded indices):
    # AU(0,1), UA(1,0), GC(2,3), CG(3,2), GU(2,1), UG(1,2)
    
    for i in range(n):
        for j in range(i + 1, n):
            # CRITICAL: Skip positions that cannot physically pair
            if (j - i) < min_separation:
                continue
            
            ref_has_pair = (reference_matrix[i, j] == 1)
            
            pred_i = pred_indices[i]
            pred_j = pred_indices[j]
            ref_i = ref_indices[i]
            ref_j = ref_indices[j]
            
            # Skip if invalid base index
            if pred_i < 0 or pred_j < 0 or ref_i < 0 or ref_j < 0:
                continue
            
            if ref_has_pair:
                paired_positions += 1
                
                # Check if prediction is a valid pair (Numba-friendly explicit checks)
                if (
                    (pred_i == 0 and pred_j == 1) or
                    (pred_i == 1 and pred_j == 0) or
                    (pred_i == 2 and pred_j == 3) or
                    (pred_i == 3 and pred_j == 2) or
                    (pred_i == 2 and pred_j == 1) or
                    (pred_i == 1 and pred_j == 2)
                ):
                    # Look up score from pre-computed array
                    base_score = score_array[pred_i, pred_j, ref_i, ref_j]
                    
                    if base_score > 0:
                        total_score += base_score + 1
                        paired_correct += 1
            else:
                # Sequence-only structural scoring:
                # Without a predicted structure matrix, treat all unpaired
                # positions in the reference as correct and award 1 point.
                unpaired_positions += 1
                total_score += 1
                unpaired_correct += 1
    
    return total_score, paired_positions, unpaired_positions, paired_correct, unpaired_correct


def load_reference_matrix(matrix_file: str) -> np.ndarray:
    """
    Load reference base pair matrix from file.
    
    Expected format: CSV, NPY, or text file with binary matrix (1 = paired, 0 = unpaired)
    
    Note: Avoid using this in training loops - load once and reuse.
    
    Args:
        matrix_file: Path to matrix file
        
    Returns:
        NxN numpy array with binary base pair information
    """
    try:
        if matrix_file.endswith('.npy'):
            mat = np.load(matrix_file)
        elif matrix_file.endswith('.csv'):
            mat = np.loadtxt(matrix_file, delimiter=',', dtype=np.int8)
        else:
            mat = np.loadtxt(matrix_file, dtype=np.int8)

        # Ensure 2D binary matrix and Numba-friendly dtype
        if mat.ndim != 2:
            raise ValueError(f"Reference matrix must be 2D, got shape {mat.shape}")
        mat = (mat > 0).astype(np.int8)
        return np.ascontiguousarray(mat)
    except Exception as e:
        raise ValueError(f"Error loading reference matrix: {e}")


def main():
    """Command line interface for standalone scoring."""
    parser = argparse.ArgumentParser(
        description='Score predicted RNA sequences against PDB reference structures',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Default (edit distance)
  python rna_scorer.py predicted.fasta reference.fasta matrix.csv
  
  # Using recovery rate
  python rna_scorer.py predicted.fasta reference.fasta matrix.csv --seq-score recovery_rate
  
  # Custom lambda and sequence score type
  python rna_scorer.py predicted.fasta reference.fasta matrix.csv --lambda_param 0.5 --seq-score recovery_rate
  
  # Disable Numba
  python rna_scorer.py predicted.fasta reference.fasta matrix.csv --no-numba
  
  # Custom minimum separation
  python rna_scorer.py predicted.fasta reference.fasta matrix.csv --min-sep 5
        """
    )
    
    parser.add_argument('predicted_fasta', help='Path to predicted sequence FASTA file')
    parser.add_argument('reference_sequence', help='Reference RNA sequence (string or FASTA file)')
    parser.add_argument('reference_matrix', help='Path to reference base pair matrix file')
    parser.add_argument('--lambda_param', type=float, default=0.5, 
                       help='Lambda parameter (0-1) for score weighting (default: 0.5)')
    parser.add_argument('--min-sep', type=int, default=4,
                       help='Minimum sequence separation |j-i| for base pairing (default: 4)')
    parser.add_argument('--seq-score', choices=['edit_distance', 'recovery_rate'], default='recovery_rate',
                       help='Sequence similarity score type: recovery_rate (default) or edit_distance')
    parser.add_argument('--no-numba', action='store_true',
                       help='Disable Numba JIT compilation')
    parser.add_argument('--unpaired-weight', type=float, default=0.0,
                       help='Weight (0-1) for unpaired positions within base-pair score; default 0 (paired-only)')
    
    args = parser.parse_args()
    
    # Initialize scorer
    scorer = RNASequenceScorer(use_numba=not args.no_numba, min_separation=args.min_sep, unpaired_weight=args.unpaired_weight)
    
    try:
        # Load predicted sequence
        predicted_seq = scorer.load_fasta_sequence(args.predicted_fasta)
        
        # Load reference sequence
        if args.reference_sequence.endswith(('.fasta', '.fa')):
            reference_seq = scorer.load_fasta_sequence(args.reference_sequence)
        else:
            reference_seq = args.reference_sequence.upper()
        
        # Load reference matrix
        reference_matrix = load_reference_matrix(args.reference_matrix)
        
        # Calculate scores
        results = scorer.calculate_combined_score(
            predicted_seq, reference_seq, reference_matrix,
            args.lambda_param, args.seq_score
        )
        
        # Print results
        scorer.print_detailed_results(results, predicted_seq, reference_seq)
        
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    # Example usage for model training:
    """
    # Initialize once (outside training loop)
    scorer = RNASequenceScorer(use_numba=True, min_separation=4)
    
    # Pre-load reference data (outside training loop)
    reference_matrices = [load_reference_matrix(f) for f in matrix_files]
    reference_seqs = [load_sequence(f) for f in ref_files]
    
    # Inside training loop (FAST!)
    for epoch in range(num_epochs):
        for pred_seq, ref_seq, ref_matrix in zip(predictions, reference_seqs, reference_matrices):
            # Using edit distance (default)
            result = scorer.calculate_combined_score(
                pred_seq, ref_seq, ref_matrix, 
                lambda_param=0.5, 
                seq_score_type='edit_distance'
            )
            
            # Or using recovery rate (more intuitive)
            result = scorer.calculate_combined_score(
                pred_seq, ref_seq, ref_matrix, 
                lambda_param=0.5, 
                seq_score_type='recovery_rate'
            )
            
            score = result['combined_score']
            # Use score for loss calculation
    
    # Note: min_separation should match the value used when generating reference matrices
    # to ensure consistency between structure generation and scoring.
    
    # Edit Distance vs Recovery Rate:
    # - Edit Distance: Considers insertions/deletions, uses dynamic programming
    # - Recovery Rate: Simple position-wise accuracy, more intuitive
    # - Both give same results for equal-length sequences with only substitutions
    """
    
    main()