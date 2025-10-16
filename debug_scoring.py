#!/usr/bin/env python3

import numpy as np

# Load the reference matrix
reference_matrix = np.loadtxt('/Users/ju/Documents/Dev/pdb_prune/pdb1aqo.mat')
sequence = "GGAGUGCUUCAACAGUGCUUGGACGCUCC"

# Base pair scoring matrix
bp_score_matrix = {
    ('AU', 'AU'): 3, ('AU', 'UA'): 2,
    ('GC', 'GC'): 3, ('GC', 'CG'): 2,
    ('GU', 'GU'): 3, ('UG', 'GU'): 2,
    ('GC', 'GU'): 1, ('GC', 'AU'): 1,
    ('AU', 'GC'): 1, ('AU', 'GU'): 1,
    ('GU', 'GC'): 1, ('GU', 'AU'): 1
}

valid_pairs = {'AU', 'UA', 'GC', 'CG', 'GU', 'UG'}

print("=== Base Pair Scoring Analysis ===")
print(f"Sequence: {sequence}")
print(f"Length: {len(sequence)}")
print()

# Find all base pairs from reference matrix
base_pairs = []
for i in range(len(sequence)):
    for j in range(i + 1, len(sequence)):
        if reference_matrix[i, j] == 1:
            base_pairs.append((i, j))

print(f"Total base pairs in reference structure: {len(base_pairs)}")
print("Base pairs:")
for i, j in base_pairs:
    pair = sequence[i] + sequence[j]
    print(f"  Position {i+1}-{j+1}: {sequence[i]}{sequence[j]} ({pair})")

print()
print("=== Scoring Analysis ===")
total_score = 0
max_possible_score = 3 * len(base_pairs)

for i, j in base_pairs:
    pred_pair = sequence[i] + sequence[j]  # Same as reference since sequences are identical
    ref_pair = sequence[i] + sequence[j]
    
    if pred_pair in valid_pairs:
        score = bp_score_matrix.get((pred_pair, ref_pair), 0)
        total_score += score
        print(f"  {pred_pair} -> {ref_pair}: score = {score}")
    else:
        print(f"  {pred_pair} -> {ref_pair}: INVALID PAIR, score = 0")

print()
print(f"Total raw score: {total_score}")
print(f"Max possible score: {max_possible_score}")
print(f"Normalized score: {total_score / max_possible_score:.4f}")
print()
print("=== Why the score is not 1.0 ===")
print("The scoring matrix gives different scores for different pair types:")
print("- Perfect match (e.g., GC->GC): 3 points")
print("- Reverse match (e.g., GC->CG): 2 points")
print("- Cross matches: 1 point")
print("- Invalid pairs: 0 points")
print()
print("Even with identical sequences, some base pairs get less than the maximum 3 points!")