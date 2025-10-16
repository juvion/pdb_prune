#!/usr/bin/env python3
"""
Example usage of the PDB base pair analysis tool.

This script demonstrates how to use the PDBBasePairAnalyzer class
and the convenience function analyze_base_pairs.
"""

import sys
import os
import numpy as np
from pdb_basepair import PDBBasePairAnalyzer, analyze_base_pairs


def example_basic_usage():
    """Demonstrate basic usage with the convenience function."""
    print("=== Basic Usage Example ===")
    
    # Example PDB file path (replace with actual file)
    pdb_file = "example_structure.pdb"
    
    if not os.path.exists(pdb_file):
        print(f"Example PDB file '{pdb_file}' not found.")
        print("Please provide a valid PDB file path.")
        return
    
    try:
        # Use the convenience function
        matrix, residue_info = analyze_base_pairs(pdb_file, distance_cutoff=10.0)
        
        print(f"Base pairing matrix shape: {matrix.shape}")
        print(f"Number of residues: {len(residue_info)}")
        print(f"Number of base pairs: {np.sum(matrix) // 2}")
        
        # Display the matrix
        print("\nBase pairing matrix:")
        print(matrix.astype(int))
        
    except Exception as e:
        print(f"Error in basic usage: {e}")


def example_advanced_usage():
    """Demonstrate advanced usage with the PDBBasePairAnalyzer class."""
    print("\n=== Advanced Usage Example ===")
    
    # Example PDB file path (replace with actual file)
    pdb_file = "example_structure.pdb"
    
    if not os.path.exists(pdb_file):
        print(f"Example PDB file '{pdb_file}' not found.")
        print("Please provide a valid PDB file path.")
        return
    
    try:
        # Initialize analyzer with custom parameters
        analyzer = PDBBasePairAnalyzer(
            distance_cutoff=12.0,  # Larger distance cutoff
            chains=['A', 'B']      # Analyze specific chains
        )
        
        # Load and analyze PDB
        analyzer.load_pdb(pdb_file)
        matrix = analyzer.get_base_pair_matrix()
        
        # Print summary
        analyzer.print_summary()
        
        # Get detailed pairing information
        pairs = analyzer.get_base_pairs()
        print(f"\nDetailed base pair information:")
        for i, pair in enumerate(pairs[:5]):  # Show first 5 pairs
            print(f"  {i+1}. {pair['res1']}({pair['base1']}) - "
                  f"{pair['res2']}({pair['base2']}) "
                  f"(distance: {pair['distance']:.2f} Å)")
        
        if len(pairs) > 5:
            print(f"  ... and {len(pairs) - 5} more pairs")
        
        # Visualize matrix
        analyzer.plot_matrix(save_path='basepair_matrix.png')
        
    except Exception as e:
        print(f"Error in advanced usage: {e}")


def example_batch_processing():
    """Demonstrate batch processing of multiple PDB files."""
    print("\n=== Batch Processing Example ===")
    
    # List of PDB files to process
    pdb_files = [
        "structure1.pdb",
        "structure2.pdb", 
        "structure3.pdb"
    ]
    
    results = []
    
    for pdb_file in pdb_files:
        if os.path.exists(pdb_file):
            try:
                print(f"Processing {pdb_file}...")
                matrix, residue_info = analyze_base_pairs(pdb_file)
                
                n_pairs = np.sum(matrix) // 2
                results.append({
                    'file': pdb_file,
                    'n_residues': len(residue_info),
                    'n_pairs': n_pairs,
                    'matrix': matrix
                })
                
                print(f"  Residues: {len(residue_info)}, Base pairs: {n_pairs}")
                
            except Exception as e:
                print(f"  Error processing {pdb_file}: {e}")
        else:
            print(f"  File {pdb_file} not found, skipping...")
    
    # Summary of batch processing
    if results:
        print(f"\nBatch processing summary:")
        total_residues = sum(r['n_residues'] for r in results)
        total_pairs = sum(r['n_pairs'] for r in results)
        print(f"  Files processed: {len(results)}")
        print(f"  Total residues: {total_residues}")
        print(f"  Total base pairs: {total_pairs}")


def example_custom_analysis():
    """Demonstrate custom analysis and filtering."""
    print("\n=== Custom Analysis Example ===")
    
    pdb_file = "example_structure.pdb"
    
    if not os.path.exists(pdb_file):
        print(f"Example PDB file '{pdb_file}' not found.")
        return
    
    try:
        analyzer = PDBBasePairAnalyzer(distance_cutoff=10.0)
        analyzer.load_pdb(pdb_file)
        matrix = analyzer.get_base_pair_matrix()
        pairs = analyzer.get_base_pairs()
        
        # Filter pairs by distance
        close_pairs = [p for p in pairs if p['distance'] < 8.0]
        print(f"Pairs with distance < 8.0 Å: {len(close_pairs)}")
        
        # Filter pairs by base type
        gc_pairs = [p for p in pairs if set([p['base1'], p['base2']]) == {'G', 'C'}]
        au_pairs = [p for p in pairs if set([p['base1'], p['base2']]) == {'A', 'U'}]
        gu_pairs = [p for p in pairs if set([p['base1'], p['base2']]) == {'G', 'U'}]
        
        print(f"G-C pairs: {len(gc_pairs)}")
        print(f"A-U pairs: {len(au_pairs)}")
        print(f"G-U pairs: {len(gu_pairs)}")
        
        # Calculate base pair density
        n_residues = len(analyzer.residues)
        if n_residues > 0:
            density = len(pairs) / n_residues
            print(f"Base pair density: {density:.3f} pairs per residue")
        
    except Exception as e:
        print(f"Error in custom analysis: {e}")


def main():
    """Run all examples."""
    print("PDB Base Pair Analysis - Example Usage")
    print("=" * 50)
    
    # Check if a PDB file is provided as command line argument
    if len(sys.argv) > 1:
        pdb_file = sys.argv[1]
        if os.path.exists(pdb_file):
            print(f"Using PDB file: {pdb_file}")
            # Update the example file paths
            globals()['pdb_file'] = pdb_file
        else:
            print(f"PDB file '{pdb_file}' not found.")
            return
    
    # Run examples
    example_basic_usage()
    example_advanced_usage()
    example_batch_processing()
    example_custom_analysis()
    
    print("\n" + "=" * 50)
    print("Examples completed!")


if __name__ == "__main__":
    main()
