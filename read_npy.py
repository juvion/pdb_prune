#!/usr/bin/env python3

import numpy as np
import sys

def print_npy_file(file_path):
    """Read and print a NumPy array from a .npy file.
    
    Args:
        file_path (str): Path to the .npy file
    """
    try:
        # Set NumPy print options to show full array
        np.set_printoptions(threshold=np.inf,  # Show all elements
                          linewidth=np.inf,    # Don't wrap lines
                          precision=3,         # Show 3 decimal places
                          suppress=True)       # Don't use scientific notation
        
        # Load the array
        array = np.load(file_path)
        
        # Print array information
        print("\nArray Information:")
        print(f"Shape: {array.shape}")
        print(f"Data type: {array.dtype}")
        print(f"Number of dimensions: {array.ndim}")
        
        # Print the array
        print("\nArray Contents:")
        print(array)
        
        # If it's a 3D array (like in our RNA case), print more detailed info
        if array.ndim == 3:
            print("\nDetailed Information:")
            print(f"Number of residues: {array.shape[0]}")
            print(f"Number of atoms per residue: {array.shape[1]}")
            print(f"Coordinates per atom: {array.shape[2]}")
            
            # Print first residue as example
            print("\nFirst residue coordinates:")
            print(array[0])
            
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

def main():
    if len(sys.argv) != 2:
        print("Usage: python read_npy.py <path_to_npy_file>")
        sys.exit(1)
    
    file_path = sys.argv[1]
    print_npy_file(file_path)

if __name__ == "__main__":
    main() 