#!/usr/bin/env python3
"""
Helper script to resume PDB search and download from where it left off.
This script can analyze the current state and suggest the best resume parameters.
"""

import os
import argparse
from utils.pdb_downloader import download_rna_pdbs

def analyze_progress(output_dir):
    """
    Analyze the current download progress and suggest resume parameters.
    
    Args:
        output_dir (str): Directory containing downloaded files
        
    Returns:
        dict: Analysis results with suggested parameters
    """
    results = {
        'found_pdb_ids': 0,
        'downloaded_files': 0,
        'failed_downloads': 0,
        'suggested_search_page': 0,
        'suggested_start_from': 0
    }
    
    # Check found PDB IDs
    pdb_ids_file = os.path.join(output_dir, "found_pdb_ids.txt")
    if os.path.exists(pdb_ids_file):
        with open(pdb_ids_file, 'r') as f:
            pdb_ids = [line.strip() for line in f if line.strip()]
            results['found_pdb_ids'] = len(pdb_ids)
            # Estimate search page based on found IDs (assuming ~100 IDs per page)
            results['suggested_search_page'] = len(pdb_ids) // 100
    
    # Count downloaded files
    if os.path.exists(output_dir):
        downloaded_files = []
        for file in os.listdir(output_dir):
            if file.endswith(('.pdb', '.ent', '.cif')):
                downloaded_files.append(file)
        results['downloaded_files'] = len(downloaded_files)
        results['suggested_start_from'] = len(downloaded_files)
    
    # Check failed downloads
    failed_file = os.path.join(output_dir, "failed_downloads.txt")
    if os.path.exists(failed_file):
        with open(failed_file, 'r') as f:
            failed_ids = [line.strip() for line in f if line.strip()]
            results['failed_downloads'] = len(failed_ids)
    
    return results

def main():
    parser = argparse.ArgumentParser(description='Resume PDB search and download')
    parser.add_argument('--output_dir', type=str, default="downloaded_rna_pdbs2", 
                       help='Output directory for downloaded files')
    parser.add_argument('--max_entries', type=int, default=7000, 
                       help='Maximum number of entries to download')
    parser.add_argument('--batch_size', type=int, default=50, 
                       help='Number of entries to process in each batch')
    parser.add_argument('--max_retries', type=int, default=3, 
                       help='Maximum retry attempts for failed downloads')
    parser.add_argument('--search_start_page', type=int, default=None, 
                       help='Starting page for PDB API search (auto-detected if not specified)')
    parser.add_argument('--start_from', type=int, default=None, 
                       help='Starting index for resuming downloads (auto-detected if not specified)')
    parser.add_argument('--analyze_only', action='store_true', 
                       help='Only analyze progress without starting download')
    
    args = parser.parse_args()
    
    # Analyze current progress
    print(f"Analyzing progress in {args.output_dir}...")
    progress = analyze_progress(args.output_dir)
    
    print("\n=== Progress Analysis ===")
    print(f"Found PDB IDs: {progress['found_pdb_ids']}")
    print(f"Downloaded files: {progress['downloaded_files']}")
    print(f"Failed downloads: {progress['failed_downloads']}")
    print(f"Suggested search start page: {progress['suggested_search_page']}")
    print(f"Suggested download start index: {progress['suggested_start_from']}")
    print("========================\n")
    
    if args.analyze_only:
        print("Analysis complete. Use the suggested parameters to resume:")
        print(f"python utils/pdb_downloader.py --output_dir {args.output_dir} \\")
        print(f"  --search_start_page {progress['suggested_search_page']} \\")
        print(f"  --start_from {progress['suggested_start_from']} \\")
        print(f"  --max_entries {args.max_entries}")
        return
    
    # Use auto-detected values if not specified
    search_start_page = args.search_start_page if args.search_start_page is not None else progress['suggested_search_page']
    start_from = args.start_from if args.start_from is not None else progress['suggested_start_from']
    
    print(f"Resuming download with:")
    print(f"  Search start page: {search_start_page}")
    print(f"  Download start index: {start_from}")
    print(f"  Max entries: {args.max_entries}")
    print(f"  Batch size: {args.batch_size}")
    print("\nStarting download...\n")
    
    try:
        download_rna_pdbs(
            download_directory=args.output_dir,
            max_entries=args.max_entries,
            start_from=start_from,
            batch_size=args.batch_size,
            max_retries=args.max_retries,
            search_start_page=search_start_page
        )
    except KeyboardInterrupt:
        print("\nDownload interrupted by user.")
        print("Progress has been saved. You can resume using the same command.")
    except Exception as e:
        print(f"\nError during download: {e}")
        print("You can try resuming with the same parameters.")

if __name__ == "__main__":
    main()