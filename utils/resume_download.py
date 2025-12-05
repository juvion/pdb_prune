#!/usr/bin/env python3
"""
Script to resume RNA PDB downloads from where they left off.
This script helps recover from interrupted downloads and provides flexible options.
"""

import os
import argparse
from utils.pdb_downloader import download_rna_pdbs

def count_existing_files(directory):
    """Count existing PDB files in the directory."""
    if not os.path.exists(directory):
        return 0
    
    count = 0
    for filename in os.listdir(directory):
        if filename.endswith(('.pdb', '.ent', '.cif')):
            count += 1
    return count

def main():
    parser = argparse.ArgumentParser(description='Resume RNA PDB downloads')
    parser.add_argument('--output_dir', type=str, default="downloaded_rna_pdbs2", 
                       help='Output directory for downloaded files')
    parser.add_argument('--target_total', type=int, default=7000, 
                       help='Target total number of files to download')
    parser.add_argument('--batch_size', type=int, default=50, 
                       help='Number of entries to process in each batch')
    parser.add_argument('--max_retries', type=int, default=3, 
                       help='Maximum retry attempts for failed downloads')
    parser.add_argument('--check_only', action='store_true', 
                       help='Only check current status without downloading')
    
    args = parser.parse_args()
    
    # Check current status
    existing_count = count_existing_files(args.output_dir)
    remaining = max(0, args.target_total - existing_count)
    
    print(f"\n=== Download Status ===")
    print(f"Output directory: {args.output_dir}")
    print(f"Existing files: {existing_count}")
    print(f"Target total: {args.target_total}")
    print(f"Remaining to download: {remaining}")
    print(f"======================\n")
    
    if args.check_only:
        return
    
    if remaining <= 0:
        print("✅ Download already complete!")
        return
    
    # Resume download
    print(f"🚀 Resuming download of {remaining} files...")
    
    try:
        download_rna_pdbs(
            download_directory=args.output_dir,
            max_entries=remaining,
            start_from=existing_count,
            batch_size=args.batch_size,
            max_retries=args.max_retries
        )
    except KeyboardInterrupt:
        print("\n⏸️  Download interrupted by user")
        final_count = count_existing_files(args.output_dir)
        print(f"Files downloaded in this session: {final_count - existing_count}")
        print(f"Total files now: {final_count}")
        print(f"Remaining: {args.target_total - final_count}")
    except Exception as e:
        print(f"\n❌ Error during download: {e}")
        final_count = count_existing_files(args.output_dir)
        print(f"Files downloaded before error: {final_count - existing_count}")
        print(f"Total files now: {final_count}")

if __name__ == "__main__":
    main()