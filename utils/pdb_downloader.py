#!/usr/bin/env python3

import os
import logging
from Bio.PDB import PDBList, PDBParser
import requests
import json
import time
from pathlib import Path
from typing import List, Dict, Any

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

class PDBDownloader:
    def __init__(self, output_dir: str = "raw_pdbs"):
        """Initialize the PDB downloader.
        
        Args:
            output_dir (str): Directory to save downloaded PDB files
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.pdb_parser = PDBParser(QUIET=True)
        
    def download_pdbs(self, pdb_ids: List[str], max_pdbs: int = None) -> None:
        """Download PDB files for a list of PDB IDs.
        
        Args:
            pdb_ids (List[str]): List of PDB IDs to download
            max_pdbs (int, optional): Maximum number of PDBs to download
        """
        if max_pdbs is not None:
            pdb_ids = pdb_ids[:max_pdbs]
            
        print(f"\nDownloading {len(pdb_ids)} PDB files...")
        for pdb_id in pdb_ids:
            try:
                self.download_single_pdb(pdb_id)
            except Exception as e:
                print(f"Error downloading {pdb_id}: {e}")
                continue
    
    def download_from_file(self, file_path: str, max_pdbs: int = None) -> None:
        """Download PDB files from a list in a text file.
        
        Args:
            file_path (str): Path to file containing PDB IDs (one per line)
            max_pdbs (int, optional): Maximum number of PDBs to download
        """
        try:
            with open(file_path, 'r') as f:
                pdb_ids = [line.strip() for line in f if line.strip()]
            
            if max_pdbs is not None:
                pdb_ids = pdb_ids[:max_pdbs]
                
            print(f"\nDownloading {len(pdb_ids)} PDB files from {file_path}...")
            self.download_pdbs(pdb_ids)
            
        except FileNotFoundError:
            print(f"Error: File not found: {file_path}")
        except Exception as e:
            print(f"Error reading file: {e}")
    
    def search_and_download(self, criteria: Dict[str, Any], max_pdbs: int = None) -> None:
        """Search for PDB files matching criteria and download them.
        
        Args:
            criteria (Dict[str, Any]): Search criteria
            max_pdbs (int, optional): Maximum number of PDBs to download
        """
        try:
            # Search for PDB files matching criteria
            pdb_ids = self.search_pdb(criteria)
            
            if max_pdbs is not None:
                pdb_ids = pdb_ids[:max_pdbs]
                
            print(f"\nFound {len(pdb_ids)} matching PDB files")
            self.download_pdbs(pdb_ids)
            
        except Exception as e:
            print(f"Error searching PDB: {e}")
    
    def search_pdb(self, criteria: Dict[str, Any]) -> List[str]:
        """Search for PDB files matching criteria.
        
        Args:
            criteria (Dict[str, Any]): Search criteria including:
                - resolution: Maximum resolution
                - rna_only: Whether to include only RNA structures
                - min_length: Minimum sequence length
                
        Returns:
            List[str]: List of matching PDB IDs
        """
        # Example implementation - replace with actual PDB search API
        # This is a placeholder that returns some example PDB IDs
        return ['1ABC', '2XYZ', '3DEF']
    
    def download_single_pdb(self, pdb_id: str) -> None:
        """Download a single PDB file.
        
        Args:
            pdb_id (str): PDB ID to download
        """
        pdb_id = pdb_id.lower()
        output_file = self.output_dir / f"{pdb_id}.pdb"
        
        if output_file.exists():
            print(f"File already exists: {output_file}")
            return
        
        try:
            # First try PDB format
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            response = requests.get(url)
            
            # If PDB format fails, try mmCIF format
            if response.status_code == 404:
                url = f"https://files.rcsb.org/download/{pdb_id}.cif"
                response = requests.get(url)
                response.raise_for_status()
                
                # Save as CIF file
                output_file = self.output_dir / f"{pdb_id}.cif"
                with open(output_file, 'w') as f:
                    f.write(response.text)
                print(f"Downloaded: {output_file}")
            else:
                response.raise_for_status()
                # Save as PDB file
                with open(output_file, 'w') as f:
                    f.write(response.text)
                print(f"Downloaded: {output_file}")
            
        except requests.exceptions.RequestException as e:
            print(f"Error downloading {pdb_id}: {e}")
            if output_file.exists():
                output_file.unlink()

def download_rna_pdbs(download_directory="data/download_pdbs/downloaded_rna_pdbs", max_entries=20, start_from=0, batch_size=100, max_retries=3, search_start_page=0):
    """
    Download PDB files containing RNA structures.
    
    Args:
        download_directory (str): Directory to save downloaded PDB files
        max_entries (int): Maximum number of PDB entries to download
        start_from (int): Starting index for resuming downloads
        batch_size (int): Number of entries to process in each batch
        max_retries (int): Maximum number of retry attempts for failed downloads
        search_start_page (int): Starting page for PDB API search (for resuming search)
    """
    # Initialize PDBList
    pdbl = PDBList()
    logging.info("Successfully initialized PDBList")
    
    # Create download directory if it doesn't exist
    os.makedirs(download_directory, exist_ok=True)
    
    # Search for RNA entries using RCSB PDB API
    logging.info("Searching PDB for entries containing RNA...")
    
    # API endpoint for search
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    # Query body
    query_body = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "entity_poly.rcsb_entity_polymer_type",
                "operator": "exact_match",
                "value": "RNA"
            }
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": 100
            }
        }
    }
    
    # Headers
    headers = {
        "Content-Type": "application/json"
    }
    
    # Collect unique PDB IDs
    unique_pdb_ids = set()
    search_start = search_start_page * 100  # Convert page to start index
    
    # Load existing PDB IDs if resuming
    pdb_ids_file = os.path.join(download_directory, "found_pdb_ids.txt")
    if os.path.exists(pdb_ids_file) and search_start_page > 0:
        logging.info(f"Loading existing PDB IDs from {pdb_ids_file}...")
        try:
            with open(pdb_ids_file, 'r') as f:
                existing_ids = {line.strip() for line in f if line.strip()}
                unique_pdb_ids.update(existing_ids)
                logging.info(f"Loaded {len(existing_ids)} existing PDB IDs")
        except Exception as e:
            logging.warning(f"Could not load existing PDB IDs: {e}")
    
    max_search_attempts = 50  # Limit total search attempts to prevent infinite loops
    search_attempts = 0
    
    while len(unique_pdb_ids) < max_entries and search_attempts < max_search_attempts:
        retry_count = 0
        success = False
        search_attempts += 1
        
        while retry_count < max_retries and not success:
            try:
                # Update pagination in the body
                query_body["request_options"]["paginate"]["start"] = search_start
                
                # Send request
                current_page = search_start//100 + 1
                logging.info(f"Sending search request to PDB API for page {current_page} (start={search_start})... (attempt {retry_count + 1})")
                response = requests.post(
                    search_url,
                    json=query_body,
                    headers=headers,
                    timeout=45  # Increased timeout
                )
                response.raise_for_status()
                
                # Parse response
                data = response.json()
                
                # Extract PDB IDs
                pdb_ids = [entry["identifier"].lower() for entry in data.get("result_set", [])]
                logging.info(f"PDB IDs from page {current_page}: {len(pdb_ids)} entries")
                
                # Add to unique set
                before_count = len(unique_pdb_ids)
                unique_pdb_ids.update(pdb_ids)
                new_ids_added = len(unique_pdb_ids) - before_count
                logging.info(f"Added {new_ids_added} new unique PDB IDs. Total: {len(unique_pdb_ids)}")
                
                # Save progress periodically
                if current_page % 5 == 0:  # Save every 5 pages
                    with open(pdb_ids_file, 'w') as f:
                        for pdb_id in sorted(unique_pdb_ids):
                            f.write(f"{pdb_id}\n")
                    logging.info(f"Progress saved to {pdb_ids_file}")
                
                # If we got fewer results than requested, we've reached the end
                if len(pdb_ids) < query_body["request_options"]["paginate"]["rows"]:
                    logging.info(f"Reached end of search results at page {current_page}")
                    success = True
                    break
                    
                success = True
                
            except (requests.exceptions.RequestException, json.JSONDecodeError) as e:
                retry_count += 1
                logging.error(f"Error during PDB search (attempt {retry_count}): {e}")
                if retry_count < max_retries:
                    wait_time = min(2 ** retry_count, 30)  # Cap wait time at 30 seconds
                    logging.info(f"Retrying in {wait_time} seconds...")
                    time.sleep(wait_time)
                else:
                    logging.error(f"Failed to fetch page {current_page} after {max_retries} attempts. Stopping search.")
                    break
            except KeyboardInterrupt:
                logging.info("Search interrupted by user. Saving progress...")
                with open(pdb_ids_file, 'w') as f:
                    for pdb_id in sorted(unique_pdb_ids):
                        f.write(f"{pdb_id}\n")
                logging.info(f"Progress saved to {pdb_ids_file}")
                raise
        
        if not success:
            logging.warning(f"Failed to fetch page {current_page}. Continuing with available {len(unique_pdb_ids)} PDB IDs.")
            break
            
        # Increment start for next page
        search_start += query_body["request_options"]["paginate"]["rows"]
        
        # Add a small delay to avoid overwhelming the API
        time.sleep(2)  # Increased delay
    
    # Save final list of found PDB IDs
    with open(pdb_ids_file, 'w') as f:
        for pdb_id in sorted(unique_pdb_ids):
            f.write(f"{pdb_id}\n")
    logging.info(f"Final PDB IDs list saved to {pdb_ids_file}")
    
    if len(unique_pdb_ids) < max_entries:
        logging.warning(f"Only found {len(unique_pdb_ids)} unique PDB entries, less than requested {max_entries}.")
    
    # Convert to list and apply start_from offset
    pdb_ids_list = list(unique_pdb_ids)
    if start_from > 0:
        pdb_ids_list = pdb_ids_list[start_from:]
        logging.info(f"Resuming download from index {start_from}")
    
    # Limit to max_entries
    if len(pdb_ids_list) > max_entries:
        pdb_ids_list = pdb_ids_list[:max_entries]
    
    # Download the PDB files in batches
    logging.info(f"Found {len(unique_pdb_ids)} unique PDB entries containing RNA. Starting download of {len(pdb_ids_list)} entries...")
    downloaded_count = 0
    failed_downloads = []
    processed_count = 0
    
    for i in range(0, len(pdb_ids_list), batch_size):
        batch = pdb_ids_list[i:i + batch_size]
        logging.info(f"Processing batch {i//batch_size + 1}: entries {i+1} to {min(i+batch_size, len(pdb_ids_list))}")
        
        for pdb_id in batch:
            processed_count += 1
            
            # Check if file already exists
            existing_files = [
                os.path.join(download_directory, f"pdb{pdb_id}.ent"),
                os.path.join(download_directory, f"{pdb_id}.pdb"),
                os.path.join(download_directory, f"{pdb_id}.cif")
            ]
            
            if any(os.path.exists(f) for f in existing_files):
                logging.info(f"File already exists for {pdb_id}, skipping...")
                downloaded_count += 1
                continue
            
            retry_count = 0
            success = False
            
            while retry_count < max_retries and not success:
                try:
                    logging.info(f"Attempting to download {pdb_id}... (attempt {retry_count + 1})")
                    pdbl.retrieve_pdb_file(
                        pdb_id,
                        pdir=download_directory,
                        file_format="pdb"
                    )
                    logging.info(f"Successfully downloaded {pdb_id}")
                    downloaded_count += 1
                    success = True
                    
                except Exception as e:
                    retry_count += 1
                    logging.error(f"Error downloading {pdb_id} (attempt {retry_count}): {e}")
                    if retry_count < max_retries:
                        wait_time = 2 ** retry_count  # Exponential backoff
                        logging.info(f"Retrying in {wait_time} seconds...")
                        time.sleep(wait_time)
                    else:
                        logging.error(f"Failed to download {pdb_id} after {max_retries} attempts")
                        failed_downloads.append(pdb_id)
            
            # Small delay between downloads to be respectful to the server
            time.sleep(0.5)
        
        # Longer delay between batches
        if i + batch_size < len(pdb_ids_list):
            logging.info(f"Completed batch {i//batch_size + 1}. Waiting 5 seconds before next batch...")
            time.sleep(5)
    
    # Print summary
    logging.info("\n--- Download Summary ---")
    logging.info(f"Total unique found: {len(unique_pdb_ids)}")
    logging.info(f"Total processed: {processed_count}")
    logging.info(f"Successfully downloaded: {downloaded_count}")
    logging.info(f"Failed downloads: {len(failed_downloads)}")
    if failed_downloads:
        logging.info(f"Failed PDB IDs: {failed_downloads[:10]}{'...' if len(failed_downloads) > 10 else ''}")
    logging.info("-----------------------")
    
    # Save failed downloads to file for retry
    if failed_downloads:
        failed_file = os.path.join(download_directory, "failed_downloads.txt")
        with open(failed_file, 'w') as f:
            for pdb_id in failed_downloads:
                f.write(f"{pdb_id}\n")
        logging.info(f"Failed downloads saved to: {failed_file}")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Download RNA PDB files')
    parser.add_argument('--max_entries', type=int, default=7000, help='Maximum number of entries to download')
    parser.add_argument('--start_from', type=int, default=0, help='Starting index for resuming downloads')
    parser.add_argument('--batch_size', type=int, default=50, help='Number of entries to process in each batch')
    parser.add_argument('--max_retries', type=int, default=3, help='Maximum retry attempts for failed downloads')
    parser.add_argument('--output_dir', type=str, default="downloaded_rna_pdbs2", help='Output directory for downloaded files')
    parser.add_argument('--search_start_page', type=int, default=0, help='Starting page for PDB API search (for resuming search)')
    
    args = parser.parse_args()
    
    download_rna_pdbs(
        download_directory=args.output_dir,
        max_entries=args.max_entries,
        start_from=args.start_from,
        batch_size=args.batch_size,
        max_retries=args.max_retries,
        search_start_page=args.search_start_page
    )