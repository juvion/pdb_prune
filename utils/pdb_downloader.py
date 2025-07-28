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

def download_rna_pdbs(download_directory="data/download_pdbs/downloaded_rna_pdbs", max_entries=20):
    """
    Download PDB files containing RNA structures.
    
    Args:
        download_directory (str): Directory to save downloaded PDB files
        max_entries (int): Maximum number of PDB entries to download
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
    start = 0
    
    while len(unique_pdb_ids) < max_entries:
        try:
            # Update pagination in the body
            query_body["request_options"]["paginate"]["start"] = start
            
            # Send request
            logging.info(f"Sending search request to PDB API for start={start}...")
            response = requests.post(
                search_url,
                json=query_body,
                headers=headers
            )
            response.raise_for_status()
            
            # Parse response
            data = response.json()
            
            # Extract PDB IDs
            pdb_ids = [entry["identifier"].lower() for entry in data.get("result_set", [])]
            logging.info(f"PDB IDs from page {start//100 + 1}: {pdb_ids}")
            
            # Add to unique set
            unique_pdb_ids.update(pdb_ids)
            logging.info(f"Current unique PDB IDs count: {len(unique_pdb_ids)}")
            
            # If we got fewer results than requested, we've reached the end
            if len(pdb_ids) < query_body["request_options"]["paginate"]["rows"]:
                break
                
            # Increment start for next page
            start += query_body["request_options"]["paginate"]["rows"]
            
            # Add a small delay to avoid overwhelming the API
            time.sleep(1)
            
        except requests.exceptions.RequestException as e:
            logging.error(f"Network error during PDB search: {e}")
            break
        except json.JSONDecodeError as e:
            logging.error(f"Error parsing API response: {e}")
            break
    
    if len(unique_pdb_ids) < max_entries:
        logging.warning(f"Only found {len(unique_pdb_ids)} unique PDB entries, less than requested {max_entries}.")
    
    # Download the PDB files
    logging.info(f"Found {len(unique_pdb_ids)} unique PDB entries containing RNA. Starting download...")
    downloaded_count = 0
    
    for pdb_id in unique_pdb_ids:
        if downloaded_count >= max_entries:
            break
            
        try:
            logging.info(f"Attempting to download {pdb_id}...")
            pdbl.retrieve_pdb_file(
                pdb_id,
                pdir=download_directory,
                file_format="pdb"
            )
            logging.info(f"Successfully downloaded {pdb_id}")
            downloaded_count += 1
        except Exception as e:
            logging.error(f"Error downloading {pdb_id}: {e}")
    
    # Print summary
    logging.info("\n--- Download Summary ---")
    logging.info(f"Total unique found: {len(unique_pdb_ids)}")
    logging.info(f"Successfully downloaded: {downloaded_count}")
    logging.info("-----------------------")

if __name__ == "__main__":
    download_rna_pdbs(max_entries=5000)