import os
import logging
from Bio.PDB import PDBList
import requests
import json
import time

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def download_rna_pdbs(download_directory="downloaded_rna_pdbs", max_entries=20):
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
    download_rna_pdbs(max_entries=20)