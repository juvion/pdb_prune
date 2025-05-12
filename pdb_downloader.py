import os
from Bio.PDB import PDBList
import time
import requests
import json
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def download_rna_pdbs(destination_dir):
    """
    Searches the PDB database for entries containing RNA and downloads
    the PDB files to the specified directory.

    Args:
        destination_dir (str): The path to the directory where PDB files
                               should be downloaded. This directory will
                               be created if it doesn't exist.
    """
    # Ensure the destination directory exists
    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir)
        logger.info(f"Created directory: {destination_dir}")

    # Initialize PDBList for downloading
    try:
        pdbl = PDBList(pdb=destination_dir)
        logger.info("Successfully initialized PDBList")
    except Exception as e:
        logger.error(f"Failed to initialize PDBList: {e}")
        return

    logger.info("Searching PDB for entries containing RNA...")

    # Use PDB REST API to search for RNA structures
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    # Updated query structure for RNA
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "entity_poly.rcsb_entity_polymer_type",
                "operator": "exact_match",
                "value": "RNA"
            }
        },
        "return_type": "entry"
    }

    try:
        # Perform the search
        logger.info("Sending search request to PDB API...")
        response = requests.post(search_url, json=query)
        response.raise_for_status()
        result = response.json()
        
        # Log the raw response for debugging
        logger.debug(f"API Response: {json.dumps(result, indent=2)}")
        
        pdb_ids = [item["identifier"].lower() for item in result.get("result_set", [])]

        if not pdb_ids:
            logger.warning("No PDB entries containing RNA found with the specified query.")
            return

        logger.info(f"Found {len(pdb_ids)} PDB entries containing RNA. Starting download...")

        # Download each PDB file
        download_count = 0
        failed_downloads = []
        for pdb_id in pdb_ids:
            try:
                logger.info(f"Attempting to download {pdb_id}...")
                pdbl.retrieve_pdb_file(pdb_id, file_format='pdb')
                logger.info(f"Successfully downloaded {pdb_id}")
                download_count += 1
                # Optional: add a small delay to be polite to the server
                time.sleep(1)
            except Exception as e:
                logger.error(f"Error downloading {pdb_id}: {str(e)}")
                failed_downloads.append(pdb_id)
                time.sleep(1)

        logger.info("\n--- Download Summary ---")
        logger.info(f"Total found: {len(pdb_ids)}")
        logger.info(f"Successfully downloaded: {download_count}")
        if failed_downloads:
            logger.warning(f"Failed to download: {len(failed_downloads)} ({', '.join(failed_downloads)})")
        logger.info("-----------------------")

    except requests.exceptions.RequestException as e:
        logger.error(f"Network error during PDB search: {e}")
    except json.JSONDecodeError as e:
        logger.error(f"Failed to parse API response: {e}")
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}")

# --- Example Usage ---
if __name__ == "__main__":
    # Define the directory where you want to save the files
    download_directory = "downloaded_rna_pdbs"

    # Run the function
    download_rna_pdbs(download_directory)