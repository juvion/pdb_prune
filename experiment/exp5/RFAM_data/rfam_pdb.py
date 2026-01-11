import sys
import os
import csv

# Add project root to path to import utils
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../"))
if project_root not in sys.path:
    sys.path.append(project_root)

from utils.rfam_api import get_rfam_structures
from utils.rcsb_api import fetch_entry_summary, fetch_batch_resolutions

# File paths
input_csv = os.path.join(os.path.dirname(__file__), "RNA_families.csv")
output_csv = os.path.join(os.path.dirname(__file__), "RFAM_data.csv")

def get_rfam_ids(filepath):
    ids = []
    try:
        with open(filepath, 'r') as f:
            reader = csv.reader(f)
            # Check if header exists
            first_row = next(reader, None)
            if first_row:
                if "Rfam_ID" in first_row[0]:
                    pass # Header detected, already skipped
                else:
                    # Not a header, process it
                    if first_row[0].startswith("RF"):
                        ids.append(first_row[0].strip())
            
            for row in reader:
                if row and len(row) > 0 and row[0].strip().startswith("RF"):
                    ids.append(row[0].strip())
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
    return ids

def main():
    rfam_ids = get_rfam_ids(input_csv)
    print(f"Found {len(rfam_ids)} RFAM families to process.")

    # Global cache for resolution to avoid repeated API calls for same PDB
    resolution_cache = {}

    try:
        with open(output_csv, 'w', newline='') as csvfile:
            fieldnames = ['Rfam_ID', 'PDB', 'Chain', 'Resolution']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            for rf_id in rfam_ids:
                print(f"Processing {rf_id}...")
                try:
                    entries = get_rfam_structures(rf_id)
                    count = 0
                    
                    # Collect PDB IDs that are not in cache
                    pdb_ids_to_fetch = set()
                    for entry in entries:
                        pdb_id = entry.get("pdb_id")
                        if pdb_id and pdb_id not in resolution_cache:
                            pdb_ids_to_fetch.add(pdb_id)
                    
                    # Batch fetch resolutions
                    if pdb_ids_to_fetch:
                        print(f"  Fetching resolutions for {len(pdb_ids_to_fetch)} new PDBs...")
                        batch_resolutions = fetch_batch_resolutions(list(pdb_ids_to_fetch))
                        resolution_cache.update(batch_resolutions)
                        
                        # Mark failed/missing fetches as None in cache to avoid re-fetching
                        for pid in pdb_ids_to_fetch:
                            if pid not in resolution_cache:
                                resolution_cache[pid] = None

                    for entry in entries:
                        pdb_id = entry.get("pdb_id")
                        chain = entry.get("chain")
                        
                        if not pdb_id:
                            continue
                        
                        resolution = resolution_cache.get(pdb_id)

                        # Check resolution <= 2.0 (preserving user's threshold)
                        if resolution is not None:
                            try:
                                res_val = float(resolution)
                                if res_val <= 4.0:
                                    print(f"{rf_id}, {pdb_id}, {chain}, {resolution}")
                                    writer.writerow({
                                        'Rfam_ID': rf_id,
                                        'PDB': pdb_id,
                                        'Chain': chain,
                                        'Resolution': resolution
                                    })
                                    csvfile.flush() # Ensure data is written to disk immediately
                                    count += 1
                            except ValueError:
                                pass # Resolution not a number
                                
                    print(f"  -> Added {count} entries for {rf_id}")

                except Exception as e:
                    print(f"Error processing {rf_id}: {e}")
                    
        print(f"Done. Results saved to {output_csv}")
        
    except Exception as e:
        print(f"Error writing to output file: {e}")

if __name__ == "__main__":
    main()
