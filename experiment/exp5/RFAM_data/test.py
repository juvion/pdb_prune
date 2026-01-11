import sys
import os
import csv

# Add project root to path to import utils
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../"))
if project_root not in sys.path:
    sys.path.append(project_root)

from utils.rfam_api import get_rfam_structures
from utils.rcsb_api import fetch_entry_summary

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

    for rf_id in rfam_ids:
        print(f"Processing {rf_id}...")


if __name__ == "__main__":
    main()
