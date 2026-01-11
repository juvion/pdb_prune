import argparse
import csv
import os
import sys
def read_rfam_ids(csv_path):
    ids = []
    with open(csv_path, newline="") as f:
        r = csv.reader(f)
        rows = list(r)
        if not rows:
            return ids
        start_idx = 1 if rows[0] and rows[0][0].strip().lower() in {"rfam_id", "rfam_acc"} else 0
        for row in rows[start_idx:]:
            if not row:
                continue
            v = row[0].strip()
            if v:
                ids.append(v)
    return ids
def filter_regions(txt_path, rfam_ids):
    rfset = set(rfam_ids)
    matched = []
    per_id = {rid: [] for rid in rfam_ids}
    with open(txt_path, "r") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            parts = s.split()
            rid = parts[0] if parts else ""
            if rid in rfset:
                matched.append(line)
                per_id[rid].append(line)
    return matched, per_id
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", default="/Users/ju/Documents/Dev/pdb_prune/experiment/exp5/RFAM_data/RNA_families.csv")
    parser.add_argument("--txt", default="/Users/ju/Documents/Dev/pdb_prune/experiment/exp5/RFAM_data/pdb_full_region.txt")
    parser.add_argument("--split", action="store_true")
    args = parser.parse_args()
    rfam_ids = read_rfam_ids(args.csv)
    matched, per_id = filter_regions(args.txt, rfam_ids)
    out_dir = os.path.dirname(args.txt)
    combined_out = os.path.join(out_dir, "pdb_full_region.filtered.txt")
    with open(combined_out, "w") as w:
        w.writelines(matched)
    print(f"Filtered {len(matched)} lines written to {combined_out}")
    if args.split:
        for rid, lines in per_id.items():
            if not lines:
                continue
            out_path = os.path.join(out_dir, f"pdb_full_region_{rid}.txt")
            with open(out_path, "w") as w:
                w.writelines(lines)
        print("Per-ID files written")
if __name__ == "__main__":
    main()
