import requests
from typing import Any, Dict, List, Optional
from .rcsb_api import fetch_entry_summary

def get_rfam_structures(rfam_acc: str) -> List[Dict[str, Any]]:
    url = f"https://rfam.org/family/{rfam_acc}/structures?content-type=application/json"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    j = r.json()
    mappings = j.get("mapping") or []
    out: List[Dict[str, Any]] = []
    for m in mappings:
        pdb_id = str(m.get("pdb_id", "")).lower()
        chain = str(m.get("chain", "")).strip()
        out.append({
            "rfam_acc": rfam_acc,
            "pdb_id": pdb_id,
            "chain": chain,
            "bit_score": m.get("bit_score"),
            "evalue_score": m.get("evalue_score"),
            "pdb_start": m.get("pdb_start"),
            "pdb_end": m.get("pdb_end"),
            "cm_start": m.get("cm_start"),
            "cm_end": m.get("cm_end"),
        })
    return out

def get_rfam_structures_with_resolution(rfam_acc: str) -> List[Dict[str, Any]]:
    entries = get_rfam_structures(rfam_acc)
    # fetch resolution once per pdb_id to avoid repeated calls
    cache: Dict[str, Optional[float]] = {}
    for e in entries:
        pid = e["pdb_id"]
        if pid not in cache:
            summary = fetch_entry_summary(pid)
            cache[pid] = summary.get("resolution")
        e["resolution"] = cache[pid]
    return entries

