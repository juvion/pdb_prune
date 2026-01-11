import requests
import json
from typing import Any, Dict, List

def parse_search_response(data: Any) -> List[str]:
    if isinstance(data, dict):
        rs = data.get("result_set")
        if isinstance(rs, list):
            out = []
            for e in rs:
                if isinstance(e, dict):
                    v = e.get("identifier")
                    if v:
                        out.append(str(v).lower())
            return out
        dl = data.get("data")
        if isinstance(dl, list):
            out = []
            for e in dl:
                if isinstance(e, dict):
                    v = e.get("identifier") or e.get("pdb_id") or e.get("id")
                    if v:
                        out.append(str(v).lower())
                else:
                    out.append(str(e).lower())
            return out
        return []
    if isinstance(data, list):
        out = []
        for e in data:
            if isinstance(e, dict):
                v = e.get("identifier") or e.get("pdb_id") or e.get("id")
                if v:
                    out.append(str(v).lower())
            else:
                out.append(str(e).lower())
        return out
    if isinstance(data, str):
        try:
            return parse_search_response(json.loads(data))
        except Exception:
            return [data.lower()]
    return []

def fetch_entry_summary(pdb_id: str, include_chains: bool = True) -> Dict[str, Any]:
    entry_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    r = requests.get(entry_url, timeout=20)
    resolution = None
    chains: List[str] = []
    if r.status_code == 200:
        j = r.json()
        res = j.get("rcsb_entry_info", {}).get("resolution_combined")
        if isinstance(res, list) and res:
            resolution = res[0]
        
        if include_chains:
            eids = j.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids") or []
            for eid in eids:
                try:
                    e_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{eid}"
                    er = requests.get(e_url, timeout=20)
                    if er.status_code == 200:
                        ej = er.json()
                        typ = ej.get("entity_poly", {}).get("rcsb_entity_polymer_type")
                        if typ == "RNA":
                            ids_list = ej.get("rcsb_polymer_entity_container_identifiers", {}).get("auth_asym_ids") or []
                            chains.extend([str(c) for c in ids_list])
                except Exception:
                    continue
    return {"pdb_id": pdb_id.lower(), "chains": sorted(set(chains)), "resolution": resolution}

def fetch_batch_resolutions(pdb_ids: List[str]) -> Dict[str, float]:
    """
    Fetches resolutions for multiple PDB IDs using RCSB GraphQL API.
    Returns a dictionary mapping PDB ID (lowercase) to resolution (float or None).
    """
    if not pdb_ids:
        return {}
        
    # GraphQL endpoint
    url = "https://data.rcsb.org/graphql"
    
    # Query template
    query = """
    query($ids: [String!]!) {
      entries(entry_ids: $ids) {
        rcsb_id
        rcsb_entry_info {
          resolution_combined
        }
      }
    }
    """
    
    results = {}
    
    # Process in chunks of 50 to avoid request limits
    chunk_size = 50
    for i in range(0, len(pdb_ids), chunk_size):
        chunk = pdb_ids[i:i + chunk_size]
        variables = {"ids": [pid.upper() for pid in chunk]}
        
        try:
            r = requests.post(url, json={"query": query, "variables": variables}, timeout=30)
            if r.status_code == 200:
                data = r.json()
                entries = data.get("data", {}).get("entries") or []
                for entry in entries:
                    pid = entry.get("rcsb_id", "").lower()
                    res_list = entry.get("rcsb_entry_info", {}).get("resolution_combined")
                    resolution = res_list[0] if isinstance(res_list, list) and res_list else None
                    if pid:
                        results[pid] = resolution
        except Exception:
            # On error, we just skip this chunk, or we could log it.
            # For robustness, we could try fallback or just leave as None
            continue
            
    return results
