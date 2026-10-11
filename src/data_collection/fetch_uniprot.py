"""
Step 1 (part B): Pull canonical protein sequence + domain/feature annotations
from UniProt for each gene in the panel.

Uses the UniProt REST API (https://rest.uniprot.org) — no authentication
required. Same sandbox caveat as fetch_clinvar.py: this domain is not
reachable from the authoring environment; run from your own machine.

Usage:
    python -m src.data_collection.fetch_uniprot
"""
import json
import time
import requests
from tqdm import tqdm

from src.config import RAW_DIR, load_gene_panel, ensure_dirs

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"
REQUEST_DELAY = 0.2

# Positional UniProt annotations kept for feature engineering. Broad
# "Region" / "Compositional bias" spans are kept too, but Step 3 treats
# them separately from structured domains and functional sites.
UNIPROT_FEATURE_TYPES = (
    "Domain",
    "Region",
    "Repeat",
    "Zinc finger",
    "Motif",
    "Coiled coil",
    "Compositional bias",
    "Active site",
    "Binding site",
    "Site",
    "Disulfide bond",
    "Transmembrane",
    "Intramembrane",
    "Topological domain",
)


def fetch_entry(uniprot_id: str) -> dict:
    """Fetch the full UniProtKB JSON entry for one accession."""
    resp = requests.get(f"{UNIPROT_BASE}/{uniprot_id}.json", timeout=30)
    resp.raise_for_status()
    return resp.json()


def extract_relevant_fields(entry: dict, gene_symbol: str) -> dict:
    """Pull out sequence + domain/feature annotations relevant to Step 3 feature engineering."""
    sequence = entry.get("sequence", {}).get("value", "")
    features = []
    for feat in entry.get("features", []):
        # Domains, active sites, binding sites, disulfide bonds, etc. — anything
        # with a positional span is potentially useful for the "domain
        # membership" feature described in the methodology.
        if feat.get("type") in UNIPROT_FEATURE_TYPES:
            loc = feat.get("location", {})
            start = loc.get("start", {}).get("value")
            end = loc.get("end", {}).get("value")
            features.append(
                {
                    "type": feat.get("type"),
                    "description": feat.get("description"),
                    "start": start,
                    "end": end,
                }
            )

    return {
        "gene_symbol": gene_symbol,
        "uniprot_id": entry.get("primaryAccession"),
        "sequence": sequence,
        "sequence_length": len(sequence),
        "features": features,
    }


def main():
    ensure_dirs()
    genes = load_gene_panel("both")
    all_entries = []
    for gene in tqdm(genes, desc="Fetching UniProt entries"):
        try:
            entry = fetch_entry(gene["uniprot_id"])
            all_entries.append(extract_relevant_fields(entry, gene["symbol"]))
        except requests.RequestException as e:
            print(f"[warn] failed to fetch {gene['symbol']} ({gene['uniprot_id']}): {e}")
        time.sleep(REQUEST_DELAY)

    out_path = RAW_DIR / "uniprot_entries.json"
    with open(out_path, "w") as f:
        json.dump(all_entries, f, indent=2)
    print(f"Saved {len(all_entries)} UniProt entries to {out_path}")


if __name__ == "__main__":
    main()
