"""
Step 1 (part A): Pull ClinVar variant records for each gene in the panel.

Uses NCBI E-utilities (esearch -> esummary) against the ClinVar database.
esearch is paged with retstart so that every missense SNV for a gene is
retrieved. The earlier version stopped at the first 500 IDs, which left
genes such as MYH7 (3,000+ missense records) with only a few dozen labelled
variants.

Set NCBI_EMAIL (and optionally NCBI_API_KEY) in the environment before
running; see src/config.py.

Usage:
    python -m src.data_collection.fetch_clinvar
"""
import time

import pandas as pd
import requests
from tqdm import tqdm

from src.config import RAW_DIR, load_gene_panel, NCBI_CONTACT_EMAIL, NCBI_API_KEY, USER_AGENT, ensure_dirs

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
# NCBI allows 3 requests/sec without an API key and 10/sec with one.
REQUEST_DELAY = 0.12 if NCBI_API_KEY else 0.36
ESEARCH_PAGE = 5000
ESUMMARY_BATCH = 200
MAX_RETRIES = 5


def _eutils_params(extra: dict) -> dict:
    params = {"retmode": "json", "email": NCBI_CONTACT_EMAIL, "tool": "cvd-xai", **extra}
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    return params


def _get(endpoint: str, params: dict) -> dict:
    """GET an E-utilities endpoint, backing off on rate limits and 5xx errors."""
    for attempt in range(MAX_RETRIES):
        try:
            resp = requests.get(
                f"{EUTILS_BASE}/{endpoint}",
                params=_eutils_params(params),
                headers={"User-Agent": USER_AGENT},
                timeout=60,
            )
            if resp.status_code == 429 or resp.status_code >= 500:
                raise requests.HTTPError(f"HTTP {resp.status_code}", response=resp)
            resp.raise_for_status()
            return resp.json()
        except (requests.RequestException, ValueError):
            if attempt == MAX_RETRIES - 1:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def search_clinvar_ids(gene_symbol: str) -> list[str]:
    """Return every ClinVar UID for missense single-nucleotide variants in a gene."""
    query = (
        f'{gene_symbol}[gene] AND "single nucleotide variant"[Type] '
        f'AND "missense variant"[molecular consequence]'
    )
    uids: list[str] = []
    retstart = 0
    while True:
        result = _get(
            "esearch.fcgi",
            {"db": "clinvar", "term": query, "retmax": ESEARCH_PAGE, "retstart": retstart},
        )["esearchresult"]
        page = result.get("idlist", [])
        uids.extend(page)
        total = int(result.get("count", 0))
        retstart += len(page)
        time.sleep(REQUEST_DELAY)
        if not page or retstart >= total:
            break
    return uids


def _summary_record(uid: str, rec: dict) -> dict:
    classification = rec.get("germline_classification") or {}
    genes = rec.get("genes") or []
    return {
        "clinvar_uid": uid,
        "title": rec.get("title"),
        "gene_sort": rec.get("gene_sort"),
        "clinical_significance": classification.get("description"),
        "review_status": classification.get("review_status"),
        "variation_set": rec.get("variation_set"),
        "protein_change": rec.get("protein_change"),
        "molecular_consequence_list": rec.get("molecular_consequence_list"),
        "clinvar_genes": [g.get("symbol") for g in genes if isinstance(g, dict)],
    }


def fetch_clinvar_summaries(uids: list[str]) -> list[dict]:
    """Fetch ClinVar esummary records for up to ESUMMARY_BATCH UIDs."""
    if not uids:
        return []
    result = _get("esummary.fcgi", {"db": "clinvar", "id": ",".join(uids)}).get("result", {})
    return [_summary_record(uid, result.get(uid, {})) for uid in result.get("uids", [])]


def fetch_gene_variants(gene_symbol: str) -> pd.DataFrame:
    uids = search_clinvar_ids(gene_symbol)
    records = []
    for i in tqdm(range(0, len(uids), ESUMMARY_BATCH), desc=gene_symbol, leave=False):
        records.extend(fetch_clinvar_summaries(uids[i : i + ESUMMARY_BATCH]))
        time.sleep(REQUEST_DELAY)
    df = pd.DataFrame(records)
    if not df.empty:
        df["gene_symbol"] = gene_symbol
    return df


def main():
    ensure_dirs()
    genes = load_gene_panel("both")
    all_records = []
    for gene in tqdm(genes, desc="Fetching ClinVar records"):
        try:
            df = fetch_gene_variants(gene["symbol"])
        except requests.RequestException as e:
            print(f"[warn] failed to fetch {gene['symbol']}: {e}")
            continue
        df["panel"] = gene["panel"]
        print(f"  {gene['symbol']}: {len(df)} records")
        all_records.append(df)

    combined = pd.concat(all_records, ignore_index=True) if all_records else pd.DataFrame()

    out_path = RAW_DIR / "clinvar_raw.csv"
    combined.to_csv(out_path, index=False)
    print(f"Saved {len(combined)} ClinVar records to {out_path}")


if __name__ == "__main__":
    main()
