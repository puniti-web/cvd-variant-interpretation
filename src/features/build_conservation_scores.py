"""
Build protein-level conservation scores for the CVD gene panel.

For each human UniProt reference protein:
    1. Retrieve vertebrate orthologs from UniProtKB: entries with the same
       gene name, excluding Homo sapiens, one sequence per species
       (reviewed entry preferred, then the longest). Swiss-Prot alone has
       only a handful of orthologs for most of these genes.
    2. Globally align each ortholog to the human sequence (BLOSUM62,
       affine gaps).
    3. For every human residue, record which amino acid each ortholog has
       at the aligned position.

Outputs:
    data/raw/conservation_scores.json
        {gene: {position: fraction of aligned orthologs that keep the
        human residue}}. Same format as before; the API reads it.
    data/raw/conservation_profiles.json
        {gene: {position: string of ortholog residues at that position}},
        used to ask whether the mutant residue already occurs in another
        species (a strong sign the substitution is tolerated).

This is a small-ortholog-set proxy, complementary to the genome-level
phyloP / GERP scores from Ensembl VEP. Positions with fewer than
MIN_HOMOLOGS aligned orthologs get no score.

Usage:
    python -m src.features.build_conservation_scores
"""
import json
import time

import requests
from Bio.Align import PairwiseAligner, substitution_matrices

from src.config import RAW_DIR, USER_AGENT, ensure_dirs

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
CONSERVATION_OUTPUT = RAW_DIR / "conservation_scores.json"
PROFILE_OUTPUT = RAW_DIR / "conservation_profiles.json"
REQUEST_TIMEOUT = 60
MIN_HOMOLOGS = 3
MAX_HOMOLOGS = 60
SEARCH_PAGE_SIZE = 500
VERTEBRATA_TAXON = 7742
# Orthologs that align over less than this fraction of the human sequence
# are usually fragments or mis-annotated paralogs.
MIN_ALIGNED_FRACTION = 0.5


def load_uniprot_entries():
    path = RAW_DIR / "uniprot_entries.json"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found. Run `python -m src.data_collection.fetch_uniprot` first.")
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def search_homologs(gene_symbol: str) -> list[dict]:
    """One UniProt sequence per non-human vertebrate species with the same gene name."""
    params = {
        "query": f"gene_exact:{gene_symbol} AND taxonomy_id:{VERTEBRATA_TAXON} NOT organism_id:9606 "
                 f"AND fragment:false",
        "format": "json",
        "size": SEARCH_PAGE_SIZE,
        "fields": "accession,reviewed,organism_name,organism_id,sequence",
    }
    resp = requests.get(UNIPROT_SEARCH_URL, params=params, headers={"User-Agent": USER_AGENT}, timeout=REQUEST_TIMEOUT)
    resp.raise_for_status()
    by_species: dict[int, dict] = {}
    for result in resp.json().get("results", []):
        sequence = (result.get("sequence") or {}).get("value")
        organism = result.get("organism") or {}
        if not sequence:
            continue
        candidate = {
            "accession": result.get("primaryAccession"),
            "organism": organism.get("scientificName"),
            "reviewed": "reviewed" in str(result.get("entryType", "")).lower()
            and "unreviewed" not in str(result.get("entryType", "")).lower(),
            "sequence": sequence,
        }
        current = by_species.get(organism.get("taxonId"))
        rank = (candidate["reviewed"], len(sequence))
        if current is None or rank > (current["reviewed"], len(current["sequence"])):
            by_species[organism.get("taxonId")] = candidate
    homologs = sorted(by_species.values(), key=lambda h: (not h["reviewed"], -len(h["sequence"])))
    return homologs[:MAX_HOMOLOGS]


def create_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    return aligner


def aligned_residues(aligner: PairwiseAligner, reference: str, homolog: str) -> dict[int, str]:
    """Map 1-based reference positions to the aligned homolog residue (gaps omitted)."""
    alignment = aligner.align(reference, homolog)[0]
    mapping = {}
    for (r_start, r_end), (h_start, _) in zip(*alignment.aligned):
        for k in range(r_end - r_start):
            mapping[r_start + k + 1] = homolog[h_start + k]
    return mapping


def build_profiles(reference: str, homologs: list[dict]) -> dict[int, str]:
    aligner = create_aligner()
    profiles: dict[int, list[str]] = {pos: [] for pos in range(1, len(reference) + 1)}
    for homolog in homologs:
        mapping = aligned_residues(aligner, reference, homolog["sequence"])
        if len(mapping) < MIN_ALIGNED_FRACTION * len(reference):
            continue
        for pos, aa in mapping.items():
            profiles[pos].append(aa)
    return {pos: "".join(aas) for pos, aas in profiles.items()}


def identity_scores(reference: str, profiles: dict[int, str]) -> dict[str, float]:
    return {
        str(pos): round(aas.count(reference[pos - 1]) / len(aas), 4)
        for pos, aas in profiles.items()
        if len(aas) >= MIN_HOMOLOGS
    }


def main():
    ensure_dirs()
    scores, profiles = {}, {}
    for entry in load_uniprot_entries():
        gene, reference = entry["gene_symbol"], entry["sequence"]
        try:
            homologs = search_homologs(gene)
            gene_profiles = build_profiles(reference, homologs) if len(homologs) >= MIN_HOMOLOGS else {}
        except requests.RequestException as exc:
            print(f"[warn] {gene}: homolog search failed: {exc}")
            gene_profiles = {}
        scores[gene] = identity_scores(reference, gene_profiles)
        profiles[gene] = {str(p): aas for p, aas in gene_profiles.items() if aas}
        depth = max((len(a) for a in gene_profiles.values()), default=0)
        print(f"  {gene}: {len(homologs) if gene_profiles else 0} orthologs, max depth {depth}, "
              f"{len(scores[gene])} scored residues")
        time.sleep(0.5)

    with open(CONSERVATION_OUTPUT, "w", encoding="utf-8") as f:
        json.dump(scores, f)
    with open(PROFILE_OUTPUT, "w", encoding="utf-8") as f:
        json.dump(profiles, f)
    print(f"Saved {CONSERVATION_OUTPUT} and {PROFILE_OUTPUT}")


if __name__ == "__main__":
    main()
