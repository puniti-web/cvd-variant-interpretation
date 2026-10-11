"""
Step 2: Data Preprocessing.

- Parses the amino-acid substitution from the ClinVar title, which carries
  the preferred (MANE Select) transcript, e.g.
  "NM_000257.4(MYH7):c.1208G>A (p.Arg403Gln)". ClinVar's protein_change
  field lists one change per isoform ("R69C, R147C, R105C") in no fixed
  order, so it is only used as a fallback.
- Validates every substitution against the UniProt canonical sequence:
  the wild-type residue must match sequence[position - 1]. Variants that
  match no candidate are dropped, so every downstream feature (conservation,
  domains, structure) indexes the right residue.
- Maps ClinVar review status to the star rating and keeps only variants
  with at least one star in the training set.
- Collapses multiple ClinVar records for the same protein change (different
  nucleotide changes, separate submissions). If those records disagree
  (pathogenic vs benign), the variant is excluded rather than guessed.
- Writes the labelled training set, the VUS set held out for qualitative
  review, and a QC report with counts at every filter step.

Usage:
    python -m src.preprocessing.clean_variants
"""
import ast
import json
import re

import pandas as pd

from src.config import RAW_DIR, PROCESSED_DIR, ensure_dirs, load_gene_panel

PATHOGENIC_LABELS = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
BENIGN_LABELS = {"Benign", "Likely benign", "Benign/Likely benign"}
VUS_LABELS = {"Uncertain significance"}
CONFLICT_LABELS = {
    "Conflicting classifications of pathogenicity",
    "Conflicting interpretations of pathogenicity",
}

# ClinVar review status -> star rating
# https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/
REVIEW_STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, single submitter": 1,
}
MIN_TRAINING_STARS = 1

# 3-letter -> 1-letter amino acid code map, used to parse ClinVar's
# p.Xxx###Yyy notation into a compact form.
AA_3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*",
}
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

PROTEIN_CHANGE_RE = re.compile(
    r"(?:p\.)?(?:(?P<wt3>[A-Za-z]{3})(?P<pos>\d+)(?P<mut3>[A-Za-z]{3})|(?P<wt1>[A-Z])(?P<pos1>\d+)(?P<mut1>[A-Z]))"
)
TITLE_RE = re.compile(
    r"\((?P<gene>[A-Za-z0-9-]+)\):c\.\S+\s+\(p\.(?P<wt>[A-Z][a-z]{2})(?P<pos>\d+)(?P<mut>[A-Z][a-z]{2})\)"
)
ONE_LETTER_RE = re.compile(r"\b([A-Z])(\d+)([A-Z])\b")


def parse_protein_change(text: str):
    """Extract (wt_aa, position, mut_aa) from ClinVar protein-change strings.

    Supports both:
        p.Arg403Gln
        Arg403Gln
        R403Q
    """
    if not isinstance(text, str):
        return None, None, None

    m = PROTEIN_CHANGE_RE.search(text)
    if not m:
        return None, None, None

    # Three-letter amino-acid notation
    if m.group("wt3"):
        wt1 = AA_3TO1.get(m.group("wt3"))
        mut1 = AA_3TO1.get(m.group("mut3"))
        pos = int(m.group("pos"))
        return wt1, pos, mut1

    # One-letter amino-acid notation
    wt1 = m.group("wt1")
    mut1 = m.group("mut1")
    pos = int(m.group("pos1"))

    return wt1, pos, mut1


def candidate_changes(title, protein_change, gene_symbol: str) -> list[tuple[str, int, str, str]]:
    """
    All (wt, pos, mut, source) substitutions a ClinVar record could refer
    to, preferred transcript first. Records whose title names a different
    gene (overlapping loci returned by the [gene] search) yield nothing.
    """
    candidates = []
    if isinstance(title, str):
        m = TITLE_RE.search(title)
        if m:
            if m.group("gene").upper() != gene_symbol.upper():
                return []
            wt, mut = AA_3TO1.get(m.group("wt")), AA_3TO1.get(m.group("mut"))
            if wt and mut:
                candidates.append((wt, int(m.group("pos")), mut, "title"))
    if isinstance(protein_change, str):
        for wt, pos, mut in ONE_LETTER_RE.findall(protein_change):
            cand = (wt, int(pos), mut, "alt_isoform")
            if cand[:3] not in {c[:3] for c in candidates}:
                candidates.append(cand)
    return candidates


def resolve_against_sequence(candidates, sequence: str | None):
    """Return the first candidate whose wild-type residue matches the canonical sequence."""
    if not sequence:
        return None
    for wt, pos, mut, source in candidates:
        if 1 <= pos <= len(sequence) and sequence[pos - 1] == wt:
            return wt, pos, mut, source
    return None


def normalize_label(raw_significance: str) -> str:
    if not isinstance(raw_significance, str):
        return "unknown"
    sig = raw_significance.strip()
    if sig in PATHOGENIC_LABELS:
        return "pathogenic"
    if sig in BENIGN_LABELS:
        return "benign"
    if sig in VUS_LABELS:
        return "vus"
    if sig in CONFLICT_LABELS:
        return "conflicting"
    return "other"  # risk factor, drug response, not provided, etc.


def review_stars(status) -> int:
    if not isinstance(status, str):
        return 0
    return REVIEW_STARS.get(status.strip().lower(), 0)


def parse_variation_ids(variation_set) -> dict:
    """Pull cDNA change, dbSNP rsID and canonical SPDI from ClinVar's variation_set."""
    out = {"cdna_change": None, "rsid": None, "spdi": None}
    if isinstance(variation_set, str):
        try:
            variation_set = ast.literal_eval(variation_set)
        except (ValueError, SyntaxError):
            return out
    if not isinstance(variation_set, list) or not variation_set:
        return out
    measure = variation_set[0] if isinstance(variation_set[0], dict) else {}
    out["cdna_change"] = measure.get("cdna_change") or None
    out["spdi"] = measure.get("canonical_spdi") or None
    for xref in measure.get("variation_xrefs") or []:
        if isinstance(xref, dict) and xref.get("db_source") == "dbSNP" and xref.get("db_id"):
            out["rsid"] = f"rs{xref['db_id']}"
            break
    return out


def load_canonical_sequences() -> dict[str, str]:
    path = RAW_DIR / "uniprot_entries.json"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found — run `python -m src.data_collection.fetch_uniprot` first."
        )
    with open(path, encoding="utf-8") as f:
        return {e["gene_symbol"]: e["sequence"] for e in json.load(f)}


def parse_and_validate(df: pd.DataFrame, sequences: dict[str, str], qc: dict) -> pd.DataFrame:
    rows = []
    reasons = {"other_gene_or_unparsed": 0, "sequence_mismatch": 0, "nonsense_or_synonymous": 0, "start_loss": 0}
    for rec in df.to_dict("records"):
        candidates = candidate_changes(rec.get("title"), rec.get("protein_change"), rec["gene_symbol"])
        if not candidates:
            reasons["other_gene_or_unparsed"] += 1
            continue
        resolved = resolve_against_sequence(candidates, sequences.get(rec["gene_symbol"]))
        if resolved is None:
            reasons["sequence_mismatch"] += 1
            continue
        wt, pos, mut, source = resolved
        if wt not in STANDARD_AA or mut not in STANDARD_AA or wt == mut:
            reasons["nonsense_or_synonymous"] += 1
            continue
        if pos == 1 and wt == "M":
            # p.Met1X removes the start codon: a loss-of-initiation variant, not a missense change.
            reasons["start_loss"] += 1
            continue
        rows.append({**rec, "wt_aa": wt, "position": pos, "mut_aa": mut, "position_source": source})
    qc["dropped_at_parsing"] = reasons
    return pd.DataFrame(rows)


def collapse_duplicates(df: pd.DataFrame, qc: dict) -> pd.DataFrame:
    """
    One row per (gene, wt, position, mut). The highest-star record wins,
    preferring a definitive classification over VUS at equal stars. If
    one-star-or-better records both call it pathogenic and benign, the
    variant is marked as conflicting.
    """
    keys = ["gene_symbol", "wt_aa", "position", "mut_aa"]
    df = df.copy()
    df["_definitive"] = df["label"].isin(["pathogenic", "benign"]).astype(int)
    df = df.sort_values(["review_stars", "_definitive", "clinvar_uid"], ascending=[False, False, True])

    confident = df[(df["review_stars"] >= MIN_TRAINING_STARS) & (df["_definitive"] == 1)]
    labels_per_variant = confident.groupby(keys)["label"].nunique()
    conflicted = set(labels_per_variant[labels_per_variant > 1].index)

    counts = df.groupby(keys).size().rename("n_clinvar_records")
    best = df.drop_duplicates(subset=keys, keep="first").set_index(keys)
    best = best.join(counts)
    best["label_conflict"] = [idx in conflicted for idx in best.index]
    best = best.reset_index().drop(columns="_definitive")

    qc["duplicate_records_collapsed"] = int(len(df) - len(best))
    qc["variants_with_label_conflict"] = len(conflicted)
    return best


def clean(df: pd.DataFrame, sequences: dict[str, str], qc: dict | None = None) -> pd.DataFrame:
    qc = {} if qc is None else qc
    qc["raw_records"] = int(len(df))

    df = parse_and_validate(df, sequences, qc)
    qc["validated_records"] = int(len(df))

    df["label"] = df["clinical_significance"].map(normalize_label)
    df["review_stars"] = df["review_status"].map(review_stars)
    ids = df["variation_set"].map(parse_variation_ids).apply(pd.Series)
    df = pd.concat([df, ids], axis=1)

    df = collapse_duplicates(df, qc)
    df["position"] = df["position"].astype(int)
    qc["unique_variants"] = int(len(df))
    qc["unique_variants_by_label"] = df["label"].value_counts().to_dict()
    return df.reset_index(drop=True)


def split_sets(cleaned: pd.DataFrame, qc: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    definitive = cleaned["label"].isin(["pathogenic", "benign"])
    enough_stars = cleaned["review_stars"] >= MIN_TRAINING_STARS
    labelled = cleaned[definitive & enough_stars & ~cleaned["label_conflict"]]
    vus = cleaned[cleaned["label"] == "vus"]

    qc["excluded_from_training"] = {
        "zero_star_definitive": int((definitive & ~enough_stars).sum()),
        "label_conflict": int((definitive & enough_stars & cleaned["label_conflict"]).sum()),
    }
    qc["labelled"] = {
        "total": int(len(labelled)),
        "pathogenic": int((labelled["label"] == "pathogenic").sum()),
        "benign": int((labelled["label"] == "benign").sum()),
    }
    qc["vus"] = int(len(vus))
    return labelled.reset_index(drop=True), vus.reset_index(drop=True)


def per_gene_table(labelled: pd.DataFrame, vus: pd.DataFrame) -> pd.DataFrame:
    panel = {g["symbol"]: g["panel"] for g in load_gene_panel("both")}
    table = (
        labelled.groupby(["gene_symbol", "label"]).size().unstack(fill_value=0)
        .reindex(columns=["pathogenic", "benign"], fill_value=0)
    )
    table = table.reindex(sorted(panel), fill_value=0)
    table["labelled_total"] = table["pathogenic"] + table["benign"]
    table["vus"] = vus["gene_symbol"].value_counts().reindex(table.index, fill_value=0)
    table.insert(0, "panel", [panel[g] for g in table.index])
    return table.rename_axis("gene_symbol").reset_index()


OUTPUT_COLUMNS = [
    "clinvar_uid", "title", "gene_sort", "clinical_significance", "review_status",
    "variation_set", "protein_change", "molecular_consequence_list", "gene_symbol",
    "panel", "wt_aa", "position", "mut_aa", "label",
    "review_stars", "position_source", "n_clinvar_records", "cdna_change", "rsid", "spdi",
]


def main():
    ensure_dirs()
    raw_path = RAW_DIR / "clinvar_raw.csv"
    if not raw_path.exists():
        raise FileNotFoundError(
            f"{raw_path} not found — run `python -m src.data_collection.fetch_clinvar` first."
        )

    qc: dict = {}
    df = pd.read_csv(raw_path)
    cleaned = clean(df, load_canonical_sequences(), qc)
    labelled, vus = split_sets(cleaned, qc)

    labelled[OUTPUT_COLUMNS].to_csv(PROCESSED_DIR / "variants_labelled.csv", index=False)
    vus[OUTPUT_COLUMNS].to_csv(PROCESSED_DIR / "variants_vus.csv", index=False)

    genes = per_gene_table(labelled, vus)
    genes.to_csv(PROCESSED_DIR / "qc_per_gene.csv", index=False)
    with open(PROCESSED_DIR / "qc_report.json", "w", encoding="utf-8") as f:
        json.dump(qc, f, indent=2)

    print(json.dumps(qc, indent=2))
    print(genes.to_string(index=False))


if __name__ == "__main__":
    main()
