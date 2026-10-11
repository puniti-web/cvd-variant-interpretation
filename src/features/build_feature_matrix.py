"""
Step 3 (part C): Assemble the feature matrices used for classification.

Joins each variant (Step 2 output) with:
  - physicochemical deltas, BLOSUM62 and Grantham distance (offline)
  - ortholog conservation (build_conservation_scores.py)
  - genome-level conservation: phyloP, phastCons, GERP (Ensembl VEP)
  - UniProt domain / site / disorder annotations (fetch_uniprot.py)
  - AlphaFold structural context: pLDDT, RSA, contacts, secondary
    structure (structure_features.py)

Published predictor scores from VEP (SIFT, PolyPhen-2, AlphaMissense,
REVEL) are carried along as columns for baseline comparison only. The
model feature lists in src/modeling/train_classifier.py do not include
them.

Outputs:
    data/features/feature_matrix.csv       labelled variants
    data/features/feature_matrix_vus.csv   VUS, for prioritisation
    data/features/feature_coverage.csv     missing-value report

Usage:
    python -m src.features.build_feature_matrix
"""
import dataclasses
import json
from functools import lru_cache

import pandas as pd

from src.config import PROCESSED_DIR, FEATURES_DIR, RAW_DIR, ensure_dirs
from src.features.physicochemical import compute_physchem_features
from src.features.conservation_and_domain import (
    compute_annotation_features,
    get_conservation_score,
    load_conservation_profiles,
    load_conservation_scores,
    mutant_frequency_in_orthologs,
)

ID_COLUMNS = ["clinvar_uid", "gene_symbol", "panel", "position", "wt_aa", "mut_aa", "label", "review_stars", "spdi"]
VEP_CONSERVATION = ["phylop", "phastcons", "gerp"]
BASELINE_SCORES = ["sift_score", "polyphen_score", "alphamissense", "revel"]
STRUCTURE_COLUMNS = ["plddt", "rsa", "contact_count", "sec_struct"]


def load_uniprot() -> dict:
    path = RAW_DIR / "uniprot_entries.json"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found — run `python -m src.data_collection.fetch_uniprot` first.")
    with open(path, encoding="utf-8") as f:
        return {e["gene_symbol"]: e for e in json.load(f)}


def load_optional_csv(path, columns: list[str]) -> pd.DataFrame:
    if path.exists():
        return pd.read_csv(path)
    print(f"[warn] {path} not found — {', '.join(columns)} will be missing.")
    return pd.DataFrame(columns=columns)


def sequence_features(variants: pd.DataFrame, uniprot: dict) -> pd.DataFrame:
    scores, profiles = load_conservation_scores(), load_conservation_profiles()
    rows = []
    for rec in variants.to_dict("records"):
        gene, pos, wt, mut = rec["gene_symbol"], int(rec["position"]), rec["wt_aa"], rec["mut_aa"]
        entry = uniprot[gene]
        rows.append({
            "relative_position": round(pos / entry["sequence_length"], 4),
            **dataclasses.asdict(compute_physchem_features(wt, mut)),
            "conservation_score": get_conservation_score(gene, pos, scores),
            "mut_aa_in_orthologs": mutant_frequency_in_orthologs(gene, pos, mut, profiles),
            **dataclasses.asdict(compute_annotation_features(pos, entry["features"])),
        })
    return pd.DataFrame(rows, index=variants.index)


def build(variants: pd.DataFrame, uniprot: dict, vep: pd.DataFrame, structure: pd.DataFrame) -> pd.DataFrame:
    df = variants[[c for c in ID_COLUMNS if c in variants.columns]].copy()
    df = pd.concat([df, sequence_features(variants, uniprot)], axis=1)

    df = df.merge(vep[["spdi"] + VEP_CONSERVATION + BASELINE_SCORES], on="spdi", how="left")
    df = df.merge(structure[["gene_symbol", "position"] + STRUCTURE_COLUMNS], on=["gene_symbol", "position"], how="left")
    df["has_structure"] = df["plddt"].notna()
    df["sec_struct"] = df["sec_struct"].fillna("unknown")
    return df


@lru_cache(maxsize=1)
def _cached_sources() -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    vep = load_optional_csv(RAW_DIR / "vep_annotations.csv", ["spdi"] + VEP_CONSERVATION + BASELINE_SCORES)
    structure = load_optional_csv(
        RAW_DIR / "alphafold_residue_features.csv", ["gene_symbol", "position"] + STRUCTURE_COLUMNS
    )
    return load_uniprot(), vep, structure


def features_for_variant(record: dict) -> dict:
    """
    Feature row for a single variant (keys gene_symbol, position, wt_aa,
    mut_aa, panel, optional spdi), computed exactly as for the training
    matrix. Used by the API for per-variant predictions.
    """
    uniprot, vep, structure = _cached_sources()
    variant = pd.DataFrame([{c: record.get(c) for c in ID_COLUMNS}])
    return build(variant, uniprot, vep, structure).iloc[0].to_dict()


def coverage_report(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.isna().mean().rename("missing_fraction").round(4)
        .to_frame().rename_axis("column").reset_index()
        .sort_values("missing_fraction", ascending=False)
    )


def main():
    ensure_dirs()
    labelled_path = PROCESSED_DIR / "variants_labelled.csv"
    if not labelled_path.exists():
        raise FileNotFoundError(f"{labelled_path} not found — run `python -m src.preprocessing.clean_variants` first.")

    uniprot = load_uniprot()
    vep = load_optional_csv(RAW_DIR / "vep_annotations.csv", ["spdi"] + VEP_CONSERVATION + BASELINE_SCORES)
    structure = load_optional_csv(
        RAW_DIR / "alphafold_residue_features.csv", ["gene_symbol", "position"] + STRUCTURE_COLUMNS
    )

    labelled = build(pd.read_csv(labelled_path), uniprot, vep, structure)
    labelled.to_csv(FEATURES_DIR / "feature_matrix.csv", index=False)
    coverage = coverage_report(labelled)
    coverage.to_csv(FEATURES_DIR / "feature_coverage.csv", index=False)
    print(f"Labelled feature matrix: {labelled.shape[0]} rows x {labelled.shape[1]} columns")
    print(coverage[coverage["missing_fraction"] > 0].to_string(index=False))

    vus_path = PROCESSED_DIR / "variants_vus.csv"
    if vus_path.exists():
        vus = build(pd.read_csv(vus_path), uniprot, vep, structure)
        vus.to_csv(FEATURES_DIR / "feature_matrix_vus.csv", index=False)
        print(f"VUS feature matrix: {vus.shape[0]} rows")


if __name__ == "__main__":
    main()
