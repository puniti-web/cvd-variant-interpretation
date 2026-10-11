"""Offline tests for ClinVar parsing, sequence validation and de-duplication."""
import pandas as pd

from src.preprocessing.clean_variants import (
    candidate_changes,
    clean,
    review_stars,
    resolve_against_sequence,
    split_sets,
)

# Toy canonical sequence: position 3 is R, position 5 is G.
SEQ = "MARKGLV"


def test_title_is_preferred_over_protein_change_list():
    cands = candidate_changes("NM_1.1(GENE1):c.8G>A (p.Arg3Gln)", "R9Q, R3Q", "GENE1")
    assert cands[0] == ("R", 3, "Q", "title")
    assert ("R", 9, "Q", "alt_isoform") in cands


def test_title_for_other_gene_is_rejected():
    assert candidate_changes("NM_1.1(OTHER):c.8G>A (p.Arg3Gln)", "R3Q", "GENE1") == []


def test_falls_back_to_isoform_that_matches_sequence():
    # Title position does not match the canonical sequence; the isoform list does.
    cands = candidate_changes("NM_1.1(GENE1):c.8G>A (p.Arg9Gln)", "R9Q, R3Q", "GENE1")
    assert resolve_against_sequence(cands, SEQ) == ("R", 3, "Q", "alt_isoform")


def test_no_candidate_matches_sequence():
    cands = candidate_changes("NM_1.1(GENE1):c.8G>A (p.Trp3Gln)", None, "GENE1")
    assert resolve_against_sequence(cands, SEQ) is None


def test_review_stars():
    assert review_stars("reviewed by expert panel") == 3
    assert review_stars("criteria provided, single submitter") == 1
    assert review_stars("no assertion criteria provided") == 0
    assert review_stars(None) == 0


def _record(uid, title, significance, status):
    return {
        "clinvar_uid": uid, "title": title, "gene_sort": "GENE1",
        "clinical_significance": significance, "review_status": status,
        "variation_set": None, "protein_change": None, "molecular_consequence_list": None,
        "gene_symbol": "GENE1", "panel": "cardiomyopathy",
    }


def test_duplicates_collapse_to_highest_star_and_conflicts_are_excluded():
    one_star = "criteria provided, single submitter"
    df = pd.DataFrame([
        # Same protein change, two records: expert panel pathogenic beats 1-star VUS.
        _record(1, "NM_1.1(GENE1):c.8G>A (p.Arg3Gln)", "Uncertain significance", one_star),
        _record(2, "NM_1.1(GENE1):c.8G>C (p.Arg3Gln)", "Pathogenic", "reviewed by expert panel"),
        # Conflicting confident labels for G5D -> excluded from training.
        _record(3, "NM_1.1(GENE1):c.14G>A (p.Gly5Asp)", "Pathogenic", one_star),
        _record(4, "NM_1.1(GENE1):c.14G>A (p.Gly5Asp)", "Benign", one_star),
        # Zero-star benign -> excluded from training.
        _record(5, "NM_1.1(GENE1):c.20T>C (p.Val7Ala)", "Benign", "no assertion criteria provided"),
        # Synonymous -> dropped.
        _record(6, "NM_1.1(GENE1):c.20T>C (p.Val7Val)", "Benign", one_star),
        # Start loss -> dropped.
        _record(7, "NM_1.1(GENE1):c.2T>C (p.Met1Thr)", "Benign", one_star),
    ])
    qc = {}
    cleaned = clean(df, {"GENE1": SEQ}, qc)
    labelled, vus = split_sets(cleaned, qc)

    assert len(cleaned) == 3
    assert list(labelled["position"]) == [3]
    assert labelled.iloc[0]["label"] == "pathogenic"
    assert labelled.iloc[0]["n_clinvar_records"] == 2
    assert qc["excluded_from_training"] == {"zero_star_definitive": 1, "label_conflict": 1}
    assert len(vus) == 0
    assert qc["dropped_at_parsing"]["start_loss"] == 1
