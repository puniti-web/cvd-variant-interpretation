"""Tests for the rule-based parts of the integrated report."""
import math

import pytest

from src.report.generate_report import (
    balanced_probability,
    conclusion,
    conservation_evidence,
    ddg_category,
    gather_evidence,
    make_report,
    stability_evidence,
    structure_evidence,
    tier,
)


def test_tiers_cover_the_unit_interval():
    assert tier(0.95) == "strong model support for pathogenic"
    assert tier(0.5) == "inconclusive"
    assert tier(0.05) == "strong model support for benign"


def test_balanced_probability_removes_training_prior():
    # At the training prevalence the calibrated value is exactly "no evidence".
    assert balanced_probability(0.85, 0.85) == pytest.approx(0.5)
    assert balanced_probability(0.5, 0.5) == pytest.approx(0.5)
    assert balanced_probability(0.95, 0.85) < 0.95


@pytest.mark.parametrize("ddg, expected", [
    (3.0, "strongly destabilising"), (1.2, "destabilising"), (0.3, "neutral"),
    (-1.5, "stabilising"), (float("nan"), "not modelled"), (None, "not modelled"),
])
def test_ddg_categories(ddg, expected):
    assert ddg_category(ddg) == expected


def test_neutral_ddg_is_not_benign_evidence():
    assert stability_evidence({"foldx_ddg": 0.1, "plddt": 90}).direction == "neutral"
    assert stability_evidence({"foldx_ddg": 2.5, "plddt": 90}).direction == "pathogenic"
    missing = stability_evidence({"foldx_ddg": math.nan, "plddt": 40})
    assert missing.direction == "unavailable" and "pLDDT 40" in missing.detail


def test_conservation_and_structure_rules():
    conserved = {"conservation_score": 1.0, "mut_aa_in_orthologs": 0.0, "phylop": 7.0}
    assert conservation_evidence(conserved).direction == "pathogenic"
    assert conservation_evidence({**conserved, "mut_aa_in_orthologs": 0.2}).direction == "benign"
    assert structure_evidence({"plddt": 92, "rsa": 0.05}).direction == "pathogenic"
    assert structure_evidence({"plddt": 40, "rsa": 0.8}).direction == "benign"


def test_conclusion_counts_agreement_and_flags_conflicts():
    row = {"conservation_score": 1.0, "mut_aa_in_orthologs": 0.0, "phylop": 7.0,
           "foldx_ddg": 3.0, "plddt": 92, "rsa": 0.05, "alphamissense": 0.1, "revel": 0.9}
    text, flags = conclusion(0.95, gather_evidence(row), "benign", gene_benign=3)
    assert "3 of 3" in text
    assert any("AlphaMissense disagrees" in f for f in flags)
    assert any("Disagrees with the ClinVar" in f for f in flags)
    assert any("Only 3 benign" in f for f in flags)


def test_legacy_make_report_mentions_ddg():
    report = make_report("MYH7", "p.Arg403Gln", "Pathogenic", 0.9, [{"feature": "phylop", "shap_value": 1.0}],
                         None, 1.5, "n/a")
    assert "+1.50 kcal/mol (destabilising)" in report.plain_language_summary
