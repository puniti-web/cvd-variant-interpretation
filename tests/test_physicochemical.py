"""
Unit tests for src/features/physicochemical.py — this module needs no
network access, so it's the easiest place to start verifying the pipeline
actually works before wiring up live API calls.

Run with:
    pytest tests/test_physicochemical.py -v
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.features.physicochemical import compute_physchem_features, grantham_distance


def test_charge_reversal_arg_to_glu():
    # R403Q is a well-known MYH7 hypertrophic cardiomyopathy mutation.
    # Arg (R) is positively charged; Gln (Q) is neutral -> delta_charge < 0.
    f = compute_physchem_features("R", "Q")
    assert f.delta_charge < 0
    assert f.polarity_change is False  # both R and Q are polar in our classification


def test_polarity_change_detected():
    # Leucine (nonpolar) -> Arginine (polar) should flip polarity class.
    f = compute_physchem_features("L", "R")
    assert f.polarity_change is True


def test_proline_introduction_flagged():
    f = compute_physchem_features("A", "P")
    assert f.mut_is_proline is True
    assert f.wt_is_proline is False


def test_identity_substitution_has_zero_deltas():
    f = compute_physchem_features("A", "A")
    assert f.delta_hydrophobicity == 0
    assert f.delta_molecular_weight == 0
    assert f.delta_volume == 0
    assert f.delta_charge == 0
    assert f.polarity_change is False


def test_grantham_matches_published_values():
    # Grantham (1974) Table 2
    assert abs(grantham_distance("R", "Q") - 43) <= 2
    assert abs(grantham_distance("L", "I") - 5) <= 2
    assert abs(grantham_distance("C", "W") - 215) <= 2
    assert grantham_distance("A", "A") == 0


def test_blosum62():
    assert compute_physchem_features("R", "K").blosum62 == 2
    assert compute_physchem_features("W", "G").blosum62 == -2


if __name__ == "__main__":
    test_charge_reversal_arg_to_glu()
    test_polarity_change_detected()
    test_proline_introduction_flagged()
    test_identity_substitution_has_zero_deltas()
    print("All tests passed.")
