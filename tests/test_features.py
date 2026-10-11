"""Offline tests for annotation, conservation-profile and structure helpers."""
from src.data_collection.fetch_ensembl_vep import spdi_to_vep_region, _first_number
from src.features.conservation_and_domain import (
    MAX_SITE_DISTANCE,
    compute_annotation_features,
    compute_domain_membership,
    mutant_frequency_in_orthologs,
)
from src.features.structure_features import smooth_sec_struct

FEATURES = [
    {"type": "Domain", "description": "Myosin motor", "start": 10, "end": 50},
    {"type": "Region", "description": "Disordered", "start": 60, "end": 80},
    {"type": "Region", "description": "Actin-binding", "start": 1, "end": 200},
    {"type": "Binding site", "description": "ATP", "start": 30, "end": 30},
    {"type": "Disulfide bond", "description": None, "start": 90, "end": 120},
    {"type": "Transmembrane", "description": "Helical", "start": 130, "end": 150},
    {"type": "Coiled coil", "description": None, "start": 160, "end": 190},
]


def test_broad_region_is_not_a_domain():
    assert compute_domain_membership(55, FEATURES) == (False, None)
    assert compute_domain_membership(20, FEATURES) == (True, "Myosin motor")


def test_annotation_features():
    f = compute_annotation_features(30, FEATURES)
    assert f.in_functional_site and f.dist_to_functional_site == 0

    f = compute_annotation_features(70, FEATURES)
    assert f.in_disordered_region and not f.in_annotated_domain
    assert f.dist_to_functional_site == 20  # nearest is the Cys at 90

    f = compute_annotation_features(140, FEATURES)
    assert f.in_transmembrane and f.in_annotated_domain

    assert compute_annotation_features(170, FEATURES).in_coiled_coil
    assert compute_annotation_features(5, []).dist_to_functional_site == MAX_SITE_DISTANCE


def test_disulfide_counts_only_the_bonded_cysteines():
    assert compute_annotation_features(120, FEATURES).in_functional_site
    assert not compute_annotation_features(105, FEATURES).in_functional_site


def test_mutant_frequency_in_orthologs():
    profiles = {"G": {"5": "RRRRQ"}}
    assert mutant_frequency_in_orthologs("G", 5, "Q", profiles) == 0.2
    assert mutant_frequency_in_orthologs("G", 5, "W", profiles) == 0.0
    assert mutant_frequency_in_orthologs("G", 6, "Q", profiles) is None


def test_secondary_structure_smoothing():
    labels = ["helix"] * 5 + ["strand"] * 2 + ["helix"] * 3 + ["strand"] * 3
    plddt = [90] * 12 + [30]
    out = smooth_sec_struct(labels, plddt)
    assert out[:5] == ["helix"] * 5
    assert out[5:10] == ["coil"] * 5          # short runs dropped
    assert out[10:12] == ["coil"] * 2         # strand run broken by low pLDDT
    assert out[12] == "coil"


def test_spdi_conversion():
    assert spdi_to_vep_region("NC_000014.9:23429278:G:A") == "14 23429279 . G A . . ."
    assert spdi_to_vep_region("NC_000023.11:100:G:A") == "X 101 . G A . . ."
    assert spdi_to_vep_region("NC_000014.9:23429278:GG:A") is None
    assert spdi_to_vep_region(None) is None


def test_dbnsfp_multi_transcript_values():
    assert _first_number(".,0.822,.") == 0.822
    assert _first_number(".,.") is None
    assert _first_number(1.5) == 1.5
