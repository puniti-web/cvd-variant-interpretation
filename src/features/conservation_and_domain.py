"""
Step 3 (part B): Conservation and UniProt annotation features.

Conservation (protein level, from build_conservation_scores.py):
  conservation_score      fraction of vertebrate orthologs keeping the
                          human residue
  mut_aa_in_orthologs     fraction of orthologs that already carry the
                          mutant residue (tolerated in another species)

UniProt positional annotations (from fetch_uniprot.py):
  in_annotated_domain     inside a structured domain: Domain, Repeat, Zinc
                          finger, Transmembrane or Intramembrane. Broad
                          "Region" and topology spans (which cover most of a
                          channel) no longer count, so this flag
                          discriminates.
  in_transmembrane        inside a transmembrane / intramembrane segment
  in_functional_site      at an active/binding site, a "Site", a short
                          motif, or a disulfide-bonded cysteine
  dist_to_functional_site residues to the nearest such site (capped)
  in_disordered_region    inside a UniProt "Disordered" region
  in_coiled_coil          inside a coiled-coil segment
"""
import json
import math
from dataclasses import dataclass

from src.config import RAW_DIR

CONSERVATION_FILE = RAW_DIR / "conservation_scores.json"
PROFILE_FILE = RAW_DIR / "conservation_profiles.json"

STRUCTURED_DOMAIN_TYPES = {"Domain", "Repeat", "Zinc finger", "Transmembrane", "Intramembrane"}
MEMBRANE_TYPES = {"Transmembrane", "Intramembrane"}
SITE_TYPES = {"Active site", "Binding site", "Site", "Motif"}
MAX_SITE_DISTANCE = 200


@dataclass
class ConservationFeatures:
    conservation_score: float | None
    in_annotated_domain: bool
    domain_type: str | None


@dataclass
class AnnotationFeatures:
    in_annotated_domain: bool
    domain_type: str | None  # description of the structured domain, e.g. "Myosin motor"
    in_transmembrane: bool
    in_functional_site: bool
    dist_to_functional_site: int
    in_disordered_region: bool
    in_coiled_coil: bool


def _span(feat: dict) -> tuple[int, int] | None:
    start, end = feat.get("start"), feat.get("end")
    if start is None or end is None:
        return None
    return int(start), int(end)


def compute_domain_membership(position: int, uniprot_features: list[dict]) -> tuple[bool, str | None]:
    """Whether a residue lies in a structured domain, and that domain's description."""
    for feat in uniprot_features:
        span = _span(feat)
        if feat.get("type") in STRUCTURED_DOMAIN_TYPES and span and span[0] <= position <= span[1]:
            return True, feat.get("description") or feat.get("type")
    return False, None


def _site_positions(uniprot_features: list[dict]) -> list[int]:
    positions = []
    for feat in uniprot_features:
        span = _span(feat)
        if not span:
            continue
        if feat.get("type") == "Disulfide bond":
            positions.extend(span)  # start/end are the two bonded cysteines
        elif feat.get("type") in SITE_TYPES:
            positions.extend(range(span[0], span[1] + 1))
    return positions


def compute_annotation_features(position: int, uniprot_features: list[dict]) -> AnnotationFeatures:
    in_domain, domain_type = compute_domain_membership(position, uniprot_features)

    def inside(predicate) -> bool:
        return any(
            predicate(feat) and (span := _span(feat)) and span[0] <= position <= span[1]
            for feat in uniprot_features
        )

    sites = _site_positions(uniprot_features)
    distance = min((abs(position - s) for s in sites), default=MAX_SITE_DISTANCE)
    return AnnotationFeatures(
        in_annotated_domain=in_domain,
        domain_type=domain_type,
        in_transmembrane=inside(lambda f: f.get("type") in MEMBRANE_TYPES),
        in_functional_site=distance == 0,
        dist_to_functional_site=min(distance, MAX_SITE_DISTANCE),
        in_disordered_region=inside(
            lambda f: f.get("type") == "Region" and "disordered" in (f.get("description") or "").lower()
        ),
        in_coiled_coil=inside(lambda f: f.get("type") == "Coiled coil"),
    )


def _load_json(path) -> dict:
    if not path.exists():
        return {}
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError):
        return {}


def load_conservation_scores() -> dict:
    """{gene: {position: identity fraction}}"""
    return _load_json(CONSERVATION_FILE)


def load_conservation_profiles() -> dict:
    """{gene: {position: ortholog residues at that position}}"""
    return _load_json(PROFILE_FILE)


def get_conservation_score(gene_symbol: str, position: int, conservation_scores: dict) -> float | None:
    score = conservation_scores.get(gene_symbol, {}).get(str(position))
    if score is None:
        return None
    try:
        score = float(score)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(score) else score


def mutant_frequency_in_orthologs(gene_symbol: str, position: int, mut_aa: str, profiles: dict) -> float | None:
    residues = profiles.get(gene_symbol, {}).get(str(position))
    if not residues:
        return None
    return round(residues.count(mut_aa) / len(residues), 4)


def compute_conservation_and_domain_features(
    gene_symbol: str,
    position: int,
    uniprot_features: list[dict],
    conservation_scores: dict | None = None,
) -> ConservationFeatures:
    """Kept for the API, which calls it per variant."""
    in_domain, domain_type = compute_domain_membership(position, uniprot_features)
    if conservation_scores is None:
        conservation_scores = load_conservation_scores()
    return ConservationFeatures(
        conservation_score=get_conservation_score(gene_symbol, position, conservation_scores),
        in_annotated_domain=in_domain,
        domain_type=domain_type,
    )
