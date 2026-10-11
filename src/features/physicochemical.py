"""
Step 3 (part A): Physicochemical property features.

Standard amino acid property tables (Kyte-Doolittle hydrophobicity,
molecular weight, side-chain charge, polarity, side-chain volume). For each
variant we compute the *change* (mutant - wild-type) in each property,
which is generally more informative for pathogenicity prediction than the
absolute values.

This module has no external dependencies and needs no network access, so
it can be developed and unit-tested entirely offline.
"""
import math
from dataclasses import dataclass

from Bio.Align import substitution_matrices

BLOSUM62 = substitution_matrices.load("BLOSUM62")

# Kyte & Doolittle (1982) hydropathy index
HYDROPHOBICITY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Average molecular weight (g/mol) of the amino acid residue
MOLECULAR_WEIGHT = {
    "A": 89.1, "R": 174.2, "N": 132.1, "D": 133.1, "C": 121.2,
    "Q": 146.2, "E": 147.1, "G": 75.1, "H": 155.2, "I": 131.2,
    "L": 131.2, "K": 146.2, "M": 149.2, "F": 165.2, "P": 115.1,
    "S": 105.1, "T": 119.1, "W": 204.2, "Y": 181.2, "V": 117.1,
}

# Side-chain volume (A^3), Zamyatnin (1972)
VOLUME = {
    "A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5,
    "Q": 143.8, "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7,
    "L": 166.7, "K": 168.6, "M": 162.9, "F": 189.9, "P": 112.7,
    "S": 89.0, "T": 116.1, "W": 227.8, "Y": 193.6, "V": 140.0,
}

# Net charge at physiological pH
CHARGE = {
    "D": -1, "E": -1, "K": 1, "R": 1, "H": 0.1,
    # all others default to 0 via .get(aa, 0)
}

# Polar vs nonpolar side chains
POLAR = set("RNDQEHKSTY")  # polar / charged
NONPOLAR = set("AVLIMFWPGC")  # nonpolar / hydrophobic

# Grantham (1974) side-chain composition, polarity and volume, used to
# compute the Grantham physicochemical distance between two residues.
GRANTHAM_PROPS = {
    #     composition, polarity, volume
    "S": (1.42, 9.2, 32.0), "R": (0.65, 10.5, 124.0), "L": (0.00, 4.9, 111.0),
    "P": (0.39, 8.0, 32.5), "T": (0.71, 8.6, 61.0), "A": (0.00, 8.1, 31.0),
    "V": (0.00, 5.9, 84.0), "G": (0.74, 9.0, 3.0), "I": (0.00, 5.2, 111.0),
    "F": (0.00, 5.2, 132.0), "Y": (0.20, 6.2, 136.0), "C": (2.75, 5.5, 55.0),
    "H": (0.58, 10.4, 96.0), "Q": (0.89, 10.5, 85.0), "N": (1.33, 11.6, 56.0),
    "K": (0.33, 11.3, 119.0), "D": (1.38, 13.0, 54.0), "E": (0.92, 12.3, 83.0),
    "M": (0.00, 5.7, 105.0), "W": (0.13, 5.4, 170.0),
}


def grantham_distance(wt_aa: str, mut_aa: str) -> float:
    """Grantham distance: 0 (identical) to ~215 (Cys <-> Trp)."""
    (c1, p1, v1), (c2, p2, v2) = GRANTHAM_PROPS[wt_aa], GRANTHAM_PROPS[mut_aa]
    d = math.sqrt(1.833 * (c1 - c2) ** 2 + 0.1018 * (p1 - p2) ** 2 + 0.000399 * (v1 - v2) ** 2)
    return round(50.723 * d, 1)


@dataclass
class VariantPhysChemFeatures:
    delta_hydrophobicity: float
    delta_molecular_weight: float
    delta_volume: float
    delta_charge: float
    polarity_change: bool  # True if polar<->nonpolar class changes
    wt_is_proline: bool
    mut_is_proline: bool  # proline introduction/removal often disrupts secondary structure
    wt_is_glycine: bool
    mut_is_glycine: bool  # glycine loss/gain affects backbone flexibility
    blosum62: float  # substitution log-odds; negative = rarely seen in evolution
    grantham: float  # physicochemical distance between the two residues


def compute_physchem_features(wt_aa: str, mut_aa: str) -> VariantPhysChemFeatures:
    wt_aa, mut_aa = wt_aa.upper(), mut_aa.upper()

    def polarity_class(aa):
        if aa in POLAR:
            return "polar"
        if aa in NONPOLAR:
            return "nonpolar"
        return "other"

    return VariantPhysChemFeatures(
        delta_hydrophobicity=HYDROPHOBICITY.get(mut_aa, 0) - HYDROPHOBICITY.get(wt_aa, 0),
        delta_molecular_weight=MOLECULAR_WEIGHT.get(mut_aa, 0) - MOLECULAR_WEIGHT.get(wt_aa, 0),
        delta_volume=VOLUME.get(mut_aa, 0) - VOLUME.get(wt_aa, 0),
        delta_charge=CHARGE.get(mut_aa, 0) - CHARGE.get(wt_aa, 0),
        polarity_change=polarity_class(wt_aa) != polarity_class(mut_aa),
        wt_is_proline=(wt_aa == "P"),
        mut_is_proline=(mut_aa == "P"),
        wt_is_glycine=(wt_aa == "G"),
        mut_is_glycine=(mut_aa == "G"),
        blosum62=float(BLOSUM62[wt_aa][mut_aa]),
        grantham=grantham_distance(wt_aa, mut_aa),
    )
