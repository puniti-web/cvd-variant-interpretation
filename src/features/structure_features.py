"""
Step 3 (part D): Residue-level structural context from AlphaFold DB models.

For each gene the predicted full-length structure is downloaded from
AlphaFold DB (via its API, so the current model version is used), checked
against the UniProt canonical sequence, and summarised per residue:

  plddt          AlphaFold confidence (0-100). Low values mark disordered
                 or poorly modelled regions.
  rsa            relative solvent accessibility: Shrake-Rupley SASA divided
                 by the residue's maximum ASA (Tien et al. 2013). Buried
                 residues (low rsa) tolerate substitutions less.
  contact_count  number of C-alpha atoms within 10 Å (|i - j| > 2), a
                 packing-density proxy.
  sec_struct     helix / strand / coil, assigned from backbone phi/psi
                 angles. DSSP would be more accurate but needs a separate
                 binary, so this is an approximation.

AlphaFold DB has no single model for proteins longer than 2,700 residues
(RYR2 in this panel). Those genes get no structural features; the model
imputes them and sees a has_structure flag.

Output: data/raw/alphafold_residue_features.csv

Usage:
    python -m src.features.structure_features
"""
import json
import math
import time

import numpy as np
import pandas as pd
import requests
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.SASA import ShrakeRupley
from tqdm import tqdm

from src.config import RAW_DIR, ensure_dirs, load_gene_panel

AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
# AlphaFold DB file hosting rejects requests without a browser-like agent.
BROWSER_HEADERS = {"User-Agent": "Mozilla/5.0 (cvd-xai research project)"}
ALPHAFOLD_DIR = RAW_DIR / "alphafold"
OUT_PATH = RAW_DIR / "alphafold_residue_features.csv"

# Theoretical maximum ASA in Å², Tien et al. (2013) PLoS ONE 8:e80635
MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
CONTACT_RADIUS = 10.0
MIN_RUN = {"helix": 4, "strand": 3}


def download_alphafold_model(uniprot_id: str) -> tuple[str | None, str | None]:
    """Download the current AlphaFold DB model. Returns (local path, model URL)."""
    ALPHAFOLD_DIR.mkdir(parents=True, exist_ok=True)
    resp = requests.get(AFDB_API.format(uniprot_id=uniprot_id), headers=BROWSER_HEADERS, timeout=30)
    if resp.status_code == 404:
        return None, None
    resp.raise_for_status()
    entries = [e for e in resp.json() if e.get("uniprotStart") == 1]
    if not entries:
        return None, None
    url = entries[0]["pdbUrl"]
    dest = ALPHAFOLD_DIR / url.rsplit("/", 1)[-1]
    if not dest.exists():
        pdb = requests.get(url, headers=BROWSER_HEADERS, timeout=120)
        pdb.raise_for_status()
        dest.write_bytes(pdb.content)
    return str(dest), url


def classify_phi_psi(phi: float | None, psi: float | None) -> str:
    """Coarse Ramachandran-region secondary structure assignment."""
    if phi is None or psi is None:
        return "coil"
    phi, psi = math.degrees(phi), math.degrees(psi)
    if -160 <= phi <= -20 and -120 <= psi <= 50:
        return "helix"
    if -180 <= phi <= -40 and (psi >= 90 or psi <= -150):
        return "strand"
    return "coil"


def smooth_sec_struct(labels: list[str], plddt: list[float]) -> list[str]:
    """
    Keep helix only in runs of >= MIN_RUN["helix"] residues and strand in
    runs of >= MIN_RUN["strand"]; isolated assignments become coil. Very
    low-confidence residues (pLDDT < 50) are treated as coil because their
    backbone geometry is not meaningful.
    """
    out = ["coil" if p < 50 else s for s, p in zip(labels, plddt)]
    i = 0
    while i < len(out):
        j = i
        while j < len(out) and out[j] == out[i]:
            j += 1
        if out[i] in MIN_RUN and j - i < MIN_RUN[out[i]]:
            out[i:j] = ["coil"] * (j - i)
        i = j
    return out


def residue_features(pdb_path: str, sequence: str) -> pd.DataFrame | None:
    structure = PDBParser(QUIET=True).get_structure("af", pdb_path)
    chain = next(structure[0].get_chains())
    residues = [r for r in chain if r.id[0] == " " and r.get_resname() in THREE_TO_ONE]

    model_seq = "".join(THREE_TO_ONE[r.get_resname()] for r in residues)
    if model_seq != sequence:
        return None  # model built on a different sequence version

    ShrakeRupley().compute(structure[0], level="R")

    phi_psi = {}
    for pp in PPBuilder().build_peptides(chain):
        for res, (phi, psi) in zip(pp, pp.get_phi_psi_list()):
            phi_psi[res.id[1]] = (phi, psi)

    ca = np.array([r["CA"].coord for r in residues])
    dist = np.sqrt(((ca[:, None, :] - ca[None, :, :]) ** 2).sum(-1))
    idx = np.arange(len(residues))
    near = (dist < CONTACT_RADIUS) & (np.abs(idx[:, None] - idx[None, :]) > 2)
    contacts = near.sum(1)

    plddt = [float(r["CA"].get_bfactor()) for r in residues]
    raw_ss = [classify_phi_psi(*phi_psi.get(r.id[1], (None, None))) for r in residues]
    sec_struct = smooth_sec_struct(raw_ss, plddt)

    rows = []
    for i, r in enumerate(residues):
        rows.append({
            "position": r.id[1],
            "plddt": plddt[i],
            "rsa": min(float(r.sasa) / MAX_ASA[r.get_resname()], 1.0),
            "contact_count": int(contacts[i]),
            "sec_struct": sec_struct[i],
        })
    return pd.DataFrame(rows)


def main():
    ensure_dirs()
    with open(RAW_DIR / "uniprot_entries.json", encoding="utf-8") as f:
        sequences = {e["gene_symbol"]: e["sequence"] for e in json.load(f)}

    frames, manifest = [], []
    for gene in tqdm(load_gene_panel("both"), desc="AlphaFold features"):
        symbol, uid = gene["symbol"], gene["uniprot_id"]
        status, url = "ok", None
        try:
            path, url = download_alphafold_model(uid)
        except requests.RequestException as e:
            path, status = None, f"download failed: {e}"
        if path is None and status == "ok":
            status = "no single full-length AlphaFold DB model"
        if path:
            feats = residue_features(path, sequences[symbol])
            if feats is None:
                status = "model sequence differs from UniProt canonical"
            else:
                feats.insert(0, "gene_symbol", symbol)
                frames.append(feats)
        manifest.append({"gene_symbol": symbol, "uniprot_id": uid, "model_url": url, "status": status})
        time.sleep(0.2)

    pd.concat(frames, ignore_index=True).to_csv(OUT_PATH, index=False)
    pd.DataFrame(manifest).to_csv(RAW_DIR / "alphafold_manifest.csv", index=False)
    print(pd.DataFrame(manifest)[["gene_symbol", "status"]].to_string(index=False))
    print(f"Saved residue features to {OUT_PATH}")


if __name__ == "__main__":
    main()
