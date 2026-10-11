"""
Step 1 (part D): Annotate every variant with Ensembl VEP.

Variants are sent to the Ensembl VEP REST API (https://rest.ensembl.org,
free, no auth) in batches of 200, using the genomic coordinates from
ClinVar's canonical SPDI. Two kinds of annotation come back:

  Conservation (used as model features in Step 3)
    - phyloP100way_vertebrate, phastCons100way_vertebrate, GERP++ RS
      (via the dbNSFP plugin)

  Published predictors (used only as BASELINES in Step 4, never as
  features: REVEL and AlphaMissense were trained or calibrated on ClinVar
  labels, so feeding them to the model would leak the label)
    - SIFT, PolyPhen-2, AlphaMissense, REVEL

Output: data/raw/vep_annotations.csv, one row per SPDI.

Usage:
    python -m src.data_collection.fetch_ensembl_vep [--labelled-only]
"""
import argparse
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import requests
from tqdm import tqdm

from src.config import PROCESSED_DIR, RAW_DIR, USER_AGENT, ensure_dirs

VEP_URL = "https://rest.ensembl.org/vep/homo_sapiens/region"
VEP_PARAMS = {
    "dbNSFP": "phyloP100way_vertebrate,phastCons100way_vertebrate,GERP++_RS,REVEL_score",
    "AlphaMissense": 1,
    "canonical": 1,
    "mane": 1,
}
BATCH_SIZE = 200  # VEP REST POST limit
REQUEST_DELAY = 0.2
# Batches with the dbNSFP plugin take 1-2 minutes server-side, so a few run
# in parallel. Ensembl allows 15 requests/sec, far above this.
WORKERS = 4
CHECKPOINT_EVERY = 5  # batches
MAX_RETRIES = 5
OUT_PATH = RAW_DIR / "vep_annotations.csv"


def spdi_to_vep_region(spdi: str) -> str | None:
    """
    Convert a ClinVar SPDI for an SNV (e.g. NC_000014.9:23429278:G:A, where
    the position is 0-based) into VEP's VCF-like region string
    ("14 23429279 . G A . . .").
    """
    if not isinstance(spdi, str):
        return None
    parts = spdi.split(":")
    if len(parts) != 4 or not parts[0].startswith("NC_"):
        return None
    accession, pos, ref, alt = parts
    if len(ref) != 1 or len(alt) != 1:
        return None
    chrom_num = int(accession[3:].split(".")[0])
    chrom = {23: "X", 24: "Y"}.get(chrom_num, str(chrom_num))
    return f"{chrom} {int(pos) + 1} . {ref} {alt} . . ."


def _post(regions: list[str]) -> list[dict]:
    for attempt in range(MAX_RETRIES):
        try:
            resp = requests.post(
                VEP_URL,
                params=VEP_PARAMS,
                json={"variants": regions},
                headers={"Content-Type": "application/json", "Accept": "application/json", "User-Agent": USER_AGENT},
                timeout=300,
            )
            if resp.status_code == 429:
                time.sleep(float(resp.headers.get("Retry-After", 2 ** attempt)))
                continue
            if resp.status_code >= 500:
                raise requests.HTTPError(f"HTTP {resp.status_code}", response=resp)
            resp.raise_for_status()
            return resp.json()
        except requests.RequestException:
            if attempt == MAX_RETRIES - 1:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError("VEP request kept hitting the rate limit")


def _first_number(value) -> float | None:
    """dbNSFP fields can be per-transcript lists like '.,0.82,.'; take the max numeric value."""
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    numbers = [float(v) for v in str(value).split(",") if v not in (".", "")]
    return max(numbers) if numbers else None


def _pick_transcript(consequences: list[dict], gene_symbol: str) -> dict | None:
    """Prefer the gene's MANE Select transcript, then its Ensembl canonical one."""
    own = [t for t in consequences if t.get("gene_symbol") == gene_symbol]
    for predicate in (lambda t: t.get("mane_select"), lambda t: t.get("canonical"), lambda t: True):
        for t in own:
            if predicate(t) and "missense_variant" in t.get("consequence_terms", []):
                return t
    return own[0] if own else None


def parse_vep_result(result: dict, gene_symbol: str) -> dict:
    t = _pick_transcript(result.get("transcript_consequences", []), gene_symbol) or {}
    am = t.get("alphamissense") or {}
    return {
        "vep_transcript": t.get("transcript_id"),
        "vep_protein_position": t.get("protein_start"),
        "vep_amino_acids": t.get("amino_acids"),
        "phylop": _first_number(t.get("phylop100way_vertebrate")),
        "phastcons": _first_number(t.get("phastcons100way_vertebrate")),
        "gerp": _first_number(t.get("gerp++_rs")),
        "sift_score": t.get("sift_score"),
        "polyphen_score": t.get("polyphen_score"),
        "alphamissense": am.get("am_pathogenicity"),
        "revel": _first_number(t.get("revel_score")),
    }


def _annotate_batch(batch: list[str], gene_by_region: dict, spdi_by_region: dict) -> list[dict]:
    rows = []
    for result in _post(batch):
        region = result.get("input")
        if region in gene_by_region:
            rows.append({"spdi": spdi_by_region[region], **parse_vep_result(result, gene_by_region[region])})
    time.sleep(REQUEST_DELAY)
    return rows


def annotate(variants: pd.DataFrame, checkpoint=None) -> pd.DataFrame:
    """
    variants needs columns spdi and gene_symbol; returns one annotation row
    per SPDI. checkpoint(rows) is called every few batches so long runs can
    resume after an interruption.
    """
    variants = variants.dropna(subset=["spdi"]).drop_duplicates("spdi").copy()
    variants["region"] = variants["spdi"].map(spdi_to_vep_region)
    variants = variants.dropna(subset=["region"])
    gene_by_region = dict(zip(variants["region"], variants["gene_symbol"]))
    spdi_by_region = dict(zip(variants["region"], variants["spdi"]))

    regions = list(gene_by_region)
    batches = [regions[i : i + BATCH_SIZE] for i in range(0, len(regions), BATCH_SIZE)]
    rows: list[dict] = []
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futures = [pool.submit(_annotate_batch, b, gene_by_region, spdi_by_region) for b in batches]
        for n, future in enumerate(tqdm(as_completed(futures), total=len(futures), desc="VEP batches"), 1):
            try:
                rows.extend(future.result())
            except requests.RequestException as e:
                print(f"[warn] VEP batch failed: {e}")
            if checkpoint and n % CHECKPOINT_EVERY == 0:
                checkpoint(rows)
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--labelled-only", action="store_true", help="skip the VUS set")
    args = parser.parse_args()

    ensure_dirs()
    names = ["variants_labelled.csv"] + ([] if args.labelled_only else ["variants_vus.csv"])
    frames = [
        pd.read_csv(PROCESSED_DIR / name, usecols=["spdi", "gene_symbol"])
        for name in names
        if (PROCESSED_DIR / name).exists()
    ]
    if not frames:
        raise FileNotFoundError("Run `python -m src.preprocessing.clean_variants` first.")
    variants = pd.concat(frames, ignore_index=True)

    # Re-use earlier annotations so an interrupted run can resume.
    done = pd.read_csv(OUT_PATH) if OUT_PATH.exists() else pd.DataFrame(columns=["spdi"])
    todo = variants[~variants["spdi"].isin(done["spdi"])]
    print(f"{len(done)} SPDIs already annotated, {todo['spdi'].nunique()} to go")

    def save(rows):
        combined = pd.concat([done, pd.DataFrame(rows)], ignore_index=True).drop_duplicates("spdi", keep="last")
        combined.to_csv(OUT_PATH, index=False)
        return combined

    combined = save(annotate(todo, checkpoint=save).to_dict("records"))
    print(f"Saved {len(combined)} VEP annotations to {OUT_PATH}")


if __name__ == "__main__":
    main()
