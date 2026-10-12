"""
Step 5: FoldX stability change (ddG) for every missense variant with a
confident AlphaFold structure.

Three subcommands, so the heavy part can run on a remote machine:

  prepare   (local)  Cut each AlphaFold model into chunks, assign variants
                     to chunks, write FoldX mutant batches. Output is one
                     self-contained work directory.
  run       (remote) RepairPDB each chunk, then BuildModel each batch, in
                     parallel. Standard library only, so the remote machine
                     needs nothing but python3 and the FoldX binary. Every
                     batch writes its own result file, so a restarted run
                     resumes where it stopped.
  package   (local)  Bundle the work directory and this script for upload.
  collect   (local)  Merge batch results into data/raw/foldx_ddg.csv.

Approximations, chosen to keep the full run to a few hours:
  - Only variants at residues with pLDDT >= 70 are modelled. FoldX energies
    on low-confidence backbone are not meaningful; other variants get NaN
    and the classifier imputes them.
  - Residues with pLDDT < 50 are removed before modelling.
  - Proteins longer than WHOLE_MAX residues are cut into chunks: a core of
    CORE_SIZE consecutive residues plus every residue whose C-alpha lies
    within SHELL_RADIUS of a core C-alpha. A variant is modelled in the
    chunk whose core contains it, so its whole local environment is
    present.
  - One BuildModel run per variant. On a 10-variant test, 3 runs changed
    ddG by at most 0.07 kcal/mol and took 3x longer.

ddG = mutant - wild type total energy in kcal/mol; positive is destabilising.

Usage:
    python -m src.structure.foldx_ddg prepare
    python3 src/structure/foldx_ddg.py run --foldx /opt/foldx/foldx_20261231 --workers 32
    python -m src.structure.foldx_ddg collect
"""
import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
WORK_DIR = PROJECT_ROOT / "data" / "raw" / "foldx"
OUT_PATH = PROJECT_ROOT / "data" / "raw" / "foldx_ddg.csv"

SITE_PLDDT_MIN = 70.0
KEEP_PLDDT_MIN = 50.0
WHOLE_MAX = 700
CORE_SIZE = 300
SHELL_RADIUS = 12.0
BATCH_SIZE = 10
NUMBER_OF_RUNS = 1
REPAIR_TIMEOUT = 4 * 3600
BUILD_TIMEOUT = 2 * 3600


# ---------------------------------------------------------------- prepare

def chunk_model(pdb_path: str, plddt_min: float = KEEP_PLDDT_MIN):
    """
    Split one AlphaFold model into chunks.
    Returns [(core_start, core_end, residue_numbers)], numbering as in the model.
    """
    import numpy as np
    from Bio.PDB import PDBParser

    chain = next(PDBParser(QUIET=True).get_structure("af", pdb_path)[0].get_chains())
    residues = [r for r in chain if r.id[0] == " " and "CA" in r and r["CA"].get_bfactor() >= plddt_min]
    numbers = np.array([r.id[1] for r in residues])
    if len(residues) <= WHOLE_MAX:
        return [(1, int(numbers.max()), numbers.tolist())]

    ca = np.array([r["CA"].coord for r in residues])
    chunks = []
    for start in range(1, int(numbers.max()) + 1, CORE_SIZE):
        end = start + CORE_SIZE - 1
        core = (numbers >= start) & (numbers <= end)
        if not core.any():
            continue
        dist = np.sqrt(((ca[:, None, :] - ca[None, core, :]) ** 2).sum(-1)).min(1)
        chunks.append((start, end, numbers[dist <= SHELL_RADIUS].tolist()))
    return chunks


def write_chunk(pdb_path: str, keep: list[int], dest: Path):
    from Bio.PDB import PDBIO, PDBParser, Select

    keep = set(keep)

    class Keep(Select):
        def accept_residue(self, residue):
            return residue.id[0] == " " and residue.id[1] in keep

    io = PDBIO()
    io.set_structure(PDBParser(QUIET=True).get_structure("af", pdb_path))
    io.save(str(dest), Keep())


def prepare(work_dir: Path = WORK_DIR, batch_size: int = BATCH_SIZE):
    import pandas as pd

    from src.config import PROCESSED_DIR, RAW_DIR, load_gene_panel

    chunk_dir = work_dir / "chunks"
    chunk_dir.mkdir(parents=True, exist_ok=True)

    variants = pd.concat([
        pd.read_csv(PROCESSED_DIR / "variants_labelled.csv"),
        pd.read_csv(PROCESSED_DIR / "variants_vus.csv"),
    ])[["gene_symbol", "position", "wt_aa", "mut_aa", "label"]]
    variants = variants.drop_duplicates(["gene_symbol", "position", "wt_aa", "mut_aa"])
    residues = pd.read_csv(RAW_DIR / "alphafold_residue_features.csv")[["gene_symbol", "position", "plddt"]]
    variants = variants.merge(residues, on=["gene_symbol", "position"], how="left")
    manifest = pd.read_csv(RAW_DIR / "alphafold_manifest.csv").set_index("gene_symbol")

    chunks, jobs = [], []
    for gene in load_gene_panel("both"):
        symbol = gene["symbol"]
        url = manifest.loc[symbol, "model_url"]
        sites = variants[(variants["gene_symbol"] == symbol) & (variants["plddt"] >= SITE_PLDDT_MIN)]
        if not isinstance(url, str) or sites.empty:
            continue
        model = RAW_DIR / "alphafold" / url.rsplit("/", 1)[-1]
        for start, end, keep in chunk_model(str(model)):
            in_core = sites[(sites["position"] >= start) & (sites["position"] <= end)]
            if in_core.empty:
                continue
            name = f"{symbol}_{start}-{end}"
            write_chunk(str(model), keep, chunk_dir / f"{name}.pdb")
            chunks.append({"chunk": name, "gene_symbol": symbol, "n_residues": len(keep),
                           "n_variants": len(in_core)})
            # Labelled variants first, so the training set finishes early.
            ordered = in_core.assign(vus=in_core["label"] == "vus").sort_values(["vus", "position"])
            rows = ordered[["position", "wt_aa", "mut_aa", "vus"]].to_dict("records")
            for i in range(0, len(rows), batch_size):
                batch = rows[i:i + batch_size]
                jobs.append({
                    "job_id": f"{name}_{i // batch_size:04d}",
                    "chunk": name,
                    "gene_symbol": symbol,
                    "priority": int(all(r["vus"] for r in batch)),
                    "mutations": [f"{r['wt_aa']}A{r['position']}{r['mut_aa']}" for r in batch],
                })

    with open(work_dir / "chunks.json", "w") as f:
        json.dump(chunks, f, indent=1)
    with open(work_dir / "jobs.json", "w") as f:
        json.dump(jobs, f, indent=1)

    n_var = sum(len(j["mutations"]) for j in jobs)
    print(f"{len(chunks)} chunks, {len(jobs)} batches, {n_var} variants "
          f"({(variants['plddt'] < SITE_PLDDT_MIN).sum() + variants['plddt'].isna().sum()} skipped: "
          f"pLDDT < {SITE_PLDDT_MIN:g} or no model)")
    print(pd.DataFrame(chunks).to_string(index=False))
    print(f"Work directory: {work_dir}")


# -------------------------------------------------------------------- run

def shard_chunks(chunks: list[dict], weights: list[float], shard: int) -> set[str]:
    """Greedy, deterministic split of chunks across machines of unequal size."""
    load = [0.0] * len(weights)
    owner = {}
    for c in sorted(chunks, key=lambda c: (-c["n_variants"] * c["n_residues"], c["chunk"])):
        i = min(range(len(weights)), key=lambda k: (load[k] + c["n_variants"] * c["n_residues"]) / weights[k])
        load[i] += c["n_variants"] * c["n_residues"]
        owner[c["chunk"]] = i
    return {name for name, i in owner.items() if i == shard}


def foldx_env(foldx: Path, cwd: Path):
    """FoldX 5.1 looks for its molecules/ folder in the working directory."""
    molecules = foldx.parent / "molecules"
    link = cwd / "molecules"
    if molecules.is_dir() and not link.exists():
        try:
            link.symlink_to(molecules, target_is_directory=True)
        except OSError:
            shutil.copytree(molecules, link)


def run_foldx(foldx: Path, cwd: Path, args: list[str], timeout: int, log: Path) -> bool:
    foldx_env(foldx, cwd)
    with open(log, "w") as out:
        try:
            proc = subprocess.run([str(foldx), *args], cwd=cwd, stdout=out, stderr=subprocess.STDOUT,
                                  timeout=timeout)
        except subprocess.TimeoutExpired:
            out.write("\nTIMEOUT\n")
            return False
    return proc.returncode == 0


def repair(foldx: Path, work_dir: Path, chunk: str) -> Path | None:
    repaired = work_dir / "repaired" / f"{chunk}_Repair.pdb"
    if repaired.exists():
        return repaired
    tmp = work_dir / "tmp" / f"repair_{chunk}"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    shutil.copy(work_dir / "chunks" / f"{chunk}.pdb", tmp / f"{chunk}.pdb")
    t0 = time.time()
    ok = run_foldx(foldx, tmp, ["--command=RepairPDB", f"--pdb={chunk}.pdb"], REPAIR_TIMEOUT,
                   tmp / "foldx.log")
    out = tmp / f"{chunk}_Repair.pdb"
    if not ok or not out.exists():
        print(f"[repair failed] {chunk}, see {tmp / 'foldx.log'}", flush=True)
        return None
    repaired.parent.mkdir(exist_ok=True)
    shutil.move(str(out), repaired)
    shutil.rmtree(tmp, ignore_errors=True)
    print(f"[repaired] {chunk} in {time.time() - t0:.0f} s", flush=True)
    return repaired


def parse_dif(path: Path, n: int, runs: int = NUMBER_OF_RUNS) -> list[float]:
    """Mean total-energy ddG per mutant, in mutant-file order (runs are consecutive rows)."""
    lines = path.read_text().splitlines()
    header = next(i for i, line in enumerate(lines) if line.startswith("Pdb\t"))
    cols = lines[header].split("\t")
    k = cols.index("total energy")
    values = [float(line.split("\t")[k]) for line in lines[header + 1:] if line.strip()]
    if len(values) != n * runs:
        raise ValueError(f"{path}: expected {n * runs} rows, found {len(values)}")
    return [sum(values[i * runs:(i + 1) * runs]) / runs for i in range(n)]


def build_batch(foldx: Path, repaired: Path, tmp: Path, mutations: list[str]) -> list[float] | None:
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    shutil.copy(repaired, tmp / repaired.name)
    # FoldX recognises the individual-list format by the file name.
    (tmp / "individual_list.txt").write_text("".join(f"{m};\n" for m in mutations))
    ok = run_foldx(foldx, tmp, [
        "--command=BuildModel", f"--pdb={repaired.name}", "--mutant-file=individual_list.txt",
        f"--numberOfRuns={NUMBER_OF_RUNS}", "--out-pdb=false",
    ], BUILD_TIMEOUT, tmp / "foldx.log")
    dif = tmp / f"Dif_{repaired.stem}.fxout"
    if not ok or not dif.exists():
        return None
    try:
        return parse_dif(dif, len(mutations))
    except (ValueError, StopIteration):
        return None


def build(foldx: Path, work_dir: Path, job: dict, repaired: Path) -> int:
    out = work_dir / "results" / f"{job['job_id']}.csv"
    if out.exists():
        return 0
    tmp = work_dir / "tmp" / job["job_id"]
    ddg = build_batch(foldx, repaired, tmp, job["mutations"])
    if ddg is None:
        # One bad mutation fails the whole batch; retry them one at a time.
        ddg = [(r or [None])[0] for r in
               (build_batch(foldx, repaired, tmp, [m]) for m in job["mutations"])]
    rows = [(job["gene_symbol"], m[0], int(m[2:-1]), m[-1], job["chunk"], "" if d is None else f"{d:.4f}")
            for m, d in zip(job["mutations"], ddg)]
    out.parent.mkdir(exist_ok=True)
    with open(out.with_suffix(".part"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene_symbol", "wt_aa", "position", "mut_aa", "chunk", "foldx_ddg"])
        w.writerows(rows)
    os.replace(out.with_suffix(".part"), out)
    shutil.rmtree(tmp, ignore_errors=True)
    return len(rows)


def run(foldx: Path, work_dir: Path, workers: int, weights: list[float], shard: int,
        only_chunks: list[str] | None, max_jobs: int | None):
    foldx = foldx.resolve()
    chunks = json.loads((work_dir / "chunks.json").read_text())
    jobs = json.loads((work_dir / "jobs.json").read_text())
    mine = shard_chunks(chunks, weights, shard)
    if only_chunks:
        mine &= set(only_chunks)
    jobs = sorted((j for j in jobs if j["chunk"] in mine), key=lambda j: (j["priority"], j["job_id"]))
    jobs = [j for j in jobs if not (work_dir / "results" / f"{j['job_id']}.csv").exists()]
    if max_jobs:
        jobs = jobs[:max_jobs]
    needed = sorted({j["chunk"] for j in jobs},
                    key=lambda c: -next(x["n_residues"] for x in chunks if x["chunk"] == c))
    n_total = sum(len(j["mutations"]) for j in jobs)
    print(f"shard {shard}: {len(needed)} chunks, {len(jobs)} batches, {n_total} variants, "
          f"{workers} workers", flush=True)

    t0, done = time.time(), 0
    with ThreadPoolExecutor(workers) as pool:
        # Largest chunks are repaired first; each chunk's batches are queued
        # as soon as its repair finishes.
        repairs = {pool.submit(repair, foldx, work_dir, c): c for c in needed}
        builds = []
        for fut in as_completed(repairs):
            chunk, repaired = repairs[fut], fut.result()
            if repaired is None:
                continue
            builds += [pool.submit(build, foldx, work_dir, j, repaired) for j in jobs if j["chunk"] == chunk]
        for fut in as_completed(builds):
            done += fut.result()
            rate = done / max(time.time() - t0, 1)
            print(f"[{time.strftime('%H:%M:%S')}] {done}/{n_total} variants, "
                  f"ETA {(n_total - done) / max(rate, 1e-9) / 60:.0f} min", flush=True)
    print(f"done in {(time.time() - t0) / 60:.1f} min", flush=True)


def package(work_dir: Path = WORK_DIR) -> Path:
    """Bundle the inputs and this script for upload: work.tar.gz next to work_dir."""
    import tarfile

    dest = work_dir.parent / "foldx_work.tar.gz"
    with tarfile.open(dest, "w:gz") as tar:
        for name in ("chunks", "chunks.json", "jobs.json"):
            tar.add(work_dir / name, arcname=f"work/{name}")
        tar.add(Path(__file__), arcname="foldx_ddg.py")
    print(f"Packaged {dest} ({dest.stat().st_size / 1e6:.1f} MB)")
    return dest


# ---------------------------------------------------------------- collect

def collect(work_dir: Path = WORK_DIR):
    import pandas as pd

    files = sorted((work_dir / "results").glob("*.csv"))
    if not files:
        raise FileNotFoundError(f"no results in {work_dir / 'results'}")
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    df = df.drop_duplicates(["gene_symbol", "position", "wt_aa", "mut_aa"])
    df.to_csv(OUT_PATH, index=False)
    jobs = json.loads((work_dir / "jobs.json").read_text())
    expected = sum(len(j["mutations"]) for j in jobs)
    print(f"{df['foldx_ddg'].notna().sum()} ddG values ({df['foldx_ddg'].isna().sum()} failed) "
          f"of {expected} planned, saved to {OUT_PATH}")
    print(df.groupby("gene_symbol")["foldx_ddg"].describe()[["count", "mean", "50%"]].round(2).to_string())


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--work-dir", type=Path, default=WORK_DIR)
    p = sub.add_parser("run")
    p.add_argument("--foldx", type=Path, default=os.environ.get("FOLDX_BIN"), required="FOLDX_BIN" not in os.environ)
    p.add_argument("--work-dir", type=Path, default=WORK_DIR)
    p.add_argument("--workers", type=int, default=os.cpu_count())
    p.add_argument("--shard-weights", default="1", help="relative machine sizes, e.g. 32,16")
    p.add_argument("--shard", type=int, default=0, help="index into --shard-weights")
    p.add_argument("--chunks", nargs="*", help="limit to these chunks (testing)")
    p.add_argument("--max-jobs", type=int, help="limit the number of batches (testing)")
    p = sub.add_parser("package")
    p.add_argument("--work-dir", type=Path, default=WORK_DIR)
    p = sub.add_parser("collect")
    p.add_argument("--work-dir", type=Path, default=WORK_DIR)
    args = parser.parse_args()

    if args.cmd == "prepare":
        prepare(args.work_dir)
    elif args.cmd == "package":
        package(args.work_dir)
    elif args.cmd == "run":
        weights = [float(w) for w in args.shard_weights.split(",")]
        run(args.foldx, args.work_dir, args.workers, weights, args.shard, args.chunks, args.max_jobs)
    else:
        collect(args.work_dir)


if __name__ == "__main__":
    if __package__ is None:
        sys.path.insert(0, str(PROJECT_ROOT))
    main()
