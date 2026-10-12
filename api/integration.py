"""Read-only access to the phase 2-6 pipeline outputs the original API predates.

Everything here reads files written by src.* stages; nothing is recomputed
except the integrated report, which calls src.report.generate_report.
Missing files give None so callers can return a reason instead of failing.
"""

from __future__ import annotations

import json
import math
from functools import lru_cache
from pathlib import Path
from typing import Any

ALPHAFOLD_CHAIN = "A"


def _paths() -> dict[str, Path]:
    from src.config import OUTPUTS_DIR, RAW_DIR, REPORTS_DIR

    return {
        "alphafold_dir": RAW_DIR / "alphafold",
        "alphafold_manifest": RAW_DIR / "alphafold_manifest.csv",
        "alphafold_residues": RAW_DIR / "alphafold_residue_features.csv",
        "foldx_ddg": RAW_DIR / "foldx_ddg.csv",
        "vus_report": REPORTS_DIR / "vus_report.csv",
        "metrics": OUTPUTS_DIR / "metrics.json",
    }


def artifact_specs() -> list[tuple[str, Path]]:
    """Files listed on /health next to the original artifacts."""
    p = _paths()
    return [(key, p[key]) for key in ("alphafold_manifest", "alphafold_residues", "foldx_ddg", "vus_report", "metrics")]


def _clean(value: Any) -> Any:
    if value is None:
        return None
    if hasattr(value, "item") and not isinstance(value, (str, bytes)):
        value = value.item()
    if isinstance(value, float) and math.isnan(value):
        return None
    return value


@lru_cache(maxsize=1)
def _manifest() -> dict[str, dict[str, Any]]:
    import pandas as pd

    path = _paths()["alphafold_manifest"]
    if not path.is_file():
        return {}
    return {row["gene_symbol"]: row for row in pd.read_csv(path).to_dict("records")}


@lru_cache(maxsize=1)
def _residues() -> dict[tuple[str, int], dict[str, Any]]:
    import pandas as pd

    path = _paths()["alphafold_residues"]
    if not path.is_file():
        return {}
    frame = pd.read_csv(path)
    return {(r["gene_symbol"], int(r["position"])): r for r in frame.to_dict("records")}


@lru_cache(maxsize=1)
def _ddg() -> dict[tuple[str, int, str, str], float]:
    import pandas as pd

    path = _paths()["foldx_ddg"]
    if not path.is_file():
        return {}
    frame = pd.read_csv(path).dropna(subset=["foldx_ddg"])
    return {(r.gene_symbol, int(r.position), r.wt_aa, r.mut_aa): float(r.foldx_ddg) for r in frame.itertuples()}


@lru_cache(maxsize=1)
def vus_ranking() -> dict[str, dict[str, Any]]:
    """variant_id ('MYH7_R403Q') -> row of outputs/reports/vus_report.csv."""
    import pandas as pd

    path = _paths()["vus_report"]
    if not path.is_file():
        return {}
    frame = pd.read_csv(path)
    return {row["variant"].replace(":", "_"): {k: _clean(v) for k, v in row.items()} for row in frame.to_dict("records")}


@lru_cache(maxsize=1)
def metrics() -> dict[str, Any] | None:
    path = _paths()["metrics"]
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def alphafold_model(gene: str) -> dict[str, Any] | None:
    """File name and source URL of the gene's AlphaFold DB model, if downloaded."""
    row = _manifest().get(gene)
    url = row and row.get("model_url")
    if not isinstance(url, str):
        return None
    name = url.rsplit("/", 1)[-1]
    if not (_paths()["alphafold_dir"] / name).is_file():
        return None
    return {"name": name, "source_url": url, "uniprot_id": row.get("uniprot_id")}


def alphafold_status(gene: str) -> str:
    row = _manifest().get(gene)
    if row is None:
        return "data/raw/alphafold_manifest.csv has no row for this gene."
    return str(row.get("status") or "No AlphaFold model was downloaded.")


def alphafold_path(name: str) -> Path | None:
    directory = _paths()["alphafold_dir"].resolve()
    path = (directory / name).resolve()
    if path.parent != directory or not path.is_file():
        return None
    return path


def residue(gene: str, position: int) -> dict[str, Any] | None:
    row = _residues().get((gene, int(position)))
    if row is None:
        return None
    return {
        "plddt": _clean(row.get("plddt")),
        "rsa": _clean(row.get("rsa")),
        "contact_count": _clean(row.get("contact_count")),
        "sec_struct": _clean(row.get("sec_struct")),
    }


def foldx_ddg(gene: str, position: int, wt: str, mut: str) -> float | None:
    return _ddg().get((gene, int(position), wt, mut))


def ddg_reason(gene: str, position: int) -> str:
    if not _paths()["foldx_ddg"].is_file():
        return "data/raw/foldx_ddg.csv is not present; run src/structure/foldx_ddg.py."
    if alphafold_model(gene) is None:
        return f"Not modelled: {alphafold_status(gene)}"
    ctx = residue(gene, position)
    if ctx and ctx["plddt"] is not None and ctx["plddt"] < 70:
        return f"Not modelled: residue pLDDT {ctx['plddt']:.0f} is below 70, so FoldX energies are not meaningful."
    return "FoldX did not return a value for this variant."


def ddg_category(ddg: float | None) -> str:
    from src.report.generate_report import ddg_category as category

    return category(ddg)


def integrated_report(variant_id: str) -> tuple[dict[str, Any] | None, str | None]:
    """Phase 6 report for a variant ID like 'MYH7_R403Q'."""
    from dataclasses import asdict

    try:
        from src.report.generate_report import report_for_variant

        gene, change = variant_id.split("_", 1)
        report = report_for_variant(f"{gene}:{change}")
    except Exception as exc:
        return None, f"Integrated report was not built: {exc}"
    body = asdict(report)
    body["markdown"] = report.to_markdown()
    return json.loads(json.dumps(body, default=_clean)), None
