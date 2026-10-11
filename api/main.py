"""Read-only FastAPI surface. Start from the repo root:

    uvicorn api.main:app --host 127.0.0.1 --port 8000
"""

from __future__ import annotations

import logging
import re
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Literal

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, PlainTextResponse

from api.load import get_catalog, load_catalog
from api.query import (
    ablation_payload,
    health_payload,
    list_variants,
    list_vus,
    metrics_payload,
    overview_payload,
    render_markdown,
    report_payload,
    resolve_pdb,
    resolve_variant,
    shap_global_payload,
    structure_payload,
    unavailable_detail,
    unavailable_report,
    unavailable_structure,
    variant_detail,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s %(message)s")


@asynccontextmanager
async def lifespan(app: FastAPI):
    app.state.catalog = load_catalog()
    yield


app = FastAPI(
    title="CVD variant interpretation API",
    version="0.1.0",
    summary="Read-only view of the cardiovascular variant pipeline.",
    lifespan=lifespan,
)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    allow_methods=["GET", "OPTIONS"],
    allow_headers=["*"],
)


def _catalog():
    try:
        return get_catalog()
    except RuntimeError as exc:
        raise HTTPException(status_code=503, detail={"value": None, "reason": str(exc)}) from exc


def _fail(status: int, reason: str) -> None:
    raise HTTPException(status_code=status, detail={"value": None, "reason": reason})


def _check_page(page: int, page_size: int) -> None:
    if page < 1:
        _fail(400, "page must be >= 1.")
    if page_size < 1 or page_size > 100:
        _fail(400, "page_size must be between 1 and 100.")


def _check_label(label: str | None) -> None:
    if label is not None and label not in {"pathogenic", "benign"}:
        _fail(400, "label must be pathogenic or benign.")


def _check_split(split: str | None) -> None:
    if split is not None and split not in {"train_gene", "test_gene"}:
        _fail(400, "split must be train_gene or test_gene.")


@app.get("/health")
def health():
    return health_payload(_catalog())


@app.get("/overview")
def overview():
    return overview_payload(_catalog())


@app.get("/metrics")
def metrics():
    return metrics_payload(_catalog())


@app.get("/shap/global")
def shap_global():
    return shap_global_payload(_catalog())


@app.get("/variants")
def variants(
    q: str | None = None,
    gene: str | None = None,
    label: str | None = None,
    split: str | None = None,
    page: int = Query(1),
    page_size: int = Query(25),
):
    _check_page(page, page_size)
    _check_label(label)
    _check_split(split)
    return list_variants(
        _catalog(),
        q=q or None,
        gene=gene or None,
        label=label,
        split=split,
        page=page,
        page_size=page_size,
    )


@app.get("/vus")
def vus(
    q: str | None = None,
    gene: str | None = None,
    label: str | None = None,
    split: str | None = None,
    page: int = Query(1),
    page_size: int = Query(25),
):
    # label and split are accepted and ignored. VUS rows are held out of the split.
    _check_page(page, page_size)
    return list_vus(
        _catalog(),
        q=q or None,
        gene=gene or None,
        page=page,
        page_size=page_size,
    )


@app.get("/variants/{variant_id}")
def variant(variant_id: str):
    catalog = _catalog()
    resolved = resolve_variant(catalog, variant_id)
    if resolved.status == 404:
        _fail(404, resolved.reason or "Variant was not found.")
    if resolved.record is None:
        return unavailable_detail(variant_id, resolved.reason or "Variant index is not loaded.")
    return variant_detail(catalog, resolved.record)


@app.get("/variants/{variant_id}/structure")
def variant_structure(variant_id: str):
    catalog = _catalog()
    resolved = resolve_variant(catalog, variant_id)
    if resolved.status == 404:
        _fail(404, resolved.reason or "Variant was not found.")
    if resolved.record is None:
        return unavailable_structure(variant_id, resolved.reason or "Variant index is not loaded.")
    return structure_payload(catalog, resolved.record)


@app.get("/variants/{variant_id}/report")
def variant_report(
    variant_id: str,
    format: Literal["json", "markdown"] = "json",
):
    catalog = _catalog()
    resolved = resolve_variant(catalog, variant_id)
    if resolved.status == 404:
        _fail(404, resolved.reason or "Variant was not found.")
    if resolved.record is None:
        body = unavailable_report(variant_id, resolved.reason or "Variant index is not loaded.")
    else:
        body = report_payload(catalog, resolved.record)
    if format == "markdown":
        return PlainTextResponse(render_markdown(body), media_type="text/markdown")
    return body


@app.get("/files/pdb/{name}")
def pdb_file(name: str):
    catalog = _catalog()
    if not re.fullmatch(r"^[A-Za-z0-9_.-]+\.pdb$", name):
        _fail(404, "PDB name is not on the outputs/*.pdb whitelist.")
    path = resolve_pdb(catalog, name)
    if path is None:
        _fail(404, f"outputs/{name} is not present.")
    return FileResponse(path, media_type="chemical/x-pdb", filename=name)


@app.get("/ablation")
def ablation():
    return ablation_payload(_catalog())


_DIST = Path(__file__).resolve().parents[1] / "web" / "dist"
_BUILD_HINT = (
    "The React dashboard has not been built.\n"
    "From the web folder run: npm install && npm run build\n"
    "Then restart this server and open http://127.0.0.1:8000/"
)


def _dist_file(relative: str) -> Path | None:
    if not relative or relative.endswith("/"):
        return None
    root = _DIST.resolve()
    candidate = (_DIST / relative).resolve()
    if candidate != root and root not in candidate.parents:
        return None
    if candidate.is_file():
        return candidate
    return None


@app.get("/", include_in_schema=False)
def dashboard():
    index = _DIST / "index.html"
    if index.is_file():
        return FileResponse(index)
    return PlainTextResponse(_BUILD_HINT, status_code=503)


@app.get("/{full_path:path}", include_in_schema=False)
def dashboard_files(full_path: str):
    asset = _dist_file(full_path)
    if asset is not None:
        return FileResponse(asset)
    index = _DIST / "index.html"
    if index.is_file():
        return FileResponse(index)
    _fail(404, "Frontend build is not present. From web/, run npm install && npm run build.")
