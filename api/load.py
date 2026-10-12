"""Load pipeline artifacts once. Missing files are recorded; they do not raise."""

from __future__ import annotations

import json
import logging
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

logger = logging.getLogger("api")

STATUS_PATH = Path(__file__).resolve().parent / "pipeline_status.json"

WT_PDB_NAME = "ACTC1_R39H_wt.pdb"
MUTANT_PDB_NAME = "ACTC1_R39H_mutant.pdb"

SPLIT_CAVEAT = (
    "The feature matrix is absent. "
    "Group membership is computed from variants_labelled.csv with gene_level_split, "
    "which holds the same pathogenic/benign rows the models train on."
)


def _rel(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except (OSError, ValueError):
        return path.as_posix()


def _readable(path: Path) -> tuple[bool, str | None]:
    try:
        with path.open("rb") as handle:
            handle.read(1)
    except OSError as exc:
        return False, str(exc)
    return True, None


@dataclass
class SplitInfo:
    source: str
    caveat: str | None
    train_genes: list[str]
    test_genes: list[str]
    train_variants: int
    test_variants: int
    train_pathogenic: int
    train_benign: int
    test_pathogenic: int
    test_benign: int
    filtered_features: Any = None
    test_frame: Any = None


@dataclass
class Catalog:
    root: Path
    outputs_dir: Path
    artifacts: list[dict[str, Any]] = field(default_factory=list)
    genes: list[dict[str, Any]] = field(default_factory=list)
    labelled: dict[str, dict[str, Any]] | None = None
    labelled_reason: str | None = None
    vus: dict[str, dict[str, Any]] | None = None
    vus_reason: str | None = None
    labelled_counts: dict[str, dict[str, int]] | None = None
    vus_counts: dict[str, int] | None = None
    features: dict[str, dict[str, Any]] | None = None
    feature_frame: Any = None
    collected_count: int | None = None
    collected_reason: str | None = None
    uniprot: dict[str, dict[str, Any]] | None = None
    uniprot_reason: str | None = None
    conservation: dict[str, dict[str, float]] | None = None
    conservation_reason: str | None = None
    manifest: list[dict[str, Any]] | None = None
    manifest_reason: str | None = None
    models: dict[str, Any] = field(default_factory=lambda: {"rf": None, "xgb": None})
    model_reasons: dict[str, str | None] = field(default_factory=dict)
    feature_cols: list[str] | None = None
    model_feature_groups: tuple[list[str], list[str]] | None = None
    gene_level_split: Any = None
    training_import_error: str | None = None
    split: SplitInfo | None = None
    split_error: str | None = None
    gene_split: dict[str, str] = field(default_factory=dict)
    shap_module: Any = None
    shap_import_error: str | None = None
    explainers: dict[str, Any] = field(default_factory=dict)
    explainer_errors: dict[str, str] = field(default_factory=dict)
    global_shap: dict[str, Any] | None = None
    conservation_median: float | None = None
    conservation_median_reason: str | None = None
    pipeline: dict[str, Any] | None = None
    pipeline_reason: str | None = None
    pdb_wt_present: bool = False
    pdb_mut_present: bool = False

    def artifact(self, key: str) -> dict[str, Any] | None:
        for item in self.artifacts:
            if item["key"] == key:
                return item
        return None

    def indexes_complete(self) -> bool:
        return self.labelled is not None and self.vus is not None


_CATALOG: Catalog | None = None


def get_catalog() -> Catalog:
    if _CATALOG is None:
        raise RuntimeError("API catalog is not loaded.")
    return _CATALOG


def load_catalog() -> Catalog:
    global _CATALOG
    catalog = _build()
    _CATALOG = catalog
    loaded = [item["key"] for item in catalog.artifacts if item["loaded"]]
    missing = [item["key"] for item in catalog.artifacts if not item["present"]]
    logger.info("Artifacts loaded: %s", ", ".join(loaded) or "(none)")
    logger.info("Artifacts missing: %s", ", ".join(missing) or "(none)")
    return catalog


def _build() -> Catalog:
    try:
        from src.config import (
            FEATURES_DIR,
            OUTPUTS_DIR,
            PROCESSED_DIR,
            PROJECT_ROOT,
            RAW_DIR,
            load_gene_panel,
        )
    except Exception as exc:
        root = Path.cwd()
        catalog = Catalog(root=root, outputs_dir=root / "outputs")
        catalog.training_import_error = f"Cannot import src.config: {exc}"
        catalog.split_error = catalog.training_import_error
        catalog.collected_reason = catalog.training_import_error
        catalog.pipeline_reason = catalog.training_import_error
        return catalog

    catalog = Catalog(root=PROJECT_ROOT, outputs_dir=OUTPUTS_DIR)
    try:
        catalog.genes = load_gene_panel("both")
    except Exception as exc:
        catalog.genes = []
        logger.warning("Gene panel failed to load: %s", exc)

    specs = [
        ("clinvar_raw", RAW_DIR / "clinvar_raw.csv"),
        ("uniprot_entries", RAW_DIR / "uniprot_entries.json"),
        ("structure_manifest", RAW_DIR / "structure_manifest.csv"),
        ("conservation_scores", RAW_DIR / "conservation_scores.json"),
        ("variants_labelled", PROCESSED_DIR / "variants_labelled.csv"),
        ("variants_vus", PROCESSED_DIR / "variants_vus.csv"),
        ("feature_matrix", FEATURES_DIR / "feature_matrix.csv"),
        ("model_random_forest", OUTPUTS_DIR / "model_random_forest.joblib"),
        ("model_xgboost", OUTPUTS_DIR / "model_xgboost.joblib"),
        ("shap_global_plot", OUTPUTS_DIR / "shap_global_summary.png"),
        ("actc1_wt_pdb", OUTPUTS_DIR / WT_PDB_NAME),
        ("actc1_mutant_pdb", OUTPUTS_DIR / MUTANT_PDB_NAME),
    ]

    try:
        import pandas as pd
    except Exception as exc:
        pd = None
        pandas_error = f"pandas is not installed: {exc}"
    else:
        pandas_error = None

    from api.integration import artifact_specs

    specs += artifact_specs()
    for key, path in specs:
        catalog.artifacts.append(_load_one(catalog, key, path, pd, pandas_error))

    _import_training(catalog)
    _import_shap(catalog)
    _compute_split(catalog)
    _conservation_median(catalog)
    _load_pipeline_status(catalog)
    wt = catalog.artifact("actc1_wt_pdb")
    mut = catalog.artifact("actc1_mutant_pdb")
    catalog.pdb_wt_present = bool(wt and wt["present"])
    catalog.pdb_mut_present = bool(mut and mut["present"])
    return catalog


def _load_one(catalog: Catalog, key: str, path: Path, pd: Any, pandas_error: str | None) -> dict[str, Any]:
    rel = _rel(path, catalog.root)
    present = path.is_file()
    record = {"key": key, "path": rel, "present": present, "loaded": False, "reason": None}
    if not present:
        record["reason"] = f"{rel} is not present."
        _remember_missing(catalog, key, record["reason"])
        return record
    if pandas_error and key.endswith((".csv",)) or (pandas_error and key in {
        "clinvar_raw",
        "structure_manifest",
        "variants_labelled",
        "variants_vus",
        "feature_matrix",
    }):
        record["reason"] = pandas_error
        _remember_missing(catalog, key, pandas_error)
        return record
    try:
        _ingest(catalog, key, path, pd)
    except Exception as exc:
        record["reason"] = f"{rel} could not be read: {exc}"
        _remember_missing(catalog, key, record["reason"])
        logger.warning("%s", record["reason"])
        return record
    record["loaded"] = True
    return record


def _remember_missing(catalog: Catalog, key: str, reason: str) -> None:
    if key == "clinvar_raw":
        catalog.collected_reason = reason
    elif key == "uniprot_entries":
        catalog.uniprot_reason = reason
    elif key == "structure_manifest":
        catalog.manifest_reason = reason
    elif key == "conservation_scores":
        catalog.conservation_reason = reason
    elif key == "variants_labelled":
        catalog.labelled_reason = reason
    elif key == "variants_vus":
        catalog.vus_reason = reason
    elif key == "model_random_forest":
        catalog.model_reasons["rf"] = reason
    elif key == "model_xgboost":
        catalog.model_reasons["xgb"] = reason


def _ingest(catalog: Catalog, key: str, path: Path, pd: Any) -> None:
    if key == "clinvar_raw":
        catalog.collected_count = int(len(pd.read_csv(path)))
        catalog.collected_reason = None
        return
    if key == "variants_labelled":
        frame = pd.read_csv(path)
        catalog.labelled, catalog.labelled_counts = _index_variants(frame, "labelled")
        catalog.labelled_reason = None
        catalog._labelled_frame = frame  # type: ignore[attr-defined]
        return
    if key == "variants_vus":
        frame = pd.read_csv(path)
        catalog.vus, counts = _index_variants(frame, "vus")
        catalog.vus_counts = {gene: sum(labels.values()) for gene, labels in counts.items()}
        catalog.vus_reason = None
        return
    if key == "feature_matrix":
        frame = pd.read_csv(path)
        catalog.feature_frame = frame
        catalog.features = _index_features(frame)
        return
    if key == "uniprot_entries":
        payload = json.loads(path.read_text(encoding="utf-8"))
        catalog.uniprot = {}
        for entry in payload:
            symbol = entry.get("gene_symbol")
            if not symbol:
                continue
            catalog.uniprot[symbol] = {
                "uniprot_id": entry.get("uniprot_id"),
                "sequence": entry.get("sequence") or "",
                "sequence_length": entry.get("sequence_length") or len(entry.get("sequence") or ""),
                "features": list(entry.get("features") or []),
            }
        catalog.uniprot_reason = None
        return
    if key == "conservation_scores":
        payload = json.loads(path.read_text(encoding="utf-8"))
        cleaned: dict[str, dict[str, float]] = {}
        for gene, scores in payload.items():
            gene_scores: dict[str, float] = {}
            for pos, score in (scores or {}).items():
                try:
                    value = float(score)
                except (TypeError, ValueError):
                    continue
                if math.isnan(value):
                    continue
                gene_scores[str(pos)] = value
            cleaned[str(gene)] = gene_scores
        catalog.conservation = cleaned
        catalog.conservation_reason = None
        return
    if key == "structure_manifest":
        frame = pd.read_csv(path)
        rows = []
        for _, row in frame.iterrows():
            rows.append(
                {
                    "gene_symbol": None if pd.isna(row.get("gene_symbol")) else str(row.get("gene_symbol")),
                    "uniprot_id": None if pd.isna(row.get("uniprot_id")) else str(row.get("uniprot_id")),
                    "source": None if pd.isna(row.get("source")) else str(row.get("source")),
                    "structure_id": None if pd.isna(row.get("structure_id")) else str(row.get("structure_id")),
                }
            )
        catalog.manifest = rows
        catalog.manifest_reason = None
        return
    if key == "model_random_forest":
        catalog.models["rf"] = _load_joblib(path)
        catalog.model_reasons["rf"] = None
        return
    if key == "model_xgboost":
        catalog.models["xgb"] = _load_joblib(path)
        catalog.model_reasons["xgb"] = None
        return
    if key in {"shap_global_plot", "actc1_wt_pdb", "actc1_mutant_pdb", "alphafold_manifest",
               "alphafold_residues", "foldx_ddg", "vus_report", "metrics"}:
        ok, err = _readable(path)
        if not ok:
            raise OSError(err or "unreadable file")
        return
    raise RuntimeError(f"No loader for {key}")


def _load_joblib(path: Path) -> Any:
    import joblib

    return joblib.load(path)


def _na(value: Any) -> bool:
    if value is None:
        return True
    try:
        return bool(value != value)  # NaN
    except Exception:
        return False


def _text(value: Any) -> str | None:
    if _na(value):
        return None
    return str(value)


def _variant_id(gene: str, wt: str, position: int, mut: str) -> str:
    return f"{gene}_{wt}{position}{mut}"


def _index_variants(frame: Any, cohort: str) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, int]]]:
    records: list[dict[str, Any]] = []
    counts: dict[str, dict[str, int]] = {}
    for _, row in frame.iterrows():
        gene = _text(row.get("gene_symbol"))
        wt = _text(row.get("wt_aa"))
        mut = _text(row.get("mut_aa"))
        position_raw = row.get("position")
        if gene is None or wt is None or mut is None or _na(position_raw):
            continue
        position = int(position_raw)
        label = _text(row.get("label")) or cohort
        uid_raw = row.get("clinvar_uid")
        record = {
            "variant_id": _variant_id(gene, wt, position, mut),
            "cohort": "vus" if cohort == "vus" else "labelled",
            "gene_symbol": gene,
            "panel": _text(row.get("panel")),
            "wt_aa": wt,
            "mut_aa": mut,
            "position": position,
            "label": label,
            "clinical_significance": _text(row.get("clinical_significance")),
            "review_status": _text(row.get("review_status")),
            "protein_change": _text(row.get("protein_change")),
            "title": _text(row.get("title")),
            "clinvar_uid": None if _na(uid_raw) else int(uid_raw),
            "variation_set": _text(row.get("variation_set")),
            "molecular_consequence_list": _text(row.get("molecular_consequence_list")),
            "spdi": _text(row.get("spdi")),
        }
        records.append(record)
        bucket = counts.setdefault(gene, {})
        bucket[label] = bucket.get(label, 0) + 1
    records.sort(key=lambda item: (item["gene_symbol"], item["position"], item["mut_aa"], item["variant_id"]))
    index: dict[str, dict[str, Any]] = {}
    for record in records:
        index.setdefault(record["variant_id"], record)
    return index, counts


def _index_features(frame: Any) -> dict[str, dict[str, Any]]:
    index: dict[str, dict[str, Any]] = {}
    for _, row in frame.iterrows():
        gene = _text(row.get("gene_symbol"))
        wt = _text(row.get("wt_aa"))
        mut = _text(row.get("mut_aa"))
        position_raw = row.get("position")
        if gene is None or wt is None or mut is None or _na(position_raw):
            continue
        vid = _variant_id(gene, wt, int(position_raw), mut)
        index.setdefault(vid, {column: row[column] for column in frame.columns})
    return index


def _import_training(catalog: Catalog) -> None:
    try:
        from src.modeling.train_classifier import (
            BOOLEAN_FEATURES,
            CATEGORICAL_FEATURES,
            NUMERIC_FEATURES,
            gene_level_split,
        )
    except Exception as exc:
        catalog.training_import_error = f"Cannot import src.modeling.train_classifier: {exc}"
        catalog.model_feature_groups = None
        catalog.feature_cols = None
        catalog.gene_level_split = None
        return
    numeric = list(NUMERIC_FEATURES)
    boolean = list(BOOLEAN_FEATURES)
    categorical = list(CATEGORICAL_FEATURES)
    catalog.model_feature_groups = (numeric + boolean, categorical)
    catalog.feature_cols = numeric + boolean + categorical
    catalog.gene_level_split = gene_level_split


def _import_shap(catalog: Catalog) -> None:
    try:
        import shap
    except Exception as exc:
        catalog.shap_import_error = f"shap is not installed: {exc}"
        return
    catalog.shap_module = shap


def _label_count(frame: Any, label: str) -> int:
    if frame is None or "label" not in frame.columns:
        return 0
    return int((frame["label"] == label).sum())


def _compute_split(catalog: Catalog) -> None:
    if catalog.gene_level_split is None:
        catalog.split_error = catalog.training_import_error or "gene_level_split is unavailable."
        return
    try:
        if catalog.feature_frame is not None and catalog.feature_cols is not None:
            missing = [
                name for name in catalog.feature_cols + ["label", "gene_symbol"]
                if name not in catalog.feature_frame.columns
            ]
            if missing:
                catalog.split_error = "feature_matrix.csv is missing columns: " + ", ".join(missing)
                return
            # Same row selection as training (src.modeling.train_classifier.training_rows):
            # missing feature values are imputed in the pipeline, so no rows are dropped.
            from src.modeling.train_classifier import training_rows

            filtered = training_rows(catalog.feature_frame)
            train_df, test_df = catalog.gene_level_split(filtered)
            catalog.split = _split_info("feature_matrix.csv", None, train_df, test_df, filtered, test_df)
        elif getattr(catalog, "_labelled_frame", None) is not None:
            frame = catalog._labelled_frame  # type: ignore[attr-defined]
            train_df, test_df = catalog.gene_level_split(frame)
            catalog.split = _split_info("variants_labelled.csv", SPLIT_CAVEAT, train_df, test_df, None, None)
        else:
            catalog.split_error = (
                "Neither data/features/feature_matrix.csv nor "
                "data/processed/variants_labelled.csv is loaded."
            )
            return
    except Exception as exc:
        catalog.split_error = f"gene_level_split failed: {exc}"
        catalog.split = None
        return
    catalog.gene_split = {gene: "train_gene" for gene in catalog.split.train_genes}
    catalog.gene_split.update({gene: "test_gene" for gene in catalog.split.test_genes})


def _split_info(source: str, caveat: str | None, train_df: Any, test_df: Any, filtered: Any, test_frame: Any) -> SplitInfo:
    return SplitInfo(
        source=source,
        caveat=caveat,
        train_genes=sorted(train_df["gene_symbol"].astype(str).unique()),
        test_genes=sorted(test_df["gene_symbol"].astype(str).unique()),
        train_variants=int(len(train_df)),
        test_variants=int(len(test_df)),
        train_pathogenic=_label_count(train_df, "pathogenic"),
        train_benign=_label_count(train_df, "benign"),
        test_pathogenic=_label_count(test_df, "pathogenic"),
        test_benign=_label_count(test_df, "benign"),
        filtered_features=filtered,
        test_frame=test_frame,
    )


def _conservation_median(catalog: Catalog) -> None:
    groups = catalog.model_feature_groups
    if groups is None:
        catalog.conservation_median_reason = catalog.training_import_error
        return
    numeric_boolean, _categorical = groups
    try:
        index = numeric_boolean.index("conservation_score")
    except ValueError:
        catalog.conservation_median_reason = "conservation_score is not in the imported model feature list."
        return
    chosen = catalog.models.get("rf") or catalog.models.get("xgb")
    if chosen is None:
        catalog.conservation_median_reason = (
            "The training pipeline median-imputes a missing conservation score. "
            "No model is loaded, so that training-fold median is not available."
        )
        return
    try:
        imputer = chosen.named_steps["preprocess"].named_transformers_["num"]
        value = float(imputer.statistics_[index])
    except Exception as exc:
        catalog.conservation_median_reason = f"Could not read the fitted conservation median: {exc}"
        return
    if math.isnan(value):
        catalog.conservation_median_reason = "The fitted imputer has no finite median for conservation_score."
        return
    catalog.conservation_median = value


def _load_pipeline_status(catalog: Catalog) -> None:
    if not STATUS_PATH.is_file():
        catalog.pipeline_reason = "api/pipeline_status.json is not present."
        return
    try:
        catalog.pipeline = json.loads(STATUS_PATH.read_text(encoding="utf-8"))
    except Exception as exc:
        catalog.pipeline_reason = f"api/pipeline_status.json could not be read: {exc}"
