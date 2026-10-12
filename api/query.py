"""Endpoint payloads. Scientific steps call into src.*; they are not reimplemented."""

from __future__ import annotations

import ast
import json
import math
import re
import threading
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

from api.load import Catalog, get_catalog
from api.schemas import absent, maybe

POC_VARIANT_ID = "ACTC1_R39H"
WT_PDB_NAME = "ACTC1_R39H_wt.pdb"
MUTANT_PDB_NAME = "ACTC1_R39H_mutant.pdb"
PDB_NAME = re.compile(r"^[A-Za-z0-9_.-]+\.pdb$")

DISCLAIMER = "Research prototype — not for clinical use"
METHOD_NOTE = (
    "The mutant file is a residue-identity swap: ARG 39 is renamed to HIS on chain A "
    "and the coordinates are unchanged (src/structure/mutate_pdb.py). No side-chain rebuild "
    "and no minimization. A Cα RMSD of 0.00 Å is a method artifact, not evidence that the "
    "substitution has no structural effect."
)
ESMFOLD_NOTE = (
    "src/structure/predict_mutant_structure.py can call the public ESMFold endpoint. "
    "That endpoint timed out during the project. This API does not call it. "
    "The proof-of-concept PDBs are not ESMFold models."
)
DDG_REASON = (
    "FoldX BuildModel is not implemented in src/structure/stability_analysis.py. "
    "ΔΔG was not computed."
)
OTHER_STRUCTURE_REASON = "Structural proof-of-concept is implemented only for ACTC1_R39H."
OTHER_STABILITY = "Structural comparison was not run for this variant. ΔΔG was not computed."
KNOWN_STRUCTURE_NOTE = (
    "fetch_structures keeps a single structure per gene (first RCSB hit, else AlphaFold DB). "
    "This is not a complete PDB list."
)
CONSERVATION_METHOD = (
    "Homolog-identity proxy from src/features/build_conservation_scores.py. "
    "Not ConSurf, phyloP, or phastCons."
)
SHAP_NOTE = (
    "Mean |SHAP| on the labelled matrix after the training row filter, train and test genes together. "
    "The saved PNG is XGBoost-only and drops rows with a missing conservation score."
)
REPORTED_SOURCE = "README.md and app.py ML page literals; not recomputed"
ABLATION_REASON = (
    "No ablation has been run. outputs/ablation_results.json is not present. "
    "This endpoint does not train."
)
ABLATION_PLAN = {
    "sets": [
        "A_physchem",
        "B_physchem_conservation",
        "C_physchem_conservation_domain",
    ],
    "position": ["with_position", "without_position"],
    "note": (
        "The current model also includes panel, which is constant per gene. "
        "domain_type is not a model input; the domain feature in set C is in_annotated_domain. "
        "Gene-level split and class weighting stay as in train_classifier.py."
    ),
}
REPORTED = {
    "source": REPORTED_SOURCE,
    "random_forest": {
        "accuracy": 0.622,
        "precision": 0.647,
        "recall": 0.887,
        "f1": 0.748,
        "roc_auc": 0.588,
        "pr_auc": 0.706,
        "confusion_matrix": [[6, 30], [7, 55]],
    },
    "xgboost": {
        "accuracy": 0.633,
        "precision": 0.686,
        "recall": 0.774,
        "f1": 0.727,
        "roc_auc": 0.559,
        "pr_auc": 0.703,
        "confusion_matrix": [[14, 22], [14, 48]],
    },
}

PHYS_KEYS = [
    "delta_hydrophobicity",
    "delta_molecular_weight",
    "delta_volume",
    "delta_charge",
    "polarity_change",
    "wt_is_proline",
    "mut_is_proline",
    "wt_is_glycine",
    "mut_is_glycine",
]
PHYS_BOOLS = {
    "polarity_change",
    "wt_is_proline",
    "mut_is_proline",
    "wt_is_glycine",
    "mut_is_glycine",
}
MODEL_KEYS = (("rf", "random_forest", "Random Forest"), ("xgb", "xgboost", "XGBoost"))
_SHAP_LOCK = threading.Lock()


def _missing(value: Any) -> bool:
    if value is None:
        return True
    try:
        return bool(value != value)
    except Exception:
        return False


def _py_scalar(value: Any) -> Any:
    if _missing(value):
        return None
    if hasattr(value, "item") and not isinstance(value, (str, bytes)):
        try:
            value = value.item()
        except Exception:
            pass
    if isinstance(value, float) and math.isnan(value):
        return None
    return value


def _as_bool(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "t", "yes"}
    return bool(value)


def _finite(value: Any) -> float | None:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def _prediction(label: str | None, prob: float | None, reason: str | None) -> dict[str, Any]:
    return {
        "label": label,
        "prob_pathogenic": prob,
        "threshold": 0.5,
        "reason": None if prob is not None else (reason or "Prediction was not computed."),
    }


def _empty_shap(reason: str) -> dict[str, Any]:
    return {
        "shap_space": "probability",
        "base_value": None,
        "output_value": None,
        "features": [],
        "reason": reason,
    }


def _legit_null() -> dict[str, Any]:
    return {"value": None, "reason": None}


@dataclass
class Resolved:
    status: int
    record: dict[str, Any] | None
    reason: str | None


def resolve_variant(catalog: Catalog, variant_id: str) -> Resolved:
    if catalog.labelled is not None and variant_id in catalog.labelled:
        return Resolved(200, catalog.labelled[variant_id], None)
    if catalog.vus is not None and variant_id in catalog.vus:
        return Resolved(200, catalog.vus[variant_id], None)
    if catalog.indexes_complete():
        return Resolved(404, None, f"No labelled or VUS row has id {variant_id}")
    parts: list[str] = []
    if catalog.labelled is None:
        parts.append(catalog.labelled_reason or "The labelled index is not loaded.")
    else:
        parts.append(f"{variant_id} is not in the labelled table.")
    if catalog.vus is None:
        parts.append(catalog.vus_reason or "The VUS index is not loaded.")
    else:
        parts.append(f"{variant_id} is not in the VUS table.")
    return Resolved(200, None, " ".join(parts))


def health_payload(catalog: Catalog) -> dict[str, Any]:
    return {"status": "ok", "artifacts": catalog.artifacts}


def overview_payload(catalog: Catalog) -> dict[str, Any]:
    labelled_reason = catalog.labelled_reason or "data/processed/variants_labelled.csv is not loaded."
    vus_reason = catalog.vus_reason or "data/processed/variants_vus.csv is not loaded."
    if catalog.labelled is None:
        labelled = pathogenic = benign = absent(labelled_reason)
    else:
        labelled = maybe(len(catalog.labelled))
        pathogenic = maybe(sum(1 for row in catalog.labelled.values() if row["label"] == "pathogenic"))
        benign = maybe(sum(1 for row in catalog.labelled.values() if row["label"] == "benign"))
    vus = absent(vus_reason) if catalog.vus is None else maybe(len(catalog.vus))
    collected = (
        maybe(catalog.collected_count)
        if catalog.collected_count is not None
        else absent(catalog.collected_reason or "data/raw/clinvar_raw.csv is not present.")
    )

    genes = []
    for gene in catalog.genes:
        symbol = gene["symbol"]
        if catalog.labelled_counts is None:
            benign_n = pathogenic_n = absent(labelled_reason)
        else:
            bucket = catalog.labelled_counts.get(symbol, {})
            benign_n = maybe(int(bucket.get("benign", 0)))
            pathogenic_n = maybe(int(bucket.get("pathogenic", 0)))
        if catalog.vus_counts is None:
            vus_n = absent(vus_reason)
        else:
            vus_n = maybe(int(catalog.vus_counts.get(symbol, 0)))
        genes.append(
            {
                "symbol": symbol,
                "panel": gene.get("panel"),
                "condition": gene.get("condition"),
                "uniprot_id": gene.get("uniprot_id"),
                "benign": benign_n,
                "pathogenic": pathogenic_n,
                "vus": vus_n,
            }
        )

    if catalog.manifest is None:
        structure_records = absent(catalog.manifest_reason or "data/raw/structure_manifest.csv is not present.")
    else:
        structure_records = maybe(len(catalog.manifest))

    split_error = catalog.split_error or "The gene-level split could not be computed."
    if catalog.split is None:
        split = {
            "source": None,
            "caveat": None,
            "train_genes": absent(split_error),
            "test_genes": absent(split_error),
            "train_variants": absent(split_error),
            "test_variants": absent(split_error),
            "train_pathogenic": absent(split_error),
            "train_benign": absent(split_error),
            "test_pathogenic": absent(split_error),
            "test_benign": absent(split_error),
        }
    else:
        info = catalog.split
        split = {
            "source": info.source,
            "caveat": info.caveat,
            "train_genes": maybe(info.train_genes),
            "test_genes": maybe(info.test_genes),
            "train_variants": maybe(info.train_variants),
            "test_variants": maybe(info.test_variants),
            "train_pathogenic": maybe(info.train_pathogenic),
            "train_benign": maybe(info.train_benign),
            "test_pathogenic": maybe(info.test_pathogenic),
            "test_benign": maybe(info.test_benign),
        }

    matrix_reason = _artifact_reason(catalog, "feature_matrix", "data/features/feature_matrix.csv is not loaded.")
    frame = catalog.feature_frame
    if frame is None or "conservation_score" not in getattr(frame, "columns", []):
        scored = missing_scores = absent(matrix_reason)
    else:
        missing_n = int(frame["conservation_score"].isna().sum())
        scored = maybe(int(len(frame) - missing_n))
        missing_scores = maybe(missing_n)

    return {
        "counts": {
            "collected": collected,
            "labelled": labelled,
            "pathogenic": pathogenic,
            "benign": benign,
            "vus": vus,
        },
        "genes": genes,
        "structure_records": structure_records,
        "split": split,
        "conservation": {
            "scored_variants": scored,
            "missing_variants": missing_scores,
            "method": "homolog_identity_proxy",
        },
        "pipeline": catalog.pipeline,
        "pipeline_reason": catalog.pipeline_reason,
    }


def _artifact_reason(catalog: Catalog, key: str, fallback: str) -> str:
    artifact = catalog.artifact(key)
    if artifact and artifact.get("reason"):
        return str(artifact["reason"])
    return fallback


def metrics_payload(catalog: Catalog) -> dict[str, Any]:
    info = catalog.split
    if info is None:
        split_block = {
            "test_genes": absent(catalog.split_error or "The gene-level split could not be computed."),
            "n_test": absent(catalog.split_error or "The gene-level split could not be computed."),
            "source": None,
        }
    else:
        split_block = {
            "test_genes": maybe(info.test_genes),
            "n_test": maybe(info.test_variants),
            "source": info.source,
        }

    blockers: list[str] = []
    if catalog.feature_frame is None or info is None or info.test_frame is None:
        blockers.append(
            _artifact_reason(catalog, "feature_matrix", "data/features/feature_matrix.csv is not loaded.")
        )
    if catalog.models.get("rf") is None:
        blockers.append(catalog.model_reasons.get("rf") or "outputs/model_random_forest.joblib is not loaded.")
    if catalog.models.get("xgb") is None:
        blockers.append(catalog.model_reasons.get("xgb") or "outputs/model_xgboost.joblib is not loaded.")

    payload: dict[str, Any] = {
        "split": split_block,
        "reported": REPORTED,
        "recomputed": None,
        "comparison": None,
        "reason": None,
    }
    if blockers:
        payload["reason"] = " ".join(blockers)
        return payload

    try:
        recomputed = {
            "random_forest": _evaluate_model(catalog, "rf", "RandomForest"),
            "xgboost": _evaluate_model(catalog, "xgb", "XGBoost"),
        }
    except Exception as exc:
        payload["reason"] = f"Metrics could not be recomputed: {exc}"
        return payload
    payload["recomputed"] = recomputed
    payload["comparison"] = {
        "random_forest": _compare(recomputed["random_forest"], REPORTED["random_forest"]),
        "xgboost": _compare(recomputed["xgboost"], REPORTED["xgboost"]),
    }
    return payload


def _evaluate_model(catalog: Catalog, model_key: str, display_name: str) -> dict[str, Any]:
    from sklearn.metrics import (
        accuracy_score,
        average_precision_score,
        confusion_matrix,
        f1_score,
        precision_recall_curve,
        precision_score,
        recall_score,
        roc_auc_score,
        roc_curve,
    )

    from src.modeling.evaluate import evaluate_classifier

    model = catalog.models[model_key]
    frame = catalog.split.test_frame
    feature_cols = catalog.feature_cols or []
    x_test = frame[feature_cols]
    y_test = (frame["label"] == "pathogenic").astype(int)
    y_pred = model.predict(x_test)
    y_proba = model.predict_proba(x_test)[:, 1]
    reason_parts: list[str] = []
    try:
        scalars = {key: float(value) for key, value in evaluate_classifier(model, x_test, y_test, name=display_name).items()}
    except Exception as exc:
        reason_parts.append(f"evaluate_classifier failed: {exc}")
        scalars = {
            "accuracy": float(accuracy_score(y_test, y_pred)),
            "precision": float(precision_score(y_test, y_pred, zero_division=0)),
            "recall": float(recall_score(y_test, y_pred, zero_division=0)),
            "f1": float(f1_score(y_test, y_pred, zero_division=0)),
            "roc_auc": None,
            "pr_auc": None,
        }
        try:
            scalars["roc_auc"] = float(roc_auc_score(y_test, y_proba))
        except Exception as roc_exc:
            reason_parts.append(f"roc_auc unavailable: {roc_exc}")
        try:
            scalars["pr_auc"] = float(average_precision_score(y_test, y_proba))
        except Exception as pr_exc:
            reason_parts.append(f"pr_auc unavailable: {pr_exc}")

    matrix = confusion_matrix(y_test, y_pred, labels=[0, 1]).tolist()
    matrix = [[int(value) for value in row] for row in matrix]
    try:
        fpr, tpr, roc_thresholds = roc_curve(y_test, y_proba)
        roc = [
            {"fpr": _finite(a), "tpr": _finite(b), "threshold": _finite(c)}
            for a, b, c in zip(fpr, tpr, roc_thresholds)
        ]
    except Exception as exc:
        roc = []
        reason_parts.append(f"roc_curve failed: {exc}")
    try:
        precision, recall, pr_thresholds = precision_recall_curve(y_test, y_proba)
        pr = []
        for index, (prec, rec) in enumerate(zip(precision, recall)):
            threshold = pr_thresholds[index] if index < len(pr_thresholds) else None
            pr.append(
                {
                    "precision": _finite(prec),
                    "recall": _finite(rec),
                    "threshold": _finite(threshold) if threshold is not None else None,
                }
            )
    except Exception as exc:
        pr = []
        reason_parts.append(f"precision_recall_curve failed: {exc}")

    return {
        "accuracy": scalars["accuracy"],
        "precision": scalars["precision"],
        "recall": scalars["recall"],
        "f1": scalars["f1"],
        "roc_auc": scalars["roc_auc"],
        "pr_auc": scalars["pr_auc"],
        "confusion_matrix": matrix,
        "roc": roc,
        "pr": pr,
        "reason": " ".join(reason_parts) or None,
    }


def _compare(got: dict[str, Any], expected: dict[str, Any]) -> dict[str, Any]:
    flags = {}
    for key in ("accuracy", "precision", "recall", "f1", "roc_auc", "pr_auc"):
        actual = got.get(key)
        flags[key] = actual is not None and round(float(actual), 3) == round(float(expected[key]), 3)
    flags["confusion_matrix"] = got.get("confusion_matrix") == expected["confusion_matrix"]
    flags["all_matched"] = all(flags.values())
    return flags


def shap_global_payload(catalog: Catalog) -> dict[str, Any]:
    if catalog.global_shap is not None:
        return catalog.global_shap
    frame, frame_reason = _global_frame(catalog)
    if frame is None:
        row_count = absent(frame_reason)
        imputed = absent(frame_reason)
    else:
        row_count = maybe(int(len(frame)))
        if "conservation_score" in frame.columns:
            imputed = maybe(int(frame["conservation_score"].isna().sum()))
        else:
            imputed = absent("feature_matrix.csv has no conservation_score column.")
    payload = {
        "shap_space": "probability",
        "row_count": row_count,
        "imputed_conservation_rows": imputed,
        "note": SHAP_NOTE,
        "random_forest": _global_for_model(catalog, "rf", frame, frame_reason),
        "xgboost": _global_for_model(catalog, "xgb", frame, frame_reason),
    }
    catalog.global_shap = payload
    return payload


def _global_frame(catalog: Catalog) -> tuple[Any, str]:
    if catalog.split is not None and catalog.split.filtered_features is not None:
        return catalog.split.filtered_features, ""
    return None, _artifact_reason(
        catalog,
        "feature_matrix",
        catalog.split_error or "data/features/feature_matrix.csv is not loaded.",
    )


def _global_for_model(catalog: Catalog, model_key: str, frame: Any, frame_reason: str) -> dict[str, Any]:
    model = catalog.models.get(model_key)
    if model is None:
        return absent(catalog.model_reasons.get(model_key) or "Model is not loaded.")
    if frame is None:
        return absent(frame_reason)
    if catalog.shap_module is None:
        return absent(catalog.shap_import_error or "shap is not installed.")
    try:
        values, names, _base = _explain_matrix(catalog, model_key, model, frame[catalog.feature_cols])
    except Exception as exc:
        return absent(f"Global SHAP failed: {exc}")
    means = values.__abs__().mean(axis=0)
    rows = []
    for index, name in enumerate(names):
        rows.append(
            {
                "feature": name,
                "raw_feature": _raw_feature_name(name),
                "mean_abs_shap": float(means[index]),
            }
        )
    rows.sort(key=lambda item: (-item["mean_abs_shap"], item["feature"]))
    return maybe(rows)


def _explain_matrix(catalog: Catalog, model_key: str, model: Any, frame: Any) -> tuple[Any, list[str], Any]:
    import numpy as np

    preprocess = model.named_steps["preprocess"]
    transformed = preprocess.transform(frame)
    if hasattr(transformed, "toarray"):
        transformed = transformed.toarray()
    transformed = np.asarray(transformed)
    names = [str(name) for name in preprocess.get_feature_names_out()]
    with _SHAP_LOCK:
        explainer = _explainer(catalog, model_key, model)
        explanation = explainer(transformed)
    values, base = _pathogenic_slice(explanation)
    if values.shape[0] != transformed.shape[0]:
        raise RuntimeError(
            f"SHAP row count {values.shape[0]} does not match the transformed matrix {transformed.shape[0]}."
        )
    if values.shape[-1] != len(names):
        raise RuntimeError(
            f"SHAP feature count {values.shape[-1]} does not match transformed names {len(names)}."
        )
    return values, names, base


SHAP_BACKGROUND_ROWS = 200


def _background(catalog: Catalog, model: Any) -> Any:
    if catalog.split is None or catalog.feature_frame is None or catalog.feature_cols is None:
        raise RuntimeError("No training rows available for the SHAP background.")
    from src.modeling.train_classifier import training_rows

    frame = training_rows(catalog.feature_frame)
    train_genes = set(catalog.split.train_genes)
    rows = frame[frame["gene_symbol"].isin(train_genes)]
    rows = rows.sample(min(SHAP_BACKGROUND_ROWS, len(rows)), random_state=0)
    data = model.named_steps["preprocess"].transform(rows[catalog.feature_cols])
    return data.toarray() if hasattr(data, "toarray") else data


def _explainer(catalog: Catalog, model_key: str, model: Any) -> Any:
    if model_key in catalog.explainer_errors:
        raise RuntimeError(catalog.explainer_errors[model_key])
    cached = catalog.explainers.get(model_key)
    if cached is not None:
        return cached
    try:
        # model_output="probability" puts the waterfall in P(pathogenic). SHAP
        # only supports it with background data (interventional), so use a
        # fixed sample of the training-gene rows the split models were fit on.
        explainer = catalog.shap_module.TreeExplainer(
            model.named_steps["clf"],
            data=_background(catalog, model),
            model_output="probability",
            feature_perturbation="interventional",
        )
    except Exception as exc:
        catalog.explainer_errors[model_key] = str(exc)
        raise
    catalog.explainers[model_key] = explainer
    return explainer


def _pathogenic_slice(explanation: Any) -> tuple[Any, Any]:
    import numpy as np

    values = np.asarray(explanation.values)
    base = np.asarray(explanation.base_values)
    class_index = None
    if values.ndim == 3:
        # Pathogenic class is index 1 when both classes are present.
        class_index = 1 if values.shape[-1] >= 2 else 0
        values = values[:, :, class_index]
    elif values.ndim == 1:
        values = values.reshape(1, -1)
    if class_index is not None and base.ndim >= 1 and base.shape[-1] >= 2 and base.shape[-1] == (2 if class_index == 1 else 1):
        base = np.take(base, class_index, axis=-1)
    if base.ndim == 0:
        base = np.array([float(base)])
    return values, base


def _raw_feature_name(transformed: str) -> str:
    if transformed.startswith("num__"):
        return transformed[len("num__") :]
    if transformed.startswith("cat__"):
        rest = transformed[len("cat__") :]
        if rest.startswith("panel_"):
            return "panel"
        return rest
    return transformed


def list_variants(
    catalog: Catalog,
    *,
    q: str | None,
    gene: str | None,
    label: str | None,
    split: str | None,
    page: int,
    page_size: int,
) -> dict[str, Any]:
    if catalog.labelled is None:
        return _empty_page(page, page_size, catalog.labelled_reason or "data/processed/variants_labelled.csv is not loaded.")
    if split and not catalog.gene_split:
        return _empty_page(page, page_size, catalog.split_error or "The gene-level split could not be computed.")
    rows = [
        row
        for row in catalog.labelled.values()
        if _matches(row, q=q, gene=gene, label=label, split=split, catalog=catalog)
    ]
    return _page(rows, page, page_size, _summary)


def list_vus(
    catalog: Catalog,
    *,
    q: str | None,
    gene: str | None,
    page: int,
    page_size: int,
) -> dict[str, Any]:
    if catalog.vus is None:
        return _empty_page(page, page_size, catalog.vus_reason or "data/processed/variants_vus.csv is not loaded.")
    rows = [row for row in catalog.vus.values() if _matches(row, q=q, gene=gene, label=None, split=None, catalog=catalog)]
    page_rows, total = _slice(rows, page, page_size)
    return {
        "page": page,
        "page_size": page_size,
        "total": total,
        "items": [_vus_summary(catalog, row) for row in page_rows],
        "reason": None,
    }


def _empty_page(page: int, page_size: int, reason: str) -> dict[str, Any]:
    return {"page": page, "page_size": page_size, "total": 0, "items": [], "reason": reason}


def _slice(rows: list[dict[str, Any]], page: int, page_size: int) -> tuple[list[dict[str, Any]], int]:
    start = (page - 1) * page_size
    return rows[start : start + page_size], len(rows)


def _page(rows: list[dict[str, Any]], page: int, page_size: int, mapper) -> dict[str, Any]:
    chunk, total = _slice(rows, page, page_size)
    return {
        "page": page,
        "page_size": page_size,
        "total": total,
        "items": [mapper(row) for row in chunk],
        "reason": None,
    }


def _matches(
    row: dict[str, Any],
    *,
    q: str | None,
    gene: str | None,
    label: str | None,
    split: str | None,
    catalog: Catalog,
) -> bool:
    if gene and row["gene_symbol"].casefold() != gene.casefold():
        return False
    if label and row["label"] != label:
        return False
    if split:
        if catalog.gene_split.get(row["gene_symbol"]) != split:
            return False
    if q:
        needle = q.casefold()
        haystack = " ".join(
            [
                row["variant_id"],
                row["gene_symbol"],
                row.get("protein_change") or "",
                row.get("title") or "",
            ]
        ).casefold()
        if needle not in haystack:
            return False
    return True


def _summary(row: dict[str, Any]) -> dict[str, Any]:
    catalog = get_catalog()
    split_name = catalog.gene_split.get(row["gene_symbol"])
    return {
        "variant_id": row["variant_id"],
        "gene_symbol": row["gene_symbol"],
        "panel": row.get("panel"),
        "wt_aa": row["wt_aa"],
        "mut_aa": row["mut_aa"],
        "position": row["position"],
        "label": row["label"],
        "clinical_significance": row.get("clinical_significance"),
        "protein_change": row.get("protein_change"),
        "split": split_name,
        "split_reason": None if split_name else (catalog.split_error or "The gene-level split could not be computed."),
        "has_structure": _has_structure(catalog, row["variant_id"]),
    }


def _vus_summary(catalog: Catalog, row: dict[str, Any]) -> dict[str, Any]:
    inputs = build_inputs(catalog, row)
    predictions, predictions_reason = _vus_predictions(catalog, inputs)
    return {
        "variant_id": row["variant_id"],
        "gene_symbol": row["gene_symbol"],
        "panel": row.get("panel"),
        "wt_aa": row["wt_aa"],
        "mut_aa": row["mut_aa"],
        "position": row["position"],
        "label": "vus",
        "clinical_significance": row.get("clinical_significance"),
        "protein_change": row.get("protein_change"),
        "predictions": predictions,
        "predictions_reason": predictions_reason,
    }


def _vus_predictions(catalog: Catalog, inputs: "Inputs") -> tuple[dict[str, Any] | None, str | None]:
    rf_loaded = catalog.models.get("rf") is not None
    xgb_loaded = catalog.models.get("xgb") is not None
    if inputs.frame is None:
        return None, inputs.frame_reason
    if not rf_loaded and not xgb_loaded:
        return None, " ".join(
            reason
            for reason in (
                catalog.model_reasons.get("rf"),
                catalog.model_reasons.get("xgb"),
            )
            if reason
        ) or "Neither model joblib is loaded."
    return {
        "rf": _predict(catalog, "rf", inputs.frame, inputs.frame_reason),
        "xgb": _predict(catalog, "xgb", inputs.frame, inputs.frame_reason),
    }, None


def _has_structure(catalog: Catalog, variant_id: str) -> bool:
    return variant_id == POC_VARIANT_ID and catalog.pdb_wt_present and catalog.pdb_mut_present


@dataclass
class Inputs:
    physchem: dict[str, Any]
    conservation_score: dict[str, Any]
    conservation_imputed: bool
    conservation_imputed_value: dict[str, Any]
    in_annotated_domain: dict[str, Any]
    domain_type: dict[str, Any]
    frame: Any
    frame_reason: str | None
    clinical: dict[str, Any]
    wt_three: str | None
    mut_three: str | None
    hgvs: str | None
    sequence: dict[str, Any]


def build_inputs(catalog: Catalog, record: dict[str, Any]) -> Inputs:
    feature_row = (catalog.features or {}).get(record["variant_id"])
    physchem, frame_bits = _physchem(record, feature_row)
    conservation_score, raw_score = _conservation(catalog, record, feature_row)
    in_domain, domain_type, domain_known = _domain(catalog, record, feature_row)
    frame, frame_reason = _frame(catalog, record, physchem, raw_score, in_domain, domain_known)
    if frame_bits and frame_reason is None:
        frame_reason = frame_bits
    imputed, imputed_value = _imputation(catalog, raw_score)
    return Inputs(
        physchem=physchem,
        conservation_score=conservation_score,
        conservation_imputed=imputed,
        conservation_imputed_value=imputed_value,
        in_annotated_domain=in_domain,
        domain_type=domain_type,
        frame=frame,
        frame_reason=frame_reason,
        clinical=_clinical(record),
        wt_three=_three(record["wt_aa"]),
        mut_three=_three(record["mut_aa"]),
        hgvs=_hgvs(record),
        sequence=_sequence(catalog, record),
    )


def _physchem(record: dict[str, Any], feature_row: dict[str, Any] | None) -> tuple[dict[str, Any], str | None]:
    if feature_row is not None and all(key in feature_row for key in PHYS_KEYS):
        values = {}
        for key in PHYS_KEYS:
            raw = _py_scalar(feature_row[key])
            values[key] = _as_bool(raw) if key in PHYS_BOOLS and raw is not None else raw
        if all(values[key] is not None for key in PHYS_KEYS):
            values["source"] = "feature_matrix.csv"
            values["reason"] = None
            return values, None
    try:
        from src.features.physicochemical import compute_physchem_features

        computed = compute_physchem_features(record["wt_aa"], record["mut_aa"])
    except Exception as exc:
        blank = {key: None for key in PHYS_KEYS}
        blank["source"] = None
        blank["reason"] = f"compute_physchem_features failed: {exc}"
        return blank, blank["reason"]
    values = {key: getattr(computed, key) for key in PHYS_KEYS}
    values["source"] = "src.features.physicochemical.compute_physchem_features"
    values["reason"] = None
    return values, None


def _conservation(
    catalog: Catalog,
    record: dict[str, Any],
    feature_row: dict[str, Any] | None,
) -> tuple[dict[str, Any], float | None]:
    if feature_row is not None and "conservation_score" in feature_row:
        raw = _py_scalar(feature_row["conservation_score"])
        if raw is None:
            return absent(_score_gap_reason(catalog, record, from_matrix=True)), None
        return maybe(float(raw)), float(raw)
    if catalog.conservation is not None:
        from src.features.conservation_and_domain import get_conservation_score

        score = get_conservation_score(record["gene_symbol"], int(record["position"]), catalog.conservation)
        if score is None:
            return absent(_score_gap_reason(catalog, record, from_matrix=False)), None
        return maybe(float(score)), float(score)
    return absent(
        "data/raw/conservation_scores.json is not present and "
        "data/features/feature_matrix.csv has no row for this variant."
    ), None


def _score_gap_reason(catalog: Catalog, record: dict[str, Any], *, from_matrix: bool) -> str:
    where = (
        f"No homolog-identity score is stored for {record['gene_symbol']} position {record['position']}."
    )
    if from_matrix and catalog.conservation is None:
        return (
            "conservation_score is empty on the feature-matrix row. "
            "data/raw/conservation_scores.json is not present."
        )
    if from_matrix:
        return where + " The feature-matrix cell is empty."
    return where


def _domain(
    catalog: Catalog,
    record: dict[str, Any],
    feature_row: dict[str, Any] | None,
) -> tuple[dict[str, Any], dict[str, Any], bool]:
    if feature_row is not None and "in_annotated_domain" in feature_row:
        in_domain = _as_bool(_py_scalar(feature_row["in_annotated_domain"]))
        domain_raw = _py_scalar(feature_row["domain_type"]) if "domain_type" in feature_row else None
        domain_value = None if domain_raw is None else str(domain_raw)
        return maybe(in_domain), {"value": domain_value, "reason": None}, True
    if catalog.uniprot is not None:
        from src.features.conservation_and_domain import compute_domain_membership

        features = catalog.uniprot.get(record["gene_symbol"], {}).get("features", [])
        in_domain, domain_type = compute_domain_membership(int(record["position"]), features)
        return maybe(bool(in_domain)), {"value": domain_type, "reason": None}, True
    reason = (
        "data/raw/uniprot_entries.json is not present, so domain membership was not checked. "
        "An empty annotation list would make the feature builder store false; "
        "this API does not treat unknown as false."
    )
    return absent(reason), absent(reason), False


def _frame(
    catalog: Catalog,
    record: dict[str, Any],
    physchem: dict[str, Any],
    raw_score: float | None,
    in_domain: dict[str, Any],
    domain_known: bool,
) -> tuple[Any, str | None]:
    if catalog.feature_cols is None:
        return None, catalog.training_import_error or "Model feature columns could not be imported."
    if physchem.get("reason"):
        return None, physchem["reason"]
    if not domain_known:
        return None, in_domain.get("reason") or "Domain membership is unknown."
    if not record.get("panel"):
        return None, "panel is missing on this row, so the model frame cannot be built."
    import pandas as pd

    # The feature set grew in pipeline phase 2 (VEP conservation, AlphaFold
    # context, UniProt annotations), so the row is built by the same src
    # function that builds the training matrix rather than by hand here.
    try:
        from src.features.build_feature_matrix import features_for_variant

        data = features_for_variant(record)
    except Exception as exc:
        return None, f"Feature row could not be built: {exc}"
    missing = [column for column in catalog.feature_cols if column not in data]
    if missing:
        return None, "Feature row is missing columns: " + ", ".join(missing)
    return pd.DataFrame([{column: data[column] for column in catalog.feature_cols}]), None


def _imputation(catalog: Catalog, raw_score: float | None) -> tuple[bool, dict[str, Any]]:
    if raw_score is not None:
        return False, _legit_null()
    model_loaded = catalog.models.get("rf") is not None or catalog.models.get("xgb") is not None
    if not model_loaded:
        return False, absent(
            catalog.conservation_median_reason
            or "No fitted pipeline is loaded, so the training-fold median is not available."
        )
    if catalog.conservation_median is None:
        return True, absent(catalog.conservation_median_reason or "The fitted conservation median is unavailable.")
    return True, maybe(catalog.conservation_median)


def _three(letter: str) -> str | None:
    try:
        from src.preprocessing.clean_variants import AA_3TO1
    except Exception:
        return None
    inverse: dict[str, str] = {}
    for three, one in AA_3TO1.items():
        inverse.setdefault(one, three)
    return inverse.get(letter)


def _hgvs(record: dict[str, Any]) -> str | None:
    wt = _three(record["wt_aa"])
    mut = _three(record["mut_aa"])
    if wt is None or mut is None:
        return None
    return f"p.{wt}{record['position']}{mut}"


def _clinical(record: dict[str, Any]) -> dict[str, Any]:
    reasons: list[str] = []
    molecular = _parse_literal(record.get("molecular_consequence_list"), "molecular_consequence_list", reasons)
    variation = _parse_literal(record.get("variation_set"), "variation_set", reasons)
    cdna = None
    dbsnp = None
    location = None
    if isinstance(molecular, list):
        molecular_out = [str(item) for item in molecular]
    else:
        molecular_out = None
        if molecular is not None:
            reasons.append("molecular_consequence_list did not parse to a list.")
    if isinstance(variation, list) and variation and isinstance(variation[0], dict):
        measure = variation[0]
        cdna = measure.get("cdna_change")
        for xref in measure.get("variation_xrefs") or []:
            if isinstance(xref, dict) and xref.get("db_source") == "dbSNP" and xref.get("db_id") is not None:
                dbsnp = str(xref["db_id"])
                break
        location = _location(measure.get("variation_loc") or [])
    elif variation is not None:
        reasons.append("variation_set did not parse to a non-empty list of records.")
    return {
        "clinvar_uid": record.get("clinvar_uid"),
        "title": record.get("title"),
        "clinical_significance": record.get("clinical_significance"),
        "review_status": record.get("review_status"),
        "protein_change": record.get("protein_change"),
        "molecular_consequence": molecular_out,
        "cdna_change": cdna,
        "dbsnp": dbsnp,
        "location": location,
        "parse_reason": " ".join(reasons) or None,
    }


def _parse_literal(raw: str | None, field: str, reasons: list[str]) -> Any:
    if raw is None:
        reasons.append(f"{field} is empty.")
        return None
    if not isinstance(raw, str):
        return raw
    try:
        return ast.literal_eval(raw)
    except (SyntaxError, ValueError) as exc:
        reasons.append(f"{field} could not be parsed: {exc}")
        return None


def _location(locs: list[Any]) -> dict[str, Any] | None:
    chosen = None
    for loc in locs:
        if isinstance(loc, dict) and loc.get("assembly_name") == "GRCh38" and loc.get("status") == "current":
            chosen = loc
            break
    if chosen is None:
        for loc in locs:
            if isinstance(loc, dict) and loc.get("assembly_name") == "GRCh38":
                chosen = loc
                break
    if chosen is None:
        return None
    return {
        "assembly": chosen.get("assembly_name"),
        "chr": None if chosen.get("chr") in (None, "") else str(chosen.get("chr")),
        "start": None if chosen.get("start") in (None, "") else str(chosen.get("start")),
        "stop": None if chosen.get("stop") in (None, "") else str(chosen.get("stop")),
    }


def _sequence(catalog: Catalog, record: dict[str, Any]) -> dict[str, Any]:
    if catalog.uniprot is None:
        return absent(catalog.uniprot_reason or "data/raw/uniprot_entries.json is not present.")
    entry = catalog.uniprot.get(record["gene_symbol"])
    if entry is None:
        return absent(f"uniprot_entries.json has no entry for {record['gene_symbol']}.")
    track = None
    if catalog.conservation is not None:
        track = {
            str(pos): float(score)
            for pos, score in (catalog.conservation.get(record["gene_symbol"]) or {}).items()
        }
    domains = []
    for feature in entry.get("features") or []:
        domains.append(
            {
                "type": feature.get("type"),
                "description": feature.get("description"),
                "start": feature.get("start"),
                "end": feature.get("end"),
            }
        )
    return maybe(
        {
            "uniprot_id": entry.get("uniprot_id"),
            "sequence": entry.get("sequence") or "",
            "sequence_length": entry.get("sequence_length") or len(entry.get("sequence") or ""),
            "domains": domains,
            "conservation_by_position": track,
        }
    )


def _predict(catalog: Catalog, model_key: str, frame: Any, frame_reason: str | None) -> dict[str, Any]:
    model = catalog.models.get(model_key)
    if model is None:
        return _prediction(None, None, catalog.model_reasons.get(model_key) or "Model is not loaded.")
    if frame is None:
        return _prediction(None, None, frame_reason)
    try:
        probability = float(model.predict_proba(frame)[0, 1])
        predicted = int(model.predict(frame)[0])
    except Exception as exc:
        return _prediction(None, None, f"Prediction failed: {exc}")
    label = "pathogenic" if predicted == 1 else "benign"
    return _prediction(label, probability, None)


def _split_maybe(catalog: Catalog, record: dict[str, Any]) -> dict[str, Any]:
    if record["cohort"] == "vus" or record["label"] == "vus":
        return absent("VUS are held out of the gene-level split.")
    assigned = catalog.gene_split.get(record["gene_symbol"])
    if assigned:
        return maybe(assigned)
    return absent(catalog.split_error or f"No split assignment is available for {record['gene_symbol']}.")


def variant_detail(catalog: Catalog, record: dict[str, Any]) -> dict[str, Any]:
    inputs = build_inputs(catalog, record)
    gene = _gene_meta(catalog, record["gene_symbol"])
    predictions = {
        "rf": _predict(catalog, "rf", inputs.frame, inputs.frame_reason),
        "xgb": _predict(catalog, "xgb", inputs.frame, inputs.frame_reason),
    }
    return {
        "variant_id": record["variant_id"],
        "found": True,
        "reason": None,
        "gene_symbol": record["gene_symbol"],
        "panel": record.get("panel"),
        "condition": gene.get("condition"),
        "uniprot_id": gene.get("uniprot_id"),
        "label": record["label"],
        "clinical": inputs.clinical,
        "wt": {"one": record["wt_aa"], "three": inputs.wt_three},
        "mut": {"one": record["mut_aa"], "three": inputs.mut_three},
        "position": record["position"],
        "hgvs_protein": inputs.hgvs,
        "physchem": inputs.physchem,
        "conservation_score": inputs.conservation_score,
        "conservation_imputed": inputs.conservation_imputed,
        "conservation_imputed_value": inputs.conservation_imputed_value,
        "in_annotated_domain": inputs.in_annotated_domain,
        "domain_type": inputs.domain_type,
        "sequence": inputs.sequence,
        "predictions": predictions,
        "shap": {
            "rf": local_shap(catalog, "rf", record["variant_id"]),
            "xgb": local_shap(catalog, "xgb", record["variant_id"]),
        },
        "split": _split_maybe(catalog, record),
        "has_structure": _has_structure(catalog, record["variant_id"]),
    }


def unavailable_detail(variant_id: str, reason: str) -> dict[str, Any]:
    empty_pred = _prediction(None, None, reason)
    return {
        "variant_id": variant_id,
        "found": False,
        "reason": reason,
        "gene_symbol": None,
        "panel": None,
        "condition": None,
        "uniprot_id": None,
        "label": None,
        "clinical": None,
        "wt": None,
        "mut": None,
        "position": None,
        "hgvs_protein": None,
        "physchem": None,
        "conservation_score": absent(reason),
        "conservation_imputed": False,
        "conservation_imputed_value": absent(reason),
        "in_annotated_domain": absent(reason),
        "domain_type": absent(reason),
        "sequence": absent(reason),
        "predictions": {"rf": empty_pred, "xgb": empty_pred},
        "shap": {"rf": _empty_shap(reason), "xgb": _empty_shap(reason)},
        "split": absent(reason),
        "has_structure": False,
    }


def _gene_meta(catalog: Catalog, symbol: str) -> dict[str, Any]:
    for gene in catalog.genes:
        if gene.get("symbol") == symbol:
            return gene
    return {}


def local_shap(catalog: Catalog, model_key: str, variant_id: str) -> dict[str, Any]:
    return json.loads(cached_local_shap(model_key, variant_id))


@lru_cache(maxsize=4096)
def cached_local_shap(model_key: str, variant_id: str) -> str:
    return json.dumps(_local_shap_uncached(model_key, variant_id))


def _local_shap_uncached(model_key: str, variant_id: str) -> dict[str, Any]:
    catalog = get_catalog()
    resolved = resolve_variant(catalog, variant_id)
    if resolved.record is None:
        return _empty_shap(resolved.reason or "Variant is not available.")
    model = catalog.models.get(model_key)
    if model is None:
        return _empty_shap(catalog.model_reasons.get(model_key) or "Model is not loaded.")
    if catalog.shap_module is None:
        return _empty_shap(catalog.shap_import_error or "shap is not installed.")
    inputs = build_inputs(catalog, resolved.record)
    if inputs.frame is None:
        return _empty_shap(inputs.frame_reason or "The model frame could not be built.")
    try:
        import numpy as np

        values, names, base = _explain_matrix(catalog, model_key, model, inputs.frame)
        row = np.asarray(values[0], dtype=float)
        base_value = float(np.asarray(base).reshape(-1)[0])
        transformed = _dense_transform(model, inputs.frame)
    except Exception as exc:
        return _empty_shap(f"Local SHAP failed: {exc}")
    features = []
    for index, name in enumerate(names):
        raw_name = _raw_feature_name(name)
        raw_value = _py_scalar(inputs.frame.iloc[0][raw_name]) if raw_name in inputs.frame.columns else None
        if isinstance(raw_value, float):
            raw_json: Any = raw_value
        elif isinstance(raw_value, (bool, int, str)) or raw_value is None:
            raw_json = raw_value
        else:
            raw_json = str(raw_value)
        features.append(
            {
                "feature": name,
                "raw_feature": raw_name,
                "shap_value": float(row[index]),
                "feature_value_transformed": float(transformed[index]),
                "feature_value_raw": raw_json,
            }
        )
    features.sort(key=lambda item: (-abs(item["shap_value"]), item["feature"]))
    return {
        "shap_space": "probability",
        "base_value": base_value,
        "output_value": float(base_value + float(row.sum())),
        "features": features,
        "reason": None,
    }


def _dense_transform(model: Any, frame: Any) -> Any:
    import numpy as np

    transformed = model.named_steps["preprocess"].transform(frame)
    if hasattr(transformed, "toarray"):
        transformed = transformed.toarray()
    return np.asarray(transformed)[0]


def structure_payload(catalog: Catalog, record: dict[str, Any]) -> dict[str, Any]:
    known = _known_structures(catalog, record["gene_symbol"])
    ddg = absent(DDG_REASON)
    if record["variant_id"] != POC_VARIANT_ID:
        return {
            "variant_id": record["variant_id"],
            "has_structure": False,
            "reason": OTHER_STRUCTURE_REASON,
            "poc": None,
            "files": None,
            "rmsd_angstrom": absent("No mutant/WT PDB pair for this variant."),
            "ddg": ddg,
            "method_note": None,
            "esmfold_note": ESMFOLD_NOTE,
            "known_structures": known,
        }
    wt_present = catalog.pdb_wt_present
    mut_present = catalog.pdb_mut_present
    files = {
        "wt": {
            "name": WT_PDB_NAME,
            "present": wt_present,
            "url": f"/files/pdb/{WT_PDB_NAME}" if wt_present else None,
        },
        "mutant": {
            "name": MUTANT_PDB_NAME,
            "present": mut_present,
            "url": f"/files/pdb/{MUTANT_PDB_NAME}" if mut_present else None,
        },
    }
    missing = []
    if not wt_present:
        missing.append(f"outputs/{WT_PDB_NAME} is not present.")
    if not mut_present:
        missing.append(f"outputs/{MUTANT_PDB_NAME} is not present.")
    has_structure = not missing
    if has_structure:
        rmsd = _compute_rmsd(catalog)
        reason = None
    else:
        rmsd = absent(" ".join(missing) + " RMSD was not computed.")
        reason = " ".join(missing)
    return {
        "variant_id": record["variant_id"],
        "has_structure": has_structure,
        "reason": reason,
        "poc": {
            "defined": True,
            "pdb_id": "9B3R",
            "chain": "A",
            "residue": 39,
            "wt_aa": "R",
            "mut_aa": "H",
            "source_note": (
                "Hardcoded proof-of-concept in src/structure/mutate_pdb.py. "
                "Not read from structure_manifest.csv."
            ),
        },
        "files": files,
        "rmsd_angstrom": rmsd,
        "ddg": ddg,
        "method_note": METHOD_NOTE,
        "esmfold_note": ESMFOLD_NOTE,
        "known_structures": known,
    }


def unavailable_structure(variant_id: str, reason: str) -> dict[str, Any]:
    return {
        "variant_id": variant_id,
        "has_structure": False,
        "reason": reason,
        "poc": None,
        "files": None,
        "rmsd_angstrom": absent(reason),
        "ddg": absent(DDG_REASON),
        "method_note": None,
        "esmfold_note": ESMFOLD_NOTE,
        "known_structures": absent(reason),
    }


def _known_structures(catalog: Catalog, gene: str) -> dict[str, Any]:
    if catalog.manifest is None:
        return absent(catalog.manifest_reason or "data/raw/structure_manifest.csv is not present.")
    rows = []
    for row in catalog.manifest:
        if row.get("gene_symbol") != gene:
            continue
        rows.append(
            {
                "gene_symbol": row.get("gene_symbol"),
                "uniprot_id": row.get("uniprot_id"),
                "source": row.get("source"),
                "structure_id": row.get("structure_id"),
                "note": KNOWN_STRUCTURE_NOTE,
            }
        )
    return maybe(rows)


def _compute_rmsd(catalog: Catalog) -> dict[str, Any]:
    wt = catalog.outputs_dir / WT_PDB_NAME
    mutant = catalog.outputs_dir / MUTANT_PDB_NAME
    try:
        from src.structure.stability_analysis import compute_rmsd

        value = float(compute_rmsd(str(wt), str(mutant)))
    except Exception as exc:
        return absent(f"compute_rmsd failed: {exc}")
    if math.isnan(value):
        return absent("compute_rmsd returned NaN.")
    return maybe(value)


def report_payload(catalog: Catalog, record: dict[str, Any]) -> dict[str, Any]:
    detail = variant_detail(catalog, record)
    structure = structure_payload(catalog, record)
    inputs = build_inputs(catalog, record)
    rmsd_value = structure["rmsd_angstrom"]["value"]
    stability = METHOD_NOTE if record["variant_id"] == POC_VARIANT_ID else OTHER_STABILITY
    summaries = {}
    top_shap = {}
    for model_key, _reported_key, _title in MODEL_KEYS:
        shap_body = detail["shap"][model_key]
        prediction = detail["predictions"][model_key]
        top = list(shap_body["features"][:5])
        top_shap[model_key] = top
        summaries[model_key] = _summary_from_report(
            record,
            inputs,
            prediction,
            top,
            rmsd_value if isinstance(rmsd_value, (int, float)) else None,
            stability,
        )
    rf_label = detail["predictions"]["rf"]["label"]
    xgb_label = detail["predictions"]["xgb"]["label"]
    if rf_label and xgb_label:
        agreement = "agree" if rf_label == xgb_label else "disagree"
    else:
        agreement = None
    missing = [{"field": "predicted_ddg", "reason": DDG_REASON}]
    if detail["conservation_score"]["value"] is None:
        missing.append(
            {"field": "conservation_score", "reason": detail["conservation_score"]["reason"]}
        )
    if structure["rmsd_angstrom"]["value"] is None:
        missing.append({"field": "structural_rmsd", "reason": structure["rmsd_angstrom"]["reason"]})
    for model_key, _reported_key, _title in MODEL_KEYS:
        prediction = detail["predictions"][model_key]
        if prediction["prob_pathogenic"] is None:
            missing.append({"field": f"predictions.{model_key}", "reason": prediction["reason"]})
        if summaries[model_key]["value"] is None:
            missing.append({"field": f"summaries.{model_key}", "reason": summaries[model_key]["reason"]})
    body = {
        "variant_id": record["variant_id"],
        "disclaimer": DISCLAIMER,
        "clinical": detail["clinical"],
        "physchem": detail["physchem"],
        "conservation": {
            "score": detail["conservation_score"],
            "imputed": detail["conservation_imputed"],
            "method": CONSERVATION_METHOD,
        },
        "predictions": detail["predictions"],
        "agreement": agreement,
        "top_shap": top_shap,
        "structure": structure,
        "summaries": summaries,
        "interpretation": _interpretation(catalog, record, detail, structure),
        "missing": missing,
    }
    return body


def _foldx_ddg(inputs: Inputs) -> float | None:
    """FoldX ddG from the feature frame, or None when the variant was not modelled."""
    frame = inputs.frame
    if frame is None or "foldx_ddg" not in getattr(frame, "columns", []):
        return None
    value = frame["foldx_ddg"].iloc[0]
    return None if value is None or math.isnan(float(value)) else float(value)


def _summary_from_report(
    record: dict[str, Any],
    inputs: Inputs,
    prediction: dict[str, Any],
    top: list[dict[str, Any]],
    rmsd: float | None,
    stability: str,
) -> dict[str, Any]:
    if prediction.get("prob_pathogenic") is None:
        return absent(prediction.get("reason") or "The model was not run, so make_report was not called.")
    if not top:
        return absent("SHAP was not computed, so make_report was not called.")
    try:
        from src.report.generate_report import make_report

        label = "Pathogenic" if prediction["label"] == "pathogenic" else "Benign"
        made = make_report(
            gene_symbol=record["gene_symbol"],
            variant=inputs.hgvs or record["variant_id"],
            predicted_label=label,
            confidence_score=float(prediction["prob_pathogenic"]),
            key_features=[
                {"feature": item["feature"], "shap_value": float(item["shap_value"])} for item in top
            ],
            structural_rmsd=rmsd,
            predicted_ddg=_foldx_ddg(inputs),
            stability_interpretation=stability,
        )
    except Exception as exc:
        return absent(f"make_report failed: {exc}")
    return maybe(made.plain_language_summary)


def _interpretation(
    catalog: Catalog,
    record: dict[str, Any],
    detail: dict[str, Any],
    structure: dict[str, Any],
) -> str:
    clinical = detail.get("clinical") or {}
    hgvs = detail.get("hgvs_protein") or record["variant_id"]
    significance = clinical.get("clinical_significance") or "unspecified significance"
    review = clinical.get("review_status") or "unspecified review status"
    paragraphs = [
        (
            f"ClinVar records {record['gene_symbol']} {hgvs} as {significance} "
            f"(normalized label: {record['label']}). Review status: {review}. "
            "This clinical label is separate from the model predictions."
        )
    ]
    if record["label"] == "vus" or record["cohort"] == "vus":
        paragraphs.append(
            "This row is a variant of uncertain significance. "
            "It was held out of training and out of the gene-level split."
        )
    else:
        assigned = catalog.gene_split.get(record["gene_symbol"])
        if assigned == "test_gene":
            paragraphs.append(
                f"{record['gene_symbol']} is in the test fold of the gene-level split "
                "(GroupShuffleSplit, test_size=0.25, random_state=42). "
                "No variant of this gene was available to the model during training."
            )
        elif assigned == "train_gene":
            paragraphs.append(
                f"{record['gene_symbol']} is in the training fold of the gene-level split "
                "(GroupShuffleSplit, test_size=0.25, random_state=42). "
                "Every labelled variant of this gene was eligible for training, "
                "so this prediction is not an unseen-gene result."
            )
        else:
            paragraphs.append(catalog.split_error or "The gene-level split could not be computed for this gene.")
        if catalog.split is not None and catalog.split.caveat:
            paragraphs.append(catalog.split.caveat)
    paragraphs.append(
        "The trained model uses panel as a feature. Panel is constant for every variant of a gene "
        f"({record.get('panel') or 'unknown'} for {record['gene_symbol']}), so it carries the "
        "disease-panel base rate rather than a property of this substitution. The model also uses "
        "position, an absolute residue index that is not comparable across proteins."
    )
    score = detail["conservation_score"]
    if score["value"] is not None:
        paragraphs.append(
            f"The homolog-identity conservation score at this position is {score['value']}. {CONSERVATION_METHOD}"
        )
    elif detail["conservation_imputed"] and detail["conservation_imputed_value"]["value"] is not None:
        paragraphs.append(
            "No homolog-identity score was available for this position. "
            f"The loaded pipeline would median-impute it with {detail['conservation_imputed_value']['value']}. "
            f"{CONSERVATION_METHOD}"
        )
    else:
        paragraphs.append(f"The conservation score was not computed. {score['reason']} {CONSERVATION_METHOD}")
    for model_key, title in (("rf", "Random Forest"), ("xgb", "XGBoost")):
        prediction = detail["predictions"][model_key]
        if prediction["prob_pathogenic"] is None:
            paragraphs.append(f"{title} was not run. {prediction['reason']}")
        else:
            paragraphs.append(
                f"{title} predicts {prediction['label']} with "
                f"P(pathogenic)={prediction['prob_pathogenic']:.3f} at the default 0.5 threshold. "
                "This probability is not calibrated."
            )
    if record["variant_id"] == POC_VARIANT_ID:
        rmsd = structure["rmsd_angstrom"]
        if rmsd["value"] is None:
            paragraphs.append(f"{METHOD_NOTE} RMSD was not recomputed. {rmsd['reason']} {DDG_REASON}")
        else:
            paragraphs.append(
                f"Recomputed Cα RMSD is {rmsd['value']} Å. {METHOD_NOTE} {DDG_REASON}"
            )
    else:
        paragraphs.append(f"{OTHER_STABILITY} {DDG_REASON}")
    paragraphs.append(DISCLAIMER)
    return "\n\n".join(paragraphs)


def unavailable_report(variant_id: str, reason: str) -> dict[str, Any]:
    empty_pred = _prediction(None, None, reason)
    return {
        "variant_id": variant_id,
        "disclaimer": DISCLAIMER,
        "clinical": None,
        "physchem": None,
        "conservation": {"score": absent(reason), "imputed": False, "method": CONSERVATION_METHOD},
        "predictions": {"rf": empty_pred, "xgb": empty_pred},
        "agreement": None,
        "top_shap": {"rf": [], "xgb": []},
        "structure": unavailable_structure(variant_id, reason),
        "summaries": {"rf": absent(reason), "xgb": absent(reason)},
        "interpretation": f"{reason}\n\n{DISCLAIMER}",
        "missing": [{"field": "variant", "reason": reason}, {"field": "predicted_ddg", "reason": DDG_REASON}],
    }


def render_markdown(body: dict[str, Any]) -> str:
    lines = [
        f"# {body['variant_id']}",
        "",
        f"**{body['disclaimer']}**",
        "",
        "## Interpretation",
        "",
        body.get("interpretation") or "",
        "",
        "## Clinical",
        "",
        _md_dump(body.get("clinical")),
        "",
        "## Physicochemical features",
        "",
        _md_dump(body.get("physchem")),
        "",
        "## Conservation",
        "",
        _md_dump(body.get("conservation")),
        "",
        "## Predictions",
        "",
        _md_dump(body.get("predictions")),
        "",
        f"Agreement: {body.get('agreement')}",
        "",
        "## Top SHAP contributions",
        "",
        _md_dump(body.get("top_shap")),
        "",
        "## Structure",
        "",
        _md_dump(body.get("structure")),
        "",
        "## Model summaries",
        "",
        _md_dump(body.get("summaries")),
        "",
        "## Missing",
        "",
        _md_dump(body.get("missing")),
        "",
    ]
    return "\n".join(lines)


def _md_dump(value: Any) -> str:
    return "```json\n" + json.dumps(value, indent=2) + "\n```"


def resolve_pdb(catalog: Catalog, name: str) -> Path | None:
    if not PDB_NAME.fullmatch(name):
        return None
    root = catalog.outputs_dir.resolve()
    path = (root / name).resolve()
    if path.parent != root:
        return None
    if not path.is_file():
        return None
    return path


def ablation_payload(catalog: Catalog) -> dict[str, Any]:
    path = catalog.outputs_dir / "ablation_results.json"
    if not path.is_file():
        return {"value": None, "reason": ABLATION_REASON, "plan": ABLATION_PLAN}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        return {
            "value": None,
            "reason": f"outputs/ablation_results.json could not be read: {exc}",
            "plan": ABLATION_PLAN,
        }
    return {"value": payload, "reason": None, "plan": ABLATION_PLAN}
