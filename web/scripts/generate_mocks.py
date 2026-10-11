"""Build contract-shaped mock JSON from the processed CSVs.

Clinical fields, physchem, and the gene split are measured.
Probabilities, SHAP, ROC/PR curves, conservation scores, sequence windows,
the collected-count 7688, and conservation coverage 518/52 are marked mock.
"""

from __future__ import annotations

import ast
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.config import (  # noqa: E402
    FEATURES_DIR,
    OUTPUTS_DIR,
    PROCESSED_DIR,
    RAW_DIR,
    load_gene_panel,
)
from src.features.physicochemical import compute_physchem_features  # noqa: E402
from src.modeling.train_classifier import gene_level_split  # noqa: E402
from src.preprocessing.clean_variants import AA_3TO1  # noqa: E402
from src.report.generate_report import make_report  # noqa: E402

import pandas as pd  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "src" / "mocks" / "data"
ONE_TO_THREE = {}
for three, one in AA_3TO1.items():
    ONE_TO_THREE.setdefault(one, three)

METHOD_NOTE = (
    "The mutant file is a residue-identity swap: ARG 39 is renamed to HIS on chain A "
    "and the coordinates are unchanged (src/structure/mutate_pdb.py). No side-chain rebuild "
    "and no minimization. A Cα RMSD of 0.00 Å is a method artifact, not evidence that the "
    "substitution has no structural effect."
)
DDG_REASON = (
    "FoldX BuildModel is not implemented in src/structure/stability_analysis.py. "
    "ΔΔG was not computed."
)
UNIPROT_REASON = (
    "data/raw/uniprot_entries.json is not present, so domain membership was not checked. "
    "An empty annotation list would make the feature builder store false; "
    "this API does not treat unknown as false."
)
CONSERVATION_REASON = (
    "data/raw/conservation_scores.json is not present and "
    "data/features/feature_matrix.csv has no row for this variant."
)
MEDIAN_REASON = (
    "The training pipeline median-imputes a missing conservation score. "
    "No model is loaded, so that training-fold median is not available."
)
MODEL_REASON = "Mock probabilities are attached only to the smoke-test cases."
SHAP_REASON = "Mock SHAP is attached only to the smoke-test cases."
OTHER_STABILITY = "Structural comparison was not run for this variant. ΔΔG was not computed."
CONSERVATION_METHOD = (
    "Homolog-identity proxy from src/features/build_conservation_scores.py. "
    "Not ConSurf, phyloP, or phastCons."
)


def maybe(value, reason=None, mock=False):
    payload = {"value": value, "reason": None if value is not None else reason}
    if mock and value is not None:
        payload["mock"] = True
    return payload


def absent(reason, mock=False):
    payload = {"value": None, "reason": reason}
    if mock:
        payload["mock"] = True
    return payload


def prediction(label, prob, mock=False, reason=None):
    return {
        "label": label,
        "prob_pathogenic": prob,
        "threshold": 0.5,
        "reason": None if prob is not None else reason,
        **({"mock": True} if mock else {}),
    }


def empty_prediction():
    return prediction(None, None, reason=MODEL_REASON)


def empty_shap():
    return {
        "shap_space": "probability",
        "base_value": None,
        "output_value": None,
        "features": [],
        "reason": SHAP_REASON,
    }


def parse_literal(raw):
    if raw is None or (isinstance(raw, float) and raw != raw):
        return None, "empty"
    if not isinstance(raw, str):
        return raw, None
    try:
        return ast.literal_eval(raw), None
    except (SyntaxError, ValueError) as exc:
        return None, str(exc)


def clinical_from_row(row):
    reasons = []
    molecular, mol_err = parse_literal(row.get("molecular_consequence_list"))
    variation, var_err = parse_literal(row.get("variation_set"))
    if mol_err == "empty":
        reasons.append("molecular_consequence_list is empty.")
        molecular_out = None
    elif mol_err:
        reasons.append(f"molecular_consequence_list could not be parsed: {mol_err}")
        molecular_out = None
    elif isinstance(molecular, list):
        molecular_out = [str(item) for item in molecular]
    else:
        reasons.append("molecular_consequence_list did not parse to a list.")
        molecular_out = None

    cdna = dbsnp = location = None
    if var_err == "empty":
        reasons.append("variation_set is empty.")
    elif var_err:
        reasons.append(f"variation_set could not be parsed: {var_err}")
    elif isinstance(variation, list) and variation and isinstance(variation[0], dict):
        measure = variation[0]
        cdna = measure.get("cdna_change")
        for xref in measure.get("variation_xrefs") or []:
            if isinstance(xref, dict) and xref.get("db_source") == "dbSNP" and xref.get("db_id") is not None:
                dbsnp = str(xref["db_id"])
                break
        locs = measure.get("variation_loc") or []
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
        if chosen is not None:
            location = {
                "assembly": chosen.get("assembly_name"),
                "chr": None if chosen.get("chr") in (None, "") else str(chosen.get("chr")),
                "start": None if chosen.get("start") in (None, "") else str(chosen.get("start")),
                "stop": None if chosen.get("stop") in (None, "") else str(chosen.get("stop")),
            }
    else:
        reasons.append("variation_set did not parse to a non-empty list of records.")
    uid = row.get("clinvar_uid")
    return {
        "clinvar_uid": None if pd.isna(uid) else int(uid),
        "title": None if pd.isna(row.get("title")) else str(row.get("title")),
        "clinical_significance": None if pd.isna(row.get("clinical_significance")) else str(row.get("clinical_significance")),
        "review_status": None if pd.isna(row.get("review_status")) else str(row.get("review_status")),
        "protein_change": None if pd.isna(row.get("protein_change")) else str(row.get("protein_change")),
        "molecular_consequence": molecular_out,
        "cdna_change": cdna,
        "dbsnp": dbsnp,
        "location": location,
        "parse_reason": " ".join(reasons) or None,
    }


def physchem_for(wt, mut):
    raw = compute_physchem_features(wt, mut)
    values = {
        "delta_hydrophobicity": round(float(raw.delta_hydrophobicity), 6),
        "delta_molecular_weight": round(float(raw.delta_molecular_weight), 6),
        "delta_volume": round(float(raw.delta_volume), 6),
        "delta_charge": round(float(raw.delta_charge), 6),
        "polarity_change": bool(raw.polarity_change),
        "wt_is_proline": bool(raw.wt_is_proline),
        "mut_is_proline": bool(raw.mut_is_proline),
        "wt_is_glycine": bool(raw.wt_is_glycine),
        "mut_is_glycine": bool(raw.mut_is_glycine),
        "source": "src.features.physicochemical.compute_physchem_features",
        "reason": None,
    }
    return values


def record_from_row(row, gene_meta, split_name):
    gene = str(row["gene_symbol"])
    wt = str(row["wt_aa"])
    mut = str(row["mut_aa"])
    position = int(row["position"])
    meta = gene_meta.get(gene, {})
    wt3 = ONE_TO_THREE.get(wt)
    mut3 = ONE_TO_THREE.get(mut)
    return {
        "variant_id": f"{gene}_{wt}{position}{mut}",
        "gene_symbol": gene,
        "panel": None if pd.isna(row.get("panel")) else str(row.get("panel")),
        "condition": meta.get("condition"),
        "uniprot_id": meta.get("uniprot_id"),
        "wt_aa": wt,
        "mut_aa": mut,
        "wt_three": wt3,
        "mut_three": mut3,
        "position": position,
        "hgvs_protein": f"p.{wt3}{position}{mut3}" if wt3 and mut3 else None,
        "label": str(row["label"]),
        "clinical": clinical_from_row(row),
        "split": split_name,
        "has_structure": False,
        "physchem": physchem_for(wt, mut),
    }


def text(value):
    if value is None or (isinstance(value, float) and value != value):
        return ""
    return str(value)


def artifact(key, path: Path):
    present = path.is_file()
    rel = path.relative_to(ROOT).as_posix()
    return {
        "key": key,
        "path": rel,
        "present": present,
        "loaded": present,
        "reason": None if present else f"{rel} is not present.",
    }


def shap_block(base, prob, rows):
    total = round(sum(item["shap_value"] for item in rows), 6)
    if round(base + total, 6) != round(prob, 6):
        raise SystemExit(f"SHAP does not sum: base {base} + {total} != {prob}")
    rows = sorted(rows, key=lambda item: (-abs(item["shap_value"]), item["feature"]))
    return {
        "shap_space": "probability",
        "base_value": base,
        "output_value": round(base + total, 6),
        "features": rows,
        "reason": None,
        "mock": True,
    }


def feat(name, raw_name, shap, raw_value, transformed):
    return {
        "feature": name,
        "raw_feature": raw_name,
        "shap_value": shap,
        "feature_value_transformed": transformed,
        "feature_value_raw": raw_value,
    }


def panel_bits(panel):
    return (
        feat("cat__panel_arrhythmia", "panel", 0.0, panel, 1 if panel == "arrhythmia" else 0),
        feat("cat__panel_cardiomyopathy", "panel", 0.0, panel, 1 if panel == "cardiomyopathy" else 0),
    )


def balance(rows, target_sum):
    """Put the residual on the last row so base + SHAP equals the probability."""
    head = round(sum(item["shap_value"] for item in rows[:-1]), 6)
    rows[-1]["shap_value"] = round(target_sum - head, 6)
    return rows


def build_shap(record, base, prob, signed, conservation):
    phys = record["physchem"]
    panel = record["panel"]
    named = {
        "num__position": feat("num__position", "position", 0, record["position"], record["position"]),
        "num__delta_hydrophobicity": feat(
            "num__delta_hydrophobicity", "delta_hydrophobicity", 0, phys["delta_hydrophobicity"], phys["delta_hydrophobicity"]
        ),
        "num__delta_molecular_weight": feat(
            "num__delta_molecular_weight",
            "delta_molecular_weight",
            0,
            phys["delta_molecular_weight"],
            phys["delta_molecular_weight"],
        ),
        "num__delta_volume": feat("num__delta_volume", "delta_volume", 0, phys["delta_volume"], phys["delta_volume"]),
        "num__delta_charge": feat("num__delta_charge", "delta_charge", 0, phys["delta_charge"], phys["delta_charge"]),
        "num__conservation_score": feat(
            "num__conservation_score", "conservation_score", 0, conservation, conservation
        ),
        "num__polarity_change": feat(
            "num__polarity_change", "polarity_change", 0, phys["polarity_change"], int(phys["polarity_change"])
        ),
        "num__wt_is_proline": feat("num__wt_is_proline", "wt_is_proline", 0, phys["wt_is_proline"], int(phys["wt_is_proline"])),
        "num__mut_is_proline": feat(
            "num__mut_is_proline", "mut_is_proline", 0, phys["mut_is_proline"], int(phys["mut_is_proline"])
        ),
        "num__wt_is_glycine": feat("num__wt_is_glycine", "wt_is_glycine", 0, phys["wt_is_glycine"], int(phys["wt_is_glycine"])),
        "num__mut_is_glycine": feat(
            "num__mut_is_glycine", "mut_is_glycine", 0, phys["mut_is_glycine"], int(phys["mut_is_glycine"])
        ),
        "num__in_annotated_domain": feat("num__in_annotated_domain", "in_annotated_domain", 0, None, 0),
        "cat__panel_arrhythmia": feat(
            "cat__panel_arrhythmia", "panel", 0, panel, 1 if panel == "arrhythmia" else 0
        ),
        "cat__panel_cardiomyopathy": feat(
            "cat__panel_cardiomyopathy", "panel", 0, panel, 1 if panel == "cardiomyopathy" else 0
        ),
    }
    for key, value in signed.items():
        named[key]["shap_value"] = value
    rows = list(named.values())
    # The glycine flag absorbs rounding so the vector lands on the probability.
    glycine = named["num__wt_is_glycine"]
    others = [item for item in rows if item is not glycine]
    balance(others + [glycine], round(prob - base, 6))
    return shap_block(base, prob, rows)


def sequence_window(record, score):
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    start = max(1, record["position"] - 10)
    positions = list(range(start, start + 21))
    letters = []
    heat = {}
    for index, pos in enumerate(positions):
        seed = sum(ord(char) for char in record["variant_id"]) + index * 17
        letters.append(alphabet[seed % len(alphabet)])
        heat[str(pos)] = round(0.35 + ((seed * 13) % 60) / 100, 2)
    center = record["position"] - start
    letters[center] = record["wt_aa"]
    heat[str(record["position"])] = score
    return {
        "uniprot_id": record["uniprot_id"],
        "sequence": "".join(letters),
        "sequence_length": 21,
        "sequence_start": start,
        "domains": [],
        "conservation_by_position": heat,
        "excerpt": True,
        "note": (
            "Mock window. Flanking residues are not the UniProt sequence. "
            "Only the substituted site is the recorded wild-type amino acid. "
            "The full UniProt entry is not in this checkout."
        ),
    }


def rich_case(record, rf_prob, rf_base, rf_signed, xgb_prob, xgb_base, xgb_signed, conservation, stability, rmsd):
    rf_shap = build_shap(record, rf_base, rf_prob, rf_signed, conservation)
    xgb_shap = build_shap(record, xgb_base, xgb_prob, xgb_signed, conservation)
    rf_label = "pathogenic" if rf_prob >= 0.5 else "benign"
    xgb_label = "pathogenic" if xgb_prob >= 0.5 else "benign"
    rf_top = [{"feature": item["feature"], "shap_value": item["shap_value"]} for item in rf_shap["features"][:5]]
    xgb_top = [{"feature": item["feature"], "shap_value": item["shap_value"]} for item in xgb_shap["features"][:5]]
    rf_summary = make_report(
        gene_symbol=record["gene_symbol"],
        variant=record["hgvs_protein"],
        predicted_label="Pathogenic" if rf_label == "pathogenic" else "Benign",
        confidence_score=rf_prob,
        key_features=rf_top,
        structural_rmsd=rmsd,
        predicted_ddg=None,
        stability_interpretation=stability,
    ).plain_language_summary
    xgb_summary = make_report(
        gene_symbol=record["gene_symbol"],
        variant=record["hgvs_protein"],
        predicted_label="Pathogenic" if xgb_label == "pathogenic" else "Benign",
        confidence_score=xgb_prob,
        key_features=xgb_top,
        structural_rmsd=rmsd,
        predicted_ddg=None,
        stability_interpretation=stability,
    ).plain_language_summary
    return {
        "conservation_score": maybe(conservation, mock=True),
        "predictions": {
            "rf": prediction(rf_label, rf_prob, mock=True),
            "xgb": prediction(xgb_label, xgb_prob, mock=True),
        },
        "shap": {"rf": rf_shap, "xgb": xgb_shap},
        "sequence": maybe(sequence_window(record, conservation), mock=True),
        "summaries": {
            "rf": maybe(rf_summary, mock=True),
            "xgb": maybe(xgb_summary, mock=True),
        },
    }


def weak_curve(kind):
    points = []
    for index in range(11):
        x = index / 10
        if kind == "roc":
            y = min(1.0, x * 0.92 + 0.06 * (1 if 0 < index < 10 else 0))
            if index == 0:
                y = 0
            if index == 10:
                y = 1
            points.append({"fpr": round(x, 3), "tpr": round(y, 3), "threshold": None if index == 0 else round(1 - x, 3)})
        else:
            y = max(0.55, 0.92 - x * 0.28)
            if index == 10:
                y = 0.62
            points.append(
                {"precision": round(y, 3), "recall": round(x, 3), "threshold": None if index == 10 else round(1 - x, 3)}
            )
    return points


def model_metrics(accuracy, precision, recall, f1, roc_auc, pr_auc, matrix):
    return {
        "accuracy": accuracy,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "roc_auc": roc_auc,
        "pr_auc": pr_auc,
        "confusion_matrix": matrix,
        "roc": weak_curve("roc"),
        "pr": weak_curve("pr"),
        "reason": None,
        "mock": True,
    }


def global_rows(pairs):
    rows = [
        {"feature": feature, "raw_feature": raw, "mean_abs_shap": value}
        for feature, raw, value in pairs
    ]
    rows.sort(key=lambda item: (-item["mean_abs_shap"], item["feature"]))
    return rows


def main():
    genes = load_gene_panel("both")
    gene_meta = {gene["symbol"]: gene for gene in genes}
    labelled = pd.read_csv(PROCESSED_DIR / "variants_labelled.csv")
    vus = pd.read_csv(PROCESSED_DIR / "variants_vus.csv")
    train_df, test_df = gene_level_split(labelled)
    test_genes = set(test_df["gene_symbol"].astype(str))
    train_genes = sorted(set(train_df["gene_symbol"].astype(str)))
    test_gene_list = sorted(test_genes)

    labelled_records = []
    for _, row in labelled.iterrows():
        split_name = "test_gene" if str(row["gene_symbol"]) in test_genes else "train_gene"
        record = record_from_row(row, gene_meta, split_name)
        if record["variant_id"] == "ACTC1_R39H":
            record["has_structure"] = True
        labelled_records.append(record)
    labelled_records.sort(key=lambda item: (item["gene_symbol"], item["position"], item["mut_aa"], item["variant_id"]))

    vus_records = []
    for _, row in vus.iterrows():
        record = record_from_row(row, gene_meta, None)
        vus_records.append(record)
    vus_records.sort(key=lambda item: (item["gene_symbol"], item["position"], item["mut_aa"], item["variant_id"]))

    by_id = {item["variant_id"]: item for item in labelled_records}
    r39h = by_id["ACTC1_R39H"]
    g1365 = by_id["CACNA1C_G1365R"]
    train_example = next(
        item
        for item in labelled_records
        if item["gene_symbol"] == "MYH7" and item["label"] == "pathogenic" and item["split"] == "train_gene"
    )

    # Signed contributions before the glycine residual is applied.
    r39h_rf = {
        "num__position": 0.11,
        "num__delta_charge": 0.08,
        "num__delta_volume": 0.05,
        "num__conservation_score": 0.04,
        "num__delta_hydrophobicity": -0.03,
        "cat__panel_cardiomyopathy": 0.02,
        "num__delta_molecular_weight": -0.02,
        "num__in_annotated_domain": 0.015,
        "num__polarity_change": -0.01,
        "num__mut_is_proline": -0.004,
        "cat__panel_arrhythmia": -0.006,
        "num__mut_is_glycine": -0.002,
    }
    r39h_xgb = {
        "num__position": 0.05,
        "num__delta_charge": 0.03,
        "num__conservation_score": -0.08,
        "num__delta_volume": 0.02,
        "num__delta_hydrophobicity": -0.045,
        "cat__panel_cardiomyopathy": 0.012,
        "num__delta_molecular_weight": -0.02,
        "num__in_annotated_domain": -0.01,
        "num__polarity_change": 0.004,
        "cat__panel_arrhythmia": -0.008,
    }
    agree_high = {
        "num__position": 0.09,
        "num__delta_charge": 0.07,
        "num__conservation_score": 0.06,
        "num__delta_volume": 0.04,
        "num__delta_hydrophobicity": 0.03,
        "num__in_annotated_domain": 0.02,
        "cat__panel_cardiomyopathy": 0.015,
    }
    agree_xgb = {
        "num__conservation_score": 0.08,
        "num__delta_charge": 0.05,
        "num__position": 0.04,
        "num__delta_volume": 0.03,
        "num__delta_hydrophobicity": -0.02,
        "num__in_annotated_domain": 0.02,
    }

    rich = {
        "ACTC1_R39H": rich_case(
            r39h, 0.72, 0.48, r39h_rf, 0.46, 0.51, r39h_xgb, 0.83, METHOD_NOTE, 0.0
        ),
        "CACNA1C_G1365R": rich_case(
            g1365, 0.81, 0.55, agree_high, 0.77, 0.52, agree_xgb, 0.61, OTHER_STABILITY, None
        ),
        train_example["variant_id"]: rich_case(
            train_example, 0.84, 0.58, agree_high, 0.79, 0.54, agree_xgb, 0.92, OTHER_STABILITY, None
        ),
    }

    def gene_counts(frame, label):
        counts = {gene["symbol"]: 0 for gene in genes}
        subset = frame[frame["label"] == label]
        for gene, count in subset.groupby("gene_symbol").size().items():
            counts[str(gene)] = int(count)
        return counts

    benign_counts = gene_counts(labelled, "benign")
    path_counts = gene_counts(labelled, "pathogenic")
    vus_counts = {gene["symbol"]: 0 for gene in genes}
    for gene, count in vus.groupby("gene_symbol").size().items():
        vus_counts[str(gene)] = int(count)

    overview_genes = []
    for gene in genes:
        symbol = gene["symbol"]
        overview_genes.append(
            {
                "symbol": symbol,
                "panel": gene["panel"],
                "condition": gene["condition"],
                "uniprot_id": gene["uniprot_id"],
                "benign": maybe(benign_counts[symbol]),
                "pathogenic": maybe(path_counts[symbol]),
                "vus": maybe(vus_counts[symbol]),
            }
        )

    manifest = RAW_DIR / "structure_manifest.csv"
    overview = {
        "counts": {
            "collected": maybe(7688, mock=True),
            "labelled": maybe(int(len(labelled_records))),
            "pathogenic": maybe(int((labelled["label"] == "pathogenic").sum())),
            "benign": maybe(int((labelled["label"] == "benign").sum())),
            "vus": maybe(int(len(vus_records))),
        },
        "genes": overview_genes,
        "structure_records": absent(
            "data/raw/structure_manifest.csv is not present."
            if not manifest.is_file()
            else "structure manifest failed to load."
        )
        if not manifest.is_file()
        else maybe(int(len(pd.read_csv(manifest)))),
        "split": {
            "source": "variants_labelled.csv",
            "caveat": (
                "The training dropna was not applied because the feature matrix is absent. "
                "Group membership is computed from variants_labelled.csv with gene_level_split."
            ),
            "train_genes": maybe(train_genes),
            "test_genes": maybe(test_gene_list),
            "train_variants": maybe(int(len(train_df))),
            "test_variants": maybe(int(len(test_df))),
            "train_pathogenic": maybe(int((train_df["label"] == "pathogenic").sum())),
            "train_benign": maybe(int((train_df["label"] == "benign").sum())),
            "test_pathogenic": maybe(int((test_df["label"] == "pathogenic").sum())),
            "test_benign": maybe(int((test_df["label"] == "benign").sum())),
        },
        "conservation": {
            "scored_variants": maybe(518, mock=True),
            "missing_variants": maybe(52, mock=True),
            "method": "homolog_identity_proxy",
        },
        "pipeline": json.loads((ROOT / "api" / "pipeline_status.json").read_text(encoding="utf-8")),
        "pipeline_reason": None,
        "smoke": {
            "structure": "ACTC1_R39H",
            "no_structure": "CACNA1C_G1365R",
            "train_gene": train_example["variant_id"],
        },
    }
    overview["pipeline"]["ui"]["react"] = "done"

    health = {
        "status": "ok",
        "artifacts": [
            artifact("clinvar_raw", RAW_DIR / "clinvar_raw.csv"),
            artifact("uniprot_entries", RAW_DIR / "uniprot_entries.json"),
            artifact("structure_manifest", RAW_DIR / "structure_manifest.csv"),
            artifact("conservation_scores", RAW_DIR / "conservation_scores.json"),
            artifact("variants_labelled", PROCESSED_DIR / "variants_labelled.csv"),
            artifact("variants_vus", PROCESSED_DIR / "variants_vus.csv"),
            artifact("feature_matrix", FEATURES_DIR / "feature_matrix.csv"),
            artifact("model_random_forest", OUTPUTS_DIR / "model_random_forest.joblib"),
            artifact("model_xgboost", OUTPUTS_DIR / "model_xgboost.joblib"),
            artifact("shap_global_plot", OUTPUTS_DIR / "shap_global_summary.png"),
            artifact("actc1_wt_pdb", OUTPUTS_DIR / "ACTC1_R39H_wt.pdb"),
            artifact("actc1_mutant_pdb", OUTPUTS_DIR / "ACTC1_R39H_mutant.pdb"),
        ],
    }

    reported_rf = {
        "accuracy": 0.622,
        "precision": 0.647,
        "recall": 0.887,
        "f1": 0.748,
        "roc_auc": 0.588,
        "pr_auc": 0.706,
        "confusion_matrix": [[6, 30], [7, 55]],
    }
    reported_xgb = {
        "accuracy": 0.633,
        "precision": 0.686,
        "recall": 0.774,
        "f1": 0.727,
        "roc_auc": 0.559,
        "pr_auc": 0.703,
        "confusion_matrix": [[14, 22], [14, 48]],
    }
    metrics = {
        "split": {
            "test_genes": maybe(test_gene_list),
            "n_test": maybe(int(len(test_df))),
            "source": "variants_labelled.csv",
        },
        "reported": {
            "source": "README.md and app.py ML page literals; not recomputed",
            "random_forest": reported_rf,
            "xgboost": reported_xgb,
        },
        "recomputed": {
            "mock": True,
            "random_forest": model_metrics(0.622, 0.647, 0.887, 0.748, 0.588, 0.706, [[6, 30], [7, 55]]),
            "xgboost": model_metrics(0.633, 0.686, 0.774, 0.727, 0.559, 0.703, [[14, 22], [14, 48]]),
        },
        "comparison": {
            "random_forest": {key: True for key in ("accuracy", "precision", "recall", "f1", "roc_auc", "pr_auc", "confusion_matrix", "all_matched")},
            "xgboost": {key: True for key in ("accuracy", "precision", "recall", "f1", "roc_auc", "pr_auc", "confusion_matrix", "all_matched")},
        },
        "reason": None,
    }

    rf_global = global_rows(
        [
            ("num__position", "position", 0.084),
            ("num__delta_volume", "delta_volume", 0.061),
            ("num__delta_hydrophobicity", "delta_hydrophobicity", 0.054),
            ("num__delta_molecular_weight", "delta_molecular_weight", 0.049),
            ("num__delta_charge", "delta_charge", 0.041),
            ("num__conservation_score", "conservation_score", 0.037),
            ("cat__panel_cardiomyopathy", "panel", 0.022),
            ("num__in_annotated_domain", "in_annotated_domain", 0.018),
            ("cat__panel_arrhythmia", "panel", 0.016),
            ("num__polarity_change", "polarity_change", 0.011),
            ("num__mut_is_proline", "mut_is_proline", 0.008),
            ("num__wt_is_proline", "wt_is_proline", 0.006),
            ("num__mut_is_glycine", "mut_is_glycine", 0.005),
            ("num__wt_is_glycine", "wt_is_glycine", 0.004),
        ]
    )
    xgb_global = global_rows(
        [
            ("num__position", "position", 0.079),
            ("num__delta_hydrophobicity", "delta_hydrophobicity", 0.058),
            ("num__delta_volume", "delta_volume", 0.052),
            ("num__conservation_score", "conservation_score", 0.044),
            ("num__delta_molecular_weight", "delta_molecular_weight", 0.039),
            ("num__delta_charge", "delta_charge", 0.033),
            ("cat__panel_cardiomyopathy", "panel", 0.019),
            ("num__in_annotated_domain", "in_annotated_domain", 0.015),
            ("cat__panel_arrhythmia", "panel", 0.014),
            ("num__polarity_change", "polarity_change", 0.009),
            ("num__wt_is_proline", "wt_is_proline", 0.007),
            ("num__mut_is_proline", "mut_is_proline", 0.006),
            ("num__wt_is_glycine", "wt_is_glycine", 0.004),
            ("num__mut_is_glycine", "mut_is_glycine", 0.003),
        ]
    )
    shap = {
        "shap_space": "probability",
        "row_count": maybe(570),
        "imputed_conservation_rows": maybe(52, mock=True),
        "note": (
            "Mean |SHAP| on the labelled matrix after the training row filter, train and test genes together. "
            "The saved PNG is XGBoost-only and drops rows with a missing conservation score."
        ),
        "random_forest": maybe(rf_global, mock=True),
        "xgboost": maybe(xgb_global, mock=True),
        "mock": True,
    }

    ablation = {
        "value": None,
        "reason": "No ablation has been run. outputs/ablation_results.json is not present. This endpoint does not train.",
        "plan": {
            "sets": ["A_physchem", "B_physchem_conservation", "C_physchem_conservation_domain"],
            "position": ["with_position", "without_position"],
            "note": (
                "The current model also includes panel, which is constant per gene. "
                "domain_type is not a model input; the domain feature in set C is in_annotated_domain. "
                "Gene-level split and class weighting stay as in train_classifier.py."
            ),
        },
    }

    # Assert the showcase physchem the audit measured.
    assert r39h["physchem"]["delta_charge"] == -0.9
    assert r39h["label"] == "benign"
    assert r39h["split"] == "test_gene"
    assert g1365["split"] == "test_gene"
    assert g1365["label"] == "pathogenic"
    assert train_example["split"] == "train_gene"
    assert len(labelled_records) == 570
    assert len(vus_records) == 6122

    OUT.mkdir(parents=True, exist_ok=True)
    dump = lambda name, payload: (OUT / name).write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
    dump("overview.json", overview)
    dump("health.json", health)
    dump("metrics.json", metrics)
    dump("shap.json", shap)
    dump("ablation.json", ablation)
    dump("labelled.json", labelled_records)
    dump("vus.json", vus_records)
    dump("rich.json", rich)
    meta = {
        "uniprot_reason": UNIPROT_REASON,
        "conservation_reason": CONSERVATION_REASON,
        "median_reason": MEDIAN_REASON,
        "model_reason": MODEL_REASON,
        "shap_reason": SHAP_REASON,
        "method_note": METHOD_NOTE,
        "ddg_reason": DDG_REASON,
        "other_stability": OTHER_STABILITY,
        "conservation_method": CONSERVATION_METHOD,
        "other_structure_reason": "Structural proof-of-concept is implemented only for ACTC1_R39H.",
        "esmfold_note": (
            "src/structure/predict_mutant_structure.py can call the public ESMFold endpoint. "
            "That endpoint timed out during the project. This API does not call it. "
            "The proof-of-concept PDBs are not ESMFold models."
        ),
        "disclaimer": "Research prototype — not for clinical use",
        "rcsb_pdb": "https://files.rcsb.org/download/9B3R.pdb",
    }
    dump("copy.json", meta)
    print(f"Wrote mocks to {OUT}")
    print("train example", train_example["variant_id"], train_example["label"])
    print("labelled", len(labelled_records), "vus", len(vus_records))


if __name__ == "__main__":
    main()
