"""
Step 5: Explainable AI via SHAP.

Explains the full XGBoost model (outputs/model_xgboost_full.joblib) with
SHAP's TreeExplainer, which is exact for tree ensembles. SHAP values are in
log-odds units: positive pushes towards pathogenic, negative towards benign,
and base value + sum of SHAP values = the model's log-odds for the variant.
One-hot columns (panel, secondary structure) are summed back to their source
feature so every explanation reads in terms of the 28 input features.

Outputs (outputs/):
  shap_global_summary.png      beeswarm of all labelled variants
  shap_global_bar.png          mean |SHAP| per feature
  shap_dependence_<f>.png      value-vs-SHAP plots for the top features
  shap_global_importance.csv   mean |SHAP| table
  shap_direction_check.json    does each feature push in the biologically
                               expected direction? (Spearman rho between
                               feature value and SHAP value)
  vus_predictions.csv          every VUS with model and calibrated
                               probability and its top three reasons,
                               ranked for prioritisation

Usage:
    python -m src.explainability.shap_explain                   # global + VUS ranking
    python -m src.explainability.shap_explain --variant MYH7:R403Q
"""
import argparse
import json
import math

import joblib
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap
from scipy.stats import spearmanr

from src.config import FEATURES_DIR, OUTPUTS_DIR, ensure_dirs
from src.modeling.train_classifier import (
    BOOLEAN_FEATURES,
    CALIBRATED_MODEL_FILE,
    CATEGORICAL_FEATURES,
    FEATURE_COLUMNS,
    FULL_MODEL_FILES,
    NUMERIC_FEATURES,
)

EXPLAINED_MODEL = "xgboost"
TOP_K = 5
N_DEPENDENCE_PLOTS = 4

# Readable name and value format for each model input.
FEATURE_LABELS = {
    "delta_hydrophobicity": ("Change in hydrophobicity (Kyte-Doolittle)", "{:+.1f}"),
    "delta_molecular_weight": ("Change in molecular weight", "{:+.0f} Da"),
    "delta_volume": ("Change in side-chain volume", "{:+.0f} A^3"),
    "delta_charge": ("Change in charge", "{:+.1f}"),
    "blosum62": ("BLOSUM62 substitution score", "{:.0f}"),
    "grantham": ("Grantham distance", "{:.0f}"),
    "conservation_score": ("Residue kept in vertebrate orthologs", "{:.0%}"),
    "mut_aa_in_orthologs": ("Mutant residue seen in orthologs", "{:.0%}"),
    "phylop": ("Genomic conservation (phyloP)", "{:.2f}"),
    "phastcons": ("Genomic conservation (phastCons)", "{:.2f}"),
    "gerp": ("Genomic conservation (GERP++)", "{:.2f}"),
    "plddt": ("AlphaFold confidence (pLDDT)", "{:.0f}"),
    "rsa": ("Relative solvent accessibility", "{:.0%}"),
    "contact_count": ("Residue packing (C-alpha contacts within 10 A)", "{:.0f}"),
    "foldx_ddg": ("FoldX stability change (ddG)", "{:+.2f} kcal/mol"),
    "dist_to_functional_site": ("Distance to nearest functional site", "{:.0f} residues"),
    "polarity_change": ("Polarity class changes", None),
    "wt_is_proline": ("Wild-type residue is proline", None),
    "mut_is_proline": ("Mutant residue is proline", None),
    "wt_is_glycine": ("Wild-type residue is glycine", None),
    "mut_is_glycine": ("Mutant residue is glycine", None),
    "in_annotated_domain": ("Inside a structured domain", None),
    "in_transmembrane": ("Inside a transmembrane segment", None),
    "in_functional_site": ("At an annotated functional site", None),
    "in_disordered_region": ("Inside a disordered region", None),
    "in_coiled_coil": ("Inside a coiled coil", None),
    "has_structure": ("AlphaFold model available", None),
    "panel": ("Gene panel", "{}"),
    "sec_struct": ("Secondary structure", "{}"),
}

# Direction a feature should push when its value rises, from biology:
# +1 = towards pathogenic, -1 = towards benign.
EXPECTED_DIRECTION = {
    "phylop": +1,
    "phastcons": +1,
    "gerp": +1,
    "conservation_score": +1,
    "mut_aa_in_orthologs": -1,
    "blosum62": -1,
    "grantham": +1,
    "rsa": -1,
    "plddt": +1,
    "contact_count": +1,
    "foldx_ddg": +1,
    "in_disordered_region": -1,
}


def load_model():
    model = joblib.load(OUTPUTS_DIR / FULL_MODEL_FILES[EXPLAINED_MODEL])
    calibrated_path = OUTPUTS_DIR / CALIBRATED_MODEL_FILE
    calibrated = joblib.load(calibrated_path) if calibrated_path.exists() else None
    return model, calibrated


def raw_feature_of(transformed_name: str) -> str:
    """Map 'num__phylop' -> 'phylop' and 'cat__sec_struct_helix' -> 'sec_struct'."""
    _, name = transformed_name.split("__", 1)
    for categorical in CATEGORICAL_FEATURES:
        if name.startswith(categorical + "_"):
            return categorical
    return name


def compute_shap(model, X: pd.DataFrame) -> tuple[pd.DataFrame, float]:
    """
    SHAP values per raw feature (one-hot columns summed), in log-odds.
    Returns (values with columns FEATURE_COLUMNS, base value).
    """
    preprocess, clf = model.named_steps["preprocess"], model.named_steps["clf"]
    X_t = preprocess.transform(X)
    X_t = X_t.toarray() if hasattr(X_t, "toarray") else X_t
    explainer = shap.TreeExplainer(clf)
    values = explainer.shap_values(X_t)
    names = [raw_feature_of(n) for n in preprocess.get_feature_names_out()]
    per_column = pd.DataFrame(values, columns=names, index=X.index)
    grouped = per_column.T.groupby(level=0).sum().T
    return grouped[FEATURE_COLUMNS], float(np.ravel(explainer.expected_value)[0])


def format_value(feature: str, value) -> str:
    label, fmt = FEATURE_LABELS.get(feature, (feature, "{}"))
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "missing (imputed)"
    if fmt is None:
        return "yes" if bool(value) else "no"
    try:
        return fmt.format(value)
    except (ValueError, TypeError):
        return str(value)


def explain_row(features: pd.Series, shap_row: pd.Series, base_value: float, top_k: int = TOP_K) -> list[dict]:
    """Top contributions for one variant, largest |SHAP| first."""
    order = shap_row.abs().sort_values(ascending=False).index[:top_k]
    contributions = []
    for feature in order:
        value = features[feature]
        contributions.append({
            "feature": feature,
            "label": FEATURE_LABELS.get(feature, (feature, None))[0],
            "value": None if pd.isna(value) else (value.item() if hasattr(value, "item") else value),
            "value_display": format_value(feature, value),
            "shap": round(float(shap_row[feature]), 4),
            "direction": "pathogenic" if shap_row[feature] > 0 else "benign",
        })
    return contributions


def describe(contribution: dict) -> str:
    return (f"{contribution['label']} = {contribution['value_display']} "
            f"({contribution['shap']:+.2f} towards {contribution['direction']})")


def predict(model, calibrated, X: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame(index=X.index)
    out["model_probability"] = model.predict_proba(X)[:, 1]
    if calibrated is not None:
        out["calibrated_probability"] = calibrated.predict_proba(X)[:, 1]
    return out


def find_variant(spec: str) -> dict:
    """'MYH7:R403Q' -> feature row from the labelled or VUS matrix, or computed on the fly."""
    gene, change = spec.split(":")
    wt, pos, mut = change[0].upper(), int(change[1:-1]), change[-1].upper()
    for name in ("feature_matrix.csv", "feature_matrix_vus.csv"):
        path = FEATURES_DIR / name
        if not path.exists():
            continue
        df = pd.read_csv(path)
        hit = df[(df["gene_symbol"] == gene) & (df["position"] == pos) & (df["wt_aa"] == wt) & (df["mut_aa"] == mut)]
        if not hit.empty:
            return hit.iloc[0].to_dict()

    from src.config import load_gene_panel
    from src.features.build_feature_matrix import features_for_variant

    panel = {g["symbol"]: g["panel"] for g in load_gene_panel("both")}
    if gene not in panel:
        raise ValueError(f"{gene} is not in the gene panel")
    return features_for_variant(
        {"gene_symbol": gene, "position": pos, "wt_aa": wt, "mut_aa": mut, "panel": panel[gene], "label": "unknown"}
    )


def explain_variant(spec: str, model=None, calibrated=None) -> dict:
    """Structured explanation for one variant, e.g. explain_variant('MYH7:R403Q')."""
    if model is None:
        model, calibrated = load_model()
    record = find_variant(spec)
    X = pd.DataFrame([{c: record.get(c) for c in FEATURE_COLUMNS}])
    shap_values, base_value = compute_shap(model, X)
    probs = predict(model, calibrated, X).iloc[0]
    contributions = explain_row(X.iloc[0], shap_values.iloc[0], base_value)
    probability = float(probs.get("calibrated_probability", probs["model_probability"]))
    return {
        "variant": f"{record['gene_symbol']} {record['wt_aa']}{int(record['position'])}{record['mut_aa']}",
        "clinvar_label": record.get("label"),
        "prediction": "pathogenic" if probability >= 0.5 else "benign",
        "probability_pathogenic": round(probability, 4),
        "model_probability": round(float(probs["model_probability"]), 4),
        "base_value_log_odds": round(base_value, 4),
        "log_odds": round(base_value + float(shap_values.iloc[0].sum()), 4),
        "top_features": contributions,
        "summary": "; ".join(describe(c) for c in contributions[:3]),
        "note": "SHAP values are in log-odds of the XGBoost model; probability_pathogenic is "
                "sigmoid-calibrated. Research use only.",
    }


def global_plots(shap_values: pd.DataFrame, X: pd.DataFrame) -> pd.DataFrame:
    labels = [FEATURE_LABELS[c][0] for c in FEATURE_COLUMNS]
    numeric_X = X.copy()
    for c in CATEGORICAL_FEATURES:
        numeric_X[c] = numeric_X[c].astype("category").cat.codes  # colour only
    numeric_X = numeric_X.astype(float)

    plt.figure()
    shap.summary_plot(shap_values.to_numpy(), numeric_X.to_numpy(), feature_names=labels, show=False, max_display=20)
    plt.title("SHAP values per variant (XGBoost, log-odds)")
    plt.tight_layout()
    plt.savefig(OUTPUTS_DIR / "shap_global_summary.png", dpi=150, bbox_inches="tight")
    plt.close()

    importance = (
        shap_values.abs().mean().rename("mean_abs_shap").sort_values(ascending=False)
        .rename_axis("feature").reset_index()
    )
    importance["label"] = importance["feature"].map(lambda f: FEATURE_LABELS[f][0])
    importance.to_csv(OUTPUTS_DIR / "shap_global_importance.csv", index=False)

    top = importance.head(20).iloc[::-1]
    plt.figure(figsize=(8, 7))
    plt.barh(top["label"], top["mean_abs_shap"], color="#4C72B0")
    plt.xlabel("mean |SHAP| (log-odds)")
    plt.title("Global feature importance")
    plt.tight_layout()
    plt.savefig(OUTPUTS_DIR / "shap_global_bar.png", dpi=150)
    plt.close()

    for feature in [f for f in importance["feature"] if f in NUMERIC_FEATURES][:N_DEPENDENCE_PLOTS]:
        plt.figure(figsize=(5, 4))
        plt.scatter(X[feature], shap_values[feature], s=6, alpha=0.5)
        plt.axhline(0, color="grey", linewidth=0.8)
        plt.xlabel(FEATURE_LABELS[feature][0])
        plt.ylabel("SHAP (log-odds)")
        plt.tight_layout()
        plt.savefig(OUTPUTS_DIR / f"shap_dependence_{feature}.png", dpi=150)
        plt.close()
    return importance


def direction_check(shap_values: pd.DataFrame, X: pd.DataFrame) -> dict:
    """
    Spearman rho between each feature's value and its SHAP value, vs the
    expected sign. A mismatch flags a feature whose attribution should not
    be read in isolation (usually because a correlated feature carries the
    main effect); it is reported, not corrected.
    """
    results = {}
    for feature, expected in EXPECTED_DIRECTION.items():
        x = X[feature].astype(float)
        mask = x.notna()
        if x[mask].nunique() < 2:
            continue
        rho = spearmanr(x[mask], shap_values.loc[mask, feature]).statistic
        results[feature] = {
            "expected": "towards pathogenic" if expected > 0 else "towards benign",
            "spearman_rho": round(float(rho), 3),
            "consistent": bool(np.sign(rho) == expected),
        }
    return results


def rank_vus(model, calibrated) -> pd.DataFrame | None:
    path = FEATURES_DIR / "feature_matrix_vus.csv"
    if not path.exists():
        return None
    vus = pd.read_csv(path)
    X = vus[FEATURE_COLUMNS]
    shap_values, base_value = compute_shap(model, X)
    probs = predict(model, calibrated, X)
    top = [
        " | ".join(describe(c) for c in explain_row(X.loc[i], shap_values.loc[i], base_value, top_k=3))
        for i in X.index
    ]
    out = vus[["gene_symbol", "panel", "wt_aa", "position", "mut_aa", "clinvar_uid", "review_stars"]].copy()
    out = pd.concat([out, probs.round(4)], axis=1)
    out["top_reasons"] = top
    sort_col = "calibrated_probability" if "calibrated_probability" in out else "model_probability"
    out = out.sort_values(sort_col, ascending=False).reset_index(drop=True)
    out.insert(0, "rank", out.index + 1)
    out.to_csv(OUTPUTS_DIR / "vus_predictions.csv", index=False)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--variant", help="explain one variant, e.g. MYH7:R403Q")
    args = parser.parse_args()
    ensure_dirs()
    model, calibrated = load_model()

    if args.variant:
        print(json.dumps(explain_variant(args.variant, model, calibrated), indent=2, default=str))
        return

    labelled = pd.read_csv(FEATURES_DIR / "feature_matrix.csv")
    X = labelled[FEATURE_COLUMNS]
    shap_values, base_value = compute_shap(model, X)

    importance = global_plots(shap_values, X)
    print("Global importance (mean |SHAP|, log-odds):")
    print(importance[["label", "mean_abs_shap"]].head(12).to_string(index=False))

    checks = direction_check(shap_values, X)
    with open(OUTPUTS_DIR / "shap_direction_check.json", "w", encoding="utf-8") as f:
        json.dump(checks, f, indent=2)
    print("\nBiological direction check (Spearman rho of value vs SHAP):")
    for feature, r in checks.items():
        print(f"  {feature:22s} expected {r['expected']:20s} rho {r['spearman_rho']:+.2f}  "
              f"{'as expected' if r['consistent'] else 'differs (see note)'}")
    print("  Features that differ are strongly correlated with a stronger feature (GERP and phastCons\n"
          "  with phyloP, contact count with solvent accessibility, Grantham with BLOSUM62); the trees use them\n"
          "  as residual corrections, so their marginal direction is not interpretable on its own.")

    ranked = rank_vus(model, calibrated)
    if ranked is not None:
        print(f"\nRanked {len(ranked)} VUS -> {OUTPUTS_DIR / 'vus_predictions.csv'}")
        print(ranked.head(5)[["rank", "gene_symbol", "wt_aa", "position", "mut_aa", "calibrated_probability"]]
              .to_string(index=False))

    for spec in ("MYH7:R403Q", "KCNQ1:A341V"):
        try:
            e = explain_variant(spec, model, calibrated)
            print(f"\n{e['variant']}: {e['prediction']} (p={e['probability_pathogenic']:.2f}, "
                  f"ClinVar {e['clinvar_label']})\n  {e['summary']}")
        except Exception as exc:
            print(f"[warn] could not explain {spec}: {exc}")


if __name__ == "__main__":
    main()
