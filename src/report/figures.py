"""
Figures for the project write-up (docs/report.md).

Reads the pipeline outputs and writes PNGs to docs/figures/, which is
tracked in git so the report renders without regenerating data.

Usage:
    python -m src.report.figures
"""
import json
import shutil

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve

from src.config import FEATURES_DIR, OUTPUTS_DIR, PROJECT_ROOT, REPORTS_DIR

FIG_DIR = PROJECT_ROOT / "docs" / "figures"
MODELS = {"logistic_regression": "Logistic regression", "random_forest": "Random forest", "xgboost": "XGBoost"}
PATHOGENIC, BENIGN, NEUTRAL = "#b2182b", "#2166ac", "#7f7f7f"


def _save(fig, name: str):
    fig.tight_layout()
    fig.savefig(FIG_DIR / name, dpi=150)
    plt.close(fig)
    print(f"wrote docs/figures/{name}")


def labelled_with_oof() -> pd.DataFrame:
    features = pd.read_csv(FEATURES_DIR / "feature_matrix.csv")
    oof = pd.read_csv(OUTPUTS_DIR / "oof_predictions.csv")
    return oof.merge(features[["clinvar_uid", "alphamissense", "revel", "foldx_ddg"]], on="clinvar_uid", how="left")


def roc_figure(df: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(5.5, 5))
    metrics = json.loads((OUTPUTS_DIR / "metrics.json").read_text())["gene_grouped_cv"]
    for key, name in MODELS.items():
        fpr, tpr, _ = roc_curve(df["y"], df[f"{key}_gene_grouped"])
        ax.plot(fpr, tpr, label=f"{name} ({metrics[key]['roc_auc']['mean']:.3f})")
    for col, name in (("alphamissense", "AlphaMissense*"), ("revel", "REVEL*")):
        sub = df.dropna(subset=[col])
        fpr, tpr, _ = roc_curve(sub["y"], sub[col])
        auc = roc_auc_score(sub["y"], sub[col])
        ax.plot(fpr, tpr, linestyle="--", label=f"{name} ({auc:.3f})")
    ax.plot([0, 1], [0, 1], color="lightgrey", linewidth=1)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("Leave-one-gene-out ROC (out-of-fold)")
    ax.legend(fontsize=8, loc="lower right")
    ax.text(0.02, 0.98, "* trained or calibrated with ClinVar", transform=ax.transAxes, fontsize=7, va="top")
    _save(fig, "fig2_roc_leave_one_gene_out.png")


def ablation_figure():
    ablation = json.loads((OUTPUTS_DIR / "ablation_results.json").read_text())["models"]
    settings = list(ablation["logistic_regression"].keys())
    x = np.arange(len(settings))
    fig, ax = plt.subplots(figsize=(7, 3.8))
    for i, model in enumerate(("logistic_regression", "xgboost")):
        values = [ablation[model][s]["roc_auc"]["mean"] for s in settings]
        ax.bar(x + (i - 0.5) * 0.38, values, width=0.38, label=MODELS[model])
    ax.set_xticks(x, [s.replace("without_", "minus ").replace("_", " ") for s in settings], rotation=20, ha="right")
    ax.set_ylim(0.75, 0.95)
    ax.set_ylabel("ROC AUC (leave one gene out)")
    ax.set_title("Feature-group ablation")
    ax.legend(fontsize=8)
    _save(fig, "fig3_ablation.png")


def per_gene_figure():
    per_gene = pd.DataFrame(json.loads((OUTPUTS_DIR / "metrics.json").read_text())["per_gene"])
    per_gene = per_gene.sort_values("n", ascending=False)
    fig, ax = plt.subplots(figsize=(8, 3.8))
    colours = [NEUTRAL if pd.isna(a) else PATHOGENIC for a in per_gene["roc_auc"]]
    ax.bar(per_gene["gene_symbol"], per_gene["roc_auc"].fillna(0), color=colours)
    for i, row in enumerate(per_gene.itertuples()):
        label = "n/a" if pd.isna(row.roc_auc) else f"{row.roc_auc:.2f}"
        ax.text(i, (0 if pd.isna(row.roc_auc) else row.roc_auc) + 0.02, f"{label}\n{row.pathogenic}/{row.benign}",
                ha="center", fontsize=6)
    ax.set_ylim(0, 1.2)
    ax.set_ylabel("ROC AUC, gene held out")
    ax.set_title("Per-gene performance (labels: AUC, pathogenic/benign); n/a = one class only")
    ax.tick_params(axis="x", rotation=45, labelsize=8)
    _save(fig, "fig4_per_gene.png")


def ddg_figure(df: pd.DataFrame):
    sub = df.dropna(subset=["foldx_ddg"])
    fig, ax = plt.subplots(figsize=(5.5, 3.8))
    bins = np.linspace(-3, 8, 45)
    for y, name, colour in ((0, "Benign", BENIGN), (1, "Pathogenic", PATHOGENIC)):
        values = sub.loc[sub["y"] == y, "foldx_ddg"].clip(-3, 8)
        ax.hist(values, bins=bins, density=True, alpha=0.55, color=colour,
                label=f"{name} (n={len(values)}, median {values.median():.2f})")
    ax.axvline(1, color="black", linewidth=0.8, linestyle=":")
    ax.set_xlabel("FoldX ΔΔG (kcal/mol, clipped to [-3, 8])")
    ax.set_ylabel("Density")
    ax.set_title("Stability change by ClinVar class")
    ax.legend(fontsize=8)
    _save(fig, "fig6_foldx_ddg.png")


def vus_figure():
    table = pd.read_csv(REPORTS_DIR / "vus_report.csv")
    order = [
        "strong model support for benign", "moderate model support for benign", "inconclusive",
        "moderate model support for pathogenic", "strong model support for pathogenic",
    ]
    counts = table["tier"].value_counts().reindex(order, fill_value=0)
    am = pd.cut(table["alphamissense"], [-1, 0.34, 0.564, 2], labels=["benign", "ambiguous", "pathogenic"])
    agree = pd.crosstab(table["tier"], am, normalize="index").reindex(order)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.8))
    labels = [t.replace("model support for ", "").replace("moderate", "mod.") for t in order]
    a1.bar(labels, counts.values, color=[BENIGN, "#92c5de", NEUTRAL, "#f4a582", PATHOGENIC])
    a1.set_title(f"{len(table):,} VUS by model-support tier")
    a1.tick_params(axis="x", rotation=30, labelsize=8)
    bottom = np.zeros(len(order))
    for col, colour in (("benign", BENIGN), ("ambiguous", NEUTRAL), ("pathogenic", PATHOGENIC)):
        a2.bar(labels, agree[col].values, bottom=bottom, color=colour, label=f"AlphaMissense {col}")
        bottom += agree[col].values
    a2.set_title("AlphaMissense class within each tier")
    a2.tick_params(axis="x", rotation=30, labelsize=8)
    a2.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=3)
    _save(fig, "fig7_vus_tiers.png")


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    df = labelled_with_oof()
    roc_figure(df)
    ablation_figure()
    per_gene_figure()
    ddg_figure(df)
    vus_figure()
    for src, dest in (("shap_global_bar.png", "fig5a_shap_importance.png"),
                      ("shap_global_summary.png", "fig5b_shap_summary.png")):
        shutil.copy(OUTPUTS_DIR / src, FIG_DIR / dest)
        print(f"copied docs/figures/{dest}")


if __name__ == "__main__":
    main()
