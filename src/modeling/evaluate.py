"""
Model evaluation utilities.

Reports threshold metrics (accuracy, balanced accuracy, precision, recall,
specificity, F1, MCC) and ranking metrics (ROC-AUC, PR-AUC), plus the Brier
score. PR-AUC and MCC are included because ROC-AUC and accuracy look
misleadingly strong when one class dominates (here ~85% pathogenic).

run_evaluation() produces everything written to outputs/metrics.json:
leave-one-gene-out CV, within-gene CV, per-gene results, published-predictor
baselines and a feature-group ablation. Each pooled metric carries a 95%
bootstrap confidence interval (variants resampled with replacement).
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    classification_report,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    precision_recall_curve,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import LeaveOneGroupOut, StratifiedKFold

from src.config import OUTPUTS_DIR

N_SPLITS = 5
SEEDS = [0, 1, 2, 3, 4]  # repeats of the within-gene CV; leave-one-gene-out is deterministic
N_BOOTSTRAP = 1000
MIN_PER_CLASS_FOR_GENE_AUC = 3

# Published predictors: column -> True if a higher score means more damaging
BASELINES = {"alphamissense": True, "revel": True, "polyphen_score": True, "sift_score": False}

# Feature groups for the ablation study (columns removed from the full set)
ABLATION_GROUPS = {
    "conservation": ["conservation_score", "mut_aa_in_orthologs", "phylop", "phastcons", "gerp"],
    "structure": ["plddt", "rsa", "contact_count", "has_structure", "sec_struct", "foldx_ddg"],
    "stability": ["foldx_ddg"],
    "annotation": [
        "in_annotated_domain", "in_transmembrane", "in_functional_site",
        "dist_to_functional_site", "in_disordered_region", "in_coiled_coil",
    ],
    "panel": ["panel"],
}


def metrics_from_scores(y_true, proba, threshold: float = 0.5) -> dict:
    y_true = np.asarray(y_true)
    proba = np.asarray(proba)
    y_pred = (proba >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    return {
        "accuracy": accuracy_score(y_true, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "precision": precision_score(y_true, y_pred, zero_division=0),
        "recall": recall_score(y_true, y_pred, zero_division=0),
        "specificity": tn / (tn + fp) if (tn + fp) else float("nan"),
        "f1": f1_score(y_true, y_pred, zero_division=0),
        "mcc": matthews_corrcoef(y_true, y_pred),
        "roc_auc": roc_auc_score(y_true, proba),
        "pr_auc": average_precision_score(y_true, proba),
        "brier": brier_score_loss(y_true, proba),
    }


def evaluate_classifier(model, X_test, y_test, name: str = "model") -> dict:
    """Score a fitted model on a held-out set (used by the API's split view)."""
    y_proba = model.predict_proba(X_test)[:, 1]
    metrics = metrics_from_scores(y_test, y_proba)
    print(f"[{name}] " + " | ".join(f"{k}={v:.3f}" for k, v in metrics.items()))
    y_pred = (y_proba >= 0.5).astype(int)
    print(confusion_matrix(y_test, y_pred, labels=[0, 1]))
    print(classification_report(y_test, y_pred, labels=[0, 1], target_names=["benign", "pathogenic"], zero_division=0))
    return metrics


def _splits(scheme: str, X, y, groups, seed: int):
    if scheme == "gene_grouped":
        return LeaveOneGroupOut().split(X, y, groups)
    return StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=seed).split(X, y)


def cross_validate(df: pd.DataFrame, feature_cols: list[str], factory, scheme: str) -> tuple[list[dict], np.ndarray]:
    """
    Leave-one-gene-out ("gene_grouped") or repeated stratified K-fold
    ("within_gene"). Metrics are computed on pooled out-of-fold predictions
    (many genes have only one class, so per-fold AUC is undefined).
    Returns per-repeat metrics and the out-of-fold probabilities averaged
    over repeats.
    """
    X, y, groups = df[feature_cols], df["y"].to_numpy(), df["gene_symbol"].to_numpy()
    per_seed, oof_sum = [], np.zeros(len(df))
    seeds = SEEDS if scheme == "within_gene" else [0]
    for seed in seeds:
        oof = np.zeros(len(df))
        for train_idx, test_idx in _splits(scheme, X, y, groups, seed):
            model = factory(y[train_idx]).fit(X.iloc[train_idx], y[train_idx])
            oof[test_idx] = model.predict_proba(X.iloc[test_idx])[:, 1]
        per_seed.append(metrics_from_scores(y, oof))
        oof_sum += oof
    return per_seed, oof_sum / len(seeds)


def bootstrap_ci(y_true, proba, n: int = N_BOOTSTRAP, seed: int = 0) -> dict:
    """95% percentile intervals for the pooled metrics, resampling variants."""
    y_true, proba = np.asarray(y_true), np.asarray(proba)
    rng = np.random.default_rng(seed)
    samples = []
    for _ in range(n):
        idx = rng.integers(0, len(y_true), len(y_true))
        if len(set(y_true[idx])) == 2:
            samples.append(metrics_from_scores(y_true[idx], proba[idx]))
    frame = pd.DataFrame(samples)
    return {k: [round(frame[k].quantile(0.025), 4), round(frame[k].quantile(0.975), 4)] for k in frame.columns}


def summarize(per_seed: list[dict], ci: dict | None = None) -> dict:
    """Mean and spread across repeats; ci95 from the bootstrap when given."""
    frame = pd.DataFrame(per_seed)
    out = {k: {"mean": round(frame[k].mean(), 4), "std": round(frame[k].std(ddof=0), 4)} for k in frame.columns}
    for k, interval in (ci or {}).items():
        out[k]["ci95"] = interval
    return out


def per_gene_metrics(df: pd.DataFrame, proba: np.ndarray) -> list[dict]:
    rows = []
    for gene, idx in df.groupby("gene_symbol").indices.items():
        y = df["y"].to_numpy()[idx]
        n_pos, n_neg = int(y.sum()), int(len(y) - y.sum())
        enough = min(n_pos, n_neg) >= MIN_PER_CLASS_FOR_GENE_AUC
        rows.append({
            "gene_symbol": gene,
            "n": len(y),
            "pathogenic": n_pos,
            "benign": n_neg,
            "roc_auc": round(roc_auc_score(y, proba[idx]), 4) if enough else None,
            "accuracy": round(accuracy_score(y, (proba[idx] >= 0.5).astype(int)), 4),
        })
    return rows


def baseline_metrics(df: pd.DataFrame, model_proba: np.ndarray) -> dict:
    """Published predictors vs the model, each on the variants that predictor covers."""
    out = {}
    y = df["y"].to_numpy()
    for column, higher_is_damaging in BASELINES.items():
        if column not in df.columns:
            continue
        score = pd.to_numeric(df[column], errors="coerce").to_numpy()
        mask = ~np.isnan(score)
        if mask.sum() == 0 or len(set(y[mask])) < 2:
            continue
        s = score[mask] if higher_is_damaging else -score[mask]
        out[column] = {
            "coverage": round(float(mask.mean()), 4),
            "roc_auc": round(roc_auc_score(y[mask], s), 4),
            "pr_auc": round(average_precision_score(y[mask], s), 4),
            "model_roc_auc_same_variants": round(roc_auc_score(y[mask], model_proba[mask]), 4),
            "model_pr_auc_same_variants": round(average_precision_score(y[mask], model_proba[mask]), 4),
        }
    return out


def ablation(df: pd.DataFrame, feature_cols: list[str], factory) -> dict:
    """Gene-grouped CV with each feature group removed in turn."""
    per_seed, _ = cross_validate(df, feature_cols, factory, "gene_grouped")
    out = {"all_features": summarize(per_seed)}
    for group, removed in ABLATION_GROUPS.items():
        cols = [c for c in feature_cols if c not in removed]
        per_seed, _ = cross_validate(df, cols, _with_columns(factory, cols), "gene_grouped")
        out[f"without_{group}"] = summarize(per_seed)
    return out


def _with_columns(factory, cols: list[str]):
    """Wrap a model factory so its preprocessor only uses the given columns."""
    def make(y_train):
        model = factory(y_train)
        pre = model.named_steps["preprocess"]
        pre.transformers = [
            (name, transformer, [c for c in columns if c in cols]) for name, transformer, columns in pre.transformers
        ]
        return model
    return make


def plot_curves(df: pd.DataFrame, oof: dict[str, np.ndarray], out_path) -> None:
    y = df["y"].to_numpy()
    fig, (ax_roc, ax_pr) = plt.subplots(1, 2, figsize=(11, 4.5))
    for name, proba in oof.items():
        fpr, tpr, _ = roc_curve(y, proba)
        prec, rec, _ = precision_recall_curve(y, proba)
        ax_roc.plot(fpr, tpr, label=f"{name} (AUC {roc_auc_score(y, proba):.3f})")
        ax_pr.plot(rec, prec, label=f"{name} (AP {average_precision_score(y, proba):.3f})")
    ax_roc.plot([0, 1], [0, 1], color="grey", linestyle="--", linewidth=0.8)
    ax_pr.axhline(y.mean(), color="grey", linestyle="--", linewidth=0.8, label=f"prevalence {y.mean():.2f}")
    ax_roc.set(xlabel="False positive rate", ylabel="True positive rate", title="ROC, leave-one-gene-out")
    ax_pr.set(xlabel="Recall", ylabel="Precision", title="Precision-recall, leave-one-gene-out")
    ax_roc.legend(loc="lower right", fontsize=8)
    ax_pr.legend(loc="lower left", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def run_evaluation(df: pd.DataFrame, feature_cols: list[str], factories: dict) -> dict:
    metrics = {
        "n_variants": len(df),
        "n_genes": int(df["gene_symbol"].nunique()),
        "class_counts": {"pathogenic": int(df["y"].sum()), "benign": int((1 - df["y"]).sum())},
        "cv": {
            "gene_grouped": "leave-one-gene-out",
            "within_gene": f"stratified {N_SPLITS}-fold, {len(SEEDS)} repeats",
            "ci95": f"bootstrap over variants, {N_BOOTSTRAP} resamples",
        },
        "gene_grouped_cv": {},
        "within_gene_cv": {},
    }
    oof = pd.DataFrame({
        "clinvar_uid": df["clinvar_uid"], "gene_symbol": df["gene_symbol"],
        "variant": df["wt_aa"] + df["position"].astype(str) + df["mut_aa"], "y": df["y"],
    })
    grouped_oof = {}
    for name, factory in factories.items():
        for scheme in ("gene_grouped", "within_gene"):
            per_seed, proba = cross_validate(df, feature_cols, factory, scheme)
            metrics[f"{scheme}_cv"][name] = summarize(per_seed, bootstrap_ci(df["y"], proba))
            oof[f"{name}_{scheme}"] = proba
            if scheme == "gene_grouped":
                grouped_oof[name] = proba
            m = metrics[f"{scheme}_cv"][name]
            lo, hi = m["roc_auc"]["ci95"]
            print(f"{name:20s} {scheme:13s} ROC-AUC {m['roc_auc']['mean']:.3f} [{lo:.3f}, {hi:.3f}]  "
                  f"PR-AUC {m['pr_auc']['mean']:.3f}  MCC {m['mcc']['mean']:.3f}  "
                  f"bal.acc {m['balanced_accuracy']['mean']:.3f}")

    best = max(grouped_oof, key=lambda n: metrics["gene_grouped_cv"][n]["roc_auc"]["mean"])
    metrics["best_model_gene_grouped"] = best
    metrics["per_gene"] = per_gene_metrics(df, grouped_oof[best])
    metrics["baselines"] = baseline_metrics(df, grouped_oof[best])
    metrics["ablation_gene_grouped"] = {
        "removed_groups": ABLATION_GROUPS,
        "models": {name: ablation(df, feature_cols, factories[name]) for name in ("logistic_regression", "xgboost")},
    }

    plot_curves(df, grouped_oof, OUTPUTS_DIR / "roc_pr_curves.png")
    print("\nBaselines (each on the variants it covers):")
    for k, v in metrics["baselines"].items():
        print(f"  {k:15s} coverage {v['coverage']:.2f}  ROC-AUC {v['roc_auc']:.3f}  "
              f"(model on same variants {v['model_roc_auc_same_variants']:.3f})")
    print("\nAblation (gene-grouped ROC-AUC):")
    for name, rows in metrics["ablation_gene_grouped"]["models"].items():
        for k, v in rows.items():
            print(f"  {name:20s} {k:25s} {v['roc_auc']['mean']:.3f}")
    return {"metrics": metrics, "oof": oof}
