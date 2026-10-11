"""
Step 4: Variant classification and evaluation.

Models: Logistic Regression, Random Forest and XGBoost, all with class
weighting (no synthetic oversampling).

Evaluation (src/modeling/evaluate.py):
  - Leave-one-gene-out cross-validation: each gene is predicted by a
    model that never saw any variant from it. This is the headline
    estimate.
  - Within-gene stratified cross-validation (rows split at random), for
    comparison. The gap between the two shows how much a model relies on
    gene identity rather than variant-level biology.
  - Published predictors (SIFT, PolyPhen-2, AlphaMissense, REVEL) scored
    on the same variants as baselines. REVEL and AlphaMissense were
    trained or calibrated with ClinVar data, so their scores here are
    optimistic.

Saved models (outputs/):
  model_{random_forest,xgboost,logistic_regression}.joblib
      trained on the training genes of gene_level_split only; the API
      evaluates them on the held-out genes.
  model_*_full.joblib
      refit on all labelled variants; used for VUS scoring and SHAP.
  model_xgboost_full_calibrated.joblib
      sigmoid-calibrated XGBoost (gene-grouped folds) for the confidence
      shown in reports.
The RF and XGBoost pipelines keep the Pipeline(preprocess, clf) layout with
a SimpleImputer as the "num" transformer, which the API depends on.

Usage:
    python -m src.modeling.train_classifier
"""
import json

import joblib
import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.calibration import CalibratedClassifierCV
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupShuffleSplit, StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from src.config import FEATURES_DIR, OUTPUTS_DIR, RANDOM_STATE, ensure_dirs
from src.modeling.evaluate import evaluate_classifier, run_evaluation

# relative_position is computed but not used: in gene-grouped CV it lowered
# ROC-AUC for all three models (it lets the model memorise where in a
# protein the pathogenic variants of the training genes cluster).
NUMERIC_FEATURES = [
    "delta_hydrophobicity",
    "delta_molecular_weight",
    "delta_volume",
    "delta_charge",
    "blosum62",
    "grantham",
    "conservation_score",
    "mut_aa_in_orthologs",
    "phylop",
    "phastcons",
    "gerp",
    "plddt",
    "rsa",
    "contact_count",
    "dist_to_functional_site",
]
BOOLEAN_FEATURES = [
    "polarity_change",
    "wt_is_proline",
    "mut_is_proline",
    "wt_is_glycine",
    "mut_is_glycine",
    "in_annotated_domain",
    "in_transmembrane",
    "in_functional_site",
    "in_disordered_region",
    "in_coiled_coil",
    "has_structure",
]
CATEGORICAL_FEATURES = ["panel", "sec_struct"]
FEATURE_COLUMNS = NUMERIC_FEATURES + BOOLEAN_FEATURES + CATEGORICAL_FEATURES

MODEL_FILES = {
    "random_forest": "model_random_forest.joblib",
    "xgboost": "model_xgboost.joblib",
    "logistic_regression": "model_logistic_regression.joblib",
}
FULL_MODEL_FILES = {name: path.replace(".joblib", "_full.joblib") for name, path in MODEL_FILES.items()}
CALIBRATED_MODEL_FILE = "model_xgboost_full_calibrated.joblib"


def load_feature_matrix() -> pd.DataFrame:
    path = FEATURES_DIR / "feature_matrix.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found — run `python -m src.features.build_feature_matrix` first."
        )
    return training_rows(pd.read_csv(path))


def training_rows(df: pd.DataFrame) -> pd.DataFrame:
    """
    Rows the models train on: every pathogenic/benign variant. Missing
    feature values are imputed inside the pipeline, so no rows are dropped.
    The API imports this so its gene split matches the one used here.
    """
    df = df[df["label"].isin(["pathogenic", "benign"])].reset_index(drop=True)
    df["y"] = (df["label"] == "pathogenic").astype(int)
    return df


def gene_level_split(df: pd.DataFrame, test_size: float = 0.25, random_state: int = RANDOM_STATE):
    """Single gene-held-out split. Kept for the API; evaluation uses grouped CV."""
    splitter = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=random_state)
    train_idx, test_idx = next(splitter.split(df, groups=df["gene_symbol"]))
    return df.iloc[train_idx].reset_index(drop=True), df.iloc[test_idx].reset_index(drop=True)


def build_preprocessor(scale: bool = False) -> ColumnTransformer:
    numeric = SimpleImputer(strategy="median")
    if scale:
        numeric = Pipeline([("impute", numeric), ("scale", StandardScaler())])
    return ColumnTransformer(
        transformers=[
            ("num", numeric, NUMERIC_FEATURES + BOOLEAN_FEATURES),
            ("cat", OneHotEncoder(handle_unknown="ignore"), CATEGORICAL_FEATURES),
        ]
    )


def make_logistic_regression(y_train=None) -> Pipeline:
    return Pipeline([
        ("preprocess", build_preprocessor(scale=True)),
        ("clf", LogisticRegression(C=0.5, class_weight="balanced", max_iter=2000)),
    ])


def make_random_forest(y_train=None) -> Pipeline:
    return Pipeline([
        ("preprocess", build_preprocessor()),
        ("clf", RandomForestClassifier(
            n_estimators=500,
            min_samples_leaf=3,
            max_features="sqrt",
            class_weight="balanced",
            random_state=RANDOM_STATE,
            n_jobs=-1,
        )),
    ])


def make_xgboost(y_train=None) -> Pipeline:
    scale_pos_weight = 1.0
    if y_train is not None:
        y_train = np.asarray(y_train)
        scale_pos_weight = (y_train == 0).sum() / max((y_train == 1).sum(), 1)
    return Pipeline([
        ("preprocess", build_preprocessor()),
        ("clf", xgb.XGBClassifier(
            n_estimators=300,
            max_depth=3,
            learning_rate=0.05,
            subsample=0.8,
            colsample_bytree=0.8,
            min_child_weight=3,
            scale_pos_weight=scale_pos_weight,
            eval_metric="logloss",
            random_state=RANDOM_STATE,
            n_jobs=-1,
        )),
    ])


MODEL_FACTORIES = {
    "logistic_regression": make_logistic_regression,
    "random_forest": make_random_forest,
    "xgboost": make_xgboost,
}


def train_random_forest(X_train, y_train) -> Pipeline:
    return make_random_forest(y_train).fit(X_train, y_train)


def train_xgboost(X_train, y_train) -> Pipeline:
    return make_xgboost(y_train).fit(X_train, y_train)


def calibrated_xgboost(X, y, groups) -> CalibratedClassifierCV:
    folds = list(StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE).split(X, y, groups))
    model = CalibratedClassifierCV(make_xgboost(y), method="sigmoid", cv=folds)
    return model.fit(X, y)


def main():
    ensure_dirs()
    df = load_feature_matrix()
    X, y, groups = df[FEATURE_COLUMNS], df["y"], df["gene_symbol"]
    print(f"{len(df)} labelled variants, {groups.nunique()} genes, "
          f"{int(y.sum())} pathogenic / {int((1 - y).sum())} benign")

    results = run_evaluation(df, FEATURE_COLUMNS, MODEL_FACTORIES)
    with open(OUTPUTS_DIR / "metrics.json", "w", encoding="utf-8") as f:
        json.dump(results["metrics"], f, indent=2)
    # Read as-is by the API's /ablation endpoint.
    with open(OUTPUTS_DIR / "ablation_results.json", "w", encoding="utf-8") as f:
        json.dump(results["metrics"]["ablation_gene_grouped"], f, indent=2)
    results["oof"].to_csv(OUTPUTS_DIR / "oof_predictions.csv", index=False)

    # Split models: trained on the training genes of gene_level_split only.
    # The API scores these on the held-out genes, so they must not see them.
    train_df, test_df = gene_level_split(df)
    print(f"\nHeld-out split for the API: test genes {sorted(test_df['gene_symbol'].unique())}")
    for name, factory in MODEL_FACTORIES.items():
        model = factory(train_df["y"]).fit(train_df[FEATURE_COLUMNS], train_df["y"])
        joblib.dump(model, OUTPUTS_DIR / MODEL_FILES[name])
        evaluate_classifier(model, test_df[FEATURE_COLUMNS], test_df["y"], name=f"{name} (held-out genes)")

    # Full models: trained on every labelled variant, used to score VUS and
    # for SHAP explanations. XGBoost is the explained model; it is within
    # ~0.02 ROC-AUC of logistic regression and SHAP's TreeExplainer is exact for it.
    for name, factory in MODEL_FACTORIES.items():
        joblib.dump(factory(y).fit(X, y), OUTPUTS_DIR / FULL_MODEL_FILES[name])
    joblib.dump(calibrated_xgboost(X, y, groups), OUTPUTS_DIR / CALIBRATED_MODEL_FILE)
    print(f"\nModels saved to {OUTPUTS_DIR}")


if __name__ == "__main__":
    main()
