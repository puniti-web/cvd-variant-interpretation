"""Fast tests for the modeling pipeline on a small synthetic feature matrix."""
import numpy as np
import pandas as pd

from src.modeling.evaluate import bootstrap_ci, cross_validate, metrics_from_scores, _with_columns
from src.modeling.train_classifier import (
    BOOLEAN_FEATURES,
    CATEGORICAL_FEATURES,
    FEATURE_COLUMNS,
    NUMERIC_FEATURES,
    make_logistic_regression,
    make_xgboost,
)


def synthetic_matrix(n_genes: int = 4, per_gene: int = 30, seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    n = n_genes * per_gene
    y = rng.integers(0, 2, n)
    df = pd.DataFrame({c: rng.normal(size=n) for c in NUMERIC_FEATURES})
    df["phylop"] = y * 3 + rng.normal(size=n)  # informative feature
    for c in BOOLEAN_FEATURES:
        df[c] = rng.integers(0, 2, n).astype(bool)
    df["panel"] = rng.choice(["cardiomyopathy", "arrhythmia"], n)
    df["sec_struct"] = rng.choice(["helix", "strand", "coil"], n)
    df.loc[::7, "conservation_score"] = np.nan  # imputation path
    df["gene_symbol"] = np.repeat([f"G{i}" for i in range(n_genes)], per_gene)
    df["y"] = y
    return df


def test_feature_lists_have_no_duplicates_or_identifiers():
    assert len(FEATURE_COLUMNS) == len(set(FEATURE_COLUMNS))
    for leaky in ("position", "relative_position", "gene_symbol", "label", "revel", "alphamissense"):
        assert leaky not in FEATURE_COLUMNS
    assert set(CATEGORICAL_FEATURES) <= set(FEATURE_COLUMNS)


def test_metrics_from_scores_perfect_ranking():
    m = metrics_from_scores([0, 0, 1, 1], [0.1, 0.2, 0.8, 0.9])
    assert m["roc_auc"] == 1.0 and m["mcc"] == 1.0 and m["specificity"] == 1.0


def test_leave_one_gene_out_learns_informative_feature():
    df = synthetic_matrix()
    per_seed, oof = cross_validate(df, FEATURE_COLUMNS, make_logistic_regression, "gene_grouped")
    assert len(per_seed) == 1  # deterministic, single run
    assert per_seed[0]["roc_auc"] > 0.8
    assert oof.shape == (len(df),)


def test_bootstrap_ci_brackets_point_estimate():
    df = synthetic_matrix()
    _, oof = cross_validate(df, FEATURE_COLUMNS, make_logistic_regression, "gene_grouped")
    ci = bootstrap_ci(df["y"], oof, n=200)
    point = metrics_from_scores(df["y"], oof)["roc_auc"]
    assert ci["roc_auc"][0] <= point <= ci["roc_auc"][1]


def test_column_subset_wrapper_drops_features():
    df = synthetic_matrix()
    cols = [c for c in FEATURE_COLUMNS if c not in ("phylop", "sec_struct")]
    model = _with_columns(make_xgboost, cols)(df["y"]).fit(df[cols], df["y"])
    names = list(model.named_steps["preprocess"].get_feature_names_out())
    assert not any("phylop" in n or "sec_struct" in n for n in names)
