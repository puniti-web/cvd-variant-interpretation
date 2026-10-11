"""SHAP helpers: name mapping, value formatting and additivity on a small model."""
import math

import numpy as np

from src.explainability.shap_explain import compute_shap, explain_row, format_value, raw_feature_of
from src.modeling.train_classifier import FEATURE_COLUMNS, make_xgboost
from tests.test_modeling import synthetic_matrix


def test_raw_feature_of_groups_one_hot_columns():
    assert raw_feature_of("num__phylop") == "phylop"
    assert raw_feature_of("cat__sec_struct_helix") == "sec_struct"
    assert raw_feature_of("cat__panel_arrhythmia") == "panel"


def test_format_value():
    assert format_value("conservation_score", 0.95) == "95%"
    assert format_value("in_transmembrane", True) == "yes"
    assert format_value("phylop", float("nan")) == "missing (imputed)"


def test_shap_values_add_up_to_model_log_odds():
    df = synthetic_matrix()
    X = df[FEATURE_COLUMNS]
    model = make_xgboost(df["y"]).fit(X, df["y"])
    values, base = compute_shap(model, X)
    assert list(values.columns) == FEATURE_COLUMNS
    proba = model.predict_proba(X)[:, 1]
    log_odds = np.log(proba / (1 - proba))
    np.testing.assert_allclose(base + values.sum(axis=1).to_numpy(), log_odds, atol=1e-3)

    top = explain_row(X.iloc[0], values.iloc[0], base, top_k=3)
    assert len(top) == 3
    assert abs(top[0]["shap"]) >= abs(top[1]["shap"]) >= abs(top[2]["shap"])
    assert all(not math.isnan(c["shap"]) for c in top)
