"""
End-to-end checks of the API fields added in phase 7. They need the
generated pipeline outputs (gitignored), so they skip on a fresh checkout.
"""
import pytest

from src.config import OUTPUTS_DIR, RAW_DIR, REPORTS_DIR

REQUIRED = [
    RAW_DIR / "foldx_ddg.csv",
    RAW_DIR / "alphafold_manifest.csv",
    REPORTS_DIR / "vus_report.csv",
    OUTPUTS_DIR / "metrics.json",
    OUTPUTS_DIR / "model_xgboost.joblib",
]
pytestmark = pytest.mark.skipif(
    not all(p.exists() for p in REQUIRED), reason="pipeline outputs not generated"
)

VARIANT = "MYH7_R403Q"


@pytest.fixture(scope="module")
def client():
    from fastapi.testclient import TestClient

    from api.main import app

    with TestClient(app) as c:
        yield c


def test_structure_serves_alphafold_model_and_ddg(client):
    body = client.get(f"/variants/{VARIANT}/structure").json()
    assert body["has_structure"] and body["model"]["residue"] == 403
    assert body["ddg"]["value"] is not None and body["ddg_category"] != "not modelled"
    assert client.get(body["files"]["wt"]["url"]).status_code == 200


def test_report_includes_integrated_interpretation(client):
    body = client.get(f"/variants/{VARIANT}/report").json()
    integrated = body["integrated"]["value"]
    assert integrated["conclusion"] and integrated["evidence"]
    assert "FoldX" in body["summaries"]["xgb"]["value"]


def test_vus_are_ranked(client):
    items = client.get("/vus", params={"page_size": 5}).json()["items"]
    ranks = [item["ranking"]["rank"] for item in items]
    assert ranks == sorted(ranks) and ranks[0] == 1


def test_metrics_match_saved_held_out_scores(client):
    body = client.get("/metrics").json()
    assert body["cross_validation"]["gene_grouped"]
    assert all(flags["all_matched"] for flags in body["comparison"].values())
