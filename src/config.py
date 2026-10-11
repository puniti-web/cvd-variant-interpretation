"""
Shared configuration loader for the pipeline.

Every stage imports PROJECT_ROOT / DATA_DIR / load_gene_panel() from here so
paths stay consistent regardless of which script is run or from where.
"""
import os
from pathlib import Path
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIR = PROJECT_ROOT / "config"
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
FEATURES_DIR = DATA_DIR / "features"
OUTPUTS_DIR = PROJECT_ROOT / "outputs"
REPORTS_DIR = OUTPUTS_DIR / "reports"

# NCBI E-utilities etiquette: set a real contact email before running at
# volume, and keep request rates under NCBI's guidelines (<=3/sec without
# an API key, <=10/sec with one). Get a free key at
# https://www.ncbi.nlm.nih.gov/account/settings/
# Both can be supplied through environment variables so they stay out of git.
NCBI_CONTACT_EMAIL = os.environ.get("NCBI_EMAIL", "your_email@example.com")
NCBI_API_KEY = os.environ.get("NCBI_API_KEY")  # optional

# Fixed seed shared by every stochastic stage (splits, models, SHAP sampling).
RANDOM_STATE = 42

USER_AGENT = "cvd-xai-project/0.1 (student research project)"


def load_gene_panel(which: str = "both") -> list[dict]:
    """
    Load the configured gene panel.

    Parameters
    ----------
    which : "cardiomyopathy" | "arrhythmia" | "both"

    Returns
    -------
    list of dicts, each with keys: symbol, condition, uniprot_id, panel
    """
    with open(CONFIG_DIR / "genes.yaml") as f:
        cfg = yaml.safe_load(f)

    genes = []
    if which in ("cardiomyopathy", "both"):
        for g in cfg["cardiomyopathy_genes"]:
            genes.append({**g, "panel": "cardiomyopathy"})
    if which in ("arrhythmia", "both"):
        for g in cfg["arrhythmia_genes"]:
            genes.append({**g, "panel": "arrhythmia"})
    return genes


def ensure_dirs():
    for d in (RAW_DIR, PROCESSED_DIR, FEATURES_DIR, OUTPUTS_DIR, REPORTS_DIR):
        d.mkdir(parents=True, exist_ok=True)
