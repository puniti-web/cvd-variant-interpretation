# Explainable AI for Cardiovascular Variant Interpretation

A research prototype that classifies missense variants in inherited cardiac
disease genes as pathogenic or benign, explains each prediction with SHAP, and
(in later phases) adds protein-structure and stability evidence.

Workflow:

```text
ClinVar / UniProt / Ensembl VEP / AlphaFold DB
        │
        ▼
Variant dataset (validated against UniProt canonical sequence)
        │
        ▼
Features: physicochemical + substitution matrices + conservation
          + domain annotation + AlphaFold structural context
        │
        ▼
Classifier (Logistic Regression / Random Forest / XGBoost), gene-grouped CV
        │
        ▼
SHAP explanations (global + per variant)
        │
        ▼
Structure + stability (FoldX ΔΔG)              ← phase 5, not started
        │
        ▼
Integrated variant report                      ← phase 6, not started
```

This is not a clinical tool. Predictions must not be used for diagnosis.

## Status

| Phase | Scope | Status |
| --- | --- | --- |
| 0 | Repo hygiene, dependencies, config | done |
| 1 | Data foundation (complete ClinVar pull, isoform-safe parsing, QC) | in progress |
| 2 | Feature set v2 | planned |
| 3 | Modeling and evaluation | planned |
| 4 | Explainability | planned |
| 5 | Structure and stability (FoldX) | to be discussed |
| 6 | Integrated report | planned |
| 7 | Frontend integration | planned |

## Setup

```bash
python -m venv .venv
.venv\Scripts\activate          # Windows
pip install -r requirements.txt
```

NCBI asks for a contact email on E-utilities requests. Set it through the
environment rather than editing code:

```bash
set NCBI_EMAIL=you@example.com      # Windows cmd
$env:NCBI_EMAIL="you@example.com"   # PowerShell
```

`NCBI_API_KEY` is optional and raises the rate limit.

## Running the pipeline

```bash
python -m src.pipeline --stage all
```

or one stage at a time: `collect`, `preprocess`, `features`, `train`, `explain`.

## Tests

```bash
pytest -q
```

## Frontends

- `web/` + `api/`: React app with a read-only FastAPI layer. See
  `web/README_FRONTEND.md` and `api/API_CONTRACT.md`. Start both with `dev.ps1`.
- `app.py`: earlier Streamlit prototype (`streamlit run app.py`).

## Project layout

```text
config/genes.yaml           gene panel (18 cardiomyopathy / arrhythmia genes)
src/data_collection/        ClinVar, UniProt, structure retrieval
src/preprocessing/          parsing, labelling, QC
src/features/               feature engineering
src/modeling/               training and evaluation
src/explainability/         SHAP
src/structure/              structure and stability (phase 5)
src/report/                 integrated report (phase 6)
data/processed/             labelled and VUS variant tables (tracked)
data/raw, data/features     generated, gitignored
outputs/                    models, plots, metrics (gitignored)
```

## Known limitations

- ClinVar labels are noisy, and benign missense variants are scarce in these
  genes, so the benign class is small.
- Evaluation holds out whole genes. Scores are lower than within-gene splits
  but closer to performance on a gene the model has not seen.
- Earlier structural proof of concept (ACTC1 R39H) renamed the residue without
  rebuilding side chains, so its 0.00 Å RMSD carries no biological meaning.
  Phase 5 replaces it.
