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
| 1 | Data foundation (complete ClinVar pull, isoform-safe parsing, QC) | done |
| 2 | Feature set v2 | in progress |
| 3 | Modeling and evaluation | planned |
| 4 | Explainability | planned |
| 5 | Structure and stability (FoldX) | to be discussed |
| 6 | Integrated report | planned |
| 7 | Frontend integration | planned |

## Data

ClinVar missense SNVs for the 18-gene panel, pulled in full through NCBI
E-utilities (`src/data_collection/fetch_clinvar.py`).

| Step | Count |
| --- | ---: |
| Raw ClinVar records | 23,039 |
| Dropped: not missense on the UniProt canonical sequence | 252 |
| Dropped: other gene / unparseable / synonymous | 15 |
| Unique protein changes after collapsing duplicates | 22,216 |
| **Labelled training set (≥1 review star, no label conflict)** | **1,902** |
| &nbsp;&nbsp;pathogenic / likely pathogenic | 1,611 |
| &nbsp;&nbsp;benign / likely benign | 291 |
| VUS held out for prioritisation | 16,883 |
| Conflicting classifications (excluded) | 2,322 |

Every substitution is checked against the UniProt canonical sequence
(`sequence[position - 1] == wt_aa`). The position comes from the ClinVar
title, which uses the MANE Select transcript; other isoform numberings are
used only when they match the canonical sequence. Per-gene counts are in
`data/processed/qc_per_gene.csv` and the full filter log in
`data/processed/qc_report.json`.

The benign class is small (about 1 benign per 5.5 pathogenic). Benign
missense variants in these genes are rarely submitted to ClinVar.

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
