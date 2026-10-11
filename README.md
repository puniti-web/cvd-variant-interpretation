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
| 2 | Feature set v2 | done |
| 3 | Modeling and evaluation | done |
| 4 | Explainability | in progress |
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
| Dropped: start-loss (p.Met1X) | 66 |
| Unique protein changes after collapsing duplicates | 22,162 |
| **Labelled training set (≥1 review star, no label conflict)** | **1,886** |
| &nbsp;&nbsp;pathogenic / likely pathogenic | 1,598 |
| &nbsp;&nbsp;benign / likely benign | 288 |
| VUS held out for prioritisation | 16,860 |
| Conflicting classifications (excluded) | 2,315 |

Every substitution is checked against the UniProt canonical sequence
(`sequence[position - 1] == wt_aa`). The position comes from the ClinVar
title, which uses the MANE Select transcript; other isoform numberings are
used only when they match the canonical sequence. Per-gene counts are in
`data/processed/qc_per_gene.csv` and the full filter log in
`data/processed/qc_report.json`.

The benign class is small (about 1 benign per 5.5 pathogenic). Benign
missense variants in these genes are rarely submitted to ClinVar.

## Features

Built by `src/features/build_feature_matrix.py` into
`data/features/feature_matrix.csv` (labelled) and `feature_matrix_vus.csv`.

| Group | Features | Source |
| --- | --- | --- |
| Substitution | Δ hydrophobicity, Δ molecular weight, Δ volume, Δ charge, polarity change, Pro/Gly gain or loss, BLOSUM62, Grantham distance | amino-acid property tables |
| Position | relative position in the protein | UniProt sequence length |
| Ortholog conservation | fraction of 60 vertebrate orthologs keeping the human residue; fraction already carrying the mutant residue | UniProtKB orthologs, BLOSUM62 global alignment |
| Genomic conservation | phyloP 100-way, phastCons 100-way, GERP++ | Ensembl VEP REST (dbNSFP) |
| Annotation | in structured domain, in transmembrane segment, at functional site, distance to nearest site, in disordered region, in coiled coil | UniProt |
| Structure | pLDDT, relative solvent accessibility, Cα contact count, secondary structure | AlphaFold DB model |

Raw residue position is no longer a feature: residue 400 means nothing
comparable across proteins. AlphaFold DB has no single full-length model
for RYR2 (4,967 aa), so its structural features are imputed and flagged by
`has_structure`.

SIFT, PolyPhen-2, AlphaMissense and REVEL scores from VEP are stored next
to the features as baselines only. REVEL and AlphaMissense were trained or
calibrated on ClinVar, so using them as inputs would leak the label.

## Model performance

`python -m src.modeling.train_classifier` writes everything below to
`outputs/metrics.json`, out-of-fold predictions to
`outputs/oof_predictions.csv` and curves to `outputs/roc_pr_curves.png`.
1886 labelled variants, 1598 pathogenic / 288 benign. Threshold 0.5.

**Leave-one-gene-out** (headline): each gene is scored by a model that never
saw a variant from it.

| Model | ROC-AUC [95% CI] | PR-AUC | MCC | Balanced acc. | Recall | Specificity |
| --- | --- | --- | --- | --- | --- | --- |
| Logistic Regression | 0.910 [0.884–0.932] | 0.976 | 0.624 | 0.849 | 0.903 | 0.795 |
| Random Forest | 0.885 [0.854–0.915] | 0.960 | 0.654 | 0.797 | 0.969 | 0.625 |
| XGBoost | 0.890 [0.861–0.918] | 0.967 | 0.614 | 0.836 | 0.911 | 0.760 |

**Within-gene** (5-fold, 5 repeats, genes shared between train and test),
for comparison:

| Model | ROC-AUC [95% CI] | PR-AUC | MCC | Balanced acc. | Recall | Specificity |
| --- | --- | --- | --- | --- | --- | --- |
| Logistic Regression | 0.926 [0.906–0.946] | 0.982 | 0.649 | 0.863 | 0.908 | 0.817 |
| Random Forest | 0.929 [0.909–0.951] | 0.982 | 0.737 | 0.853 | 0.970 | 0.735 |
| XGBoost | 0.931 [0.913–0.954] | 0.983 | 0.687 | 0.873 | 0.926 | 0.819 |

The small gap between the two tables means the models are not leaning on
gene identity. Logistic regression generalises best to unseen genes;
XGBoost is within about 0.02 ROC-AUC and is the model explained with SHAP.

**Published predictors** on the same variants (best leave-one-gene-out
model alongside):

| Predictor | Coverage | ROC-AUC | Our model, same variants |
| --- | --- | --- | --- |
| alphamissense | 93% | 0.970 | 0.915 |
| revel | 99% | 0.958 | 0.918 |
| polyphen_score | 82% | 0.891 | 0.929 |
| sift_score | 84% | 0.759 | 0.921 |

The model beats SIFT and PolyPhen-2 but not REVEL or AlphaMissense. Both of
those used ClinVar labels in training or calibration, so their scores on
ClinVar variants are optimistic and not a fair comparison.

**Feature-group ablation** (leave-one-gene-out ROC-AUC):

| Removed | Logistic regression | XGBoost |
| --- | --- | --- |
| all features | 0.910 | 0.890 |
| without conservation | 0.803 | 0.800 |
| without structure | 0.888 | 0.861 |
| without annotation | 0.914 | 0.898 |
| without panel | 0.904 | 0.899 |

Conservation carries most of the signal, with structure second. UniProt
annotation features add little once conservation and structure are present.
`relative_position` was dropped after this kind of comparison showed it
lowered leave-one-gene-out ROC-AUC for all three models. Choosing features
on the same CV adds slight optimism.

Saved models (`outputs/`):

- `model_*.joblib`: trained on the training genes of `gene_level_split`
  only, so the API can score them on genes they never saw.
- `model_*_full.joblib`: refit on all labelled variants, used for VUS
  scoring and explanations.
- `model_xgboost_full_calibrated.joblib`: sigmoid-calibrated XGBoost
  (gene-grouped folds) for report confidence values.

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

or one stage at a time: `collect`, `preprocess`, `annotate`, `features`,
`train`, `explain`. `annotate` calls Ensembl VEP for every variant and takes
about an hour for the full VUS set; it resumes from
`data/raw/vep_annotations.csv` if interrupted. Use
`python -m src.data_collection.fetch_ensembl_vep --labelled-only` for the
training set alone (a few minutes).

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
