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
Structure + stability (FoldX ΔΔG, run on AWS EC2)
        │
        ▼
Integrated variant report (rule-based conclusion, ranked VUS)
```

This is not a clinical tool. Predictions must not be used for diagnosis.

## Status

| Phase | Scope | Status |
| --- | --- | --- |
| 0 | Repo hygiene, dependencies, config | done |
| 1 | Data foundation (complete ClinVar pull, isoform-safe parsing, QC) | done |
| 2 | Feature set v2 | done |
| 3 | Modeling and evaluation | done |
| 4 | Explainability | done |
| 5 | Structure and stability (FoldX ΔΔG) | done |
| 6 | Integrated report | done |
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
| Stability | FoldX ΔΔG (mutant − wild type, kcal/mol) | FoldX 5.1 on the AlphaFold model, see below |

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
| Logistic Regression | 0.910 [0.884–0.932] | 0.976 | 0.621 | 0.848 | 0.902 | 0.795 |
| Random Forest | 0.889 [0.859–0.918] | 0.962 | 0.636 | 0.787 | 0.968 | 0.608 |
| XGBoost | 0.893 [0.865–0.919] | 0.969 | 0.642 | 0.858 | 0.907 | 0.809 |

**Within-gene** (5-fold, 5 repeats, genes shared between train and test),
for comparison:

| Model | ROC-AUC [95% CI] | PR-AUC | MCC | Balanced acc. | Recall | Specificity |
| --- | --- | --- | --- | --- | --- | --- |
| Logistic Regression | 0.925 [0.906–0.946] | 0.982 | 0.647 | 0.862 | 0.907 | 0.816 |
| Random Forest | 0.929 [0.908–0.952] | 0.982 | 0.734 | 0.852 | 0.969 | 0.734 |
| XGBoost | 0.932 [0.914–0.954] | 0.983 | 0.689 | 0.873 | 0.927 | 0.819 |

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
| all features | 0.910 | 0.893 |
| without conservation | 0.810 | 0.806 |
| without structure (incl. ΔΔG) | 0.888 | 0.861 |
| without FoldX ΔΔG only | 0.910 | 0.890 |
| without annotation | 0.915 | 0.902 |
| without panel | 0.904 | 0.903 |

Conservation carries most of the signal, with structure second. FoldX ΔΔG
adds little on top of the other structural features (see Stability). UniProt
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

## Explainability

`python -m src.explainability.shap_explain` explains the full XGBoost model
with SHAP's TreeExplainer (exact for trees). Values are in log-odds:
positive pushes towards pathogenic. One-hot columns are summed back to
their source feature.

Global importance (mean |SHAP|, `outputs/shap_global_bar.png` and
`shap_global_summary.png`):

| Rank | Feature | mean \|SHAP\| |
| --- | --- | --- |
| 1 | Residue kept in vertebrate orthologs | 0.95 |
| 2 | Genomic conservation (phyloP) | 0.72 |
| 3 | AlphaFold confidence (pLDDT) | 0.45 |
| 4 | Genomic conservation (GERP++) | 0.44 |
| 5 | FoldX stability change (ΔΔG) | 0.40 |
| 6 | Relative solvent accessibility | 0.35 |
| 7 | Gene panel | 0.29 |
| 8 | Change in molecular weight | 0.24 |

**Direction check.** For features whose effect is known from biology, the
Spearman correlation between feature value and SHAP value should have the
expected sign (`outputs/shap_direction_check.json`):

| Feature | Expected | ρ | Matches |
| --- | --- | --- | --- |
| phylop | towards pathogenic | +0.86 | yes |
| phastcons | towards pathogenic | -0.36 | no |
| gerp | towards pathogenic | -0.59 | no |
| conservation_score | towards pathogenic | +0.88 | yes |
| mut_aa_in_orthologs | towards benign | -0.52 | yes |
| blosum62 | towards benign | -0.90 | yes |
| grantham | towards pathogenic | -0.17 | no |
| rsa | towards benign | -0.85 | yes |
| plddt | towards pathogenic | +0.87 | yes |
| contact_count | towards pathogenic | -0.78 | no |
| foldx_ddg | towards pathogenic | +0.54 | yes |
| in_disordered_region | towards benign | -0.51 | yes |
| plddt | towards pathogenic | +0.90 | yes |
| contact_count | towards pathogenic | -0.72 | no |
| in_disordered_region | towards benign | -0.52 | yes |

The three mismatches are each strongly correlated with a stronger feature:
GERP with phyloP, contact count with solvent accessibility (ρ = −0.85),
Grantham with BLOSUM62 (ρ = −0.82). Their raw pathogenic rates are also not
monotonic. The trees use them as residual corrections, so their individual
attributions should not be read on their own. An XGBoost with monotonic
constraints forcing the expected directions was tested: it lost about 0.01
ROC-AUC and 0.035 MCC (leave-one-gene-out), so the unconstrained model was
kept. "Inside a structured domain" also leans slightly benign once
conservation and structure are in the model.

**Per-variant explanations:**

```bash
python -m src.explainability.shap_explain --variant MYH7:R403Q
```

returns JSON with the calibrated probability, the model log-odds and the
top five contributions, for example:

```text
MYH7 R403Q: pathogenic (p=0.98, ClinVar pathogenic)
  Genomic conservation (phyloP) = 7.77 (+1.05 towards pathogenic)
  Residue kept in vertebrate orthologs = 100% (+1.00 towards pathogenic)
  Genomic conservation (GERP++) = 4.04 (+0.92 towards pathogenic)
```

Variants not in either matrix are featurised on the fly, but genomic
conservation needs a VEP lookup and is imputed for them.

**VUS prioritisation.** All 16,860 VUS are scored and ranked in
`outputs/vus_predictions.csv` with their top three reasons. These are
research hypotheses, not reclassifications.

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

## Stability (FoldX ΔΔG)

`src/structure/foldx_ddg.py` computes the FoldX 5.1 stability change for
every variant at a confidently modelled residue:

- Only residues with AlphaFold pLDDT ≥ 70 are modelled (8,984 of 18,746
  variants, 1,268 of them labelled). Others get no ΔΔG and the classifier
  imputes it. RYR2 has no single AlphaFold model.
- Residues with pLDDT < 50 are removed. Proteins over 700 residues are cut
  into 300-residue cores plus every residue within 12 Å (Cα) of the core,
  so each variant keeps its full local environment.
- `RepairPDB` once per chunk (39 chunks), then `BuildModel` with one run per
  variant. On a 10-variant test, three runs changed ΔΔG by at most
  0.07 kcal/mol.

The run takes about 15 core-hours. It ran on one `m7a.8xlarge` Spot instance
(32 cores) in 30 minutes, using `scripts/aws/foldx_instance.sh` as user data:
the instance pulls FoldX and the inputs from a private S3 bucket, syncs
results back every 5 minutes (a restarted run resumes) and terminates
itself. All 8,984 variants returned a value.

```bash
python -m src.structure.foldx_ddg prepare     # chunks + mutation batches
python -m src.structure.foldx_ddg package     # data/raw/foldx_work.tar.gz for upload
python3 foldx_ddg.py run --foldx <binary>     # on the compute machine
python -m src.structure.foldx_ddg collect     # -> data/raw/foldx_ddg.csv
```

FoldX needs an academic licence; the binary is not in this repository.

Result: median ΔΔG is 1.13 kcal/mol for pathogenic and 0.41 for benign
variants. On its own ΔΔG reaches ROC-AUC 0.68 on the labelled variants it
covers. In the model it is the fifth most important feature by SHAP and
points the expected way (higher ΔΔG towards pathogenic), but leave-one-gene-out
ROC-AUC moves by only +0.003 (XGBoost, RF) and 0.000 (logistic regression),
well inside the confidence intervals. Most of what ΔΔG captures is already
present in solvent accessibility, packing and conservation.

## Integrated report

`python -m src.report.generate_report` (or `--stage report`) combines every
earlier stage into one report per variant (`outputs/reports/`):

- ClinVar record and disease from `config/genes.yaml`
- model probability: for labelled variants the out-of-fold XGBoost score
  with the variant's gene held out; for VUS the calibrated XGBoost
- top SHAP contributions
- FoldX ΔΔG and its category, AlphaFold context, UniProt domain
- conservation, and AlphaMissense / REVEL as reference only

The conclusion is assembled from fixed rules in `src/report/generate_report.py`,
so every sentence traces to a threshold:

| Element | Rule |
| --- | --- |
| Tier | P ≥ 0.9 strong pathogenic, ≥ 0.7 moderate, 0.3–0.7 inconclusive, ≤ 0.3 moderate benign, < 0.1 strong benign |
| Probability scale | Calibrated P re-expressed for a 50/50 prior. The calibrated value assumes the training mix (85% pathogenic) and would call most VUS pathogenic |
| Conservation | pathogenic if ≥ 90% of orthologs keep the residue and phyloP ≥ 2; benign if the mutant residue occurs in ≥ 5% of orthologs or < 60% keep the residue |
| Stability | ΔΔG ≥ 1 kcal/mol destabilising (≥ 2 strongly). Neutral or stabilising ΔΔG is *not* counted as benign evidence |
| Structure | pathogenic if buried (RSA < 25%) at pLDDT ≥ 70; benign if pLDDT < 50 or in a disordered region |
| Reference predictors | AlphaMissense 0.34 / 0.564, REVEL 0.290 / 0.644 (ClinGen). Reported, not counted, since both used ClinVar |
| Flags | disagreement with an evidence line, a reference predictor or ClinVar; fewer than 10 benign training variants in the gene; no ΔΔG |

Example (`MYH7_R403Q.md`): *Moderate model support for pathogenic (P = 0.74).
Independent evidence agreeing: 2 of 2 (conservation, stability).* The
out-of-fold P is lower than the full model's 0.97 because MYH7 was held out.

**Checks.** On labelled variants (out-of-fold), the share that ClinVar calls
pathogenic is 97% in the strong-pathogenic tier and 12% in the strong-benign
tier. On VUS, 86% of strong-pathogenic calls are AlphaMissense-pathogenic and
89% of strong-benign calls are AlphaMissense-benign.

**VUS.** All 16,860 VUS are scored in `outputs/reports/vus_report.csv`
(tier, ΔΔG, evidence agreement, top reasons) and the top 50 get full reports
in `vus_top50.md`. Tier counts: 1,654 strong pathogenic, 5,850 moderate
pathogenic, 2,857 inconclusive, 2,342 moderate benign, 4,157 strong benign.
Probabilities saturate near 0.91 at the top, so ranks within the top tier
are effectively ties.

Tiers are model-support levels, not ACMG/AMP classes.

The API's report endpoint now passes ΔΔG into its summary, and its
per-variant SHAP works again (`TreeExplainer` with `model_output="probability"`
needs background data; it now uses 200 training-gene rows).

## Project layout

```text
config/genes.yaml           gene panel (18 cardiomyopathy / arrhythmia genes)
src/data_collection/        ClinVar, UniProt, structure retrieval
src/preprocessing/          parsing, labelling, QC
src/features/               feature engineering
src/modeling/               training and evaluation
src/explainability/         SHAP
src/structure/              FoldX ΔΔG pipeline (phase 5)
scripts/aws/                EC2 user data for the FoldX run
src/report/                 integrated variant report (phase 6)
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
  FoldX ΔΔG (phase 5) replaces it as the stability measure; the API still
  shows the old proof of concept until phase 7.
- FoldX ΔΔG is computed on predicted structures, single chains, and chunks
  of large proteins. It ignores partner proteins and membranes, so variants
  at interfaces (e.g. myosin-actin, channel tetramers) are underestimated.
