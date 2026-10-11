# Phase 0 audit

Read-only inspection of this checkout on 2026-10-11. No existing file was modified. Nothing was retrained. SHAP was not run.

The API must import `src.*` from the repo root and read artifacts. It must not reimplement feature engineering, conservation, SHAP, or RMSD, and it must not train.

## What is actually on disk

| Path | Present | Notes |
| --- | --- | --- |
| `config/genes.yaml` | yes | 18 genes. Source for the panel list. |
| `data/processed/variants_labelled.csv` | yes | 570 × 14. Not gitignored. |
| `data/processed/variants_vus.csv` | yes | 6,122 × 14. Not gitignored. |
| `data/raw/clinvar_raw.csv` | no | Gitignored. README's 7,688 cannot be recounted. |
| `data/raw/uniprot_entries.json` | no | Gitignored. |
| `data/raw/structure_manifest.csv` | no | Gitignored. |
| `data/raw/conservation_scores.json` | no | Gitignored. |
| `data/features/feature_matrix.csv` | no | Gitignored. |
| `outputs/model_random_forest.joblib` | no | Directory absent. |
| `outputs/model_xgboost.joblib` | no | Directory absent. |
| `outputs/shap_global_summary.png` | no | Only SHAP artifact the script writes. |
| `outputs/ACTC1_R39H_wt.pdb` | no | Input to the mutator, not an output of it. |
| `outputs/ACTC1_R39H_mutant.pdb` | no | |
| `models/` | no | Gitignored, and unused. Models are written under `outputs/`. |
| Ablation results | no | No script and no results file anywhere in the repo. |

`data/processed/` is the only generated table available. Counts below that come from those two CSVs were measured. Counts that need a missing file are not filled in from the README.

The project `venv/` exists and has none of the pipeline dependencies installed. This audit's split reproduction used the system interpreter's scikit-learn **1.5.2**. `xgboost` is not importable there. The sklearn version that originally trained the (absent) joblibs is unknown. A future mismatch against the published metrics should be reported by the API, not patched.

## Measured labelled cohort

`variants_labelled.csv` and `variants_vus.csv` share this header:

`clinvar_uid`, `title`, `gene_sort`, `clinical_significance`, `review_status`, `variation_set`, `protein_change`, `molecular_consequence_list`, `gene_symbol`, `panel`, `wt_aa`, `position`, `mut_aa`, `label`

Labelled: **570** rows, **0** duplicate ids on `(gene_symbol, wt_aa, position, mut_aa)`, **0** nulls. **424** pathogenic, **146** benign. Panel: cardiomyopathy **388**, arrhythmia **182**. All **18** panel genes appear. Amino acids are the 20 standard one-letter codes (no stops).

Raw `clinical_significance` values that `normalize_label` collapsed:

| Raw ClinVar string | Rows | Normalized `label` |
| --- | ---: | --- |
| Likely pathogenic | 280 | pathogenic |
| Pathogenic | 88 | pathogenic |
| Pathogenic/Likely pathogenic | 56 | pathogenic |
| Likely benign | 133 | benign |
| Benign | 8 | benign |
| Benign/Likely benign | 5 | benign |

VUS file: **6,122** rows, label `vus` on every row, raw significance `Uncertain significance`, same 18 genes, **0** duplicate ids. One null `protein_change` (the parsed `wt_aa` / `position` / `mut_aa` columns are still filled).

Rows with any other ClinVar significance (conflicting, risk factor, and so on) are written to neither file. `570 + 6,122 = 6,692`. The README's 7,688 raw records are **not** in this checkout, so the **996**-record gap cannot be verified here. The API must count `clinvar_raw.csv` when it exists and return null until then.

Per-gene labelled counts (measured):

| Gene | Panel | Benign | Pathogenic | Total |
| --- | --- | ---: | ---: | ---: |
| ACTC1 | cardiomyopathy | 1 | 20 | 21 |
| MYBPC3 | cardiomyopathy | 8 | 3 | 11 |
| MYH7 | cardiomyopathy | 5 | 38 | 43 |
| MYL2 | cardiomyopathy | 1 | 12 | 13 |
| MYL3 | cardiomyopathy | 0 | 3 | 3 |
| TNNI3 | cardiomyopathy | 2 | 56 | 58 |
| TNNT2 | cardiomyopathy | 3 | 48 | 51 |
| TPM1 | cardiomyopathy | 27 | 44 | 71 |
| LMNA | cardiomyopathy | 25 | 63 | 88 |
| PRKAG2 | cardiomyopathy | 17 | 12 | 29 |
| CACNA1C | arrhythmia | 17 | 20 | 37 |
| CASQ2 | arrhythmia | 7 | 9 | 16 |
| KCNE1 | arrhythmia | 1 | 0 | 1 |
| KCNE2 | arrhythmia | 3 | 1 | 4 |
| KCNH2 | arrhythmia | 9 | 19 | 28 |
| KCNQ1 | arrhythmia | 9 | 48 | 57 |
| RYR2 | arrhythmia | 5 | 15 | 20 |
| SCN5A | arrhythmia | 6 | 13 | 19 |

VUS counts by gene (measured): ACTC1 308, CACNA1C 451, CASQ2 272, KCNE1 20, KCNE2 68, KCNH2 456, KCNQ1 415, LMNA 372, MYBPC3 481, MYH7 451, MYL2 237, MYL3 232, PRKAG2 430, RYR2 479, SCN5A 424, TNNI3 315, TNNT2 349, TPM1 362.

### ACTC1 R39H is a labelled benign test-gene variant

One matching row:

| Field | Value |
| --- | --- |
| `clinvar_uid` | 2681695 |
| `title` | `NM_005159.5(ACTC1):c.116G>A (p.Arg39His)` |
| `clinical_significance` | Likely benign |
| `review_status` | no assertion criteria provided |
| `protein_change` | R39H |
| `panel` | cardiomyopathy |
| `label` | **benign** |
| Variant id | `ACTC1_R39H` |

`variation_set` is a Python literal (not JSON), about 1,033 characters. For this row it contains `cdna_change` `c.116G>A`, dbSNP `2504185412`, and GRCh38 chr 15 position 34794693. `molecular_consequence_list` is the literal `['missense variant']`. Parse with `ast.literal_eval` only. Never `eval`.

The structure proof-of-concept is this same variant. Its clinical label is Likely benign, on a weak review status. The model prediction, when artifacts exist, is a separate field. Do not present R39H as a pathogenic case.

A real test-gene variant with no structure wiring, for the later smoke test: **`CACNA1C_G1365R`**, label pathogenic, raw significance `Pathogenic/Likely pathogenic`.

`compute_physchem_features("R", "H")` returns, and these are the only physchem numbers that may be shown for R39H until the feature matrix exists:

- `delta_hydrophobicity` 1.3 (the float is `1.2999999999999998`)
- `delta_molecular_weight` -19.0
- `delta_volume` -20.2 (the float is `-20.200000000000017`)
- `delta_charge` -0.9
- `polarity_change` false (both R and H are in the polar set)
- all four proline/glycine flags false

Conservation, domain membership, probabilities, and SHAP for this variant are **not computable** in this checkout.

## 1. Models: path, format, columns, encoding

Written by `src/modeling/train_classifier.py` `main()`:

- `outputs/model_random_forest.joblib`
- `outputs/model_xgboost.joblib`

`joblib.dump` of a full `sklearn.pipeline.Pipeline`, not a bare estimator and not a file under `models/`. Steps, in order:

1. `preprocess` — `ColumnTransformer`
2. `clf` — `RandomForestClassifier` or `XGBClassifier`

There is no sidecar JSON for columns, medians, categories, the split, or metrics. `evaluate_classifier` prints metrics and returns a dict that `main()` discards.

### Estimator settings

Random Forest (`train_random_forest`): `n_estimators=500`, `max_depth=None`, `class_weight="balanced"`, `random_state=42`, `n_jobs=-1`.

XGBoost (`train_xgboost`): `n_estimators=400`, `max_depth=5`, `learning_rate=0.05`, `eval_metric="logloss"`, `random_state=42`, `n_jobs=-1`, `scale_pos_weight = n_benign / n_pathogenic` computed on the **training** labels only.

SMOTE is not used. `imbalanced-learn` is in `requirements.txt`. The module docstring says a SMOTE path is commented out; there is no such block in the file.

Decision threshold is scikit-learn's default **0.5**. Nothing is calibrated. Positive class is `label == "pathogenic"` (`y = 1`).

### Columns that enter the model, in this order

`feature_cols = NUMERIC_FEATURES + BOOLEAN_FEATURES + CATEGORICAL_FEATURES` in `train_classifier.py`:

| # | Column | Transformer |
| --- | --- | --- |
| 1 | `position` | median imputer |
| 2 | `delta_hydrophobicity` | median imputer |
| 3 | `delta_molecular_weight` | median imputer |
| 4 | `delta_volume` | median imputer |
| 5 | `delta_charge` | median imputer |
| 6 | `conservation_score` | median imputer |
| 7 | `polarity_change` | median imputer (treated as numeric) |
| 8 | `wt_is_proline` | median imputer |
| 9 | `mut_is_proline` | median imputer |
| 10 | `wt_is_glycine` | median imputer |
| 11 | `mut_is_glycine` | median imputer |
| 12 | `in_annotated_domain` | median imputer |
| 13 | `panel` | `OneHotEncoder(handle_unknown="ignore")` |

The imputer is `SimpleImputer(strategy="median")` and is **fit on the training fold only**, inside the pipeline. A missing conservation score is replaced by the training-fold median, not a median computed in the API.

The comment under `CATEGORICAL_FEATURES` (lines 50–52 of `train_classifier.py`) says `conservation_score` is left out of `NUMERIC_FEATURES`. That comment is stale. The column **is** in `NUMERIC_FEATURES`. `build_feature_matrix.py`'s docstring is also stale: conservation is no longer a placeholder returning `None` for every row. The live path is `conservation_and_domain.py`.

### Columns in the 18-wide matrix that do **not** enter the model

Built by `build_feature_row` in this order:

`gene_symbol`, `panel`, `position`, `wt_aa`, `mut_aa`, `label`, `delta_hydrophobicity`, `delta_molecular_weight`, `delta_volume`, `delta_charge`, `polarity_change`, `wt_is_proline`, `mut_is_proline`, `wt_is_glycine`, `mut_is_glycine`, `conservation_score`, `in_annotated_domain`, `domain_type`

Not model inputs: `gene_symbol`, `wt_aa`, `mut_aa`, `label`, `domain_type`.

`domain_type` is stored and displayed, then dropped. Only the boolean `in_annotated_domain` is trained on. Domain type is the **first** UniProt feature whose span covers the position, among `Domain`, `Region`, `Active site`, `Binding site`, `Disulfide bond`, `Transmembrane`, `Topological domain`.

### Encoded feature names

Confirmed with scikit-learn 1.5.2 (`get_feature_names_out`), which matches `requirements.txt` (`scikit-learn>=1.3`, dense `OneHotEncoder` output). A saved pipeline from an older sklearn could differ; the API should read names from the loaded preprocessor rather than hardcode them.

```text
num__position
num__delta_hydrophobicity
num__delta_molecular_weight
num__delta_volume
num__delta_charge
num__conservation_score
num__polarity_change
num__wt_is_proline
num__mut_is_proline
num__wt_is_glycine
num__mut_is_glycine
num__in_annotated_domain
cat__panel_arrhythmia
cat__panel_cardiomyopathy
```

14 columns after transform. Panel categories sort alphabetically. Both panels are guaranteed to appear in training if all 18 genes are present: the test fold is 5 genes, and the smaller panel (arrhythmia) has 8 genes, so neither panel can be entirely held out.

`ColumnTransformer` was fit on a DataFrame, so prediction must pass a DataFrame with those 13 column names. Booleans should be real bools, `conservation_score` a float or NaN, `panel` a string.

Training `dropna` keeps rows with a null `conservation_score` and drops rows missing any other model feature or `label`. Physchem and the boolean flags are total functions of the two amino acids, and `in_annotated_domain` is always a bool, so on the current labelled table that dropna removes nothing. Unverified until `feature_matrix.csv` exists.

## 2. Leakage check

`gene_symbol` is the split key only. It is **not** a model feature. Putting it in the matrix does not leak gene identity into the trees. Under this split an unseen gene would have been all zeros anyway (`handle_unknown="ignore"`), which is why excluding it is the correct choice. Do not add it in the API.

`panel` **is** a model feature. It is constant for every variant of a gene (`genes.yaml` assigns each gene to one panel). It is not a property of the substitution. With a gene-level split it still transfers the cardiomyopathy-vs-arrhythmia base rate onto unseen genes. Both panels occur on both sides of the reproduced split, so the one-hot is not a dead column on the test set. **Flag `panel` as a gene-level covariate with no variant-level signal.** The published metrics include it.

`position` **is** a model feature. Residue numbers are not comparable across proteins (ACTC1 39 vs MYH7 403 vs a CACNA1C site in the thousands; labelled positions run from 1 to 4880). The README says the global SHAP plot is driven substantially by `position`, `delta_volume`, `delta_hydrophobicity`, `delta_molecular_weight`, `delta_charge`, and `conservation_score`. That plot is not in this checkout, so the ranking is **not** remeasured here. A large SHAP share for `position` should be treated as a possible dataset artifact, not as evidence that absolute residue index is biological mechanism. This is exactly what the ablation "with vs without `position`" is for. That ablation has not been run.

`wt_aa` and `mut_aa` are not one-hot encoded. Their effect reaches the model only through the physchem deltas and the proline/glycine flags.

## 3. How the 13/5 split is chosen

The docstring says GroupKFold. The code is one `GroupShuffleSplit`:

```python
GroupShuffleSplit(n_splits=1, test_size=0.25, random_state=42)
# groups = gene_symbol
```

`gene_level_split` in `train_classifier.py`. No gene list is saved.

scikit-learn takes the unique group labels (`np.unique`, so alphabetical), then holds out `ceil(0.25 * n_genes)` groups. For 18 genes that is **5 test genes and 13 train genes**. Membership does not depend on row order. It does depend on the set of genes that survive `dropna`, and on the sklearn/numpy permutation used by `ShuffleSplit`.

Reproduced here with scikit-learn 1.5.2 on `variants_labelled.csv` (groups = `gene_symbol`):

**Test genes (98 variants: 62 pathogenic, 36 benign)**

`ACTC1`, `CACNA1C`, `KCNE1`, `KCNH2`, `MYBPC3`

**Train genes (472 variants: 362 pathogenic, 110 benign)**

`CASQ2`, `KCNE2`, `KCNQ1`, `LMNA`, `MYH7`, `MYL2`, `MYL3`, `PRKAG2`, `RYR2`, `SCN5A`, `TNNI3`, `TNNT2`, `TPM1`

This matches the README's 13/472 and 5/98. It is the split the API should recompute by calling `gene_level_split` after the training `dropna`. Do not hardcode the lists. If a future feature matrix drops rows or a different sklearn changes the permutation, the comparison flag against the README metrics is how that shows up.

`ACTC1` is a **test** gene. The structure proof-of-concept is an unseen-gene example, not an in-sample one. The UI chip for R39H is "test gene", not the train-gene warning. The train-gene warning is for the other 13 genes, whose variants were in the training fold (the split is by gene, so every variant of a train gene was eligible for training).

`KCNE1` contributes a single labelled variant, and it is benign, and the gene is in the test fold.

The test-set confusion-matrix margins in the README match these label counts: 36 benign and 62 pathogenic.

`python -m src.pipeline` never runs structure, stability, conservation, or the report. Its stages are only `collect`, `preprocess`, `features`, `train`, `explain`.

## 4. SHAP is not saved per variant

`src/explainability/shap_explain.py`:

- Default model is **`outputs/model_xgboost.joblib` only**. Random Forest SHAP is not produced unless the caller changes `load_model_and_data`.
- The only file written is `outputs/shap_global_summary.png`.
- No per-variant SHAP array, no `expected_value`, no mean-|SHAP| table.

`compute_shap_values` builds `shap.TreeExplainer` on `pipeline.named_steps["clf"]` and explains `pipeline.named_steps["preprocess"].transform(X)`. `explain_single_variant` returns the top-k rows of that explanation. Both are callable. The API should use them (or `TreeExplainer` the same way) at request time and cache per variant. It must not invent SHAP values when the model is missing.

### Divergences the API has to handle, without editing `src/`

1. **Different row filter.** Training `dropna` **keeps** null conservation scores and lets the imputer fill them. `load_model_and_data` **drops** any row with a null in `NUMERIC_FEATURES + BOOLEAN_FEATURES`, which includes `conservation_score`. The PNG is therefore fit on the scored subset (README: 518 scored, 52 missing), and it mixes train-gene and test-gene rows. The API's global and local SHAP should transform through the **saved pipeline** so the 52 missing scores are imputed the way the model saw them, and each local explanation should set `conservation_imputed: true` when the raw score is null. Do not copy the SHAP script's `dropna` if the goal is "what the model did".

2. **Output space.** `TreeExplainer` is constructed with no `model_output`. For XGBoost that default is the raw margin (log-odds), not `predict_proba`. A waterfall in that space will not land on P(pathogenic). Random Forest explanations are in probability (or a two-class array, depending on the SHAP version). `explain_single_variant` assumes `shap_values.values[index]` is 1-D; a `(n_features, n_classes)` row will be ranked wrong. The API should request `model_output="probability"`, and if the array still has a class axis, take the pathogenic class (index 1) before ranking or before calling `explain_single_variant`. Record `shap_space: "probability"` on the payload. Base value is `explainer.expected_value` for that class. This is a call-site choice. It is not a change to `src/explainability/shap_explain.py`. The PNG remains the script's default-space plot.

3. **`feature_value` is post-transform.** `explain_single_variant` zips SHAP against the imputed / one-hot row, under names like `num__conservation_score`. The dossier should also show the raw feature value from the variant row. One-hot `panel` columns do not have a single raw number; report the raw `panel` string beside them.

4. Global mean |SHAP| is not a function in `src/`. Reduce the matrix returned by `compute_shap_values`. Do not read it off the PNG.

Per-variant SHAP is feasible live for a labelled row or a VUS row once the corresponding joblib loads and the 13-column frame can be built. It is not feasible in this checkout.

## 5. Scoring one VUS without changing `src/`

Feasible, by composition. There is no single "score this variant" function.

| Function | Role |
| --- | --- |
| `build_feature_row(row, uniprot_features_by_gene)` | One-row assembler. Needs a Series with `gene_symbol`, `wt_aa`, `mut_aa`, `position`, and ideally `panel` and `label`. |
| `compute_physchem_features(wt_aa, mut_aa)` | Works now, from the amino-acid tables only. |
| `compute_conservation_and_domain_features(gene_symbol, position, uniprot_features, conservation_scores=...)` | Domain + conservation. Pass a preloaded dict. |
| Saved `Pipeline.predict` / `predict_proba` | Score. Label is not an input. |

`build_feature_row` does **not** accept a preloaded conservation dict. Each call re-reads `conservation_scores.json`. For a 6,122-row VUS pass, call `compute_physchem_features` and `compute_conservation_and_domain_features` with the dict already loaded. That is still `src` code.

VUS rows already have `panel` and `label="vus"`. The model never sees `label`.

What comes back null, and why:

| Missing artifact | Effect |
| --- | --- |
| `uniprot_entries.json` | `in_annotated_domain` false, `domain_type` null, because the domain function returns `(False, None)` on an empty feature list. Say so. Do not pretend the residue was checked against UniProt. |
| `conservation_scores.json` | `conservation_score` null. If the model is loaded, the imputer still fills the training median and `conservation_imputed` is true. |
| Either joblib | That model's `label` and `prob_pathogenic` are null. |
| `variants_vus.csv` | The VUS endpoint has nothing to list. |

Do not score VUS by silently appending them to the training frame. They stay out of `/metrics` and out of global SHAP.

List endpoints must not compute SHAP for every row. SHAP stays on the detail endpoint, behind an `lru_cache` keyed by variant id and model name. A process restart is how newly generated artifacts get picked up; the lifespan loader reads files once.

## 6. What `src/report/` exposes

One module, `src/report/generate_report.py`. Not called by the pipeline.

Callable:

- `VariantReport` dataclass: `gene_symbol`, `variant`, `predicted_pathogenicity`, `confidence_score`, `key_features`, `structural_rmsd`, `predicted_ddg`, `stability_interpretation`, `plain_language_summary`
- `make_report(...)` fills `plain_language_summary` via `build_plain_language_summary`
- `save_report(report, path)` writes that dataclass as JSON

The summary sentence is a template: the model "classifies {gene} {variant} as {label} with {confidence:.0%} confidence", names the first key feature, then appends `stability_interpretation`. It has one prediction, not two. It does not know about clinical review status, conservation, physchem, or the RMSD method artifact. It does not include a clinical disclaimer.

`build_plain_language_summary` says the variant "increases" when `confidence_score >= 0.5` and "decreases" otherwise, without saying increases what. Pass the pathogenic probability in `confidence_score` so the sentence matches the dial. Call `make_report` twice, once per model, when that model actually produced a probability. If the model is missing, do not call it with a made-up score.

Running the module as `__main__` writes `outputs/example_variant_report.json` for a **fictional** MYH7 p.Arg403Gln example (`confidence_score=0.91`, `structural_rmsd=1.8`). That file is a smoke fixture. The API must never serve it as a variant report.

`interpret_stability` in `stability_analysis.py` is a different helper. For RMSD 0.00 and ddG `None` it returns "modest structural deviation (RMSD=0.00 Å); ddG not available". That sentence contradicts the README and `app.py`, which correctly call 0.00 Å a method artifact. **Do not pass `interpret_stability`'s RMSD-only string into `make_report` for the identity-swap proof-of-concept.** Pass an interpretation that states the artifact. Leave `predicted_ddg=None`.

## 7. Structure, RMSD, ΔΔG

Only one variant is wired, and it is hardcoded.

`src/structure/mutate_pdb.py` is a script with **module-level side effects**. Importing it parses `outputs/ACTC1_R39H_wt.pdb`, renames residue 39 from ARG to HIS on the first chain where the wild-type residue is R, and writes `outputs/ACTC1_R39H_mutant.pdb`. Coordinates are unchanged. The API must **not** import this module.

Constants in that file and in the README / `app.py`: PDB **9B3R**, chain **A**, residue **39**, WT **R**, mutant **H**. `9B3R` is not read from `structure_manifest.csv`. The manifest, when it exists, has one row per gene (`gene_symbol`, `uniprot_id`, `source`, `structure_id`) and keeps only the **first** RCSB hit, else an AlphaFold DB id. It is not an exhaustive PDB list, and its `structure_id` for ACTC1 is not guaranteed to be 9B3R.

Nothing in `src/` downloads 9B3R or writes `outputs/ACTC1_R39H_wt.pdb`. `fetch_structures` saves `data/raw/structures/{gene}_{first_pdb_id}.pdb`. The WT proof-of-concept file is a manual input.

`compute_rmsd(wt, mut)` in `stability_analysis.py` superimposes Cα atoms of chain **A** only. Call it when both PDB files exist. Do not return the README's 0.00 unless that call returns 0.00. The method note is still required when the value is 0.00, because that is the expected artifact of the residue rename.

`run_foldx_ddg` returns `None` when the `foldx` binary is absent. When the binary **is** on `PATH`, it raises `NotImplementedError`. The API must not call it. ΔΔG is null, reason: FoldX BuildModel is a stub.

`predict_mutant_structure.py` can POST a sequence to the public ESMFold endpoint (`https://api.esmatlas.com/foldSequence/v1/pdb/`). The README says that endpoint timed out. The API must not call it. The proof-of-concept PDBs did not come from ESMFold.

`has_structure` is true only for variant id `ACTC1_R39H` when both output PDB files exist. Other variants get `has_structure: false` plus whatever single manifest row exists for the gene.

## 8. Commands that regenerate the missing artifacts

From the repo root, after `NCBI_CONTACT_EMAIL` in `src/config.py` is set. Do not use `--stage all` if the published conservation-aware metrics are the target: `build_conservation_scores` is **not** a pipeline stage, and `build_feature_matrix` will silently write an all-null `conservation_score` column when the JSON is missing. Those nulls are then median-imputed to a constant.

```text
python -m src.pipeline --stage collect
python -m src.pipeline --stage preprocess
python -m src.features.build_conservation_scores
python -m src.pipeline --stage features
python -m src.pipeline --stage train
python -m src.pipeline --stage explain
```

Equivalent entry points: `fetch_clinvar`, `fetch_uniprot`, `fetch_structures`, `clean_variants`, `build_feature_matrix`, `train_classifier`, `shap_explain`.

`fetch_ensembl_vep.py` is not part of the pipeline and is not required for any endpoint.

Structure proof-of-concept, after placing the experimental PDB yourself:

```text
# save https://files.rcsb.org/download/9B3R.pdb as outputs/ACTC1_R39H_wt.pdb
python -m src.structure.mutate_pdb
python -m src.structure.stability_analysis --wt outputs/ACTC1_R39H_wt.pdb --mut outputs/ACTC1_R39H_mutant.pdb
```

`stability_analysis` prints `interpret_stability` and does not write a results file. There is no ablation command. Do not add one inside `src/`.

## 9. Streamlit behavior the API must not copy

`app.py` is left untouched. These are traps if the API is written by imitating it:

- The dashboard counts pathogenic/benign with `label == 1` and `label == 0`. The feature matrix stores the strings `pathogenic` and `benign`, so those metrics fall through to zero once the matrix exists.
- The ML page displays the README metrics as literals. It does not load the joblibs.
- The structures page displays `0.00 Å` as a literal even when the PDB files are missing. The file-status table does check that the files exist.
- SHAP is an image lookup for `shap_summary.png`, `shap_summary_plot.png`, or `shap_global_summary.png`. Only the last name is what `shap_explain.py` writes.

## 10. Backend changes I did not make

These are real, and they stay in `src/` until you ask:

1. SHAP `dropna` drops the 52 unscored variants; training imputes them. Global explanations and the training frame are not the same row set.
2. `interpret_stability` describes RMSD 0.00 Å as a modest structural deviation.
3. `mutate_pdb.py` runs on import and is hardcoded to one variant.
4. Stale comments: GroupKFold, conservation "left out", SMOTE "commented out", conservation "placeholder".
5. `--stage all` never builds conservation scores, so a one-shot run does not reproduce the conservation-aware models.
6. `app.py` label comparison against `0`/`1`.

No ablation module exists. The `/ablation` endpoint stays empty until a results file appears under `outputs/`. It must not train.
