# API contract

Base URL: `http://127.0.0.1:8000`. Process must be started with the **repo root** as the working directory so `import src...` resolves. CORS allows `http://localhost:5173`.

Phase 0 writes this contract only. The FastAPI app is Phase 1. The frontend mock client must implement this document exactly, including nulls.

## Rules

1. Missing artifacts do not crash the process. Lifespan load records a reason and continues.
2. A null measurement is `{ "value": null, "reason": "<why>" }`. When `value` is present, `reason` is null.
3. Do not substitute README numbers for a missing file. The one exception is `reported` inside `GET /metrics`, which is explicitly the published comparison baseline and is never placed in `recomputed`.
4. Do not import `src.structure.mutate_pdb` (it mutates a PDB at import time).
5. Do not call `run_foldx_ddg` (it raises `NotImplementedError` if FoldX is on `PATH`).
6. Do not call the ESMFold endpoint.
7. Do not read `outputs/example_variant_report.json` as a real report.
8. Do not train, and do not write model or SHAP artifacts.
9. Variant id is `{gene}_{wt}{position}{mut}` using the one-letter codes in the CSV, for example `ACTC1_R39H`. Ids are unique inside each of the two processed tables and do not collide across them (labels are exclusive).
10. Restart the process to pick up artifacts that appear later. Loaders and `lru_cache` do not watch the disk.

Shared error body for a genuine unknown id (indexes loaded, id in neither table):

```json
{ "detail": { "value": null, "reason": "No labelled or VUS row has id CACNA1C_Z1Z" } }
```

HTTP 404. If the indexes themselves failed to load, return HTTP 200 with null sections and reasons, not 404. 404 means "we looked and it is absent".

Query validation failures (page < 1, unknown `label` / `split`) are HTTP 400 with the same `detail` shape.

## Artifact flags

`GET /health` reports these paths.

| Key | Path |
| --- | --- |
| `clinvar_raw` | `data/raw/clinvar_raw.csv` |
| `uniprot_entries` | `data/raw/uniprot_entries.json` |
| `structure_manifest` | `data/raw/structure_manifest.csv` |
| `conservation_scores` | `data/raw/conservation_scores.json` |
| `variants_labelled` | `data/processed/variants_labelled.csv` |
| `variants_vus` | `data/processed/variants_vus.csv` |
| `feature_matrix` | `data/features/feature_matrix.csv` |
| `model_random_forest` | `outputs/model_random_forest.joblib` |
| `model_xgboost` | `outputs/model_xgboost.joblib` |
| `shap_global_plot` | `outputs/shap_global_summary.png` |
| `actc1_wt_pdb` | `outputs/ACTC1_R39H_wt.pdb` |
| `actc1_mutant_pdb` | `outputs/ACTC1_R39H_mutant.pdb` |

A file that exists but fails to load (`joblib` version skew, bad CSV) is `present: true`, `loaded: false`, with `reason` set to the exception message. `present: false` leaves `loaded: false`.

## Pipeline status file

Phase 1 adds `api/pipeline_status.json` (new file, not an `src/` change). `GET /overview` returns it under `pipeline`. If the file is missing, `pipeline` is null and `reason` says so. This file describes the project's implementation status. It is not a claim that artifacts exist on this machine; artifact presence stays on `/health`.

```json
{
  "steps": [
    { "id": 1, "name": "Data collection", "status": "done", "note": "ClinVar, UniProt, PDB manifest." },
    { "id": 2, "name": "Preprocessing", "status": "done", "note": "Labelled vs VUS split. VUS held out." },
    { "id": 3, "name": "Feature engineering", "status": "done", "note": "Physchem, homolog-identity conservation, domain flag." },
    { "id": 4, "name": "Classification", "status": "done", "note": "Random Forest and XGBoost. Gene-level split. Class weighting, no SMOTE." },
    { "id": 5, "name": "SHAP", "status": "done", "note": "Global plot is XGBoost-only. Per-variant SHAP is computed live by the API." },
    { "id": 6, "name": "Structure", "status": "done_poc", "note": "ACTC1 R39H, PDB 9B3R, chain A. Residue-identity swap." },
    {
      "id": 7,
      "name": "Stability",
      "status": "partial",
      "note": "Cα RMSD is implemented. FoldX ΔΔG is a stub.",
      "parts": [
        { "name": "Cα RMSD", "status": "done" },
        { "name": "FoldX ΔΔG", "status": "stub" }
      ]
    },
    { "id": 8, "name": "Integrated report", "status": "partial", "note": "make_report exists and is not wired to the pipeline." }
  ],
  "ui": {
    "react": "done"
  }
}
```

`status` is one of `done`, `done_poc`, `partial`, `stub`, `not_started`. The helix has these 8 steps. `done_poc`, `partial`, and `stub` render dashed / half-lit. `ui` is the status board, not a ninth rung.

## Types

```ts
type Maybe<T> = { value: T | null; reason: string | null };

type Label = "pathogenic" | "benign" | "vus";
type Split = "train_gene" | "test_gene";
type Panel = "cardiomyopathy" | "arrhythmia";

interface Artifact {
  key: string;
  path: string;
  present: boolean;
  loaded: boolean;
  reason: string | null;
}

interface Prediction {
  label: "pathogenic" | "benign" | null;
  prob_pathogenic: number | null;
  threshold: 0.5;
  reason: string | null;
}

interface ShapFeature {
  feature: string;            // transformed name, e.g. num__delta_charge
  raw_feature: string;        // column name before the transformer
  shap_value: number;         // probability space, pathogenic class
  feature_value_transformed: number;
  feature_value_raw: number | string | boolean | null;
}

interface LocalShap {
  shap_space: "probability";
  base_value: number;
  output_value: number;       // base_value + sum(shap_value); should match prob_pathogenic within 1e-3
  features: ShapFeature[];    // all transformed features, sorted by |shap_value| descending
  reason: string | null;
}
```

`Prediction.reason` is null when `prob_pathogenic` is a number. `LocalShap.reason` is null when `features` is non-empty. A missing model uses `features: []` and a reason, with `base_value` and `output_value` null.

## `GET /health`

```json
{
  "status": "ok",
  "artifacts": [Artifact]
}
```

`status` is `ok` even when every artifact is missing. The process is up. Clients read `artifacts`.

## `GET /overview`

```json
{
  "counts": {
    "collected": Maybe<number>,
    "labelled": Maybe<number>,
    "pathogenic": Maybe<number>,
    "benign": Maybe<number>,
    "vus": Maybe<number>
  },
  "genes": [
    {
      "symbol": "ACTC1",
      "panel": "cardiomyopathy",
      "condition": "Hypertrophic / Dilated Cardiomyopathy",
      "uniprot_id": "P68032",
      "benign": Maybe<number>,
      "pathogenic": Maybe<number>,
      "vus": Maybe<number>
    }
  ],
  "structure_records": Maybe<number>,
  "split": {
    "source": "feature_matrix.csv" | "variants_labelled.csv" | null,
    "caveat": string | null,
    "train_genes": Maybe<string[]>,
    "test_genes": Maybe<string[]>,
    "train_variants": Maybe<number>,
    "test_variants": Maybe<number>,
    "train_pathogenic": Maybe<number>,
    "train_benign": Maybe<number>,
    "test_pathogenic": Maybe<number>,
    "test_benign": Maybe<number>
  },
  "conservation": {
    "scored_variants": Maybe<number>,
    "missing_variants": Maybe<number>,
    "method": "homolog_identity_proxy"
  },
  "pipeline": object | null,
  "pipeline_reason": string | null
}
```

Sources:

- `genes[]` identity fields always come from `src.config.load_gene_panel("both")` (`config/genes.yaml`). 18 genes.
- Class counts come from the processed CSVs. If a CSV is missing, that count is null.
- `collected` is `len` of `clinvar_raw.csv` only. Do not return 7688 unless that file has 7688 rows.
- `structure_records` is the row count of `structure_manifest.csv`.
- Split: if `feature_matrix.csv` is loaded, apply the same `dropna` as `train_classifier.main` (every model feature except `conservation_score`, plus `label`) and call `gene_level_split`. `source` is `feature_matrix.csv`, `caveat` is null. If the matrix is missing and `variants_labelled.csv` is loaded, call `gene_level_split` on that table and set `source` to `variants_labelled.csv` with caveat: the training dropna was not applied because the feature matrix is absent. If neither table is loaded, every split field is null.
- `conservation.scored_variants` / `missing_variants` count non-null vs null `conservation_score` on the **labelled feature matrix** only. Do not return 518 / 52 unless the matrix says so.

On this checkout the labelled and VUS counts, the per-gene table, and the labelled-CSV split are computable. Collected, structure records, and conservation coverage are not. See `AUDIT.md` for the measured tables. Expected split if the matrix does not drop rows: test `ACTC1`, `CACNA1C`, `KCNE1`, `KCNH2`, `MYBPC3` (98); the other 13 genes train (472).

## `GET /metrics`

Recompute on the test fold of the same split used for training. Call `src.modeling.evaluate.evaluate_classifier` for the scalar metrics (it prints to stdout; that is existing behavior). It does **not** return the confusion matrix or curves. Also call `sklearn.metrics.confusion_matrix`, `roc_curve`, and `precision_recall_curve` on that same test fold. Positive class is pathogenic (`1`). Matrix labels are `["benign", "pathogenic"]`, so `[[tn, fp], [fn, tp]]`.

```json
{
  "split": {
    "test_genes": Maybe<string[]>,
    "n_test": Maybe<number>,
    "source": string | null
  },
  "reported": {
    "source": "README.md and app.py ML page literals; not recomputed",
    "random_forest": {
      "accuracy": 0.622, "precision": 0.647, "recall": 0.887,
      "f1": 0.748, "roc_auc": 0.588, "pr_auc": 0.706,
      "confusion_matrix": [[6, 30], [7, 55]]
    },
    "xgboost": {
      "accuracy": 0.633, "precision": 0.686, "recall": 0.774,
      "f1": 0.727, "roc_auc": 0.559, "pr_auc": 0.703,
      "confusion_matrix": [[14, 22], [14, 48]]
    }
  },
  "recomputed": null | {
    "random_forest": ModelMetrics,
    "xgboost": ModelMetrics
  },
  "comparison": null | {
    "random_forest": MetricMatch,
    "xgboost": MetricMatch
  },
  "reason": string | null
}
```

```ts
interface ModelMetrics {
  accuracy: number;
  precision: number;
  recall: number;
  f1: number;
  roc_auc: number | null;
  pr_auc: number | null;
  confusion_matrix: [[number, number], [number, number]];
  roc: { fpr: number; tpr: number; threshold: number }[];
  pr: { precision: number; recall: number; threshold: number }[];
  reason: string | null;
}

interface MetricMatch {
  accuracy: boolean;
  precision: boolean;
  recall: boolean;
  f1: boolean;
  roc_auc: boolean;
  pr_auc: boolean;
  confusion_matrix: boolean;
  all_matched: boolean;
}
```

A scalar matches when both sides round to the same 3 decimal places. The matrix matches exactly. `roc_auc` / `pr_auc` are null, with `reason`, if the test fold has a single class. `recomputed` and `comparison` are null when either model failed to load or the feature matrix is missing; `reason` says which. `reported` is always present and must be labelled as the baseline in the UI, not as a live score.

Curves use the test fold only (98 points at the published split size). Thresholds may include non-finite sentinels from scikit-learn; serialize those as null.

## `GET /shap/global`

Mean absolute SHAP per **transformed** feature, both models, probability space, pathogenic class. Row set: the labelled matrix after the training `dropna` (null conservation scores kept), passed through the saved preprocessor. Includes train-gene and test-gene rows. This is wider than `shap_explain.main`, which drops null conservation and explains XGBoost only.

```json
{
  "shap_space": "probability",
  "row_count": Maybe<number>,
  "imputed_conservation_rows": Maybe<number>,
  "note": "Mean |SHAP| on the labelled matrix after the training row filter, train and test genes together. The saved PNG is XGBoost-only and drops rows with a missing conservation score.",
  "random_forest": Maybe<GlobalFeature[]>,
  "xgboost": Maybe<GlobalFeature[]>
}
```

```ts
interface GlobalFeature {
  feature: string;
  raw_feature: string;
  mean_abs_shap: number;
}
```

Sort each list by `mean_abs_shap` descending. If a model is missing, that side is null. Do not parse `shap_global_summary.png`.

## `GET /variants`

Query: `q`, `gene`, `label` (`pathogenic` | `benign`), `split` (`train_gene` | `test_gene`), `page` (default 1), `page_size` (default 25, max 100).

`q` is a case-insensitive substring over `variant_id`, `gene_symbol`, `protein_change`, and `title`. Filters AND together. Sort is `gene_symbol`, `position`, `mut_aa`, `variant_id` so pages stay stable.

Labelled rows only.

```json
{
  "page": 1,
  "page_size": 25,
  "total": 0,
  "items": [VariantSummary],
  "reason": string | null
}
```

`reason` is set when the labelled table is missing (then `items` is `[]` and `total` is 0). A filter that matches nothing is `reason: null` and `total: 0`.

```ts
interface VariantSummary {
  variant_id: string;
  gene_symbol: string;
  panel: Panel;
  wt_aa: string;
  mut_aa: string;
  position: number;
  label: "pathogenic" | "benign";
  clinical_significance: string;
  protein_change: string | null;
  split: Split | null;
  split_reason: string | null;
  has_structure: boolean;
}
```

`has_structure` follows the structure rule below, using files loaded at startup. `split` is null when the split could not be computed.

## `GET /vus`

Same query string as `/variants` except `label` and `split` are ignored (VUS are held out of the split; every row is `vus`).

```ts
interface VusSummary {
  variant_id: string;
  gene_symbol: string;
  panel: Panel;
  wt_aa: string;
  mut_aa: string;
  position: number;
  label: "vus";
  clinical_significance: string;
  protein_change: string | null;
  predictions: { rf: Prediction; xgb: Prediction } | null;
  predictions_reason: string | null;
}
```

`predictions` is populated only when **both** of these are true: the joblib for that model loaded, and a 13-column frame could be built. Physchem always can. `panel` is on the VUS CSV. Conservation and domain features follow the null rules in `AUDIT.md` section 5. One side can be present while the other model failed: then `predictions` is non-null and the failed side carries the reason.

Do not attach SHAP here.

Scoring uses `compute_physchem_features` and `compute_conservation_and_domain_features` with conservation preloaded. It does not use `build_feature_row` in a loop (that re-reads the JSON per row).

## `GET /variants/{id}`

Resolve a labelled id first, then a VUS id. This is the case-file payload for both explorer tabs.

```json
{
  "variant_id": "ACTC1_R39H",
  "found": true,
  "gene_symbol": "ACTC1",
  "panel": "cardiomyopathy",
  "condition": "Hypertrophic / Dilated Cardiomyopathy",
  "uniprot_id": "P68032",
  "label": "benign",
  "clinical": {
    "clinvar_uid": 2681695,
    "title": "NM_005159.5(ACTC1):c.116G>A (p.Arg39His)",
    "clinical_significance": "Likely benign",
    "review_status": "no assertion criteria provided",
    "protein_change": "R39H",
    "molecular_consequence": ["missense variant"],
    "cdna_change": "c.116G>A",
    "dbsnp": "2504185412",
    "location": {
      "assembly": "GRCh38",
      "chr": "15",
      "start": "34794693",
      "stop": "34794693"
    },
    "parse_reason": null
  },
  "wt": { "one": "R", "three": "Arg" },
  "mut": { "one": "H", "three": "His" },
  "position": 39,
  "hgvs_protein": "p.Arg39His",
  "physchem": {
    "delta_hydrophobicity": 1.3,
    "delta_molecular_weight": -19.0,
    "delta_volume": -20.2,
    "delta_charge": -0.9,
    "polarity_change": false,
    "wt_is_proline": false,
    "mut_is_proline": false,
    "wt_is_glycine": false,
    "mut_is_glycine": false,
    "source": "feature_matrix.csv" | "src.features.physicochemical.compute_physchem_features",
    "reason": null
  },
  "conservation_score": Maybe<number>,
  "conservation_imputed": boolean,
  "conservation_imputed_value": Maybe<number>,
  "in_annotated_domain": Maybe<boolean>,
  "domain_type": Maybe<string | null>,
  "sequence": Maybe<{
    "uniprot_id": string,
    "sequence": string,
    "sequence_length": number,
    "domains": { "type": string, "description": string | null, "start": number | null, "end": number | null }[],
    "conservation_by_position": Record<string, number> | null
  }>,
  "predictions": { "rf": Prediction, "xgb": Prediction },
  "shap": { "rf": LocalShap, "xgb": LocalShap },
  "split": Maybe<Split>,
  "has_structure": boolean
}
```

Rules:

- Three-letter codes are the inverse of `src.preprocessing.clean_variants.AA_3TO1` (`Arg`, `His`, …). `hgvs_protein` is `p.{three}{position}{three}`.
- Prefer the feature-matrix row when it exists. Otherwise compute physchem with `compute_physchem_features`. Those floats are deterministic; they are not mock data.
- `conservation_score.value` is the raw score, null when missing. `conservation_imputed` is true when the raw score is null **and** a loaded pipeline would impute it. `conservation_imputed_value` is that pipeline's training-fold median for `conservation_score` (read it from the fitted `SimpleImputer`, do not recompute a median). If the model is not loaded, the imputed value is null and the reason says the median lives in the unfitted-from-here pipeline.
- `in_annotated_domain` is null, not false, when UniProt features were not loaded. False is a real "checked, not in a span" answer.
- `sequence.conservation_by_position` is the gene's slice of `conservation_scores.json` (`{"39": 0.91, ...}`). Null when the file is missing, even if the sequence loaded.
- `split` for a VUS is null, reason `VUS are held out of the gene-level split`.
- `clinical.parse_reason` is set when `ast.literal_eval` fails on `variation_set` or `molecular_consequence_list`. Failed fields are null. Do not `eval`.
- Local SHAP uses `TreeExplainer(..., model_output="probability")` on the single transformed row, pathogenic class, cached with `lru_cache`. `features` is the full vector sorted by absolute value, not top-5 (the report endpoint slices top-5).
- The R39H example embedded above is the measured clinical row plus physchem from `compute_physchem_features`. Conservation, sequence, predictions, and SHAP are null in this checkout. The mock layer may not copy the nulls; it must still use this clinical label (benign / Likely benign) and this physchem. Mock probabilities and SHAP are invented and must be marked `mock: true` on those objects only.

## `GET /variants/{id}/structure`

```json
{
  "variant_id": "ACTC1_R39H",
  "has_structure": false,
  "reason": "outputs/ACTC1_R39H_wt.pdb is not present",
  "poc": {
    "defined": true,
    "pdb_id": "9B3R",
    "chain": "A",
    "residue": 39,
    "wt_aa": "R",
    "mut_aa": "H",
    "source_note": "Hardcoded proof-of-concept in src/structure/mutate_pdb.py. Not read from structure_manifest.csv."
  },
  "files": {
    "wt": { "name": "ACTC1_R39H_wt.pdb", "present": false, "url": null },
    "mutant": { "name": "ACTC1_R39H_mutant.pdb", "present": false, "url": null }
  },
  "rmsd_angstrom": Maybe<number>,
  "ddg": Maybe<number>,
  "method_note": "The mutant file is a residue-identity swap: ARG 39 is renamed to HIS on chain A and the coordinates are unchanged (src/structure/mutate_pdb.py). No side-chain rebuild and no minimization. A Cα RMSD of 0.00 Å is a method artifact, not evidence that the substitution has no structural effect.",
  "esmfold_note": "src/structure/predict_mutant_structure.py can call the public ESMFold endpoint. That endpoint timed out during the project. This API does not call it. The proof-of-concept PDBs are not ESMFold models.",
  "known_structures": Maybe<KnownStructure[]>
}
```

```ts
interface KnownStructure {
  gene_symbol: string;
  uniprot_id: string;
  source: "PDB" | "AlphaFoldDB" | null;
  structure_id: string | null;
  note: "fetch_structures keeps a single structure per gene (first RCSB hit, else AlphaFold DB). This is not a complete PDB list.";
}
```

- `has_structure` is true only for `ACTC1_R39H` when both PDB files exist. Then `files.*.url` is `/files/pdb/ACTC1_R39H_wt.pdb` and `/files/pdb/ACTC1_R39H_mutant.pdb`, and `rmsd_angstrom` is `compute_rmsd` from `src.structure.stability_analysis` (chain A). If that call throws, the value is null and `reason` is the exception.
- For every other id, `poc` is null, `files` is null, `has_structure` is false, `reason` is `Structural proof-of-concept is implemented only for ACTC1_R39H`, and `method_note` is null. `esmfold_note` stays.
- `ddg` is always `{ "value": null, "reason": "FoldX BuildModel is not implemented in src/structure/stability_analysis.py. ΔΔG was not computed." }`.
- `known_structures` is the manifest rows for this variant's gene, or null when the manifest is missing.
- Mock mode, in the frontend only: when `has_structure` is false and the id is `ACTC1_R39H`, the viewer loads `https://files.rcsb.org/download/9B3R.pdb` for both WT and mutant. The real API does not proxy RCSB and does not put that URL in `files`.

## `GET /variants/{id}/report`

Assemble a dossier. Use `src.report.generate_report.make_report` once per model that has a probability and a SHAP vector. `variant` is the HGVS protein string. `confidence_score` is `prob_pathogenic`. `key_features` is the first 5 local-SHAP rows (`feature`, `shap_value`). `structural_rmsd` is the recomputed RMSD or `None`. `predicted_ddg` is `None`. `stability_interpretation` for `ACTC1_R39H` is the method-artifact note, not the string from `interpret_stability`. For every other variant it is `Structural comparison was not run for this variant. ΔΔG was not computed.`

```json
{
  "variant_id": "ACTC1_R39H",
  "disclaimer": "Research prototype — not for clinical use",
  "clinical": {},
  "physchem": {},
  "conservation": {
    "score": Maybe<number>,
    "imputed": boolean,
    "method": "Homolog-identity proxy from src/features/build_conservation_scores.py. Not ConSurf, phyloP, or phastCons."
  },
  "predictions": { "rf": Prediction, "xgb": Prediction },
  "agreement": "agree" | "disagree" | null,
  "top_shap": { "rf": ShapFeature[], "xgb": ShapFeature[] },
  "structure": {},
  "summaries": {
    "rf": Maybe<string>,
    "xgb": Maybe<string>
  },
  "interpretation": string,
  "missing": [{ "field": "predicted_ddg", "reason": "..." }]
}
```

`clinical`, `physchem`, and `structure` are the same objects as on the detail and structure endpoints.

`agreement` compares the two predicted labels when both exist, else null.

`interpretation` is template prose, conservative, and must include:

- the clinical label and the review status, kept separate from the model labels
- whether the gene was in the train fold or the test fold (or that the row is VUS)
- that `panel` is a gene-level covariate in the model and `position` is an absolute index
- the conservation method, or that the score was not computed / was imputed
- the structure method note when this is the proof-of-concept, otherwise that structure was not run
- the disclaimer sentence again

`summaries.*.value` is `plain_language_summary` from `make_report`. When a model did not run, that summary is null and `interpretation` says the model was not run. Do not invent a confidence.

`GET /report` also accepts `?format=markdown` and `?format=json` as the raw dossier (same JSON, or a Markdown rendering of it). The default is JSON. The frontend downloads those two renderings; it does not assemble a second report from scratch.

## `GET /files/pdb/{name}`

Whitelist:

- `name` matches `^[A-Za-z0-9_.-]+\.pdb$`
- resolved path's parent is exactly `outputs/`
- the file exists

Otherwise HTTP 404 with a reason. No `..`, no absolute paths, no other suffix. Response is the file bytes with `Content-Type: chemical/x-pdb`.

## `GET /ablation`

No ablation script exists. Do not train.

Look only for `outputs/ablation_results.json`. If it is absent:

```json
{
  "value": null,
  "reason": "No ablation has been run. outputs/ablation_results.json is not present. This endpoint does not train.",
  "plan": {
    "sets": [
      "A_physchem",
      "B_physchem_conservation",
      "C_physchem_conservation_domain"
    ],
    "position": ["with_position", "without_position"],
    "note": "The current model also includes panel, which is constant per gene. domain_type is not a model input; the domain feature in set C is in_annotated_domain. Gene-level split and class weighting stay as in train_classifier.py."
  }
}
```

If a future file appears, return `{ "value": <parsed json>, "reason": null, "plan": { ...same plan... } }` without validating it into fake metrics. The UI shows the plan either way, and shows results only when `value` is an object.

## Mock client

`VITE_USE_MOCKS` defaults to `true`. Mocks live in `web/src/mocks/` and return this schema.

Mock-only additions, which the real API omits:

- Every page shows a visible `MOCK DATA` ribbon.
- Objects that are not measured (probabilities, SHAP, ROC/PR points, global mean |SHAP|, conservation scores, sequences, RMSD) carry `"mock": true`.
- Clinical fields, gene panel, the measured 570/424/146 and 6,122 counts, the per-gene table, and the reproduced 13/5 gene lists may be served from the processed CSVs' measured values **without** `mock: true`, because they were read from disk during this audit. Copy them into mock JSON; do not have the browser read the CSV unless the mock layer does it itself.
- `ACTC1_R39H` mock predictions and SHAP must be plausible and labelled mock. Clinical label stays benign. Physchem stays the `compute_physchem_features("R","H")` values. Structure viewer uses RCSB `9B3R` for both WT and mutant. `rmsd_angstrom` may be the documented method-artifact `0.00` with `mock: true` and the method note always visible. `ddg.value` stays null.
- `CACNA1C_G1365R` is the test-gene variant without structure.
- Collected count 7,688 and conservation 518/52 may appear in mocks because the README states them, and they must be tagged `mock: true` until the raw files exist. The real API returns null for those until the files exist.
- `pipeline.ui.react` in mocks may be `done` once the UI is built. The status file in `api/` should be updated in the same phase, not before.

## Out of scope for the API

- Writing or retraining models
- Editing `src/`, `app.py`, `config/`, `data/`, `tests/`, `README.md`, `requirements.txt`, `.gitignore`
- Proxying NCBI, UniProt, RCSB, or ESMFold
- Serving `example_variant_report.json`
