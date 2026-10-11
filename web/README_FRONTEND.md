# Specimen Archive frontend

React bench for the cardiovascular variant project. It reads the read-only API in `api/` or, by default, a mock client with the same JSON shapes.

## Run

From the repo root, with mocks (no API required):

```powershell
cd web
npm install
npm run dev
```

Open http://127.0.0.1:5173. Or from the repo root, `./dev.ps1` starts Uvicorn on port 8000 and Vite on port 5173.

`VITE_USE_MOCKS` defaults to on. It is off only when the variable is exactly `false`. Copy `web/.env.example` to `web/.env.local` and set:

```
VITE_USE_MOCKS=false
VITE_API_URL=http://127.0.0.1:8000
```

Start the API from the repo root so `import src` resolves:

```powershell
python -m uvicorn api.main:app --host 127.0.0.1 --port 8000
```

The mock ribbon is visible on every page while mocks are on. Clinical labels, physchem deltas, gene identity, and the labelled/VUS counts are measured from the processed tables. Probabilities, SHAP, ROC/PR points, global mean |SHAP|, conservation scores, sequence flanks, the collected count 7688, and conservation coverage 518/52 are marked mock.

## Pages

| Path | What it shows |
| --- | --- |
| `/` | Sealed drawer. Open it to enter the bench |
| `/bench` | Counts, gene split, pipeline helix, per-gene stacks, metrics, global SHAP, artifact health |
| `/explore?variant=` | Labelled and VUS lists, filters, case file |
| `/structure/:id` | Molstar panes, RMSD, empty ΔΔG |
| `/report/:id` | Dossier, Markdown and JSON downloads |
| `/models` | Metrics, leakage notes, ablation plan |
| `/about` | Sources, split, limits |

Smoke cases in the mock catalog: `ACTC1_R39H` (structure), `CACNA1C_G1365R` (test gene, no structure), `MYH7_F252Y` (train gene). Find them from the bench chips or Ctrl/Cmd+K.

The structure viewer loads the Molstar browser bundle (`molstar/build/viewer/molstar.js`). `pdbe-molstar` is installed for the same viewer family; the prebuilt Molstar script is what the page executes, because the plugin's CommonJS graph does not bundle cleanly through Vite on Node 20. In mock mode, both wild-type and mutant panes load deposited PDB 9B3R from RCSB. The real API does not proxy RCSB. When local PDBs exist, the viewer uses `/files/pdb/...` from the API.

## Regenerate mocks

From the repo root, with the project dependencies importable (`src` on the path, including pandas and scikit-learn):

```powershell
python web/scripts/generate_mocks.py
```

The script reads `data/processed/variants_labelled.csv` and `variants_vus.csv`, computes physchem with `compute_physchem_features`, and reproduces the gene split with `gene_level_split`. It does not train.

To fill the real API, the pipeline has to be run from the repo root. `python -m src.pipeline --stage all` does **not** build conservation scores. The order that matches the training code is:

```powershell
python -m src.pipeline --stage collect
python -m src.pipeline --stage preprocess
python -m src.features.build_conservation_scores
python -m src.pipeline --stage features
python -m src.pipeline --stage train
python -m src.pipeline --stage explain
```

Structure is manual: place `outputs/ACTC1_R39H_wt.pdb` (PDB 9B3R), then run `python -m src.structure.mutate_pdb` and the stability step. Do not import `mutate_pdb` from the API. FoldX ΔΔG stays unimplemented. Do not call ESMFold from the API.

## Checks

```powershell
cd web
npm run lint
npx tsc -b
npm run build
```
