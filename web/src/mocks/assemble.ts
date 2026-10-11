import type {
  Ablation,
  Clinical,
  GlobalShap,
  Health,
  Label,
  ListQuery,
  LocalShap,
  Maybe,
  Metrics,
  Overview,
  Page,
  Panel,
  Physchem,
  Prediction,
  Report,
  SequenceWindow,
  Split,
  StructurePayload,
  SuggestHit,
  VariantDetail,
  VariantSummary,
  VusSummary,
} from '../api/types.ts'

export interface CopyFields {
  uniprot_reason: string
  conservation_reason: string
  median_reason: string
  model_reason: string
  shap_reason: string
  method_note: string
  ddg_reason: string
  other_stability: string
  conservation_method: string
  other_structure_reason: string
  esmfold_note: string
  disclaimer: string
  rcsb_pdb: string
}

interface RecordRow {
  variant_id: string
  gene_symbol: string
  panel: Panel
  condition: string | null
  uniprot_id: string | null
  wt_aa: string
  mut_aa: string
  wt_three: string | null
  mut_three: string | null
  position: number
  hgvs_protein: string | null
  label: Label
  clinical: Clinical
  split: Split | null
  has_structure: boolean
  physchem: Physchem
}

interface RichCase {
  conservation_score: Maybe<number>
  predictions: { rf: Prediction; xgb: Prediction }
  shap: { rf: LocalShap; xgb: LocalShap }
  sequence: Maybe<SequenceWindow>
  summaries: { rf: Maybe<string>; xgb: Maybe<string> }
}

interface Catalog {
  labelled: RecordRow[]
  vus: RecordRow[]
  rich: Record<string, RichCase>
  overview: Overview
  health: Health
  metrics: Metrics
  shap: GlobalShap
  ablation: Ablation
  copy: CopyFields
  byId: Map<string, { row: RecordRow; kind: 'labelled' | 'vus' }>
}

let pending: Promise<Catalog> | null = null

function as<T>(value: unknown): T {
  return value as T
}

export function loadCatalog() {
  if (!pending) {
    pending = (async () => {
      await new Promise((resolve) => setTimeout(resolve, 160))
      const [labelledMod, vusMod, richMod, overviewMod, healthMod, metricsMod, shapMod, ablationMod, copyMod] =
        await Promise.all([
          import('./data/labelled.json'),
          import('./data/vus.json'),
          import('./data/rich.json'),
          import('./data/overview.json'),
          import('./data/health.json'),
          import('./data/metrics.json'),
          import('./data/shap.json'),
          import('./data/ablation.json'),
          import('./data/copy.json'),
        ])
      const labelled = as<RecordRow[]>(labelledMod.default)
      const vus = as<RecordRow[]>(vusMod.default)
      const byId = new Map<string, { row: RecordRow; kind: 'labelled' | 'vus' }>()
      for (const row of labelled) byId.set(row.variant_id, { row, kind: 'labelled' })
      for (const row of vus) byId.set(row.variant_id, { row, kind: 'vus' })
      return {
        labelled,
        vus,
        rich: as<Record<string, RichCase>>(richMod.default),
        overview: as<Overview>(overviewMod.default),
        health: as<Health>(healthMod.default),
        metrics: as<Metrics>(metricsMod.default),
        shap: as<GlobalShap>(shapMod.default),
        ablation: as<Ablation>(ablationMod.default),
        copy: as<CopyFields>(copyMod.default),
        byId,
      }
    })().catch((error: unknown) => {
      pending = null
      throw error
    })
  }
  return pending
}

function absent<T>(reason: string): Maybe<T> {
  return { value: null, reason }
}

function emptyPrediction(reason: string): Prediction {
  return { label: null, prob_pathogenic: null, threshold: 0.5, reason }
}

function emptyShap(reason: string): LocalShap {
  return {
    shap_space: 'probability',
    base_value: null,
    output_value: null,
    features: [],
    reason,
  }
}

function assertPage(query: ListQuery) {
  const page = query.page ?? 1
  const pageSize = query.page_size ?? 25
  if (!Number.isInteger(page) || page < 1) {
    throw Object.assign(new Error('page must be an integer >= 1'), { status: 400 })
  }
  if (!Number.isInteger(pageSize) || pageSize < 1 || pageSize > 100) {
    throw Object.assign(new Error('page_size must be an integer from 1 to 100'), { status: 400 })
  }
  if (query.label && query.label !== 'pathogenic' && query.label !== 'benign') {
    throw Object.assign(new Error('label must be pathogenic or benign'), { status: 400 })
  }
  if (query.split && query.split !== 'train_gene' && query.split !== 'test_gene') {
    throw Object.assign(new Error('split must be train_gene or test_gene'), { status: 400 })
  }
  return { page, pageSize }
}

function matchesQuery(row: RecordRow, query: ListQuery, ignoreLabelSplit: boolean) {
  if (query.gene && row.gene_symbol !== query.gene) return false
  if (!ignoreLabelSplit && query.label && row.label !== query.label) return false
  if (!ignoreLabelSplit && query.split && row.split !== query.split) return false
  const needle = query.q?.trim().toLowerCase()
  if (!needle) return true
  const haystack = [row.variant_id, row.gene_symbol, row.clinical.protein_change ?? '', row.clinical.title ?? '']
  return haystack.some((part) => part.toLowerCase().includes(needle))
}

function slicePage<T>(items: T[], page: number, pageSize: number): Page<T> {
  const start = (page - 1) * pageSize
  return {
    page,
    page_size: pageSize,
    total: items.length,
    items: items.slice(start, start + pageSize),
    reason: null,
  }
}

export async function overview() {
  return (await loadCatalog()).overview
}

export async function health() {
  return (await loadCatalog()).health
}

export async function metrics() {
  return (await loadCatalog()).metrics
}

export async function globalShap() {
  return (await loadCatalog()).shap
}

export async function ablation() {
  return (await loadCatalog()).ablation
}

export async function variants(query: ListQuery): Promise<Page<VariantSummary>> {
  const { page, pageSize } = assertPage(query)
  const catalog = await loadCatalog()
  const items = catalog.labelled.filter((row) => matchesQuery(row, query, false)).map(toSummary)
  return slicePage(items, page, pageSize)
}

export async function vus(query: ListQuery): Promise<Page<VusSummary>> {
  const { page, pageSize } = assertPage(query)
  const catalog = await loadCatalog()
  const items = catalog.vus.filter((row) => matchesQuery(row, query, true)).map((row) => toVus(row, catalog.copy))
  return slicePage(items, page, pageSize)
}

function toSummary(row: RecordRow): VariantSummary {
  return {
    variant_id: row.variant_id,
    gene_symbol: row.gene_symbol,
    panel: row.panel,
    wt_aa: row.wt_aa,
    mut_aa: row.mut_aa,
    position: row.position,
    label: row.label === 'benign' ? 'benign' : 'pathogenic',
    clinical_significance: row.clinical.clinical_significance ?? '',
    protein_change: row.clinical.protein_change,
    split: row.split,
    split_reason: row.split ? null : 'Split was not computed.',
    has_structure: row.has_structure,
  }
}

function toVus(row: RecordRow, copy: CopyFields): VusSummary {
  return {
    variant_id: row.variant_id,
    gene_symbol: row.gene_symbol,
    panel: row.panel,
    wt_aa: row.wt_aa,
    mut_aa: row.mut_aa,
    position: row.position,
    label: 'vus',
    clinical_significance: row.clinical.clinical_significance ?? '',
    protein_change: row.clinical.protein_change,
    predictions: null,
    predictions_reason: copy.uniprot_reason,
  }
}

function requireRow(catalog: Catalog, id: string) {
  const hit = catalog.byId.get(id)
  if (!hit) {
    throw Object.assign(new Error(`No labelled or VUS row has id ${id}`), { status: 404 })
  }
  return hit
}

export async function variant(id: string): Promise<VariantDetail> {
  const catalog = await loadCatalog()
  const { row, kind } = requireRow(catalog, id)
  const rich = catalog.rich[id]
  const copy = catalog.copy
  const score = rich?.conservation_score ?? absent(copy.conservation_reason)
  return {
    variant_id: row.variant_id,
    found: true,
    gene_symbol: row.gene_symbol,
    panel: row.panel,
    condition: row.condition,
    uniprot_id: row.uniprot_id,
    label: row.label,
    clinical: row.clinical,
    wt: { one: row.wt_aa, three: row.wt_three },
    mut: { one: row.mut_aa, three: row.mut_three },
    position: row.position,
    hgvs_protein: row.hgvs_protein,
    physchem: row.physchem,
    conservation_score: score,
    conservation_imputed: false,
    conservation_imputed_value: score.value == null ? absent(copy.median_reason) : { value: null, reason: null },
    in_annotated_domain: absent(copy.uniprot_reason),
    domain_type: { value: null, reason: null },
    sequence: rich?.sequence ?? absent(copy.uniprot_reason),
    predictions: rich?.predictions ?? {
      rf: emptyPrediction(copy.model_reason),
      xgb: emptyPrediction(copy.model_reason),
    },
    shap: rich?.shap ?? {
      rf: emptyShap(copy.shap_reason),
      xgb: emptyShap(copy.shap_reason),
    },
    split:
      kind === 'vus'
        ? { value: null, reason: 'VUS are held out of the gene-level split.' }
        : { value: row.split, reason: row.split ? null : 'Split was not computed.' },
    has_structure: row.has_structure,
  }
}

export async function structure(id: string): Promise<StructurePayload> {
  const catalog = await loadCatalog()
  const { row } = requireRow(catalog, id)
  const copy = catalog.copy
  const known: StructurePayload['known_structures'] = {
    value: null,
    reason:
      catalog.overview.structure_records.value == null
        ? (catalog.overview.structure_records.reason ?? 'Structure manifest is not present.')
        : 'Per-gene manifest rows are not in the mock catalog.',
  }
  const ddg = absent<number>(copy.ddg_reason)
  if (id !== 'ACTC1_R39H') {
    return {
      variant_id: id,
      has_structure: false,
      reason: copy.other_structure_reason,
      poc: null,
      files: null,
      rmsd_angstrom: absent(copy.other_stability),
      ddg,
      method_note: null,
      esmfold_note: copy.esmfold_note,
      known_structures: known,
      mock_viewer: null,
    }
  }
  return {
    variant_id: id,
    has_structure: row.has_structure,
    reason: row.has_structure ? null : 'outputs/ACTC1_R39H_wt.pdb is not present',
    poc: {
      defined: true,
      pdb_id: '9B3R',
      chain: 'A',
      residue: row.position,
      wt_aa: row.wt_aa,
      mut_aa: row.mut_aa,
      source_note: 'Hardcoded proof-of-concept in src/structure/mutate_pdb.py. Not read from structure_manifest.csv.',
    },
    files: {
      wt: { name: 'ACTC1_R39H_wt.pdb', present: false, url: null },
      mutant: { name: 'ACTC1_R39H_mutant.pdb', present: false, url: null },
    },
    rmsd_angstrom: { value: 0, reason: null, mock: true },
    ddg,
    method_note: copy.method_note,
    esmfold_note: copy.esmfold_note,
    known_structures: known,
    mock_viewer: { wt: copy.rcsb_pdb, mutant: copy.rcsb_pdb, mock: true },
  }
}

function agreementOf(predictions: VariantDetail['predictions']): Report['agreement'] {
  const left = predictions.rf.label
  const right = predictions.xgb.label
  if (!left || !right) return null
  return left === right ? 'agree' : 'disagree'
}

function interpretation(detail: VariantDetail, struct: StructurePayload, copy: CopyFields) {
  const clinical = detail.clinical.clinical_significance ?? 'unavailable'
  const review = detail.clinical.review_status ?? 'review status unavailable'
  const rf = detail.predictions.rf
  const xgb = detail.predictions.xgb
  const rfText = rf.label ? `Random Forest ${rf.label}` : 'Random Forest was not run'
  const xgbText = xgb.label ? `XGBoost ${xgb.label}` : 'XGBoost was not run'
  const fold =
    detail.label === 'vus'
      ? 'This row is VUS and was held out of the gene-level split.'
      : detail.split.value === 'train_gene'
        ? 'This gene was in the training fold.'
        : detail.split.value === 'test_gene'
          ? 'This gene was in the test fold.'
          : 'Fold membership was not computed.'
  const conservation =
    detail.conservation_score.value == null
      ? 'The conservation score was not computed.'
      : `Conservation score ${detail.conservation_score.value}. Method: ${copy.conservation_method}`
  const structureText = struct.method_note ?? 'Structural comparison was not run for this variant.'
  return [
    `Clinical record: ${clinical} (${review}). Model labels are separate: ${rfText}; ${xgbText}.`,
    fold,
    'The model includes panel as a gene-level covariate. Position is an absolute residue index and is not comparable across proteins.',
    conservation,
    structureText,
    copy.disclaimer,
  ].join(' ')
}

export async function report(id: string, format: 'json' | 'markdown' = 'json'): Promise<Report | string> {
  const catalog = await loadCatalog()
  const detail = await variant(id)
  const struct = await structure(id)
  const rich = catalog.rich[id]
  const copy = catalog.copy
  const missing: Report['missing'] = []
  if (detail.conservation_score.value == null) {
    missing.push({ field: 'conservation_score', reason: detail.conservation_score.reason })
  }
  if (detail.in_annotated_domain.value == null) {
    missing.push({ field: 'in_annotated_domain', reason: detail.in_annotated_domain.reason })
  }
  if (detail.predictions.rf.prob_pathogenic == null) {
    missing.push({ field: 'rf_probability', reason: detail.predictions.rf.reason })
  }
  if (detail.predictions.xgb.prob_pathogenic == null) {
    missing.push({ field: 'xgb_probability', reason: detail.predictions.xgb.reason })
  }
  if (struct.rmsd_angstrom.value == null) {
    missing.push({ field: 'structural_rmsd', reason: struct.rmsd_angstrom.reason })
  }
  if (struct.ddg.value == null) {
    missing.push({ field: 'predicted_ddg', reason: struct.ddg.reason })
  }
  const payload: Report = {
    variant_id: id,
    disclaimer: copy.disclaimer,
    clinical: detail.clinical,
    physchem: detail.physchem,
    conservation: {
      score: detail.conservation_score,
      imputed: detail.conservation_imputed,
      method: copy.conservation_method,
    },
    predictions: detail.predictions,
    agreement: agreementOf(detail.predictions),
    top_shap: {
      rf: detail.shap.rf.features.slice(0, 5),
      xgb: detail.shap.xgb.features.slice(0, 5),
    },
    structure: struct,
    summaries: rich?.summaries ?? {
      rf: absent(detail.predictions.rf.reason ?? 'The model was not run.'),
      xgb: absent(detail.predictions.xgb.reason ?? 'The model was not run.'),
    },
    interpretation: interpretation(detail, struct, copy),
    missing,
  }
  if (format === 'markdown') return renderMarkdown(payload, detail)
  return payload
}

function renderMarkdown(report: Report, detail: VariantDetail) {
  const lines = [
    `# ${detail.hgvs_protein ?? report.variant_id}`,
    '',
    `> ${report.disclaimer}`,
    '',
    `Variant: ${report.variant_id}`,
    `Gene: ${detail.gene_symbol} (${detail.panel})`,
    `Clinical: ${detail.clinical.clinical_significance ?? 'Not computed'}`,
    `Review: ${detail.clinical.review_status ?? 'Not computed'}`,
    `ClinVar title: ${detail.clinical.title ?? 'Not computed'}`,
    '',
    '## Models',
    modelLine('Random Forest', report.predictions.rf),
    modelLine('XGBoost', report.predictions.xgb),
    `Agreement: ${report.agreement ?? 'Not computed'}`,
    '',
    '## Interpretation',
    report.interpretation,
    '',
    '## Missing',
    ...report.missing.map((item) => `- ${item.field}: ${item.reason ?? 'Not computed'}`),
  ]
  return lines.join('\n')
}

function modelLine(name: string, prediction: Prediction) {
  if (prediction.prob_pathogenic == null) return `${name}: Not computed. ${prediction.reason ?? ''}`
  return `${name}: ${prediction.label} (P pathogenic ${prediction.prob_pathogenic})`
}

export async function suggest(q: string): Promise<SuggestHit[]> {
  const catalog = await loadCatalog()
  const needle = q.trim().toLowerCase()
  const pool = needle
    ? [...catalog.labelled, ...catalog.vus].filter((row) =>
        matchesQuery(row, { q: needle }, true),
      )
    : catalog.overview.smoke
      ? [
          catalog.byId.get(catalog.overview.smoke.structure)?.row,
          catalog.byId.get(catalog.overview.smoke.no_structure)?.row,
          catalog.byId.get(catalog.overview.smoke.train_gene)?.row,
        ].filter((row): row is RecordRow => Boolean(row))
      : catalog.labelled.slice(0, 8)
  return pool.slice(0, 12).map((row) => ({
    variant_id: row.variant_id,
    gene_symbol: row.gene_symbol,
    label: row.label,
    protein_change: row.clinical.protein_change,
    position: row.position,
  }))
}
