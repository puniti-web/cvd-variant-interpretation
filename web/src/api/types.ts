export type Label = 'pathogenic' | 'benign' | 'vus'
export type ModelLabel = 'pathogenic' | 'benign'
export type Split = 'train_gene' | 'test_gene'
export type Panel = 'cardiomyopathy' | 'arrhythmia'

export interface Maybe<T> {
  value: T | null
  reason: string | null
  mock?: boolean
}

export interface Artifact {
  key: string
  path: string
  present: boolean
  loaded: boolean
  reason: string | null
}

export interface Health {
  status: string
  artifacts: Artifact[]
}

export interface Prediction {
  label: ModelLabel | null
  prob_pathogenic: number | null
  threshold: number
  reason: string | null
  mock?: boolean
}

export interface ShapFeature {
  feature: string
  raw_feature: string
  shap_value: number
  feature_value_transformed: number
  feature_value_raw: number | string | boolean | null
}

export interface LocalShap {
  shap_space: 'probability'
  base_value: number | null
  output_value: number | null
  features: ShapFeature[]
  reason: string | null
  mock?: boolean
}

export interface PipelinePart {
  name: string
  status: string
}

export interface PipelineStep {
  id: number
  name: string
  status: string
  note: string
  parts?: PipelinePart[]
}

export interface PipelineStatus {
  steps: PipelineStep[]
  ui: { streamlit: string; react: string }
}

export interface GeneCount {
  symbol: string
  panel: Panel
  condition: string
  uniprot_id: string
  benign: Maybe<number>
  pathogenic: Maybe<number>
  vus: Maybe<number>
}

export interface Overview {
  counts: {
    collected: Maybe<number>
    labelled: Maybe<number>
    pathogenic: Maybe<number>
    benign: Maybe<number>
    vus: Maybe<number>
  }
  genes: GeneCount[]
  structure_records: Maybe<number>
  split: {
    source: string | null
    caveat: string | null
    train_genes: Maybe<string[]>
    test_genes: Maybe<string[]>
    train_variants: Maybe<number>
    test_variants: Maybe<number>
    train_pathogenic: Maybe<number>
    train_benign: Maybe<number>
    test_pathogenic: Maybe<number>
    test_benign: Maybe<number>
  }
  conservation: {
    scored_variants: Maybe<number>
    missing_variants: Maybe<number>
    method: string
  }
  pipeline: PipelineStatus | null
  pipeline_reason: string | null
  smoke?: {
    structure: string
    no_structure: string
    train_gene: string
  }
}

export interface CurvePoint {
  fpr?: number
  tpr?: number
  precision?: number
  recall?: number
  threshold: number | null
}

export interface ModelMetrics {
  accuracy: number
  precision: number
  recall: number
  f1: number
  roc_auc: number | null
  pr_auc: number | null
  confusion_matrix: [[number, number], [number, number]]
  roc?: CurvePoint[]
  pr?: CurvePoint[]
  reason: string | null
  mock?: boolean
}

export interface MetricMatch {
  accuracy: boolean
  precision: boolean
  recall: boolean
  f1: boolean
  roc_auc: boolean
  pr_auc: boolean
  confusion_matrix: boolean
  all_matched: boolean
}

export interface CvStat {
  mean: number | null
  ci95: [number, number] | null
}

export interface Baseline {
  coverage: number
  roc_auc: number
  pr_auc: number
  model_roc_auc_same_variants: number
  model_pr_auc_same_variants: number
}

export interface PerGene {
  gene_symbol: string
  n: number
  pathogenic: number
  benign: number
  roc_auc: number | null
  accuracy: number
}

export interface CrossValidation {
  scheme: Record<string, string> | null
  n_variants: number | null
  n_genes: number | null
  best_model: string | null
  gene_grouped: Record<string, Record<string, CvStat>>
  within_gene: Record<string, Record<string, CvStat>>
  baselines: Record<string, Baseline> | null
  per_gene: PerGene[] | null
}

export interface Metrics {
  split: {
    test_genes: Maybe<string[]>
    n_test: Maybe<number>
    source: string | null
  }
  reported: {
    source: string
    random_forest: ModelMetrics
    xgboost: ModelMetrics
  } | null
  cross_validation?: CrossValidation | null
  recomputed: null | {
    mock?: boolean
    random_forest: ModelMetrics
    xgboost: ModelMetrics
  }
  comparison: null | {
    random_forest: MetricMatch
    xgboost: MetricMatch
  }
  reason: string | null
}

export interface GlobalFeature {
  feature: string
  raw_feature: string
  mean_abs_shap: number
}

export interface GlobalShap {
  shap_space: string
  row_count: Maybe<number>
  imputed_conservation_rows: Maybe<number>
  note: string
  random_forest: Maybe<GlobalFeature[]>
  xgboost: Maybe<GlobalFeature[]>
  mock?: boolean
}

export interface Ablation {
  value: unknown
  reason: string | null
  plan: {
    sets: string[]
    position: string[]
    note: string
  }
}

export interface VariantSummary {
  variant_id: string
  gene_symbol: string
  panel: Panel
  wt_aa: string
  mut_aa: string
  position: number
  label: 'pathogenic' | 'benign'
  clinical_significance: string
  protein_change: string | null
  split: Split | null
  split_reason: string | null
  has_structure: boolean
}

export interface VusSummary {
  variant_id: string
  gene_symbol: string
  panel: Panel
  wt_aa: string
  mut_aa: string
  position: number
  label: 'vus'
  clinical_significance: string
  protein_change: string | null
  predictions: { rf: Prediction; xgb: Prediction } | null
  predictions_reason: string | null
  ranking?: VusRanking | null
}

export interface VusRanking {
  rank: number | null
  probability: number | null
  calibrated_probability: number | null
  tier: string | null
  foldx_ddg: number | null
  ddg_category: string | null
  evidence_agree: number | null
  evidence_informative: number | null
  top_reasons: string | null
}

export interface Page<T> {
  page: number
  page_size: number
  total: number
  items: T[]
  reason: string | null
}

export interface ListQuery {
  q?: string
  gene?: string
  label?: string
  split?: string
  page?: number
  page_size?: number
}

export interface Clinical {
  clinvar_uid: number | null
  title: string | null
  clinical_significance: string | null
  review_status: string | null
  protein_change: string | null
  molecular_consequence: string[] | null
  cdna_change: string | null
  dbsnp: string | null
  location: {
    assembly: string | null
    chr: string | null
    start: string | null
    stop: string | null
  } | null
  parse_reason: string | null
}

export interface Physchem {
  delta_hydrophobicity: number
  delta_molecular_weight: number
  delta_volume: number
  delta_charge: number
  polarity_change: boolean
  wt_is_proline: boolean
  mut_is_proline: boolean
  wt_is_glycine: boolean
  mut_is_glycine: boolean
  source: string
  reason: string | null
}

export interface SequenceWindow {
  uniprot_id: string
  sequence: string
  sequence_length: number
  sequence_start?: number
  domains: {
    type: string
    description: string | null
    start: number | null
    end: number | null
  }[]
  conservation_by_position: Record<string, number> | null
  excerpt?: boolean
  note?: string
}

export interface VariantDetail {
  variant_id: string
  found: true
  gene_symbol: string
  panel: Panel
  condition: string | null
  uniprot_id: string | null
  label: Label
  clinical: Clinical
  wt: { one: string; three: string | null }
  mut: { one: string; three: string | null }
  position: number
  hgvs_protein: string | null
  physchem: Physchem
  conservation_score: Maybe<number>
  conservation_imputed: boolean
  conservation_imputed_value: Maybe<number>
  in_annotated_domain: Maybe<boolean>
  domain_type: Maybe<string | null>
  sequence: Maybe<SequenceWindow>
  predictions: { rf: Prediction; xgb: Prediction }
  shap: { rf: LocalShap; xgb: LocalShap }
  split: Maybe<Split>
  has_structure: boolean
}

export interface PdbFile {
  name: string
  present: boolean
  url: string | null
}

export interface KnownStructure {
  gene_symbol: string
  uniprot_id: string
  source: 'PDB' | 'AlphaFoldDB' | null
  structure_id: string | null
  note: string
}

export interface StructurePayload {
  variant_id: string
  has_structure: boolean
  reason: string | null
  poc: null | {
    defined: boolean
    pdb_id: string
    chain: string
    residue: number
    wt_aa: string
    mut_aa: string
    source_note: string
  }
  files: null | { wt: PdbFile; mutant: PdbFile }
  rmsd_angstrom: Maybe<number>
  ddg: Maybe<number>
  method_note: string | null
  esmfold_note: string
  known_structures: Maybe<KnownStructure[]>
  mock_viewer?: null | { wt: string; mutant: string; mock: true }
  ddg_category?: string
  model?: StructureModel | null
}

export interface ResidueContext {
  plddt: number | null
  rsa: number | null
  contact_count: number | null
  sec_struct: string | null
}

export interface StructureModel {
  source: string
  name: string
  source_url: string
  uniprot_id: string | null
  chain: string
  residue: number
  wt_aa: string
  mut_aa: string
  context: ResidueContext | null
}

export interface EvidenceLine {
  name: string
  direction: 'pathogenic' | 'benign' | 'neutral' | 'unavailable'
  detail: string
  independent: boolean
}

export interface IntegratedContribution {
  feature: string
  label: string
  value_display: string
  shap: number
  direction: 'pathogenic' | 'benign'
}

export interface IntegratedReport {
  variant: string
  gene_symbol: string
  protein_change: string
  condition: string
  panel: string
  clinvar: { label: string | null; review_stars: number; clinvar_uid: number | null }
  prediction: {
    probability: number
    source: string
    label: 'pathogenic' | 'benign'
    calibrated_at_training_prevalence: number | null
  }
  tier: string
  key_features: IntegratedContribution[]
  stability: { foldx_ddg: number | null; category: string }
  structure: Record<string, number | string | null>
  evidence: EvidenceLine[]
  conclusion: string
  flags: string[]
  disclaimer: string
  markdown: string
}

export interface Report {
  variant_id: string
  disclaimer: string
  clinical: Clinical
  physchem: Physchem
  conservation: {
    score: Maybe<number>
    imputed: boolean
    method: string
  }
  predictions: { rf: Prediction; xgb: Prediction }
  agreement: 'agree' | 'disagree' | null
  top_shap: { rf: ShapFeature[]; xgb: ShapFeature[] }
  structure: StructurePayload
  summaries: { rf: Maybe<string>; xgb: Maybe<string> }
  interpretation: string
  integrated?: Maybe<IntegratedReport>
  missing: { field: string; reason: string | null }[]
}

export interface SuggestHit {
  variant_id: string
  gene_symbol: string
  label: Label
  protein_change: string | null
  position: number
}
