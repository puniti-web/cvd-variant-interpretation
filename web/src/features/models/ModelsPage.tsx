import { useQuery } from '@tanstack/react-query'
import { api } from '../../api/client.ts'
import type { Ablation, CrossValidation, CvStat } from '../../api/types.ts'
import { EmptyState, PageIntro, QueryGate } from '../../components/ui.tsx'
import { formatCount } from '../../lib/format.ts'

const MODEL_NAMES: Record<string, string> = {
  logistic_regression: 'Logistic regression',
  random_forest: 'Random forest',
  xgboost: 'XGBoost',
}

const COLUMNS: [string, string][] = [
  ['roc_auc', 'ROC AUC'],
  ['pr_auc', 'PR AUC'],
  ['mcc', 'MCC'],
  ['balanced_accuracy', 'Balanced acc.'],
  ['recall', 'Recall'],
  ['specificity', 'Specificity'],
  ['brier', 'Brier'],
]

export function ModelsPage() {
  const metrics = useQuery({ queryKey: ['metrics'], queryFn: api.metrics })
  const ablation = useQuery({ queryKey: ['ablation'], queryFn: api.ablation })
  return (
    <div className="space-y-8">
      <PageIntro
        kicker="Models"
        title="How well the classifier generalises."
        lede="Headline numbers hold out a whole gene at a time, so every score comes from a model that never saw a variant of that gene. Brackets are bootstrap 95% intervals over variants."
      />
      <QueryGate query={metrics}>
        {(data) =>
          data.cross_validation ? (
            <CvSections cv={data.cross_validation} />
          ) : (
            <EmptyState
              title="No cross-validation results"
              detail="outputs/metrics.json is missing or predates the cross-validation block. Run python -m src.modeling.train_classifier."
            />
          )
        }
      </QueryGate>
      <section className="space-y-3">
        <h2 className="font-display text-2xl">Feature-group ablation</h2>
        <QueryGate query={ablation}>{(data) => <AblationTable data={data} />}</QueryGate>
      </section>
    </div>
  )
}

function CvSections({ cv }: { cv: CrossValidation }) {
  return (
    <>
      <section className="space-y-3">
        <h2 className="font-display text-2xl">Leave one gene out</h2>
        <p className="text-sm text-mute">
          {cv.n_variants != null ? `${formatCount(cv.n_variants)} labelled variants` : 'Labelled variants'}
          {cv.n_genes != null ? ` across ${cv.n_genes} genes` : ''}. Threshold 0.5.
          {cv.best_model ? ` Best model: ${MODEL_NAMES[cv.best_model] ?? cv.best_model}.` : ''}
        </p>
        <MetricTable block={cv.gene_grouped} />
      </section>
      <section className="space-y-3">
        <h2 className="font-display text-2xl">Within gene, for comparison</h2>
        <p className="text-sm text-mute">
          Stratified 5-fold, 5 repeats. Genes are shared between training and test folds, so these scores are
          optimistic for a gene the model has not seen.
        </p>
        <MetricTable block={cv.within_gene} />
      </section>
      {cv.baselines ? (
        <section className="space-y-3">
          <h2 className="font-display text-2xl">Published predictors on the same variants</h2>
          <p className="text-sm text-mute">
            REVEL and AlphaMissense used ClinVar labels in training or calibration, so their scores here are
            optimistic.
          </p>
          <div className="overflow-x-auto">
            <table className="w-full text-sm">
              <thead>
                <tr className="text-left text-mute">
                  <th className="py-1 pr-4 font-normal">Predictor</th>
                  <th className="py-1 pr-4 font-normal">Coverage</th>
                  <th className="py-1 pr-4 font-normal">ROC AUC</th>
                  <th className="py-1 pr-4 font-normal">Our model, same variants</th>
                </tr>
              </thead>
              <tbody className="font-mono">
                {Object.entries(cv.baselines).map(([name, row]) => (
                  <tr key={name} className="border-t rule">
                    <td className="py-1 pr-4">{name}</td>
                    <td className="py-1 pr-4">{Math.round(row.coverage * 100)}%</td>
                    <td className="py-1 pr-4">{row.roc_auc.toFixed(3)}</td>
                    <td className="py-1 pr-4">{row.model_roc_auc_same_variants.toFixed(3)}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </section>
      ) : null}
      {cv.per_gene ? (
        <section className="space-y-3">
          <h2 className="font-display text-2xl">Per held-out gene</h2>
          <p className="text-sm text-mute">ROC AUC is undefined for genes with only one class.</p>
          <div className="overflow-x-auto">
            <table className="w-full text-sm">
              <thead>
                <tr className="text-left text-mute">
                  <th className="py-1 pr-4 font-normal">Gene</th>
                  <th className="py-1 pr-4 font-normal">Variants</th>
                  <th className="py-1 pr-4 font-normal">Pathogenic / benign</th>
                  <th className="py-1 pr-4 font-normal">ROC AUC</th>
                  <th className="py-1 pr-4 font-normal">Accuracy</th>
                </tr>
              </thead>
              <tbody className="font-mono">
                {cv.per_gene.map((row) => (
                  <tr key={row.gene_symbol} className="border-t rule">
                    <td className="py-1 pr-4">{row.gene_symbol}</td>
                    <td className="py-1 pr-4">{row.n}</td>
                    <td className="py-1 pr-4">
                      {row.pathogenic} / {row.benign}
                    </td>
                    <td className="py-1 pr-4">{row.roc_auc == null ? '—' : row.roc_auc.toFixed(3)}</td>
                    <td className="py-1 pr-4">{row.accuracy.toFixed(3)}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </section>
      ) : null}
    </>
  )
}

function cell(stat: CvStat | undefined, withCi: boolean) {
  if (!stat || stat.mean == null) return '—'
  const value = stat.mean.toFixed(3)
  if (!withCi || !stat.ci95) return value
  return `${value} [${stat.ci95[0].toFixed(3)}–${stat.ci95[1].toFixed(3)}]`
}

function MetricTable({ block }: { block: Record<string, Record<string, CvStat>> }) {
  return (
    <div className="overflow-x-auto">
      <table className="w-full text-sm">
        <thead>
          <tr className="text-left text-mute">
            <th className="py-1 pr-4 font-normal">Model</th>
            {COLUMNS.map(([key, label]) => (
              <th key={key} className="py-1 pr-4 font-normal">
                {label}
              </th>
            ))}
          </tr>
        </thead>
        <tbody className="font-mono">
          {Object.entries(block).map(([model, stats]) => (
            <tr key={model} className="border-t rule">
              <td className="py-1 pr-4 font-sans">{MODEL_NAMES[model] ?? model}</td>
              {COLUMNS.map(([key]) => (
                <td key={key} className="py-1 pr-4 whitespace-nowrap">
                  {cell(stats[key], key === 'roc_auc')}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  )
}

interface AblationValue {
  removed_groups: Record<string, string[]>
  models: Record<string, Record<string, Record<string, { mean: number }>>>
}

function isAblationValue(value: unknown): value is AblationValue {
  return Boolean(value && typeof value === 'object' && 'models' in value && 'removed_groups' in value)
}

function AblationTable({ data }: { data: Ablation }) {
  if (!isAblationValue(data.value)) {
    return <EmptyState title="No ablation results" detail={data.reason ?? 'outputs/ablation_results.json has an unexpected shape.'} />
  }
  const value = data.value
  const models = Object.keys(value.models)
  const settings = Object.keys(value.models[models[0]] ?? {})
  return (
    <>
      <p className="text-sm text-mute">{data.plan.note}</p>
      <div className="overflow-x-auto">
        <table className="w-full text-sm">
          <thead>
            <tr className="text-left text-mute">
              <th className="py-1 pr-4 font-normal">Feature set</th>
              {models.map((model) => (
                <th key={model} className="py-1 pr-4 font-normal">
                  {MODEL_NAMES[model] ?? model} ROC AUC
                </th>
              ))}
            </tr>
          </thead>
          <tbody className="font-mono">
            {settings.map((setting) => {
              const group = setting.replace(/^without_/, '')
              const removed = value.removed_groups[group]
              return (
                <tr key={setting} className="border-t rule">
                  <td className="py-1 pr-4 font-sans">
                    {setting.replaceAll('_', ' ')}
                    {removed ? <span className="block text-xs text-mute">{removed.join(', ')}</span> : null}
                  </td>
                  {models.map((model) => {
                    const score = value.models[model][setting]?.roc_auc?.mean
                    return (
                      <td key={model} className="py-1 pr-4">
                        {score == null ? '—' : score.toFixed(3)}
                      </td>
                    )
                  })}
                </tr>
              )
            })}
          </tbody>
        </table>
      </div>
    </>
  )
}
