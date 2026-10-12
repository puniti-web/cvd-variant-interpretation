import { Link } from 'react-router-dom'
import { useQuery } from '@tanstack/react-query'
import { api } from '../../api/client.ts'
import type { GeneCount, Overview } from '../../api/types.ts'
import { ConfusionTable, MeanShapBars, RocPr } from '../../charts/Charts.tsx'
import { PipelineHelix } from '../../charts/Helix.tsx'
import { Missing, PageIntro, QueryGate, Readout, Swatch } from '../../components/ui.tsx'
import { caseLabel, formatCount } from '../../lib/format.ts'

export function BenchPage() {
  const overview = useQuery({ queryKey: ['overview'], queryFn: api.overview })
  const metrics = useQuery({ queryKey: ['metrics'], queryFn: api.metrics })
  const shap = useQuery({ queryKey: ['shap'], queryFn: api.globalShap })
  const health = useQuery({ queryKey: ['health'], queryFn: api.health })
  return (
    <div className="space-y-8">
      <PageIntro
        kicker="Bench"
        title="The labelled drawer, then the held-out VUS."
        lede="Eighteen cardiomyopathy and arrhythmia genes. Counts below come from the catalog. A dashed readout means the file is not on this machine."
      />
      <QueryGate query={overview}>
        {(data) => (
          <>
            <SmokeLinks overview={data} />
            <section className="grid gap-4 sm:grid-cols-2 lg:grid-cols-5">
              <Readout label="Collected" maybe={data.counts.collected} />
              <Readout label="Labelled" maybe={data.counts.labelled} />
              <Readout label="Pathogenic" maybe={data.counts.pathogenic} />
              <Readout label="Benign" maybe={data.counts.benign} />
              <Readout label="VUS held out" maybe={data.counts.vus} />
            </section>
            {data.split.caveat ? <p className="text-sm text-mute">{data.split.caveat}</p> : null}
            <section className="grid gap-4 md:grid-cols-2">
              <GeneList title="Train genes" maybe={data.split.train_genes} count={data.split.train_variants} />
              <GeneList title="Test genes" maybe={data.split.test_genes} count={data.split.test_variants} dashed />
            </section>
            <section>
              <h2 className="font-display text-2xl">Pipeline</h2>
              <PipelineHelix pipeline={data.pipeline} reason={data.pipeline_reason} />
            </section>
            <GeneStacks genes={data.genes} />
          </>
        )}
      </QueryGate>
      <section className="space-y-4">
        <h2 className="font-display text-2xl">Test-fold scores</h2>
        <QueryGate query={metrics}>
          {(data) => {
            const live = data.recomputed
            const reported = data.reported
            if (!live && !reported) return <Missing reason={data.reason} />
            const shown = live ?? {
              random_forest: reported!.random_forest,
              xgboost: reported!.xgboost,
            }
            const baseline = reported ?? shown
            const mock = Boolean(live?.mock || live?.random_forest.mock)
            return (
              <div className="space-y-4">
                <p className="text-sm text-mute">
                  {live
                    ? `Split models scored on the held-out genes${data.split.test_genes.value ? ` (${data.split.test_genes.value.join(', ')})` : ''}. The line under each number is the value saved at training time. The cross-validated headline numbers are on the Models page.`
                    : data.reason ?? 'Recomputed metrics were not returned. The cards are the values saved at training time.'}
                </p>
                <div className="grid gap-3 sm:grid-cols-2">
                  <MetricCard
                    name="Random Forest"
                    metrics={shown.random_forest}
                    baseline={baseline.random_forest}
                    mock={mock}
                    matched={data.comparison?.random_forest.all_matched}
                  />
                  <MetricCard
                    name="XGBoost"
                    metrics={shown.xgboost}
                    baseline={baseline.xgboost}
                    mock={mock}
                    matched={data.comparison?.xgboost.all_matched}
                  />
                </div>
                <RocPr roc={live?.random_forest.roc} pr={live?.random_forest.pr} mock={mock} />
                <div className="grid gap-4 lg:grid-cols-2">
                  <ConfusionTable matrix={shown.random_forest.confusion_matrix} label="Random Forest" mock={mock} />
                  <ConfusionTable matrix={shown.xgboost.confusion_matrix} label="XGBoost" mock={mock} />
                </div>
              </div>
            )
          }}
        </QueryGate>
      </section>
      <section className="space-y-3">
        <h2 className="font-display text-2xl">Global SHAP</h2>
        <QueryGate query={shap}>
          {(data) => (
            <>
              <p className="text-sm text-mute">{data.note}</p>
              <p className="text-sm">
                Rows:{' '}
                {data.row_count.value == null ? <Missing reason={data.row_count.reason} /> : formatCount(data.row_count.value)}
                . Imputed conservation rows:{' '}
                {data.imputed_conservation_rows.value == null ? (
                  <Missing reason={data.imputed_conservation_rows.reason} />
                ) : (
                  formatCount(data.imputed_conservation_rows.value)
                )}
                {data.imputed_conservation_rows.mock ? ' (mock)' : ''}.
              </p>
              <MeanShapBars features={data.xgboost.value} mock={data.xgboost.mock || data.mock} label="XGBoost" />
            </>
          )}
        </QueryGate>
      </section>
      <section>
        <h2 className="font-display text-2xl">Artifacts on disk</h2>
        <QueryGate query={health}>
          {(data) => (
            <ul className="mt-3 grid gap-2 sm:grid-cols-2 lg:grid-cols-3">
              {data.artifacts.map((artifact) => (
                <li
                  key={artifact.key}
                  className={`p-3 text-sm ${artifact.loaded ? 'border border-ink' : 'border border-dashed border-graphite'}`}
                >
                  <p className="font-mono text-xs">{artifact.key}</p>
                  <p>{artifact.loaded ? 'Loaded' : artifact.present ? 'Present, not loaded' : 'Absent'}</p>
                  {artifact.reason ? <p className="text-xs text-mute">{artifact.reason}</p> : null}
                </li>
              ))}
            </ul>
          )}
        </QueryGate>
      </section>
    </div>
  )
}

function SmokeLinks({ overview }: { overview: Overview }) {
  const smoke = overview.smoke
  if (!smoke) return null
  const links = [
    ['Structure proof of concept', smoke.structure],
    ['Test gene, no structure', smoke.no_structure],
    ['Train gene', smoke.train_gene],
  ] as const
  return (
    <ul className="flex flex-wrap gap-2">
      {links.map(([label, id]) => (
        <li key={id}>
          <Link className="control inline-flex items-center font-mono text-sm" to={`/explore?variant=${id}`}>
            {label}: {caseLabel(id)}
          </Link>
        </li>
      ))}
    </ul>
  )
}

function GeneList({
  title,
  maybe,
  count,
  dashed,
}: {
  title: string
  maybe: Overview['split']['train_genes']
  count: Overview['split']['train_variants']
  dashed?: boolean
}) {
  return (
    <section className={`p-3 ${dashed ? 'border border-dashed border-graphite' : 'border border-ink'}`}>
      <h2 className="font-display text-xl">{title}</h2>
      {count.value == null ? <Missing reason={count.reason} /> : <p className="font-mono text-sm">{formatCount(count.value)} variants</p>}
      {maybe.value == null ? (
        <Missing reason={maybe.reason} />
      ) : (
        <ul className="mt-2 flex flex-wrap gap-2">
          {maybe.value.map((gene) => (
            <li key={gene}>
              <Link className="inline-flex min-h-11 items-center font-mono text-sm underline" to={`/explore?gene=${gene}`}>
                {gene}
              </Link>
            </li>
          ))}
        </ul>
      )}
    </section>
  )
}

function GeneStacks({ genes }: { genes: GeneCount[] }) {
  const totals = genes.map((gene) => (gene.pathogenic.value ?? 0) + (gene.benign.value ?? 0) + (gene.vus.value ?? 0))
  const max = Math.max(...totals, 1)
  const panels = ['cardiomyopathy', 'arrhythmia'] as const
  return (
    <section className="space-y-4">
      <div className="flex flex-wrap gap-4">
        <Swatch tone="pathogenic">Pathogenic</Swatch>
        <Swatch tone="benign">Benign</Swatch>
        <Swatch tone="vus">VUS</Swatch>
      </div>
      {panels.map((panel) => (
        <div key={panel}>
          <h2 className="font-display text-2xl capitalize">{panel}</h2>
          <ul className="mt-2 space-y-2">
            {genes
              .filter((gene) => gene.panel === panel)
              .map((gene) => {
                const path = gene.pathogenic.value ?? 0
                const benign = gene.benign.value ?? 0
                const vus = gene.vus.value ?? 0
                const width = ((path + benign + vus) / max) * 100
                const missing = gene.pathogenic.value == null || gene.benign.value == null || gene.vus.value == null
                return (
                  <li key={gene.symbol} className="grid grid-cols-[5.5rem_1fr] items-center gap-2">
                    <Link className="font-mono text-sm underline" to={`/explore?gene=${gene.symbol}`}>
                      {gene.symbol}
                    </Link>
                    {missing ? (
                      <Missing reason="A class count is missing for this gene." />
                    ) : (
                      <span className="flex min-w-0 items-center gap-2">
                        <span
                          className="flex h-6 min-w-0 flex-1"
                          aria-hidden
                        >
                          <span className="flex h-6" style={{ width: `${width}%` }}>
                            <span className="bg-pathogenic" style={{ width: `${(path / (path + benign + vus || 1)) * 100}%` }} />
                            <span className="bg-benign" style={{ width: `${(benign / (path + benign + vus || 1)) * 100}%` }} />
                            <span className="bg-vus" style={{ width: `${(vus / (path + benign + vus || 1)) * 100}%` }} />
                          </span>
                        </span>
                        <span className="shrink-0 font-mono text-[11px]">
                          P {path} · B {benign} · V {vus}
                        </span>
                      </span>
                    )}
                  </li>
                )
              })}
          </ul>
        </div>
      ))}
    </section>
  )
}

function MetricCard({
  name,
  metrics,
  baseline,
  mock,
  matched,
}: {
  name: string
  metrics: { accuracy: number; f1: number; roc_auc: number | null; pr_auc: number | null }
  baseline: { accuracy: number; f1: number; roc_auc: number | null; pr_auc: number | null }
  mock: boolean
  matched?: boolean
}) {
  const rows = [
    ['Accuracy', metrics.accuracy, baseline.accuracy],
    ['F1', metrics.f1, baseline.f1],
    ['ROC AUC', metrics.roc_auc, baseline.roc_auc],
    ['PR AUC', metrics.pr_auc, baseline.pr_auc],
  ] as const
  return (
    <section className="border rule p-3">
      <h3 className="font-display text-xl">
        {name}
        {mock ? <span className="ml-2 font-mono text-xs tracking-widest uppercase">Mock</span> : null}
      </h3>
      {matched != null ? (
        <p className="text-xs text-mute">{matched ? 'Matches the saved training-time values at 3 decimals.' : 'Does not match the saved training-time values.'}</p>
      ) : null}
      <dl className="mt-2 space-y-1">
        {rows.map(([label, value, base]) => (
          <div key={label} className="flex items-baseline justify-between gap-3">
            <dt>{label}</dt>
            <dd className="text-right font-mono">
              {value == null ? 'Not computed' : value.toFixed(3)}
              <span className="ml-2 text-xs text-mute">saved {base == null ? '—' : base.toFixed(3)}</span>
            </dd>
          </div>
        ))}
      </dl>
    </section>
  )
}
