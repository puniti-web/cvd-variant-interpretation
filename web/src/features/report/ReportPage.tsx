import { useState } from 'react'
import { Link, useParams } from 'react-router-dom'
import { useQuery } from '@tanstack/react-query'
import { api } from '../../api/client.ts'
import type { IntegratedReport, Report, ShapFeature } from '../../api/types.ts'
import { ModelDials } from '../../charts/Dials.tsx'
import { PhyschemBars } from '../../charts/Physchem.tsx'
import { EmptyState, Missing, MockTag, PageIntro, QueryGate, Stamp } from '../../components/ui.tsx'
import { caseLabel, featureLabel, formatNumber, formatProb } from '../../lib/format.ts'

export function ReportPage() {
  const params = useParams()
  const id = params.id ?? ''
  const query = useQuery({
    queryKey: ['report', id],
    queryFn: async () => {
      const body = await api.report(id, 'json')
      if (typeof body === 'string') throw new Error('The report endpoint returned markdown for a JSON request.')
      return body
    },
    enabled: id.length > 0,
  })
  return (
    <div className="space-y-5">
      <PageIntro kicker="Report" title={id ? caseLabel(id) : 'No case'} />
      {!id ? <EmptyState title="Missing case id" /> : <QueryGate query={query}>{(report) => <Dossier report={report} />}</QueryGate>}
    </div>
  )
}

function Dossier({ report }: { report: Report }) {
  const [note, setNote] = useState<string | null>(null)
  async function download(format: 'json' | 'markdown') {
    setNote(null)
    try {
      const body = await api.report(report.variant_id, format)
      const text = typeof body === 'string' ? body : JSON.stringify(body, null, 2)
      const blob = new Blob([text], { type: format === 'json' ? 'application/json' : 'text/markdown' })
      const url = URL.createObjectURL(blob)
      const anchor = document.createElement('a')
      anchor.href = url
      anchor.download = `${report.variant_id}.${format === 'json' ? 'json' : 'md'}`
      anchor.click()
      URL.revokeObjectURL(url)
    } catch (error) {
      setNote(error instanceof Error ? error.message : 'Download failed')
    }
  }
  const clinicalLabel =
    report.clinical.clinical_significance?.toLowerCase().includes('pathogenic') &&
    !report.clinical.clinical_significance.toLowerCase().includes('benign')
      ? 'pathogenic'
      : report.clinical.clinical_significance?.toLowerCase().includes('benign')
        ? 'benign'
        : 'vus'
  return (
    <article className="space-y-5">
      <div className="no-print flex flex-wrap gap-2">
        <button type="button" className="control" onClick={() => void download('markdown')}>
          Download Markdown
        </button>
        <button type="button" className="control" onClick={() => void download('json')}>
          Download JSON
        </button>
        <button type="button" className="control" onClick={() => window.print()}>
          Print
        </button>
        <Link className="control inline-flex items-center" to={`/explore?variant=${report.variant_id}`}>
          Case file
        </Link>
      </div>
      {note ? <p role="alert">{note}</p> : null}
      <Stamp variantId={report.variant_id} label={clinicalLabel} />
      <p className="border border-ink p-3 font-display text-xl">{report.disclaimer}</p>
      {report.integrated?.value ? (
        <IntegratedPanel report={report.integrated.value} />
      ) : report.integrated ? (
        <section>
          <h2 className="font-display text-xl">Integrated interpretation</h2>
          <Missing reason={report.integrated.reason} />
        </section>
      ) : null}
      <details>
        <summary className="cursor-pointer font-display text-xl">Model-by-model notes</summary>
        <p className="mt-2 whitespace-pre-line">{report.interpretation}</p>
      </details>
      <ModelDials rf={report.predictions.rf} xgb={report.predictions.xgb} agreement={report.agreement} />
      <section className="grid gap-3 md:grid-cols-2">
        <Summary title="Random Forest summary" text={report.summaries.rf} />
        <Summary title="XGBoost summary" text={report.summaries.xgb} />
      </section>
      <div className="grid gap-4 lg:grid-cols-2">
        <TopList title="Random Forest top SHAP" rows={report.top_shap.rf} />
        <TopList title="XGBoost top SHAP" rows={report.top_shap.xgb} />
      </div>
      <section>
        <h2 className="font-display text-xl">
          Conservation
          {report.conservation.score.mock ? <MockTag /> : null}
        </h2>
        {report.conservation.score.value == null ? (
          <Missing reason={report.conservation.score.reason} />
        ) : (
          <p className="font-mono text-2xl">{formatProb(report.conservation.score.value)}</p>
        )}
        <p className="text-sm text-mute">{report.conservation.method}</p>
        <p className="text-sm">{report.conservation.imputed ? 'Score was imputed.' : 'Score was not imputed.'}</p>
      </section>
      <PhyschemBars physchem={report.physchem} />
      <section>
        <h2 className="font-display text-xl">Structure</h2>
        {report.structure.method_note ? <p>{report.structure.method_note}</p> : <p>{report.structure.reason ?? 'Structure was not run.'}</p>}
        <p className="mt-2 text-sm">
          RMSD:{' '}
          {report.structure.rmsd_angstrom.value == null ? (
            <Missing reason={report.structure.rmsd_angstrom.reason} />
          ) : (
            <span className="font-mono">
              {formatProb(report.structure.rmsd_angstrom.value)} Å{report.structure.rmsd_angstrom.mock ? ' (mock)' : ''}
            </span>
          )}
        </p>
        <p className="mt-2 text-sm">
          FoldX ΔΔG:{' '}
          {report.structure.ddg.value == null ? (
            <Missing reason={report.structure.ddg.reason} />
          ) : (
            <span className="font-mono">
              {report.structure.ddg.value > 0 ? '+' : ''}
              {formatNumber(report.structure.ddg.value, 2)} kcal/mol ({report.structure.ddg_category ?? 'uncategorised'})
            </span>
          )}
        </p>
      </section>
      <section>
        <h2 className="font-display text-xl">Clinical record</h2>
        <dl className="grid gap-2 text-sm sm:grid-cols-2">
          <div>
            <dt className="text-mute">Significance</dt>
            <dd>{report.clinical.clinical_significance ?? 'Not computed'}</dd>
          </div>
          <div>
            <dt className="text-mute">Review</dt>
            <dd>{report.clinical.review_status ?? 'Not computed'}</dd>
          </div>
          <div>
            <dt className="text-mute">Title</dt>
            <dd>{report.clinical.title ?? 'Not computed'}</dd>
          </div>
          <div>
            <dt className="text-mute">ClinVar UID</dt>
            <dd className="font-mono">{report.clinical.clinvar_uid ?? 'Not computed'}</dd>
          </div>
        </dl>
      </section>
      <section>
        <h2 className="font-display text-xl">Not computed</h2>
        {report.missing.length === 0 ? (
          <p>Every dossier field in this rendering has a value.</p>
        ) : (
          <ul className="space-y-2">
            {report.missing.map((item) => (
              <li key={item.field}>
                <span className="font-mono text-sm">{item.field}</span>
                <Missing reason={item.reason} />
              </li>
            ))}
          </ul>
        )}
      </section>
    </article>
  )
}

const DIRECTION_TONE: Record<string, string> = {
  pathogenic: 'bg-pathogenic',
  benign: 'bg-benign',
  neutral: 'bg-graphite',
  unavailable: 'bg-graphite',
}

function IntegratedPanel({ report }: { report: IntegratedReport }) {
  const p = report.prediction
  return (
    <section className="space-y-4 border border-ink p-4">
      <div>
        <p className="font-mono text-xs tracking-[0.18em] text-mute uppercase">Integrated interpretation</p>
        <p className="mt-1 font-display text-2xl">{report.conclusion}</p>
        <p className="mt-1 text-sm text-mute">
          P(pathogenic) {formatProb(p.probability)} from {p.source}
          {p.calibrated_at_training_prevalence != null
            ? `; ${formatProb(p.calibrated_at_training_prevalence)} at the training prevalence`
            : ''}
          . {report.condition ? `Disease context: ${report.condition}.` : ''}
        </p>
      </div>
      <div className="overflow-x-auto">
        <table className="w-full text-sm">
          <thead>
            <tr className="text-left text-mute">
              <th className="py-1 pr-4 font-normal">Evidence line</th>
              <th className="py-1 pr-4 font-normal">Points to</th>
              <th className="py-1 pr-4 font-normal">Detail</th>
            </tr>
          </thead>
          <tbody>
            {report.evidence.map((line) => (
              <tr key={line.name} className="border-t rule align-top">
                <td className="py-1 pr-4">
                  {line.name}
                  {line.independent ? null : <span className="block text-xs text-mute">reference only</span>}
                </td>
                <td className="py-1 pr-4 whitespace-nowrap">
                  <span className="inline-flex items-center gap-2">
                    <span className={`inline-block h-3 w-3 ${DIRECTION_TONE[line.direction] ?? 'bg-graphite'}`} aria-hidden />
                    {line.direction}
                  </span>
                </td>
                <td className="py-1 pr-4">{line.detail}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      <div>
        <h3 className="font-display text-lg">Why the model decided this</h3>
        <ol className="mt-1 space-y-1 text-sm">
          {report.key_features.map((item) => (
            <li key={item.feature} className="flex flex-wrap justify-between gap-3">
              <span>
                {item.label} <span className="font-mono text-mute">= {item.value_display}</span>
              </span>
              <span className="font-mono">
                {item.shap >= 0 ? '+' : '−'}
                {formatProb(Math.abs(item.shap))} towards {item.direction}
              </span>
            </li>
          ))}
        </ol>
        <p className="mt-1 text-xs text-mute">SHAP values of the full XGBoost model, in log-odds.</p>
      </div>
      {report.flags.length > 0 ? (
        <div>
          <h3 className="font-display text-lg">Flags</h3>
          <ul className="mt-1 list-disc space-y-1 pl-5 text-sm">
            {report.flags.map((flag) => (
              <li key={flag}>{flag}</li>
            ))}
          </ul>
        </div>
      ) : null}
      <p className="text-xs text-mute">{report.disclaimer}</p>
    </section>
  )
}

function Summary({ title, text }: { title: string; text: Report['summaries']['rf'] }) {
  return (
    <section className="border rule p-3">
      <h2 className="font-display text-xl">
        {title}
        {text.mock ? <MockTag /> : null}
      </h2>
      {text.value == null ? <Missing reason={text.reason} /> : <p className="mt-2 text-sm">{text.value}</p>}
    </section>
  )
}

function TopList({ title, rows }: { title: string; rows: ShapFeature[] }) {
  return (
    <section className="border rule p-3">
      <h2 className="font-display text-xl">{title}</h2>
      {rows.length === 0 ? (
        <Missing reason="SHAP was not computed for this model." />
      ) : (
        <ol className="mt-2 space-y-1 font-mono text-sm">
          {rows.map((row) => (
            <li key={row.feature} className="flex justify-between gap-3">
              <span>{featureLabel(row.feature)}</span>
              <span>
                {row.shap_value >= 0 ? '+' : '−'}
                {formatProb(Math.abs(row.shap_value))}
              </span>
            </li>
          ))}
        </ol>
      )}
    </section>
  )
}
