import { useQuery } from '@tanstack/react-query'
import { useSearchParams } from 'react-router-dom'
import { api } from '../../api/client.ts'
import type { Page, VariantSummary, VusRanking, VusSummary } from '../../api/types.ts'
import { CaseFile } from '../../components/CaseFile.tsx'
import { ClassMark, EmptyState, PageIntro, QueryGate } from '../../components/ui.tsx'
import { caseLabel, splitWord } from '../../lib/format.ts'

export function ExplorePage() {
  const [params, setParams] = useSearchParams()
  const tab = params.get('tab') === 'vus' ? 'vus' : 'labelled'
  const q = params.get('q') ?? ''
  const gene = params.get('gene') ?? ''
  const label = params.get('label') ?? ''
  const split = params.get('split') ?? ''
  const page = Math.max(1, Number(params.get('page') ?? '1') || 1)
  const variantId = params.get('variant') ?? ''
  const genes = useQuery({ queryKey: ['overview'], queryFn: api.overview })
  const labelled = useQuery({
    queryKey: ['labelled', q, gene, label, split, page],
    queryFn: () => api.variants({ q, gene, label, split, page, page_size: 25 }),
    enabled: tab === 'labelled',
  })
  const held = useQuery({
    queryKey: ['vus', q, gene, page],
    queryFn: () => api.vus({ q, gene, page, page_size: 25 }),
    enabled: tab === 'vus',
  })
  const detail = useQuery({
    queryKey: ['variant', variantId],
    queryFn: () => api.variant(variantId),
    enabled: variantId.length > 0,
  })

  function update(next: Record<string, string | null>, resetPage = false) {
    const following = new URLSearchParams(params)
    for (const [key, value] of Object.entries(next)) {
      if (!value) following.delete(key)
      else following.set(key, value)
    }
    if (resetPage) following.delete('page')
    setParams(following)
  }

  return (
    <div className="space-y-5">
      <PageIntro
        kicker="Explorer"
        title="Open a case file."
        lede="Labelled variants sit in the training split. VUS stay held out and are listed by model support for pathogenicity, strongest first. Filters combine."
      />
      <div className="flex flex-wrap gap-2" role="tablist" aria-label="Catalog">
        <button
          type="button"
          role="tab"
          aria-selected={tab === 'labelled'}
          className={`control ${tab === 'labelled' ? 'bg-ink text-paper' : ''}`}
          onClick={() => update({ tab: null }, true)}
        >
          Labelled
        </button>
        <button
          type="button"
          role="tab"
          aria-selected={tab === 'vus'}
          className={`control ${tab === 'vus' ? 'bg-ink text-paper' : ''}`}
          onClick={() => update({ tab: 'vus', label: null, split: null }, true)}
        >
          VUS
        </button>
      </div>
      <form
        className="grid gap-2 md:grid-cols-4"
        onSubmit={(event) => {
          event.preventDefault()
          const data = new FormData(event.currentTarget)
          update(
            {
              q: String(data.get('q') ?? ''),
              gene: String(data.get('gene') ?? ''),
              label: String(data.get('label') ?? ''),
              split: String(data.get('split') ?? ''),
            },
            true,
          )
        }}
      >
        <input className="control font-mono" name="q" defaultValue={q} aria-label="Search" placeholder="Search id, gene, change, title" key={q} />
        <select
          className="control"
          name="gene"
          defaultValue={gene}
          aria-label="Gene"
          key={`gene-${gene}`}
          onChange={(event) => update({ gene: event.target.value }, true)}
        >
          <option value="">All genes</option>
          {genes.data?.genes.map((item) => (
            <option key={item.symbol} value={item.symbol}>
              {item.symbol}
            </option>
          ))}
        </select>
        {tab === 'labelled' ? (
          <>
            <select
              className="control"
              name="label"
              defaultValue={label}
              aria-label="Clinical class"
              key={`label-${label}`}
              onChange={(event) => update({ label: event.target.value }, true)}
            >
              <option value="">Any class</option>
              <option value="pathogenic">Pathogenic</option>
              <option value="benign">Benign</option>
            </select>
            <select
              className="control"
              name="split"
              defaultValue={split}
              aria-label="Fold"
              key={`split-${split}`}
              onChange={(event) => update({ split: event.target.value }, true)}
            >
              <option value="">Any fold</option>
              <option value="train_gene">Train gene</option>
              <option value="test_gene">Test gene</option>
            </select>
          </>
        ) : null}
        <button type="submit" className="control md:col-span-4">
          Apply search
        </button>
      </form>
      <div className="grid gap-4 lg:grid-cols-[minmax(0,1fr)_minmax(0,1.1fr)]">
        <div>
          {tab === 'vus' ? (
            <QueryGate query={held}>
              {(data) => <ResultList data={data} vus variantId={variantId} page={page} update={update} />}
            </QueryGate>
          ) : (
            <QueryGate query={labelled}>
              {(data) => <ResultList data={data} vus={false} variantId={variantId} page={page} update={update} />}
            </QueryGate>
          )}
        </div>
        <div>
          {!variantId ? <EmptyState title="Select a case" detail="The file opens beside the list." /> : null}
          {variantId ? (
            <QueryGate query={detail}>{(row) => <CaseFile detail={row} />}</QueryGate>
          ) : null}
        </div>
      </div>
    </div>
  )
}

function ResultList({
  data,
  vus,
  variantId,
  page,
  update,
}: {
  data: Page<VariantSummary> | Page<VusSummary>
  vus: boolean
  variantId: string
  page: number
  update: (next: Record<string, string | null>, resetPage?: boolean) => void
}) {
  const pageCount = Math.max(1, Math.ceil(data.total / data.page_size))
  if (data.items.length === 0) {
    return <EmptyState title="No variants in this drawer" detail={data.reason ?? 'Nothing matches these filters.'} />
  }
  return (
    <>
      <p className="mb-2 font-mono text-xs text-mute">
        {data.total} matches · page {data.page} of {pageCount}
      </p>
      <ul className="space-y-2">
        {data.items.map((item) => (
          <li key={item.variant_id}>
            <button
              type="button"
              className={`w-full border p-3 text-left ${item.variant_id === variantId ? 'border-ink bg-panel' : 'rule'}`}
              aria-current={item.variant_id === variantId ? 'true' : undefined}
              onClick={() => update({ variant: item.variant_id })}
            >
              <Card item={item} vus={vus} />
            </button>
          </li>
        ))}
      </ul>
      <div className="mt-3 flex gap-2">
        <button
          type="button"
          className="control disabled:opacity-40"
          disabled={page <= 1}
          onClick={() => update({ page: page <= 2 ? null : String(page - 1) })}
        >
          Previous
        </button>
        <button
          type="button"
          className="control disabled:opacity-40"
          disabled={page >= pageCount}
          onClick={() => update({ page: String(page + 1) })}
        >
          Next
        </button>
      </div>
    </>
  )
}

function RankingLine({ ranking }: { ranking: VusRanking }) {
  const ddg =
    ranking.foldx_ddg == null
      ? 'no ΔΔG'
      : `ΔΔG ${ranking.foldx_ddg > 0 ? '+' : ''}${ranking.foldx_ddg.toFixed(2)} (${ranking.ddg_category})`
  return (
    <span className="mt-1 block text-xs">
      <span className="font-mono">#{ranking.rank}</span> · {ranking.tier}
      {ranking.probability != null ? ` · P ${ranking.probability.toFixed(2)}` : ''} · {ddg}
      {ranking.evidence_informative ? ` · evidence ${ranking.evidence_agree}/${ranking.evidence_informative}` : ''}
    </span>
  )
}

function Card({ item, vus }: { item: VariantSummary | VusSummary; vus: boolean }) {
  return (
    <span className="flex flex-wrap items-start justify-between gap-2">
      <span>
        <span className="block font-mono text-sm">{caseLabel(item.variant_id)}</span>
        <span className="block text-sm text-mute">
          {item.gene_symbol} {item.protein_change ?? `${item.wt_aa}${item.position}${item.mut_aa}`} · {item.clinical_significance}
        </span>
        {vus && 'predictions_reason' in item && item.predictions == null ? (
          <span className="mt-1 block text-xs text-mute">Scores not computed. {item.predictions_reason}</span>
        ) : null}
        {vus && 'ranking' in item && item.ranking ? <RankingLine ranking={item.ranking} /> : null}
      </span>
      <span className="flex flex-col items-end gap-1 text-sm">
        <ClassMark label={item.label} />
        {!vus && 'split' in item ? <span>{item.split ? splitWord(item.split) : 'Split not computed'}</span> : null}
        {!vus && 'has_structure' in item && item.has_structure ? <span>AlphaFold model</span> : null}
      </span>
    </span>
  )
}
