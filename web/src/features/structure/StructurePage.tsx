import { useEffect, useRef, useState } from 'react'
import { useParams } from 'react-router-dom'
import { useQuery } from '@tanstack/react-query'
import cssUrl from 'molstar/build/viewer/molstar.css?url'
import jsUrl from 'molstar/build/viewer/molstar.js?url'
import { api, fileUrl } from '../../api/client.ts'
import type { StructurePayload } from '../../api/types.ts'
import { EmptyState, Missing, MockTag, PageIntro, QueryGate } from '../../components/ui.tsx'
import { caseLabel, formatNumber } from '../../lib/format.ts'
import { useTheme } from '../../lib/theme.ts'

type Pane = 'wt' | 'mutant' | 'overlay'

interface MolstarViewer {
  loadStructureFromUrl: (url: string, format?: string, isBinary?: boolean) => Promise<void>
  structureInteractivity?: (args: { action: 'focus'; expression?: (builder: MolBuilder) => unknown }) => void
  dispose: () => void
  handleResize?: () => void
}

interface MolBuilder {
  struct: {
    generator: { atomGroups: (query: Record<string, unknown>) => unknown }
    atomProperty?: { macromolecular?: { auth_asym_id?: () => unknown; auth_seq_id?: () => unknown } }
  }
  core: { rel: { eq: (pair: unknown[]) => unknown } }
}

interface MolstarApi {
  Viewer: {
    create: (element: HTMLElement, options?: Record<string, unknown>) => Promise<MolstarViewer>
  }
}

let molstarPromise: Promise<MolstarApi> | null = null

function loadMolstar() {
  if (!molstarPromise) {
    molstarPromise = new Promise((resolve, reject) => {
      if (!document.querySelector('link[data-molstar-css]')) {
        const link = document.createElement('link')
        link.rel = 'stylesheet'
        link.href = cssUrl
        link.dataset.molstarCss = '1'
        document.head.appendChild(link)
      }
      const existing = document.querySelector('script[data-molstar]')
      const ready = () => {
        const apiGlobal = (window as unknown as { molstar?: MolstarApi }).molstar
        if (!apiGlobal?.Viewer) {
          molstarPromise = null
          reject(new Error('Molstar loaded without a viewer.'))
          return
        }
        resolve(apiGlobal)
      }
      if (existing) {
        ready()
        return
      }
      const script = document.createElement('script')
      script.src = jsUrl
      script.async = true
      script.dataset.molstar = '1'
      script.onload = () => ready()
      script.onerror = () => {
        molstarPromise = null
        reject(new Error('The Molstar viewer script did not load.'))
      }
      document.body.appendChild(script)
    })
  }
  return molstarPromise
}

function viewerUrls(structure: StructurePayload) {
  const wt = fileUrl(structure.files?.wt.url ?? null) ?? structure.mock_viewer?.wt ?? null
  const mutant = fileUrl(structure.files?.mutant.url ?? null) ?? structure.mock_viewer?.mutant ?? null
  return { wt, mutant, mock: Boolean(structure.mock_viewer) }
}

export function StructurePage() {
  const params = useParams()
  const id = params.id ?? ''
  const query = useQuery({
    queryKey: ['structure', id],
    queryFn: () => api.structure(id),
    enabled: id.length > 0,
  })
  return (
    <div className="space-y-5">
      <PageIntro kicker="Structure" title={id ? caseLabel(id) : 'No case'} lede="Wild type, mutant, and an overlay. ΔΔG stays empty until FoldX exists." />
      {!id ? <EmptyState title="Missing case id" /> : <QueryGate query={query}>{(data) => <StructureBody data={data} />}</QueryGate>}
    </div>
  )
}

function StructureBody({ data }: { data: StructurePayload }) {
  const urls = viewerUrls(data)
  const canView = Boolean(urls.wt || urls.mutant)
  return (
    <div className="space-y-4">
      {data.poc ? (
        <p className="text-sm">
          PDB {data.poc.pdb_id}, chain {data.poc.chain}, residue {data.poc.residue} {data.poc.wt_aa}→{data.poc.mut_aa}.{' '}
          {data.poc.source_note}
        </p>
      ) : null}
      {data.method_note ? <p className="border border-dashed border-graphite p-3 text-sm">{data.method_note}</p> : null}
      <div className="grid gap-3 sm:grid-cols-2">
        <Measure label="Cα RMSD (Å)" value={data.rmsd_angstrom.value} reason={data.rmsd_angstrom.reason} mock={data.rmsd_angstrom.mock} />
        <Measure label="ΔΔG" value={data.ddg.value} reason={data.ddg.reason} mock={data.ddg.mock} />
      </div>
      {canView ? (
        <Viewer urls={urls} chain={data.poc?.chain ?? null} residue={data.poc?.residue ?? null} />
      ) : (
        <EmptyState title="No coordinates to load" detail={data.reason ?? 'This variant has no structure files.'} />
      )}
      <p className="text-sm text-mute">{data.esmfold_note}</p>
      <section>
        <h2 className="font-display text-xl">Known structures for this gene</h2>
        {data.known_structures.value == null ? (
          <Missing reason={data.known_structures.reason} />
        ) : data.known_structures.value.length === 0 ? (
          <EmptyState title="No manifest rows" />
        ) : (
          <ul className="mt-2 space-y-2">
            {data.known_structures.value.map((row) => (
              <li key={`${row.structure_id}-${row.source}`} className="border rule p-3 text-sm">
                <p className="font-mono">
                  {row.structure_id ?? 'Not computed'} · {row.source ?? 'source not recorded'}
                </p>
                <p>{row.note}</p>
              </li>
            ))}
          </ul>
        )}
      </section>
    </div>
  )
}

function Measure({
  label,
  value,
  reason,
  mock,
}: {
  label: string
  value: number | null
  reason: string | null
  mock?: boolean
}) {
  return (
    <section className="border rule p-3">
      <h2 className="text-xs tracking-wide text-mute uppercase">
        {label}
        {mock ? <MockTag /> : null}
      </h2>
      {value == null ? <Missing reason={reason} /> : <p className="font-mono text-3xl">{formatNumber(value, 2)}</p>}
    </section>
  )
}

function Viewer({
  urls,
  chain,
  residue,
}: {
  urls: { wt: string | null; mutant: string | null; mock: boolean }
  chain: string | null
  residue: number | null
}) {
  const theme = useTheme()
  const host = useRef<HTMLDivElement>(null)
  const [pane, setPane] = useState<Pane>('wt')
  const target = pane === 'mutant' ? urls.mutant : urls.wt
  const extra = pane === 'overlay' ? urls.mutant : null
  const viewKey = `${pane}|${target ?? ''}|${extra ?? ''}|${theme.theme}`
  const [ticket, setTicket] = useState({ key: '', error: null as string | null, ready: false })
  if (ticket.key !== viewKey) {
    setTicket({ key: viewKey, error: null, ready: false })
  }
  const missing = !target
  const pending = !missing && ticket.key === viewKey && !ticket.ready && !ticket.error
  const error = missing ? 'This pane has no coordinate file.' : ticket.error

  useEffect(() => {
    const node = host.current
    if (!node || !target) return
    let viewer: MolstarViewer | null = null
    let cancelled = false
    void loadMolstar()
      .then(async (molstar) => {
        if (cancelled) return
        viewer = await molstar.Viewer.create(node, {
          layoutIsExpanded: false,
          layoutShowControls: false,
          layoutShowSequence: true,
          layoutShowLog: false,
          layoutShowLeftPanel: false,
          collapseRightPanel: true,
          viewportBackgroundColor: theme.theme === 'dark' ? '#0B0D10' : '#F3EFE6',
        })
        await viewer.loadStructureFromUrl(target, 'pdb', false)
        if (extra) await viewer.loadStructureFromUrl(extra, 'pdb', false)
        if (chain && residue) focusResidue(viewer, chain, residue)
        if (!cancelled) setTicket({ key: viewKey, error: null, ready: true })
      })
      .catch((loadError: unknown) => {
        if (cancelled) return
        setTicket({
          key: viewKey,
          error: loadError instanceof Error ? loadError.message : 'The viewer failed.',
          ready: false,
        })
      })
    return () => {
      cancelled = true
      viewer?.dispose()
    }
  }, [viewKey, target, extra, chain, residue, theme.theme])

  return (
    <section>
      <div className="mb-2 flex flex-wrap gap-2" role="tablist" aria-label="Coordinate pane">
        {(['wt', 'mutant', 'overlay'] as const).map((item) => (
          <button
            key={item}
            type="button"
            role="tab"
            aria-selected={pane === item}
            className={`control ${pane === item ? 'bg-ink text-paper' : ''}`}
            onClick={() => setPane(item)}
          >
            {item === 'wt' ? 'Wild type' : item === 'mutant' ? 'Mutant' : 'Overlay'}
          </button>
        ))}
      </div>
      {urls.mock ? (
        <p className="mb-2 text-sm">
          <MockTag /> Both panes load deposited 9B3R from RCSB. The local mutant PDB is not on disk, so the overlay is
          the same coordinates twice.
        </p>
      ) : null}
      <div className="relative h-[28rem] border rule bg-panel">
        <div ref={host} className="h-full w-full" />
        {pending ? <p className="absolute top-3 left-3 bg-paper px-2 py-1 text-sm">Loading coordinates…</p> : null}
        {error ? (
          <div className="absolute inset-4">
            <EmptyState title="Viewer unavailable" detail={error} />
          </div>
        ) : null}
      </div>
      {chain && residue ? (
        <p className="mt-2 text-sm">
          Focus request: chain {chain}, author residue {residue}.
        </p>
      ) : null}
    </section>
  )
}

function focusResidue(viewer: MolstarViewer, chain: string, residue: number) {
  if (!viewer.structureInteractivity) return
  try {
    viewer.structureInteractivity({
      action: 'focus',
      expression: (builder) => {
        const asym = builder.struct.atomProperty?.macromolecular?.auth_asym_id?.()
        const seq = builder.struct.atomProperty?.macromolecular?.auth_seq_id?.()
        if (!asym || !seq) return undefined
        return builder.struct.generator.atomGroups({
          'chain-test': builder.core.rel.eq([asym, chain]),
          'residue-test': builder.core.rel.eq([seq, residue]),
        })
      },
    })
  } catch {
    /* The residue is still named in the page caption. */
  }
}
