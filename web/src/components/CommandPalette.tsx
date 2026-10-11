import { useEffect, useId, useState } from 'react'
import { useNavigate } from 'react-router-dom'
import { useQuery } from '@tanstack/react-query'
import { api } from '../api/client.ts'
import { ClassMark } from './ui.tsx'
import { caseLabel } from '../lib/format.ts'

export function CommandPalette({ onClose }: { onClose: () => void }) {
  const navigate = useNavigate()
  const titleId = useId()
  const [text, setText] = useState('')
  const [debounced, setDebounced] = useState('')
  const query = useQuery({
    queryKey: ['suggest', debounced],
    queryFn: () => api.suggest(debounced),
  })
  const hits = query.data ?? []
  const [cursor, setCursor] = useState({ q: '', index: 0 })
  const active = cursor.q === debounced ? cursor.index : 0

  useEffect(() => {
    const timer = window.setTimeout(() => setDebounced(text), 150)
    return () => window.clearTimeout(timer)
  }, [text])

  function choose(index: number) {
    setCursor({ q: debounced, index })
  }

  function openHit(id: string) {
    const params = new URLSearchParams({ variant: id })
    void navigate(`/explore?${params.toString()}`)
    onClose()
  }

  return (
    <div className="no-print fixed inset-0 z-50 flex items-start justify-center bg-ink/40 p-4" onMouseDown={onClose}>
      <div
        role="dialog"
        aria-modal="true"
        aria-labelledby={titleId}
        className="mt-16 w-full max-w-xl border border-ink bg-paper text-ink shadow-none"
        onMouseDown={(event) => event.stopPropagation()}
      >
        <h2 id={titleId} className="sr-only">
          Find a case
        </h2>
        <input
          autoFocus
          className="control w-full border-0 border-b font-mono"
          placeholder="Gene, protein change, or case id"
          value={text}
          aria-label="Find a case"
          onChange={(event) => setText(event.target.value)}
          onKeyDown={(event) => {
            if (event.key === 'Escape') onClose()
            if (event.key === 'ArrowDown') {
              event.preventDefault()
              setCursor({ q: debounced, index: Math.min(active + 1, Math.max(hits.length - 1, 0)) })
            }
            if (event.key === 'ArrowUp') {
              event.preventDefault()
              setCursor({ q: debounced, index: Math.max(active - 1, 0) })
            }
            if (event.key === 'Enter' && hits[active]) openHit(hits[active].variant_id)
          }}
        />
        <ul className="max-h-80 overflow-auto">
          {query.isPending ? <li className="px-3 py-3 text-sm text-mute">Opening the index…</li> : null}
          {query.isError ? <li className="px-3 py-3 text-sm">The index did not load.</li> : null}
          {hits.map((hit, index) => (
            <li key={hit.variant_id}>
              <button
                type="button"
                className={`flex min-h-11 w-full items-center justify-between gap-3 px-3 text-left ${index === active ? 'bg-panel' : ''}`}
                onMouseEnter={() => choose(index)}
                onClick={() => openHit(hit.variant_id)}
              >
                <span>
                  <span className="font-mono text-sm">{caseLabel(hit.variant_id)}</span>
                  <span className="mt-0.5 block text-xs text-mute">
                    {hit.gene_symbol} {hit.protein_change ?? ''} · residue {hit.position}
                  </span>
                </span>
                <ClassMark label={hit.label} />
              </button>
            </li>
          ))}
          {!query.isPending && hits.length === 0 ? (
            <li className="px-3 py-3 text-sm text-mute">No case matches that text.</li>
          ) : null}
        </ul>
      </div>
    </div>
  )
}
