import type { Maybe, SequenceWindow } from '../api/types.ts'
import { Missing, MockTag } from '../components/ui.tsx'

export function SequenceRibbon({
  sequence,
  position,
  mut,
}: {
  sequence: Maybe<SequenceWindow>
  position: number
  mut: string
}) {
  const windowed = sequence.value
  if (!windowed) return <Missing reason={sequence.reason} />
  const start = windowed.sequence_start ?? 1
  return (
    <section>
      <h3 className="font-display text-xl">
        Sequence window
        {sequence.mock ? <MockTag /> : null}
      </h3>
      <p className="text-xs text-mute">
        UniProt {windowed.uniprot_id}. Violet depth is the conservation value at that residue. The marked site is the
        recorded wild-type residue.
      </p>
      <div className="mt-3 flex flex-wrap gap-1" aria-label="Sequence window">
        {windowed.sequence.split('').map((letter, index) => {
          const residue = start + index
          const heat = windowed.conservation_by_position?.[String(residue)]
          const site = residue === position
          return (
            <span key={`${residue}-${letter}`} className="inline-flex w-7 flex-col items-center">
              <span className="h-4 font-mono text-[10px] text-ink">{site ? mut : ''}</span>
              <span
                className={`flex h-8 w-7 items-center justify-center font-mono text-sm text-ink ${site ? 'outline outline-2 outline-ink' : ''}`}
                style={{ background: `rgba(139, 124, 255, ${heat == null ? 0 : Math.min(1, heat)})` }}
                title={heat == null ? `Residue ${residue}` : `Residue ${residue}, conservation ${heat}`}
              >
                {letter}
              </span>
              <span className="font-mono text-[9px] text-mute">{site ? residue : ''}</span>
            </span>
          )
        })}
      </div>
      {windowed.note ? <p className="mt-2 text-xs text-mute">{windowed.note}</p> : null}
      {windowed.domains.length === 0 ? (
        <p className="mt-1 text-xs text-mute">No domain spans are attached to this window.</p>
      ) : (
        <ul className="mt-2 text-sm">
          {windowed.domains.map((domain) => (
            <li key={`${domain.type}-${domain.start}`}>
              {domain.type}
              {domain.description ? ` · ${domain.description}` : ''} ({domain.start ?? '—'}–{domain.end ?? '—'})
            </li>
          ))}
        </ul>
      )}
    </section>
  )
}
