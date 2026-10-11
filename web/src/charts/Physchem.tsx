import type { Physchem } from '../api/types.ts'
import { formatNumber } from '../lib/format.ts'
import { Check, Minus } from 'lucide-react'

const AXES = [
  { key: 'delta_hydrophobicity', label: 'Δ hydrophobicity', limit: 5 },
  { key: 'delta_molecular_weight', label: 'Δ mass', limit: 150 },
  { key: 'delta_volume', label: 'Δ volume', limit: 150 },
  { key: 'delta_charge', label: 'Δ charge', limit: 2 },
] as const

export function PhyschemBars({ physchem }: { physchem: Physchem }) {
  return (
    <section>
      <h3 className="font-display text-xl">Physicochemical deltas</h3>
      <p className="text-xs text-mute">
        Source: {physchem.source}. Axis limits are display scales, not clinical cutoffs. Bars grow left for a negative
        delta and right for a positive delta.
      </p>
      <ul className="mt-3 space-y-2">
        {AXES.map((axis) => {
          const value = physchem[axis.key]
          const fraction = Math.max(-1, Math.min(1, value / axis.limit))
          return (
            <li key={axis.key} className="grid grid-cols-[8.5rem_1fr_4.5rem] items-center gap-2">
              <span className="text-sm">{axis.label}</span>
              <span className="relative h-6 bg-panel">
                <span className="absolute top-0 bottom-0 left-1/2 w-px bg-ink" aria-hidden />
                <span
                  className="absolute top-1 h-4 bg-ink"
                  style={{
                    width: `${Math.abs(fraction) * 50}%`,
                    left: fraction >= 0 ? '50%' : `${50 - Math.abs(fraction) * 50}%`,
                  }}
                />
              </span>
              <span className="text-right font-mono text-sm">{formatNumber(value, 3)}</span>
            </li>
          )
        })}
      </ul>
      <ul className="mt-3 flex flex-wrap gap-2 text-sm">
        <Flag on={physchem.polarity_change} label="Polarity change" />
        <Flag on={physchem.wt_is_proline} label="Wild type proline" />
        <Flag on={physchem.mut_is_proline} label="Mutant proline" />
        <Flag on={physchem.wt_is_glycine} label="Wild type glycine" />
        <Flag on={physchem.mut_is_glycine} label="Mutant glycine" />
      </ul>
    </section>
  )
}

function Flag({ on, label }: { on: boolean; label: string }) {
  return (
    <li className="inline-flex min-h-11 items-center gap-1 border rule px-2">
      {on ? <Check size={14} aria-hidden /> : <Minus size={14} aria-hidden />}
      <span>
        {label}: {on ? 'yes' : 'no'}
      </span>
    </li>
  )
}
