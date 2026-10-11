import type { Prediction } from '../api/types.ts'
import { formatProb } from '../lib/format.ts'
import { ClassMark, Missing, MockTag } from '../components/ui.tsx'

export function ModelDials({
  rf,
  xgb,
  agreement,
}: {
  rf: Prediction
  xgb: Prediction
  agreement: 'agree' | 'disagree' | null
}) {
  return (
    <section className="border rule p-3">
      <div className="mb-2 flex flex-wrap items-center justify-between gap-2">
        <h3 className="font-display text-xl">Model dials</h3>
        <Agreement agreement={agreement} />
      </div>
      <div className="grid gap-4 sm:grid-cols-2">
        <Dial name="Random Forest" prediction={rf} />
        <Dial name="XGBoost" prediction={xgb} />
      </div>
      <p className="mt-2 text-xs text-mute">
        Each dial is P(pathogenic). The call is pathogenic at 0.50 and above. The word under the dial is the call; the
        arc is only a mark.
      </p>
    </section>
  )
}

function Agreement({ agreement }: { agreement: 'agree' | 'disagree' | null }) {
  if (!agreement) return <Missing reason="Both model labels are required before agreement can be called." />
  const disagree = agreement === 'disagree'
  return (
    <span className={`border px-2 py-1 text-sm ${disagree ? 'border-dashed border-graphite' : 'border-ink'}`}>
      Models {disagree ? 'disagree' : 'agree'}
    </span>
  )
}

function Dial({ name, prediction }: { name: string; prediction: Prediction }) {
  const prob = prediction.prob_pathogenic
  const label =
    prob == null ? `${name} pathogenic probability not computed` : `${name} pathogenic probability ${formatProb(prob)}`
  return (
    <figure className="m-0">
      <svg viewBox="0 0 120 78" role="img" aria-label={label} className="mx-auto w-full max-w-56">
        <path d="M12 64 A48 48 0 0 1 108 64" fill="none" stroke="var(--graphite)" strokeWidth="8" strokeDasharray="4 4" />
        {prob != null ? <ValueArc probability={prob} /> : null}
        {prob != null ? <Needle probability={prob} /> : null}
        <text x="60" y="58" textAnchor="middle" fill="currentColor" fontFamily="JetBrains Mono, monospace" fontSize="16">
          {prob == null ? '—' : formatProb(prob)}
        </text>
      </svg>
      <figcaption className="text-center">
        <p className="font-medium">{name}</p>
        {prediction.mock ? <MockTag /> : null}
        <div className="mt-1 flex justify-center">
          {prob == null ? <Missing reason={prediction.reason} /> : <ClassMark label={prediction.label} />}
        </div>
      </figcaption>
    </figure>
  )
}

function ValueArc({ probability }: { probability: number }) {
  const clamped = Math.min(1, Math.max(0, probability))
  const start = Math.PI
  const end = Math.PI * (1 - clamped)
  const large = start - end > Math.PI ? 1 : 0
  const x1 = 60 + 48 * Math.cos(start)
  const y1 = 64 - 48 * Math.sin(start)
  const x2 = 60 + 48 * Math.cos(end)
  const y2 = 64 - 48 * Math.sin(end)
  const tone = clamped >= 0.5 ? '#ff7a1a' : '#3ddc84'
  return (
    <path
      d={`M ${x1} ${y1} A 48 48 0 ${large} 1 ${x2} ${y2}`}
      fill="none"
      stroke={tone}
      strokeWidth="8"
      strokeLinecap="butt"
    />
  )
}

function Needle({ probability }: { probability: number }) {
  const angle = Math.PI * (1 - Math.min(1, Math.max(0, probability)))
  const x = 60 + 36 * Math.cos(angle)
  const y = 64 - 36 * Math.sin(angle)
  return <line x1="60" y1="64" x2={x} y2={y} stroke="currentColor" strokeWidth="2" />
}
