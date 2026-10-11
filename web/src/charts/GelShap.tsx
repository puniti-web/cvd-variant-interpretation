import { useState } from 'react'
import type { LocalShap, ShapFeature } from '../api/types.ts'
import { featureLabel, formatProb } from '../lib/format.ts'
import { Missing, MockTag } from '../components/ui.tsx'

export function GelShap({ shap, title }: { shap: LocalShap; title: string }) {
  const [mode, setMode] = useState<'gel' | 'waterfall'>('gel')
  const empty = shap.features.length === 0 || shap.base_value == null
  return (
    <section className="border rule p-3">
      <div className="mb-3 flex flex-wrap items-center justify-between gap-2">
        <h3 className="font-display text-xl">
          {title}
          {shap.mock ? <MockTag /> : null}
        </h3>
        <div className="flex" role="group" aria-label={`${title} view`}>
          <button
            type="button"
            className={`control ${mode === 'gel' ? 'bg-ink text-paper' : ''}`}
            aria-pressed={mode === 'gel'}
            onClick={() => setMode('gel')}
          >
            Gel lanes
          </button>
          <button
            type="button"
            className={`control border-l-0 ${mode === 'waterfall' ? 'bg-ink text-paper' : ''}`}
            aria-pressed={mode === 'waterfall'}
            onClick={() => setMode('waterfall')}
          >
            Waterfall
          </button>
        </div>
      </div>
      {empty ? (
        <Missing reason={shap.reason} />
      ) : mode === 'gel' ? (
        <GelLanes shap={shap} />
      ) : (
        <Waterfall shap={shap} />
      )}
      <p className="mt-3 text-xs text-mute">
        Left of the hairline lowers pathogenic probability. Right of it raises pathogenic probability. Bar length and
        brightness follow the absolute SHAP value. The number is the probability-space contribution.
      </p>
    </section>
  )
}

function GelLanes({ shap }: { shap: LocalShap }) {
  const maxAbs = Math.max(...shap.features.map((feature) => Math.abs(feature.shap_value)), 1e-6)
  return (
    <div className="space-y-1">
      <div className="grid grid-cols-[7rem_1fr_4.5rem] gap-2 font-mono text-[10px] tracking-wide text-mute uppercase">
        <span>Benign</span>
        <span className="text-center">Lane</span>
        <span className="text-right">Pathogenic</span>
      </div>
      {shap.features.map((feature) => {
        const magnitude = Math.abs(feature.shap_value) / maxAbs
        const positive = feature.shap_value >= 0
        return (
          <div key={feature.feature} className="grid grid-cols-[minmax(0,9rem)_1fr_4.5rem] items-center gap-2">
            <span className="truncate font-mono text-xs" title={feature.feature}>
              {featureLabel(feature.feature)}
            </span>
            <div className="relative h-7">
              <div className="absolute top-0 bottom-0 left-1/2 w-px bg-ink" aria-hidden />
              <div
                className={positive ? 'bg-pathogenic' : 'bg-benign'}
                style={{
                  position: 'absolute',
                  top: 4,
                  height: 20,
                  width: `${magnitude * 48}%`,
                  left: positive ? '50%' : undefined,
                  right: positive ? undefined : '50%',
                  opacity: 0.4 + magnitude * 0.6,
                }}
                title={`${feature.feature} ${feature.shap_value}`}
              />
            </div>
            <span className="text-right font-mono text-xs">
              {positive ? '+' : '−'}
              {formatProb(Math.abs(feature.shap_value))}
            </span>
          </div>
        )
      })}
    </div>
  )
}

function Waterfall({ shap }: { shap: LocalShap }) {
  const base = shap.base_value ?? 0
  const { steps, cursor } = walkShap(shap.features, base)
  const output = shap.output_value ?? cursor
  const samples = [base, output, ...steps.flatMap((step) => [step.start, step.end])]
  const min = Math.min(...samples)
  const max = Math.max(...samples)
  const span = max - min || 1
  const xOf = (value: number) => ((value - min) / span) * 100
  return (
    <div className="space-y-1">
      <p className="font-mono text-xs">
        Base {formatProb(base)} → output {formatProb(output)}
      </p>
      <Bar label="Base" start={xOf(Math.min(0, base))} end={xOf(Math.max(0, base))} tone="graphite" value={base} />
      {steps.map((step) => (
        <Bar
          key={step.feature.feature}
          label={featureLabel(step.feature.feature)}
          start={xOf(Math.min(step.start, step.end))}
          end={xOf(Math.max(step.start, step.end))}
          tone={step.feature.shap_value >= 0 ? 'pathogenic' : 'benign'}
          value={step.feature.shap_value}
        />
      ))}
      <Bar label="Output" start={xOf(Math.min(0, output))} end={xOf(Math.max(0, output))} tone="graphite" value={output} />
    </div>
  )
}

function walkShap(features: ShapFeature[], base: number) {
  const steps: { feature: ShapFeature; start: number; end: number }[] = []
  let cursor = base
  for (const feature of features) {
    const start = cursor
    cursor += feature.shap_value
    steps.push({ feature, start, end: cursor })
  }
  return { steps, cursor }
}

function Bar({
  label,
  start,
  end,
  tone,
  value,
}: {
  label: string
  start: number
  end: number
  tone: 'pathogenic' | 'benign' | 'graphite'
  value: number
}) {
  const fill = tone === 'pathogenic' ? 'bg-pathogenic' : tone === 'benign' ? 'bg-benign' : 'bg-graphite'
  return (
    <div className="grid grid-cols-[minmax(0,9rem)_1fr_4.5rem] items-center gap-2">
      <span className="truncate font-mono text-xs" title={label}>
        {label}
      </span>
      <div className="relative h-6 bg-panel">
        <div
          className={fill}
          style={{ position: 'absolute', top: 4, height: 16, left: `${start}%`, width: `${Math.max(end - start, 0.8)}%` }}
        />
      </div>
      <span className="text-right font-mono text-xs">
        {value >= 0 ? '+' : '−'}
        {formatProb(Math.abs(value))}
      </span>
    </div>
  )
}
