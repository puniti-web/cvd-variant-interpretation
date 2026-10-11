import type { ReactNode } from 'react'
import { AxisBottom, AxisLeft } from '@visx/axis'
import { Group } from '@visx/group'
import { ParentSize } from '@visx/responsive'
import { scaleLinear } from '@visx/scale'
import { Line, LinePath } from '@visx/shape'
import type { CurvePoint, GlobalFeature } from '../api/types.ts'
import { featureLabel, formatNumber } from '../lib/format.ts'
import { Missing, MockTag } from '../components/ui.tsx'

export function RocPr({
  roc,
  pr,
  mock,
}: {
  roc: CurvePoint[] | undefined
  pr: CurvePoint[] | undefined
  mock?: boolean
}) {
  return (
    <div className="grid gap-4 lg:grid-cols-2">
      <ChartFrame title="ROC" mock={mock}>
        {roc && roc.length > 1 ? (
          <ParentSize>
            {({ width }) => (
              <Curve
                width={width}
                points={roc.map((point) => ({ x: point.fpr ?? 0, y: point.tpr ?? 0 }))}
                xLabel="False positive rate"
                yLabel="True positive rate"
                diagonal
              />
            )}
          </ParentSize>
        ) : (
          <Missing reason="ROC points were not returned." />
        )}
      </ChartFrame>
      <ChartFrame title="Precision–recall" mock={mock}>
        {pr && pr.length > 1 ? (
          <ParentSize>
            {({ width }) => (
              <Curve
                width={width}
                points={pr.map((point) => ({ x: point.recall ?? 0, y: point.precision ?? 0 }))}
                xLabel="Recall"
                yLabel="Precision"
              />
            )}
          </ParentSize>
        ) : (
          <Missing reason="Precision–recall points were not returned." />
        )}
      </ChartFrame>
    </div>
  )
}

function ChartFrame({ title, mock, children }: { title: string; mock?: boolean; children: ReactNode }) {
  return (
    <section className="border rule p-3">
      <h3 className="font-display text-xl">
        {title}
        {mock ? <MockTag /> : null}
      </h3>
      <div className="mt-2 h-56">{children}</div>
    </section>
  )
}

function Curve({
  width,
  points,
  xLabel,
  yLabel,
  diagonal,
}: {
  width: number
  points: { x: number; y: number }[]
  xLabel: string
  yLabel: string
  diagonal?: boolean
}) {
  if (width < 40) return null
  const height = 210
  const margin = { top: 8, right: 12, bottom: 36, left: 40 }
  const innerWidth = width - margin.left - margin.right
  const innerHeight = height - margin.top - margin.bottom
  const x = scaleLinear({ domain: [0, 1], range: [0, innerWidth] })
  const y = scaleLinear({ domain: [0, 1], range: [innerHeight, 0] })
  return (
    <svg width={width} height={height} role="img" aria-label={`${yLabel} against ${xLabel}`}>
      <Group left={margin.left} top={margin.top}>
        {diagonal ? (
          <Line from={{ x: x(0), y: y(0) }} to={{ x: x(1), y: y(1) }} stroke="var(--graphite)" strokeDasharray="4 4" />
        ) : null}
        <LinePath data={points} x={(point) => x(point.x)} y={(point) => y(point.y)} stroke="currentColor" strokeWidth={2} />
        <AxisBottom
          scale={x}
          top={innerHeight}
          stroke="currentColor"
          tickStroke="currentColor"
          tickLabelProps={{ fill: 'currentColor', fontSize: 10, fontFamily: 'JetBrains Mono, monospace' }}
          label={xLabel}
          labelProps={{ fill: 'currentColor', fontSize: 11 }}
        />
        <AxisLeft
          scale={y}
          stroke="currentColor"
          tickStroke="currentColor"
          tickLabelProps={{ fill: 'currentColor', fontSize: 10, fontFamily: 'JetBrains Mono, monospace' }}
          label={yLabel}
          labelProps={{ fill: 'currentColor', fontSize: 11 }}
        />
      </Group>
    </svg>
  )
}

export function MeanShapBars({
  features,
  mock,
  label,
}: {
  features: GlobalFeature[] | null
  mock?: boolean
  label: string
}) {
  if (!features || features.length === 0) return <Missing reason={`${label} global SHAP was not returned.`} />
  const max = Math.max(...features.map((feature) => feature.mean_abs_shap), 1e-6)
  return (
    <section className="border rule p-3">
      <h3 className="font-display text-xl">
        {label} mean |SHAP|
        {mock ? <MockTag /> : null}
      </h3>
      <ul className="mt-3 space-y-1">
        {features.map((feature) => (
          <li key={feature.feature} className="grid grid-cols-[minmax(0,10rem)_1fr_3.5rem] items-center gap-2">
            <span className="truncate font-mono text-xs" title={feature.feature}>
              {featureLabel(feature.feature)}
            </span>
            <span className="h-3 bg-panel">
              <span className="block h-3 bg-graphite" style={{ width: `${(feature.mean_abs_shap / max) * 100}%` }} />
            </span>
            <span className="text-right font-mono text-xs">{formatNumber(feature.mean_abs_shap, 3)}</span>
          </li>
        ))}
      </ul>
    </section>
  )
}

export function ConfusionTable({
  matrix,
  mock,
  label,
}: {
  matrix: [[number, number], [number, number]]
  mock?: boolean
  label: string
}) {
  const [[tn, fp], [fn, tp]] = matrix
  const cells = [
    ['True negative', tn, 'Benign called benign'],
    ['False positive', fp, 'Benign called pathogenic'],
    ['False negative', fn, 'Pathogenic called benign'],
    ['True positive', tp, 'Pathogenic called pathogenic'],
  ] as const
  return (
    <section>
      <h3 className="font-display text-xl">
        {label} confusion
        {mock ? <MockTag /> : null}
      </h3>
      <div className="mt-2 grid grid-cols-2 gap-2">
        {cells.map(([name, count, detail]) => (
          <div key={name} className="border rule p-3">
            <p className="text-xs tracking-wide text-mute uppercase">{name}</p>
            <p className="font-mono text-2xl">{count}</p>
            <p className="text-xs text-mute">{detail}</p>
          </div>
        ))}
      </div>
    </section>
  )
}
