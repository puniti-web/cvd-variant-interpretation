import { useRef } from 'react'
import { motion, useReducedMotion, useScroll, useTransform } from 'framer-motion'
import type { PipelineStatus } from '../api/types.ts'
import { Missing } from '../components/ui.tsx'

function lit(status: string) {
  return status === 'done' || status === 'done_poc' || status === 'partial'
}

export function PipelineHelix({ pipeline, reason }: { pipeline: PipelineStatus | null; reason: string | null }) {
  const ref = useRef<HTMLDivElement>(null)
  const reduce = useReducedMotion()
  const { scrollYProgress } = useScroll({ target: ref, offset: ['start 85%', 'end 35%'] })
  const travel = useTransform(scrollYProgress, [0, 1], ['0%', '92%'])
  if (!pipeline) return <Missing reason={reason} />
  const lastLit = Math.max(
    0,
    pipeline.steps.reduce((index, step, stepIndex) => (lit(step.status) ? stepIndex : index), 0),
  )
  const staticTop = pipeline.steps.length <= 1 ? '0%' : `${(lastLit / (pipeline.steps.length - 1)) * 92}%`
  return (
    <div ref={ref} className="relative">
      <ol className="relative space-y-3 border-l border-dashed border-graphite pl-6">
        <motion.span
          aria-hidden
          className="absolute -left-[5px] h-2.5 w-2.5 rounded-full bg-ink"
          style={{ top: reduce ? staticTop : travel }}
        />
        {pipeline.steps.map((step) => (
          <li key={step.id}>
            <div className="flex flex-wrap items-baseline gap-2">
              <span className="font-mono text-xs text-mute">{String(step.id).padStart(2, '0')}</span>
              <span className="font-display text-lg">{step.name}</span>
            </div>
            <p className="text-sm text-mute">{step.note}</p>
          </li>
        ))}
      </ol>
    </div>
  )
}
