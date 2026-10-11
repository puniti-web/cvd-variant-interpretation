import type { ReactNode } from 'react'
import type { UseQueryResult } from '@tanstack/react-query'
import { Check, CircleQuestionMark, Minus, TriangleAlert } from 'lucide-react'
import type { Label, Maybe } from '../api/types.ts'
import { ApiError } from '../api/client.ts'
import { caseLabel, classWord, formatCount, formatNumber } from '../lib/format.ts'

export function ClassMark({ label }: { label: string | null | undefined }) {
  const word = classWord(label)
  const icon =
    label === 'pathogenic' ? (
      <TriangleAlert size={16} aria-hidden />
    ) : label === 'benign' ? (
      <Check size={16} aria-hidden />
    ) : label === 'vus' ? (
      <CircleQuestionMark size={16} aria-hidden />
    ) : (
      <Minus size={16} aria-hidden />
    )
  return (
    <span className="inline-flex items-center gap-1 text-ink">
      {icon}
      <span>{word}</span>
    </span>
  )
}

export function Swatch({
  tone,
  children,
}: {
  tone: 'pathogenic' | 'benign' | 'vus' | 'graphite'
  children: ReactNode
}) {
  const bg =
    tone === 'pathogenic'
      ? 'bg-pathogenic'
      : tone === 'benign'
        ? 'bg-benign'
        : tone === 'vus'
          ? 'bg-vus'
          : 'bg-graphite'
  return (
    <span className="inline-flex items-center gap-2 text-sm">
      <span className={`inline-block h-3 w-3 ${bg}`} aria-hidden />
      <span>{children}</span>
    </span>
  )
}

export function MockTag() {
  return (
    <span className="ml-2 border border-dashed border-graphite px-1 font-mono text-[10px] tracking-widest text-ink uppercase">
      Mock
    </span>
  )
}

export function Missing({ reason }: { reason: string | null }) {
  return (
    <span className="inline-flex max-w-full flex-col gap-1">
      <span className="inline-flex items-center gap-1 border border-dashed border-graphite px-2 py-1 text-sm text-ink">
        <Minus size={14} aria-hidden />
        Not computed
      </span>
      {reason ? <span className="text-xs text-mute">{reason}</span> : null}
    </span>
  )
}

export function Readout({
  label,
  maybe,
  digits = 0,
}: {
  label: string
  maybe: Maybe<number>
  digits?: number
}) {
  return (
    <div className="min-w-0">
      <p className="text-xs tracking-wide text-mute uppercase">{label}</p>
      {maybe.value == null ? (
        <Missing reason={maybe.reason} />
      ) : (
        <p className="font-mono text-2xl text-ink">
          {digits > 0 ? formatNumber(maybe.value, digits) : formatCount(maybe.value)}
          {maybe.mock ? <MockTag /> : null}
        </p>
      )}
    </div>
  )
}

export function Stamp({ variantId, label }: { variantId: string; label: Label }) {
  const fill = label === 'pathogenic' ? 'bg-pathogenic' : label === 'benign' ? 'bg-benign' : 'bg-vus'
  return (
    <div className={`stamp ${fill}`}>
      <ClassMark label={label} />
      <span>
        <span className="block text-[10px] tracking-[0.18em] uppercase">Specimen</span>
        <span className="font-mono text-sm">{caseLabel(variantId)}</span>
      </span>
    </div>
  )
}

export function GelSkeleton({ rows = 6 }: { rows?: number }) {
  return (
    <div className="space-y-2" aria-hidden>
      {Array.from({ length: rows }, (_, index) => (
        <div key={index} className="gel-shimmer h-8" style={{ width: `${88 - (index % 4) * 12}%` }} />
      ))}
      <p className="sr-only">Loading</p>
    </div>
  )
}

export function ErrorState({ reason, onRetry }: { reason: string; onRetry?: () => void }) {
  return (
    <div role="alert" className="border border-dashed border-graphite p-4">
      <p className="font-display text-xl">This drawer did not open</p>
      <p className="mt-2 text-sm text-mute">{reason}</p>
      {onRetry ? (
        <button type="button" className="control mt-4" onClick={onRetry}>
          Try again
        </button>
      ) : null}
    </div>
  )
}

export function EmptyState({ title, detail }: { title: string; detail?: string }) {
  return (
    <div className="border border-dashed border-graphite p-4">
      <p className="font-display text-xl">{title}</p>
      {detail ? <p className="mt-2 text-sm text-mute">{detail}</p> : null}
    </div>
  )
}

export function QueryGate<T>({
  query,
  children,
}: {
  query: UseQueryResult<T>
  children: (data: T) => ReactNode
}) {
  if (query.isPending) return <GelSkeleton />
  if (query.isError) {
    const reason = query.error instanceof ApiError ? query.error.reason : query.error.message
    return <ErrorState reason={reason} onRetry={() => void query.refetch()} />
  }
  if (query.data === undefined) return <EmptyState title="Nothing returned" />
  return children(query.data)
}

export function PageIntro({ kicker, title, lede }: { kicker: string; title: string; lede?: string }) {
  return (
    <header className="mb-6 max-w-3xl">
      <p className="font-mono text-xs tracking-[0.18em] text-mute uppercase">{kicker}</p>
      <h1 className="font-display text-4xl leading-tight">{title}</h1>
      {lede ? <p className="mt-2 text-mute">{lede}</p> : null}
    </header>
  )
}
