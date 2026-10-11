import { useEffect, useRef, useState } from 'react'
import { useNavigate } from 'react-router-dom'
import { useQuery } from '@tanstack/react-query'
import { motion, useReducedMotion } from 'framer-motion'
import { Moon, Sun } from 'lucide-react'
import { api } from '../../api/client.ts'
import { Swatch } from '../../components/ui.tsx'
import { useTheme } from '../../lib/theme.ts'

const ENTRY_KEY = 'specimen-unsealed'

const EMPTY_LANES = Array.from({ length: 8 }, (_, index) => ({
  id: index + 1,
  name: `Lane ${String(index + 1).padStart(2, '0')}`,
  status: 'not_started',
}))

function rememberEntry() {
  try {
    sessionStorage.setItem(ENTRY_KEY, '1')
  } catch {
    /* private mode */
  }
}

export function LandingPage() {
  const navigate = useNavigate()
  const reduce = useReducedMotion()
  const theme = useTheme()
  const overview = useQuery({ queryKey: ['overview'], queryFn: api.overview })
  const [opening, setOpening] = useState(false)
  const [lit, setLit] = useState(-1)
  const timer = useRef<number | null>(null)
  const tick = useRef<number | null>(null)
  const steps = overview.data?.pipeline?.steps?.slice(0, 8) ?? EMPTY_LANES

  useEffect(() => {
    return () => {
      if (timer.current != null) window.clearTimeout(timer.current)
      if (tick.current != null) window.clearInterval(tick.current)
    }
  }, [])

  function enter() {
    if (tick.current != null) window.clearInterval(tick.current)
    if (timer.current != null) window.clearTimeout(timer.current)
    rememberEntry()
    void navigate('/bench')
  }

  function openDrawer() {
    if (opening) return
    if (reduce) {
      enter()
      return
    }
    setOpening(true)
    setLit(0)
    let index = 0
    tick.current = window.setInterval(() => {
      index += 1
      if (index >= steps.length) {
        if (tick.current != null) window.clearInterval(tick.current)
        return
      }
      setLit(index)
    }, 140)
    timer.current = window.setTimeout(() => {
      if (tick.current != null) window.clearInterval(tick.current)
      enter()
    }, 1500)
  }

  const caption = opening
    ? lit >= 0 && steps[lit]
      ? steps[lit].name
      : 'The drawer is opening.'
    : 'Eight lanes. They stay dark until the drawer opens.'

  return (
    <div className="grain flex min-h-dvh flex-col bg-paper text-ink">
      <a className="skip" href="#open-drawer">
        Skip to open
      </a>
      {api.mocks ? (
        <div className="border-b border-ink bg-pathogenic px-4 py-2 text-sm text-[#14161a]">
          <strong className="font-mono tracking-[0.16em]">MOCK DATA</strong>
          <span className="ml-2">
            Probabilities, SHAP, curves, conservation coverage, and sequence flanks are invented. Clinical labels and
            physchem deltas are measured.
          </span>
        </div>
      ) : null}
      <header className="flex items-center justify-between px-4 py-3">
        <p className="font-mono text-[10px] tracking-[0.22em] uppercase">Cardiovascular variants</p>
        <button
          type="button"
          className="control inline-flex items-center gap-2"
          onClick={theme.toggle}
          aria-pressed={theme.theme === 'dark'}
        >
          {theme.theme === 'dark' ? <Sun size={16} aria-hidden /> : <Moon size={16} aria-hidden />}
          <span>{theme.theme === 'dark' ? 'Paper' : 'UV'}</span>
        </button>
      </header>
      <main className="flex flex-1 items-center justify-center px-4 py-8">
        <motion.section
          className="w-full max-w-xl border-2 border-ink bg-paper p-6 md:p-8"
          style={{ boxShadow: '8px 8px 0 var(--fg)' }}
          animate={opening ? { y: -8 } : { y: 0 }}
          transition={{ duration: 0.6, ease: 'easeOut' }}
        >
          <motion.p
            className="inline-block border-2 border-ink px-2 py-1 font-mono text-[11px] tracking-[0.22em] text-ink uppercase"
            animate={opening ? { rotate: -16, opacity: 0.35 } : { rotate: -8, opacity: 1 }}
          >
            Sealed
          </motion.p>
          <h1 className="mt-5 font-display text-5xl leading-[0.95] md:text-7xl">Specimen Archive</h1>
          <p className="mt-4 max-w-md text-mute">
            A research bench for labelled heart-gene variants, the VUS held beside them, and two classifiers. Not a
            clinical report.
          </p>
          <div className="mt-4 flex flex-wrap gap-4">
            <Swatch tone="pathogenic">Pathogenic</Swatch>
            <Swatch tone="benign">Benign</Swatch>
            <Swatch tone="vus">VUS</Swatch>
          </div>
          <div className="mt-8 border border-ink bg-panel p-3" aria-hidden={false}>
            <div className="flex h-28 items-end gap-1.5">
              {steps.map((step, index) => {
                const dashed = step.status === 'done_poc' || step.status === 'partial' || step.status === 'stub'
                const active = opening && index <= lit
                return (
                  <div key={step.id} className="flex h-full min-w-0 flex-1 flex-col justify-end">
                    <div className="relative h-full overflow-hidden border rule bg-paper">
                      <motion.div
                        className={`absolute inset-x-0 bottom-0 ${dashed ? 'border border-dashed border-ink bg-graphite/50' : 'bg-ink'}`}
                        initial={false}
                        animate={{ height: active ? '100%' : '0%' }}
                        transition={{ duration: reduce ? 0 : 0.35, ease: 'easeOut' }}
                      />
                    </div>
                  </div>
                )
              })}
            </div>
            <div className="mt-2 flex justify-between font-mono text-[10px] tracking-wide text-mute">
              {steps.map((step) => (
                <span key={step.id}>{String(step.id).padStart(2, '0')}</span>
              ))}
            </div>
          </div>
          <p className="mt-3 min-h-6 font-mono text-xs" aria-live="polite">
            {caption}
          </p>
          <div className="mt-4 flex flex-wrap items-center gap-3">
            <button
              id="open-drawer"
              type="button"
              className="control bg-ink px-5 text-paper disabled:opacity-60"
              onClick={openDrawer}
              disabled={opening}
              aria-busy={opening}
            >
              {opening ? 'Opening the drawer' : 'Open the drawer'}
            </button>
            {opening ? (
              <button type="button" className="min-h-11 px-2 text-sm underline" onClick={enter}>
                Skip ahead
              </button>
            ) : null}
          </div>
        </motion.section>
      </main>
      <footer className="px-4 py-4 text-xs text-mute">Research prototype — not for clinical use</footer>
    </div>
  )
}
