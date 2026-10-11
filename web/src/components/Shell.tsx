import { useEffect, useState } from 'react'
import { NavLink, Outlet } from 'react-router-dom'
import { motion, useReducedMotion } from 'framer-motion'
import { Moon, Search, Sun } from 'lucide-react'
import { api } from '../api/client.ts'
import { useTheme } from '../lib/theme.ts'
import { CommandPalette } from './CommandPalette.tsx'

const ENTRY_KEY = 'specimen-unsealed'

const NAV = [
  { to: '/bench', label: 'Bench', end: true },
  { to: '/explore', label: 'Explore', end: false },
  { to: '/models', label: 'Models', end: false },
  { to: '/about', label: 'About', end: false },
]

function NavItems({ compact }: { compact?: boolean }) {
  return (
    <ul className={compact ? 'grid grid-cols-4' : 'space-y-1'}>
      {NAV.map((item) => (
        <li key={item.to}>
          <NavLink
            to={item.to}
            end={item.end}
            className={({ isActive }) =>
              `flex min-h-11 items-center px-3 ${compact ? 'justify-center text-xs' : ''} ${isActive ? 'bg-panel font-medium' : ''}`
            }
          >
            {item.label}
          </NavLink>
        </li>
      ))}
    </ul>
  )
}

export function Shell() {
  const theme = useTheme()
  const reduce = useReducedMotion()
  const [palette, setPalette] = useState(false)
  const [arrive] = useState(() => {
    try {
      return sessionStorage.getItem(ENTRY_KEY) === '1'
    } catch {
      return false
    }
  })

  useEffect(() => {
    try {
      sessionStorage.removeItem(ENTRY_KEY)
    } catch {
      /* private mode */
    }
  }, [])

  useEffect(() => {
    const onKey = (event: KeyboardEvent) => {
      if ((event.ctrlKey || event.metaKey) && event.key.toLowerCase() === 'k') {
        event.preventDefault()
        setPalette(true)
      }
    }
    window.addEventListener('keydown', onKey)
    return () => window.removeEventListener('keydown', onKey)
  }, [])

  return (
    <motion.div
      className="grain min-h-screen bg-paper text-ink"
      initial={arrive && !reduce ? { opacity: 0, y: 36 } : false}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.55, ease: 'easeOut' }}
    >
      <a className="skip" href="#content">
        Skip to content
      </a>
      {api.mocks ? (
        <div className="no-print sticky top-0 z-40 border-b border-ink bg-pathogenic px-4 py-2 text-sm text-[#14161a]">
          <strong className="font-mono tracking-[0.16em]">MOCK DATA</strong>
          <span className="ml-2">
            Probabilities, SHAP, curves, conservation coverage, and sequence flanks are invented. Clinical labels and
            physchem deltas are measured.
          </span>
        </div>
      ) : null}
      <div className="md:grid md:grid-cols-[14.5rem_minmax(0,1fr)]">
        <aside className="no-print sticky top-0 z-20 hidden h-screen flex-col border-r rule bg-paper md:flex">
          <div className="border-b rule px-4 py-5">
            <p className="font-mono text-[10px] tracking-[0.2em] uppercase">CVD variants</p>
            <p className="font-display text-2xl leading-none">Specimen Archive</p>
          </div>
          <nav aria-label="Primary" className="flex-1 py-3">
            <NavItems />
          </nav>
          <p className="border-t rule px-4 py-3 text-xs text-mute">Research prototype — not for clinical use</p>
        </aside>
        <div className="min-w-0">
          <header className="no-print flex items-center justify-between gap-3 border-b rule px-4 py-3">
            <p className="font-display text-xl md:hidden">Specimen Archive</p>
            <p className="hidden text-sm text-mute md:block">Explainable cardiovascular variant bench</p>
            <div className="flex gap-2">
              <button type="button" className="control inline-flex items-center gap-2" onClick={() => setPalette(true)}>
                <Search size={16} aria-hidden />
                <span>Find case</span>
                <span className="hidden font-mono text-xs text-mute sm:inline">Ctrl K</span>
              </button>
              <button
                type="button"
                className="control inline-flex items-center gap-2"
                onClick={theme.toggle}
                aria-pressed={theme.theme === 'dark'}
              >
                {theme.theme === 'dark' ? <Sun size={16} aria-hidden /> : <Moon size={16} aria-hidden />}
                <span>{theme.theme === 'dark' ? 'Paper' : 'UV'}</span>
              </button>
            </div>
          </header>
          <main id="content" className="px-4 py-6 pb-24 md:pb-12">
            <Outlet />
          </main>
          <p className="no-print px-4 pb-24 text-xs text-mute md:pb-8">Research prototype — not for clinical use</p>
        </div>
      </div>
      <nav aria-label="Primary" className="no-print fixed inset-x-0 bottom-0 z-40 border-t rule bg-paper md:hidden">
        <NavItems compact />
      </nav>
      {palette ? <CommandPalette onClose={() => setPalette(false)} /> : null}
    </motion.div>
  )
}
