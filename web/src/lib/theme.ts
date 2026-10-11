import { create } from 'zustand'

export type ThemeName = 'light' | 'dark'

const STORAGE_KEY = 'specimen-theme'

function preferredTheme(): ThemeName {
  try {
    const stored = localStorage.getItem(STORAGE_KEY)
    if (stored === 'light' || stored === 'dark') return stored
  } catch {
    /* private mode */
  }
  return window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light'
}

function applyTheme(theme: ThemeName) {
  document.documentElement.classList.toggle('dark', theme === 'dark')
}

interface ThemeState {
  theme: ThemeName
  toggle: () => void
}

export const useTheme = create<ThemeState>((set, get) => ({
  theme: preferredTheme(),
  toggle: () => {
    const next: ThemeName = get().theme === 'light' ? 'dark' : 'light'
    try {
      localStorage.setItem(STORAGE_KEY, next)
    } catch {
      /* ignore */
    }
    applyTheme(next)
    set({ theme: next })
  },
}))

export function bootTheme() {
  applyTheme(preferredTheme())
}
