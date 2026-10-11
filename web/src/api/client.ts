import type {
  Ablation,
  GlobalShap,
  Health,
  ListQuery,
  Metrics,
  Overview,
  Page,
  Report,
  StructurePayload,
  SuggestHit,
  VariantDetail,
  VariantSummary,
  VusSummary,
} from './types.ts'
import * as mocks from '../mocks/assemble.ts'

const mockFlag = import.meta.env.VITE_USE_MOCKS
export const useMocks = mockFlag === 'true' ? true : mockFlag === 'false' ? false : import.meta.env.DEV

const configuredApi = import.meta.env.VITE_API_URL?.replace(/\/$/, '') ?? ''
const apiBase = configuredApi || (import.meta.env.DEV ? 'http://127.0.0.1:8000' : '')

export class ApiError extends Error {
  status: number
  reason: string

  constructor(status: number, reason: string) {
    super(reason)
    this.name = 'ApiError'
    this.status = status
    this.reason = reason
  }
}

function throwStatus(error: unknown): never {
  if (error instanceof ApiError) throw error
  if (error instanceof Error && 'status' in error && typeof error.status === 'number') {
    throw new ApiError(error.status, error.message)
  }
  throw error
}

async function request(path: string) {
  let response: Response
  try {
    response = await fetch(`${apiBase}${path}`)
  } catch (error) {
    const message = error instanceof Error ? error.message : 'Network request failed'
    throw new ApiError(0, message)
  }
  if (!response.ok) {
    let reason = response.statusText || 'Request failed'
    try {
      const body: unknown = await response.json()
      if (body && typeof body === 'object' && 'detail' in body) {
        const detail = body.detail
        if (detail && typeof detail === 'object' && 'reason' in detail && typeof detail.reason === 'string') {
          reason = detail.reason
        }
      }
    } catch {
      /* keep the status text */
    }
    throw new ApiError(response.status, reason)
  }
  return response
}

function queryString(query: ListQuery) {
  const params = new URLSearchParams()
  if (query.q) params.set('q', query.q)
  if (query.gene) params.set('gene', query.gene)
  if (query.label) params.set('label', query.label)
  if (query.split) params.set('split', query.split)
  params.set('page', String(query.page ?? 1))
  params.set('page_size', String(query.page_size ?? 25))
  const text = params.toString()
  return text ? `?${text}` : ''
}

async function getJson<T>(path: string): Promise<T> {
  const response = await request(path)
  return (await response.json()) as T
}

export function fileUrl(url: string | null) {
  if (!url) return null
  if (url.startsWith('http://') || url.startsWith('https://')) return url
  return `${apiBase}${url.startsWith('/') ? url : `/${url}`}`
}

export const api = {
  mocks: useMocks,
  overview: () => (useMocks ? mocks.overview() : getJson<Overview>('/overview')),
  health: () => (useMocks ? mocks.health() : getJson<Health>('/health')),
  metrics: () => (useMocks ? mocks.metrics() : getJson<Metrics>('/metrics')),
  globalShap: () => (useMocks ? mocks.globalShap() : getJson<GlobalShap>('/shap/global')),
  ablation: () => (useMocks ? mocks.ablation() : getJson<Ablation>('/ablation')),
  variants: (query: ListQuery) =>
    useMocks
      ? mocks.variants(query).catch(throwStatus)
      : getJson<Page<VariantSummary>>(`/variants${queryString(query)}`),
  vus: (query: ListQuery) =>
    useMocks ? mocks.vus(query).catch(throwStatus) : getJson<Page<VusSummary>>(`/vus${queryString(query)}`),
  variant: (id: string) =>
    useMocks ? mocks.variant(id).catch(throwStatus) : getJson<VariantDetail>(`/variants/${encodeURIComponent(id)}`),
  structure: (id: string) =>
    useMocks
      ? mocks.structure(id).catch(throwStatus)
      : getJson<StructurePayload>(`/variants/${encodeURIComponent(id)}/structure`),
  report: async (id: string, format: 'json' | 'markdown' = 'json') => {
    if (useMocks) return mocks.report(id, format).catch(throwStatus)
    const response = await request(`/variants/${encodeURIComponent(id)}/report?format=${format}`)
    if (format === 'markdown') return response.text()
    return (await response.json()) as Report
  },
  suggest: async (q: string): Promise<SuggestHit[]> => {
    if (useMocks) return mocks.suggest(q)
    const [labelled, held] = await Promise.all([
      getJson<Page<VariantSummary>>(`/variants${queryString({ q, page: 1, page_size: 8 })}`),
      getJson<Page<VusSummary>>(`/vus${queryString({ q, page: 1, page_size: 8 })}`),
    ])
    return [
      ...labelled.items.map((item) => ({
        variant_id: item.variant_id,
        gene_symbol: item.gene_symbol,
        label: item.label,
        protein_change: item.protein_change,
        position: item.position,
      })),
      ...held.items.map((item) => ({
        variant_id: item.variant_id,
        gene_symbol: item.gene_symbol,
        label: item.label,
        protein_change: item.protein_change,
        position: item.position,
      })),
    ].slice(0, 12)
  },
}
