export function caseLabel(variantId: string) {
  const breakAt = variantId.indexOf('_')
  const pretty = breakAt === -1 ? variantId : `${variantId.slice(0, breakAt)}-${variantId.slice(breakAt + 1)}`
  return `CASE № ${pretty}`
}

export function featureLabel(name: string) {
  return name.replace(/^(num|cat)__/, '').replaceAll('_', ' ')
}

export function formatNumber(value: number, digits = 3) {
  return value.toLocaleString('en-US', {
    minimumFractionDigits: 0,
    maximumFractionDigits: digits,
  })
}

export function formatProb(value: number) {
  return value.toLocaleString('en-US', {
    minimumFractionDigits: 2,
    maximumFractionDigits: 2,
  })
}

export function formatCount(value: number) {
  return value.toLocaleString('en-US')
}

export function classWord(label: string | null | undefined) {
  if (label === 'pathogenic') return 'Pathogenic'
  if (label === 'benign') return 'Benign'
  if (label === 'vus') return 'VUS'
  return 'Not computed'
}

export function splitWord(split: string | null | undefined) {
  if (split === 'train_gene') return 'Train gene'
  if (split === 'test_gene') return 'Test gene'
  return 'Not in split'
}
