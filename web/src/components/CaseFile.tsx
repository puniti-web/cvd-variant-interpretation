import { Link } from 'react-router-dom'
import type { VariantDetail } from '../api/types.ts'
import { ModelDials } from '../charts/Dials.tsx'
import { GelShap } from '../charts/GelShap.tsx'
import { PhyschemBars } from '../charts/Physchem.tsx'
import { SequenceRibbon } from '../charts/Ribbon.tsx'
import { Missing, MockTag, Stamp } from './ui.tsx'
import { formatProb, splitWord } from '../lib/format.ts'

export function CaseFile({ detail }: { detail: VariantDetail }) {
  const agreement =
    detail.predictions.rf.label && detail.predictions.xgb.label
      ? detail.predictions.rf.label === detail.predictions.xgb.label
        ? 'agree'
        : 'disagree'
      : null
  return (
    <article className="space-y-5 border rule p-4">
      <div className="flex flex-wrap items-start justify-between gap-4">
        <Stamp variantId={detail.variant_id} label={detail.label} />
        <div className="min-w-0 text-sm">
          <p className="font-display text-2xl">{detail.hgvs_protein ?? detail.variant_id}</p>
          <p>
            {detail.gene_symbol} · {detail.panel}
            {detail.condition ? ` · ${detail.condition}` : ''}
          </p>
          <p className="font-mono text-xs text-mute">UniProt {detail.uniprot_id ?? 'Not computed'}</p>
        </div>
      </div>
      {detail.split.value === 'train_gene' ? (
        <p className="border border-dashed border-graphite p-3 text-sm" role="note">
          This gene was in the training fold. A probability here is not an independent test of the model.
        </p>
      ) : null}
      <dl className="grid gap-3 text-sm sm:grid-cols-2">
        <Fact term="Clinical significance" value={detail.clinical.clinical_significance} />
        <Fact term="Review status" value={detail.clinical.review_status} />
        <Fact term="ClinVar title" value={detail.clinical.title} />
        <Fact term="cDNA" value={detail.clinical.cdna_change} />
        <Fact term="dbSNP" value={detail.clinical.dbsnp} />
        <Fact
          term="GRCh38"
          value={
            detail.clinical.location
              ? `${detail.clinical.location.chr ?? '—'}:${detail.clinical.location.start ?? '—'}`
              : null
          }
        />
        <div>
          <dt className="text-xs tracking-wide text-mute uppercase">Fold</dt>
          <dd>{detail.split.value ? splitWord(detail.split.value) : <Missing reason={detail.split.reason} />}</dd>
        </div>
        <div>
          <dt className="text-xs tracking-wide text-mute uppercase">Structure</dt>
          <dd>{detail.has_structure ? 'Proof-of-concept PDB pair' : 'No proof-of-concept structure'}</dd>
        </div>
      </dl>
      {detail.clinical.parse_reason ? <p className="text-xs text-mute">{detail.clinical.parse_reason}</p> : null}
      <div className="flex flex-wrap gap-2">
        <Link className="control inline-flex items-center" to={`/structure/${detail.variant_id}`}>
          Structure
        </Link>
        <Link className="control inline-flex items-center" to={`/report/${detail.variant_id}`}>
          Report
        </Link>
      </div>
      <PhyschemBars physchem={detail.physchem} />
      <section>
        <h3 className="font-display text-xl">
          Conservation
          {detail.conservation_score.mock ? <MockTag /> : null}
        </h3>
        {detail.conservation_score.value == null ? (
          <Missing reason={detail.conservation_score.reason} />
        ) : (
          <p className="font-mono text-2xl">{formatProb(detail.conservation_score.value)}</p>
        )}
        {detail.conservation_imputed ? (
          <p className="mt-1 text-sm">
            Imputed
            {detail.conservation_imputed_value.value == null ? (
              <Missing reason={detail.conservation_imputed_value.reason} />
            ) : (
              <span className="font-mono"> {formatProb(detail.conservation_imputed_value.value)}</span>
            )}
          </p>
        ) : (
          <p className="text-xs text-mute">Not imputed.</p>
        )}
        <p className="mt-2 text-sm">
          Annotated domain:{' '}
          {detail.in_annotated_domain.value == null ? (
            <Missing reason={detail.in_annotated_domain.reason} />
          ) : detail.in_annotated_domain.value ? (
            'In an annotated domain'
          ) : (
            'Checked: not in an annotated domain'
          )}
          {detail.domain_type.value ? ` (${detail.domain_type.value})` : ''}
        </p>
      </section>
      <SequenceRibbon sequence={detail.sequence} position={detail.position} mut={detail.mut.one} />
      <ModelDials rf={detail.predictions.rf} xgb={detail.predictions.xgb} agreement={agreement} />
      <GelShap shap={detail.shap.rf} title="Random Forest gel" />
      <GelShap shap={detail.shap.xgb} title="XGBoost gel" />
    </article>
  )
}

function Fact({ term, value }: { term: string; value: string | null | undefined }) {
  return (
    <div>
      <dt className="text-xs tracking-wide text-mute uppercase">{term}</dt>
      <dd>{value ? value : <Missing reason={null} />}</dd>
    </div>
  )
}
