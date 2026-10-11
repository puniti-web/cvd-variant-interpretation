import { useQuery } from '@tanstack/react-query'
import { api } from '../../api/client.ts'
import { PageIntro, QueryGate, Readout } from '../../components/ui.tsx'

export function AboutPage() {
  const overview = useQuery({ queryKey: ['overview'], queryFn: api.overview })
  return (
    <div className="max-w-3xl space-y-6">
      <PageIntro
        kicker="About"
        title="A research bench, not a clinical report."
        lede="The archive reads ClinVar labels, UniProt context, and one proof-of-concept PDB. It does not write a diagnosis."
      />
      <QueryGate query={overview}>
        {(data) => (
          <section className="grid gap-4 sm:grid-cols-3">
            <Readout label="Labelled" maybe={data.counts.labelled} />
            <Readout label="VUS held out" maybe={data.counts.vus} />
            <Readout label="Genes" maybe={{ value: data.genes.length, reason: null }} />
          </section>
        )}
      </QueryGate>
      <section className="space-y-3 text-sm">
        <h2 className="font-display text-2xl">Sources</h2>
        <p>Clinical labels and titles come from ClinVar. Gene identity, conditions, and UniProt accessions come from the project gene panel. Coordinates, when present, come from a local PDB pair rather than a live RCSB proxy.</p>
        <p>The labelled table is split by gene with GroupShuffleSplit, test size 0.25, random state 42. The same gene does not appear in both folds. A comment elsewhere in the training module still says GroupKFold; the function that runs is the shuffle split.</p>
        <p>Training uses class weights. SMOTE is not applied. Variants of uncertain significance are held out of the feature matrix used for fitting and are not assigned a fold.</p>
        <p>Conservation is a homolog-identity proxy computed by the project. It is not ConSurf, phyloP, or phastCons. Domain membership comes from UniProt feature spans. When that file is missing, an unknown domain is left null instead of being stored as false.</p>
      </section>
      <section className="space-y-3 text-sm">
        <h2 className="font-display text-2xl">Limits</h2>
        <ul className="list-disc space-y-2 pl-5">
          <li>Panel is a gene-level covariate inside the model, so it can stand in for the gene.</li>
          <li>Position is an absolute residue index and is not comparable across proteins.</li>
          <li>The structure proof of concept renames one residue and leaves coordinates unchanged. A Cα RMSD of 0.00 Å on that pair is a method artifact.</li>
          <li>FoldX ΔΔG is not implemented. The stability number stays empty.</li>
          <li>The public ESMFold endpoint is not called. The proof-of-concept files are not ESMFold models.</li>
          <li>Per-variant SHAP is computed when a saved model is loaded. A global plot file is not parsed into numbers.</li>
        </ul>
      </section>
      <p className="border border-ink p-3 font-display text-2xl">Research prototype — not for clinical use</p>
    </div>
  )
}
