"""
Step 6: Integrated variant interpretation report.

For one variant, combines everything earlier stages produced:
  - ClinVar record and disease context (config/genes.yaml)
  - model prediction: for labelled variants the out-of-fold XGBoost
    score from leave-one-gene-out CV (the gene was held out, so the model
    never saw the variant); for VUS the calibrated full XGBoost model.
    Tiers use a balanced-prior probability (see balanced_probability), so
    they are not inflated by the 85% pathogenic training set. The class-
    weighted XGBoost score is already on that scale.
  - top SHAP contributions
  - FoldX stability change (src/structure/foldx_ddg.py)
  - structural context (AlphaFold pLDDT, solvent accessibility, secondary
    structure, UniProt domain)
  - conservation and published predictors (AlphaMissense, REVEL) for
    reference

The conclusion is built from fixed rules, not free text, so every sentence
can be traced to a threshold below. Each evidence line is scored
independently; the report says how many agree with the model and flags
conflicts and known weak spots (genes with few benign training variants,
variants without a structural model).

This is a research tool. The tiers are model-support levels, not ACMG
classifications, and must not be used for clinical decisions.

Outputs (outputs/reports/):
    <GENE>_<WT><POS><MUT>.md / .json   single-variant reports
    vus_report.csv                     every VUS with tier and evidence summary
    vus_top<N>.md                      ranked report for the top N VUS

Usage:
    python -m src.report.generate_report --variant MYH7:R403Q
    python -m src.report.generate_report --top-vus 50
    python -m src.report.generate_report          # examples + VUS report
"""
import argparse
import json
import math
from dataclasses import asdict, dataclass, field
from functools import lru_cache

import pandas as pd

from src.config import FEATURES_DIR, OUTPUTS_DIR, PROCESSED_DIR, REPORTS_DIR, ensure_dirs, load_gene_panel

# ---------------------------------------------------------------- rules

# Model support tiers on P(pathogenic). Not ACMG classes.
TIERS = [
    (0.90, "strong model support for pathogenic"),
    (0.70, "moderate model support for pathogenic"),
    (0.30, "inconclusive"),
    (0.10, "moderate model support for benign"),
    (0.00, "strong model support for benign"),
]

# FoldX ddG (kcal/mol). FoldX's own error is about 1 kcal/mol, so smaller
# changes are treated as neutral.
DDG_STRONG = 2.0
DDG_DESTABILISING = 1.0
DDG_STABILISING = -1.0

BURIED_RSA = 0.25
CONFIDENT_PLDDT = 70.0
LOW_PLDDT = 50.0
CONSERVED_FRACTION = 0.9
VARIABLE_FRACTION = 0.6
PHYLOP_CONSERVED = 2.0
TOLERATED_IN_ORTHOLOGS = 0.05

# Published thresholds: AlphaMissense (Cheng et al. 2023) and the ClinGen
# REVEL calibration (Pejaver et al. 2022, supporting level).
ALPHAMISSENSE_PATHOGENIC = 0.564
ALPHAMISSENSE_BENIGN = 0.34
REVEL_PATHOGENIC = 0.644
REVEL_BENIGN = 0.290

FEW_BENIGN = 10

DISCLAIMER = ("Research use only. Model-support tiers are not ACMG/AMP classifications "
              "and must not be used for diagnosis or clinical decisions.")


def balanced_probability(p: float, prevalence: float) -> float:
    """
    Re-express a probability calibrated at the training prevalence (85%
    pathogenic) as if pathogenic and benign were equally likely a priori.
    The calibrated value would call most VUS pathogenic simply because most
    training variants were; the balanced value reflects only the evidence.
    """
    p = min(max(p, 1e-6), 1 - 1e-6)
    odds = p / (1 - p) * (1 - prevalence) / prevalence
    return odds / (1 + odds)


def tier(probability: float) -> str:
    return next(name for cutoff, name in TIERS if probability >= cutoff)


def _missing(value) -> bool:
    return value is None or (isinstance(value, float) and math.isnan(value))


def ddg_category(ddg) -> str:
    if _missing(ddg):
        return "not modelled"
    if ddg >= DDG_STRONG:
        return "strongly destabilising"
    if ddg >= DDG_DESTABILISING:
        return "destabilising"
    if ddg <= DDG_STABILISING:
        return "stabilising"
    return "neutral"


@dataclass
class Evidence:
    name: str
    direction: str  # "pathogenic" | "benign" | "neutral" | "unavailable"
    detail: str
    independent: bool = True  # False for predictors trained on ClinVar


def conservation_evidence(r: dict) -> Evidence:
    frac, mut_seen, phylop = r.get("conservation_score"), r.get("mut_aa_in_orthologs"), r.get("phylop")
    if _missing(frac):
        return Evidence("Conservation", "unavailable", "no ortholog alignment")
    detail = f"wild-type residue kept in {frac:.0%} of vertebrate orthologs"
    if not _missing(mut_seen) and mut_seen > 0:
        detail += f", mutant residue seen in {mut_seen:.0%}"
    if not _missing(phylop):
        detail += f"; phyloP {phylop:.2f}"
    if not _missing(mut_seen) and mut_seen >= TOLERATED_IN_ORTHOLOGS:
        return Evidence("Conservation", "benign", detail + " (mutant residue tolerated in other species)")
    if frac >= CONSERVED_FRACTION and (_missing(phylop) or phylop >= PHYLOP_CONSERVED):
        return Evidence("Conservation", "pathogenic", detail + " (highly conserved)")
    if frac < VARIABLE_FRACTION:
        return Evidence("Conservation", "benign", detail + " (variable position)")
    return Evidence("Conservation", "neutral", detail)


def stability_evidence(r: dict) -> Evidence:
    ddg = r.get("foldx_ddg")
    category = ddg_category(ddg)
    if category == "not modelled":
        plddt = r.get("plddt")
        why = ("no AlphaFold model for this protein" if _missing(plddt)
               else f"residue pLDDT {plddt:.0f} is below {CONFIDENT_PLDDT:.0f}")
        return Evidence("Stability (FoldX)", "unavailable", f"not modelled: {why}")
    detail = f"ddG {ddg:+.2f} kcal/mol, {category}"
    if category in ("strongly destabilising", "destabilising"):
        return Evidence("Stability (FoldX)", "pathogenic", detail)
    # A neutral or stabilising ddG does not argue for benign: many
    # pathogenic variants act through function, not folding.
    return Evidence("Stability (FoldX)", "neutral", detail)


def structure_evidence(r: dict) -> Evidence:
    plddt, rsa = r.get("plddt"), r.get("rsa")
    if _missing(plddt):
        return Evidence("Structural context", "unavailable", "no AlphaFold model")
    parts = [f"pLDDT {plddt:.0f}"]
    if not _missing(rsa):
        parts.append(f"{rsa:.0%} solvent accessible")
    if isinstance(r.get("sec_struct"), str) and r["sec_struct"] != "unknown":
        parts.append(r["sec_struct"])
    if isinstance(r.get("domain_type"), str) and r["domain_type"]:
        parts.append(f"in {r['domain_type']} domain")
    if r.get("in_functional_site") in (True, "True"):
        parts.append("at an annotated functional site")
    detail = ", ".join(parts)
    if plddt < LOW_PLDDT or r.get("in_disordered_region") in (True, "True"):
        return Evidence("Structural context", "benign", detail + " (disordered or low-confidence region)")
    if plddt >= CONFIDENT_PLDDT and not _missing(rsa) and rsa < BURIED_RSA:
        return Evidence("Structural context", "pathogenic", detail + " (buried in a well-folded region)")
    return Evidence("Structural context", "neutral", detail)


def predictor_evidence(name: str, score, high: float, low: float) -> Evidence:
    if _missing(score):
        return Evidence(name, "unavailable", "no score", independent=False)
    direction = "pathogenic" if score >= high else "benign" if score <= low else "neutral"
    return Evidence(name, direction, f"{score:.3f} (thresholds {low} / {high})", independent=False)


def gather_evidence(r: dict) -> list[Evidence]:
    return [
        conservation_evidence(r),
        stability_evidence(r),
        structure_evidence(r),
        predictor_evidence("AlphaMissense", r.get("alphamissense"), ALPHAMISSENSE_PATHOGENIC, ALPHAMISSENSE_BENIGN),
        predictor_evidence("REVEL", r.get("revel"), REVEL_PATHOGENIC, REVEL_BENIGN),
    ]


def conclusion(probability: float, evidence: list[Evidence], clinvar_label: str | None,
               gene_benign: int) -> tuple[str, list[str]]:
    """Rule-based conclusion sentence plus a list of flags."""
    predicted = "pathogenic" if probability >= 0.5 else "benign"
    own = [e for e in evidence if e.independent and e.direction in ("pathogenic", "benign")]
    agree = [e for e in own if e.direction == predicted]
    text = f"{tier(probability).capitalize()} (P = {probability:.2f})."
    if own:
        text += (f" Independent evidence agreeing: {len(agree)} of {len(own)}"
                 f" ({', '.join(e.name.lower() for e in agree) or 'none'}).")
    flags = []
    for e in own:
        if e.direction != predicted:
            flags.append(f"{e.name} points towards {e.direction}: {e.detail}.")
    for e in evidence:
        if not e.independent and e.direction in ("pathogenic", "benign") and e.direction != predicted:
            flags.append(f"{e.name} disagrees ({e.direction}, {e.detail}).")
    if gene_benign < FEW_BENIGN:
        flags.append(f"Only {gene_benign} benign training variants in this gene, so the model has "
                     f"little evidence of what tolerated changes look like here.")
    if any(e.name == "Stability (FoldX)" and e.direction == "unavailable" for e in evidence):
        flags.append("No stability estimate; the structural part of the prediction rests on imputed values.")
    if clinvar_label in ("pathogenic", "benign") and clinvar_label != predicted:
        flags.append(f"Disagrees with the ClinVar classification ({clinvar_label}).")
    return text, flags


# ------------------------------------------------------- report object

@dataclass
class IntegratedReport:
    variant: str
    gene_symbol: str
    protein_change: str
    condition: str
    panel: str
    clinvar: dict
    prediction: dict
    tier: str
    key_features: list
    stability: dict
    structure: dict
    evidence: list
    conclusion: str
    flags: list = field(default_factory=list)
    disclaimer: str = DISCLAIMER

    def to_markdown(self) -> str:
        p = self.prediction
        lines = [
            f"# {self.gene_symbol} {self.protein_change}",
            "",
            f"**Conclusion:** {self.conclusion}",
            "",
            f"- Gene: {self.gene_symbol} ({self.panel}; {self.condition})",
            f"- ClinVar: {self.clinvar['label']} ({self.clinvar['review_stars']} star(s), "
            f"ClinVar UID {self.clinvar['clinvar_uid']})",
            f"- P(pathogenic): {p['probability']:.3f} ({p['source']})"
            + ("" if p.get("calibrated_at_training_prevalence") is None else
               f"; {p['calibrated_at_training_prevalence']:.3f} at the training prevalence"),
            "",
            "## Why the model decided this",
            "",
            "| Feature | Value | SHAP (log-odds) |",
            "| --- | --- | --- |",
        ]
        lines += [f"| {c['label']} | {c['value_display']} | {c['shap']:+.2f} towards {c['direction']} |"
                  for c in self.key_features]
        lines += ["", "## Evidence lines", "", "| Line | Points to | Detail |", "| --- | --- | --- |"]
        lines += [f"| {e['name']}{'' if e['independent'] else ' (reference)'} | {e['direction']} | {e['detail']} |"
                  for e in self.evidence]
        if self.flags:
            lines += ["", "## Flags", ""] + [f"- {f}" for f in self.flags]
        lines += ["", f"_{self.disclaimer}_", ""]
        return "\n".join(lines)


# ------------------------------------------------------------- sources

@lru_cache(maxsize=1)
def _sources():
    labelled = pd.read_csv(FEATURES_DIR / "feature_matrix.csv")
    vus_path = FEATURES_DIR / "feature_matrix_vus.csv"
    vus = pd.read_csv(vus_path) if vus_path.exists() else labelled.iloc[0:0]
    oof = pd.read_csv(OUTPUTS_DIR / "oof_predictions.csv")[["clinvar_uid", "xgboost_gene_grouped"]]
    labelled = labelled.merge(oof, on="clinvar_uid", how="left")
    qc = pd.read_csv(PROCESSED_DIR / "qc_per_gene.csv").set_index("gene_symbol")
    genes = {g["symbol"]: g for g in load_gene_panel("both")}
    return labelled, vus, qc, genes


def training_prevalence() -> float:
    qc = _sources()[2]
    return qc["pathogenic"].sum() / qc["labelled_total"].sum()


@lru_cache(maxsize=1)
def _models():
    from src.explainability.shap_explain import load_model

    return load_model()


def _scored(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, float]:
    """Probabilities and SHAP values for a block of feature rows."""
    from src.explainability.shap_explain import compute_shap, predict
    from src.modeling.train_classifier import FEATURE_COLUMNS

    model, calibrated = _models()
    X = df[FEATURE_COLUMNS]
    shap_values, base = compute_shap(model, X)
    return predict(model, calibrated, X), shap_values, base


def _build(row: dict, probability: float, source: str, contributions: list, qc, genes,
           calibrated: float | None = None) -> IntegratedReport:
    gene = row["gene_symbol"]
    change = f"{row['wt_aa']}{int(row['position'])}{row['mut_aa']}"
    evidence = gather_evidence(row)
    gene_benign = int(qc.loc[gene, "benign"]) if gene in qc.index else 0
    label = row.get("label")
    text, flags = conclusion(probability, evidence, label, gene_benign)
    ddg = row.get("foldx_ddg")
    return IntegratedReport(
        variant=f"{gene}:{change}",
        gene_symbol=gene,
        protein_change=f"p.{change}",
        condition=genes.get(gene, {}).get("condition", ""),
        panel=row.get("panel", ""),
        clinvar={"label": label, "review_stars": int(row.get("review_stars") or 0),
                 "clinvar_uid": row.get("clinvar_uid")},
        prediction={"probability": round(probability, 4), "source": source,
                    "label": "pathogenic" if probability >= 0.5 else "benign",
                    "calibrated_at_training_prevalence": None if calibrated is None else round(calibrated, 4)},
        tier=tier(probability),
        key_features=contributions,
        stability={"foldx_ddg": None if _missing(ddg) else round(float(ddg), 3), "category": ddg_category(ddg)},
        structure={k: (None if _missing(row.get(k)) else row.get(k))
                   for k in ("plddt", "rsa", "contact_count", "sec_struct", "domain_type")},
        evidence=[asdict(e) for e in evidence],
        conclusion=text,
        flags=flags,
    )


def report_for_variant(spec: str) -> IntegratedReport:
    """Integrated report for one variant, e.g. report_for_variant('MYH7:R403Q')."""
    from src.explainability.shap_explain import explain_row

    gene, change = spec.split(":")
    wt, pos, mut = change[0].upper(), int(change[1:-1]), change[-1].upper()
    labelled, vus, qc, genes = _sources()
    for df, is_labelled in ((labelled, True), (vus, False)):
        hit = df[(df["gene_symbol"] == gene) & (df["position"] == pos) & (df["wt_aa"] == wt) & (df["mut_aa"] == mut)]
        if not hit.empty:
            break
    else:
        raise ValueError(f"{spec} is not in the labelled or VUS tables")
    hit = hit.iloc[[0]].reset_index(drop=True)
    probs, shap_values, base = _scored(hit)
    row = hit.iloc[0].to_dict()
    from src.modeling.train_classifier import FEATURE_COLUMNS
    contributions = explain_row(hit[FEATURE_COLUMNS].iloc[0], shap_values.iloc[0], base)
    if is_labelled and not _missing(row.get("xgboost_gene_grouped")):
        return _build(row, float(row["xgboost_gene_grouped"]),
                      "class-weighted XGBoost, out-of-fold with this gene held out", contributions, qc, genes)
    calibrated = float(probs.iloc[0]["calibrated_probability"])
    return _build(row, balanced_probability(calibrated, training_prevalence()),
                  "calibrated XGBoost, balanced prior", contributions, qc, genes, calibrated)


def save(report: IntegratedReport) -> None:
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    stem = report.variant.replace(":", "_")
    (REPORTS_DIR / f"{stem}.md").write_text(report.to_markdown(), encoding="utf-8")
    with open(REPORTS_DIR / f"{stem}.json", "w", encoding="utf-8") as f:
        json.dump(asdict(report), f, indent=2, default=str)


def vus_report(top_n: int = 50) -> pd.DataFrame:
    """Score every VUS, write vus_report.csv and a ranked Markdown report of the top N."""
    from src.explainability.shap_explain import explain_row
    from src.modeling.train_classifier import FEATURE_COLUMNS

    _, vus, qc, genes = _sources()
    probs, shap_values, base = _scored(vus)
    prevalence = training_prevalence()
    rows = []
    reports = {}
    for i in vus.index:
        row = vus.loc[i].to_dict()
        contributions = explain_row(vus.loc[i, FEATURE_COLUMNS], shap_values.loc[i], base, top_k=3)
        calibrated = float(probs.loc[i, "calibrated_probability"])
        report = _build(row, balanced_probability(calibrated, prevalence), "calibrated XGBoost, balanced prior",
                        contributions, qc, genes, calibrated)
        own = [e for e in report.evidence if e["independent"] and e["direction"] in ("pathogenic", "benign")]
        rows.append({
            "variant": report.variant,
            "gene_symbol": report.gene_symbol,
            "clinvar_uid": row.get("clinvar_uid"),
            "probability": report.prediction["probability"],
            "calibrated_probability": round(calibrated, 4),
            "model_probability": round(float(probs.loc[i, "model_probability"]), 4),
            "tier": report.tier,
            "foldx_ddg": report.stability["foldx_ddg"],
            "ddg_category": report.stability["category"],
            "evidence_agree": sum(e["direction"] == report.prediction["label"] for e in own),
            "evidence_informative": len(own),
            "alphamissense": row.get("alphamissense"),
            "top_reasons": " | ".join(f"{c['label']} ({c['shap']:+.2f})" for c in contributions),
            "n_flags": len(report.flags),
        })
        reports[report.variant] = report
    table = pd.DataFrame(rows).sort_values(["probability", "evidence_agree", "model_probability"], ascending=False)
    table.insert(0, "rank", range(1, len(table) + 1))
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    table.to_csv(REPORTS_DIR / "vus_report.csv", index=False)

    top = table.head(top_n)
    counts = table["tier"].value_counts()
    lines = [
        f"# Top {top_n} VUS by model support for pathogenicity",
        "",
        f"{len(table)} VUS scored. Tier counts: "
        + ", ".join(f"{name}: {counts.get(name, 0)}" for _, name in TIERS) + ".",
        "",
        "P is the calibrated probability re-expressed for a 50/50 prior; the calibrated value "
        f"itself assumes the training mix ({prevalence:.0%} pathogenic) and would overstate "
        "pathogenicity for VUS. Ranked by P, then by the number of independent evidence lines "
        "(conservation, FoldX stability, structural context) that agree. Many top variants share "
        "the calibration ceiling, so treat ranks within the top tier as ties.",
        "",
        "| Rank | Variant | P | FoldX ddG | Agree | AlphaMissense | Top reasons |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for r in top.itertuples():
        ddg = "n/a" if _missing(r.foldx_ddg) else f"{r.foldx_ddg:+.2f} ({r.ddg_category})"
        am = "n/a" if _missing(r.alphamissense) else f"{r.alphamissense:.2f}"
        lines.append(f"| {r.rank} | {r.variant} | {r.probability:.2f} | {ddg} | "
                     f"{r.evidence_agree}/{r.evidence_informative} | {am} | {r.top_reasons} |")
    lines += ["", "---", ""]
    for r in top.itertuples():
        lines.append("\n".join("#" + line if line.startswith("#") else line
                               for line in reports[r.variant].to_markdown().splitlines()))
        lines += ["---", ""]
    lines.append(f"_{DISCLAIMER}_")
    (REPORTS_DIR / f"vus_top{top_n}.md").write_text("\n".join(lines), encoding="utf-8")
    return table


EXAMPLES = ["MYH7:R403Q", "KCNQ1:A341V", "LMNA:R190W", "SCN5A:R1623Q"]


# ------------------------------------------- legacy summary (used by api/ and web/)

@dataclass
class VariantReport:
    gene_symbol: str
    variant: str  # e.g. "p.Arg403Gln"
    predicted_pathogenicity: str  # "Pathogenic" | "Benign" | "Uncertain (VUS)"
    confidence_score: float
    key_features: list  # from SHAP, e.g. [{"feature": ..., "shap_value": ...}, ...]
    structural_rmsd: float | None
    predicted_ddg: float | None
    stability_interpretation: str
    plain_language_summary: str


def build_plain_language_summary(report: "VariantReport") -> str:
    top_feature = report.key_features[0]["feature"] if report.key_features else "unavailable"
    text = (
        f"The model classifies {report.gene_symbol} {report.variant} as "
        f"{report.predicted_pathogenicity} with {report.confidence_score:.0%} confidence "
        f"({tier(report.confidence_score)}). "
        f"The strongest contributing feature was '{top_feature}'."
    )
    if report.predicted_ddg is not None:
        text += (f" FoldX predicts a stability change of {report.predicted_ddg:+.2f} kcal/mol "
                 f"({ddg_category(report.predicted_ddg)}).")
    else:
        text += f" Structural analysis indicates: {report.stability_interpretation}."
    return text


def make_report(
    gene_symbol: str,
    variant: str,
    predicted_label: str,
    confidence_score: float,
    key_features: list,
    structural_rmsd: float | None,
    predicted_ddg: float | None,
    stability_interpretation: str,
) -> VariantReport:
    report = VariantReport(
        gene_symbol=gene_symbol,
        variant=variant,
        predicted_pathogenicity=predicted_label,
        confidence_score=confidence_score,
        key_features=key_features,
        structural_rmsd=structural_rmsd,
        predicted_ddg=predicted_ddg,
        stability_interpretation=stability_interpretation,
        plain_language_summary="",
    )
    report.plain_language_summary = build_plain_language_summary(report)
    return report


def main():
    parser = argparse.ArgumentParser(description="Integrated variant interpretation reports")
    parser.add_argument("--variant", action="append", help="e.g. MYH7:R403Q (repeatable)")
    parser.add_argument("--top-vus", type=int, help="write the ranked VUS report for the top N")
    args = parser.parse_args()
    ensure_dirs()

    variants = args.variant or ([] if args.top_vus else EXAMPLES)
    for spec in variants:
        try:
            report = report_for_variant(spec)
        except ValueError as e:
            print(f"[skip] {e}")
            continue
        save(report)
        print(report.to_markdown())
    if args.top_vus or not args.variant:
        n = args.top_vus or 50
        table = vus_report(n)
        print(f"\nVUS tiers:\n{table['tier'].value_counts().to_string()}")
        print(f"Saved {REPORTS_DIR / 'vus_report.csv'} and {REPORTS_DIR / f'vus_top{n}.md'}")


if __name__ == "__main__":
    main()
