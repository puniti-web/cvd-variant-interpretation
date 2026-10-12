"""
Top-level orchestrator. Runs pipeline stages in order via CLI flags, so you
can run the whole thing or resume from any stage during development.

Usage:
    python -m src.pipeline --stage all
    python -m src.pipeline --stage collect
    python -m src.pipeline --stage preprocess
    python -m src.pipeline --stage annotate
    python -m src.pipeline --stage features
    python -m src.pipeline --stage train
    python -m src.pipeline --stage explain
    python -m src.pipeline --stage report

FoldX ddG (src/structure/foldx_ddg.py) runs on separate compute and is not
a stage here; build_feature_matrix picks up data/raw/foldx_ddg.csv when present.
"""
import argparse

from src.config import ensure_dirs


def run_collect():
    from src.data_collection import fetch_clinvar, fetch_uniprot, fetch_structures

    fetch_clinvar.main()
    fetch_uniprot.main()
    fetch_structures.main()


def run_preprocess():
    from src.preprocessing import clean_variants

    clean_variants.main()


def run_annotate():
    import sys

    from src.data_collection import fetch_ensembl_vep
    from src.features import build_conservation_scores, structure_features

    build_conservation_scores.main()
    structure_features.main()
    # VEP reads its own CLI flags; run it on everything (labelled + VUS).
    argv, sys.argv = sys.argv, [sys.argv[0]]
    try:
        fetch_ensembl_vep.main()
    finally:
        sys.argv = argv


def run_features():
    from src.features import build_feature_matrix

    build_feature_matrix.main()


def run_train():
    from src.modeling import train_classifier

    train_classifier.main()


def run_explain():
    import sys

    from src.explainability import shap_explain

    # shap_explain parses its own CLI flags.
    argv, sys.argv = sys.argv, [sys.argv[0]]
    try:
        shap_explain.main()
    finally:
        sys.argv = argv


def run_report():
    import sys

    from src.report import generate_report

    argv, sys.argv = sys.argv, [sys.argv[0]]
    try:
        generate_report.main()
    finally:
        sys.argv = argv


STAGES = {
    "collect": run_collect,
    "preprocess": run_preprocess,
    "annotate": run_annotate,
    "features": run_features,
    "train": run_train,
    "explain": run_explain,
    "report": run_report,
}


def main():
    parser = argparse.ArgumentParser(description="Cardiovascular Variant XAI pipeline")
    parser.add_argument(
        "--stage",
        choices=list(STAGES.keys()) + ["all"],
        default="all",
        help="Which stage to run (default: all, in order)",
    )
    args = parser.parse_args()

    ensure_dirs()

    if args.stage == "all":
        for name, fn in STAGES.items():
            print(f"\n{'=' * 60}\nSTAGE: {name}\n{'=' * 60}")
            fn()
    else:
        STAGES[args.stage]()


if __name__ == "__main__":
    main()
