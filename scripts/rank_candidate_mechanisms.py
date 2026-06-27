#!/usr/bin/env python
"""Rank candidate epithelial-myeloid mechanisms from multi-cohort evidence."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
import yaml


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from luad_niche.mechanism_ranking import rank_candidate_axes, summarize_perturbation_candidates  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--axes",
        type=Path,
        default=PROJECT_ROOT / "config" / "candidate_mechanisms.yaml",
    )
    parser.add_argument(
        "--gene-top-celltype",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse131907_selected_gene_top_celltype.csv",
    )
    parser.add_argument(
        "--bulk-trends",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse282617_candidate_marker_trends.csv",
    )
    parser.add_argument(
        "--spatial-summary",
        type=Path,
        default=PROJECT_ROOT
        / "results"
        / "tables"
        / "specificity_filtered"
        / "gse307534_refined_signature_adjacency_by_stage.csv",
    )
    parser.add_argument(
        "--snrna-summary",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse308103_snrna_refined_signature_stage_summary.csv",
    )
    parser.add_argument(
        "--scrna-summary",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse164789_scrna_refined_signature_stage_summary.csv",
    )
    parser.add_argument(
        "--axis-spatial-summary",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse307534_candidate_axis_spatial_by_stage.csv",
    )
    parser.add_argument("--table-dir", type=Path, default=PROJECT_ROOT / "results" / "tables")
    return parser.parse_args()


def load_axes(path: Path) -> list[dict]:
    with Path(path).open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    return list(data["axes"])


def main() -> int:
    args = parse_args()
    axes = load_axes(args.axes)
    gene_top = pd.read_csv(args.gene_top_celltype)
    bulk = pd.read_csv(args.bulk_trends)
    spatial = pd.read_csv(args.spatial_summary)
    snrna = pd.read_csv(args.snrna_summary)
    scrna = pd.read_csv(args.scrna_summary)
    axis_spatial = pd.read_csv(args.axis_spatial_summary) if args.axis_spatial_summary.exists() else None

    ranked = rank_candidate_axes(axes, gene_top, bulk, spatial, snrna, scrna, axis_spatial)
    perturbations = summarize_perturbation_candidates(axes, ranked)

    args.table_dir.mkdir(parents=True, exist_ok=True)
    ranked_output = args.table_dir / "candidate_mechanism_axis_ranking.csv"
    perturbation_output = args.table_dir / "candidate_perturbation_gene_ranking.csv"
    weights_output = args.table_dir / "candidate_mechanism_ranking_weights.json"
    ranked.to_csv(ranked_output, index=False, encoding="utf-8-sig")
    perturbations.to_csv(perturbation_output, index=False, encoding="utf-8-sig")
    weights_output.write_text(
        json.dumps(
            {
                "priority_score": {
                    "spatial_score": 0.20,
                    "axis_spatial_score": 0.15,
                    "source_specificity_score": 0.17,
                    "target_specificity_score": 0.12,
                    "bulk_score": 0.13,
                    "snrna_score_component": 0.13,
                    "scrna_score_component": 0.10,
                },
                "legacy_priority_score_without_axis_spatial": {
                    "spatial_score": 0.25,
                    "source_specificity_score": 0.20,
                    "target_specificity_score": 0.15,
                    "bulk_score": 0.15,
                    "snrna_score_component": 0.15,
                    "scrna_score_component": 0.10,
                },
                "normalization": {
                    "spatial_score": "positive spatial MIA/LUAD mean delta divided by 0.10, clipped to 0-1",
                    "axis_spatial_score": "positive candidate-axis source/target MIA/LUAD adjacency mean delta divided by 0.20, clipped to 0-1",
                    "bulk_score": "positive mean bulk IAC-vs-normal delta divided by 20, clipped to 0-1",
                    "snrna_score_component": "positive late-vs-early score delta divided by 0.10, clipped to 0-1",
                    "scrna_score_component": "positive tumor-vs-adjacent score delta divided by 0.05, clipped to 0-1",
                },
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"Wrote axis ranking: {ranked_output}")
    print(f"Wrote perturbation ranking: {perturbation_output}")
    print(f"Wrote scoring weights: {weights_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
