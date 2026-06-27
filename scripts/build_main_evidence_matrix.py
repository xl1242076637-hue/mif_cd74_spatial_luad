#!/usr/bin/env python
"""Build a compact manuscript-facing evidence matrix for candidate LUAD axes."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from luad_niche.evidence_matrix import (  # noqa: E402
    best_continuous_perturbation_by_axis,
    mean_axis_spatial_by_evidence,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--axis-ranking",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "candidate_mechanism_axis_ranking.csv",
    )
    parser.add_argument(
        "--axis-spatial",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse307534_candidate_axis_spatial_by_stage.csv",
    )
    parser.add_argument(
        "--continuous-perturbation",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse307534_continuous_perturbation_mia_luad_ranking.csv",
    )
    parser.add_argument("--table-dir", type=Path, default=PROJECT_ROOT / "results" / "tables")
    parser.add_argument("--docs-dir", type=Path, default=PROJECT_ROOT / "docs")
    return parser.parse_args()


def evidence_grade(row: pd.Series) -> str:
    score = float(row.get("priority_score", 0) or 0)
    perturb = float(row.get("top_continuous_priority_score", 0) or 0)
    if score >= 0.75 and perturb >= 0.7:
        return "lead"
    if score >= 0.70:
        return "strong"
    if score >= 0.55:
        return "supporting"
    return "benchmark_or_secondary"


def short_interpretation(row: pd.Series) -> str:
    axis_id = row["axis_id"]
    if axis_id == "mif_cd74_cxcr4":
        return "Lead epithelial-to-myeloid communication hypothesis; prioritize MIF and CD74."
    if axis_id == "spp1_trem2_macrophage_epithelial":
        return "Macrophage-state readout axis; prioritize SPP1, with TREM2/PLA2G7 secondary."
    if axis_id == "inflammatory_il1_tnf_cxcl8":
        return "Published inflammatory-niche benchmark/positive control, not primary novelty."
    if axis_id == "c1q_apoe_trem2_lgals3":
        return "Myeloid immunoregulatory supporting axis with mixed late/tumor support."
    if axis_id == "cxcl9_cxcl10_cxcr3":
        return "Immune-recruitment supporting axis, less central to epithelial-myeloid mechanism."
    return "Candidate axis."


def to_markdown_table(table: pd.DataFrame) -> str:
    columns = [
        "rank",
        "axis_id",
        "evidence_grade",
        "priority_score",
        "source_spatial_delta",
        "target_spatial_delta",
        "top_perturbed_genes",
        "top_coupling_relative_delta",
        "interpretation",
    ]
    display = table[columns].copy()
    for column in [
        "priority_score",
        "source_spatial_delta",
        "target_spatial_delta",
        "top_coupling_relative_delta",
    ]:
        display[column] = pd.to_numeric(display[column], errors="coerce").round(3)
    display["top_coupling_relative_delta"] = display["top_coupling_relative_delta"].astype(object)
    display.loc[display["top_perturbed_genes"].eq("not_tested"), "top_coupling_relative_delta"] = "not_tested"
    return display.to_markdown(index=False)


def main() -> int:
    args = parse_args()
    ranking = pd.read_csv(args.axis_ranking)
    spatial = pd.read_csv(args.axis_spatial)
    perturbation = pd.read_csv(args.continuous_perturbation)

    source_spatial = mean_axis_spatial_by_evidence(spatial)
    source_spatial = source_spatial.pivot_table(
        index="axis_id",
        columns="evidence_type",
        values="mean_axis_spatial_delta",
        aggfunc="first",
    ).reset_index()
    source_spatial = source_spatial.rename(
        columns={
            "source_near_epithelial_progenitor": "source_spatial_delta",
            "target_near_epithelial_progenitor": "target_spatial_delta",
        }
    )
    best_perturbation = best_continuous_perturbation_by_axis(perturbation)

    matrix = ranking.merge(source_spatial, on="axis_id", how="left").merge(
        best_perturbation,
        on="axis_id",
        how="left",
    )
    matrix = matrix.sort_values("priority_score", ascending=False).reset_index(drop=True)
    matrix.insert(0, "rank", range(1, len(matrix) + 1))
    matrix["evidence_grade"] = matrix.apply(evidence_grade, axis=1)
    matrix["interpretation"] = matrix.apply(short_interpretation, axis=1)
    for column in ["top_perturbed_genes", "top_perturbation_factor", "top_perturbation_evidence"]:
        if column in matrix.columns:
            matrix[column] = matrix[column].fillna("not_tested")

    output_columns = [
        "rank",
        "axis_id",
        "axis_label",
        "evidence_grade",
        "priority_score",
        "spatial_mia_luad_delta",
        "axis_spatial_mia_luad_delta",
        "source_spatial_delta",
        "target_spatial_delta",
        "source_specificity_fraction",
        "target_specificity_fraction",
        "bulk_delta_mean",
        "snrna_program_delta",
        "scrna_tumor_adjacent_delta",
        "top_perturbed_genes",
        "top_perturbation_factor",
        "top_perturbation_evidence",
        "top_coupling_relative_delta",
        "top_continuous_priority_score",
        "perturbation_genes",
        "interpretation",
    ]
    matrix = matrix[[column for column in output_columns if column in matrix.columns]]

    args.table_dir.mkdir(parents=True, exist_ok=True)
    args.docs_dir.mkdir(parents=True, exist_ok=True)
    table_output = args.table_dir / "main_axis_evidence_matrix.csv"
    markdown_output = args.docs_dir / "main_axis_evidence_matrix.md"
    matrix.to_csv(table_output, index=False, encoding="utf-8-sig")
    markdown = "\n".join(
        [
            "# Main Axis Evidence Matrix",
            "",
            "This table is generated from the current multi-cohort evidence outputs.",
            "",
            to_markdown_table(matrix),
            "",
            "Interpretation: lead the manuscript with `MIF-CD74`, keep `SPP1/TREM2/PLA2G7` as the macrophage-state readout, and treat `IL1B/TNF/CXCL8` as a published inflammatory benchmark.",
            "",
        ]
    )
    markdown_output.write_text(markdown, encoding="utf-8")
    print(f"Wrote evidence matrix: {table_output}")
    print(f"Wrote markdown summary: {markdown_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
