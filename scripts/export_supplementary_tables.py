#!/usr/bin/env python
"""Package manuscript-facing supplementary CSV tables with a reproducible manifest."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]

TABLE_SPECS = [
    (
        "ST01",
        "Dataset composition and evidence roles",
        "figure1_dataset_composition_source.csv",
    ),
    (
        "ST02",
        "Specificity audit of refined GSE189357 discovery signatures",
        "gse189357_refined_signature_gse131907_specificity_audit.csv",
    ),
    (
        "ST03",
        "Per-sample GSE307534 candidate-axis spatial adjacency",
        "gse307534_candidate_axis_spatial_with_patient.csv",
    ),
    (
        "ST04",
        "Sample-level GSE307534 late-versus-precursor spatial contrasts",
        "gse307534_candidate_axis_late_vs_precursor_sample_stats.csv",
    ),
    (
        "ST05",
        "Patient-aggregated GSE307534 late-versus-precursor spatial contrasts",
        "gse307534_candidate_axis_late_vs_precursor_patient_stats.csv",
    ),
    (
        "ST06",
        "Paired-patient GSE307534 late-minus-precursor differences",
        "gse307534_candidate_axis_paired_patient_differences.csv",
    ),
    (
        "ST07",
        "Paired-patient GSE307534 spatial-axis statistics",
        "gse307534_candidate_axis_paired_patient_stats.csv",
    ),
    (
        "ST08",
        "Integrated multi-cohort mechanism-axis evidence matrix",
        "main_axis_evidence_matrix.csv",
    ),
    (
        "ST09",
        "Top-quantile score-level in-silico target-prioritization ranking",
        "gse307534_virtual_perturbation_mia_luad_ranking.csv",
    ),
    (
        "ST10",
        "Continuous-coupling perturbation ranking",
        "gse307534_continuous_perturbation_mia_luad_ranking.csv",
    ),
    (
        "ST11",
        "GSE282617 bulk candidate-marker trends",
        "gse282617_candidate_marker_trends.csv",
    ),
    (
        "ST12",
        "GSE308103 single-nucleus macrophage-state stage summary",
        "gse308103_snrna_macrophage_state_stage_summary.csv",
    ),
    (
        "ST13",
        "Expression-matched single-gene controls selected for source-side MIF",
        "gse307534_mif_expression_matched_controls.csv",
    ),
    (
        "ST14",
        "Paired-patient statistics for MIF and expression-matched controls",
        "gse307534_mif_random_control_paired_stats.csv",
    ),
    (
        "ST15",
        "MIF versus expression-matched control distribution summary",
        "gse307534_mif_random_control_summary.csv",
    ),
    (
        "ST16",
        "MIF tissue-density control associations",
        "gse307534_mif_density_control_summary.csv",
    ),
    (
        "ST17",
        "Broad macrophage-signature paired-patient control statistics",
        "gse307534_mif_broad_signature_control_paired_stats.csv",
    ),
    (
        "ST18",
        "Per-sample GSE308103 candidate-gene expression by broad compartment",
        "gse308103_snrna_candidate_gene_sample_summary.csv",
    ),
    (
        "ST19",
        "Stage-level GSE308103 candidate-gene expression by broad compartment",
        "gse308103_snrna_candidate_gene_stage_summary.csv",
    ),
    (
        "ST20",
        "Focused orthogonal-validation snRNA source data",
        "supplementary_figure_focused_orthogonal_validation_snrna_source.csv",
    ),
    (
        "ST21",
        "Focused orthogonal-validation bulk source data",
        "supplementary_figure_focused_orthogonal_validation_bulk_source.csv",
    ),
    (
        "ST22",
        "Covariate-adjusted source-side MIF spatial sensitivity models",
        "gse307534_mif_covariate_sensitivity_models.csv",
    ),
    (
        "ST23",
        "Paired-patient source-side MIF covariate changes",
        "gse307534_mif_covariate_paired_changes.csv",
    ),
    (
        "ST24",
        "Gene-level source data for SPP1 signature specificity refinement",
        "supplementary_figure_spp1_signature_refinement_gene_source.csv",
    ),
    (
        "ST25",
        "Specificity-status source data for SPP1 signature refinement",
        "supplementary_figure_spp1_signature_refinement_status_source.csv",
    ),
    (
        "ST26",
        "Top-celltype source data for SPP1 signature refinement",
        "supplementary_figure_spp1_signature_refinement_celltype_source.csv",
    ),
    (
        "ST27",
        "GSE308103 GRN-level virtual perturbation target ranking",
        "gse308103_grn_virtual_perturbation_target_ranking.csv",
    ),
    (
        "ST28",
        "GSE308103 GRN-level virtual perturbation robustness summary",
        "gse308103_grn_virtual_perturbation_robustness_summary.csv",
    ),
    (
        "ST29",
        "Cross-dataset validation rows for GRN-prioritized signatures",
        "grn_cross_dataset_signature_validation.csv",
    ),
    (
        "ST30",
        "Cross-dataset validation summary for GRN-prioritized signatures",
        "grn_cross_dataset_signature_validation_summary.csv",
    ),
    (
        "ST31",
        "Original scTenifoldKnk expanded sensitivity and boundary summary",
        "sctenifoldknk_original_expanded_interpretation.csv",
    ),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table-dir", type=Path, default=PROJECT_ROOT / "results" / "tables")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "results" / "supplementary_tables",
    )
    parser.add_argument("--docs-dir", type=Path, default=PROJECT_ROOT / "docs")
    return parser.parse_args()


def safe_filename(table_id: str, source_name: str) -> str:
    """Prefix a source CSV filename with its stable supplementary-table identifier."""
    return f"{table_id}_{source_name}"


def package_tables(table_dir: Path, output_dir: Path) -> pd.DataFrame:
    """Copy selected source tables and return a manifest with row/column counts."""
    output_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for table_id, description, source_name in TABLE_SPECS:
        source_path = table_dir / source_name
        if not source_path.exists():
            raise FileNotFoundError(f"Missing supplementary-table source: {source_path}")
        output_name = safe_filename(table_id, source_name)
        output_path = output_dir / output_name
        shutil.copyfile(source_path, output_path)
        table = pd.read_csv(output_path)
        records.append(
            {
                "table_id": table_id,
                "description": description,
                "source_path": source_path.relative_to(PROJECT_ROOT).as_posix(),
                "output_path": output_path.relative_to(PROJECT_ROOT).as_posix(),
                "n_rows": len(table),
                "n_columns": len(table.columns),
            }
        )
    manifest = pd.DataFrame(records)
    manifest.to_csv(output_dir / "supplementary_table_manifest.csv", index=False, encoding="utf-8-sig")
    return manifest


def write_index(path: Path, manifest: pd.DataFrame) -> None:
    """Write a human-readable supplementary-table index."""
    lines = [
        "# Supplementary Table Index",
        "",
        "This index is generated by `scripts/export_supplementary_tables.py`.",
        "",
        "| ID | Description | Rows | Columns | Packaged file |",
        "|---|---|---:|---:|---|",
    ]
    for row in manifest.itertuples(index=False):
        lines.append(
            f"| `{row.table_id}` | {row.description} | {row.n_rows} | {row.n_columns} | `{row.output_path}` |"
        )
    lines.extend(
        [
            "",
            "## Interpretation Note",
            "",
            "The perturbation tables report score-level in-silico target prioritization and GRN-level virtual perturbation prioritization. They do not represent wet-lab perturbation experiments or causal proof.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    manifest = package_tables(args.table_dir, args.output_dir)
    args.docs_dir.mkdir(parents=True, exist_ok=True)
    write_index(args.docs_dir / "supplementary_tables_index.md", manifest)
    print(f"Packaged {len(manifest)} supplementary tables under {args.output_dir}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
