#!/usr/bin/env python
"""Build sample-level and patient-aware GSE307534 spatial-axis statistics."""

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

from luad_niche.spatial_statistics import (  # noqa: E402
    add_patient_phase,
    build_paired_patient_differences,
    summarize_late_vs_precursor,
    summarize_paired_patient_differences,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--adjacency",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse307534_candidate_axis_spatial_adjacency.csv",
    )
    parser.add_argument("--table-dir", type=Path, default=PROJECT_ROOT / "results" / "tables")
    parser.add_argument("--docs-dir", type=Path, default=PROJECT_ROOT / "docs")
    parser.add_argument("--bootstrap-iterations", type=int, default=2_000)
    parser.add_argument("--seed", type=int, default=307534)
    return parser.parse_args()


def _format_float(value: float) -> str:
    return "NA" if pd.isna(value) else f"{value:.3f}"


def write_summary_markdown(
    path: Path,
    sample_summary: pd.DataFrame,
    patient_summary: pd.DataFrame,
    paired_summary: pd.DataFrame,
) -> None:
    focus_axes = ["mif_cd74_cxcr4", "spp1_trem2_macrophage_epithelial"]
    focus = paired_summary[paired_summary["axis_id"].isin(focus_axes)].copy()
    lines = [
        "# GSE307534 Spatial Statistics Summary",
        "",
        "Date: 2026-05-30",
        "",
        "## Design",
        "",
        "- Precursor group: `AAH` and `AIS`.",
        "- Late group: `MIA` and `LUAD`.",
        "- Sample-level analysis treats each Visium section as one observation.",
        "- Patient-aggregated analysis averages repeated lesions within each patient and phase.",
        "- Paired-patient sensitivity analysis uses patients with both precursor and late lesions.",
        "- P-values are exploratory and are accompanied by Benjamini-Hochberg adjusted q-values.",
        "",
        "## Focused Paired-Patient Results",
        "",
        "| Axis | Evidence type | Paired patients | Mean late-minus-precursor delta | 95% bootstrap CI | Positive fraction | Wilcoxon q (BH) |",
        "|---|---|---:|---:|---|---:|---:|",
    ]
    for row in focus.itertuples(index=False):
        lines.append(
            "| "
            f"`{row.axis_id}` | `{row.evidence_type}` | {row.n_paired_patients} | "
            f"{_format_float(row.paired_difference_mean)} | "
            f"[{_format_float(row.ci_95_low)}, {_format_float(row.ci_95_high)}] | "
            f"{_format_float(row.positive_fraction)} | {_format_float(row.wilcoxon_q_bh)} |"
        )
    lines.extend(
        [
            "",
            "## Output Tables",
            "",
            "- `results/tables/gse307534_candidate_axis_spatial_with_patient.csv`",
            "- `results/tables/gse307534_candidate_axis_late_vs_precursor_sample_stats.csv`",
            "- `results/tables/gse307534_candidate_axis_late_vs_precursor_patient_stats.csv`",
            "- `results/tables/gse307534_candidate_axis_paired_patient_differences.csv`",
            "- `results/tables/gse307534_candidate_axis_paired_patient_stats.csv`",
            "",
            "## Interpretation Boundary",
            "",
            "These statistics strengthen the evidence that selected spatial enrichments recur across lesions and paired patients. They remain spot-level spatial-coupling summaries and do not prove direct cell-cell contact or causal pathway activity.",
            "",
            "## Full Summary Tables",
            "",
            f"- Sample-level contrasts: {len(sample_summary)} rows.",
            f"- Patient-aggregated contrasts: {len(patient_summary)} rows.",
            f"- Paired-patient contrasts: {len(paired_summary)} rows.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    adjacency = pd.read_csv(args.adjacency)
    annotated = add_patient_phase(adjacency)
    sample_summary = summarize_late_vs_precursor(
        adjacency,
        aggregate_by_patient=False,
        bootstrap_iterations=args.bootstrap_iterations,
        seed=args.seed,
    )
    patient_summary = summarize_late_vs_precursor(
        adjacency,
        aggregate_by_patient=True,
        bootstrap_iterations=args.bootstrap_iterations,
        seed=args.seed,
    )
    paired_differences = build_paired_patient_differences(adjacency)
    paired_summary = summarize_paired_patient_differences(
        paired_differences,
        bootstrap_iterations=args.bootstrap_iterations,
        seed=args.seed,
    )

    args.table_dir.mkdir(parents=True, exist_ok=True)
    args.docs_dir.mkdir(parents=True, exist_ok=True)
    annotated.to_csv(args.table_dir / "gse307534_candidate_axis_spatial_with_patient.csv", index=False)
    sample_summary.to_csv(
        args.table_dir / "gse307534_candidate_axis_late_vs_precursor_sample_stats.csv",
        index=False,
    )
    patient_summary.to_csv(
        args.table_dir / "gse307534_candidate_axis_late_vs_precursor_patient_stats.csv",
        index=False,
    )
    paired_differences.to_csv(
        args.table_dir / "gse307534_candidate_axis_paired_patient_differences.csv",
        index=False,
    )
    paired_summary.to_csv(
        args.table_dir / "gse307534_candidate_axis_paired_patient_stats.csv",
        index=False,
    )
    write_summary_markdown(
        args.docs_dir / "gse307534_spatial_statistics_summary.md",
        sample_summary,
        patient_summary,
        paired_summary,
    )
    print(
        "Wrote GSE307534 spatial statistics: "
        f"{len(sample_summary)} sample contrasts, "
        f"{len(patient_summary)} patient contrasts, "
        f"{len(paired_summary)} paired-patient contrasts."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

