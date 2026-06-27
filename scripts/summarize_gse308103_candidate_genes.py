#!/usr/bin/env python
"""Summarize focused candidate-gene expression across GSE308103 snRNA compartments."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from luad_niche.expression import normalize_log1p_counts_by_totals  # noqa: E402
from luad_niche.orthogonal_validation import (  # noqa: E402
    FOCUSED_GENE_CONTEXTS,
    summarize_compartment_gene_expression,
    summarize_stage_expression,
)
from luad_niche.sn_matrix import parse_gse308103_filename, read_tabular_selected_genes  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=PROJECT_ROOT / "data" / "interim" / "GSE308103" / "raw_counts",
    )
    parser.add_argument(
        "--assignments",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse308103_snrna_cell_state_assignments.csv",
    )
    parser.add_argument("--table-dir", type=Path, default=PROJECT_ROOT / "results" / "tables")
    parser.add_argument("--scale-factor", type=float, default=10_000.0)
    parser.add_argument("--limit-samples", type=int, default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    genes = list(FOCUSED_GENE_CONTEXTS)
    assignments = pd.read_csv(
        args.assignments,
        usecols=["sample", "sample_name", "stage", "cell_barcode", "total_counts", "broad_class"],
    )
    assignments_by_sample = {
        sample: table.set_index("cell_barcode")
        for sample, table in assignments.groupby("sample", sort=False)
    }

    files = sorted(args.input_dir.glob("GSM*.raw_counts.mtx.txt.gz"))
    if args.limit_samples:
        files = files[: args.limit_samples]
    if not files:
        raise SystemExit(f"No GSE308103 raw-count files found under {args.input_dir}")

    summaries = []
    genes_used_by_sample: dict[str, list[str]] = {}
    for index, path in enumerate(files, start=1):
        parsed = parse_gse308103_filename(path)
        sample = parsed["sample_accession"]
        sample_name = parsed["sample_name"]
        annotation = assignments_by_sample.get(sample)
        if annotation is None or annotation.empty:
            raise ValueError(f"No cell-state assignments found for {sample}")

        counts = read_tabular_selected_genes(path, genes)
        totals = pd.to_numeric(annotation["total_counts"], errors="coerce").reindex(counts.index)
        normalized = normalize_log1p_counts_by_totals(counts, totals, scale_factor=args.scale_factor)
        working = normalized.join(annotation[["broad_class"]], how="left")
        if working["broad_class"].isna().any():
            raise ValueError(f"Cell-state assignments do not cover every barcode in {sample}")
        stages = annotation["stage"].dropna().astype(str).unique().tolist()
        if len(stages) != 1:
            raise ValueError(f"Expected exactly one stage for {sample}; found {stages}")
        genes_used_by_sample[sample] = normalized.columns.astype(str).tolist()
        summaries.append(
            summarize_compartment_gene_expression(
                working,
                sample=sample,
                sample_name=sample_name,
                stage=stages[0],
                genes=genes,
            )
        )
        print(f"[{index}/{len(files)}] {sample}: stage={stages[0]}; cells={len(working)}", flush=True)

    sample_summary = pd.concat(summaries, ignore_index=True)
    stage_summary = summarize_stage_expression(sample_summary)
    args.table_dir.mkdir(parents=True, exist_ok=True)
    sample_output = args.table_dir / "gse308103_snrna_candidate_gene_sample_summary.csv"
    stage_output = args.table_dir / "gse308103_snrna_candidate_gene_stage_summary.csv"
    genes_output = args.table_dir / "gse308103_snrna_candidate_gene_genes_used.json"
    sample_summary.to_csv(sample_output, index=False, encoding="utf-8-sig")
    stage_summary.to_csv(stage_output, index=False, encoding="utf-8-sig")
    genes_output.write_text(json.dumps(genes_used_by_sample, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Wrote sample summary: {sample_output}")
    print(f"Wrote stage summary: {stage_output}")
    print(f"Wrote gene inventory: {genes_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
