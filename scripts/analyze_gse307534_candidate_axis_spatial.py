#!/usr/bin/env python
"""Test candidate mechanism-axis spatial adjacency in GSE307534 Visium samples."""

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
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
for path in (SRC_DIR, SCRIPTS_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from luad_niche.expression import normalize_log1p_counts_by_totals  # noqa: E402
from luad_niche.h5ad import compute_panel_scores  # noqa: E402
from luad_niche.spatial_axis import build_axis_gene_panels, summarize_axis_adjacency_by_stage  # noqa: E402
from luad_niche.tenx import discover_10x_samples, read_10x_obs, read_10x_selected_genes_with_totals  # noqa: E402
from score_gse307534_refined_signatures import (  # noqa: E402
    adjacency_for_target,
    flatten_genes,
    load_signatures,
    sample_accessions_from_filelist,
    sample_metadata,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=PROJECT_ROOT / "data" / "interim" / "GSE307534" / "raw_visium",
    )
    parser.add_argument(
        "--filelist",
        type=Path,
        default=PROJECT_ROOT / "data" / "raw" / "GSE307534" / "filelist.txt",
    )
    parser.add_argument(
        "--mechanisms",
        type=Path,
        default=PROJECT_ROOT / "config" / "candidate_mechanisms.yaml",
    )
    parser.add_argument(
        "--signatures",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse189357_refined_signature_genes_gse131907_specificity_filtered.json",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "geo_sample_metadata_annotated.csv",
    )
    parser.add_argument("--table-dir", type=Path, default=PROJECT_ROOT / "results" / "tables")
    parser.add_argument("--scale-factor", type=float, default=10_000.0)
    parser.add_argument("--quantile", type=float, default=0.75)
    parser.add_argument("--radius-multiplier", type=float, default=1.0)
    parser.add_argument("--permutations", type=int, default=100)
    parser.add_argument("--seed", type=int, default=307534)
    parser.add_argument("--limit-samples", type=int, default=None)
    return parser.parse_args()


def load_axes(path: Path) -> list[dict]:
    with Path(path).open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    return list(config.get("axes", []))


def summarize_scores_by_stage(scores: pd.DataFrame, score_columns: list[str]) -> pd.DataFrame:
    records = []
    for column in score_columns:
        sample_means = (
            scores.groupby(["sample", "sample_name", "stage"])[column]
            .mean()
            .reset_index(name="mean_score")
        )
        sample_means["panel"] = column.removesuffix("_score")
        records.append(sample_means)
    sample_long = pd.concat(records, ignore_index=True)
    return (
        sample_long.groupby(["stage", "panel"])
        .agg(n_samples=("sample", "nunique"), mean_score=("mean_score", "mean"), sd_score=("mean_score", "std"))
        .reset_index()
    )


def main() -> int:
    args = parse_args()
    axes = load_axes(args.mechanisms)
    signatures = load_signatures(args.signatures)
    if "epithelial_progenitor_like" not in signatures:
        raise SystemExit("epithelial_progenitor_like signature is required.")
    panels = {"epithelial_progenitor_like": signatures["epithelial_progenitor_like"]}
    panels.update(build_axis_gene_panels(axes))
    all_genes = flatten_genes(panels)

    metadata = pd.read_csv(args.metadata)
    accession_by_sample = sample_accessions_from_filelist(args.filelist)
    samples = list(discover_10x_samples(args.input_dir, require_tissue_positions=True).values())
    samples = sorted(samples, key=lambda sample: sample.sample_name)
    if args.limit_samples:
        samples = samples[: args.limit_samples]
    if not samples:
        raise SystemExit(f"No GSE307534 Visium samples found under {args.input_dir}")

    args.table_dir.mkdir(parents=True, exist_ok=True)
    score_tables = []
    adjacency_rows = []
    genes_used = {}
    axis_by_id = {axis.get("id", ""): axis for axis in axes}
    evidence_specs = [
        ("source", "source_near_epithelial_progenitor", "source_genes"),
        ("target", "target_near_epithelial_progenitor", "target_genes"),
    ]

    for sample_index, sample in enumerate(samples):
        sample_accession = accession_by_sample.get(sample.sample_name, sample.sample_accession)
        meta = sample_metadata(metadata, sample_accession)
        obs = read_10x_obs(sample)
        counts, total_counts = read_10x_selected_genes_with_totals(sample, all_genes)
        normalized = normalize_log1p_counts_by_totals(counts, total_counts, scale_factor=args.scale_factor)
        scores = compute_panel_scores(normalized, panels)
        genes_used[sample_accession] = scores.attrs["panel_genes_used"]
        table = obs.merge(total_counts.rename_axis("spot").reset_index(), on="spot", how="left")
        table = table.merge(scores.reset_index(names="spot"), on="spot", how="left")
        table = table[table["in_tissue"] == 1].copy()
        table.insert(0, "sample_name", sample.sample_name)
        table.insert(0, "sample", sample_accession)
        for key, value in meta.items():
            table[key] = value
        score_tables.append(table)

        for axis_index, axis_id in enumerate(axis_by_id):
            axis = axis_by_id[axis_id]
            for evidence_index, (suffix, evidence_type, gene_field) in enumerate(evidence_specs):
                panel_column = f"{axis_id}_{suffix}_score"
                if panel_column not in table.columns:
                    continue
                stats = adjacency_for_target(
                    table,
                    source_column="epithelial_progenitor_like_score",
                    target_column=panel_column,
                    quantile=args.quantile,
                    radius_multiplier=args.radius_multiplier,
                    permutations=args.permutations,
                    seed=args.seed + sample_index * 100 + axis_index * 10 + evidence_index,
                )
                stats.update(
                    {
                        "axis_id": axis_id,
                        "axis_label": axis.get("label", axis_id),
                        "evidence_type": evidence_type,
                        "gene_panel": suffix,
                        "genes": ",".join(axis.get(gene_field, []) or []),
                    }
                )
                adjacency_rows.append(stats)

    scores = pd.concat(score_tables, ignore_index=True)
    adjacency = pd.DataFrame(adjacency_rows)
    stage_summary = summarize_axis_adjacency_by_stage(adjacency)
    score_columns = ["epithelial_progenitor_like_score"] + [
        f"{axis.get('id')}_{suffix}_score"
        for axis in axes
        for suffix in ("source", "target")
        if f"{axis.get('id')}_{suffix}_score" in scores.columns
    ]
    score_stage = summarize_scores_by_stage(scores, score_columns)

    adjacency.to_csv(
        args.table_dir / "gse307534_candidate_axis_spatial_adjacency.csv",
        index=False,
        encoding="utf-8-sig",
    )
    stage_summary.to_csv(
        args.table_dir / "gse307534_candidate_axis_spatial_by_stage.csv",
        index=False,
        encoding="utf-8-sig",
    )
    score_stage.to_csv(
        args.table_dir / "gse307534_candidate_axis_spatial_score_by_stage.csv",
        index=False,
        encoding="utf-8-sig",
    )
    (args.table_dir / "gse307534_candidate_axis_spatial_genes_used.json").write_text(
        json.dumps(genes_used, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(f"Wrote candidate-axis spatial adjacency for {len(samples)} samples and {len(axes)} axes.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
