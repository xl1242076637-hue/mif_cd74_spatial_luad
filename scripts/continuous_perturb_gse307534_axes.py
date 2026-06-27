#!/usr/bin/env python
"""Run continuous spatial-coupling perturbations for candidate axes in GSE307534."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]
SRC_DIR = PROJECT_ROOT / "src"
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
for path in (SRC_DIR, SCRIPTS_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from luad_niche.expression import normalize_log1p_counts_by_totals  # noqa: E402
from luad_niche.h5ad import compute_panel_scores  # noqa: E402
from luad_niche.perturbation import apply_expression_perturbation  # noqa: E402
from luad_niche.spatial_axis import build_axis_gene_panels  # noqa: E402
from luad_niche.spatial_coupling import continuous_spatial_coupling, summarize_continuous_effects  # noqa: E402
from luad_niche.spatial_niche import median_nearest_neighbor_distance  # noqa: E402
from luad_niche.tenx import discover_10x_samples, read_10x_obs, read_10x_selected_genes_with_totals  # noqa: E402
from score_gse307534_refined_signatures import (  # noqa: E402
    finite_coordinate_mask,
    flatten_genes,
    load_signatures,
    sample_accessions_from_filelist,
    sample_metadata,
)
from virtual_perturb_gse307534_axes import (  # noqa: E402
    DEFAULT_FOCUS_GENES,
    build_perturbation_specs,
    load_axes,
    parse_csv_values,
    parse_factors,
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
    parser.add_argument("--radius-multiplier", type=float, default=1.0)
    parser.add_argument("--factors", default="0,0.5")
    parser.add_argument("--focus-genes", default=",".join(DEFAULT_FOCUS_GENES))
    parser.add_argument("--limit-samples", type=int, default=None)
    return parser.parse_args()


def relative_delta(perturbed: float, baseline: float) -> float:
    if math.isnan(perturbed) or math.isnan(baseline) or baseline == 0:
        return float("nan")
    return (perturbed - baseline) / baseline


def main() -> int:
    args = parse_args()
    axes = load_axes(args.mechanisms)
    focus_genes = parse_csv_values(args.focus_genes)
    factors = parse_factors(args.factors)
    perturbation_specs = build_perturbation_specs(axes, focus_genes, factors)
    signatures = load_signatures(args.signatures)
    if "epithelial_progenitor_like" not in signatures:
        raise SystemExit("epithelial_progenitor_like signature is required.")
    panels = {"epithelial_progenitor_like": signatures["epithelial_progenitor_like"]}
    panels.update(build_axis_gene_panels(axes))
    all_genes = flatten_genes(panels)
    for spec in perturbation_specs:
        all_genes.extend(gene for gene in spec["genes"] if gene not in all_genes)

    metadata = pd.read_csv(args.metadata)
    accession_by_sample = sample_accessions_from_filelist(args.filelist)
    samples = list(discover_10x_samples(args.input_dir, require_tissue_positions=True).values())
    samples = sorted(samples, key=lambda sample: sample.sample_name)
    if args.limit_samples:
        samples = samples[: args.limit_samples]
    if not samples:
        raise SystemExit(f"No GSE307534 Visium samples found under {args.input_dir}")

    axis_by_id = {axis.get("id", ""): axis for axis in axes}
    evidence_specs = [
        ("source", "source_near_epithelial_progenitor", "source_genes"),
        ("target", "target_near_epithelial_progenitor", "target_genes"),
    ]
    effect_rows = []
    genes_used = {}
    for sample_index, sample in enumerate(samples, start=1):
        sample_accession = accession_by_sample.get(sample.sample_name, sample.sample_accession)
        meta = sample_metadata(metadata, sample_accession)
        obs = read_10x_obs(sample)
        counts, total_counts = read_10x_selected_genes_with_totals(sample, all_genes)
        normalized = normalize_log1p_counts_by_totals(counts, total_counts, scale_factor=args.scale_factor)
        baseline_scores = compute_panel_scores(normalized, panels)
        genes_used[sample_accession] = baseline_scores.attrs["panel_genes_used"]
        baseline_table = obs.merge(baseline_scores.reset_index(names="spot"), on="spot", how="left")
        baseline_table = baseline_table[baseline_table["in_tissue"] == 1].copy()
        finite = finite_coordinate_mask(baseline_table[["x", "y"]].to_numpy(dtype=float))
        baseline_table = baseline_table.loc[finite].copy()
        coords = baseline_table[["x", "y"]].to_numpy(dtype=float)
        if len(coords) < 2:
            continue
        radius = args.radius_multiplier * median_nearest_neighbor_distance(coords)

        for spec in perturbation_specs:
            axis = axis_by_id.get(str(spec["axis_id"]), {})
            perturbed_expr, perturbation_metadata = apply_expression_perturbation(
                normalized,
                list(spec["genes"]),
                factor=float(spec["factor"]),
            )
            perturbed_scores = compute_panel_scores(perturbed_expr, panels)
            perturbed_table = baseline_table[["spot", "x", "y", "in_tissue"]].merge(
                perturbed_scores.reset_index(names="spot"),
                on="spot",
                how="left",
            )
            for suffix, evidence_type, gene_field in evidence_specs:
                panel_name = f"{spec['axis_id']}_{suffix}"
                panel_column = f"{panel_name}_score"
                panel_genes = set(panels.get(panel_name, []))
                if panel_column not in baseline_table.columns or panel_column not in perturbed_table.columns:
                    continue
                if not set(spec["genes"]) & panel_genes:
                    continue
                baseline = continuous_spatial_coupling(
                    baseline_table,
                    source_column="epithelial_progenitor_like_score",
                    target_column=panel_column,
                    radius=radius,
                )
                perturbed = continuous_spatial_coupling(
                    perturbed_table,
                    source_column="epithelial_progenitor_like_score",
                    target_column=panel_column,
                    radius=radius,
                )
                coupling_delta = perturbed["coupling_score"] - baseline["coupling_score"]
                weighted_delta = (
                    perturbed["source_weighted_neighbor_target_mean"]
                    - baseline["source_weighted_neighbor_target_mean"]
                )
                effect_rows.append(
                    {
                        "sample": sample_accession,
                        "sample_name": sample.sample_name,
                        "stage": meta["stage"],
                        "axis_id": spec["axis_id"],
                        "axis_label": spec["axis_label"],
                        "evidence_type": evidence_type,
                        "gene_panel": suffix,
                        "panel_genes": ",".join(axis.get(gene_field, []) or []),
                        "perturbation_id": spec["perturbation_id"],
                        "perturbation_type": spec["perturbation_type"],
                        "perturbed_genes": ",".join(spec["genes"]),
                        "perturbation_factor": spec["factor"],
                        "present_perturbed_genes": ",".join(perturbation_metadata["present_genes"]),
                        "missing_perturbed_genes": ",".join(perturbation_metadata["missing_genes"]),
                        "radius": radius,
                        "baseline_coupling_score": baseline["coupling_score"],
                        "perturbed_coupling_score": perturbed["coupling_score"],
                        "coupling_delta": coupling_delta,
                        "coupling_relative_delta": relative_delta(
                            perturbed["coupling_score"],
                            baseline["coupling_score"],
                        ),
                        "baseline_source_weighted_neighbor_target_mean": baseline[
                            "source_weighted_neighbor_target_mean"
                        ],
                        "perturbed_source_weighted_neighbor_target_mean": perturbed[
                            "source_weighted_neighbor_target_mean"
                        ],
                        "source_weighted_neighbor_target_delta": weighted_delta,
                        "source_weighted_neighbor_target_relative_delta": relative_delta(
                            perturbed["source_weighted_neighbor_target_mean"],
                            baseline["source_weighted_neighbor_target_mean"],
                        ),
                        "baseline_source_mean": baseline["source_mean"],
                        "perturbed_source_mean": perturbed["source_mean"],
                        "baseline_target_mean": baseline["target_mean"],
                        "perturbed_target_mean": perturbed["target_mean"],
                        "n_spots": baseline["n_spots"],
                        "n_spots_with_neighbors": baseline["n_spots_with_neighbors"],
                    }
                )
        print(f"{sample_index}/{len(samples)} {sample_accession} {sample.sample_name} {meta['stage']}")

    args.table_dir.mkdir(parents=True, exist_ok=True)
    effects = pd.DataFrame(effect_rows)
    summary = summarize_continuous_effects(effects)
    metadata_columns = effects[
        ["perturbation_id", "perturbation_type", "perturbed_genes", "perturbation_factor"]
    ].drop_duplicates()
    summary = summary.merge(metadata_columns, on="perturbation_id", how="left")
    leading_columns = [
        "perturbation_id",
        "perturbation_type",
        "perturbed_genes",
        "perturbation_factor",
        "axis_id",
        "evidence_type",
    ]
    summary = summary[leading_columns + [column for column in summary.columns if column not in leading_columns]]
    effects.to_csv(
        args.table_dir / "gse307534_continuous_perturbation_effects.csv",
        index=False,
        encoding="utf-8-sig",
    )
    summary.to_csv(
        args.table_dir / "gse307534_continuous_perturbation_mia_luad_ranking.csv",
        index=False,
        encoding="utf-8-sig",
    )
    (args.table_dir / "gse307534_continuous_perturbation_genes_used.json").write_text(
        json.dumps({"genes_used_by_sample": genes_used, "perturbation_specs": perturbation_specs}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(f"Wrote {len(effects)} continuous perturbation effect rows.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
