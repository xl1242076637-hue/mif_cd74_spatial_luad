#!/usr/bin/env python
"""Run score-level in-silico target-prioritization perturbations for GSE307534 axes."""

from __future__ import annotations

import argparse
import json
import math
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
from luad_niche.perturbation import apply_expression_perturbation, summarize_perturbation_effects  # noqa: E402
from luad_niche.spatial_axis import build_axis_gene_panels  # noqa: E402
from luad_niche.spatial_niche import (  # noqa: E402
    finite_coordinate_mask,
    median_nearest_neighbor_distance,
    nearest_target_fraction,
    top_quantile_mask,
)
from luad_niche.tenx import discover_10x_samples, read_10x_obs, read_10x_selected_genes_with_totals  # noqa: E402
from score_gse307534_refined_signatures import (  # noqa: E402
    flatten_genes,
    load_signatures,
    sample_accessions_from_filelist,
    sample_metadata,
)


DEFAULT_FOCUS_GENES = ["MIF", "CD74", "CD44", "CXCR4", "SPP1", "TREM2", "PLA2G7"]


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
    parser.add_argument("--factors", default="0,0.5")
    parser.add_argument("--focus-genes", default=",".join(DEFAULT_FOCUS_GENES))
    parser.add_argument("--limit-samples", type=int, default=None)
    return parser.parse_args()


def load_axes(path: Path) -> list[dict]:
    with Path(path).open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    return list(config.get("axes", []))


def parse_csv_values(text: str) -> list[str]:
    return [value.strip() for value in text.split(",") if value.strip()]


def parse_factors(text: str) -> list[float]:
    return [float(value) for value in parse_csv_values(text)]


def build_perturbation_specs(axes: list[dict], focus_genes: list[str], factors: list[float]) -> list[dict[str, object]]:
    focus = set(focus_genes)
    specs: list[dict[str, object]] = []
    seen: set[tuple] = set()
    for axis in axes:
        axis_id = axis.get("id", "")
        perturbation_genes = [gene for gene in axis.get("perturbation_genes", []) if gene in focus]
        for genes, perturbation_type in (
            *[([gene], "gene") for gene in perturbation_genes],
            (perturbation_genes, "axis"),
        ):
            if not genes:
                continue
            if perturbation_type == "axis" and len(genes) < 2:
                continue
            for factor in factors:
                gene_label = "_".join(genes)
                factor_label = f"{factor:g}".replace(".", "p")
                key = (axis_id, perturbation_type, tuple(genes), factor)
                if key in seen:
                    continue
                seen.add(key)
                specs.append(
                    {
                        "perturbation_id": f"{axis_id}_{perturbation_type}_{gene_label}_x{factor_label}",
                        "axis_id": axis_id,
                        "axis_label": axis.get("label", axis_id),
                        "perturbation_type": perturbation_type,
                        "genes": genes,
                        "factor": factor,
                    }
                )
    return specs


def observed_adjacency(
    table: pd.DataFrame,
    source_column: str,
    target_column: str,
    quantile: float,
    radius_multiplier: float,
) -> dict[str, object]:
    usable = table[table[["x", "y", source_column, target_column]].notna().all(axis=1)].copy()
    finite_mask = finite_coordinate_mask(usable[["x", "y"]].to_numpy(dtype=float))
    usable = usable.loc[finite_mask].reset_index(drop=True)
    if len(usable) < 2:
        return {"observed_fraction": float("nan"), "n_source_high": 0, "n_target_high": 0, "status": "insufficient_spots"}
    coords = usable[["x", "y"]].to_numpy(dtype=float)
    radius = radius_multiplier * median_nearest_neighbor_distance(coords)
    source_high = top_quantile_mask(usable[source_column], quantile).to_numpy()
    target_high = top_quantile_mask(usable[target_column], quantile).to_numpy()
    n_source_high = int(source_high.sum())
    n_target_high = int(target_high.sum())
    if n_source_high == 0:
        observed = float("nan")
        status = "insufficient_source_high"
    elif n_target_high == 0:
        observed = 0.0
        status = "insufficient_target_high"
    else:
        observed = nearest_target_fraction(coords, source_high, target_high, radius=radius)
        status = "ok"
    return {
        "observed_fraction": observed,
        "n_source_high": n_source_high,
        "n_target_high": n_target_high,
        "status": status,
        "radius": radius,
    }


def rank_perturbations(stage_summary: pd.DataFrame, stages: tuple[str, ...] = ("MIA", "LUAD")) -> pd.DataFrame:
    subset = stage_summary[stage_summary["stage"].isin(stages)].copy()
    if subset.empty:
        return pd.DataFrame()
    subset["panel_relative_delta"] = subset["panel_mean_delta_mean"] / subset["baseline_panel_mean"]
    ranked = (
        subset.groupby(["perturbation_id", "axis_id", "evidence_type"])
        .agg(
            n_stage_rows=("stage", "size"),
            n_samples=("n_samples", "sum"),
            baseline_panel_mean=("baseline_panel_mean", "mean"),
            panel_mean_delta=("panel_mean_delta_mean", "mean"),
            panel_relative_delta=("panel_relative_delta", "mean"),
            observed_fraction_delta=("observed_fraction_delta_mean", "mean"),
        )
        .reset_index()
    )
    ranked["dropout_priority_score"] = (
        -ranked["panel_relative_delta"].clip(upper=0).fillna(0)
        + -ranked["observed_fraction_delta"].clip(upper=0).fillna(0)
    )
    ranked = ranked.sort_values(
        ["dropout_priority_score", "panel_relative_delta", "observed_fraction_delta"],
        ascending=[False, True, True],
    ).reset_index(drop=True)
    numeric_columns = ranked.select_dtypes(include="number").columns
    ranked[numeric_columns] = ranked[numeric_columns].round(12)
    return ranked


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
        baseline_table.insert(0, "sample_name", sample.sample_name)
        baseline_table.insert(0, "sample", sample_accession)
        for key, value in meta.items():
            baseline_table[key] = value

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
                baseline_adj = observed_adjacency(
                    baseline_table,
                    source_column="epithelial_progenitor_like_score",
                    target_column=panel_column,
                    quantile=args.quantile,
                    radius_multiplier=args.radius_multiplier,
                )
                perturbed_adj = observed_adjacency(
                    perturbed_table,
                    source_column="epithelial_progenitor_like_score",
                    target_column=panel_column,
                    quantile=args.quantile,
                    radius_multiplier=args.radius_multiplier,
                )
                baseline_panel_mean = float(baseline_table[panel_column].mean())
                perturbed_panel_mean = float(perturbed_table[panel_column].mean())
                baseline_observed = float(baseline_adj["observed_fraction"])
                perturbed_observed = float(perturbed_adj["observed_fraction"])
                observed_delta = (
                    perturbed_observed - baseline_observed
                    if not (math.isnan(perturbed_observed) or math.isnan(baseline_observed))
                    else float("nan")
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
                        "baseline_panel_mean": baseline_panel_mean,
                        "perturbed_panel_mean": perturbed_panel_mean,
                        "panel_mean_delta": perturbed_panel_mean - baseline_panel_mean,
                        "baseline_observed_fraction": baseline_observed,
                        "perturbed_observed_fraction": perturbed_observed,
                        "observed_fraction_delta": observed_delta,
                        "baseline_n_target_high": baseline_adj["n_target_high"],
                        "perturbed_n_target_high": perturbed_adj["n_target_high"],
                        "baseline_status": baseline_adj["status"],
                        "perturbed_status": perturbed_adj["status"],
                    }
                )
        print(f"{sample_index}/{len(samples)} {sample_accession} {sample.sample_name} {meta['stage']}")

    args.table_dir.mkdir(parents=True, exist_ok=True)
    effects = pd.DataFrame(effect_rows)
    stage_summary = summarize_perturbation_effects(effects)
    ranking = rank_perturbations(stage_summary)
    ranking_metadata = effects[
        ["perturbation_id", "perturbation_type", "perturbed_genes", "perturbation_factor"]
    ].drop_duplicates()
    ranking = ranking.merge(ranking_metadata, on="perturbation_id", how="left")
    leading_columns = [
        "perturbation_id",
        "perturbation_type",
        "perturbed_genes",
        "perturbation_factor",
        "axis_id",
        "evidence_type",
    ]
    ranking = ranking[leading_columns + [column for column in ranking.columns if column not in leading_columns]]
    effects.to_csv(args.table_dir / "gse307534_virtual_perturbation_effects.csv", index=False, encoding="utf-8-sig")
    stage_summary.to_csv(
        args.table_dir / "gse307534_virtual_perturbation_by_stage.csv",
        index=False,
        encoding="utf-8-sig",
    )
    ranking.to_csv(
        args.table_dir / "gse307534_virtual_perturbation_mia_luad_ranking.csv",
        index=False,
        encoding="utf-8-sig",
    )
    (args.table_dir / "gse307534_virtual_perturbation_genes_used.json").write_text(
        json.dumps({"genes_used_by_sample": genes_used, "perturbation_specs": perturbation_specs}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(f"Wrote {len(effects)} perturbation effect rows.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
