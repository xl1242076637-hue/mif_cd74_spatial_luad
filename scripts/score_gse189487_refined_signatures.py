#!/usr/bin/env python
"""Map GSE189357 refined signatures back to GSE189487 spatial samples."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from luad_niche.expression import normalize_log1p_counts_by_totals  # noqa: E402
from luad_niche.h5ad import compute_panel_scores  # noqa: E402
from luad_niche.spatial_niche import (  # noqa: E402
    finite_coordinate_mask,
    median_nearest_neighbor_distance,
    nearest_target_fraction,
    permutation_adjacency_test,
    top_quantile_mask,
)
from luad_niche.tenx import discover_10x_samples, read_10x_obs, read_10x_selected_genes_with_totals  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=PROJECT_ROOT / "data" / "interim" / "GSE189487" / "raw_10x",
        help="Directory containing extracted GSE189487 spatial 10x files.",
    )
    parser.add_argument(
        "--signatures",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse189357_refined_state_signature_genes.json",
        help="Signature JSON from extract_gse189357_state_markers.py.",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "geo_sample_metadata_annotated.csv",
        help="Annotated GEO sample metadata.",
    )
    parser.add_argument(
        "--table-dir",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables",
        help="Output table directory.",
    )
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=PROJECT_ROOT / "results" / "figures",
        help="Output figure directory.",
    )
    parser.add_argument("--scale-factor", type=float, default=10_000.0)
    parser.add_argument("--quantile", type=float, default=0.75)
    parser.add_argument("--radius-multiplier", type=float, default=1.0)
    parser.add_argument("--permutations", type=int, default=500)
    parser.add_argument("--seed", type=int, default=29)
    return parser.parse_args()


def clean_signature_name(name: str) -> str:
    return (
        name.replace("_vs_other_epithelial", "")
        .replace("_vs_other_macrophage", "")
        .replace("epithelial_progenitor_like", "epithelial_progenitor_like")
    )


def load_signatures(path: Path) -> dict[str, list[str]]:
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    return {clean_signature_name(name): list(genes) for name, genes in raw.items()}


def sample_metadata(metadata: pd.DataFrame, sample_accession: str) -> dict[str, object]:
    rows = metadata[metadata["sample_accession"] == sample_accession]
    if rows.empty:
        return {"stage": "", "histological_type": "", "radiological_type": "", "gender": "", "title": ""}
    row = rows.iloc[0]
    return {
        "stage": row.get("interpreted_stage", ""),
        "histological_type": row.get("histological_type", ""),
        "radiological_type": row.get("radiological_type", ""),
        "gender": row.get("gender", ""),
        "title": row.get("title", ""),
    }


def flatten_genes(panels: dict[str, list[str]]) -> list[str]:
    genes: list[str] = []
    for panel_genes in panels.values():
        genes.extend(gene for gene in panel_genes if gene not in genes)
    return genes


def adjacency_for_target(
    table: pd.DataFrame,
    source_column: str,
    target_column: str,
    quantile: float,
    radius_multiplier: float,
    permutations: int,
    seed: int,
) -> dict[str, object]:
    usable = table[table[["x", "y", source_column, target_column]].notna().all(axis=1)].copy()
    finite_mask = finite_coordinate_mask(usable[["x", "y"]].to_numpy(dtype=float))
    usable = usable.loc[finite_mask].reset_index(drop=True)
    coords = usable[["x", "y"]].to_numpy(dtype=float)
    radius = radius_multiplier * median_nearest_neighbor_distance(coords)
    source_high = top_quantile_mask(usable[source_column], quantile).to_numpy()
    target_high = top_quantile_mask(usable[target_column], quantile).to_numpy()
    overlap = source_high & target_high
    stats = permutation_adjacency_test(
        coords,
        source_high,
        target_high,
        radius=radius,
        n_permutations=permutations,
        seed=seed,
    )
    stats.update(
        {
            "sample": usable["sample"].iloc[0],
            "sample_name": usable["sample_name"].iloc[0],
            "stage": usable["stage"].iloc[0],
            "source": source_column.removesuffix("_score"),
            "target": target_column.removesuffix("_score"),
            "quantile": quantile,
            "radius_multiplier": radius_multiplier,
            "n_spots": len(usable),
            "n_source_high": int(source_high.sum()),
            "n_target_high": int(target_high.sum()),
            "n_overlap_high": int(overlap.sum()),
            "overlap_fraction_of_source": float(overlap.sum() / source_high.sum()),
            "target_to_source_fraction": nearest_target_fraction(coords, target_high, source_high, radius=radius),
        }
    )
    return stats


def plot_stage_target_summary(stage_summary: pd.DataFrame, output: Path) -> None:
    pivot = stage_summary.pivot_table(
        index="stage",
        columns="target",
        values="enrichment_delta_mean",
        fill_value=0,
    )
    order = [stage for stage in ["AIS", "MIA", "IAC"] if stage in pivot.index]
    pivot = pivot.loc[order]
    ax = pivot.plot(kind="bar", figsize=(8, 4.8), width=0.8)
    ax.axhline(0, color="#333333", linewidth=0.8)
    ax.set_ylabel("Observed minus permuted-null adjacency")
    ax.set_xlabel("")
    ax.set_title("GSE189487 refined signature adjacency by stage")
    ax.legend(frameon=False, fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    ax.figure.savefig(output, dpi=220)
    plt.close(ax.figure)


def main() -> int:
    args = parse_args()
    signatures = load_signatures(args.signatures)
    metadata = pd.read_csv(args.metadata)
    samples = discover_10x_samples(args.input_dir, require_tissue_positions=True)
    all_genes = flatten_genes(signatures)
    score_tables = []
    adjacency_rows = []
    genes_used = {}

    for sample_index, sample in enumerate(samples.values()):
        meta = sample_metadata(metadata, sample.sample_accession)
        obs = read_10x_obs(sample)
        counts, total_counts = read_10x_selected_genes_with_totals(sample, all_genes)
        normalized = normalize_log1p_counts_by_totals(counts, total_counts, scale_factor=args.scale_factor)
        scores = compute_panel_scores(normalized, signatures)
        genes_used[sample.sample_accession] = scores.attrs["panel_genes_used"]
        table = obs.merge(total_counts.rename_axis("spot").reset_index(), on="spot", how="left")
        table = table.merge(scores.reset_index(names="spot"), on="spot", how="left")
        table = table[table["in_tissue"] == 1].copy()
        table.insert(0, "sample_name", sample.sample_name)
        table.insert(0, "sample", sample.sample_accession)
        for key, value in meta.items():
            table[key] = value
        score_tables.append(table)

        source_column = "epithelial_progenitor_like_score"
        target_columns = [
            column
            for column in table.columns
            if column.endswith("_macrophage_score") and column != source_column
        ]
        for target_index, target_column in enumerate(target_columns):
            adjacency = adjacency_for_target(
                table,
                source_column=source_column,
                target_column=target_column,
                quantile=args.quantile,
                radius_multiplier=args.radius_multiplier,
                permutations=args.permutations,
                seed=args.seed + sample_index * 10 + target_index,
            )
            adjacency_rows.append(adjacency)
            print(
                f"{sample.sample_accession} {meta['stage']} {target_column}: "
                f"delta={adjacency['enrichment_delta']:.3f}; "
                f"p_greater={adjacency['empirical_p_greater']:.4f}; "
                f"p_less={adjacency['empirical_p_less']:.4f}"
            )

    args.table_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    scores = pd.concat(score_tables, ignore_index=True)
    adjacency = pd.DataFrame(adjacency_rows)
    stage_summary = (
        adjacency.groupby(["target", "stage"])
        .agg(
            n_samples=("sample", "nunique"),
            observed_fraction_mean=("observed_fraction", "mean"),
            null_mean_mean=("null_mean", "mean"),
            enrichment_delta_mean=("enrichment_delta", "mean"),
            empirical_p_greater_median=("empirical_p_greater", "median"),
            empirical_p_less_median=("empirical_p_less", "median"),
        )
        .reset_index()
    )

    scores.to_csv(args.table_dir / "gse189487_refined_signature_spatial_scores.csv", index=False, encoding="utf-8-sig")
    adjacency.to_csv(args.table_dir / "gse189487_refined_signature_adjacency.csv", index=False, encoding="utf-8-sig")
    stage_summary.to_csv(
        args.table_dir / "gse189487_refined_signature_adjacency_by_stage.csv",
        index=False,
        encoding="utf-8-sig",
    )
    (args.table_dir / "gse189487_refined_signature_genes_used.json").write_text(
        json.dumps(genes_used, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    plot_stage_target_summary(
        stage_summary,
        args.figure_dir / "gse189487_refined_signature_adjacency_by_stage.png",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
