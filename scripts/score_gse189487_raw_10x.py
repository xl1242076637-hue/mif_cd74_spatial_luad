#!/usr/bin/env python
"""Score and test epithelial-macrophage spatial niches in GSE189487 raw 10x data."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import yaml


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
    summarize_adjacency_by_stage,
    top_quantile_mask,
)
from luad_niche.tenx import (  # noqa: E402
    TenXSampleFiles,
    discover_10x_samples,
    read_10x_obs,
    read_10x_selected_genes_with_totals,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=PROJECT_ROOT / "data" / "interim" / "GSE189487" / "raw_10x",
        help="Directory containing extracted GSE189487 10x files.",
    )
    parser.add_argument(
        "--genes",
        type=Path,
        default=PROJECT_ROOT / "config" / "candidate_genes.yaml",
        help="Candidate marker gene YAML.",
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
    parser.add_argument(
        "--scale-factor",
        type=float,
        default=10_000.0,
        help="Per-spot total-count scaling factor before log1p.",
    )
    parser.add_argument(
        "--quantile",
        type=float,
        default=0.75,
        help="High-score quantile threshold.",
    )
    parser.add_argument(
        "--radius-multiplier",
        type=float,
        default=1.0,
        help="Adjacency radius as a multiplier of median nearest-neighbor distance.",
    )
    parser.add_argument(
        "--permutations",
        type=int,
        default=500,
        help="Number of target-label permutations per sample.",
    )
    parser.add_argument("--seed", type=int, default=13, help="Base random seed.")
    parser.add_argument(
        "--include-off-tissue",
        action="store_true",
        help="Keep off-tissue Visium positions in score tables and plots.",
    )
    return parser.parse_args()


def load_panels(path: Path) -> dict[str, list[str]]:
    with Path(path).open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    return {name: list(genes) for name, genes in config["panels"].items()}


def flatten_panel_genes(panels: dict[str, list[str]]) -> list[str]:
    genes: list[str] = []
    for panel_genes in panels.values():
        genes.extend(gene for gene in panel_genes if gene not in genes)
    return genes


def sample_metadata(metadata: pd.DataFrame, sample_accession: str) -> dict[str, object]:
    rows = metadata[metadata["sample_accession"] == sample_accession]
    if rows.empty:
        return {
            "sample_accession": sample_accession,
            "stage": "",
            "histological_type": "",
            "radiological_type": "",
            "gender": "",
            "title": "",
        }
    row = rows.iloc[0]
    return {
        "sample_accession": sample_accession,
        "stage": row.get("interpreted_stage", ""),
        "histological_type": row.get("histological_type", ""),
        "radiological_type": row.get("radiological_type", ""),
        "gender": row.get("gender", ""),
        "title": row.get("title", ""),
    }


def score_sample(
    sample: TenXSampleFiles,
    panels: dict[str, list[str]],
    metadata_row: dict[str, object],
    scale_factor: float,
    include_off_tissue: bool,
) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    all_genes = flatten_panel_genes(panels)
    obs = read_10x_obs(sample)
    counts, total_counts = read_10x_selected_genes_with_totals(sample, all_genes)
    normalized = normalize_log1p_counts_by_totals(counts, total_counts, scale_factor=scale_factor)
    scores = compute_panel_scores(normalized, panels)

    table = obs.merge(total_counts.rename_axis("spot").reset_index(), on="spot", how="left")
    table = table.merge(scores.reset_index(names="spot"), on="spot", how="left")
    if not include_off_tissue and "in_tissue" in table:
        table = table[table["in_tissue"] == 1].copy()

    table.insert(0, "sample_name", sample.sample_name)
    table.insert(0, "sample", sample.sample_accession)
    for key in ["stage", "histological_type", "radiological_type", "gender", "title"]:
        table[key] = metadata_row.get(key, "")
    return table.reset_index(drop=True), scores.attrs["panel_genes_used"]


def plot_score_panels(table: pd.DataFrame, output: Path) -> None:
    score_columns = [column for column in table.columns if column.endswith("_score")]
    ncols = 2
    nrows = (len(score_columns) + 1) // ncols
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 4.5 * nrows), squeeze=False)
    for ax, score_column in zip(axes.ravel(), score_columns, strict=False):
        plot_df = table[table[score_column].notna()]
        sc = ax.scatter(
            plot_df["x"],
            -plot_df["y"],
            c=plot_df[score_column],
            s=8,
            cmap="viridis",
            linewidths=0,
        )
        ax.set_title(score_column)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    for ax in axes.ravel()[len(score_columns) :]:
        ax.axis("off")
    sample = table["sample"].iloc[0]
    stage = table["stage"].iloc[0]
    fig.suptitle(f"{sample} {stage}: raw 10x candidate panel scores", y=0.995)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def adjacency_for_sample(
    table: pd.DataFrame,
    quantile: float,
    radius_multiplier: float,
    permutations: int,
    seed: int,
) -> dict[str, object]:
    usable = table[
        table[["x", "y", "epithelial_progenitor_score", "macrophage_niche_score"]]
        .notna()
        .all(axis=1)
    ].copy()
    coords_all = usable[["x", "y"]].to_numpy(dtype=float)
    finite_mask = finite_coordinate_mask(coords_all)
    dropped_nonfinite = int((~finite_mask).sum())
    usable = usable.loc[finite_mask].reset_index(drop=True)
    coords = usable[["x", "y"]].to_numpy(dtype=float)
    radius = radius_multiplier * median_nearest_neighbor_distance(coords)
    epithelial_high = top_quantile_mask(usable["epithelial_progenitor_score"], quantile).to_numpy()
    macrophage_high = top_quantile_mask(usable["macrophage_niche_score"], quantile).to_numpy()
    overlap = epithelial_high & macrophage_high

    stats = permutation_adjacency_test(
        coords,
        epithelial_high,
        macrophage_high,
        radius=radius,
        n_permutations=permutations,
        seed=seed,
    )
    stats.update(
        {
            "sample": usable["sample"].iloc[0],
            "sample_name": usable["sample_name"].iloc[0],
            "stage": usable["stage"].iloc[0],
            "histological_type": usable["histological_type"].iloc[0],
            "radiological_type": usable["radiological_type"].iloc[0],
            "source": "epithelial_progenitor_high",
            "target": "macrophage_niche_high",
            "quantile": quantile,
            "radius_multiplier": radius_multiplier,
            "n_spots": len(usable),
            "n_dropped_nonfinite_coordinates": dropped_nonfinite,
            "n_source_high": int(epithelial_high.sum()),
            "n_target_high": int(macrophage_high.sum()),
            "n_overlap_high": int(overlap.sum()),
            "overlap_fraction_of_source": float(overlap.sum() / epithelial_high.sum()),
            "target_to_source_fraction": nearest_target_fraction(
                coords, macrophage_high, epithelial_high, radius=radius
            ),
        }
    )
    return stats


def plot_adjacency_overlay(table: pd.DataFrame, output: Path, quantile: float) -> None:
    plot_df = table[
        table[["x", "y", "epithelial_progenitor_score", "macrophage_niche_score"]]
        .notna()
        .all(axis=1)
    ].copy()
    finite_mask = finite_coordinate_mask(plot_df[["x", "y"]].to_numpy(dtype=float))
    plot_df = plot_df.loc[finite_mask].reset_index(drop=True)
    epithelial_high = top_quantile_mask(plot_df["epithelial_progenitor_score"], quantile).to_numpy()
    macrophage_high = top_quantile_mask(plot_df["macrophage_niche_score"], quantile).to_numpy()
    overlap = epithelial_high & macrophage_high

    plot_df["class"] = "other"
    plot_df.loc[epithelial_high, "class"] = "epithelial_high"
    plot_df.loc[macrophage_high, "class"] = "macrophage_high"
    plot_df.loc[overlap, "class"] = "both_high"
    colors = {
        "other": "#d0d0d0",
        "epithelial_high": "#d95f02",
        "macrophage_high": "#1b9e77",
        "both_high": "#7570b3",
    }
    fig, ax = plt.subplots(figsize=(7, 7))
    for label, color in colors.items():
        subset = plot_df[plot_df["class"] == label]
        ax.scatter(subset["x"], -subset["y"], s=9, c=color, label=label, linewidths=0)
    sample = plot_df["sample"].iloc[0]
    stage = plot_df["stage"].iloc[0]
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"{sample} {stage}: epithelial-high vs macrophage-high")
    ax.legend(loc="upper right", frameon=False, markerscale=2)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_stage_summary(summary: pd.DataFrame, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    x = range(len(summary))
    width = 0.36
    ax.bar(
        [value - width / 2 for value in x],
        summary["observed_fraction_mean"],
        width=width,
        label="observed",
        color="#d95f02",
    )
    ax.bar(
        [value + width / 2 for value in x],
        summary["null_mean_mean"],
        width=width,
        label="permuted null",
        color="#7570b3",
    )
    ax.set_xticks(list(x))
    ax.set_xticklabels(summary["stage"])
    ax.set_ylabel("Fraction of epithelial-high spots near macrophage-high spots")
    ax.set_title("GSE189487 epithelial-macrophage adjacency by stage")
    ax.legend(frameon=False)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    panels = load_panels(args.genes)
    metadata = pd.read_csv(args.metadata)
    samples = discover_10x_samples(args.input_dir)

    all_scores = []
    all_adjacency = []
    genes_used_by_sample = {}
    for index, sample in enumerate(samples.values()):
        meta = sample_metadata(metadata, sample.sample_accession)
        score_table, genes_used = score_sample(
            sample,
            panels,
            meta,
            scale_factor=args.scale_factor,
            include_off_tissue=args.include_off_tissue,
        )
        all_scores.append(score_table)
        genes_used_by_sample[sample.sample_accession] = genes_used

        per_sample_scores = args.table_dir / f"gse189487_{sample.sample_accession.lower()}_raw_spatial_panel_scores.csv"
        score_table.to_csv(per_sample_scores, index=False, encoding="utf-8-sig")
        plot_score_panels(
            score_table,
            args.figure_dir / f"gse189487_{sample.sample_accession.lower()}_raw_spatial_panel_scores.png",
        )

        adjacency = adjacency_for_sample(
            score_table,
            quantile=args.quantile,
            radius_multiplier=args.radius_multiplier,
            permutations=args.permutations,
            seed=args.seed + index,
        )
        all_adjacency.append(adjacency)
        plot_adjacency_overlay(
            score_table,
            args.figure_dir / f"gse189487_{sample.sample_accession.lower()}_raw_epithelial_macrophage_adjacency.png",
            quantile=args.quantile,
        )
        print(
            f"{sample.sample_accession}: stage={meta['stage']}; "
            f"spots={len(score_table)}; observed={adjacency['observed_fraction']:.3f}; "
            f"null={adjacency['null_mean']:.3f}; delta={adjacency['enrichment_delta']:.3f}; "
            f"p_greater={adjacency['empirical_p_greater']:.4f}; "
            f"p_less={adjacency['empirical_p_less']:.4f}"
        )

    args.table_dir.mkdir(parents=True, exist_ok=True)
    combined_scores = pd.concat(all_scores, ignore_index=True)
    combined_scores_output = args.table_dir / "gse189487_raw_spatial_panel_scores.csv"
    combined_scores.to_csv(combined_scores_output, index=False, encoding="utf-8-sig")

    genes_used_output = args.table_dir / "gse189487_raw_panel_genes_used.json"
    genes_used_output.write_text(
        json.dumps(genes_used_by_sample, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    adjacency_df = pd.DataFrame(all_adjacency)
    adjacency_output = args.table_dir / "gse189487_raw_epithelial_macrophage_adjacency.csv"
    adjacency_df.to_csv(adjacency_output, index=False, encoding="utf-8-sig")

    stage_summary = summarize_adjacency_by_stage(adjacency_df)
    stage_summary_output = args.table_dir / "gse189487_raw_epithelial_macrophage_adjacency_by_stage.csv"
    stage_summary.to_csv(stage_summary_output, index=False, encoding="utf-8-sig")
    plot_stage_summary(
        stage_summary,
        args.figure_dir / "gse189487_raw_epithelial_macrophage_adjacency_by_stage.png",
    )

    print(f"Wrote combined scores: {combined_scores_output}")
    print(f"Wrote genes used: {genes_used_output}")
    print(f"Wrote adjacency: {adjacency_output}")
    print(f"Wrote stage summary: {stage_summary_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
