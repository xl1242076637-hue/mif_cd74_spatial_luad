#!/usr/bin/env python
"""Parameter robustness for the GSE189487 epithelial progenitor-SPP1 macrophage niche."""

from __future__ import annotations

import argparse
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
from luad_niche.signatures import build_top_n_signature_panels  # noqa: E402
from luad_niche.spatial_niche import (  # noqa: E402
    finite_coordinate_mask,
    median_nearest_neighbor_distance,
    nearest_target_fraction,
    permutation_adjacency_test,
    top_quantile_mask,
)
from luad_niche.tenx import discover_10x_samples, read_10x_obs, read_10x_selected_genes_with_totals  # noqa: E402


SOURCE_CONTRAST = "epithelial_progenitor_like_vs_other_epithelial"
TARGET_CONTRAST = "spp1_macrophage_vs_other_macrophage"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=PROJECT_ROOT / "data" / "interim" / "GSE189487" / "raw_10x",
    )
    parser.add_argument(
        "--ranked-markers",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gse189357_refined_state_markers.csv",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "geo_sample_metadata_annotated.csv",
    )
    parser.add_argument("--top-ns", default="10,20,30")
    parser.add_argument("--quantiles", default="0.70,0.75,0.80")
    parser.add_argument("--radius-multipliers", default="0.75,1.00,1.25")
    parser.add_argument("--permutations", type=int, default=200)
    parser.add_argument("--seed", type=int, default=101)
    parser.add_argument("--table-dir", type=Path, default=PROJECT_ROOT / "results" / "tables")
    parser.add_argument("--figure-dir", type=Path, default=PROJECT_ROOT / "results" / "figures")
    return parser.parse_args()


def parse_number_list(value: str, cast):
    return [cast(item.strip()) for item in value.split(",") if item.strip()]


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


def adjacency_stats(
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
            "n_spots": len(usable),
            "n_source_high": int(source_high.sum()),
            "n_target_high": int(target_high.sum()),
            "n_overlap_high": int(overlap.sum()),
            "overlap_fraction_of_source": float(overlap.sum() / source_high.sum()),
            "target_to_source_fraction": nearest_target_fraction(coords, target_high, source_high, radius),
        }
    )
    return stats


def plot_mia_robustness(stage_summary: pd.DataFrame, output: Path) -> None:
    mia = stage_summary[stage_summary["stage"].eq("MIA")]
    top_ns = sorted(mia["top_n"].unique())
    fig, axes = plt.subplots(
        1,
        len(top_ns),
        figsize=(4.6 * len(top_ns), 4.0),
        squeeze=False,
        constrained_layout=True,
    )
    image = None
    for panel_index, (ax, top_n) in enumerate(zip(axes.ravel(), top_ns, strict=False)):
        subset = mia[mia["top_n"].eq(top_n)]
        pivot = subset.pivot_table(
            index="radius_multiplier",
            columns="quantile",
            values="enrichment_delta_mean",
        ).sort_index(ascending=False)
        image = ax.imshow(pivot.to_numpy(), cmap="coolwarm", vmin=-0.15, vmax=0.15, aspect="auto")
        ax.set_title(f"top{top_n}")
        ax.set_xticks(range(len(pivot.columns)))
        ax.set_xticklabels([f"{value:.2f}" for value in pivot.columns])
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels([f"{value:.2f}" for value in pivot.index])
        ax.set_xlabel("quantile")
        ax.set_ylabel("radius" if panel_index == 0 else "")
        for y, radius in enumerate(pivot.index):
            for x, quantile in enumerate(pivot.columns):
                value = pivot.loc[radius, quantile]
                ax.text(x, y, f"{value:.2f}", ha="center", va="center", fontsize=8)
    if image is not None:
        fig.colorbar(image, ax=axes.ravel().tolist(), shrink=0.88, pad=0.02, label="MIA mean delta")
    fig.suptitle("MIA epithelial progenitor-SPP1 macrophage niche robustness")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    top_ns = parse_number_list(args.top_ns, int)
    quantiles = parse_number_list(args.quantiles, float)
    radius_multipliers = parse_number_list(args.radius_multipliers, float)
    ranked = pd.read_csv(args.ranked_markers)
    panels = build_top_n_signature_panels(
        ranked,
        contrasts=[SOURCE_CONTRAST, TARGET_CONTRAST],
        top_ns=top_ns,
        min_pct_group=0.05,
    )
    all_genes = flatten_genes(panels)
    metadata = pd.read_csv(args.metadata)
    samples = discover_10x_samples(args.input_dir, require_tissue_positions=True)

    rows = []
    for sample_index, sample in enumerate(samples.values()):
        meta = sample_metadata(metadata, sample.sample_accession)
        obs = read_10x_obs(sample)
        counts, totals = read_10x_selected_genes_with_totals(sample, all_genes)
        normalized = normalize_log1p_counts_by_totals(counts, totals)
        scores = compute_panel_scores(normalized, panels)
        table = obs.merge(totals.rename_axis("spot").reset_index(), on="spot", how="left")
        table = table.merge(scores.reset_index(names="spot"), on="spot", how="left")
        table = table[table["in_tissue"] == 1].copy()
        table["sample"] = sample.sample_accession
        table["sample_name"] = sample.sample_name
        table["stage"] = meta["stage"]

        for top_n in top_ns:
            source_column = f"{SOURCE_CONTRAST}_top{top_n}_score"
            target_column = f"{TARGET_CONTRAST}_top{top_n}_score"
            for quantile in quantiles:
                for radius_multiplier in radius_multipliers:
                    stats = adjacency_stats(
                        table,
                        source_column,
                        target_column,
                        quantile=quantile,
                        radius_multiplier=radius_multiplier,
                        permutations=args.permutations,
                        seed=args.seed + sample_index * 1000 + top_n * 10 + int(quantile * 100) + int(radius_multiplier * 100),
                    )
                    stats.update(
                        {
                            "sample": sample.sample_accession,
                            "sample_name": sample.sample_name,
                            "stage": meta["stage"],
                            "top_n": top_n,
                            "quantile": quantile,
                            "radius_multiplier": radius_multiplier,
                        }
                    )
                    rows.append(stats)
        print(f"Robustness complete for {sample.sample_accession} ({meta['stage']})")

    results = pd.DataFrame(rows)
    stage_summary = (
        results.groupby(["top_n", "quantile", "radius_multiplier", "stage"])
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
    mia_summary = stage_summary[stage_summary["stage"].eq("MIA")].copy()
    mia_summary["positive_and_p_greater_lt_0_05"] = (
        (mia_summary["enrichment_delta_mean"] > 0)
        & (mia_summary["empirical_p_greater_median"] < 0.05)
    )

    args.table_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.table_dir / "gse189487_spp1_niche_robustness.csv", index=False, encoding="utf-8-sig")
    stage_summary.to_csv(
        args.table_dir / "gse189487_spp1_niche_robustness_by_stage.csv",
        index=False,
        encoding="utf-8-sig",
    )
    mia_summary.to_csv(
        args.table_dir / "gse189487_spp1_niche_robustness_mia_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )
    plot_mia_robustness(
        stage_summary,
        args.figure_dir / "gse189487_spp1_niche_robustness_mia.png",
    )
    print(
        "MIA positive parameter sets: "
        f"{int(mia_summary['positive_and_p_greater_lt_0_05'].sum())}/{len(mia_summary)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
