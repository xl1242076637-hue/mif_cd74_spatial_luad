#!/usr/bin/env python
"""Test epithelial-high and macrophage-high spatial adjacency in one scored sample."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from luad_niche.spatial_niche import (  # noqa: E402
    finite_coordinate_mask,
    median_nearest_neighbor_distance,
    nearest_target_fraction,
    permutation_adjacency_test,
    top_quantile_mask,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scores",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gsm5702474_spatial_panel_scores.csv",
        help="Per-spot panel score table.",
    )
    parser.add_argument(
        "--sample",
        default="GSM5702474",
        help="Sample identifier.",
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
        default=2.0,
        help="Radius as a multiplier of median nearest-neighbor distance.",
    )
    parser.add_argument(
        "--permutations",
        type=int,
        default=500,
        help="Number of target-label permutations.",
    )
    parser.add_argument(
        "--table-output",
        type=Path,
        default=PROJECT_ROOT / "results" / "tables" / "gsm5702474_epithelial_macrophage_adjacency.csv",
        help="Output summary CSV.",
    )
    parser.add_argument(
        "--figure-output",
        type=Path,
        default=PROJECT_ROOT / "results" / "figures" / "gsm5702474_epithelial_macrophage_adjacency.png",
        help="Output overlay figure.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    df = pd.read_csv(args.scores)
    coords_all = df[["x", "y"]].to_numpy(dtype=float)
    finite_mask = finite_coordinate_mask(coords_all)
    dropped_nonfinite = int((~finite_mask).sum())
    df = df.loc[finite_mask].reset_index(drop=True)
    coords = df[["x", "y"]].to_numpy(dtype=float)
    radius = args.radius_multiplier * median_nearest_neighbor_distance(coords)
    epithelial_high = top_quantile_mask(df["epithelial_progenitor_score"], args.quantile).to_numpy()
    macrophage_high = top_quantile_mask(df["macrophage_niche_score"], args.quantile).to_numpy()
    overlap = epithelial_high & macrophage_high

    stats = permutation_adjacency_test(
        coords,
        epithelial_high,
        macrophage_high,
        radius=radius,
        n_permutations=args.permutations,
    )
    stats.update(
        {
            "sample": args.sample,
            "source": "epithelial_progenitor_high",
            "target": "macrophage_niche_high",
            "quantile": args.quantile,
            "n_spots": len(df),
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

    args.table_output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([stats]).to_csv(args.table_output, index=False, encoding="utf-8-sig")

    plot_df = df.copy()
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
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"{args.sample}: epithelial-high vs macrophage-high spots")
    ax.legend(loc="upper right", frameon=False, markerscale=2)
    fig.tight_layout()
    args.figure_output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.figure_output, dpi=220)
    plt.close(fig)

    print(f"Wrote table: {args.table_output}")
    print(f"Wrote figure: {args.figure_output}")
    print(pd.DataFrame([stats]).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
