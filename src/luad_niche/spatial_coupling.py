"""Continuous spatial coupling scores for source-target expression programs."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

from luad_niche.spatial_niche import finite_coordinate_mask


def continuous_spatial_coupling(
    table: pd.DataFrame,
    source_column: str,
    target_column: str,
    radius: float,
) -> dict[str, float | int]:
    """Compute a continuous source-to-neighbor-target coupling score."""
    usable = table[table[["x", "y", source_column, target_column]].notna().all(axis=1)].copy()
    finite = finite_coordinate_mask(usable[["x", "y"]].to_numpy(dtype=float))
    usable = usable.loc[finite].reset_index(drop=True)
    if len(usable) < 2 or radius <= 0:
        return {
            "n_spots": int(len(usable)),
            "n_spots_with_neighbors": 0,
            "source_mean": float("nan"),
            "target_mean": float("nan"),
            "neighbor_target_mean": float("nan"),
            "coupling_score": float("nan"),
            "source_weighted_neighbor_target_mean": float("nan"),
        }

    coords = usable[["x", "y"]].to_numpy(dtype=float)
    source = usable[source_column].to_numpy(dtype=float)
    target = usable[target_column].to_numpy(dtype=float)
    neighbors = NearestNeighbors(radius=radius)
    neighbors.fit(coords)
    neighbor_indices = neighbors.radius_neighbors(coords, return_distance=False)

    neighbor_target_means = np.full(len(usable), np.nan, dtype=float)
    for index, indices in enumerate(neighbor_indices):
        other_indices = [neighbor_index for neighbor_index in indices if neighbor_index != index]
        if other_indices:
            neighbor_target_means[index] = float(np.mean(target[other_indices]))

    has_neighbors = np.isfinite(neighbor_target_means)
    coupling_values = source[has_neighbors] * neighbor_target_means[has_neighbors]
    source_with_neighbors = source[has_neighbors]
    source_sum = float(np.sum(source_with_neighbors))
    weighted_neighbor = (
        float(np.sum(coupling_values) / source_sum)
        if source_sum > 0 and len(coupling_values)
        else float("nan")
    )
    return {
        "n_spots": int(len(usable)),
        "n_spots_with_neighbors": int(has_neighbors.sum()),
        "source_mean": float(np.mean(source)),
        "target_mean": float(np.mean(target)),
        "neighbor_target_mean": float(np.nanmean(neighbor_target_means)) if has_neighbors.any() else float("nan"),
        "coupling_score": float(np.mean(coupling_values)) if len(coupling_values) else float("nan"),
        "source_weighted_neighbor_target_mean": weighted_neighbor,
    }


def summarize_continuous_effects(
    effects: pd.DataFrame,
    stages: tuple[str, ...] = ("MIA", "LUAD"),
) -> pd.DataFrame:
    """Summarize continuous perturbation effects for selected late/progressive stages."""
    subset = effects[effects["stage"].isin(stages)].copy()
    if subset.empty:
        return pd.DataFrame()
    group_columns = ["perturbation_id", "axis_id", "evidence_type"]
    summary = (
        subset.groupby(group_columns, dropna=False)
        .agg(
            n_samples=("sample", "nunique"),
            baseline_coupling_score=("baseline_coupling_score", "mean"),
            perturbed_coupling_score=("perturbed_coupling_score", "mean"),
            coupling_delta_mean=("coupling_delta", "mean"),
            coupling_relative_delta_mean=("coupling_relative_delta", "mean"),
        )
        .reset_index()
    )
    summary["continuous_priority_score"] = summary["coupling_relative_delta_mean"].apply(
        lambda value: -float(value) if not math.isnan(float(value)) and float(value) < 0 else 0.0
    )
    numeric_columns = summary.select_dtypes(include="number").columns
    summary[numeric_columns] = summary[numeric_columns].round(12)
    return summary.sort_values(
        ["continuous_priority_score", "coupling_relative_delta_mean"],
        ascending=[False, True],
    ).reset_index(drop=True)
