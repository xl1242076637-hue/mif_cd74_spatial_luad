"""Spatial adjacency and niche enrichment utilities."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors


def top_quantile_mask(values: pd.Series, quantile: float = 0.75) -> pd.Series:
    """Return a boolean mask for values strictly above the requested quantile."""
    threshold = values.quantile(quantile)
    return values > threshold


def finite_coordinate_mask(coords: np.ndarray) -> np.ndarray:
    """Return rows with finite x/y coordinates."""
    coords = np.asarray(coords, dtype=float)
    return np.isfinite(coords).all(axis=1)


def median_nearest_neighbor_distance(coords: np.ndarray) -> float:
    """Return the median nearest-neighbor distance across all coordinates."""
    coords = np.asarray(coords, dtype=float)
    if len(coords) < 2:
        raise ValueError("At least two coordinates are required.")
    nn = NearestNeighbors(n_neighbors=2)
    nn.fit(coords)
    distances, _ = nn.kneighbors(coords)
    return float(np.median(distances[:, 1]))


def nearest_target_fraction(
    coords: np.ndarray,
    source_mask: np.ndarray,
    target_mask: np.ndarray,
    radius: float,
) -> float:
    """Fraction of source points with at least one target point within radius."""
    coords = np.asarray(coords, dtype=float)
    source_mask = np.asarray(source_mask, dtype=bool)
    target_mask = np.asarray(target_mask, dtype=bool)
    source_coords = coords[source_mask]
    target_coords = coords[target_mask]
    if len(source_coords) == 0 or len(target_coords) == 0:
        return float("nan")
    nn = NearestNeighbors(n_neighbors=1)
    nn.fit(target_coords)
    distances, _ = nn.kneighbors(source_coords)
    return float(np.mean(distances[:, 0] <= radius))


def permutation_adjacency_test(
    coords: np.ndarray,
    source_mask: np.ndarray,
    target_mask: np.ndarray,
    radius: float,
    n_permutations: int = 500,
    seed: int = 13,
) -> dict:
    """Compare observed source-target adjacency to random target-label permutations."""
    observed = nearest_target_fraction(coords, source_mask, target_mask, radius)
    rng = np.random.default_rng(seed)
    target_count = int(np.asarray(target_mask, dtype=bool).sum())
    null = []
    for _ in range(n_permutations):
        shuffled = np.zeros(len(target_mask), dtype=bool)
        shuffled[rng.choice(len(target_mask), size=target_count, replace=False)] = True
        null.append(nearest_target_fraction(coords, source_mask, shuffled, radius))
    null_array = np.asarray(null, dtype=float)
    null_mean = float(np.nanmean(null_array))
    empirical_p = (float(np.sum(null_array >= observed)) + 1.0) / (len(null_array) + 1.0)
    empirical_p_less = (float(np.sum(null_array <= observed)) + 1.0) / (len(null_array) + 1.0)
    return {
        "observed_fraction": observed,
        "null_mean": null_mean,
        "null_sd": float(np.nanstd(null_array, ddof=1)),
        "enrichment_delta": observed - null_mean,
        "empirical_p_greater": empirical_p,
        "empirical_p_less": empirical_p_less,
        "n_permutations": n_permutations,
        "radius": radius,
    }


def summarize_adjacency_by_stage(
    adjacency: pd.DataFrame,
    stage_order: tuple[str, ...] = ("AIS", "MIA", "IAC"),
) -> pd.DataFrame:
    """Summarize per-sample epithelial-macrophage adjacency statistics by stage."""
    df = adjacency.copy()
    df["enrichment_delta"] = df["observed_fraction"] - df["null_mean"]
    summary = (
        df.groupby("stage", dropna=False)
        .agg(
            n_samples=("sample", "nunique"),
            observed_fraction_mean=("observed_fraction", "mean"),
            observed_fraction_sd=("observed_fraction", "std"),
            null_mean_mean=("null_mean", "mean"),
            enrichment_delta_mean=("enrichment_delta", "mean"),
            enrichment_delta_sd=("enrichment_delta", "std"),
            empirical_p_greater_median=("empirical_p_greater", "median"),
            empirical_p_less_median=("empirical_p_less", "median"),
        )
        .reset_index()
    )
    order = {stage: index for index, stage in enumerate(stage_order)}
    summary["_stage_order"] = summary["stage"].map(order).fillna(len(order))
    summary = summary.sort_values(["_stage_order", "stage"]).drop(columns="_stage_order")
    numeric_columns = summary.select_dtypes(include="number").columns
    summary[numeric_columns] = summary[numeric_columns].round(12)
    return summary
