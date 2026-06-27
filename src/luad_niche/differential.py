"""Sparse aggregate helpers for lightweight marker ranking."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.sparse import spmatrix


def aggregate_contrast(
    matrix: spmatrix,
    group_mask: np.ndarray,
    background_mask: np.ndarray,
) -> pd.DataFrame:
    """Aggregate feature-by-cell sparse counts for a group/background contrast."""
    group_mask = np.asarray(group_mask, dtype=bool)
    background_mask = np.asarray(background_mask, dtype=bool)
    group_matrix = matrix[:, group_mask]
    background_matrix = matrix[:, background_mask]
    aggregate = pd.DataFrame(
        {
            "group_sum": np.asarray(group_matrix.sum(axis=1)).ravel(),
            "background_sum": np.asarray(background_matrix.sum(axis=1)).ravel(),
            "group_detected": np.asarray((group_matrix > 0).sum(axis=1)).ravel(),
            "background_detected": np.asarray((background_matrix > 0).sum(axis=1)).ravel(),
        }
    )
    aggregate.attrs["n_group"] = int(group_mask.sum())
    aggregate.attrs["n_background"] = int(background_mask.sum())
    return aggregate


def rank_markers(
    aggregate: pd.DataFrame,
    pseudocount: float = 1e-6,
) -> pd.DataFrame:
    """Rank marker genes from aggregate sums and detection counts."""
    ranked = aggregate.copy()
    ranked["mean_group"] = ranked["group_sum"] / ranked["n_group"]
    ranked["mean_background"] = ranked["background_sum"] / ranked["n_background"]
    ranked["pct_group"] = ranked["group_detected"] / ranked["n_group"]
    ranked["pct_background"] = ranked["background_detected"] / ranked["n_background"]
    ranked["log2fc"] = np.log2(
        (ranked["mean_group"] + pseudocount) / (ranked["mean_background"] + pseudocount)
    )
    ranked["delta_pct"] = ranked["pct_group"] - ranked["pct_background"]
    ranked["marker_score"] = ranked["log2fc"] * ranked["delta_pct"]
    return ranked.sort_values(
        ["marker_score", "log2fc", "pct_group"],
        ascending=[False, False, False],
    ).reset_index(drop=True)
