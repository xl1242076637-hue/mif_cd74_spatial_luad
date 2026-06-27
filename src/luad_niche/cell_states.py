"""Helpers for lightweight marker-score cell-state assignment."""

from __future__ import annotations

import numpy as np
import pandas as pd


def assign_dominant_labels(
    scores: pd.DataFrame,
    label_to_column: dict[str, str],
    min_score: float = 0.05,
    min_margin: float = 0.05,
    unassigned_label: str = "unassigned",
    ambiguous_label: str = "ambiguous",
) -> pd.Series:
    """Assign each row to the label with the highest score."""
    missing = [column for column in label_to_column.values() if column not in scores.columns]
    if missing:
        raise ValueError(f"Missing score columns: {missing}")

    columns = list(label_to_column.values())
    labels = list(label_to_column.keys())
    values = scores[columns].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy()
    assigned = []
    for row in values:
        order = np.argsort(row)[::-1]
        best_index = int(order[0])
        best_score = float(row[best_index])
        second_score = float(row[int(order[1])]) if len(order) > 1 else float("-inf")
        if best_score < min_score:
            assigned.append(unassigned_label)
        elif best_score - second_score < min_margin:
            assigned.append(ambiguous_label)
        else:
            assigned.append(labels[best_index])
    return pd.Series(assigned, index=scores.index, name="dominant_label")


def add_within_group_high_flags(
    df: pd.DataFrame,
    group_column: str,
    score_column: str,
    group_value: str,
    quantile: float,
    output_column: str,
) -> pd.DataFrame:
    """Mark rows above a score quantile computed inside one group."""
    result = df.copy()
    result[output_column] = False
    group_mask = result[group_column] == group_value
    if not group_mask.any():
        return result
    threshold = result.loc[group_mask, score_column].quantile(quantile)
    result.loc[group_mask, output_column] = result.loc[group_mask, score_column] > threshold
    return result
