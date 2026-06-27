"""Patient-aware statistical summaries for spatial candidate-axis results."""

from __future__ import annotations

import re
from collections.abc import Iterable

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon


PHASE_BY_STAGE = {
    "AAH": "precursor",
    "AIS": "precursor",
    "MIA": "late",
    "LUAD": "late",
}


def extract_patient_id(sample_name: str) -> str | None:
    """Extract a stable P<number> patient identifier from a spatial sample name."""
    match = re.search(r"(?:^|[^A-Za-z0-9])P(\d+)(?!\d)", str(sample_name), flags=re.IGNORECASE)
    if not match:
        return None
    return f"P{int(match.group(1))}"


def add_patient_phase(adjacency: pd.DataFrame) -> pd.DataFrame:
    """Annotate spatial rows with patient IDs and a precursor-versus-late phase."""
    result = adjacency.copy()
    if "status" not in result.columns:
        result["status"] = "ok"
    result["patient_id"] = result["sample_name"].map(extract_patient_id)
    result["phase"] = result["stage"].map(PHASE_BY_STAGE)
    return result


def benjamini_hochberg(p_values: Iterable[float]) -> np.ndarray:
    """Return Benjamini-Hochberg adjusted p-values while preserving NaNs."""
    values = np.asarray(list(p_values), dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    valid_indexes = np.flatnonzero(np.isfinite(values))
    if not len(valid_indexes):
        return adjusted
    valid_values = values[valid_indexes]
    order = np.argsort(valid_values)
    ranked = valid_values[order]
    corrected = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    corrected = np.minimum.accumulate(corrected[::-1])[::-1]
    corrected = np.clip(corrected, 0.0, 1.0)
    adjusted[valid_indexes[order]] = corrected
    return adjusted


def _bootstrap_mean_difference(
    late_values: np.ndarray,
    precursor_values: np.ndarray,
    *,
    iterations: int,
    seed: int,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    differences = np.empty(iterations, dtype=float)
    for index in range(iterations):
        late_sample = rng.choice(late_values, size=len(late_values), replace=True)
        precursor_sample = rng.choice(precursor_values, size=len(precursor_values), replace=True)
        differences[index] = late_sample.mean() - precursor_sample.mean()
    return tuple(np.quantile(differences, [0.025, 0.975]))


def _bootstrap_mean(values: np.ndarray, *, iterations: int, seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    means = np.empty(iterations, dtype=float)
    for index in range(iterations):
        means[index] = rng.choice(values, size=len(values), replace=True).mean()
    return tuple(np.quantile(means, [0.025, 0.975]))


def _valid_phase_rows(adjacency: pd.DataFrame) -> pd.DataFrame:
    annotated = add_patient_phase(adjacency)
    valid = annotated[
        annotated["status"].eq("ok")
        & annotated["phase"].isin(["precursor", "late"])
        & annotated["enrichment_delta"].notna()
        & annotated["patient_id"].notna()
    ].copy()
    return valid


def summarize_late_vs_precursor(
    adjacency: pd.DataFrame,
    *,
    aggregate_by_patient: bool,
    bootstrap_iterations: int = 2_000,
    seed: int = 307534,
) -> pd.DataFrame:
    """Compare MIA/LUAD versus AAH/AIS enrichment at sample or patient level."""
    valid = _valid_phase_rows(adjacency)
    group_columns = ["axis_id", "axis_label", "evidence_type"]
    if aggregate_by_patient:
        valid = (
            valid.groupby(group_columns + ["patient_id", "phase"], as_index=False)
            .agg(enrichment_delta=("enrichment_delta", "mean"))
        )
    records = []
    for row_index, (keys, group) in enumerate(valid.groupby(group_columns, sort=True)):
        late = group.loc[group["phase"].eq("late"), "enrichment_delta"].to_numpy(dtype=float)
        precursor = group.loc[group["phase"].eq("precursor"), "enrichment_delta"].to_numpy(dtype=float)
        if not len(late) or not len(precursor):
            continue
        ci_low, ci_high = _bootstrap_mean_difference(
            late,
            precursor,
            iterations=bootstrap_iterations,
            seed=seed + row_index,
        )
        records.append(
            {
                "axis_id": keys[0],
                "axis_label": keys[1],
                "evidence_type": keys[2],
                "analysis_level": "patient_aggregated" if aggregate_by_patient else "sample",
                "n_late": len(late),
                "n_precursor": len(precursor),
                "late_mean": late.mean(),
                "precursor_mean": precursor.mean(),
                "mean_difference": late.mean() - precursor.mean(),
                "ci_95_low": ci_low,
                "ci_95_high": ci_high,
                "late_median": np.median(late),
                "precursor_median": np.median(precursor),
                "mannwhitney_p": mannwhitneyu(late, precursor, alternative="two-sided").pvalue,
            }
        )
    summary = pd.DataFrame(records)
    if summary.empty:
        return summary
    summary["mannwhitney_q_bh"] = benjamini_hochberg(summary["mannwhitney_p"])
    return summary.sort_values(["mean_difference", "axis_id", "evidence_type"], ascending=[False, True, True])


def build_paired_patient_differences(adjacency: pd.DataFrame) -> pd.DataFrame:
    """Return paired patient-level late-minus-precursor enrichment differences."""
    valid = _valid_phase_rows(adjacency)
    group_columns = ["axis_id", "axis_label", "evidence_type"]
    patient_means = (
        valid.groupby(group_columns + ["patient_id", "phase"], as_index=False)
        .agg(enrichment_delta=("enrichment_delta", "mean"))
    )
    paired = (
        patient_means.pivot(
            index=group_columns + ["patient_id"],
            columns="phase",
            values="enrichment_delta",
        )
        .reset_index()
        .dropna(subset=["late", "precursor"])
    )
    paired["late_minus_precursor"] = paired["late"] - paired["precursor"]
    return paired.sort_values(group_columns + ["patient_id"]).reset_index(drop=True)


def summarize_paired_patient_differences(
    paired: pd.DataFrame,
    *,
    bootstrap_iterations: int = 2_000,
    seed: int = 307534,
) -> pd.DataFrame:
    """Summarize paired patient late-minus-precursor spatial enrichment changes."""
    group_columns = ["axis_id", "axis_label", "evidence_type"]
    records = []
    for row_index, (keys, group) in enumerate(paired.groupby(group_columns, sort=True)):
        differences = group["late_minus_precursor"].to_numpy(dtype=float)
        ci_low, ci_high = _bootstrap_mean(
            differences,
            iterations=bootstrap_iterations,
            seed=seed + row_index,
        )
        try:
            p_value = wilcoxon(differences, alternative="two-sided").pvalue
        except ValueError:
            p_value = np.nan
        records.append(
            {
                "axis_id": keys[0],
                "axis_label": keys[1],
                "evidence_type": keys[2],
                "n_paired_patients": len(differences),
                "paired_difference_mean": differences.mean(),
                "paired_difference_median": np.median(differences),
                "ci_95_low": ci_low,
                "ci_95_high": ci_high,
                "positive_fraction": (differences > 0).mean(),
                "wilcoxon_p": p_value,
            }
        )
    summary = pd.DataFrame(records)
    if summary.empty:
        return summary
    summary["wilcoxon_q_bh"] = benjamini_hochberg(summary["wilcoxon_p"])
    return summary.sort_values(
        ["paired_difference_mean", "axis_id", "evidence_type"],
        ascending=[False, True, True],
    )
