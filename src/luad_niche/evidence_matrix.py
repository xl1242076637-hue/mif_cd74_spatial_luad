"""Build compact evidence matrices from multi-cohort LUAD niche outputs."""

from __future__ import annotations

import pandas as pd


def mean_axis_spatial_by_evidence(
    spatial: pd.DataFrame,
    stages: tuple[str, ...] = ("MIA", "LUAD"),
) -> pd.DataFrame:
    """Average candidate-axis spatial effects across selected progression stages."""
    subset = spatial[spatial["stage"].isin(stages)].copy()
    subset["enrichment_delta_mean"] = pd.to_numeric(subset["enrichment_delta_mean"], errors="coerce")
    summary = (
        subset.groupby(["axis_id", "evidence_type"], dropna=False)
        .agg(
            mean_axis_spatial_delta=("enrichment_delta_mean", "mean"),
            n_axis_spatial_stage_rows=("stage", "size"),
        )
        .reset_index()
    )
    summary["mean_axis_spatial_delta"] = summary["mean_axis_spatial_delta"].round(12)
    return summary


def best_continuous_perturbation_by_axis(ranking: pd.DataFrame) -> pd.DataFrame:
    """Select the strongest continuous perturbation result for each mechanism axis."""
    working = ranking.copy()
    working["continuous_priority_score"] = pd.to_numeric(
        working["continuous_priority_score"],
        errors="coerce",
    )
    working["coupling_relative_delta_mean"] = pd.to_numeric(
        working["coupling_relative_delta_mean"],
        errors="coerce",
    )
    working = working.sort_values(
        ["axis_id", "continuous_priority_score", "coupling_relative_delta_mean"],
        ascending=[True, False, True],
    )
    best = working.groupby("axis_id", as_index=False).head(1).reset_index(drop=True)
    return best.rename(
        columns={
            "perturbed_genes": "top_perturbed_genes",
            "perturbation_factor": "top_perturbation_factor",
            "evidence_type": "top_perturbation_evidence",
            "coupling_relative_delta_mean": "top_coupling_relative_delta",
            "continuous_priority_score": "top_continuous_priority_score",
        }
    )[
        [
            "axis_id",
            "top_perturbed_genes",
            "top_perturbation_factor",
            "top_perturbation_evidence",
            "top_coupling_relative_delta",
            "top_continuous_priority_score",
        ]
    ]
