"""Utilities for in-silico gene and axis perturbation analyses."""

from __future__ import annotations

import pandas as pd


def apply_expression_perturbation(
    expr: pd.DataFrame,
    genes: list[str],
    factor: float,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Scale selected gene-expression columns and report present/missing genes."""
    requested = list(dict.fromkeys(genes))
    present = [gene for gene in requested if gene in expr.columns]
    missing = [gene for gene in requested if gene not in expr.columns]
    perturbed = expr.copy()
    for gene in present:
        perturbed[gene] = perturbed[gene] * factor
    return perturbed, {"present_genes": present, "missing_genes": missing, "factor": factor}


def summarize_perturbation_effects(effects: pd.DataFrame) -> pd.DataFrame:
    """Summarize per-sample perturbation effects by stage and evidence type."""
    group_columns = ["perturbation_id", "axis_id", "evidence_type", "stage"]
    summary = (
        effects.groupby(group_columns, dropna=False)
        .agg(
            n_samples=("sample", "nunique"),
            baseline_panel_mean=("baseline_panel_mean", "mean"),
            perturbed_panel_mean=("perturbed_panel_mean", "mean"),
            panel_mean_delta_mean=("panel_mean_delta", "mean"),
            baseline_observed_fraction=("baseline_observed_fraction", "mean"),
            perturbed_observed_fraction=("perturbed_observed_fraction", "mean"),
            observed_fraction_delta_mean=("observed_fraction_delta", "mean"),
        )
        .reset_index()
    )
    numeric_columns = summary.select_dtypes(include="number").columns
    summary[numeric_columns] = summary[numeric_columns].round(12)
    return summary.sort_values(group_columns).reset_index(drop=True)
