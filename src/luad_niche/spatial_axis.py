"""Spatial summaries for candidate ligand-receptor mechanism axes."""

from __future__ import annotations

import pandas as pd


def _deduplicate(values: list[str]) -> list[str]:
    result: list[str] = []
    for value in values:
        if value and value not in result:
            result.append(value)
    return result


def build_axis_gene_panels(axes: list[dict]) -> dict[str, list[str]]:
    """Return source and target gene panels for each candidate mechanism axis."""
    panels: dict[str, list[str]] = {}
    for axis in axes:
        axis_id = axis.get("id", "")
        if not axis_id:
            continue
        source_genes = _deduplicate(list(axis.get("source_genes", []) or []))
        target_genes = _deduplicate(list(axis.get("target_genes", []) or []))
        if source_genes:
            panels[f"{axis_id}_source"] = source_genes
        if target_genes:
            panels[f"{axis_id}_target"] = target_genes
    return panels


def summarize_axis_adjacency_by_stage(adjacency: pd.DataFrame) -> pd.DataFrame:
    """Summarize per-sample axis adjacency while excluding invalid tests from effect means."""
    summary_input = adjacency.copy()
    if "status" not in summary_input.columns:
        summary_input["status"] = "ok"
    group_columns = ["axis_id", "evidence_type", "stage"]
    all_counts = (
        summary_input.groupby(group_columns, dropna=False)
        .agg(
            n_samples=("sample", "nunique"),
            n_tests=("sample", "size"),
            n_invalid_tests=("status", lambda values: int((values != "ok").sum())),
        )
        .reset_index()
    )
    valid = summary_input[summary_input["status"].eq("ok")].copy()
    valid_summary = (
        valid.groupby(group_columns, dropna=False)
        .agg(
            n_valid_samples=("sample", "nunique"),
            n_valid_tests=("sample", "size"),
            observed_fraction_mean=("observed_fraction", "mean"),
            null_mean_mean=("null_mean", "mean"),
            enrichment_delta_mean=("enrichment_delta", "mean"),
            empirical_p_greater_median=("empirical_p_greater", "median"),
            empirical_p_less_median=("empirical_p_less", "median"),
        )
        .reset_index()
    )
    summary = all_counts.merge(valid_summary, on=group_columns, how="left")
    for column in ["n_valid_samples", "n_valid_tests"]:
        summary[column] = summary[column].fillna(0).astype(int)
    return summary.sort_values(group_columns).reset_index(drop=True)
