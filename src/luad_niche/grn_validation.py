"""Helpers for external validation of GRN-prioritized signatures."""

from __future__ import annotations

import pandas as pd


SIGNATURE_STATE_MAP = {
    "c1q_macrophage_vs_other_macrophage": "c1q_macrophage",
    "spp1_macrophage_vs_other_macrophage": "spp1_macrophage",
    "inflammatory_macrophage_vs_other_macrophage": "inflammatory_macrophage",
    "resident_macrophage_vs_other_macrophage": "resident_macrophage",
    "proliferating_epithelial_vs_other_epithelial": "proliferating_epithelial",
    "epithelial_progenitor_like_vs_other_epithelial": "epithelial_progenitor_like",
}

STATE_EXPECTED_CELLTYPE = {
    "c1q_macrophage": "Myeloid cells",
    "spp1_macrophage": "Myeloid cells",
    "inflammatory_macrophage": "Myeloid cells",
    "resident_macrophage": "Myeloid cells",
    "proliferating_epithelial": "Epithelial cells",
    "epithelial_progenitor_like": "Epithelial cells",
}


def canonical_signature_state(signature: str) -> str:
    """Map GRN signature names to state names used by external summaries."""
    text = str(signature)
    return SIGNATURE_STATE_MAP.get(text, text)


def expected_celltype_for_state(state: str) -> str:
    """Return the expected GSE131907 broad reference cell type for a state."""
    return STATE_EXPECTED_CELLTYPE.get(str(state), "")


def summarize_group_delta(
    table: pd.DataFrame,
    *,
    state_column: str,
    state_value: str,
    group_column: str,
    value_column: str,
    baseline_groups: list[str],
    comparison_groups: list[str],
    n_samples_column: str = "n_samples",
) -> dict[str, object]:
    """Summarize comparison-minus-baseline state abundance or score deltas."""
    required = {state_column, group_column, value_column}
    missing = required.difference(table.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"validation table is missing required columns: {missing_text}")

    working = table[table[state_column].astype(str).eq(str(state_value))].copy()
    working[value_column] = pd.to_numeric(working[value_column], errors="coerce")
    baseline = working[working[group_column].astype(str).isin([str(item) for item in baseline_groups])]
    comparison = working[working[group_column].astype(str).isin([str(item) for item in comparison_groups])]
    baseline_values = baseline[value_column].dropna()
    comparison_values = comparison[value_column].dropna()
    baseline_mean = float(baseline_values.mean()) if len(baseline_values) else float("nan")
    comparison_mean = float(comparison_values.mean()) if len(comparison_values) else float("nan")
    if n_samples_column in working.columns:
        baseline_n = int(pd.to_numeric(baseline[n_samples_column], errors="coerce").fillna(0).sum())
        comparison_n = int(pd.to_numeric(comparison[n_samples_column], errors="coerce").fillna(0).sum())
    else:
        baseline_n = int(len(baseline))
        comparison_n = int(len(comparison))
    return {
        "state": state_value,
        "baseline_groups": ",".join(baseline_groups),
        "comparison_groups": ",".join(comparison_groups),
        "baseline_mean": baseline_mean,
        "comparison_mean": comparison_mean,
        "delta": comparison_mean - baseline_mean,
        "baseline_n_samples": baseline_n,
        "comparison_n_samples": comparison_n,
    }
