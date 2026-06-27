"""Focused orthogonal-validation helpers for the early-LUAD niche project."""

from __future__ import annotations

from collections.abc import Iterable, Mapping

import numpy as np
import pandas as pd


SNRNA_STAGE_ORDER = ("Normal", "AAH", "AIS", "MIA", "LUAD")
BULK_STAGE_ORDER = ("Normal", "AIS", "MIA", "IAC")

FOCUSED_GENE_CONTEXTS = {
    "MIF": "epithelial",
    "CD74": "macrophage",
    "CD44": "macrophage",
    "CXCR4": "macrophage",
    "SPP1": "macrophage",
    "TREM2": "macrophage",
    "PLA2G7": "macrophage",
    "IL1B": "macrophage",
    "TNF": "macrophage",
    "CXCL8": "macrophage",
}


def _ordered_existing(values: Iterable[str], preferred: Iterable[str]) -> list[str]:
    observed = {str(value) for value in values}
    return [value for value in preferred if value in observed]


def _row_zscores(table: pd.DataFrame, value_column: str) -> pd.Series:
    """Return population row z-scores, keeping constant rows at zero."""
    numeric = pd.to_numeric(table[value_column], errors="coerce")
    means = numeric.groupby(table["gene"]).transform("mean")
    stds = numeric.groupby(table["gene"]).transform(lambda values: values.std(ddof=0))
    return ((numeric - means) / stds.replace(0, np.nan)).fillna(0.0)


def summarize_compartment_gene_expression(
    table: pd.DataFrame,
    *,
    sample: str,
    sample_name: str,
    stage: str,
    genes: Iterable[str],
) -> pd.DataFrame:
    """Summarize normalized selected-gene expression within broad compartments."""
    if "broad_class" not in table.columns:
        raise ValueError("table is missing required column: broad_class")

    contexts = {
        "all_cells": pd.Series(True, index=table.index),
        "epithelial": table["broad_class"].eq("epithelial"),
        "macrophage": table["broad_class"].eq("macrophage"),
    }
    records: list[dict[str, object]] = []
    for context, mask in contexts.items():
        subset = table.loc[mask]
        for gene in genes:
            present = gene in subset.columns
            values = pd.to_numeric(subset[gene], errors="coerce").fillna(0.0) if present else None
            records.append(
                {
                    "sample": sample,
                    "sample_name": sample_name,
                    "stage": stage,
                    "context": context,
                    "gene": gene,
                    "n_cells": len(subset),
                    "mean_expression": float(values.mean()) if present and len(subset) else None,
                    "detection_fraction": float(values.gt(0).mean()) if present and len(subset) else None,
                    "status": "present" if present else "missing",
                }
            )
    return pd.DataFrame(records)


def summarize_stage_expression(sample_summary: pd.DataFrame) -> pd.DataFrame:
    """Average compartment-level selected-gene expression across samples per stage."""
    required = {
        "sample",
        "stage",
        "context",
        "gene",
        "n_cells",
        "mean_expression",
        "detection_fraction",
    }
    missing = required.difference(sample_summary.columns)
    if missing:
        raise ValueError(f"sample summary is missing required columns: {', '.join(sorted(missing))}")

    working = sample_summary.copy()
    working["mean_expression"] = pd.to_numeric(working["mean_expression"], errors="coerce")
    working["detection_fraction"] = pd.to_numeric(working["detection_fraction"], errors="coerce")
    result = (
        working.groupby(["stage", "context", "gene"], dropna=False)
        .agg(
            n_samples=("sample", "nunique"),
            mean_expression=("mean_expression", "mean"),
            sd_expression=("mean_expression", "std"),
            mean_detection_fraction=("detection_fraction", "mean"),
            total_cells=("n_cells", "sum"),
        )
        .reset_index()
    )
    return result


def prepare_snrna_focused_stage_source(
    stage_summary: pd.DataFrame,
    *,
    gene_contexts: Mapping[str, str] = FOCUSED_GENE_CONTEXTS,
    stages: Iterable[str] = SNRNA_STAGE_ORDER,
) -> pd.DataFrame:
    """Select biologically matched snRNA compartments and add display z-scores."""
    required = {"stage", "context", "gene", "mean_expression"}
    missing = required.difference(stage_summary.columns)
    if missing:
        raise ValueError(f"snRNA stage summary is missing required columns: {', '.join(sorted(missing))}")

    stage_order = _ordered_existing(stage_summary["stage"], stages)
    records = []
    for gene, context in gene_contexts.items():
        subset = stage_summary[
            stage_summary["gene"].eq(gene) & stage_summary["context"].eq(context)
        ].copy()
        subset = subset[subset["stage"].isin(stage_order)].copy()
        subset["gene"] = pd.Categorical(subset["gene"], [gene], ordered=True)
        subset["stage"] = pd.Categorical(subset["stage"], stage_order, ordered=True)
        records.append(subset)
    if not records:
        return pd.DataFrame()

    result = pd.concat(records, ignore_index=True)
    result["gene"] = result["gene"].astype(str)
    result["stage"] = result["stage"].astype(str)
    gene_order = list(gene_contexts)
    result["_gene_order"] = result["gene"].map({gene: index for index, gene in enumerate(gene_order)})
    result["_stage_order"] = result["stage"].map({stage: index for index, stage in enumerate(stage_order)})
    result = result.sort_values(["_gene_order", "_stage_order"]).drop(columns=["_gene_order", "_stage_order"])
    result["row_zscore"] = _row_zscores(result, "mean_expression")
    return result.reset_index(drop=True)


def prepare_bulk_focused_stage_source(
    marker_table: pd.DataFrame,
    *,
    genes: Iterable[str] = tuple(FOCUSED_GENE_CONTEXTS),
    stages: Iterable[str] = BULK_STAGE_ORDER,
) -> pd.DataFrame:
    """Prepare focused bulk stage means, row z-scores, and late-versus-normal deltas."""
    required = {"gene", "stage", "mean_expression", "status"}
    missing = required.difference(marker_table.columns)
    if missing:
        raise ValueError(f"bulk marker table is missing required columns: {', '.join(sorted(missing))}")

    genes = list(genes)
    stage_order = list(stages)
    result = marker_table[
        marker_table["gene"].isin(genes)
        & marker_table["stage"].isin(stage_order)
        & marker_table["status"].eq("present")
    ].copy()
    result["mean_expression"] = pd.to_numeric(result["mean_expression"], errors="coerce")
    result["_gene_order"] = result["gene"].map({gene: index for index, gene in enumerate(genes)})
    result["_stage_order"] = result["stage"].map({stage: index for index, stage in enumerate(stage_order)})
    result = result.sort_values(["_gene_order", "_stage_order"]).drop(columns=["_gene_order", "_stage_order"])
    result["row_zscore"] = _row_zscores(result, "mean_expression")

    wide = result.pivot_table(index="gene", columns="stage", values="mean_expression", aggfunc="mean")
    first_stage = stage_order[0]
    last_stage = stage_order[-1]
    delta = wide.get(last_stage, pd.Series(dtype=float)) - wide.get(first_stage, pd.Series(dtype=float))
    result["delta_late_vs_normal"] = result["gene"].map(delta.to_dict())
    return result.reset_index(drop=True)
