"""Rank candidate epithelial-myeloid mechanisms from multi-cohort evidence."""

from __future__ import annotations

import math
from typing import Iterable

import pandas as pd


def _as_float(value) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return result


def positive_score(value, scale: float) -> float:
    """Map a positive effect size to 0-1 with clipping."""
    numeric = _as_float(value)
    if math.isnan(numeric) or scale <= 0:
        return 0.0
    return max(0.0, min(1.0, numeric / scale))


def gene_specificity_fraction(
    genes: list[str],
    top_celltypes: pd.DataFrame,
    expected_celltype: str,
) -> dict[str, object]:
    """Return the fraction of present genes whose top cell type matches expectation."""
    if not genes or not expected_celltype:
        return {"n_present": 0, "n_expected": 0, "fraction": float("nan"), "missing_genes": ",".join(genes)}
    by_gene = top_celltypes.set_index("gene", drop=False)
    present = []
    expected = []
    missing = []
    for gene in genes:
        if gene not in by_gene.index:
            missing.append(gene)
            continue
        row = by_gene.loc[gene]
        present.append(gene)
        if row.get("Cell_type.refined", "") == expected_celltype:
            expected.append(gene)
    fraction = len(expected) / len(present) if present else float("nan")
    return {
        "n_present": len(present),
        "n_expected": len(expected),
        "fraction": fraction,
        "missing_genes": ",".join(missing),
    }


def mean_bulk_delta(genes: list[str], trends: pd.DataFrame) -> dict[str, object]:
    by_gene = trends.set_index("gene", drop=False)
    values = []
    missing = []
    for gene in genes:
        if gene not in by_gene.index:
            missing.append(gene)
            continue
        value = _as_float(by_gene.loc[gene].get("delta_iac_vs_normal"))
        if not math.isnan(value):
            values.append(value)
    return {
        "n_present": len(values),
        "mean_delta": sum(values) / len(values) if values else float("nan"),
        "missing_genes": ",".join(missing),
    }


def mean_spatial_delta(
    spatial_summary: pd.DataFrame,
    target: str,
    stages: Iterable[str] = ("MIA", "LUAD"),
) -> float:
    subset = spatial_summary[
        spatial_summary["target"].eq(target) & spatial_summary["stage"].isin(list(stages))
    ].copy()
    if subset.empty:
        return float("nan")
    values = pd.to_numeric(subset["enrichment_delta_mean"], errors="coerce").dropna()
    return float(values.mean()) if len(values) else float("nan")


def mean_axis_spatial_delta(
    axis_spatial_summary: pd.DataFrame | None,
    axis_id: str,
    stages: Iterable[str] = ("MIA", "LUAD"),
) -> float:
    if axis_spatial_summary is None or axis_spatial_summary.empty:
        return float("nan")
    required = {"axis_id", "stage", "enrichment_delta_mean"}
    if not required.issubset(axis_spatial_summary.columns):
        return float("nan")
    subset = axis_spatial_summary[
        axis_spatial_summary["axis_id"].eq(axis_id) & axis_spatial_summary["stage"].isin(list(stages))
    ].copy()
    if subset.empty:
        return float("nan")
    values = pd.to_numeric(subset["enrichment_delta_mean"], errors="coerce").dropna()
    return float(values.mean()) if len(values) else float("nan")


def program_delta(
    summary: pd.DataFrame,
    context: str,
    score: str,
    early_stage: str,
    late_stage: str,
) -> float:
    subset = summary[summary["context"].eq(context) & summary["score"].eq(score)].copy()
    if subset.empty:
        return float("nan")
    by_stage = subset.set_index("stage", drop=False)
    if early_stage not in by_stage.index or late_stage not in by_stage.index:
        return float("nan")
    late = _as_float(by_stage.loc[late_stage].get("mean_score"))
    early = _as_float(by_stage.loc[early_stage].get("mean_score"))
    return late - early if not math.isnan(late) and not math.isnan(early) else float("nan")


def rank_candidate_axes(
    axes: list[dict],
    gene_top_celltypes: pd.DataFrame,
    bulk_trends: pd.DataFrame,
    spatial_summary: pd.DataFrame,
    snrna_summary: pd.DataFrame,
    scrna_summary: pd.DataFrame,
    axis_spatial_summary: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Rank candidate axes using transparent weighted evidence components."""
    rows = []
    use_axis_spatial = axis_spatial_summary is not None and not axis_spatial_summary.empty
    for axis in axes:
        source_specificity = gene_specificity_fraction(
            axis.get("source_genes", []),
            gene_top_celltypes,
            axis.get("source_expected_celltype", ""),
        )
        target_specificity = gene_specificity_fraction(
            axis.get("target_genes", []),
            gene_top_celltypes,
            axis.get("target_expected_celltype", ""),
        )
        bulk = mean_bulk_delta(axis.get("bulk_genes", []), bulk_trends)
        spatial_delta = mean_spatial_delta(spatial_summary, axis.get("spatial_target", ""))
        axis_spatial_delta = mean_axis_spatial_delta(axis_spatial_summary, axis.get("id", ""))
        snrna_delta = program_delta(
            snrna_summary,
            axis.get("snrna_context", ""),
            axis.get("snrna_score", ""),
            early_stage=axis.get("snrna_early_stage", "Normal"),
            late_stage=axis.get("snrna_late_stage", "LUAD"),
        )
        scrna_delta = program_delta(
            scrna_summary,
            axis.get("scrna_context", ""),
            axis.get("scrna_score", ""),
            early_stage=axis.get("scrna_early_stage", "Adjacent"),
            late_stage=axis.get("scrna_late_stage", "Tumor"),
        )
        source_score = 0.0 if math.isnan(source_specificity["fraction"]) else source_specificity["fraction"]
        target_score = 0.0 if math.isnan(target_specificity["fraction"]) else target_specificity["fraction"]
        bulk_score = positive_score(bulk["mean_delta"], scale=20.0)
        spatial_score = positive_score(spatial_delta, scale=0.10)
        axis_spatial_score = positive_score(axis_spatial_delta, scale=0.20)
        snrna_score = positive_score(snrna_delta, scale=0.10)
        scrna_score = positive_score(scrna_delta, scale=0.05)
        if use_axis_spatial:
            priority = (
                0.20 * spatial_score
                + 0.15 * axis_spatial_score
                + 0.17 * source_score
                + 0.12 * target_score
                + 0.13 * bulk_score
                + 0.13 * snrna_score
                + 0.10 * scrna_score
            )
        else:
            priority = (
                0.25 * spatial_score
                + 0.20 * source_score
                + 0.15 * target_score
                + 0.15 * bulk_score
                + 0.15 * snrna_score
                + 0.10 * scrna_score
            )
        rows.append(
            {
                "axis_id": axis.get("id", ""),
                "axis_label": axis.get("label", axis.get("id", "")),
                "priority_score": round(priority, 6),
                "spatial_mia_luad_delta": spatial_delta,
                "spatial_score": spatial_score,
                "axis_spatial_mia_luad_delta": axis_spatial_delta,
                "axis_spatial_score": axis_spatial_score,
                "source_specificity_fraction": source_specificity["fraction"],
                "source_specificity_score": source_score,
                "source_genes_present": source_specificity["n_present"],
                "source_genes_expected": source_specificity["n_expected"],
                "source_missing_genes": source_specificity["missing_genes"],
                "target_specificity_fraction": target_specificity["fraction"],
                "target_specificity_score": target_score,
                "target_genes_present": target_specificity["n_present"],
                "target_genes_expected": target_specificity["n_expected"],
                "target_missing_genes": target_specificity["missing_genes"],
                "bulk_delta_mean": bulk["mean_delta"],
                "bulk_score": bulk_score,
                "bulk_genes_present": bulk["n_present"],
                "bulk_missing_genes": bulk["missing_genes"],
                "snrna_program_delta": snrna_delta,
                "snrna_score_component": snrna_score,
                "scrna_tumor_adjacent_delta": scrna_delta,
                "scrna_score_component": scrna_score,
                "perturbation_genes": ",".join(axis.get("perturbation_genes", [])),
                "notes": axis.get("notes", ""),
            }
        )
    return pd.DataFrame(rows).sort_values("priority_score", ascending=False).reset_index(drop=True)


def summarize_perturbation_candidates(axes: list[dict], ranked_axes: pd.DataFrame) -> pd.DataFrame:
    """Summarize perturbation genes by the priorities of axes they appear in."""
    priority_by_axis = dict(zip(ranked_axes["axis_id"], ranked_axes["priority_score"], strict=False))
    records = []
    for axis in axes:
        axis_id = axis.get("id", "")
        priority = priority_by_axis.get(axis_id, 0.0)
        for gene in axis.get("perturbation_genes", []):
            records.append({"gene": gene, "axis_id": axis_id, "axis_priority": priority})
    if not records:
        return pd.DataFrame(columns=["gene", "n_axes", "max_axis_priority", "mean_axis_priority", "axis_ids"])
    table = pd.DataFrame(records)
    summary = (
        table.groupby("gene")
        .agg(
            n_axes=("axis_id", "nunique"),
            max_axis_priority=("axis_priority", "max"),
            mean_axis_priority=("axis_priority", "mean"),
            axis_ids=("axis_id", lambda values: ",".join(sorted(set(values)))),
        )
        .reset_index()
        .sort_values(["max_axis_priority", "n_axes", "gene"], ascending=[False, False, True])
        .reset_index(drop=True)
    )
    summary["max_axis_priority"] = summary["max_axis_priority"].round(6)
    summary["mean_axis_priority"] = summary["mean_axis_priority"].round(6)
    return summary
