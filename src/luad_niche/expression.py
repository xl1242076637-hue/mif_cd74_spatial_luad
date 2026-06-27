"""Expression matrix helpers for bulk validation analyses."""

from __future__ import annotations

import numpy as np
import pandas as pd


STAGE_ORDER = ["Normal", "AIS", "MIA", "IAC"]


def expression_by_gene(df: pd.DataFrame, value_prefix: str = "FPKM.") -> pd.DataFrame:
    """Extract a gene-name-indexed expression matrix from a mixed GEO table."""
    value_columns = [column for column in df.columns if column.startswith(value_prefix)]
    if "gene_name" not in df.columns:
        raise ValueError("Input table must contain a 'gene_name' column.")
    if not value_columns:
        raise ValueError(f"No columns start with {value_prefix!r}.")

    expr = df[["gene_name", *value_columns]].copy()
    expr = expr.dropna(subset=["gene_name"])
    expr = expr.drop_duplicates(subset=["gene_name"], keep="first")
    expr = expr.set_index("gene_name")
    expr.columns = [column.removeprefix(value_prefix) for column in value_columns]
    return expr.apply(pd.to_numeric, errors="coerce")


def normalize_log1p_counts(expr: pd.DataFrame, scale_factor: float = 10_000.0) -> pd.DataFrame:
    """Normalize spot-by-gene counts by total counts per spot and log1p transform."""
    numeric = expr.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    totals = numeric.sum(axis=1)
    return normalize_log1p_counts_by_totals(numeric, totals, scale_factor=scale_factor)


def normalize_log1p_counts_by_totals(
    expr: pd.DataFrame,
    total_counts: pd.Series,
    scale_factor: float = 10_000.0,
) -> pd.DataFrame:
    """Normalize selected genes by externally supplied per-spot library sizes."""
    numeric = expr.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    totals = pd.Series(total_counts, index=numeric.index).astype(float)
    scaled = numeric.div(totals.replace(0, np.nan), axis=0) * scale_factor
    return np.log1p(scaled).fillna(0.0)


def progression_metadata(metadata: pd.DataFrame) -> pd.DataFrame:
    """Return metadata rows included in the LUAD progression trend."""
    meta = metadata.copy()
    include = meta["include_in_luad_progression"].astype(str).str.lower().isin(
        {"true", "1", "yes"}
    )
    meta = meta[include]
    meta = meta[meta["interpreted_stage"].isin(STAGE_ORDER)]
    return meta


def compute_stage_means(expr: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """Compute mean expression per stage for samples present in both inputs."""
    meta = progression_metadata(metadata)
    sample_to_stage = dict(zip(meta["title"], meta["interpreted_stage"], strict=False))
    samples = [sample for sample in expr.columns if sample in sample_to_stage]
    if not samples:
        raise ValueError("No expression columns matched metadata titles.")

    stage_means = {}
    for stage in STAGE_ORDER:
        stage_samples = [sample for sample in samples if sample_to_stage[sample] == stage]
        if stage_samples:
            stage_means[stage] = expr[stage_samples].mean(axis=1)
    return pd.DataFrame(stage_means)


def marker_stage_table(
    expr: pd.DataFrame,
    metadata: pd.DataFrame,
    markers: list[str],
) -> pd.DataFrame:
    """Return a long table of marker stage means and missing marker flags."""
    means = compute_stage_means(expr, metadata)
    records: list[dict] = []
    for marker in markers:
        if marker not in means.index:
            records.append({"gene": marker, "stage": "", "mean_expression": None, "status": "missing"})
            continue
        for stage in means.columns:
            records.append(
                {
                    "gene": marker,
                    "stage": stage,
                    "mean_expression": means.loc[marker, stage],
                    "status": "present",
                }
            )
    return pd.DataFrame(records)


def marker_trend_summary(marker_table: pd.DataFrame) -> pd.DataFrame:
    """Summarize Normal-to-IAC changes for present markers."""
    present = marker_table[marker_table["status"] == "present"].copy()
    present["mean_expression"] = pd.to_numeric(present["mean_expression"], errors="coerce")
    wide = present.pivot_table(
        index="gene",
        columns="stage",
        values="mean_expression",
        aggfunc="mean",
    )
    records: list[dict] = []
    for gene, row in wide.iterrows():
        stage_values = {stage: row.get(stage) for stage in STAGE_ORDER if stage in wide.columns}
        max_stage = max(
            (stage for stage, value in stage_values.items() if pd.notna(value)),
            key=lambda stage: stage_values[stage],
        )
        normal = row.get("Normal")
        iac = row.get("IAC")
        records.append(
            {
                "gene": gene,
                "normal_mean": normal,
                "iac_mean": iac,
                "delta_iac_vs_normal": iac - normal
                if pd.notna(normal) and pd.notna(iac)
                else None,
                "max_stage": max_stage,
            }
        )
    return pd.DataFrame(records).sort_values("delta_iac_vs_normal", ascending=False)
