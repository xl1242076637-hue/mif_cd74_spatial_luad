"""Lightweight GRN-level virtual perturbation utilities.

These functions implement a transparent, scTenifoldKnk-inspired analysis layer:
build a gene-gene coexpression network, remove outgoing edges from selected
genes, and quantify how much propagated network signal is lost.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def build_correlation_network(
    expr: pd.DataFrame,
    *,
    min_abs_correlation: float = 0.05,
) -> pd.DataFrame:
    """Build a row-normalized positive correlation network from cells x genes data."""
    numeric = expr.apply(pd.to_numeric, errors="coerce")
    numeric = numeric.dropna(axis=1, how="all").fillna(0.0)
    variable = numeric.loc[:, numeric.var(axis=0) > 0]
    if variable.shape[1] < 2:
        return pd.DataFrame()

    corr = variable.corr(method="pearson").fillna(0.0)
    corr.values[np.diag_indices_from(corr.values)] = 0.0
    corr = corr.where(corr >= min_abs_correlation, 0.0)
    row_sums = corr.sum(axis=1)
    corr = corr.loc[row_sums > 0, row_sums > 0]
    row_sums = corr.sum(axis=1)
    if corr.empty:
        return corr
    return corr.div(row_sums, axis=0)


def _propagate(
    network: pd.DataFrame,
    seed: pd.Series,
    *,
    n_steps: int,
    restart: float,
    decay: float,
) -> pd.Series:
    current = seed.astype(float).reindex(network.index).fillna(0.0)
    cumulative = pd.Series(0.0, index=network.columns)
    for step in range(1, n_steps + 1):
        current = restart * seed.reindex(network.index).fillna(0.0) + (1 - restart) * current.dot(network)
        cumulative = cumulative.add(current.reindex(cumulative.index).fillna(0.0) * (decay**step), fill_value=0.0)
    return cumulative


def virtual_outgoing_knockout(
    network: pd.DataFrame,
    ko_genes: list[str],
    *,
    n_steps: int = 3,
    restart: float = 0.15,
    decay: float = 0.5,
) -> pd.DataFrame:
    """Remove outgoing edges from KO genes and rank propagated signal losses."""
    if network.empty:
        return pd.DataFrame(
            columns=["gene", "baseline_signal", "knockout_signal", "impact_score", "impact_rank"]
        )
    ko_present = [gene for gene in dict.fromkeys(ko_genes) if gene in network.index]
    if not ko_present:
        return pd.DataFrame(
            columns=["gene", "baseline_signal", "knockout_signal", "impact_score", "impact_rank"]
        )

    seed = pd.Series(0.0, index=network.index)
    seed.loc[ko_present] = 1.0 / len(ko_present)
    baseline = _propagate(network, seed, n_steps=n_steps, restart=restart, decay=decay)
    knockout_network = network.copy()
    knockout_network.loc[ko_present, :] = 0.0
    knockout = _propagate(knockout_network, seed, n_steps=n_steps, restart=restart, decay=decay)
    genes = [gene for gene in network.columns if gene not in ko_present]
    effects = pd.DataFrame(
        {
            "gene": genes,
            "baseline_signal": baseline.reindex(genes).fillna(0.0).to_numpy(),
            "knockout_signal": knockout.reindex(genes).fillna(0.0).to_numpy(),
        }
    )
    effects["impact_score"] = (effects["baseline_signal"] - effects["knockout_signal"]).clip(lower=0.0)
    effects = effects.sort_values(["impact_score", "gene"], ascending=[False, True]).reset_index(drop=True)
    effects["impact_rank"] = np.arange(1, len(effects) + 1)
    return effects


def summarize_signature_impacts(
    effects: pd.DataFrame,
    signatures: dict[str, list[str]],
) -> pd.DataFrame:
    """Summarize perturbation impact scores over predefined gene signatures."""
    rows: list[dict[str, object]] = []
    impact_by_gene = effects.set_index("gene")["impact_score"] if not effects.empty else pd.Series(dtype=float)
    for signature, genes in signatures.items():
        present = [gene for gene in dict.fromkeys(genes) if gene in impact_by_gene.index]
        values = impact_by_gene.reindex(present).dropna()
        rows.append(
            {
                "signature": signature,
                "n_present_genes": int(len(values)),
                "present_genes": ",".join(present),
                "mean_impact_score": float(values.mean()) if len(values) else 0.0,
                "max_impact_score": float(values.max()) if len(values) else 0.0,
                "sum_impact_score": float(values.sum()) if len(values) else 0.0,
            }
        )
    return pd.DataFrame(rows)


def summarize_target_ranking_stability(detail: pd.DataFrame) -> pd.DataFrame:
    """Summarize target-rank stability across GRN perturbation parameter runs."""
    required = {
        "run_id",
        "broad_class",
        "target_gene",
        "display_rank",
        "top_impacted_signature",
        "top_signature_mean_impact",
    }
    missing = required.difference(detail.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"GRN robustness detail is missing required columns: {missing_text}")

    table = detail.copy()
    table["display_rank"] = pd.to_numeric(table["display_rank"], errors="coerce")
    table["top_signature_mean_impact"] = pd.to_numeric(
        table["top_signature_mean_impact"],
        errors="coerce",
    ).fillna(0.0)
    rows: list[dict[str, object]] = []
    for (broad_class, target_gene), group in table.groupby(["broad_class", "target_gene"], sort=False):
        ranks = group["display_rank"].dropna()
        signatures = group["top_impacted_signature"].fillna("").astype(str)
        counts = signatures.value_counts(dropna=False)
        top_signature = str(counts.index[0]) if not counts.empty else ""
        rows.append(
            {
                "broad_class": broad_class,
                "target_gene": target_gene,
                "n_runs": int(group["run_id"].nunique()),
                "median_rank": float(ranks.median()) if len(ranks) else float("nan"),
                "mean_rank": float(ranks.mean()) if len(ranks) else float("nan"),
                "rank_min": int(ranks.min()) if len(ranks) else 0,
                "rank_max": int(ranks.max()) if len(ranks) else 0,
                "rank_iqr": float(ranks.quantile(0.75) - ranks.quantile(0.25)) if len(ranks) else 0.0,
                "top_signature_mode": top_signature,
                "top_signature_mode_fraction": float(counts.iloc[0] / len(signatures)) if len(signatures) else 0.0,
                "n_unique_top_signatures": int(signatures.nunique()) if len(signatures) else 0,
                "median_top_signature_mean_impact": float(group["top_signature_mean_impact"].median()),
            }
        )
    summary = pd.DataFrame(rows)
    if summary.empty:
        return summary
    return summary.sort_values(
        ["median_rank", "rank_iqr", "target_gene"],
        ascending=[True, True, True],
    ).reset_index(drop=True)
