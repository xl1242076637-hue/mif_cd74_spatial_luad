"""Control-analysis helpers for spatial candidate-axis prioritization."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, t


DEFAULT_EXCLUDED_PREFIXES = ("MT-", "RPL", "RPS")


def select_expression_matched_controls(
    gene_summary: pd.DataFrame,
    *,
    target_gene: str,
    n_controls: int,
    min_reference_samples: int,
    excluded_genes: Iterable[str] = (),
    excluded_prefixes: tuple[str, ...] = DEFAULT_EXCLUDED_PREFIXES,
) -> pd.DataFrame:
    """Select deterministic expression-matched control genes for a target gene."""
    required = {"gene", "mean_expression", "n_reference_samples"}
    missing = required.difference(gene_summary.columns)
    if missing:
        raise ValueError(f"gene summary is missing required columns: {', '.join(sorted(missing))}")
    working = gene_summary.copy()
    working["mean_expression"] = pd.to_numeric(working["mean_expression"], errors="coerce")
    working["n_reference_samples"] = pd.to_numeric(working["n_reference_samples"], errors="coerce")
    target_rows = working[working["gene"].eq(target_gene)]
    if target_rows.empty:
        raise ValueError(f"target gene is absent from gene summary: {target_gene}")
    target_expression = float(target_rows.iloc[0]["mean_expression"])
    excluded = set(excluded_genes) | {target_gene}
    eligible = working[
        working["mean_expression"].gt(0)
        & working["n_reference_samples"].ge(min_reference_samples)
        & ~working["gene"].isin(excluded)
    ].copy()
    for prefix in excluded_prefixes:
        eligible = eligible[~eligible["gene"].astype(str).str.startswith(prefix)]
    eligible["target_gene"] = target_gene
    eligible["target_mean_expression"] = target_expression
    eligible["expression_log_distance"] = (
        np.log1p(eligible["mean_expression"]) - np.log1p(target_expression)
    ).abs()
    return (
        eligible.sort_values(["expression_log_distance", "gene"])
        .head(n_controls)
        .reset_index(drop=True)
    )


def summarize_random_control_distribution(
    paired_stats: pd.DataFrame,
    *,
    target_gene: str,
) -> pd.DataFrame:
    """Compare a target paired effect against expression-matched control genes."""
    required = {"axis_id", "paired_difference_mean"}
    missing = required.difference(paired_stats.columns)
    if missing:
        raise ValueError(f"paired stats are missing required columns: {', '.join(sorted(missing))}")
    working = paired_stats.copy()
    working["paired_difference_mean"] = pd.to_numeric(working["paired_difference_mean"], errors="coerce")
    target_rows = working[working["axis_id"].eq(target_gene)]
    if target_rows.empty:
        raise ValueError(f"target gene is absent from paired stats: {target_gene}")
    target_value = float(target_rows.iloc[0]["paired_difference_mean"])
    controls = working[~working["axis_id"].eq(target_gene)]["paired_difference_mean"].dropna().to_numpy(dtype=float)
    if not len(controls):
        raise ValueError("at least one random control gene is required")
    return pd.DataFrame(
        [
            {
                "target_gene": target_gene,
                "target_paired_difference_mean": target_value,
                "n_control_genes": len(controls),
                "control_mean": controls.mean(),
                "control_median": np.median(controls),
                "control_sd": controls.std(ddof=1) if len(controls) > 1 else 0.0,
                "target_percentile": (controls < target_value).mean(),
                "empirical_upper_p": (1 + (controls >= target_value).sum()) / (1 + len(controls)),
            }
        ]
    )


def spearman_association(table: pd.DataFrame, x_column: str, y_column: str) -> dict[str, object]:
    """Return a compact Spearman association summary for two numeric columns."""
    usable = table[[x_column, y_column]].apply(pd.to_numeric, errors="coerce").dropna()
    if len(usable) < 3:
        return {
            "x": x_column,
            "y": y_column,
            "n": len(usable),
            "spearman_r": float("nan"),
            "spearman_p": float("nan"),
        }
    result = spearmanr(usable[x_column], usable[y_column])
    return {
        "x": x_column,
        "y": y_column,
        "n": len(usable),
        "spearman_r": float(result.statistic),
        "spearman_p": float(result.pvalue),
    }


def _standardized_covariate(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    sd = numeric.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return numeric * 0.0
    return (numeric - numeric.mean()) / sd


def fit_ols_sensitivity(
    table: pd.DataFrame,
    *,
    outcome_column: str,
    effect_column: str,
    covariate_columns: Iterable[str] = (),
    model_id: str,
) -> dict[str, object]:
    """Fit a compact OLS sensitivity model and return the effect estimate.

    Covariates are z-scored for numerical stability. The primary effect column
    is kept on its original scale, so a binary late-versus-precursor coefficient
    remains directly interpretable as an adjusted mean difference.
    """
    covariate_columns = list(covariate_columns)
    required = {outcome_column, *covariate_columns}
    if effect_column != "intercept":
        required.add(effect_column)
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"sensitivity table is missing required columns: {', '.join(sorted(missing))}")

    columns = [outcome_column, *covariate_columns]
    if effect_column != "intercept":
        columns.append(effect_column)
    working = table[columns].copy()
    for column in columns:
        working[column] = pd.to_numeric(working[column], errors="coerce")
    working = working.dropna()

    y = working[outcome_column].to_numpy(dtype=float)
    design_parts = [np.ones(len(working), dtype=float)]
    design_names = ["intercept"]
    if effect_column != "intercept":
        design_parts.append(working[effect_column].to_numpy(dtype=float))
        design_names.append(effect_column)
    for column in covariate_columns:
        design_parts.append(_standardized_covariate(working[column]).to_numpy(dtype=float))
        design_names.append(column)
    x = np.column_stack(design_parts)
    beta = np.linalg.pinv(x) @ y
    fitted = x @ beta
    residuals = y - fitted
    n_obs = len(y)
    n_parameters = x.shape[1]
    residual_df = n_obs - n_parameters
    effect_index = design_names.index(effect_column)
    if residual_df > 0:
        sigma2 = float((residuals @ residuals) / residual_df)
        covariance = sigma2 * np.linalg.pinv(x.T @ x)
        effect_se = float(np.sqrt(max(covariance[effect_index, effect_index], 0.0)))
        effect_t = float(beta[effect_index] / effect_se) if effect_se > 0 else float("nan")
        effect_p = float(2 * t.sf(abs(effect_t), residual_df)) if np.isfinite(effect_t) else float("nan")
    else:
        effect_se = float("nan")
        effect_t = float("nan")
        effect_p = float("nan")

    return {
        "model_id": model_id,
        "outcome": outcome_column,
        "effect": effect_column,
        "covariates": ",".join(covariate_columns) if covariate_columns else "none",
        "n_observations": n_obs,
        "n_covariates": len(covariate_columns),
        "residual_df": residual_df,
        "effect_estimate": float(beta[effect_index]),
        "effect_se": effect_se,
        "effect_t": effect_t,
        "effect_p": effect_p,
    }


def build_paired_change_table(
    table: pd.DataFrame,
    *,
    value_column: str,
    covariate_columns: Iterable[str],
) -> pd.DataFrame:
    """Build paired patient late-minus-precursor outcome and covariate changes."""
    covariate_columns = list(covariate_columns)
    required = {"patient_id", "phase", value_column, *covariate_columns}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"paired-change table is missing required columns: {', '.join(sorted(missing))}")

    columns = ["patient_id", "phase", value_column, *covariate_columns]
    working = table[columns].copy()
    working = working[working["phase"].isin(["precursor", "late"])].copy()
    for column in [value_column, *covariate_columns]:
        working[column] = pd.to_numeric(working[column], errors="coerce")
    grouped = working.groupby(["patient_id", "phase"], as_index=False).mean(numeric_only=True)
    wide = grouped.pivot(index="patient_id", columns="phase", values=[value_column, *covariate_columns])
    wide.columns = [f"{metric}_{phase}" for metric, phase in wide.columns]
    required_wide = [f"{value_column}_late", f"{value_column}_precursor"]
    required_wide.extend(f"{column}_{phase}" for column in covariate_columns for phase in ("late", "precursor"))
    wide = wide.dropna(subset=required_wide).reset_index()
    result = pd.DataFrame({"patient_id": wide["patient_id"]})
    result[f"{value_column}_change"] = wide[f"{value_column}_late"] - wide[f"{value_column}_precursor"]
    for column in covariate_columns:
        result[f"{column}_change"] = wide[f"{column}_late"] - wide[f"{column}_precursor"]
    return result
